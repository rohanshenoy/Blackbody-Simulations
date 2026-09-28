"""Per-angle field extraction, replicating legacy extract_waveguide_data / extract_far_field_data.

Scratch files live in the job's scratch directory and are deleted after parsing.
Every wait is bounded by ``ctx.timeout_s``.
"""
from __future__ import annotations

import logging
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from bbsim.config import ExitFieldConfig
from bbsim.geometry import FaceInfo
from bbsim.naming import frequency_label
from bbsim.readers import FieldFileError, read_exit_field_fld, read_far_field_ffd, wait_for_file
from bbsim.sampling import FarFieldGrid
from bbsim.schema import FAR_FIELD_COLUMNS, WAVEGUIDE_COLUMNS, order_columns

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class ExtractionContext:
    frequency_ghz: float
    solution: str                 # e.g. "500GHz : LastAdaptive"
    sphere_name: str
    phi_values: np.ndarray
    theta_values: np.ndarray
    grid: FarFieldGrid
    exit_field: ExitFieldConfig
    exit_face: FaceInfo
    exit_cs_name: str
    exit_cs_x: Sequence[float]
    exit_cs_y: Sequence[float]
    scratch_dir: Path
    timeout_s: float
    incoming_power_w: float

    @property
    def freq_label(self) -> str:
        return frequency_label(self.frequency_ghz)


def intrinsics(ephi: int, freq_label: str, phi_deg: float, theta_deg: float) -> list:
    return ["Ephi:=", ephi, "Freq:=", freq_label, "IWavePhi:=", f"{phi_deg}deg", "IWaveTheta:=", f"{theta_deg}deg"]


def intrinsic_variation_key(freq_label: str, phi_deg: float, theta_deg: float) -> str:
    return f"Freq='{freq_label}' IWavePhi='{phi_deg}deg' IWaveTheta='{theta_deg}deg'"


def local_bounds_mm(vertices_mm: Sequence[Sequence[float]], origin_mm: Sequence[float], cs_x: Sequence[float],
                    cs_y: Sequence[float], units: str) -> tuple[list[str], list[str]]:
    """Bounding box of the exit face expressed in the exit coordinate system (legacy vertex transform)."""
    x = np.asarray(cs_x, dtype=float)
    y = np.asarray(cs_y, dtype=float)
    basis = np.column_stack((x, y, np.cross(x, y)))
    to_local = np.linalg.inv(basis)
    origin = np.asarray(origin_mm, dtype=float)
    local = np.array([to_local @ (np.asarray(v, dtype=float) - origin) for v in vertices_mm])

    def fmt(value: float) -> str:
        return f"{round(float(value), 12) + 0.0}{units}"

    return [fmt(v) for v in local.min(axis=0)], [fmt(v) for v in local.max(axis=0)]


def evaluate_outgoing_power(hfss: Any, ctx: ExtractionContext, ephi: int) -> dict[tuple[float, float], float]:
    fields = hfss.odesign.GetModule("FieldsReporter")
    fields.CalcStack("clear")
    fields.CopyNamedExprToStack("outgoing_power")
    power: dict[tuple[float, float], float] = {}
    total = len(ctx.phi_values) * len(ctx.theta_values)
    for phi in ctx.phi_values:
        for theta in ctx.theta_values:
            args = intrinsics(ephi, ctx.freq_label, float(phi), float(theta))
            fields.ClcEval(ctx.solution, args, "Fields")
            result = fields.GetTopEntryValue(ctx.solution, args)
            fields.CalcStack("pop")
            key = (float(phi), float(theta))
            power[key] = float(result[0])
            log.info("[power %d/%d] phi=%s theta=%s Ephi=%s -> %.6e W", len(power), total, phi, theta, ephi, power[key])
    return power


def exit_field_grid_bounds(hfss: Any, ctx: ExtractionContext) -> tuple[list[str], list[str]]:
    units = hfss.modeler.model_units
    if ctx.exit_field.boundary_mm is not None:
        lo, hi = ctx.exit_field.boundary_mm
        return [f"{v}{units}" for v in lo], [f"{v}{units}" for v in hi]
    oeditor = hfss.oeditor
    vertex_ids = oeditor.GetVertexIDsFromFace(ctx.exit_face.id)
    vertices = [[float(c) for c in oeditor.GetVertexPosition(i)] for i in vertex_ids]
    return local_bounds_mm(vertices, ctx.exit_face.center_mm, ctx.exit_cs_x, ctx.exit_cs_y, units)


def extract_waveguide(hfss: Any, ctx: ExtractionContext, ephi: int) -> pd.DataFrame:
    power = evaluate_outgoing_power(hfss, ctx, ephi)

    fields = hfss.odesign.GetModule("FieldsReporter")
    fields.CalcStack("clear")
    fields.EnterQty("E")
    fields.CalcOp("Smooth")
    if ctx.exit_field.manual:
        range_min, range_max = exit_field_grid_bounds(hfss, ctx)
        resolution = [f"{v}mm" for v in ctx.exit_field.resolution_mm]
        log.info("exit-field grid in %s: min %s max %s step %s", ctx.exit_cs_name, range_min, range_max, resolution)
    else:
        fields.EnterSurf("outgoing")
        fields.CalcOp("Value")

    path = ctx.scratch_dir / f"exitfield_{ctx.freq_label}_Ephi{ephi}.fld"
    frames = []
    for phi in ctx.phi_values:
        for theta in ctx.theta_values:
            if path.exists():
                path.unlink()
            args = intrinsics(ephi, ctx.freq_label, float(phi), float(theta))
            if ctx.exit_field.manual:
                fields.ExportOnGrid(
                    str(path), range_min, range_max, resolution, ctx.solution, args + ["Phase:=", "0deg"],
                    ["NAME:ExportOption", "IncludePtInOutput:=", True, "RefCSName:=", ctx.exit_cs_name,
                     "PtInSI:=", True, "FieldInRefCS:=", False],
                    "Cartesian", ["0mm", "0mm", "0mm"], False,
                )
            else:
                fields.CalculatorWrite(str(path), ["Solution:=", ctx.solution], args)
            wait_for_file(path, ctx.timeout_s)
            df = read_exit_field_fld(path)
            path.unlink()
            if frames and len(df) != len(frames[0]):
                # Geant4 pairs exit points across polarizations by row index; every angle must have the same grid.
                raise FieldFileError(
                    f"exit field for phi={phi} theta={theta} Ephi={ephi} has {len(df)} points, "
                    f"but the first angle had {len(frames[0])}"
                )
            df["Freq"] = ctx.freq_label
            df["Ephi"] = ephi
            df["IWavePhi"] = float(phi)
            df["IWaveTheta"] = float(theta)
            df["OutgoingPower"] = power[(float(phi), float(theta))]
            frames.append(df)
            log.info("[exit %d/%d] phi=%s theta=%s Ephi=%s: %d points", len(frames), len(power), phi, theta, ephi, len(df))

    table = pd.concat(frames, ignore_index=True)
    table["IngoingPower"] = ctx.incoming_power_w
    return order_columns(table, WAVEGUIDE_COLUMNS)


def extract_far_field(hfss: Any, ctx: ExtractionContext, ephi: int) -> pd.DataFrame:
    radiation = hfss.odesign.GetModule("RadField")
    path = ctx.scratch_dir / f"farfield_{ctx.freq_label}_Ephi{ephi}.ffd"
    frames = []
    total = len(ctx.phi_values) * len(ctx.theta_values)
    for phi in ctx.phi_values:
        for theta in ctx.theta_values:
            if path.exists():
                path.unlink()
            radiation.ExportFieldsToFile([
                "ExportFileName:=", str(path),
                "SetupName:=", ctx.sphere_name,
                "IntrinsicVariationKey:=", intrinsic_variation_key(ctx.freq_label, float(phi), float(theta)),
                "DesignVariationKey:=", f"Ephi='{ephi}'",
                "SolutionName:=", ctx.solution,
                "Quantity:=", "",
            ])
            wait_for_file(path, ctx.timeout_s)
            df = read_far_field_ffd(path, ctx.grid.theta_step_deg, ctx.grid.phi_step_deg)
            path.unlink()
            df["Freq"] = ctx.freq_label
            df["Ephi"] = int(ephi)
            df["IWavePhi"] = float(phi)
            df["IWaveTheta"] = float(theta)
            frames.append(df)
            log.info("[far %d/%d] phi=%s theta=%s Ephi=%s: %d points", len(frames), total, phi, theta, ephi, len(df))
    return order_columns(pd.concat(frames, ignore_index=True), FAR_FIELD_COLUMNS)
