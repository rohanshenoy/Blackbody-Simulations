"""Per-angle field extraction, replicating legacy extract_waveguide_data / extract_far_field_data.

Scratch files live in the job's scratch directory and are deleted after parsing.
Every wait is bounded by ``ctx.timeout_s``. The exit face is described once per run
(``describe_exit_face``): its outline in the exit frame fixes the export lattice and the
cross-section declared in the dataset sidecar.
"""
from __future__ import annotations

import logging
import math
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from bbsim.config import ExitFieldConfig
from bbsim.geometry import FaceInfo, GeometryError
from bbsim.naming import frequency_label
from bbsim.readers import FieldFileError, read_calculator_scalar, read_exit_field_fld, read_far_field_ffd, wait_for_file
from bbsim.sampling import FarFieldGrid
from bbsim.schema import FAR_FIELD_COLUMNS, WAVEGUIDE_COLUMNS, order_columns

log = logging.getLogger(__name__)

GUARD_MM = 1e-9           # floating-noise guard when rounding lattice bounds outward
CONTAIN_REL_TOL = 1e-6    # retained exit points may exceed the declared cross-section by this fraction
PLANE_TOL_M = 1e-9        # retained exit points must lie on the exit plane x_e = 0 to this
RIM_BAND_REL = 0.05       # a disc's faceted rim: export disagreements tolerated within this fraction of the radius
FIELD_COLUMNS = ["Ex_real", "Ey_real", "Ez_real", "Ex_imag", "Ey_imag", "Ez_imag"]


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

    @property
    def swept_incidence(self) -> bool:
        """True when the plane wave has more than one incident direction; only then does AEDT define the
        incident-wave variables IWavePhi and IWaveTheta (see ``field_variation``)."""
        return len(self.phi_values) * len(self.theta_values) > 1


@dataclass(frozen=True)
class ExitFaceGeometry:
    """The exit-face outline in the exit frame, and what it implies for the export lattice and the sidecar."""

    local_points_mm: np.ndarray   # (N, 3): the face vertices, or edge samples when it has fewer than 3
    vertex_count: int
    bounds_method: str            # "vertices" | "edge_samples"
    cross_section: dict           # sidecar exit_field.cross_section, metres
    sagitta_mm: float             # chord error bound of the edge sampling; 0 for vertices


@dataclass(frozen=True)
class ExitGrid:
    """The ExportOnGrid lattice in the exit frame, in mm."""

    lo_mm: tuple[float, float, float]
    hi_mm: tuple[float, float, float]
    resolution_mm: tuple[float, float, float]

    @property
    def counts(self) -> tuple[int, ...]:
        return tuple(1 if step == 0 else int(math.floor((hi - lo) / step + 1e-9)) + 1
                     for lo, hi, step in zip(self.lo_mm, self.hi_mm, self.resolution_mm))

    @property
    def lattice_points(self) -> int:
        return math.prod(self.counts)

    def range_strings(self, units: str) -> tuple[list[str], list[str]]:
        return [_fmt_length(v, units) for v in self.lo_mm], [_fmt_length(v, units) for v in self.hi_mm]

    def axis_record(self) -> dict[str, dict[str, float | int]]:
        """Sidecar exit_field.grid: the lattice along exit Y and Z in metres (the CSV unit)."""
        axes = zip(("y_e", "z_e"), self.lo_mm[1:], self.hi_mm[1:], self.counts[1:])
        return {name: {"min": round(lo, 12) * 1e-3 + 0.0, "max": round(hi, 12) * 1e-3 + 0.0, "count": n}
                for name, lo, hi, n in axes}


def intrinsics(ephi: int, freq_label: str, phi_deg: float, theta_deg: float) -> list:
    return ["Ephi:=", ephi, "Freq:=", freq_label, "IWavePhi:=", f"{phi_deg}deg", "IWaveTheta:=", f"{theta_deg}deg"]


def field_variation(ctx: ExtractionContext, ephi: int, phi_deg: float, theta_deg: float) -> list:
    """Variation arguments of a calculator evaluation or a grid export at one incidence key and polarization.

    The polarization goes to AEDT as text ("1"), as PyAEDT passes design variables: over gRPC, CalculatorWrite
    ignored an integer Ephi and silently evaluated the nominal Ephi=0, so the Ephi=1 table of HPC step 5
    (job 4069263) carried Ephi=0's OutgoingPower while ExportOnGrid honoured the same integer.

    AEDT defines the incident-wave variables IWavePhi and IWaveTheta only when the plane wave sweeps more
    than one direction. With a single direction (PhiPoints = ThetaPoints = 1) naming them fails every
    evaluation and export, a constant included (HPC step 4a with --solve, 2026-10-06); the one direction is
    then implied. The legacy script always swept 15 directions. A sweep along one axis only is unverified.
    """
    args = ["Ephi:=", str(int(ephi)), "Freq:=", ctx.freq_label]
    if ctx.swept_incidence:
        args += ["IWavePhi:=", f"{phi_deg}deg", "IWaveTheta:=", f"{theta_deg}deg"]
    return args


def intrinsic_variation_key(freq_label: str, phi_deg: float, theta_deg: float) -> str:
    return f"Freq='{freq_label}' IWavePhi='{phi_deg}deg' IWaveTheta='{theta_deg}deg'"


def _fmt_length(value: float, units: str) -> str:
    return f"{round(float(value), 12) + 0.0}{units}"


def to_exit_frame(points_mm: Sequence[Sequence[float]], origin_mm: Sequence[float], cs_x: Sequence[float],
                  cs_y: Sequence[float]) -> np.ndarray:
    """Global points (mm) expressed in the exit frame (x = cs_x, y = cs_y, z = cs_x x cs_y), as (N, 3) rows."""
    x = np.asarray(cs_x, dtype=float)
    y = np.asarray(cs_y, dtype=float)
    to_local = np.linalg.inv(np.column_stack((x, y, np.cross(x, y))))
    origin = np.asarray(origin_mm, dtype=float)
    return np.array([to_local @ (np.asarray(p, dtype=float) - origin) for p in points_mm], dtype=float).reshape(-1, 3)


def local_bounds_mm(vertices_mm: Sequence[Sequence[float]], origin_mm: Sequence[float], cs_x: Sequence[float],
                    cs_y: Sequence[float], units: str) -> tuple[list[str], list[str]]:
    """Bounding box of the exit face expressed in the exit coordinate system (legacy vertex transform)."""
    local = to_exit_frame(vertices_mm, origin_mm, cs_x, cs_y)
    return [_fmt_length(v, units) for v in local.min(axis=0)], [_fmt_length(v, units) for v in local.max(axis=0)]


def round_bounds_outward(lo: Sequence[float], hi: Sequence[float], steps: Sequence[float],
                         guard_mm: float = GUARD_MM) -> tuple[list[float], list[float]]:
    """Widen each bound to the next multiple of its step; axes with step 0 keep their bounds.

    Never insets: points on the face boundary stay on the lattice, where they carry weight in BBRsim.
    The guard keeps floating noise (0.05000000000000002, or 0.025 / 0.001 = 25.000000000000004)
    from adding a whole step.
    """
    out_lo, out_hi = [], []
    for a, b, step in zip(lo, hi, steps):
        a, b, step = float(a), float(b), float(step)
        if step == 0:
            out_lo.append(a)
            out_hi.append(b)
        else:
            out_lo.append(math.floor((a + guard_mm) / step) * step)
            out_hi.append(math.ceil((b - guard_mm) / step) * step)
    return out_lo, out_hi


def classify_cross_section(local_points_mm: np.ndarray, vertex_count: int, rel_tol: float = 1e-6) -> dict:
    """Sidecar cross-section (metres) from the exit-face outline in the exit frame.

    rectangle: 4 vertices with two distinct Y and two distinct Z values; disc: fewer than 3 vertices and
    every sample equidistant from the face centre; polygon otherwise.
    """
    yz = np.asarray(local_points_mm, dtype=float)[:, 1:3]
    if vertex_count == 4:
        ys = np.unique(np.round(yz[:, 0], 9))
        zs = np.unique(np.round(yz[:, 1], 9))
        if len(ys) == 2 and len(zs) == 2:
            return {"shape": "rectangle", "y_e_half_m": float(ys[1] - ys[0]) / 2 * 1e-3,
                    "z_e_half_m": float(zs[1] - zs[0]) / 2 * 1e-3}
    if vertex_count < 3 and len(yz):
        r = np.hypot(yz[:, 0], yz[:, 1])
        mean = float(r.mean())
        if mean > 0 and bool(np.all(np.abs(r - mean) <= rel_tol * mean)):
            return {"shape": "disc", "radius_m": mean * 1e-3}
    return {"shape": "polygon", "vertices_m": [[float(y) * 1e-3, float(z) * 1e-3] for y, z in yz]}


def chord_sagitta_mm(local_points_mm: np.ndarray, samples_per_edge: int) -> float:
    """Upper bound on how far a chord between neighbouring samples of a circular edge falls inside the edge."""
    yz = np.asarray(local_points_mm, dtype=float)[:, 1:3]
    d_max = float(np.max(np.hypot(yz[:, 0], yz[:, 1]))) if len(yz) else 0.0
    return d_max * (1.0 - math.cos(math.pi / samples_per_edge))


def exit_grid(face: ExitFaceGeometry, exit_field: ExitFieldConfig) -> ExitGrid:
    """The export lattice. Manual bounds win; vertex bounds are used as they are (legacy, so the crack
    exports stay byte-identical); edge-sample bounds are rounded outward to the resolution."""
    resolution = tuple(float(v) for v in exit_field.resolution_mm)
    if exit_field.boundary_mm is not None:
        lo, hi = exit_field.boundary_mm
        return ExitGrid(tuple(float(v) for v in lo), tuple(float(v) for v in hi), resolution)
    lo, hi = face.local_points_mm.min(axis=0), face.local_points_mm.max(axis=0)
    if face.bounds_method == "edge_samples":
        lo, hi = round_bounds_outward(lo, hi, resolution)
    return ExitGrid(tuple(float(v) for v in lo), tuple(float(v) for v in hi), resolution)


def _outside_outline_hull(y: np.ndarray, z: np.ndarray, vertices_m: Sequence[Sequence[float]],
                          tol_m: float) -> np.ndarray:
    """True where a point lies farther than ``tol_m`` outside the convex hull of a polygon outline.

    The hull equals the face for convex outlines and contains it otherwise, so a point on the face is
    never rejected; the outline's point order does not matter.
    """
    from scipy.spatial import ConvexHull

    try:
        hull = ConvexHull(np.asarray(vertices_m, dtype=float))
    except Exception as exc:  # noqa: BLE001 - degenerate outline (fewer than 3 points, or collinear)
        raise FieldFileError(f"cannot test containment in a degenerate polygon outline: {exc}") from None
    # Rows of hull.equations are [unit outward normal, offset]; inside means normal . p + offset <= 0.
    distance = np.column_stack((y, z)) @ hull.equations[:, :2].T + hull.equations[:, 2]
    return np.any(distance > tol_m, axis=1)


def check_points_within(df: pd.DataFrame, cross_section: dict, where: str,
                        rel_tol: float = CONTAIN_REL_TOL, plane_tol_m: float = PLANE_TOL_M,
                        polygon_tol_m: float = 0.0) -> None:
    """Every retained exit point (metres, exit frame) must lie on x_e = 0 and inside the declared cross-section.

    Catches an HFSS that writes values (for example zeros) instead of nan for lattice points outside the solid.
    ``polygon_tol_m`` widens a polygon outline built from edge samples, whose chords cut inside the true edge.
    """
    x, y, z = (df[c].to_numpy(dtype=float) for c in ("X", "Y", "Z"))
    if len(x) and float(np.max(np.abs(x))) > plane_tol_m:
        raise FieldFileError(f"{where}: exit points off the plane x_e = 0 (max |X| = {np.max(np.abs(x)):.3e} m)")
    shape = cross_section["shape"]
    if shape == "disc":
        outside = y * y + z * z > cross_section["radius_m"] ** 2 * (1 + rel_tol)
    elif shape == "rectangle":
        outside = ((np.abs(y) > cross_section["y_e_half_m"] * (1 + rel_tol))
                   | (np.abs(z) > cross_section["z_e_half_m"] * (1 + rel_tol)))
    elif shape == "polygon":
        outline = np.asarray(cross_section["vertices_m"], dtype=float)
        size = float(np.max(np.ptp(outline, axis=0))) if len(outline) else 0.0
        outside = _outside_outline_hull(y, z, outline, rel_tol * size + polygon_tol_m)
    else:
        raise FieldFileError(f"{where}: unknown cross-section shape {shape!r}")
    if np.any(outside):
        raise FieldFileError(f"{where}: {int(np.sum(outside))} exit points lie outside the declared {shape} "
                             "cross-section; HFSS must write nan, not values, outside the solid")


def retain_disc_lattice(df: pd.DataFrame, cross_section: dict, where: str, band_rel: float = RIM_BAND_REL,
                        rel_tol: float = CONTAIN_REL_TOL) -> tuple[pd.DataFrame, dict]:
    """Every lattice point inside a declared disc, in export order, with zero field where HFSS gave none.

    ``df`` is the whole lattice as read with ``keep_unsolved``. HFSS meshes a circle as a polygon (the runner
    keeps curvilinear elements off, as in the crack reference designs), so at the rim the export disagrees with
    the declared circle, and differently for each polarization's adaptive mesh: HPC step 7a (2026-10-06, R = 50
    um at 2000 GHz) left about 60 points inside the disc without field, up to 0.76 um in, and gave field to
    about 50 outside it, up to 1.08 um out. Rim rule (Rohan, 2026-10-06, ledger HF-029): every lattice point
    inside the disc is kept, so all keys and both polarizations share one grid, which Geant4 pairs by row;
    valued points outside are dropped; inside points without field get zero field, which BBRsim never
    samples. A valued point more than ``band_rel`` of the radius beyond the rim, or an unsolved one more than
    that inside it, is an error: values written outside the solid, or a cross-section that does not match it.
    Returns the kept rows and the counts and extreme distances, in metres, for the manifest.
    """
    x, y, z = (df[c].to_numpy(dtype=float) for c in ("X", "Y", "Z"))
    if len(x) and float(np.max(np.abs(x))) > PLANE_TOL_M:
        raise FieldFileError(f"{where}: exit points off the plane x_e = 0 (max |X| = {np.max(np.abs(x)):.3e} m)")
    radius = float(cross_section["radius_m"])
    band = band_rel * radius
    r = np.hypot(y, z)
    inside = y * y + z * z <= radius ** 2 * (1 + rel_tol)
    solved = df[FIELD_COLUMNS].notna().all(axis=1).to_numpy()
    spill, holes = solved & ~inside, ~solved & inside
    beyond, depth = r[spill] - radius, np.maximum(radius - r[holes], 0.0)
    if np.any(beyond > band):
        raise FieldFileError(
            f"{where}: {int(np.sum(beyond > band))} points with field lie outside the declared disc more than "
            f"{band:.3g} m ({band_rel:.0%} of the radius) beyond its rim, the farthest {beyond.max():.3g} m; "
            "HFSS must write nan, not values, outside the solid")
    if np.any(depth > band):
        raise FieldFileError(
            f"{where}: {int(np.sum(depth > band))} points inside the declared disc have no field more than "
            f"{band:.3g} m ({band_rel:.0%} of the radius) inside its rim, the deepest {depth.max():.3g} m; "
            "the export does not cover the exit face")
    kept = df[inside].reset_index(drop=True)
    kept.loc[holes[inside], FIELD_COLUMNS] = 0.0
    return kept, {"unsolved_inside_zeroed": int(holes.sum()), "valued_outside_dropped": int(spill.sum()),
                  "deepest_unsolved_inside_m": float(depth.max()) if depth.size else 0.0,
                  "farthest_valued_outside_m": float(beyond.max()) if beyond.size else 0.0}


def describe_exit_face(hfss: Any, ctx: ExtractionContext) -> ExitFaceGeometry:
    """Outline of the exit face in the exit frame: its vertices, or samples along its edges when it has
    fewer than 3 vertices (a circle has none or one)."""
    oeditor = hfss.oeditor
    face_id = ctx.exit_face.id
    vertex_ids = list(oeditor.GetVertexIDsFromFace(face_id))
    if len(vertex_ids) >= 3:
        points = [[float(c) for c in oeditor.GetVertexPosition(v)] for v in vertex_ids]
        method = "vertices"
    else:
        n = ctx.exit_field.edge_samples
        points = [[float(c) for c in oeditor.GetEdgePositionAtNormalizedParameter(int(edge), k / n)]
                  for edge in oeditor.GetEdgeIDsFromFace(face_id) for k in range(n)]
        method = "edge_samples"
    if not points:
        raise GeometryError(f"exit face {face_id} has neither vertices nor edges to bound the export grid")
    local = to_exit_frame(points, ctx.exit_face.center_mm, ctx.exit_cs_x, ctx.exit_cs_y)
    sagitta = chord_sagitta_mm(local, ctx.exit_field.edge_samples) if method == "edge_samples" else 0.0
    positive = [s for s in ctx.exit_field.resolution_mm if s > 0]
    if positive and sagitta > 0.5 * min(positive):
        log.warning("edge-sampling sagitta %.3g mm exceeds half the smallest grid step %.3g mm; "
                    "raise exit_field.edge_samples", sagitta, min(positive))
    cross_section = classify_cross_section(local, len(vertex_ids))
    sampling = f", chord sagitta {sagitta:.3g} mm" if method == "edge_samples" else ""
    log.info("exit face %s: %d vertices, bounds from %s, cross-section %s%s", face_id, len(vertex_ids), method,
             cross_section["shape"], sampling)
    return ExitFaceGeometry(local, len(vertex_ids), method, cross_section, sagitta)


def evaluate_outgoing_power(hfss: Any, ctx: ExtractionContext, ephi: int) -> dict[tuple[float, float], float]:
    """The named expression ``outgoing_power`` (real Poynting flux through the exit face) per incident direction.

    Each value is written to a scratch file with ``CalculatorWrite`` and read back, as PyAEDT 1.7.0 evaluates
    expressions, with ``Phase`` among the intrinsics. The legacy script read the calculator stack instead
    (``ClcEval`` then ``GetTopEntryValue``); over gRPC on HPC that failed in every form tried, each time
    with the incident-angle keys of a single-direction plane wave, which ``field_variation`` now omits.
    """
    fields = hfss.odesign.GetModule("FieldsReporter")
    path = ctx.scratch_dir / f"outgoing_power_{ctx.freq_label}_Ephi{ephi}.fld"
    power: dict[tuple[float, float], float] = {}
    total = len(ctx.phi_values) * len(ctx.theta_values)
    for phi in ctx.phi_values:
        for theta in ctx.theta_values:
            if path.exists():
                path.unlink()
            args = field_variation(ctx, ephi, float(phi), float(theta)) + ["Phase:=", "0deg"]
            fields.CalcStack("clear")
            fields.CopyNamedExprToStack("outgoing_power")
            fields.CalculatorWrite(str(path), ["Solution:=", ctx.solution], args)
            wait_for_file(path, ctx.timeout_s)
            value = read_calculator_scalar(path)
            path.unlink()
            key = (float(phi), float(theta))
            power[key] = value
            log.info("[power %d/%d] phi=%s theta=%s Ephi=%s -> %.6e W", len(power), total, phi, theta, ephi, value)
    fields.CalcStack("clear")
    return power


def extract_waveguide(hfss: Any, ctx: ExtractionContext, ephi: int, face: ExitFaceGeometry,
                      rim_log: dict | None = None) -> pd.DataFrame:
    """The waveguide table of one polarization. For a disc, ``retain_disc_lattice`` keeps the whole lattice inside
    it; its per-key counts are appended to ``rim_log[str(ephi)]`` for the manifest."""
    power = evaluate_outgoing_power(hfss, ctx, ephi)

    fields = hfss.odesign.GetModule("FieldsReporter")
    fields.CalcStack("clear")
    fields.EnterQty("E")
    fields.CalcOp("Smooth")
    if ctx.exit_field.manual:
        grid = exit_grid(face, ctx.exit_field)
        range_min, range_max = grid.range_strings(hfss.modeler.model_units)
        resolution = [f"{v}mm" for v in ctx.exit_field.resolution_mm]
        log.info("exit-field grid in %s: min %s max %s step %s (%d lattice points, %s bounds)", ctx.exit_cs_name,
                 range_min, range_max, resolution, grid.lattice_points, face.bounds_method)
    else:
        fields.EnterSurf("outgoing")
        fields.CalcOp("Value")

    path = ctx.scratch_dir / f"exitfield_{ctx.freq_label}_Ephi{ephi}.fld"
    disc = ctx.exit_field.manual and face.cross_section["shape"] == "disc"
    frames = []
    for phi in ctx.phi_values:
        for theta in ctx.theta_values:
            if path.exists():
                path.unlink()
            args = field_variation(ctx, ephi, float(phi), float(theta))
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
            df = read_exit_field_fld(path, keep_unsolved=disc)
            path.unlink()
            where = f"exit field for phi={phi} theta={theta} Ephi={ephi}"
            if disc:
                df, rim = retain_disc_lattice(df, face.cross_section, where)
                log.info("%s: rim: %d unsolved inside points set to zero (deepest %.3g m), %d valued outside points "
                         "dropped (farthest %.3g m)", where, rim["unsolved_inside_zeroed"], rim["deepest_unsolved_inside_m"],
                         rim["valued_outside_dropped"], rim["farthest_valued_outside_m"])
                if rim_log is not None:
                    rim_log.setdefault(str(ephi), []).append({"key": [float(phi), float(theta)], **rim})
            check_points_within(df, face.cross_section, where, polygon_tol_m=face.sagitta_mm * 1e-3)
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
