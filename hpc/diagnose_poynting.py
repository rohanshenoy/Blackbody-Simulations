#!/usr/bin/env python3
"""Power diagnostic on a solved job, without a new solve (Codex queue CDX-009, 013, 014, 015; ledger G4-020, OR-019).

The tables give T = OutgoingPower / IngoingPower, where IngoingPower = A Ei^2 / (2 c mu0) carries no angle and
OutgoingPower integrates |Re S| over the exit face. This script measures, on a copy of a solved job's project,
the quantities that decide what that ratio means, per incidence key and polarization:

  exit face (total fields):      the integral of |Re S| (must reproduce the recorded OutgoingPower), of the axial
                                 flux Re S_z, of |Re S_z|, and of Re S . n with AEDT's own face normal;
  entrance face:                 the axial flux Re S_z with total fields (the power that enters the gap; it must
                                 equal the exit flux, the walls being PEC) and with scattered fields (the power the
                                 scattered field carries back out), with |Re S| and Re S . n alongside;
  incident field:                |E_inc| = |E_total - E_scattered| on a short line just inside the entrance, which
                                 IngoingPower assumes is Ei = 1 V/m at every angle (step 5b checked its direction);
  far field:                     the power radiated by the exit face's equivalent currents over the whole sphere,
                                 outward (x_e > 0) and backward hemispheres separately, with total and with scattered
                                 fields, on a new infinite sphere that covers both hemispheres.

The canonical pose is assumed and checked: propagation along global +z, entrance outward normal -z, exit +z, so
S_z is the flux into the gap at the entrance and out of it at the exit.

Run on HPC in a batch job after `source /home/rshenoy/BBRSim/bb_env.sh` (an AEDT session, no solver licence):

    python hpc/diagnose_poynting.py --job /home/rshenoy/BBRSim/outputs/InfParallelPlate_crack1Rohan_500GHz/job_tight1

It copies the job's project and results folder to <job>/diagnose_poynting/project/, opens the copy, never saves it,
writes <job>/diagnose_poynting/poynting.json, prints a table and a line starting "POYNTING:", and removes the copy
after a successful run (--keep-copy keeps it). The job's own files are not touched. Exit 0 when every primary
quantity was measured and the exit |Re S| integral reproduces the recorded OutgoingPower of every key to 1e-4;
1 otherwise (the JSON says what is missing); 2 when nothing was started.
"""
from __future__ import annotations

import argparse
import inspect
import json
import math
import os
import shutil
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bbsim.extract import field_variation, intrinsic_variation_key  # noqa: E402
from bbsim.hfss_setup import EXIT_CS, EXIT_FACE_LIST, far_field_sphere_kwargs  # noqa: E402
from bbsim.manifest import git_commit, sha256_file, write_json  # noqa: E402
from bbsim.readers import read_calculator_scalar, read_exit_field_fld, read_far_field_ffd, wait_for_file  # noqa: E402
from bbsim.sampling import SPEED_OF_LIGHT_M_PER_S, VACUUM_PERMEABILITY_H_PER_M, FarFieldGrid  # noqa: E402

ETA0_OHM = SPEED_OF_LIGHT_M_PER_S * VACUUM_PERMEABILITY_H_PER_M
WORK_DIR = "diagnose_poynting"
ENTRANCE_FACE_LIST = "incoming"
SPHERE = "dg_full_sphere"
CSV_TOLERANCE = 1e-4
INSIDE_MM = 0.001        # the incident-field line lies this far inside the entrance face
FIELD_COLUMNS = ["Ex_real", "Ey_real", "Ez_real", "Ex_imag", "Ey_imag", "Ez_imag"]

# Named calculator expressions: (calculator operations on Vector_RealPoynting, surface, primary). The order of the
# operations mirrors the runner's outgoing_power (named vector, operations, surface, Integrate); PyAEDT's own
# Power_Flow catalogue entry takes Re S . n with Normal then Dot.
EXIT_EXPRESSIONS = {
    "dg_exit_mag": (("Mag",), EXIT_FACE_LIST, True),
    "dg_exit_flux_z": (("ScalarZ",), EXIT_FACE_LIST, True),
    "dg_exit_abs_flux_z": (("ScalarZ", "Abs"), EXIT_FACE_LIST, True),
    "dg_exit_normal_flux": (("Normal", "Dot"), EXIT_FACE_LIST, False),
}
ENTRANCE_EXPRESSIONS = {
    "dg_entrance_flux_z": (("ScalarZ",), ENTRANCE_FACE_LIST, True),
    "dg_entrance_normal_flux": (("Normal", "Dot"), ENTRANCE_FACE_LIST, False),
    "dg_entrance_mag": (("Mag",), ENTRANCE_FACE_LIST, False),
}
EXIT_FIELDS = {"dg_exit_mag": "mag_w", "dg_exit_flux_z": "flux_z_w", "dg_exit_abs_flux_z": "abs_flux_z_w",
               "dg_exit_normal_flux": "flux_normal_w"}
ENTRANCE_FIELDS = {"dg_entrance_flux_z": "flux_z_w", "dg_entrance_normal_flux": "flux_normal_w",
                   "dg_entrance_mag": "mag_w"}
FIELD_TYPES = ("TotalFields", "ScatteredFields")


class Refused(RuntimeError):
    """Nothing was started: bad arguments, a missing results folder, or an existing work directory."""


# --- pure helpers ---------------------------------------------------------------------------------------------------

def parse_keys(items: list[str]) -> list[tuple[float, float]]:
    """``["0,180", "90,135"]`` -> ``[(0.0, 180.0), (90.0, 135.0)]`` (IWavePhi, IWaveTheta in degrees)."""
    keys = []
    for item in items:
        parts = item.split(",")
        if len(parts) != 2:
            raise ValueError(f"a key is 'phi,theta' in degrees, got {item!r}")
        keys.append((float(parts[0]), float(parts[1])))
    return keys


def full_sphere_grid(theta_step_deg: float, phi_step_deg: float) -> FarFieldGrid:
    """A whole-sphere grid in the exit frame no coarser than the run's own, with binary-exact steps.

    Theta covers 0..180 with both ends; phi covers -180 up to 180 minus one step (the circle is periodic), so the
    point counts are exact powers of two and phi = +-90 deg, the boundary between the hemispheres, lies on the grid.
    """
    n_theta = max(4, 2 ** math.ceil(math.log2(180.0 / theta_step_deg)))
    n_phi = max(4, 2 ** math.ceil(math.log2(360.0 / phi_step_deg)))
    d_theta, d_phi = 180.0 / n_theta, 360.0 / n_phi
    return FarFieldGrid(0.0, 180.0, d_theta, -180.0, 180.0 - d_phi, d_phi)


def sphere_power(df: pd.DataFrame, grid: FarFieldGrid) -> dict[str, float]:
    """Radiated power (W) of an exported far field: the integral of (|rE_theta|^2 + |rE_phi|^2) / (2 eta0) dOmega.

    Trapezoid in theta, the periodic rectangle rule in phi. Forward is the outward hemisphere x_e > 0 (|phi| < 90 deg
    in the exit frame, whose polar axis is the gap direction); points at |phi| = 90 deg count half to each side.
    """
    thetas, phis = grid.theta_values, grid.phi_values
    n = len(thetas) * len(phis)
    if len(df) != n:
        raise ValueError(f"far-field export has {len(df)} points, the grid has {n}")
    if not (np.allclose(df["Theta"].to_numpy(), np.repeat(thetas, len(phis)))
            and np.allclose(df["Phi"].to_numpy(), np.tile(phis, len(thetas)))):
        raise ValueError("far-field export points are not the grid's, theta outer and phi inner")
    intensity = sum(df[c].to_numpy() ** 2 for c in ("rEtheta_real", "rEtheta_imag", "rEphi_real", "rEphi_imag"))
    intensity = intensity.reshape(len(thetas), len(phis)) / (2 * ETA0_OHM)
    t = np.radians(thetas)
    w_theta = np.full(len(t), np.radians(grid.theta_step_deg))
    w_theta[0] *= 0.5
    w_theta[-1] *= 0.5
    w_theta *= np.sin(t)
    w_phi = np.radians(grid.phi_step_deg)
    forward_weight = np.where(np.isclose(np.abs(phis), 90.0), 0.5, (np.abs(phis) < 90.0).astype(float))
    per_phi = (w_theta[:, None] * intensity).sum(axis=0) * w_phi
    forward = float((per_phi * forward_weight).sum())
    total = float(per_phi.sum())
    return {"forward_w": forward, "backward_w": total - forward, "total_w": total}


def incident_amplitude(total: pd.DataFrame, scattered: pd.DataFrame) -> dict[str, float]:
    """|E_total - E_scattered| per point (peak V/m) on the same points; the incident wave's amplitude."""
    if len(total) != len(scattered):
        raise ValueError(f"total and scattered exports differ in length: {len(total)} and {len(scattered)}")

    def phasors(df: pd.DataFrame) -> np.ndarray:
        return np.column_stack([df[f"E{c}_real"].to_numpy() + 1j * df[f"E{c}_imag"].to_numpy() for c in "xyz"])

    magnitude = np.sqrt((np.abs(phasors(total) - phasors(scattered)) ** 2).sum(axis=1))
    return {"mean_v_per_m": float(magnitude.mean()), "min_v_per_m": float(magnitude.min()),
            "max_v_per_m": float(magnitude.max()), "points": int(len(magnitude))}


def _ratio(num: float | None, den: float | None) -> float | None:
    if num is None or den is None or den == 0 or not math.isfinite(den):
        return None
    return num / den


def derived(rec: dict[str, Any], ingoing_w: float) -> dict[str, float | None]:
    """Ratios that answer the questions; None where a quantity is missing or the geometric flux is zero."""
    exit_, ent = rec["exit"], rec["entrance"]
    e_inc = (rec["incident_e"] or {}).get("mean_v_per_m")
    cos_a = rec["cos_alpha"]
    geometric = ingoing_w * cos_a * e_inc ** 2 if e_inc is not None and cos_a > 1e-9 else None
    scattered_out = -ent["scattered"]["flux_z_w"] if ent["scattered"]["flux_z_w"] is not None else None
    ff = rec["far_field"]
    flux_out = exit_["flux_z_w"]
    return {
        "T_table": _ratio(exit_["mag_w"], ingoing_w),
        "T_flux": _ratio(flux_out, ingoing_w),
        "mag_over_flux": _ratio(exit_["mag_w"], flux_out),
        "abs_flux_over_flux": _ratio(exit_["abs_flux_z_w"], flux_out),
        "exit_over_entrance_flux": _ratio(flux_out, ent["total"]["flux_z_w"]),
        "entrance_flux_over_geometric": _ratio(ent["total"]["flux_z_w"], geometric),
        "scattered_out_over_geometric": _ratio(scattered_out, geometric),
        "far_total_over_exit_flux": _ratio((ff.get("TotalFields") or {}).get("total_w"), flux_out),
        "far_forward_over_exit_flux": _ratio((ff.get("TotalFields") or {}).get("forward_w"), flux_out),
        "far_scattered_total_over_exit_flux": _ratio((ff.get("ScatteredFields") or {}).get("total_w"), flux_out),
        "csv_relative_difference": (None if exit_["mag_w"] is None or not rec["csv_outgoing_power_w"]
                                    else exit_["mag_w"] / rec["csv_outgoing_power_w"] - 1.0),
    }


# --- AEDT side --------------------------------------------------------------------------------------------------------

def define_expression(fields: Any, name: str, ops: tuple[str, ...], surface: str) -> None:
    fields.CalcStack("clear")
    fields.CopyNamedExprToStack("Vector_RealPoynting")
    for op in ops:
        fields.CalcOp(op)
    fields.EnterSurf(surface)
    fields.CalcOp("Integrate")
    fields.AddNamedExpression(name, "Fields")


def evaluate(fields: Any, name: str, path: Path, solution: str, args: list, timeout_s: float) -> float:
    if path.exists():
        path.unlink()
    fields.CalcStack("clear")
    fields.CopyNamedExprToStack(name)
    fields.CalculatorWrite(str(path), ["Solution:=", solution], args)
    wait_for_file(path, timeout_s)
    value = read_calculator_scalar(path)
    path.unlink()
    return value


def set_field_type(hfss: Any, field_type: str) -> bool:
    hfss.odesign.GetModule("Solutions").EditSources(["FieldType:=", field_type])
    return True


def checked(result: Any, what: str) -> bool:
    """PyAEDT reports some failures by returning False rather than raising."""
    if not result:
        raise RuntimeError(f"{what} returned {result!r}")
    return True


def export_incident_line(hfss: Any, path: Path, solution: str, args: list, entrance: dict, timeout_s: float) -> pd.DataFrame:
    """E on a short line along global y through the entrance centre, INSIDE_MM inside the gap (global frame, SI)."""
    cx, cy, cz = (float(v) for v in entrance["center_mm"])
    half = 0.2 * math.sqrt(float(entrance["area_mm2"]))

    def mm(v: float) -> str:
        return f"{round(v, 9) + 0.0}mm"

    fields = hfss.odesign.GetModule("FieldsReporter")
    fields.CalcStack("clear")
    fields.EnterQty("E")
    if path.exists():
        path.unlink()
    fields.ExportOnGrid(
        str(path), [mm(cx), mm(cy - half), mm(cz + INSIDE_MM)], [mm(cx), mm(cy + half), mm(cz + INSIDE_MM)],
        ["0mm", mm(half / 5), "0mm"], solution, args,
        ["NAME:ExportOption", "IncludePtInOutput:=", True, "RefCSName:=", "Global", "PtInSI:=", True,
         "FieldInRefCS:=", False],
        "Cartesian", ["0mm", "0mm", "0mm"], False)
    wait_for_file(path, timeout_s)
    df = read_exit_field_fld(path)
    path.unlink()
    return df


def export_far_field(hfss: Any, path: Path, solution: str, freq_label: str, ephi: int, phi: float, theta: float,
                     grid: FarFieldGrid, timeout_s: float) -> dict[str, float]:
    if path.exists():
        path.unlink()
    hfss.odesign.GetModule("RadField").ExportFieldsToFile([
        "ExportFileName:=", str(path),
        "SetupName:=", SPHERE,
        "IntrinsicVariationKey:=", intrinsic_variation_key(freq_label, phi, theta),
        "DesignVariationKey:=", f"Ephi='{ephi}'",
        "SolutionName:=", solution,
        "Quantity:=", "",
    ])
    wait_for_file(path, timeout_s)
    df = read_far_field_ffd(path, grid.theta_step_deg, grid.phi_step_deg)
    path.unlink()
    return sphere_power(df, grid)


# --- the run ----------------------------------------------------------------------------------------------------------

def parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--job", type=Path, required=True, help="A finished job directory (manifest.json, project/)")
    p.add_argument("--keys", nargs="+", default=None, help="'phi,theta' pairs in degrees; default every key of the job")
    p.add_argument("--polarizations", nargs="+", type=int, default=None, help="default the job's polarizations")
    p.add_argument("--skip-far-field", action="store_true")
    p.add_argument("--skip-entrance", action="store_true")
    p.add_argument("--keep-copy", action="store_true", help="keep <job>/diagnose_poynting/project after success")
    p.add_argument("--timeout", type=float, default=None, help="per-file wait in seconds; default the job's")
    return p.parse_args(argv)


def prepare(args: argparse.Namespace) -> tuple[dict, Path, Path]:
    """Validate the job and copy its project; raise Refused before anything is changed."""
    job = args.job.resolve()
    manifest_path = job / "manifest.json"
    if not manifest_path.is_file():
        raise Refused(f"no manifest.json in {job}")
    manifest = json.loads(manifest_path.read_text())
    geo = manifest["geometry"]
    if [float(v) for v in geo["entrance_face"]["outward_normal"]] != [0.0, 0.0, -1.0] or \
            [float(v) for v in geo["exit_face"]["outward_normal"]] != [0.0, 0.0, 1.0]:
        raise Refused("this diagnostic assumes the canonical pose (entrance outward normal -z, exit +z); "
                      f"the job has {geo['entrance_face']['outward_normal']} and {geo['exit_face']['outward_normal']}")
    project = job / "project" / Path(manifest["project"]["job_copy"]).name
    results = project.with_suffix(".aedtresults")
    if not project.is_file() or not results.is_dir():
        raise Refused(f"the job's solved project is incomplete: {project} and {results} (.aedtresults) must both exist")
    work = job / WORK_DIR
    if work.exists():
        raise Refused(f"{work} exists from an earlier attempt; move or remove it first")
    copy_dir = work / "project"
    copy_dir.mkdir(parents=True)
    shutil.copy2(project, copy_dir / project.name)
    shutil.copytree(results, copy_dir / results.name)
    return manifest, work, copy_dir / project.name


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        manifest, work, copy = prepare(args)
        keys = parse_keys(args.keys) if args.keys else [
            (float(p), float(t)) for p in manifest["excitation"]["incident_phi_deg"]
            for t in manifest["excitation"]["incident_theta_deg"]]
    except (Refused, ValueError) as exc:
        print(f"diagnose_poynting: nothing started: {exc}", file=sys.stderr)
        return 2
    stale = sorted(k for k in os.environ if k.startswith("PYAEDT_"))
    for k in stale:  # as hpc/run_frequency.sbatch: a stale port or pid would attach PyAEDT to another desktop
        del os.environ[k]

    cfg = manifest["config"]
    exc_m, solver, ff_m = manifest["excitation"], manifest["solver"], manifest["far_field"]
    polarizations = args.polarizations or [int(p) for p in exc_m["polarizations"]]
    timeout = args.timeout or float(cfg["output"]["export_timeout_s"])
    freq_label = solver["frequency_label"]
    solution = f"{solver['setup']} : LastAdaptive"
    ingoing = float(exc_m["incoming_power_w"])
    ctx = SimpleNamespace(freq_label=freq_label,
                          swept_incidence=len(exc_m["incident_phi_deg"]) * len(exc_m["incident_theta_deg"]) > 1)
    grid = full_sphere_grid(float(ff_m["theta_step_deg"]), float(ff_m["phi_step_deg"]))
    scratch = work / "scratch"
    scratch.mkdir()
    entrance = manifest["geometry"]["entrance_face"]

    records: dict[tuple[int, float, float], dict] = {}
    for ephi in polarizations:
        csv = pd.read_csv(Path(manifest["outputs"][str(ephi)]) / "waveguide.csv",
                          usecols=["IWavePhi", "IWaveTheta", "OutgoingPower"]).groupby(["IWavePhi", "IWaveTheta"]).first()
        for phi, theta in keys:
            cos_a = abs(math.cos(math.radians(theta)))
            records[(ephi, phi, theta)] = {
                "ephi": ephi, "phi_deg": phi, "theta_deg": theta, "alpha_deg": 180.0 - theta, "cos_alpha": cos_a,
                "ingoing_power_w": ingoing, "csv_outgoing_power_w": float(csv.loc[(phi, theta), "OutgoingPower"]),
                "exit": {v: None for v in EXIT_FIELDS.values()},
                "entrance": {"total": {v: None for v in ENTRANCE_FIELDS.values()},
                             "scattered": {v: None for v in ENTRANCE_FIELDS.values()}},
                "incident_e": None, "far_field": {}, "derived": None,
            }
    errors: dict[str, str] = {}
    failures: list[str] = []

    def variation(ephi: int, phi: float, theta: float) -> list:
        return field_variation(ctx, ephi, phi, theta) + ["Phase:=", "0deg"]

    def attempt(label: str, primary: bool, fn):
        """Run one step; record a failure under its label and go on, so one refusal costs only its own result."""
        try:
            return fn()
        except Exception as exc:  # noqa: BLE001 - recorded in the JSON; a primary step fails the run
            errors.setdefault(label, f"{type(exc).__name__}: {exc}")
            if primary:
                failures.append(f"{label}: {exc}")
            return None

    def define_all(fields: Any, table: dict) -> list[str]:
        defined = []
        for name, (ops, surface, primary) in table.items():
            if attempt(name, primary, lambda: define_expression(fields, name, ops, surface) or True):
                defined.append(name)
        return defined

    def evaluate_all(fields: Any, names: list[str], table: dict, store) -> None:
        for name in names:
            for key, rec in records.items():
                value = attempt(name, table[name][2],
                                lambda: evaluate(fields, name, scratch / f"{name}.fld", solution, variation(*key), timeout))
                store(rec)[EXIT_FIELDS.get(name) or ENTRANCE_FIELDS[name]] = value

    from bbsim import session  # lazy: PyAEDT is imported only inside the session

    totals_e: dict[tuple[int, float, float], pd.DataFrame] = {}
    try:
        with session.aedt_session(copy, manifest["project"]["run_design"], cfg["project"]["aedt_version"]) as hfss:
            if manifest["project"]["run_design"] not in hfss.design_list:
                raise RuntimeError(f"design {manifest['project']['run_design']!r} not in the project copy")
            fields = hfss.odesign.GetModule("FieldsReporter")
            # 1. Everything that changes nothing in the modeler, with total fields.
            attempt("field type TotalFields", True, lambda: set_field_type(hfss, "TotalFields"))
            evaluate_all(fields, define_all(fields, EXIT_EXPRESSIONS), EXIT_EXPRESSIONS, lambda rec: rec["exit"])
            for key in records:
                df = attempt("incident line TotalFields", True, lambda: export_incident_line(
                    hfss, scratch / "e_total.fld", solution, variation(*key), entrance, timeout))
                if df is not None:
                    totals_e[key] = df
            sphere = None if args.skip_far_field else attempt("far-field sphere", True, lambda: checked(
                hfss.insert_infinite_sphere(**far_field_sphere_kwargs(
                    inspect.signature(hfss.insert_infinite_sphere).parameters, grid, SPHERE)), "insert_infinite_sphere"))
            if sphere is not None:
                for key, rec in records.items():
                    rec["far_field"]["TotalFields"] = attempt("far_field TotalFields", True, lambda: export_far_field(
                        hfss, scratch / "far.ffd", solution, freq_label, *key, grid, timeout))
            # 2. Scattered fields: the incident wave as total minus scattered, and the scattered far field.
            if attempt("field type ScatteredFields", True, lambda: set_field_type(hfss, "ScatteredFields")):
                for key, rec in records.items():
                    scattered = attempt("incident line ScatteredFields", True, lambda: export_incident_line(
                        hfss, scratch / "e_scattered.fld", solution, variation(*key), entrance, timeout))
                    if scattered is not None and key in totals_e:
                        rec["incident_e"] = incident_amplitude(totals_e[key], scattered)
                    if sphere is not None:
                        rec["far_field"]["ScatteredFields"] = attempt("far_field ScatteredFields", False, lambda: export_far_field(
                            hfss, scratch / "far.ffd", solution, freq_label, *key, grid, timeout))
            # 3. The entrance face list is the one change to the modeler, so it comes last.
            if not args.skip_entrance and attempt("entrance face list", True, lambda: checked(hfss.modeler.create_face_list(
                    [int(entrance["id"])], name=ENTRANCE_FACE_LIST), "create_face_list")):
                defined = define_all(fields, ENTRANCE_EXPRESSIONS)
                for field_type, label in (("ScatteredFields", "scattered"), ("TotalFields", "total")):
                    if attempt(f"field type {field_type}", True, lambda: set_field_type(hfss, field_type)):
                        evaluate_all(fields, defined, ENTRANCE_EXPRESSIONS, lambda rec: rec["entrance"][label])
            attempt("field type TotalFields", True, lambda: set_field_type(hfss, "TotalFields"))
    except Exception as exc:  # noqa: BLE001 - the JSON records how far the run got
        failures.append(f"{type(exc).__name__}: {exc}")

    mismatches = []
    for (ephi, phi, theta), rec in records.items():
        rec["derived"] = derived(rec, ingoing)
        diff = rec["derived"]["csv_relative_difference"]
        if diff is not None and abs(diff) > CSV_TOLERANCE:
            mismatches.append(f"Ephi={ephi} phi={phi} theta={theta}: exit |Re S| integral {rec['exit']['mag_w']:.6e} W "
                              f"against the recorded OutgoingPower {rec['csv_outgoing_power_w']:.6e} W ({diff:+.2e})")
    status = "failed" if failures else ("csv mismatch" if mismatches else "ok")
    write_json(work / "poynting.json", {
        "status": status, "failures": failures, "csv_mismatches": mismatches, "expression_errors": errors,
        "job": str(args.job.resolve()), "manifest_sha256": sha256_file(args.job.resolve() / "manifest.json"),
        "code_commit": git_commit(), "solution": solution, "ingoing_power_w": ingoing,
        "far_field_grid": {"theta_step_deg": grid.theta_step_deg, "phi_step_deg": grid.phi_step_deg,
                           "points": grid.points_per_angle, "coordinate_system": EXIT_CS},
        "keys": list(records.values()),
    })
    report(records, status, failures, mismatches, work)
    if status == "ok":
        shutil.rmtree(scratch, ignore_errors=True)
        if not args.keep_copy:
            shutil.rmtree(work / "project")
        return 0
    return 1


def _fmt(v: float | None, spec: str = ".4f") -> str:
    return "-" if v is None else format(v, spec)


def report(records: dict, status: str, failures: list[str], mismatches: list[str], work: Path) -> None:
    print(f"== diagnose_poynting: {status}; record {work / 'poynting.json'}")
    for line in failures + mismatches:
        print(f"   {line}")
    print("   Ephi  phi  theta  T_table  T_flux  |S|/S_z  out/in  in/geom  scat/geom  |E_inc|  ff/flux  ff_fwd/flux  ff_scat/flux")
    for (ephi, phi, theta), rec in records.items():
        d, e = rec["derived"], rec["incident_e"] or {}
        print(f"   {ephi:>4} {phi:>4.0f} {theta:>6.0f}  {_fmt(d['T_table'])}  {_fmt(d['T_flux'])}  "
              f"{_fmt(d['mag_over_flux'])}  {_fmt(d['exit_over_entrance_flux'])}  {_fmt(d['entrance_flux_over_geometric'])}  "
              f"{_fmt(d['scattered_out_over_geometric'])}  {_fmt(e.get('mean_v_per_m'))}  "
              f"{_fmt(d['far_total_over_exit_flux'])}  {_fmt(d['far_forward_over_exit_flux'])}  "
              f"{_fmt(d['far_scattered_total_over_exit_flux'])}")
    parts = []
    for (ephi, phi, theta), rec in records.items():
        d = rec["derived"]
        if d["T_table"] is not None and d["T_table"] > 1e-3:
            parts.append(f"Ephi={ephi} ({phi:g},{theta:g}) T={_fmt(d['T_table'])} |S|/S_z={_fmt(d['mag_over_flux'])} "
                         f"out/in={_fmt(d['exit_over_entrance_flux'])} in/geom={_fmt(d['entrance_flux_over_geometric'])} "
                         f"ff/flux={_fmt(d['far_total_over_exit_flux'])}")
    print(f"POYNTING: {status}; " + "; ".join(parts))


if __name__ == "__main__":
    sys.exit(main())
