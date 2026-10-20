#!/usr/bin/env python3
"""OR-035 check on copies of the solved projects, with no solve (HPC step 8d; HF-045, CDX-023, OR-035).

Two checks of what OR-035 changed, each on real AEDT output and each printed with its verdict:

(a) HF-045: a rectangle whose extent is not a whole number of steps gets the fitted step, so both walls are lattice
    nodes. On a copy of crack1's solved project, the runner's own lattice (``extract.exit_grid``, which calls
    ``fit_axis``) and export call (``extract.extract_waveguide``) run over the exit face's own y lattice and
    z in [-0.01875, 0.01875] mm at the configured z resolution: 37.5 steps of 1 um, so 38 intervals of 0.0375/38 mm
    and 39 nodes. Every key's export must hold 101 x 39 rows, with the z end nodes at -1.875e-5 and 1.875e-5 m to
    1e-12 m (``extract.NODE_TOL_M``) and every point on x_e = 0. Ephi=0 only; no power is evaluated, so the tables'
    OutgoingPower column is a placeholder 0 and is not read.
(b) CDX-023: the exit-face outline (its vertices, or samples along its edges for a face with fewer than 3) of each
    of the three projects, in the exit coordinate system, through the runner's own ``extract.describe_exit_face``,
    must lie on x_e = 0 to 1e-12 m (``extract.check_exit_plane_bounds``): max |x| <= 1e-9 mm.

Run on HPC in a batch job (hpc/check_or035.sbatch), one AEDT session per project and no solver licence:

    python hpc/check_or035.py [--crack1 JOB] [--crack2 JOB] [--round-gap JOB]

The defaults are the three solved jobs: crack1 and crack2 ``job_tight1`` (step 5d) and the round gap
``job_reference1`` (step 8). Each job's project is copied into ``<job>/check_or035/project/``, as
hpc/diagnose_poynting.py copies a solved job: crack1's with its ``.aedtresults``, which the export needs, the other
two as the ``.aedt`` alone, which their geometry needs. AEDT opens only the copies and nothing is saved; the source
jobs are never written. The copies are removed after a passing run and kept after a failure, with the raw export of
a failed key in ``check_or035/scratch/``. The last line is ``== RESULT: PASS`` (exit 0) or ``== RESULT: FAIL <what
failed>`` (exit 1). A usage error, a job that cannot be read, or an existing ``check_or035/`` directory ends
``== RESULT: FAIL nothing started: <reason>`` (exit 2) before anything is copied.
"""
from __future__ import annotations

import argparse
import dataclasses
import logging
import math
import shutil
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bbsim.config import ExitFieldConfig  # noqa: E402
from bbsim.dataset_v2 import key_text  # noqa: E402
from bbsim.extract import (NODE_TOL_M, PLANE_TOL_M, ExitFaceGeometry, ExitGrid, check_exit_plane_bounds,  # noqa: E402
                           exit_grid, fit_axis)
from bbsim.geometry import GeometryError  # noqa: E402

WORK_DIR = "check_or035"
Z_HALF_MM = 0.01875       # HF-045's test extent: 37.5 um, 37.5 steps of the configured 1 um, never a whole number
PLANE_TOL_MM = PLANE_TOL_M * 1e3
OUTPUTS = Path("/home/rshenoy/BBRSim/outputs")
DEFAULTS = {"crack1": OUTPUTS / "InfParallelPlate_crack1Rohan_500GHz" / "job_tight1",
            "crack2": OUTPUTS / "InfParallelPlate_crack2_500GHz" / "job_tight1",
            "round gap": OUTPUTS / "RoundGap_r50um_2000GHz" / "job_reference1"}
KEYS = ["IWavePhi", "IWaveTheta"]


class Refused(RuntimeError):
    """Nothing was started: a job that cannot be read, or a copy from an earlier attempt."""


# --- pure parts -------------------------------------------------------------------------------------------------------

def hf045_grid(face: ExitFaceGeometry, exit_field: ExitFieldConfig, z_half_mm: float = Z_HALF_MM
               ) -> tuple[ExitFieldConfig, ExitGrid]:
    """Check (a)'s exit-field settings and lattice: the face's own x and y bounds, z in [-z_half_mm, z_half_mm], at
    the configured resolution, through the runner's ``exit_grid``, which fits the z step (``fit_axis``). A z extent
    that is a whole number of steps would not test HF-045 and is refused (ValueError)."""
    own = exit_grid(face, exit_field)
    boundary = ((own.lo_mm[0], own.lo_mm[1], -z_half_mm), (own.hi_mm[0], own.hi_mm[1], z_half_mm))
    settings = dataclasses.replace(exit_field, boundary_mm=boundary)
    grid = exit_grid(face, settings)
    if not grid.fitted[2]:
        raise ValueError(f"z in [{-z_half_mm}, {z_half_mm}] mm is a whole number of {exit_field.resolution_mm[2]} mm "
                         f"steps on a {face.cross_section['shape']} exit face; HF-045's check needs a fitted axis")
    return settings, grid


def hf045_expectation(grid: ExitGrid) -> dict[str, Any]:
    """What every key's export must hold, from the lattice alone: the node counts along y and z, the rows per key,
    the declared bounds in metres (``ExitGrid.axis_record``, as the sidecar declares them), the steps and the
    ExportOnGrid spacing strings."""
    axes = grid.axis_record()
    counts = tuple(grid.counts[1:])
    return {"counts": counts, "rows_per_key": math.prod(counts),
            "y_m": (axes["y_e"]["min"], axes["y_e"]["max"]), "z_m": (axes["z_e"]["min"], axes["z_e"]["max"]),
            "step_mm": tuple(grid.step_mm), "step_strings": grid.step_strings()}


def _nodes(values: pd.Series, lo_m: float, step_m: float) -> int:
    """How many distinct lattice nodes the values sit on, each value taken to its nearest node index."""
    return int(len(np.unique(np.rint((values.to_numpy(dtype=float) - lo_m) / step_m))))


def hf045_problems(table: pd.DataFrame, grid: ExitGrid) -> list[str]:
    """What is wrong with an export of ``grid`` (exit frame, metres, one block per key): a key whose block does not
    hold the lattice's rows, y nodes or z nodes; a z end node off its declared wall by more than ``NODE_TOL_M``; a
    point off the plane x_e = 0 by more than ``PLANE_TOL_M``. Empty when the export is the lattice."""
    if table.empty:
        return ["the export holds no rows"]
    if not np.all(np.isfinite(table[KEYS + ["X", "Y", "Z"]].to_numpy(dtype=float))):
        return ["the export holds non-finite keys or coordinates"]
    expected = hf045_expectation(grid)
    ny, nz = expected["counts"]
    step_y, step_z = grid.step_mm[1] * 1e-3, grid.step_mm[2] * 1e-3
    problems = []
    for (phi, theta), block in table.groupby(KEYS, sort=False):
        where = key_text((phi, theta))
        if len(block) != expected["rows_per_key"]:
            problems.append(f"{where}: {len(block)} rows, expected {ny} x {nz} = {expected['rows_per_key']}")
        found = (_nodes(block.Y, expected["y_m"][0], step_y), _nodes(block.Z, expected["z_m"][0], step_z))
        for name, n, want in (("y", found[0], ny), ("z", found[1], nz)):
            if n != want:
                problems.append(f"{where}: {n} {name} nodes, expected {want}")
        for side, value, wall in (("low", float(block.Z.min()), expected["z_m"][0]),
                                  ("high", float(block.Z.max()), expected["z_m"][1])):
            if not abs(value - wall) <= NODE_TOL_M:
                problems.append(f"{where}: z end node {value!r} m on the {side} wall, expected {wall!r} m to "
                                f"{NODE_TOL_M:g} m")
        x = float(np.max(np.abs(block.X.to_numpy(dtype=float))))
        if not x <= PLANE_TOL_M:
            problems.append(f"{where}: points off the exit plane x_e = 0, max |X| = {x!r} m (limit {PLANE_TOL_M:g} m)")
    return problems


def plane_problems(local_points_mm: np.ndarray, label: str) -> tuple[float, list[str]]:
    """Max |x| (mm) of the exit-face outline in the exit frame, and the runner's refusal if it is off the exit plane
    (``check_exit_plane_bounds``: x_e = 0 to 1e-12 m on both bounds), named by ``label``."""
    points = np.asarray(local_points_mm, dtype=float).reshape(-1, 3)
    if not points.size or not np.all(np.isfinite(points)):
        return math.nan, [f"{label}: the exit-face outline is empty or not finite"]
    max_x = float(np.max(np.abs(points[:, 0])))
    try:
        check_exit_plane_bounds(points.min(axis=0), points.max(axis=0))
    except GeometryError:
        return max_x, [f"{label}: the exit-face outline is off the exit plane, max |x| = {max_x:.3g} mm (limit "
                       f"{PLANE_TOL_MM:g} mm)"]
    return max_x, []


# --- the run ---------------------------------------------------------------------------------------------------------

@dataclasses.dataclass(frozen=True)
class Job:
    label: str
    source: Any          # bbsim.run_frequency.SourceJob
    work: Path
    export: bool         # check (a) runs on this job's copy, which then needs the .aedtresults too


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--crack1", type=Path, default=DEFAULTS["crack1"], help="crack1's solved job; checks (a) and (b)")
    p.add_argument("--crack2", type=Path, default=DEFAULTS["crack2"], help="crack2's solved job; check (b)")
    p.add_argument("--round-gap", type=Path, default=DEFAULTS["round gap"],
                   help="the round gap's solved job; check (b)")
    return p.parse_args(argv)


def prepare(paths: dict[str, Path]) -> list[Job]:
    """Read every job as the runner reads a --from-solved source, before anything is copied; refuse a job that
    cannot be read, a work directory left by an earlier attempt, and a crack1 resolution that would not fit."""
    from bbsim.run_frequency import load_source_job

    jobs = []
    for label, path in paths.items():
        try:
            source = load_source_job(path)
        except Exception as exc:  # noqa: BLE001 - any unreadable job stops the run before anything starts
            raise Refused(f"the {label} job {path} cannot be read: {exc}") from exc
        work = source.job_dir / WORK_DIR
        if work.exists():
            raise Refused(f"{work} exists from an earlier attempt; move or remove it first")
        export = label == "crack1"
        resolution = source.config.exit_field.resolution_mm[2]
        if export and not fit_axis(-Z_HALF_MM, Z_HALF_MM, resolution)[1]:
            raise Refused(f"crack1's configured z resolution {resolution} mm makes [{-Z_HALF_MM}, {Z_HALF_MM}] mm a "
                          "whole number of steps; HF-045's check needs a fitted axis")
        jobs.append(Job(label, source, work, export))
    return jobs


def copy_project(job: Job) -> Path:
    """The job's project copied under its work directory (with its .aedtresults when the export needs them)."""
    copy_dir = job.work / "project"
    copy_dir.mkdir(parents=True)
    (job.work / "scratch").mkdir()
    project = copy_dir / job.source.project.name
    shutil.copy2(job.source.project, project)
    if job.export:
        shutil.copytree(job.source.results, copy_dir / job.source.results.name)
    return project


def inspect(job: Job, project: Path) -> dict[str, list[str]]:
    """Open the copy, describe its exit face as the runner does, judge (b), and on crack1 export and judge (a).
    Returns the problems per check, keyed "a" and "b"."""
    from bbsim import session
    from bbsim.extract import ExtractionContext, describe_exit_face
    from bbsim.geometry import select_faces
    from bbsim.hfss_setup import EXIT_CS, far_field_sphere
    from bbsim.naming import last_adaptive_solution
    from bbsim.sampling import incident_angle_values, incoming_power_w

    cfg = job.source.config
    run_design = job.source.manifest["project"]["run_design"]
    freq = cfg.solver.frequency_ghz
    exc = cfg.excitation
    found: dict[str, list[str]] = {}
    with session.aedt_session(project, run_design, cfg.project.aedt_version) as hfss:
        if run_design not in hfss.design_list:
            raise RuntimeError(f"design {run_design!r} not in the copy {project}; designs: {list(hfss.design_list)}")
        if hfss.design_name != run_design:
            hfss.set_active_design(run_design)
        entrance, exit_face = select_faces(hfss, cfg.geometry)
        sphere, ff_grid = far_field_sphere(cfg.far_field, freq)
        ctx = ExtractionContext(
            frequency_ghz=freq, solution=last_adaptive_solution(freq), sphere_name=sphere,
            phi_values=incident_angle_values(exc.phi_lower_deg, exc.phi_upper_deg, exc.phi_step_deg),
            theta_values=incident_angle_values(exc.theta_lower_deg, exc.theta_upper_deg, exc.theta_step_deg),
            grid=ff_grid, exit_field=cfg.exit_field, exit_face=exit_face, exit_cs_name=EXIT_CS,
            exit_cs_x=cfg.geometry.exit_cs_x, exit_cs_y=cfg.geometry.exit_cs_y, scratch_dir=job.work / "scratch",
            timeout_s=cfg.output.export_timeout_s, incoming_power_w=incoming_power_w(exc.ei_v_per_m, entrance.area_mm2),
        )
        face = describe_exit_face(hfss, ctx)
        max_x, found["b"] = plane_problems(face.local_points_mm, job.label)
        outline = (f"{face.vertex_count} vertices" if face.bounds_method == "vertices"
                   else f"{len(face.local_points_mm)} edge samples ({face.vertex_count} vertices)")
        print(f"   (b) {job.label}: exit face {exit_face.id}, {outline}, {face.cross_section['shape']}: max |x| = "
              f"{max_x:.3g} mm (limit {PLANE_TOL_MM:g} mm): {'FAIL' if found['b'] else 'PASS'}")
        for problem in found["b"]:
            print(f"       {problem}")
        if job.export:
            found["a"] = export_check(hfss, ctx, face)
    return found


def export_check(hfss: Any, ctx: Any, face: ExitFaceGeometry) -> list[str]:
    """Check (a) on the open copy of crack1: the runner's export over HF-045's fitted lattice, judged per key."""
    from bbsim.extract import extract_waveguide
    from bbsim.hfss_setup import set_total_fields

    settings, grid = hf045_grid(face, ctx.exit_field)
    expected = hf045_expectation(grid)
    ny, nz = expected["counts"]
    zlo, zhi = expected["z_m"]
    print(f"   (a) HF-045 on crack1: z in [{-Z_HALF_MM}, {Z_HALF_MM}] mm at the configured "
          f"{grid.resolution_mm[2]} mm: fitted step {grid.step_mm[2]!r} mm (fit_axis); ExportOnGrid spacing "
          f"{expected['step_strings']}")
    print(f"       expect {ny} x {nz} = {expected['rows_per_key']} rows per key, z end nodes at {zlo!r} and {zhi!r} m "
          f"to {NODE_TOL_M:g} m")
    set_total_fields(hfss)   # as prepare_from_solved: the field type is a post-processing setting of the copy
    keys = [(float(p), float(t)) for p in ctx.phi_values for t in ctx.theta_values]
    try:
        table = extract_waveguide(hfss, dataclasses.replace(ctx, exit_field=settings), 0, face,
                                  exit_powers=({k: 0.0 for k in keys}, {}))
    except Exception as exc:  # noqa: BLE001 - the runner's refusal (a short export, a missing wall node) is the verdict
        problems = [f"the runner's export refused it: {type(exc).__name__}: {exc}"]
    else:
        sizes = sorted(set(table.groupby(KEYS, sort=False).size()))
        print(f"       {table.groupby(KEYS).ngroups} keys exported (Ephi=0): rows per key {sizes}; z from "
              f"{float(table.Z.min())!r} to {float(table.Z.max())!r} m; max |X| "
              f"{float(np.max(np.abs(table.X.to_numpy(dtype=float))))!r} m")
        problems = hf045_problems(table, grid)
    for problem in problems:
        print(f"       {problem}")
    print(f"   HF-045: {'FAIL' if problems else 'PASS'}")
    return problems


def main(argv: list[str] | None = None) -> int:
    try:
        args = parse_args(argv)
    except SystemExit as exc:                 # argparse has printed its message (or the help)
        if not exc.code:
            return 0
        print("== RESULT: FAIL nothing started: usage (see above)")
        return 2
    paths = {"crack1": args.crack1, "crack2": args.crack2, "round gap": args.round_gap}
    try:
        jobs = prepare(paths)
    except Refused as exc:
        print(f"== RESULT: FAIL nothing started: {exc}")
        return 2
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    print("== check_or035: " + "; ".join(f"{job.label} {job.source.job_dir}" for job in jobs))
    failed: list[str] = []
    for job in jobs:
        print(f"== {job.label}: a copy of {job.source.project}" + (" and its .aedtresults" if job.export else ""))
        try:
            found = inspect(job, copy_project(job))
        except Exception as exc:  # noqa: BLE001 - one project's failure is its checks' failure; the others still run
            print(f"   {job.label}: stopped: {type(exc).__name__}: {exc}")
            failed.append(f"{job.label} stopped ({type(exc).__name__})")
            continue
        if found.get("a"):
            failed.append("HF-045 (a) on crack1")
        if found.get("b"):
            failed.append(f"CDX-023 (b) {job.label} exit face off the plane")
    if failed:
        kept = [str(job.work) for job in jobs if job.work.exists()]
        print(f"   copies kept for inspection (remove them by hand before a rerun): {', '.join(kept) or 'none'}")
        print(f"== RESULT: FAIL {'; '.join(failed)}")
        return 1
    for job in jobs:
        try:
            shutil.rmtree(job.work)
        except OSError as exc:  # AEDT or NFS may still hold files; the verdict is already known
            print(f"   note: could not remove {job.work} ({exc}); remove it by hand")
    print("== RESULT: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
