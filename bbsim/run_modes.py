"""CLI: the ports-only mode run of a hollow round fibre bore with a lossy silver wall (ledger OR-053).

One job builds the bore and solves its ports, in its own directory ``<root>/<dataset_id>_modes/job_<id>/``:

1. The config (``bbsim.modes_config``) is loaded and printed; ``--dry-run`` stops there.
2. The licence guard is checked (OR-045, HF-057), the code checkout read, and the job directory made; an existing one
   is refused and left as it is.
3. AEDT session 1, ``build_project.build(..., reopen=False)``, builds the bore and saves ``project/<dataset_id>.aedt``.
   Its launch lease is released by its own verified cleanup before session 2 starts, since a one-slot reservation
   refuses a second lease while one is held (HF-058(2)). It sets no temp directory (mode off). The build verifies
   nothing the run relies on: with the reopen skipped its outcome only looks verified.
4. AEDT session 2 (``run_frequency.open_session``, under the CLI's temp mode) opens that project, sets driven modal,
   and is the verification: the opened design is the spec's bore (``build_project.verify_build``: name, units, one
   bare vacuum object, the bounding box and both end faces to 1e-6 mm). port1 goes on the end face at z = 0, port2 on
   the one at z = L, the wall on the remaining face, curved. Its area as AEDT gives it is recorded beside 2 pi r L, as
   information: the faces are identified without it.
5. The wall, the mesh, two N-mode wave ports and one ports-only setup per frequency (``bbsim.modes_setup``), then one
   solve of every setup with the config's cores; AEDT's messages are logged and kept in ``logs/aedt_messages.txt``.
   ``analyze`` returning True proves nothing, so every setup must then hold its one solution before the project is
   saved.
6. Per setup, every port mode's Gamma, and a Touchstone export into ``scratch/`` as evidence; a failed export is a
   warning row, not a failure, since the acceptance check reads ``modes.csv`` alone. AEDT's messages are kept again
   after the reads, so an export's reason stays in the record.
7. AEDT's temp directory is restored, ``modes.csv`` is written and self-checked (``bbsim.modes``), and
   ``modes_manifest.json`` records the run. Its top-level ``acceptance`` block, {bore_um: the bore's diameter in um,
   rho_ohm_m: 1/sigma, null for a pec wall}, is what the bridge script hands to the acceptance check.

Exit codes: 0, every self-check row of fail severity passed; 1, a self-check failed (the CSV is written), or a read, a
write or the self-check itself failed after a finished solve (the solved project, AEDT's messages and a manifest whose
``passed`` is false are kept); 2, refused before AEDT (arguments, config, the guard, the checkout, an existing job
directory); 3, the build or an AEDT session failed, the design would not take driven modal, the opened design is not
the spec's bore, the solve left a setup without its solution, or the project could not be saved after the solve, since
the solution may then not be on disk (the project and a manifest are kept).

Once the run log is open and ``effective_config.json`` written, every run writes the manifest however it ends,
SIGTERM aside, with a null sha256 for a project it cannot read; a failure to open the log or write that file stops
the run before, with a traceback on stderr and exit 1. SIGTERM raises ``SystemExit(143)`` as in ``run_frequency``, and
its cleanup runs on the way out. The last line on stdout is ``== modes: PASS`` or ``== modes: FAIL <what>``.
"""
from __future__ import annotations

import argparse
import json
import logging
import math
import sys
import time
from collections.abc import Mapping, Sequence
from contextlib import ExitStack
from dataclasses import asdict
from pathlib import Path
from typing import Any

from bbsim import aedt_temp as tempdir
from bbsim import licence_guard
from bbsim.build_project import TOL_MM, build, verify_build
from bbsim.config import FaceSelector
from bbsim.geometry import GeometryError, select_face_on_bounding_plane
from bbsim.geometry_spec import GeometrySpec
from bbsim.hfss_setup import aedt_messages, solve
from bbsim.inventory import inventory_design
from bbsim.manifest import (REPO_ROOT, EnvironmentState, EnvironmentStateError, build_manifest, environment_state,
                            sha256_file, write_json)
from bbsim.modes import CheckRow, GammaRow, self_check, write_modes_csv
from bbsim.modes_config import ModesConfig, load_modes_config
from bbsim.modes_setup import (assign_mesh, assign_ports, assign_wall, export_touchstone_with_gamma,
                               insert_ports_only_setup, read_gamma, set_driven_modal, solution_of)
from bbsim.paths import JobDirs, create_modes_job_dirs, default_job_id
from bbsim.run_frequency import (_environment_log_line, _guard_log_line, _handle_sigterm, _log_system_exit,
                                 _require_vacuum, _setup_logging, allocation_warning, check_guard_before_job,
                                 open_session)
from bbsim.session import DEFAULT_LAUNCH_TIMEOUT_S, aedt_versions

log = logging.getLogger("bbsim.run_modes")

MANIFEST_NAME = "modes_manifest.json"
CSV_NAME = "modes.csv"
MESSAGES_NAME = "aedt_messages.txt"
TEMP_MODES = ("off", "task")   # env is not offered: AEDT 2025.2 ignored TMPDIR on the cluster (OR-047 amended)
BUILD_TEMP = tempdir.TempSettings()   # the build's launch sets no temp directory
_NAMED = 8


class _Stopped(Exception):
    """The run stopped at ``stage`` with exit ``code`` (1 after a finished solve, 3 before it)."""

    def __init__(self, code: int, stage: str, error: str) -> None:
        super().__init__(error)
        self.code, self.stage, self.error = code, stage, error


def _text(exc: BaseException) -> str:
    return f"{type(exc).__name__}: {exc}"


def parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="The fibre bore's ports-only mode run (OR-053): build the bore in its own "
                                            "job, solve its two wave ports at every frequency, write modes.csv.")
    p.add_argument("--config", type=Path, required=True, help="modes TOML (configs/fibre_modes_*.toml)")
    p.add_argument("--job-id", default=None,
                   help="Defaults to slurm<SLURM_JOB_ID> or a UTC timestamp; an existing job directory is refused")
    p.add_argument("--aedt-temp", default="off", choices=TEMP_MODES,
                   help="AEDT's temp directory for the run's session (HF-057): off (the default) leaves it; task "
                        "points it at <job>/aedt_temp and sets it to --aedt-temp-restore at the end, for a job alone "
                        "on its node. The build's session sets nothing")
    p.add_argument("--aedt-temp-restore", default=None, metavar="PATH",
                   help="task mode: the absolute TempDirectory AEDT is set back to at the end")
    p.add_argument("--launch-timeout-s", type=float, default=DEFAULT_LAUNCH_TIMEOUT_S,
                   help="seconds PyAEDT waits for the run's AEDT to start (its settings.desktop_launch_timeout); the "
                        "build's session keeps the default")
    p.add_argument("--dry-run", action="store_true", help="load and print the config; no job directory, no AEDT")
    args = p.parse_args(argv)
    if not (math.isfinite(args.launch_timeout_s) and args.launch_timeout_s > 0):
        p.error(f"--launch-timeout-s must be a finite number of seconds above 0, got {args.launch_timeout_s!r}")
    try:
        mode, template = tempdir.parse_mode(args.aedt_temp)
        args.temp_settings = tempdir.TempSettings(mode, restore=args.aedt_temp_restore, template=template)
    except ValueError as exc:
        p.error(str(exc))
    return args


# --- what the run hands on ---------------------------------------------------------------------------------------------

def bore_diameter_um(cfg: ModesConfig) -> float:
    return 2 * cfg.radius_m * 1e6


def acceptance_block(cfg: ModesConfig) -> dict[str, float | None]:
    """The manifest's ``acceptance`` block, which the bridge script reads: the bore's diameter in um and the wall's
    resistivity 1/sigma in ohm m (None for a pec wall)."""
    rho = None if cfg.wall.kind == "pec" else 1.0 / float(cfg.wall.conductivity_s_per_m)
    return {"bore_um": bore_diameter_um(cfg), "rho_ohm_m": rho}


def csv_provenance(cfg: ModesConfig, environment: EnvironmentState, versions: Mapping[str, str]) -> dict[str, str]:
    """``modes.csv``'s ``#`` lines: the config and spec, the checkout as read before the job, and the run session's
    AEDT and PyAEDT versions. An unset conductivity or thickness is empty."""
    wall = cfg.wall

    def number(value: float | None) -> str:
        return "" if value is None else repr(float(value))

    return {"dataset_id": cfg.dataset_id, "bore_um": repr(bore_diameter_um(cfg)), "wall_kind": wall.kind,
            "conductivity_s_per_m": number(wall.conductivity_s_per_m), "thickness_um": number(wall.thickness_um),
            "modes_per_port": str(cfg.ports.modes), "bbsim_commit": environment.commit,
            "tree_clean": "false" if environment.dirty else "true",
            "aedt_version": str(versions.get("aedt", "")), "pyaedt_version": str(versions.get("pyaedt", ""))}


def config_record(cfg: ModesConfig) -> dict[str, Any]:
    return {"geometry_spec": str(cfg.geometry_spec), "wall": asdict(cfg.wall), "ports": asdict(cfg.ports),
            "solver": asdict(cfg.solver), "output": {"root": str(cfg.output.root)}}


def describe(cfg: ModesConfig, config_path: Path) -> list[str]:
    """The config as the run reads it, printed before anything starts and repeated in the run log."""
    s, wall, ports, solver = cfg.spec, cfg.wall, cfg.ports, cfg.solver
    if wall.kind == "pec":
        wall_text = "pec (the default boundary: nothing is assigned)"
    else:
        sigma = float(wall.conductivity_s_per_m)
        wall_text = (f"finite_conductivity, {sigma:g} S/m (rho {1 / sigma:.6g} ohm m), "
                     + ("bulk (no thickness)" if wall.thickness_um is None else f"thickness {wall.thickness_um:g} um"))
    length = solver.port_face_max_length_mm
    return [
        f"Fibre modes run (OR-053): config {config_path}",
        f"  geometry spec {cfg.geometry_spec}: design {s.project.design}, object {s.project.object}, dataset_id "
        f"{cfg.dataset_id}, AEDT {s.project.aedt_version}",
        f"  bore: {s.solid.material} cylinder, diameter {bore_diameter_um(cfg):g} um (radius {s.solid.radius_mm:g} mm), "
        f"length {s.solid.length_mm:g} mm, axis +z from z = 0",
        f"  wall: {wall_text}",
        f"  ports: port1 on the end face at z = 0, port2 on the one at z = {s.solid.length_mm:g} mm, {ports.modes} "
        f"modes each",
        f"  setups: {len(ports.frequencies_ghz)} ports-only, port accuracy {solver.port_accuracy_percent} %, at "
        f"{', '.join(f'{f:g}' for f in ports.frequencies_ghz)} GHz",
        f"  solver: {solver.cores} cores, curvilinear {solver.curvilinear}, wall facets "
        f"{solver.wall_normal_deviation_deg:g} deg, port-face length "
        + ("none" if length is None else f"{length:g} mm"),
        f"  output: {cfg.output.root}/{cfg.dataset_id}_modes/job_<id>",
    ]


# --- the opened design against the spec ------------------------------------------------------------------------------

class BoreMismatch(RuntimeError):
    """The opened design is not the spec's bore; ``differences`` names each way it is not."""

    def __init__(self, differences: Sequence[str]) -> None:
        self.differences = list(differences)
        super().__init__("the opened design is not the spec's bore, so nothing is set up: " + "; ".join(differences))


def bore_faces(design: Mapping[str, Any], spec: GeometrySpec) -> dict[str, Any]:
    """The run session's verification of the opened design (the build's is not one), and its faces.

    The design must be the spec's bare bore (``verify_build``: name, units, one object, material, the bounding box and
    both planar end faces on the axis to 1e-6 mm, their areas) and the object vacuum (CDX-029). The entrance, port1, is
    the planar face on the spec's plane z = 0; the exit, port2, the one on z = L; the side, which takes the wall, the
    one face left, curved. Returns {entrance, exit: FaceInfo, side: its inventory record, side_area: {area_mm2 as AEDT
    gave it, expected_mm2 = 2 pi r L, ratio}}, the area information only (null when AEDT gave none);
    ``BoreMismatch`` otherwise, and nothing is set up."""
    diffs = verify_build(design, spec)
    if diffs:
        raise BoreMismatch(diffs)
    (obj,) = design["objects"]
    try:
        _require_vacuum(obj)
        box = spec.solid.expected_bounding_box_mm
        entrance = select_face_on_bounding_plane(obj["faces"], box, FaceSelector("z", "min"), tol_mm=TOL_MM)
        exit_face = select_face_on_bounding_plane(obj["faces"], box, FaceSelector("z", "max"), tol_mm=TOL_MM)
    except (RuntimeError, GeometryError) as exc:
        raise BoreMismatch([str(exc)]) from exc
    problems = []
    for what, face, z in (("port1's entrance", entrance, 0.0), ("port2's exit", exit_face, float(spec.solid.length_mm))):
        if abs(float(face.center_mm[2]) - z) > TOL_MM:
            problems.append(f"{what} face {face.id} is at z = {face.center_mm[2]} mm, not {z} mm")
    rest = [f for f in obj["faces"] if f["id"] not in (entrance.id, exit_face.id)]
    side = rest[0] if len(rest) == 1 and rest[0].get("is_planar") is not True else None
    if side is None:
        problems.append(f"expected one curved side face besides the end faces {entrance.id} and {exit_face.id}, "
                        f"found {[(f['id'], f.get('is_planar')) for f in rest]}")
    if problems:
        raise BoreMismatch(problems)
    expected = 2 * math.pi * float(spec.solid.radius_mm) * float(spec.solid.length_mm)
    area = side.get("area_mm2")
    area = float(area) if area is not None else None
    return {"entrance": entrance, "exit": exit_face, "side": side,
            "side_area": {"area_mm2": area, "expected_mm2": expected,
                          "ratio": area / expected if area is not None else None}}


# --- the job -----------------------------------------------------------------------------------------------------------

def solved_project_kept_message(project: Path) -> str:
    """The run log's line when a step after a finished solve fails: the solved project and its solution folder stay in
    the job, the evidence of what AEDT solved."""
    return (f"the solve finished but a step after it failed: the solved project {project} and its solution folder "
            f"{project.with_suffix('.aedtresults')} are kept in the job; move the job aside rather than delete it")


def _read_json(path: Path) -> Any:
    try:
        return json.loads(Path(path).read_text())
    except (OSError, ValueError) as exc:
        log.warning("could not read %s: %s", path, exc)
        return None


def _keep_messages(hfss: Any, dirs: JobDirs, record: dict) -> None:
    """AEDT's message manager into the record and ``logs/aedt_messages.txt``; entries not seen before are logged, errors
    and warnings as warnings (AEDT's reason for a failure lives only there). Never raises."""
    try:
        messages = aedt_messages(hfss)
        for message in (m for m in messages if m not in record["aedt_messages"]):
            (log.warning if "[error]" in message or "[warning]" in message else log.info)("aedt: %s", message)
        record["aedt_messages"] = messages
        (dirs.logs / MESSAGES_NAME).write_text("".join(f"{m}\n" for m in messages))
    except Exception:  # noqa: BLE001 - diagnostics must not replace the run's own error
        log.exception("AEDT's messages could not be kept")


def _build(cfg: ModesConfig, project: Path, record: dict) -> None:
    """AEDT session 1: the bore, saved to ``project``. A failed session or a build whose own inventory differs from the
    spec stops the run (exit 3); a clean one verifies nothing (the reopen is skipped), so the run checks again."""
    t = time.perf_counter()
    try:
        outcome = build(cfg.spec, cfg.geometry_spec, project, cfg.spec.project.aedt_version, reopen=False)
    except Exception as exc:  # noqa: BLE001 - build() records a failed session itself; this is anything else
        log.exception("the build failed")
        raise _Stopped(3, "build", _text(exc)) from exc
    finally:
        record["timings"]["build_s"] = time.perf_counter() - t
    record["build"] = {"record": _read_json(project.with_suffix(".build.json")), "aedt_temp": asdict(BUILD_TEMP),
                       "launch_timeout_s": DEFAULT_LAUNCH_TIMEOUT_S, "verified_by_build": False,
                       "verification": "the run session's (the build skips its reopen)"}
    if outcome.versions:      # read first thing in the build's session: it opened
        record["launches"] += 1
    if outcome.error:
        raise _Stopped(3, "build", outcome.error)
    if outcome.diffs:
        raise _Stopped(3, "build", "the build's own inventory is not the spec's bore: " + "; ".join(outcome.diffs))
    log.info("built %s in %.1f s; the build's session is released and its lease settled", project,
             record["timings"]["build_s"])


def _setup(cfg: ModesConfig, hfss: Any, record: dict) -> None:
    t = time.perf_counter()
    record["versions"] = aedt_versions(hfss)
    log.info("versions: %s", record["versions"])
    record["solution_type"] = set_driven_modal(hfss)
    design = inventory_design(hfss)
    record["verification"] = {"design": design, "differences": None, "faces": None, "side_area": None}
    try:
        faces = bore_faces(design, cfg.spec)
    except BoreMismatch as exc:
        record["verification"]["differences"] = exc.differences
        raise
    entrance, exit_face, side = faces["entrance"], faces["exit"], faces["side"]
    record["verification"].update(differences=[], faces={"entrance": asdict(entrance), "exit": asdict(exit_face),
                                                         "side": side}, side_area=faces["side_area"])
    log.info("the opened design is the spec's bore: entrance face %d (z = 0, port1), exit face %d (z = %g mm, port2), "
             "side face %d (the wall)", entrance.id, exit_face.id, cfg.spec.solid.length_mm, side["id"])
    record["wall"] = assign_wall(hfss, side["id"], cfg.wall)
    record["mesh"] = assign_mesh(hfss, cfg, [entrance.id, exit_face.id])
    names = assign_ports(hfss, entrance.id, exit_face.id, cfg.ports.modes)
    record["ports"] = {"names": names, "modes": cfg.ports.modes,
                       "faces": {names[0]: entrance.id, names[1]: exit_face.id}}
    for f in cfg.ports.frequencies_ghz:
        setup = insert_ports_only_setup(hfss, f, cfg.solver.port_accuracy_percent)
        record["setups"].append({"name": setup, "frequency_ghz": f, "ports_only": True,
                                 "port_accuracy_percent": cfg.solver.port_accuracy_percent, "solution": None,
                                 "touchstone": None})
    hfss.save_project()
    record["timings"]["setup_s"] = time.perf_counter() - t
    log.info("wall %s, %d-mode ports %s, %d ports-only setups", record["wall"]["kind"], cfg.ports.modes, names,
             len(record["setups"]))


def _solve(cfg: ModesConfig, hfss: Any, dirs: JobDirs, record: dict) -> None:
    log.info("solving %d ports-only setups with %d cores", len(record["setups"]), cfg.solver.cores)
    t = time.perf_counter()
    solved = solve(hfss, cfg.solver.cores)
    record["timings"]["solve_s"] = time.perf_counter() - t
    log.info("solve finished in %.1f s", record["timings"]["solve_s"])
    _keep_messages(hfss, dirs, record)
    if not solved:
        raise RuntimeError("hfss.analyze returned False")
    # analyze() returning True proves nothing (HPC step 4, baseline1: a licence checkout failure solved nothing): every
    # setup must hold its one solution, found by the lookup the reads use, or the session failed (exit 3).
    missing = []
    for entry in record["setups"]:
        try:
            entry["solution"] = solution_of(hfss, entry["name"])
        except RuntimeError as exc:
            missing.append(f"{entry['name']}: {exc}")
        else:
            record["solutions"].append(entry["solution"])
    if missing:
        raise RuntimeError(f"the solve left {len(missing)} of {len(record['setups'])} setups without their solution: "
                           + "; ".join(missing) + ". AEDT messages: "
                           + (" | ".join(record["aedt_messages"]) or "(none)"))
    # AEDT writes the solution into the project's results folder only when the project is saved; a failed save is
    # exit 3, since the solution may then not be on disk.
    hfss.save_project()


def _read(cfg: ModesConfig, hfss: Any, dirs: JobDirs, record: dict) -> list[GammaRow]:
    t = time.perf_counter()
    rows: list[GammaRow] = []
    names = record["ports"]["names"]
    for entry in record["setups"]:
        setup = entry["name"]     # its solution was found in the solve stage
        got = read_gamma(hfss, setup, names, cfg.ports.modes, entry["frequency_ghz"])
        rows += got
        path = dirs.scratch / f"{setup}.s{len(names) * cfg.ports.modes}p"
        error = export_touchstone_with_gamma(hfss, setup, path)
        entry["touchstone"] = {"path": str(path), "error": error}
        if error is not None:
            log.warning("%s: the Touchstone export failed (evidence only, a warning): %s", setup, error)
        log.info("%s: %d port modes read from %s", setup, len(got), entry["solution"])
    record["timings"]["read_s"] = time.perf_counter() - t
    return rows


def _build_and_solve(cfg: ModesConfig, dirs: JobDirs, project: Path, temp: tempdir.TempSettings,
                     launch_timeout_s: float, record: dict) -> list[GammaRow]:
    """Both AEDT sessions; every port mode's Gamma. ``_Stopped`` with exit 3 before the solve finished, 1 after."""
    _build(cfg, project, record)
    stage = "session"
    try:
        with open_session(project, cfg.design, cfg.spec.project.aedt_version, dirs, temp,
                          launch_timeout_s) as session:
            record["launches"] += 1
            record["aedt"], record["aedt_temp"] = session.aedt, session.aedt_temp   # restore completes the latter
            try:
                stage = "setup"
                _setup(cfg, session.hfss, record)
                stage = "solve"
                _solve(cfg, session.hfss, dirs, record)
                stage = "read"
                rows = _read(cfg, session.hfss, dirs, record)
                _keep_messages(session.hfss, dirs, record)   # a failed export's reason is only in AEDT's messages
                session.before_write()    # the run's AEDT calls are done: AEDT's temp directory goes back now
                stage = "done"
            except Exception:
                _keep_messages(session.hfss, dirs, record)
                raise
    except Exception as exc:  # noqa: BLE001 - mapped to the exit code below; the session is already released
        code = 1 if stage in ("read", "done") else 3
        log.exception("the modes run failed at %s (exit %d)", stage, code)
        if code == 1:
            log.error(solved_project_kept_message(project))
        raise _Stopped(code, stage, _text(exc)) from exc
    return rows


def _touchstone_row(setups: Sequence[Mapping[str, Any]]) -> CheckRow:
    failed = [f"{s['name']}: {s['touchstone']['error']}" for s in setups if (s["touchstone"] or {}).get("error")]
    if not failed:
        return CheckRow("export.touchstone", True, "warn",
                        f"{len(setups)} Touchstone exports with Gamma in scratch/ (evidence)")
    more = len(failed) - _NAMED
    return CheckRow("export.touchstone", False, "warn",
                    "; ".join(failed[:_NAMED]) + (f"; and {more} more" if more > 0 else ""))


def _write(cfg: ModesConfig, dirs: JobDirs, rows: list[GammaRow], record: dict, environment: EnvironmentState) -> None:
    """``modes.csv``, recorded in the manifest's ``outputs`` as soon as it is on disk."""
    path = dirs.root / CSV_NAME
    write_modes_csv(path, rows, csv_provenance(cfg, environment, record["versions"] or {}))
    record["outputs"] = {"modes_csv": {"path": str(path), "sha256": sha256_file(path), "rows": len(rows)}}
    log.info("wrote %s (%d rows)", path, len(rows))


def _check(cfg: ModesConfig, rows: list[GammaRow], record: dict) -> list[CheckRow]:
    checks = self_check(rows, cfg.radius_m, cfg.ports.modes, cfg.ports.frequencies_ghz, cfg.wall.kind)
    checks.append(_touchstone_row(record["setups"]))
    for row in checks:
        verdict = "PASS" if row.passed else ("FAIL" if row.severity == "fail" else "WARN")
        (log.info if row.passed else log.warning)("check %s: %s (%s)", row.name, verdict, row.detail)
    return checks


def _project_sha256(project: Path) -> str | None:
    """The project's sha256, or None when it is absent or cannot be read (a manifest is written either way)."""
    try:
        return sha256_file(project) if project.is_file() else None
    except Exception:  # noqa: BLE001 - recorded as null; the run's own failure keeps its exit code
        log.exception("the project %s could not be hashed", project)
        return None


def _write_manifest(cfg: ModesConfig, inputs: Mapping[str, Any], dirs: JobDirs, job_id: str, project: Path,
                    record: dict, checks: Sequence[CheckRow], failure: _Stopped | None, code: int,
                    launch_timeout_s: float, environment: EnvironmentState, guard: Mapping[str, Any]) -> bool:
    """Build and write ``modes_manifest.json``; False, logged, if either fails (the caller keeps the exit code the run
    earned, or makes a pass a failure)."""
    try:
        build_record = (record["build"] or {}).get("record") or {}
        manifest = build_manifest(
            environment=environment,
            acceptance=acceptance_block(cfg),
            job={"id": job_id, "root": str(dirs.root)},
            **inputs,
            config=config_record(cfg),
            spec=asdict(cfg.spec),
            code={"commit": environment.commit, "clean": not environment.dirty},
            guard=dict(guard),
            launches=record["launches"],
            build=record["build"],
            project={"path": str(project), "sha256": _project_sha256(project),
                     "sha256_after_build": build_record.get("output_sha256")},
            verification=record["verification"], solution_type=record["solution_type"], wall=record["wall"],
            mesh=record["mesh"], ports=record["ports"], setups=record["setups"], solutions=record["solutions"],
            versions=record["versions"], aedt=record["aedt"], aedt_temp=record["aedt_temp"],
            launch_timeout_s=launch_timeout_s, aedt_messages=record["aedt_messages"], timings=record["timings"],
            outputs=record["outputs"], self_check=[row._asdict() for row in checks], passed=code == 0,
            exit_code=code,
            failure=None if failure is None else {"stage": failure.stage, "exit_code": failure.code,
                                                  "error": failure.error},
        )
        write_json(dirs.root / MANIFEST_NAME, manifest)
    except Exception:  # noqa: BLE001 - a run without its manifest cannot pass
        log.exception("%s could not be written", MANIFEST_NAME)
        return False
    return True


def run_job(cfg: ModesConfig, inputs: Mapping[str, Any], dirs: JobDirs, job_id: str, *, temp: tempdir.TempSettings,
            launch_timeout_s: float, environment: EnvironmentState, guard: Mapping[str, Any]) -> int:
    """The job in its directory: both sessions, the CSV and its self-checks, the manifest (written however the run
    ends, SIGTERM aside), and the last line. Returns the exit code."""
    t_start = time.perf_counter()
    project = dirs.project / f"{cfg.dataset_id}.aedt"
    record: dict[str, Any] = {
        "launches": 0, "build": None, "verification": None, "solution_type": None, "wall": None, "mesh": None,
        "ports": None, "setups": [], "solutions": [], "versions": None, "aedt": None, "aedt_temp": None,
        "aedt_messages": [], "outputs": None,
        "timings": {"wall_s": None, "build_s": None, "setup_s": None, "solve_s": None, "read_s": None}}
    failure: _Stopped | None = None
    checks: list[CheckRow] = []
    try:
        rows = _build_and_solve(cfg, dirs, project, temp, launch_timeout_s, record)
    except _Stopped as stop:
        failure = stop
    else:
        stage = "write"
        try:
            _write(cfg, dirs, rows, record, environment)
            stage = "check"
            checks = _check(cfg, rows, record)
        except Exception as exc:  # noqa: BLE001 - the solve is kept; exit 1
            log.exception("%s %s failed", "writing" if stage == "write" else "self-checking", CSV_NAME)
            log.error(solved_project_kept_message(project))
            failure = _Stopped(1, stage, _text(exc))
    failed_rows = [c.name for c in checks if c.severity == "fail" and not c.passed]
    code = failure.code if failure is not None else (1 if failed_rows else 0)
    record["timings"]["wall_s"] = time.perf_counter() - t_start
    if _write_manifest(cfg, inputs, dirs, job_id, project, record, checks, failure, code, launch_timeout_s,
                       environment, guard):
        log.info("wrote %s", dirs.root / MANIFEST_NAME)
    else:
        code = code or 1
    if code == 0:
        line = "== modes: PASS"
    elif failure is not None:
        line = f"== modes: FAIL {failure.stage} (exit {failure.code}): " + " ".join(failure.error.split())
    else:
        line = "== modes: FAIL " + (", ".join(failed_rows) or f"the manifest ({MANIFEST_NAME}) was not written")
    (log.info if code == 0 else log.error)("%s", line)
    print(line)
    return code


def main(argv: Sequence[str] | None = None) -> int:
    """The CLI (see the module docstring). For the length of the run SIGTERM raises ``SystemExit(143)``, and the run's
    log handlers are removed and closed on the way out, so ``main`` can run twice in one process (HF-057)."""
    with ExitStack() as cleanup:
        cleanup.callback(_handle_sigterm())
        try:
            args = parse_args(argv)
        except SystemExit as exc:   # argparse: --help (0) or a refused argument (2), before anything started
            return exc.code if isinstance(exc.code, int) else 2
        try:
            cfg = load_modes_config(args.config)
            inputs = {"config_file": {"path": str(args.config.resolve()), "sha256": sha256_file(args.config)},
                      "spec_file": {"path": str(cfg.geometry_spec.resolve()), "sha256": sha256_file(cfg.geometry_spec)}}
        except (OSError, ValueError) as exc:
            print(f"refused before the job started: {args.config}: {exc}", file=sys.stderr)
            return 2
        summary = describe(cfg, args.config)
        print("\n".join(summary))
        if args.dry_run:
            print("\nDry run: no job directory, no AEDT session.")
            return 0
        job_id = args.job_id or default_job_id()
        try:
            guard = check_guard_before_job()
            environment = environment_state(REPO_ROOT)
            dirs = create_modes_job_dirs(cfg.output.root, cfg.dataset_id, job_id)
        except (OSError, ValueError, ImportError, licence_guard.GuardError, EnvironmentStateError) as exc:
            print(f"refused before the job started: {exc}", file=sys.stderr)
            return 2
        cleanup.callback(_setup_logging(dirs.logs / "run.log", "INFO"))
        write_json(dirs.root / "effective_config.json", {**inputs, "config": config_record(cfg),
                                                          "spec": asdict(cfg.spec)})
        log.info("job %s -> %s: the fibre bore's ports-only mode run (OR-053)", job_id, dirs.root)
        for line in summary:
            log.info("%s", line)
        log.info(_guard_log_line(guard))
        (log.warning if environment.dirty else log.info)(_environment_log_line(environment))
        log.info("AEDT temp mode for the run's session: %s; the build's session sets none", args.temp_settings.mode)
        warning = allocation_warning(cfg.solver.cores)
        if warning:
            log.warning(warning)
        try:
            return run_job(cfg, inputs, dirs, job_id, temp=args.temp_settings, launch_timeout_s=args.launch_timeout_s,
                           environment=environment, guard=guard)
        except SystemExit as exc:
            _log_system_exit("modes job", exc)
            raise
        except Exception:  # noqa: BLE001 - report and exit non-zero; the sessions are already released
            log.exception("modes job failed")
            print("== modes: FAIL the run itself raised; see logs/run.log")
            return 1


if __name__ == "__main__":
    raise SystemExit(main())
