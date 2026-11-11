"""CLI: solve one frequency headlessly from a config file and export the two CSV tables.

Two phases share one AEDT session: ``prepare_and_solve`` sets up the run design in the job's own project copy and
solves it; ``extract_and_write`` extracts both polarizations into memory, checks them, and only then writes the
CSVs, the manifest and the dataset sidecar (HF-050, HF-052). ``output.schema`` selects the
writer: 2 by default (HF-050), 1 for the 1.0 tables until the ``_v2`` trees pass BBRsim's tree check.

``--from-solved JOB_DIR`` replaces the first phase: a finished job's project and ``.aedtresults`` are copied into a
new job beside it, the solved run design is opened, and ``extract_and_write`` writes schema 2 from the copy with no
solve, under the source manifest's configuration (HF-050).

Every run (HF-057): the code checkout's state is read once before the job starts; the session opens with PyAEDT's
launch timeout, must have opened the job's own copy, and records AEDT's process id and port; AEDT's temp directory is
recorded and, by ``--aedt-temp``, pointed at the job (``bbsim.aedt_temp``); ``main`` owns its log handlers and turns
SIGTERM into ``SystemExit(143)`` so the cleanup on the way out runs.

Plan mode, one job-array task (HF-057, OR-037, OR-045, HF-058): ``--plan PLAN --plan-sha256 HEX [--plan-index N]
[--log-level L]`` and no other option. ``plan_preflight`` checks everything before a job directory or AEDT exists, and
any refusal exits 2: the plan's bytes have the given sha256; the index (by default ``SLURM_ARRAY_TASK_ID``, which must
agree with ``--plan-index``) is a task of the plan; the plan is an ``array_run`` plan of its policy; its code is the
checkout this module was imported from, clean at the plan's commit, with the plan's PyAEDT; the config, project,
mapping record and (production) column-order reference copies have the plan's sha256, the record beside the project
copy confirms the plan's dataset, and the reference fits this checkout's exporter and PyAEDT; the task's config
(``array_plan.task_config`` of the captured config bytes) is the policy's and passes ``check_before_aedt``; its
effective step along l is at most c/(4f) (OR-048) and its recorded lattice is the frozen rectangle's HF-045 lattice
(``array_plan.predict_exit_lattice``); ``SLURM_CPUS_PER_TASK``, when set, is at least the plan's cores;
``PYAEDT_DESKTOP_PORT`` and ``PYAEDT_PROCESS_ID`` are unset; the plan's guard is the shared one, named by
``BBSIM_LICENCE_GUARD``, and an array task runs under its own ``array`` reservation; the job directory does not exist.
The job is ``a<SLURM_ARRAY_JOB_ID>_<index>`` in an array, else ``<default_job_id>_<index>``. ``start.json`` is
written before AEDT starts; ``run_job`` then compares the exit face's lattice with the plan's prediction before the
solve and records the plan block in the manifest; ``bbsim.run_acceptance`` judges the run against its plan; once its
report is written, ``done.json`` records the completion and its verdict, PASS or a deterministic FAIL. The log ends
``done; manifest at <path>`` (exit 0) or ``completed; acceptance FAIL (<rows>); manifest at <path>`` (exit 1); a run
that stops before the report writes no ``done.json`` (exit 1). A launch lease the guard refuses for want of capacity
is named as a refusal; any other guard failure before AEDT starts is named as a fault.
"""
from __future__ import annotations

import argparse
import copy
import functools
import hashlib
import json
import logging
import math
import os
import re
import shutil
import signal
import socket
import subprocess
import sys
import threading
import time
import tomllib
from collections.abc import Callable, Iterator, Mapping, Sequence
from contextlib import ExitStack, contextmanager
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from bbsim import aedt_temp as tempdir
from bbsim import array_plan, dataset_v2, frequency_grid, licence_guard
from bbsim.config import RunConfig, build_config, config_to_dict, load_config
from bbsim.dataset_record import REEXPORT_TOOL, build_dataset_record, build_dataset_record_v2, modes_record_v2
from bbsim.dataset_v2 import DatasetV2Error
from bbsim.geometry_spec import is_vacuum, vacuum_refusal
from bbsim.manifest import (REPO_ROOT, EnvironmentState, EnvironmentStateError, build_manifest, environment_state,
                            sha256_file, write_json)
from bbsim.naming import dataset_record_name, frequency_design_name, frequency_label, last_adaptive_solution
from bbsim.paths import JobDirs, create_job_dirs, default_job_id, job_root
from bbsim.sampling import incident_angle_values, incoming_power_w
from bbsim.schema import KEY_COLUMNS
from bbsim.session import DEFAULT_LAUNCH_TIMEOUT_S, launch_timeout_setting

log = logging.getLogger("bbsim.run")

LOG_FORMAT = "%(asctime)s %(levelname)s %(name)s: %(message)s"
EXIT_SIGTERM = 143   # 128 + SIGTERM, as a shell reports a process ended by it

# Plan mode (HF-057): the options it takes, the variables that would attach PyAEDT to another desktop, its records
PLAN_MODE_OPTIONS = ("plan", "plan_sha256", "plan_index", "log_level")
ATTACH_VARIABLES = ("PYAEDT_DESKTOP_PORT", "PYAEDT_PROCESS_ID")
START_NAME, DONE_NAME = "start.json", "done.json"
START_SLURM_VARIABLES = ("SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID", "SLURM_JOB_ID", "SLURMD_NODENAME",
                         "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE", "SLURM_RESTART_COUNT")
SCONTROL_TIMEOUT_S = 30
_SHA256 = re.compile(r"[0-9a-f]{64}")
_DIGITS = re.compile(r"[0-9]+")
_UNSET = object()


def _task_index_arg(text: str) -> int:
    if not _DIGITS.fullmatch(text):
        raise argparse.ArgumentTypeError(f"a task index is a nonnegative integer, got {text!r}")
    return int(text)


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Headless HFSS gap simulation for one frequency.")
    p.add_argument("--plan", type=Path, default=None, metavar="PLAN",
                   help="Plan mode (HF-057): run one task of this frozen job-array plan.json; takes --plan-sha256, "
                        "--plan-index and --log-level only, since the plan supplies everything else")
    p.add_argument("--plan-sha256", default=None, metavar="HEX",
                   help="plan mode: the sha256 the submission recorded for the plan (plan.sha256); the plan's bytes "
                        "must have it")
    p.add_argument("--plan-index", type=_task_index_arg, default=None, metavar="N",
                   help="plan mode: the task's index in the plan; by default SLURM_ARRAY_TASK_ID, which must agree "
                        "when both are given")
    p.add_argument("--config", type=Path, default=None,
                   help="TOML run configuration (see configs/); required unless --from-solved")
    p.add_argument("--from-solved", type=Path, default=None, metavar="JOB_DIR",
                   help="Schema 2 re-export with no solve: copy this finished job's project and .aedtresults into a "
                        "new job beside it (--job-id required) and extract from the copy; the configuration is the "
                        "job's manifest's, so --config and the overrides below are refused")
    p.add_argument("--frequency-ghz", type=float, default=None)
    p.add_argument("--cores", type=int, default=None)
    p.add_argument("--output-root", default=None)
    p.add_argument("--project", default=None, help="Prepared .aedt project path")
    p.add_argument("--polarizations", type=int, nargs="+", default=None, choices=[0, 1])
    p.add_argument("--aedt-version", default=None)
    p.add_argument("--job-id", default=None, help="Defaults to slurm<SLURM_JOB_ID> or a UTC timestamp")
    p.add_argument("--column-order-reference", type=Path, default=None, metavar="JSON",
                   help="Schema 2: the exporter column-order reference (HF-041) that the crack1 500 GHz reference "
                        "job keeps as scratch/column_order_reference.json; required by every other schema-2 job")
    p.add_argument("--aedt-temp", default="off", metavar="{off,task,env,TEMPLATE}",
                   help="AEDT's temp directory (HF-057): off leaves it; task points it at <job>/aedt_temp and sets it "
                        "to --aedt-temp-restore at the end; a path holding {job} (the job directory) is task mode "
                        "there; env exports TMPDIR and TMP = <job>/aedt_temp before AEDT starts. Every run records "
                        "what AEDT reported at session start")
    p.add_argument("--aedt-temp-restore", default=None, metavar="PATH",
                   help="task mode: the TempDirectory AEDT is set back to at the end, the host's TempDirectory "
                        "before any plan-mode run, read from AEDT's .cfg by the cluster probe")
    p.add_argument("--aedt-temp-strict", action="store_true",
                   help="task mode: fail the run instead of warning when AEDT does not report the requested "
                        "directory or a temp-directory call raises")
    p.add_argument("--launch-timeout-s", type=float, default=DEFAULT_LAUNCH_TIMEOUT_S,
                   help="seconds PyAEDT waits for AEDT to start (its settings.desktop_launch_timeout)")
    p.add_argument("--log-level", default="INFO")
    args = p.parse_args(argv)
    if not isinstance(logging.getLevelName(args.log_level.upper()), int):
        p.error(f"--log-level {args.log_level!r} is not a logging level")
    if args.plan is not None:
        # Parsed again into a namespace pre-filled with a marker, so an option given at its default value, or
        # abbreviated, is seen as given.
        given = vars(p.parse_args(argv, namespace=argparse.Namespace(**dict.fromkeys(vars(args), _UNSET))))
        others = sorted("--" + dest.replace("_", "-") for dest, value in given.items()
                        if value is not _UNSET and dest not in PLAN_MODE_OPTIONS)
        if others:
            p.error(f"--plan takes only --plan-sha256, --plan-index and --log-level, since the plan supplies "
                    f"everything else (HF-057); drop {', '.join(others)}")
        if args.plan_sha256 is None:
            p.error("--plan needs --plan-sha256, the plan's sha256 the submission recorded")
        if not _SHA256.fullmatch(args.plan_sha256):
            p.error(f"--plan-sha256 must be 64 lowercase hex digits, got {args.plan_sha256!r}")
    elif args.plan_sha256 is not None or args.plan_index is not None:
        p.error("--plan-sha256 and --plan-index belong to --plan")
    if not (math.isfinite(args.launch_timeout_s) and args.launch_timeout_s > 0):
        p.error(f"--launch-timeout-s must be a finite number of seconds above 0, got {args.launch_timeout_s!r}")
    try:
        mode, template = tempdir.parse_mode(args.aedt_temp)
        args.temp_settings = tempdir.TempSettings(mode, restore=args.aedt_temp_restore, strict=args.aedt_temp_strict,
                                                  template=template)
    except ValueError as exc:
        p.error(str(exc))
    return args


def overrides_from_args(args: argparse.Namespace) -> dict[str, object]:
    return {
        "solver.frequency_ghz": args.frequency_ghz,
        "solver.cores": args.cores,
        "output.root": args.output_root,
        "project.path": args.project,
        "output.polarizations": args.polarizations,
        "project.aedt_version": args.aedt_version,
    }


def allocation_warning(cores: int, env: Mapping[str, str] | None = None) -> str | None:
    """Message when the configured solver cores disagree with the Slurm allocation, else None."""
    env = os.environ if env is None else env
    allocated = env.get("SLURM_CPUS_PER_TASK")
    if allocated and allocated.isdigit() and int(allocated) != cores:
        return f"solver.cores = {cores} but SLURM_CPUS_PER_TASK = {allocated}; pass --cores to match the allocation"
    return None


def _axis_summary(values, step: float) -> dict[str, float | int]:
    """The emitted far-field grid along one angle: first and last value, point count, step (degrees)."""
    return {"min": float(values[0]), "max": float(values[-1]), "count": int(len(values)), "step": float(step)}



def _setup_logging(logfile: Path, level: str) -> Callable[[], None]:
    """Attach this run's handlers, ``logfile`` (the job's ``logs/run.log``) and the console, to the root logger at
    ``level``, and return the function that removes and closes them and puts the root level back (``main`` calls it on
    the way out). ``logging.basicConfig`` does nothing once the root logger has a handler: under pytest, and in a second
    ``main()`` of one process."""
    root = logging.getLogger()
    previous_level = root.level
    handlers = [logging.FileHandler(logfile), logging.StreamHandler()]
    for handler in handlers:
        handler.setFormatter(logging.Formatter(LOG_FORMAT))
        root.addHandler(handler)
    root.setLevel(level.upper())

    def remove() -> None:
        for handler in handlers:
            root.removeHandler(handler)
            handler.close()
        root.setLevel(previous_level)

    return remove


def _exit_on_sigterm(signum: int, frame: Any) -> None:
    raise SystemExit(EXIT_SIGTERM)


def _handle_sigterm() -> Callable[[], None]:
    """For the length of ``main``, SIGTERM (Slurm's time limit, scancel) raises ``SystemExit(143)``, so the blocks on
    the way out run when Python regains control: the temp-directory restore, the session's release and launch-lease
    wait, the log handlers' removal. Best effort: a handler runs only between bytecodes, so during a blocking AEDT call
    it waits, and Slurm's SIGKILL follows after KillWait. Returns the function that puts the previous handler back;
    outside the main thread, where Python allows no handler, nothing changes."""
    if threading.current_thread() is not threading.main_thread():
        return lambda: None
    previous = signal.signal(signal.SIGTERM, _exit_on_sigterm)
    return lambda: signal.signal(signal.SIGTERM, signal.SIG_DFL if previous is None else previous)


@contextmanager
def aedt_messages_on_error(hfss: Any, seen: set[str]) -> Iterator[None]:
    """On an exception in the block, log AEDT's message-manager entries not in ``seen``, then re-raise.

    PyAEDT turns a failed AEDT command into a bare "Failed to execute gRPC AEDT command"; AEDT's reason is
    only in its message manager (HPC step 4, baseline2).
    """
    from bbsim.hfss_setup import aedt_messages

    try:
        yield
    except Exception:
        new = [m for m in aedt_messages(hfss) if m not in seen]
        for message in new or ["(no new AEDT messages)"]:
            log.warning("aedt: %s", message)
        raise


@dataclass(frozen=True)
class StagedProject:
    """The job's own project copy, made before any AEDT session, and what was checked and hashed on the way."""

    source: Path
    source_sha: str
    job_project: Path
    job_sha: str
    recorded_id: str | None
    record: Path | None


@dataclass(frozen=True)
class SolvedContext:
    """What ``prepare_and_solve`` (or a reuse of a solved job) hands to ``extract_and_write``: the extraction context,
    the exit face and its lattice, the names of everything created in the run design, the solve's evidence, and the
    project provenance. ``solve_s`` is None when no solve ran."""

    ctx: Any
    face_geometry: Any
    lattice: Any
    entrance: Any
    exit_face: Any
    solid: dict
    setup: str
    parametric: str
    sphere: str
    grid: Any
    plane_wave: str
    wall_facets: dict | None
    mesh_inherited: Any
    mesh_requested: dict
    variations: list
    messages: list
    solve_s: float | None
    source: Path
    source_sha: str
    job_project: Path
    job_sha: str
    recorded_id: str | None
    record: Path | None
    run_design: str
    versions: dict


@dataclass(frozen=True)
class OpenSession:
    """A run's open AEDT session: the app, the manifest's ``aedt`` {pid, port} and ``aedt_temp`` records, and
    ``before_write``, which restores AEDT's temp directory once the run's last AEDT call is made, so the manifest
    records the outcome (HF-057)."""

    hfss: Any
    aedt: dict
    aedt_temp: dict
    before_write: Callable[[], None]


def _positive_int(value: Any) -> int | None:
    return value if isinstance(value, int) and not isinstance(value, bool) and value > 0 else None


def check_open_project(hfss: Any, job_project: Path) -> dict[str, int | None]:
    """After the session opens, the open project must be the job's own copy (``os.path.samefile``: ``/home`` resolves
    to ``/resnick/home`` on the cluster), or the run stops before anything is set up. Returns the ``aedt`` record:
    AEDT's process id and gRPC port as PyAEDT reports them (``hfss.desktop_class``), each null when PyAEDT gives no
    positive integer (HF-057)."""
    opened = getattr(hfss, "project_file", None)
    try:
        same = bool(opened) and os.path.samefile(opened, job_project)
    except OSError:
        same = False
    if not same:
        raise RuntimeError(f"AEDT opened {opened!r}, not the job's own copy {job_project}; nothing was set up")
    desktop = getattr(hfss, "desktop_class", None)
    pid, port = getattr(desktop, "aedt_process_id", None), getattr(desktop, "port", None)
    log.info("AEDT process %r, gRPC port %r, open project %s", pid, port, opened)
    return {"pid": _positive_int(pid), "port": _positive_int(port)}


@contextmanager
def open_session(job_project: Path, design: str, version: str, dirs: JobDirs, temp: tempdir.TempSettings,
                 launch_timeout_s: float) -> Iterator[OpenSession]:
    """The run's one AEDT session (HF-057): in env mode TMPDIR and TMP are exported first; the session starts with
    the launch timeout; the open project is checked to be the job's copy; AEDT's temp directory is recorded and, in
    task mode, pointed at the job and restored on the way out, however the block ends."""
    from bbsim.session import aedt_session

    with tempdir.exported_env(temp, dirs.root):
        with aedt_session(job_project, design, version, launch_timeout_s=launch_timeout_s) as hfss:
            aedt = check_open_project(hfss, job_project)
            with tempdir.around_session(hfss, temp, dirs.root) as record:
                yield OpenSession(hfss=hfss, aedt=aedt, aedt_temp=record,
                                  before_write=functools.partial(tempdir.restore, hfss, temp, record))


def is_reference_job(cfg: RunConfig) -> bool:
    """True for the one dataset and frequency that measures the exporter column-order reference itself (HF-041)."""
    return (cfg.project.dataset_id == dataset_v2.REFERENCE_DATASET
            and float(cfg.solver.frequency_ghz) == dataset_v2.REFERENCE_GHZ)


def load_reference(path: Path) -> dict:
    """The column-order reference object (HF-041) from a JSON file, such as the reference job's
    ``scratch/column_order_reference.json``. Its content is validated against the run's versions later."""
    reference = json.loads(Path(path).read_text())
    if not isinstance(reference, dict):
        raise ValueError(f"{path} must hold the column-order reference object, got {type(reference).__name__}")
    return reference


def check_before_aedt(cfg: RunConfig, reference: Mapping | None) -> None:
    """Refuse, from the configuration alone and before any job copy or AEDT session, a run whose output cannot be
    written: under schema 1 a column-order reference (schema 1 neither needs nor applies it); under schema 2
    anything but both polarizations, a declared theta grid without a 90 deg key and one above it (``cell_90`` is
    defined by that band, HF-049 and OR-036), and a missing reference on any job but the reference job (HF-041).
    Under both schemas, first a run with no frequency (solver.frequency_ghz may be absent from a config, for a
    job-array plan to supply; HF-057, OR-047), then an acceptance reference_job that is not an existing directory: the
    job that is to judge this one must exist before a licensed solve (CDX-039)."""
    if cfg.solver.frequency_ghz is None:
        raise ValueError("no frequency: give solver.frequency_ghz in the config or --frequency-ghz; a job-array task "
                         "takes its frequency from its plan (--plan)")
    reference_job = cfg.acceptance.reference_job
    if reference_job is not None and not Path(reference_job).is_dir():
        raise ValueError(f"[acceptance] reference_job {reference_job} is not an existing job directory: the "
                         f"{cfg.acceptance.kind} that is to judge this job needs it (CDX-039)")
    if cfg.output.schema == 1:
        if reference is not None:
            raise ValueError("a column-order reference applies to schema 2 only; this run writes schema 1 "
                             "(output.schema = 1), which neither requires nor applies it")
        return
    if sorted(cfg.output.polarizations) != [0, 1]:
        raise DatasetV2Error(f"schema 2 needs both polarizations (output.polarizations = [0, 1]), got "
                             f"{list(cfg.output.polarizations)}; a one-polarization diagnostic writes schema 1")
    exc = cfg.excitation
    theta_values = incident_angle_values(exc.theta_lower_deg, exc.theta_upper_deg, exc.theta_step_deg)
    try:
        dataset_v2.cell_90_band_deg(theta_values)
    except DatasetV2Error as err:
        raise DatasetV2Error(f"schema 2: {err}; a diagnostic without that key writes schema 1") from err
    if reference is None and not is_reference_job(cfg):
        raise DatasetV2Error(
            f"schema 2 needs the column-order reference (--column-order-reference JSON): only "
            f"{dataset_v2.REFERENCE_DATASET} at {dataset_v2.REFERENCE_GHZ:g} GHz measures its own, and "
            f"{cfg.project.dataset_id} at {float(cfg.solver.frequency_ghz):g} GHz must embed that one (HF-041)")
    if reference is None:
        # The reference job measures its two cases at fixed keys; a grid without them would fail after the solve.
        phi_values = incident_angle_values(exc.phi_lower_deg, exc.phi_upper_deg, exc.phi_step_deg)
        missing = [(phi, theta) for _, phi, theta in dataset_v2.REFERENCE_CASES
                   if not (np.any(np.abs(phi_values - phi) <= 1e-9) and np.any(np.abs(theta_values - theta) <= 1e-9))]
        if missing:
            raise DatasetV2Error(f"the column-order reference job needs the keys (phi, theta) = "
                                 f"{[c[1:] for c in dataset_v2.REFERENCE_CASES]} deg in its grid; missing {missing}")
    if reference is not None and not isinstance(reference, Mapping):
        raise DatasetV2Error(f"the column-order reference must be an object, got {type(reference).__name__}")


def stage_project(cfg: RunConfig, dirs: JobDirs, *, require_record: bool = False) -> StagedProject:
    """Check the config's dataset_id against the record beside the prepared project, then copy the project into the
    job directory; both before any AEDT session, which opens the copy. A single run only warns when no record
    confirms the dataset; ``require_record`` (a plan's task) refuses then (HF-057)."""
    from bbsim.mapping import recorded_dataset_id

    source = cfg.project.path
    if not source.is_file():
        raise FileNotFoundError(f"prepared project not found: {source}")
    # The run config names the dataset; the record beside the project says which dataset its design is.
    recorded_id, record = recorded_dataset_id(source, cfg.project.design, cfg.geometry.object)
    if recorded_id is not None and recorded_id != cfg.project.dataset_id:
        raise ValueError(f"config dataset_id {cfg.project.dataset_id!r} disagrees with {record}, which gives "
                         f"{recorded_id!r} for design {cfg.project.design!r} object {cfg.geometry.object!r}")
    if recorded_id is None:
        why = (f"{record} has no entry for {cfg.project.design}/{cfg.geometry.object}" if record
               else f"no .mapping.json or .build.json beside {source}")
        if require_record:
            raise ValueError(f"dataset_id {cfg.project.dataset_id!r} not confirmed: {why}; a job-array task runs only "
                             "a project its frozen record confirms (HF-057)")
        log.warning("dataset_id %r not confirmed: %s", cfg.project.dataset_id, why)
    job_project = dirs.project / source.name
    shutil.copy2(source, job_project)
    return StagedProject(source=source, source_sha=sha256_file(source), job_project=job_project,
                         job_sha=sha256_file(job_project), recorded_id=recorded_id, record=record)


def run_job(cfg: RunConfig, dirs: JobDirs, job_id: str, reference: dict | None = None, *,
            aedt_temp: tempdir.TempSettings | None = None, launch_timeout_s: float = DEFAULT_LAUNCH_TIMEOUT_S,
            environment: EnvironmentState | None = None, expected_lattice: Mapping | None = None,
            plan_record: Mapping | None = None) -> Path:
    """One frequency: stage the project, open one AEDT session, solve, extract both polarizations, write.

    ``reference`` is the column-order reference object every schema-2 job but the reference job embeds. ``aedt_temp``
    (default off) and ``launch_timeout_s`` set the session (HF-057); ``environment`` is the code checkout's state,
    read here before anything is copied when the caller has not read it. A job-array task (plan mode) also gives
    ``expected_lattice``, its plan's predicted exit lattice, which the exit face must give before the solve, and
    ``plan_record`` {path, sha256, index}, the manifest's ``plan`` block; it runs only a project its record
    confirms."""
    from bbsim.session import aedt_versions

    t_start = time.perf_counter()
    check_before_aedt(cfg, reference)
    environment = environment_state(REPO_ROOT) if environment is None else environment
    temp = tempdir.TempSettings() if aedt_temp is None else aedt_temp
    staged = stage_project(cfg, dirs, require_record=plan_record is not None)
    with open_session(staged.job_project, cfg.project.design, cfg.project.aedt_version, dirs, temp,
                      launch_timeout_s) as session:
        hfss = session.hfss
        if cfg.output.schema == 2:
            # A stale or unreadable version would stop the run after its licensed solve; it costs seconds here.
            _check_versions_for_reference(reference, aedt_versions(hfss))
        solved = prepare_and_solve(cfg, dirs, job_id, hfss, staged, expected_lattice=expected_lattice)
        try:
            return extract_and_write(cfg, dirs, job_id, hfss, solved, t_start, reference, None,
                                     environment=environment, session=session, plan_record=plan_record)
        except Exception:
            # The licensed solve is saved in the job's copy (prepare_and_solve) and nothing here removes it; say so,
            # so that a failed check after it never ends in deleting the job and solving again.
            log.error(solved_project_kept_message(solved.job_project))
            raise


def solved_project_kept_message(job_project: Path) -> str:
    """The run log's line when a check after a finished solve fails: it names the solved project and its solution
    folder, which stay in the job, and says why they must not be deleted."""
    return (f"the solve finished but a check after it failed: the solved project {job_project} and its solution "
            f"folder {job_project.with_suffix('.aedtresults')} are kept; do not delete this job (move it aside), "
            f"because a re-extraction can reuse the solve instead of solving again")


def _require_vacuum(solid: Mapping[str, Any]) -> None:
    """CDX-029: the sidecar declares vacuum-guide cutoffs and BBRsim refuses any other material (BBR025), so a gap
    solid that is not vacuum is refused under both schemas before the solve or any extraction. The material is
    never relabelled."""
    if not is_vacuum(solid["material"]):
        raise RuntimeError(f"gap solid {solid['name']!r}: {vacuum_refusal(solid['material'])}")


def _require_v2_exit_normal(cfg: RunConfig, exit_face: Any) -> None:
    """Schema 2: exit_flux integrates ScalarZ of Re S, which is the flux Re S.n (HF-039(1), CDX-009) only on an exit
    face whose outward normal is +z; any other pose is refused, a pose the config leaves unchecked included."""
    if tuple(exit_face.outward_normal) != (0.0, 0.0, 1.0):
        raise DatasetV2Error(
            f"schema 2 stores the exit flux Integrate(ScalarZ(Re S)), which is Re S.n only when the exit face's "
            f"outward normal is (0, 0, 1) (pose canonical-z); exit face {exit_face.id} has outward normal "
            f"{list(exit_face.outward_normal)} (config pose {cfg.geometry.pose!r})")


def _require_v2_lattice(cross_section: Mapping, lattice: Any) -> None:
    """Schema 2: the quadrature weights and the modes need a rectangle or a disc on the manual exit lattice."""
    shape = cross_section["shape"]
    if shape not in dataset_v2.QUADRATURE_BY_SHAPE or lattice is None:
        raise DatasetV2Error(f"schema 2 needs a rectangle or disc exit cross-section on the manual exit lattice "
                             f"(quadrature weights, modes); got {shape!r}"
                             + ("" if lattice is not None else " without a lattice"))


def _require_v2_exit_face(cfg: RunConfig, solved: SolvedContext) -> None:
    """Both schema-2 exit-face checks on a solved context, whichever path built it (fresh solve or --from-solved)."""
    _require_v2_exit_normal(cfg, solved.exit_face)
    _require_v2_lattice(solved.face_geometry.cross_section, solved.lattice)


def _check_versions_for_reference(reference: Mapping | None, versions: Mapping[str, Any]) -> None:
    """Schema 2: the reference job's versions must be evidence; a supplied reference must match them (HF-041)."""
    aedt, pyaedt = versions.get("aedt"), versions.get("pyaedt")
    if reference is None:
        dataset_v2.check_versions(aedt, pyaedt)
    else:
        dataset_v2.check_reference_compatible(dict(reference), aedt, pyaedt, REPO_ROOT)


def lattice_record(expected: Mapping) -> tuple[list[float], list[float], list[float], list[int]]:
    """A plan task's lattice record (``array_plan.predict_exit_lattice``'s) as its bounds, steps and node counts along
    x_e, y_e and z_e; ``ValueError`` for a record that is not a lattice."""
    try:
        lo, hi, step = ([float(v) for v in expected[name]] for name in ("lo_mm", "hi_mm", "step_mm"))
        counts = list(expected["counts"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"the plan's lattice {expected!r} is not a lattice record ({exc})") from None
    if not (len(lo) == len(hi) == len(step) == len(counts) == 3 and all(type(c) is int for c in counts)):
        raise ValueError(f"the plan's lattice {expected!r} does not give three bounds, steps and node counts")
    return lo, hi, step, counts


def _require_expected_lattice(face_geometry: Any, lattice: Any, expected: Mapping) -> None:
    """Plan mode (HF-057): the exit face as AEDT describes it, before the solve, must be a rectangle whose HF-045
    lattice is the one the plan predicted from the frozen geometry; otherwise nothing is solved."""
    problems = []
    shape = face_geometry.cross_section.get("shape")
    if shape != "rectangle":
        problems.append(f"the exit cross-section is {shape!r}, and the plan predicts a rectangle")
    if lattice is None:
        problems.append("the run exports no manual lattice")
    else:
        try:
            plan = lattice_record(expected)
        except ValueError as exc:
            problems.append(str(exc))
        else:   # the array_run check's exit_grid row compares the declared lattice with the plan's the same way
            problems += array_plan.lattice_problems("the exit face's lattice", lattice.lo_mm, lattice.hi_mm,
                                                    lattice.step_mm, lattice.counts, *plan)
    if problems:
        raise RuntimeError("the exit face's lattice is not the plan's prediction from the frozen geometry, so nothing "
                           "is solved (HF-057): " + "; ".join(problems))
    log.info("exit lattice %s is the plan's prediction", " x ".join(str(n) for n in lattice.counts))


def prepare_and_solve(cfg: RunConfig, dirs: JobDirs, job_id: str, hfss: Any, staged: StagedProject, *,
                      expected_lattice: Mapping | None = None) -> SolvedContext:
    """Set up the run design in the open job copy as the legacy script did, check what can be checked before the
    solve, solve, and verify that every requested polarization was solved. Returns what extraction needs.

    ``staged`` is ``stage_project``'s result: the copy that ``hfss`` opened, made and hashed before the session.
    ``expected_lattice`` (plan mode) is the plan's predicted exit lattice, which the face's must equal before the
    solve (HF-057)."""
    from bbsim.extract import ExtractionContext, describe_exit_face, exit_grid
    from bbsim.geometry import create_exit_coordinate_system, create_exit_face_list, select_faces
    from bbsim.hfss_setup import (
        EXIT_CS,
        EXIT_FACE_LIST,
        add_exit_flux_expression,
        add_outgoing_power_expression,
        aedt_messages,
        assign_radiation_boundaries,
        assign_wall_facets,
        clear_boundaries_and_excitations,
        create_plane_wave,
        create_setup_and_polarization_sweep,
        initialize_variables,
        insert_far_field_sphere,
        missing_polarizations,
        reset_initial_mesh_settings,
        save_project,
        set_total_fields,
        solve,
        solved_variations,
    )
    from bbsim.inventory import object_record
    from bbsim.session import aedt_versions

    job_project = staged.job_project
    freq = cfg.solver.frequency_ghz
    run_design = frequency_design_name(cfg.project.design, freq)
    exc = cfg.excitation
    phi_values = incident_angle_values(exc.phi_lower_deg, exc.phi_upper_deg, exc.phi_step_deg)
    theta_values = incident_angle_values(exc.theta_lower_deg, exc.theta_upper_deg, exc.theta_step_deg)

    versions = aedt_versions(hfss)
    log.info("versions: %s", versions)
    if cfg.project.design not in hfss.design_list:
        raise RuntimeError(f"design {cfg.project.design!r} not in {job_project}; designs: {list(hfss.design_list)}")
    if run_design in hfss.design_list:
        raise RuntimeError(f"design {run_design!r} already exists in the job copy")
    hfss.set_active_design(cfg.project.design)
    hfss.duplicate_design(run_design)
    if hfss.design_name != run_design:
        hfss.set_active_design(run_design)
    if hfss.modeler.model_units != "mm":
        raise RuntimeError(f"model units are {hfss.modeler.model_units!r}; the configuration assumes mm")
    # duplicate_design copies the base design's global mesh settings (manual 0.1 mm model resolution in the
    # GUI-made base designs); the legacy script solved a fresh design with Auto. Restore the reference settings.
    mesh_inherited, mesh_requested = reset_initial_mesh_settings(hfss)
    log.info("initial mesh settings reset to %s (inherited: %s)", mesh_requested, mesh_inherited)
    wall_facets = None
    if cfg.solver.wall_normal_deviation_deg is not None:
        wall_facets = assign_wall_facets(hfss, cfg.geometry.object, cfg.solver.wall_normal_deviation_deg)
        log.info("wall facets: %s", wall_facets)

    # Same order as legacy main(): variables, faces/boundaries, exit CS, calculator, setup, sphere, plane wave.
    clear_boundaries_and_excitations(hfss)
    entrance, exit_face = select_faces(hfss, cfg.geometry)
    if cfg.output.schema == 2:
        _require_v2_exit_normal(cfg, exit_face)   # before the solve; extract_and_write checks again
    solid = object_record(hfss, cfg.geometry.object)
    _require_vacuum(solid)
    log.info("entrance face %s, exit face %s, material %s", entrance, exit_face, solid["material"])
    initialize_variables(hfss, exc.ei_v_per_m)
    create_exit_face_list(hfss, exit_face, name=EXIT_FACE_LIST)
    assign_radiation_boundaries(hfss, entrance, exit_face)
    create_exit_coordinate_system(hfss, exit_face, cfg.geometry.exit_cs_x, cfg.geometry.exit_cs_y, name=EXIT_CS)
    # add_outgoing_power_expression clears every named expression, so the flux is added after it.
    add_outgoing_power_expression(hfss)
    add_exit_flux_expression(hfss)
    setup, parametric = create_setup_and_polarization_sweep(hfss, freq, cfg.solver.max_delta_e, cfg.solver.max_passes)
    sphere, grid = insert_far_field_sphere(hfss, cfg.far_field, freq)
    plane_wave = create_plane_wave(hfss, entrance, exc, freq)
    set_total_fields(hfss)

    # The export lattice needs only geometry: describe the exit face before the licensed solve, so a
    # geometry problem costs seconds rather than a solve.
    ctx = ExtractionContext(
        frequency_ghz=freq, solution=last_adaptive_solution(freq), sphere_name=sphere,
        phi_values=phi_values, theta_values=theta_values, grid=grid, exit_field=cfg.exit_field,
        exit_face=exit_face, exit_cs_name=EXIT_CS, exit_cs_x=cfg.geometry.exit_cs_x, exit_cs_y=cfg.geometry.exit_cs_y,
        scratch_dir=dirs.scratch, timeout_s=cfg.output.export_timeout_s,
        incoming_power_w=incoming_power_w(exc.ei_v_per_m, entrance.area_mm2),
    )
    face_geometry = describe_exit_face(hfss, ctx)
    lattice = exit_grid(face_geometry, cfg.exit_field) if cfg.exit_field.manual else None
    if lattice is not None:
        log.info("exit lattice: %d points per incident angle (%s)", lattice.lattice_points, face_geometry.bounds_method)
    if expected_lattice is not None:
        _require_expected_lattice(face_geometry, lattice, expected_lattice)
    if cfg.output.schema == 2:
        _require_v2_lattice(face_geometry.cross_section, lattice)
    # Checked (CDX-046): an AEDT that will not save stops the run before the licensed solve, not after it.
    save_project(hfss, "the setup")

    log.info("solving %s with %d cores (%d incident angles x 2 polarizations)", setup, cfg.solver.cores,
             len(phi_values) * len(theta_values))
    t_solve = time.perf_counter()
    if not solve(hfss, cfg.solver.cores):
        raise RuntimeError("hfss.analyze returned False")
    solve_s = time.perf_counter() - t_solve
    log.info("solve finished in %.1f s", solve_s)
    # analyze() returning True proves nothing (HPC step 4, baseline1: 66 s, no adaptive pass, no fields).
    # Ask AEDT what it solved, and keep its own messages: the solver's reason lives only there.
    messages = aedt_messages(hfss)
    for message in messages:
        # AEDT tags its entries [info], [warning], [error]; its errors are warnings here because the
        # runner decides below what is fatal.
        (log.warning if "[error]" in message or "[warning]" in message else log.info)("aedt: %s", message)
    variations = solved_variations(hfss, ctx.solution)
    log.info("solved variations of %s: %s", ctx.solution, variations)
    missing = missing_polarizations(variations, cfg.output.polarizations)
    if missing:
        raise RuntimeError(
            f"no solved variations for Ephi={missing} in {ctx.solution!r}: AEDT holds {variations}; "
            f"the solve ran for {solve_s:.1f} s without producing them. AEDT messages: " + " | ".join(messages)
        )
    # AEDT writes solution data into the project's results folder only when the project is saved; a run
    # that fails during extraction closes the project unsaved (HPC steps 4 and 4a lost their solves that way).
    # A save AEDT refuses (PyAEDT's False, or an exception) stops the run here, before any export, and before the
    # message that the solved project is kept (CDX-046).
    save_project(hfss, "the solve")
    return SolvedContext(
        ctx=ctx, face_geometry=face_geometry, lattice=lattice, entrance=entrance, exit_face=exit_face, solid=solid,
        setup=setup, parametric=parametric, sphere=sphere, grid=grid, plane_wave=plane_wave, wall_facets=wall_facets,
        mesh_inherited=mesh_inherited, mesh_requested=mesh_requested, variations=variations, messages=messages,
        solve_s=solve_s, source=staged.source, source_sha=staged.source_sha, job_project=job_project,
        job_sha=staged.job_sha, recorded_id=staged.recorded_id, record=staged.record, run_design=run_design,
        versions=versions,
    )


# BBRsim pairs the two polarizations' exit points by row and requires equal X, Y and Z to this (BBR012 at load).
POSITION_MATCH_M = 1e-12


def _require_shared_exit_points(waveguide0: Any, waveguide1: Any) -> None:
    """Row by row, the Ephi=0 and Ephi=1 tables must hold the same key and the same exit point to
    ``POSITION_MATCH_M`` in X, Y and Z; ``dataset_v2.build`` compares lattice nodes only. Names the first key."""
    keys0 = waveguide0[KEY_COLUMNS].to_numpy(dtype=float)
    keys1 = waveguide1[KEY_COLUMNS].to_numpy(dtype=float)
    if keys0.shape != keys1.shape or not np.array_equal(keys0, keys1):
        raise DatasetV2Error("the Ephi=0 and Ephi=1 waveguide tables do not list the same keys row by row: BBRsim "
                             "pairs the polarizations' exit points by row (BBR012)")
    xyz0 = waveguide0[["X", "Y", "Z"]].to_numpy(dtype=float)
    xyz1 = waveguide1[["X", "Y", "Z"]].to_numpy(dtype=float)
    diff = np.abs(xyz0 - xyz1)
    bad = ~(diff <= POSITION_MATCH_M)
    if bad.any():
        i = int(np.flatnonzero(bad.any(axis=1))[0])
        j = int(np.flatnonzero(bad[i])[0])
        phi, theta = keys0[i]
        row = i - int(np.flatnonzero(np.all(keys0 == keys0[i], axis=1))[0])
        raise DatasetV2Error(
            f"at (phi, theta) = ({phi:g}, {theta:g}) deg the Ephi=0 and Ephi=1 exit points differ in {'XYZ'[j]} by "
            f"{diff[i, j]!r} m at row {row} of the key (limit {POSITION_MATCH_M:g} m): BBRsim pairs them by row (BBR012)")


def _column_order_check(cfg: RunConfig, dirs: JobDirs, solved: SolvedContext,
                        extracted: Mapping[int, dataset_v2.Extracted], reference: Mapping | None, commit: str
                        ) -> tuple[dict, dict, dict | None]:
    """(block, reference, raw tables to keep as evidence or None). The reference job measures the reference on its
    own raw tables (HF-041); every other job validates the supplied one against its versions and exporter first."""
    aedt, pyaedt = solved.versions.get("aedt"), solved.versions.get("pyaedt")
    if reference is None:
        if not is_reference_job(cfg):
            raise DatasetV2Error(f"schema 2 needs the column-order reference for {cfg.project.dataset_id} at "
                                 f"{float(cfg.solver.frequency_ghz):g} GHz (--column-order-reference JSON; HF-041)")
        raw = {e: extracted[e].waveguide for e in (0, 1)}
        measured = dataset_v2.build_reference(raw, cfg.project.dataset_id, cfg.solver.frequency_ghz, str(dirs.root),
                                              commit, aedt, pyaedt, REPO_ROOT)
        return dataset_v2.order_check_block("reference", measured), measured, raw
    reference = dict(reference)
    dataset_v2.check_reference_compatible(reference, aedt, pyaedt, REPO_ROOT)
    return dataset_v2.order_check_block("shared_exporter", reference), reference, None


def _write_reference_evidence(scratch: Path, reference: Mapping, raw: Mapping[int, Any]) -> None:
    """The reference job keeps the two raw tables it hashed, byte for byte, and the reference object, so its
    ``raw_table_sha256`` can be checked against files on disk (HF-041)."""
    for e in (0, 1):
        path = scratch / f"column_order_reference_Ephi{e}.csv"
        path.write_bytes(raw[e].to_csv(index=False).encode())
        if sha256_file(path) != reference["raw_table_sha256"][str(e)]:
            raise RuntimeError(f"{path} does not hash to the reference's raw_table_sha256[{e!r}]")
    write_json(scratch / "column_order_reference.json", reference)


def extract_and_write(cfg: RunConfig, dirs: JobDirs, job_id: str, hfss: Any, solved: SolvedContext, t_start: float,
                      reference: dict | None = None, reexport: dict | None = None, *,
                      environment: EnvironmentState, session: OpenSession | None = None,
                      plan_record: Mapping | None = None) -> Path:
    """Extract every polarization into memory, check, then write the CSVs, the manifest and the sidecar.

    Order (HF-052): both polarizations are extracted and checked for finite values; under
    schema 2 the column-order reference is measured (the reference job) or validated (every other job), the 2.0
    tables and sidecar blocks are built (``dataset_v2.build``); then the pairing and lattice checks. Only then is
    anything written, so a failure leaves no CSV, manifest or sidecar. Schema 1 writes the 1.0 tables and sidecar
    exactly as before. ``reexport`` (schema 2 only) is the --from-solved provenance for the sidecar; ``t_start`` is
    the job's start on ``time.perf_counter``. ``environment`` is the code checkout's state, which the caller read
    before the session: this runs after the licensed solve, so it never reads git itself. ``session``'s records go
    into the manifest, and its ``before_write`` runs after the last AEDT call (HF-057). ``plan_record`` (plan mode)
    is the manifest's ``plan`` block {path, sha256, index}; a single run records none. Returns the manifest path.
    """
    from bbsim.convergence import qualify
    from bbsim.extract import RIM_BAND_REL, evaluate_exit_powers, extract_far_field, extract_waveguide
    from bbsim.hfss_setup import EXIT_CS, EXIT_FACE_LIST, export_convergence_text, save_project
    from bbsim.readers import require_finite_table

    schema = cfg.output.schema
    if schema == 1 and reexport is not None:
        raise ValueError("a re-export record applies to schema 2 only; this run writes schema 1")
    check_before_aedt(cfg, reference)   # pure; repeated here for callers that did not come through run_job
    if schema == 2:
        _require_v2_exit_face(cfg, solved)   # every path in, a --from-solved context included
    ctx = solved.ctx
    freq = cfg.solver.frequency_ghz
    dataset_id = cfg.project.dataset_id
    extracted: dict[int, dataset_v2.Extracted] = {}
    convergence: dict[str, str] = {}
    rim_log: dict[str, list] = {}

    with aedt_messages_on_error(hfss, set(solved.messages)):
        for ephi in cfg.output.polarizations:
            if schema == 2:
                abs_w, flux_w = evaluate_exit_powers(hfss, ctx, ephi)
                waveguide = extract_waveguide(hfss, ctx, ephi, solved.face_geometry, rim_log,
                                              exit_powers=(abs_w, flux_w))
            else:   # schema 1 evaluates |Re S| inside extract_waveguide and needs no flux
                abs_w, flux_w = {}, {}
                waveguide = extract_waveguide(hfss, ctx, ephi, solved.face_geometry, rim_log)
            far_field = extract_far_field(hfss, ctx, ephi)
            for name, table in (("waveguide.csv", waveguide), ("far_field.csv", far_field)):
                require_finite_table(table, f"{dataset_id} Ephi={ephi} {name}")
            extracted[ephi] = dataset_v2.Extracted(waveguide=waveguide, far_field=far_field, flux_w=flux_w, abs_w=abs_w)
            convergence[str(ephi)] = export_convergence_text(hfss, solved.setup, f"Ephi='{ephi}'",
                                                             dirs.logs / f"convergence_Ephi{ephi}.txt")
            log.info("extracted Ephi=%s (%d waveguide rows, %d far-field rows)", ephi, len(waveguide), len(far_field))
        save_project(hfss, "the extraction")      # checked (CDX-046): a refused save publishes nothing
    if session is not None:
        session.before_write()          # the run's AEDT calls are done: AEDT's temp directory goes back now (HF-057)

    # CDX-025: whether each polarization met the configured target. Recorded, not yet a gate (OR-027): the
    # job array's assembly will refuse unqualified datasets.
    qualification = {e: qualify(text, cfg.solver.max_delta_e) for e, text in convergence.items()}
    qualified = bool(qualification) and all(q["qualified"] for q in qualification.values())
    if not qualified:
        log.warning("convergence target %g not met or not shown (CDX-025; recorded in the manifest, not yet a "
                    "gate): %s", cfg.solver.max_delta_e,
                    "; ".join(f"Ephi={e}: {', '.join(q['reasons'])}" for e, q in qualification.items()
                              if not q["qualified"]))

    result = None
    evidence = None
    if schema == 2:
        order_check, reference, evidence = _column_order_check(cfg, dirs, solved, extracted, reference,
                                                               environment.commit)
        modes = modes_record_v2(solved.face_geometry.cross_section)
        if modes is None:
            raise DatasetV2Error(f"schema 2 has no modes for a {solved.face_geometry.cross_section['shape']!r} "
                                 "cross-section")
        result = dataset_v2.build(
            extracted, grid=solved.lattice, cross_section=solved.face_geometry.cross_section,
            incoming_w=ctx.incoming_power_w, phi_values=ctx.phi_values, theta_values=ctx.theta_values,
            ff_steps_deg=(solved.grid.theta_step_deg, solved.grid.phi_step_deg),
            polarization_filter_limit_ghz=modes["polarization_filter_limit_ghz"], order_check=order_check)
        _require_shared_exit_points(extracted[0].waveguide, extracted[1].waveguide)

    # Geant4 pairs the Ephi=0 and Ephi=1 exit points of each incident angle by row index.
    exit_point_counts = {e: ex.waveguide.groupby(KEY_COLUMNS).size().to_dict() for e, ex in extracted.items()}
    first_ephi, first_counts = next(iter(exit_point_counts.items()))
    for ephi, counts in exit_point_counts.items():
        if counts != first_counts:
            raise RuntimeError(
                f"exit-point counts per incident angle differ between polarization Ephi={first_ephi} and "
                f"Ephi={ephi}; Geant4 pairs them by row index. Ephi={first_ephi}: {first_counts}; Ephi={ephi}: {counts}"
            )
    retained = int(next(iter(first_counts.values())))
    lattice = solved.lattice
    lattice_points = lattice.lattice_points if lattice is not None else None
    if lattice_points is not None and retained > lattice_points:
        # HFSS can omit lattice points outside the solid but never add any: the lattice model is wrong,
        # and so would be the grid recorded in the manifest and the sidecar.
        raise RuntimeError(
            f"exit field has {retained} points per incident angle, more than the {lattice_points} points of "
            f"the ExportOnGrid lattice {lattice.counts}; the lattice model disagrees with the HFSS export"
        )
    outside_points = "omitted" if lattice_points is not None and retained < lattice_points else "none"

    # --- write phase: nothing above has written a CSV, manifest or sidecar --------------------------------------
    if evidence is not None:
        _write_reference_evidence(dirs.scratch, reference, evidence)
    outputs: dict[str, str] = {}
    files: dict[str, dict] = {}
    for ephi in cfg.output.polarizations:
        tables = ((result.waveguide[ephi], result.far_field[ephi]) if result is not None
                  else (extracted[ephi].waveguide, extracted[ephi].far_field))
        out_dir = dirs.dataset_dir(dataset_id, freq, ephi)
        out_dir.mkdir(parents=True, exist_ok=False)
        for name, table in zip(("waveguide.csv", "far_field.csv"), tables):
            path = out_dir / name
            table.to_csv(path, index=False)
            files[f"{out_dir.name}/{name}"] = {"sha256": sha256_file(path), "bytes": path.stat().st_size,
                                               "rows": len(table)}
        outputs[str(ephi)] = str(out_dir)
        log.info("wrote %s (%d waveguide rows, %d far-field rows)", out_dir, len(tables[0]), len(tables[1]))

    exc = cfg.excitation
    entrance, exit_face, face_geometry, grid = solved.entrance, solved.exit_face, solved.face_geometry, solved.grid
    exit_cs_z = [float(v) + 0.0 for v in np.cross(cfg.geometry.exit_cs_x, cfg.geometry.exit_cs_y)]
    solver = {"setup": solved.setup, "parametric_sweep": solved.parametric, "frequency_label": frequency_label(freq),
              "formulation": "TotalFields", "max_delta_e": cfg.solver.max_delta_e,
              "max_passes": cfg.solver.max_passes, "cores": cfg.solver.cores,
              "initial_mesh_settings": {"requested": solved.mesh_requested,
                                        "inherited_from_base_design": solved.mesh_inherited},
              "wall_facets": solved.wall_facets,
              "boundaries": {"rbin": "Radiation on entrance", "rbout": "Radiation on exit", "other faces": "default PEC"},
              "solved_variations": solved.variations, "convergence": convergence, "qualified": qualified,
              "convergence_qualification": {e: {k: q[k] for k in ("achieved", "target", "passes_completed",
                                                                  "max_passes", "converged", "qualified", "reasons")}
                                            for e, q in qualification.items()}}
    if reexport is not None:
        solver["reused_solution"] = reexport.get("source_job")
    exit_field = {"manual": cfg.exit_field.manual, "resolution_mm": list(cfg.exit_field.resolution_mm),
                  "step_mm": list(lattice.step_mm) if lattice is not None else None,
                  "points_in_si": True, "field_in_ref_cs": False, "coordinate_system": EXIT_CS,
                  "grid": lattice.axis_record() if lattice is not None else None, "lattice_points": lattice_points,
                  "retained_points_per_key": retained, "outside_points": outside_points,
                  "rim": {"rule": "every lattice point inside the disc kept; inside points without field set to zero; "
                                  "points outside dropped (ledger HF-029)",
                          "band_rel": RIM_BAND_REL, "per_key": rim_log} if rim_log else None}
    extra: dict[str, Any] = {}
    if result is not None:
        # Fresh pre-clamp powers (manifest-only provenance; reexport_check derives the dark flags from them), and the
        # diagnostics that are not sidecar fields.
        exit_field["raw_powers_per_key"] = result.raw_powers_per_key
        extra["transmittance_checks"] = {"flux_over_abs_re_s": result.flux_over_abs_re_s,
                                         "rho_per_key": result.rho_per_key}
    if plan_record is not None:
        extra["plan"] = {"path": str(plan_record["path"]), "sha256": plan_record["sha256"],
                         "index": plan_record["index"]}
    manifest = build_manifest(
        environment=environment,
        config=config_to_dict(cfg),
        job={"id": job_id, "root": str(dirs.root), "wall_time_s": time.perf_counter() - t_start,
             "solve_time_s": solved.solve_s},
        project={"source_path": str(solved.source), "source_sha256": solved.source_sha,
                 "job_copy": str(solved.job_project),
                 "dataset_id_record": {"path": str(solved.record) if solved.record else None,
                                       "dataset_id": solved.recorded_id},
                 "job_copy_sha256_before_run": solved.job_sha, "base_design": cfg.project.design,
                 "run_design": solved.run_design},
        geometry={"object": cfg.geometry.object, "units": "mm", "material": solved.solid["material"],
                  "bounding_box_mm": solved.solid["bounding_box_mm"],
                  "entrance_face": asdict(entrance), "exit_face": asdict(exit_face),
                  "exit_cs": {"name": EXIT_CS, "x": list(cfg.geometry.exit_cs_x), "y": list(cfg.geometry.exit_cs_y),
                              "z": exit_cs_z},
                  "pose": cfg.geometry.pose, "symmetry": list(cfg.geometry.symmetry),
                  "cross_section": face_geometry.cross_section, "bounds_method": face_geometry.bounds_method,
                  "exit_face_vertex_count": face_geometry.vertex_count,
                  "edge_sampling_sagitta_mm": face_geometry.sagitta_mm},
        excitation={"plane_wave": solved.plane_wave, "ei_v_per_m": exc.ei_v_per_m,
                    "incident_phi_deg": ctx.phi_values.tolist(), "incident_theta_deg": ctx.theta_values.tolist(),
                    "polarizations": list(cfg.output.polarizations),
                    "polarization_convention": "Ephi=0: E_theta=1; Ephi=1: E_phi=1",
                    "incoming_power_w": ctx.incoming_power_w},
        solver=solver,
        far_field={"sphere": solved.sphere, "theta_step_deg": grid.theta_step_deg, "phi_step_deg": grid.phi_step_deg,
                   "points_per_angle": grid.points_per_angle, "coordinate_system": EXIT_CS,
                   "radiation_surface": EXIT_FACE_LIST,
                   "theta_deg": _axis_summary(grid.theta_values, grid.theta_step_deg),
                   "phi_deg": _axis_summary(grid.phi_values, grid.phi_step_deg)},
        exit_field=exit_field,
        outputs=outputs,
        versions=solved.versions,
        aedt=session.aedt if session is not None else None,
        aedt_temp=session.aedt_temp if session is not None else None,
        **extra,
    )
    manifest_path = dirs.root / "manifest.json"
    write_json(manifest_path, manifest)
    for out_dir in outputs.values():
        shutil.copy2(manifest_path, Path(out_dir) / "manifest.json")
    # Built from the manifest exactly as written, so the sidecar cannot disagree with it. Never edited afterwards:
    # the record shares objects with the parsed manifest.
    written, manifest_sha = json.loads(manifest_path.read_text()), sha256_file(manifest_path)
    if result is None:
        record = build_dataset_record(written, manifest_sha, files)
    else:
        record = build_dataset_record_v2(written, manifest_sha, files, result.blocks, reexport)
    write_json(dirs.root / dataset_record_name(dataset_id, freq), record)
    log.info("wrote the schema-%d sidecar %s", schema, dataset_record_name(dataset_id, freq))
    return manifest_path


# --- --from-solved: a fresh schema-2 extraction on a copy of a solved job, no solve (HF-050) ------------------------

@dataclass(frozen=True)
class SourceJob:
    """A finished job to re-export: its directory, its manifest exactly as read (never modified), the run
    configuration built from that manifest's ``config`` block with ``output.schema`` set to 2, its solved project
    and ``.aedtresults`` folder, and the manifest's sha256."""

    job_dir: Path
    manifest: dict
    config: RunConfig
    project: Path
    results: Path
    manifest_sha256: str


def _manifest_field(manifest: Mapping, section: str, key: str) -> Any:
    block = manifest.get(section)
    if not isinstance(block, Mapping) or key not in block:
        raise ValueError(f"the source manifest has no {section}.{key}")
    return block[key]


def _check_source_records(cfg: RunConfig, manifest: Mapping) -> None:
    """The re-export asks AEDT for the source's run design, solution, far-field sphere and incidence keys by the names
    and grids this code derives from the source's config; refuse a source that recorded different ones, before any
    copy, rather than halfway through the extraction."""
    from bbsim.hfss_setup import far_field_sphere
    from bbsim.naming import setup_name

    freq = cfg.solver.frequency_ghz
    exc = cfg.excitation
    sphere, grid = far_field_sphere(cfg.far_field, freq)
    derived = {
        ("project", "run_design"): frequency_design_name(cfg.project.design, freq),
        ("solver", "setup"): setup_name(freq),
        ("far_field", "sphere"): sphere,
        ("far_field", "theta_step_deg"): grid.theta_step_deg,
        ("far_field", "phi_step_deg"): grid.phi_step_deg,
        ("far_field", "points_per_angle"): grid.points_per_angle,
        ("excitation", "incident_phi_deg"): incident_angle_values(exc.phi_lower_deg, exc.phi_upper_deg,
                                                                  exc.phi_step_deg).tolist(),
        ("excitation", "incident_theta_deg"): incident_angle_values(exc.theta_lower_deg, exc.theta_upper_deg,
                                                                    exc.theta_step_deg).tolist(),
    }

    def same(recorded: Any, value: Any) -> bool:
        if isinstance(value, str) or isinstance(recorded, str):
            return recorded == value
        try:
            return bool(np.shape(recorded) == np.shape(value) and np.allclose(recorded, value, rtol=1e-12, atol=1e-9))
        except (TypeError, ValueError):
            return False

    differ = [f"{section}.{key}: recorded {_manifest_field(manifest, section, key)!r}, derived {value!r}"
              for (section, key), value in derived.items() if not same(_manifest_field(manifest, section, key), value)]
    if differ:
        raise ValueError("the source job's manifest records what this code would not reproduce from its config: "
                         + "; ".join(differ))


def load_source_job(job_dir: Path) -> SourceJob:
    """Read and check a finished job for --from-solved, without changing it: its ``manifest.json``, its solved
    project (``project/<name>.aedt``, the manifest's ``job_copy``) and that project's ``.aedtresults`` folder must
    exist. The configuration is a deep copy of the manifest's ``config`` block with
    ``output.schema`` = 2, whether the source recorded 1 or nothing (HF-050); every other setting is the source's.
    The job must sit in its dataset folder, where the fresh job goes."""
    job_dir = Path(os.path.abspath(job_dir))
    manifest_path = job_dir / "manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"no manifest.json in {job_dir}: --from-solved needs a finished job")
    manifest = json.loads(manifest_path.read_text())
    if not isinstance(manifest, dict):
        raise ValueError(f"{manifest_path} must hold a manifest object, got {type(manifest).__name__}")
    project = job_dir / "project" / Path(_manifest_field(manifest, "project", "job_copy")).name
    if not project.is_file():
        raise FileNotFoundError(f"the source job's solved project {project} does not exist")
    results = project.with_suffix(".aedtresults")
    if not results.is_dir():
        raise FileNotFoundError(f"the source job's solution folder {results} (.aedtresults) does not exist")
    _manifest_field(manifest, "project", "run_design")
    raw = manifest.get("config")
    if not isinstance(raw, Mapping):
        raise ValueError(f"{manifest_path} has no config block")
    raw = copy.deepcopy(dict(raw))
    output = raw.get("output")
    raw["output"] = {**(output if isinstance(output, Mapping) else {}), "schema": 2}
    # A re-export is judged against its source, whatever the source recorded (a source made before CDX-039 records
    # nothing); the fresh manifest records it, equal to the sidecar's provenance.reexport.source_job.
    raw["acceptance"] = {"kind": "reexport_check", "reference_job": str(job_dir)}
    cfg = build_config(raw)
    if cfg.solver.frequency_ghz is None:
        # A config may leave it to a plan (HF-057, OR-047), but a finished job always ran at one and records it.
        raise ValueError(f"{manifest_path}'s config records no solver.frequency_ghz: a finished job records its own")
    folder = f"{cfg.project.dataset_id}_{frequency_label(cfg.solver.frequency_ghz)}"
    if job_dir.parent.name != folder:
        raise ValueError(f"the source job {job_dir} is not in its dataset folder {folder}/, where the fresh job goes")
    record = manifest["project"].get("dataset_id_record") or {}
    if record.get("dataset_id") is not None and record["dataset_id"] != cfg.project.dataset_id:
        raise ValueError(f"the source manifest's config dataset_id {cfg.project.dataset_id!r} disagrees with its "
                         f"dataset_id record {record!r}")
    _check_source_records(cfg, manifest)
    return SourceJob(job_dir=job_dir, manifest=manifest, config=cfg, project=project, results=results,
                     manifest_sha256=sha256_file(manifest_path))


def reexport_job_dirs(source: SourceJob, job_id: str) -> JobDirs:
    """The fresh job's directories, in the source's dataset folder (``<dataset>_<label>/job_<job_id>``); an existing
    job directory, the source's own included, is refused."""
    return create_job_dirs(source.job_dir.parent.parent, source.config.project.dataset_id,
                           source.config.solver.frequency_ghz, job_id)


def check_before_copy(source: SourceJob, reference: Mapping | None) -> None:
    """--from-solved: refuse, before the source's project and ``.aedtresults`` (about 2.5 GB) are copied or AEDT
    starts, what the configuration and this checkout can refuse: ``check_before_aedt``; a source whose manifest
    records a gap material that is not vacuum (CDX-029); and a supplied column-order reference measured with other
    exporter files or another PyAEDT (HF-041). The AEDT-side checks after the copy is opened stay the authority:
    the AEDT version and the solid's material as AEDT reports it."""
    from bbsim.session import pyaedt_version

    check_before_aedt(source.config, reference)
    geometry = source.manifest.get("geometry")
    if isinstance(geometry, Mapping) and "material" in geometry and not is_vacuum(geometry["material"]):
        raise ValueError(f"the source job's manifest records gap {vacuum_refusal(geometry['material'])}")
    if reference is not None:
        try:
            installed = pyaedt_version()
        except ImportError as exc:
            raise DatasetV2Error(f"PyAEDT cannot be imported ({exc}), so the column-order reference cannot be "
                                 "checked against it") from exc
        dataset_v2.check_reference_before_aedt(dict(reference), installed, REPO_ROOT)


def stage_solved_job(source: SourceJob, dirs: JobDirs) -> StagedProject:
    """Copy the source's solved project and its ``.aedtresults`` into the fresh job, before any AEDT session, as
    step 8b's ``diagnose_poynting.py`` (archived in bbridge) copied a solved job; AEDT opens only the copy, so the
    source is never written."""
    job_project = dirs.project / source.project.name
    shutil.copy2(source.project, job_project)
    shutil.copytree(source.results, dirs.project / source.results.name)
    source_sha, job_sha = sha256_file(source.project), sha256_file(job_project)
    if job_sha != source_sha:
        raise RuntimeError(f"the copy {job_project} does not hash to its source {source.project}")
    record = source.manifest["project"].get("dataset_id_record") or {}
    return StagedProject(source=source.project, source_sha=source_sha, job_project=job_project, job_sha=job_sha,
                         recorded_id=record.get("dataset_id"),
                         record=Path(record["path"]) if record.get("path") else None)


def prepare_from_solved(source: SourceJob, dirs: JobDirs, job_id: str, hfss: Any, staged: StagedProject
                        ) -> SolvedContext:
    """Open the solved run design in the job's copy and rebuild the extraction context from its geometry, as
    ``prepare_and_solve`` builds it, without setting up or solving anything. Both polarizations must be solved
    already. ``staged`` is ``stage_solved_job``'s copy, which ``hfss`` opened. ``solve_s`` is None; the solve's mesh
    and facet settings are the source manifest's."""
    from bbsim.extract import ExtractionContext, describe_exit_face, exit_grid
    from bbsim.geometry import select_faces
    from bbsim.hfss_setup import (
        EXIT_CS,
        add_exit_flux_expression,
        add_outgoing_power_expression,
        aedt_messages,
        far_field_sphere,
        missing_polarizations,
        set_total_fields,
        solved_variations,
    )
    from bbsim.inventory import object_record
    from bbsim.naming import parametric_setup_name, plane_wave_name, setup_name
    from bbsim.session import aedt_versions

    cfg = source.config
    freq = cfg.solver.frequency_ghz
    run_design = source.manifest["project"]["run_design"]
    exc = cfg.excitation
    phi_values = incident_angle_values(exc.phi_lower_deg, exc.phi_upper_deg, exc.phi_step_deg)
    theta_values = incident_angle_values(exc.theta_lower_deg, exc.theta_upper_deg, exc.theta_step_deg)

    versions = aedt_versions(hfss)
    log.info("versions: %s", versions)
    if run_design not in hfss.design_list:
        raise RuntimeError(f"design {run_design!r} not in {staged.job_project}; designs: {list(hfss.design_list)}")
    if hfss.design_name != run_design:
        hfss.set_active_design(run_design)
    if hfss.modeler.model_units != "mm":
        raise RuntimeError(f"model units are {hfss.modeler.model_units!r}; the configuration assumes mm")
    solution = last_adaptive_solution(freq)
    variations = solved_variations(hfss, solution)
    log.info("solved variations of %s: %s", solution, variations)
    missing = missing_polarizations(variations, cfg.output.polarizations)
    if missing:
        raise RuntimeError(f"no solved variations for Ephi={missing} in {solution!r} of {staged.job_project} (copied "
                           f"from {source.job_dir}): AEDT holds {variations}, and a re-export never solves")

    entrance, exit_face = select_faces(hfss, cfg.geometry)
    _require_v2_exit_normal(cfg, exit_face)
    solid = object_record(hfss, cfg.geometry.object)
    _require_vacuum(solid)
    log.info("entrance face %s, exit face %s, material %s", entrance, exit_face, solid["material"])
    # The field type is a post-processing setting: HPC step 8b set it on copies of these solved projects and then
    # evaluated their fields. The source run set the same before its solve.
    set_total_fields(hfss)
    # HF-052: both named expressions are re-defined exactly as a fresh run defines them, the clear in
    # add_outgoing_power_expression first, not exit_flux alone: every
    # solve since the schema-2 runner (b2af196) saves exit_flux in its run design, and what AEDT does when a name
    # is added twice is unverified. Named expressions are calculator definitions; they do not touch the solution.
    add_outgoing_power_expression(hfss)
    add_exit_flux_expression(hfss)
    sphere, grid = far_field_sphere(cfg.far_field, freq)   # the solved design holds this sphere (load_source_job)

    ctx = ExtractionContext(
        frequency_ghz=freq, solution=solution, sphere_name=sphere,
        phi_values=phi_values, theta_values=theta_values, grid=grid, exit_field=cfg.exit_field,
        exit_face=exit_face, exit_cs_name=EXIT_CS, exit_cs_x=cfg.geometry.exit_cs_x, exit_cs_y=cfg.geometry.exit_cs_y,
        scratch_dir=dirs.scratch, timeout_s=cfg.output.export_timeout_s,
        incoming_power_w=incoming_power_w(exc.ei_v_per_m, entrance.area_mm2),
    )
    face_geometry = describe_exit_face(hfss, ctx)
    lattice = exit_grid(face_geometry, cfg.exit_field) if cfg.exit_field.manual else None
    _require_v2_lattice(face_geometry.cross_section, lattice)
    log.info("exit lattice: %d points per incident angle (%s)", lattice.lattice_points, face_geometry.bounds_method)
    solver = source.manifest.get("solver") or {}
    mesh = solver.get("initial_mesh_settings") or {}
    return SolvedContext(
        ctx=ctx, face_geometry=face_geometry, lattice=lattice, entrance=entrance, exit_face=exit_face, solid=solid,
        setup=setup_name(freq), parametric=parametric_setup_name(freq), sphere=sphere, grid=grid,
        plane_wave=plane_wave_name(freq), wall_facets=solver.get("wall_facets"),
        mesh_inherited=mesh.get("inherited_from_base_design"), mesh_requested=mesh.get("requested"),
        variations=variations, messages=aedt_messages(hfss), solve_s=None, source=staged.source,
        source_sha=staged.source_sha, job_project=staged.job_project, job_sha=staged.job_sha,
        recorded_id=staged.recorded_id, record=staged.record, run_design=run_design, versions=versions,
    )


def run_from_solved(source: SourceJob, dirs: JobDirs, job_id: str, reference: dict | None = None, *,
                    aedt_temp: tempdir.TempSettings | None = None, launch_timeout_s: float = DEFAULT_LAUNCH_TIMEOUT_S,
                    environment: EnvironmentState | None = None) -> Path:
    """--from-solved: copy the source's solved project and results, open the copy, extract and write schema 2 with
    no solve. The convergence tables are exported fresh from the opened design and qualified against the source's
    target (CDX-025); the sidecar records ``provenance.reexport`` and the manifest ``solver.reused_solution``.
    ``aedt_temp``, ``launch_timeout_s`` and ``environment`` as for ``run_job`` (HF-057)."""
    from bbsim.session import aedt_versions

    t_start = time.perf_counter()
    cfg = source.config
    check_before_copy(source, reference)
    environment = environment_state(REPO_ROOT) if environment is None else environment
    temp = tempdir.TempSettings() if aedt_temp is None else aedt_temp
    staged = stage_solved_job(source, dirs)
    reexport = {"source_job": str(source.job_dir), "source_manifest_sha256": source.manifest_sha256,
                "tool": REEXPORT_TOOL, "git_commit": environment.commit}
    with open_session(staged.job_project, source.manifest["project"]["run_design"], cfg.project.aedt_version, dirs,
                      temp, launch_timeout_s) as session:
        _check_versions_for_reference(reference, aedt_versions(session.hfss))
        solved = prepare_from_solved(source, dirs, job_id, session.hfss, staged)
        return extract_and_write(cfg, dirs, job_id, session.hfss, solved, t_start, reference, reexport,
                                 environment=environment, session=session)


def check_guard_before_job() -> dict:
    """The licence guard this run starts AEDT under (OR-045, HF-057), checked before any job directory exists:
    ``BBSIM_LICENCE_GUARD`` naming the shared guard and, when set, ``BBSIM_GUARD_RESERVATION`` naming this job's own
    reservation (``licence_guard.check_job_guard``). Raises ``licence_guard.GuardError``; the session checks again
    when it takes its launch lease."""
    return licence_guard.check_job_guard()


def _guard_log_line(found: Mapping[str, Any]) -> str:
    reservation = found.get("reservation") or "none (the session takes a direct one-slot reservation)"
    return f"licence guard {found.get('root')}, reservation {reservation}"


_FROM_SOLVED_EXCLUDES = (("--config", "config"), ("--frequency-ghz", "frequency_ghz"), ("--cores", "cores"),
                         ("--output-root", "output_root"), ("--project", "project"), ("--polarizations", "polarizations"),
                         ("--aedt-version", "aedt_version"))


def _environment_log_line(environment: EnvironmentState) -> str:
    dirty = " with uncommitted changes to tracked files" if environment.dirty else ""
    return f"code {REPO_ROOT} at {environment.commit}{dirty}; PyAEDT {environment.pyaedt_version}"


def _run_logged(what: str, run: Callable[[], Path]) -> int:
    """Run the job and log how it ended: 0 with the manifest's path last, 1 on an error, and ``SystemExit`` (SIGTERM's
    143) logged and raised again once the session's own cleanup has run."""
    try:
        manifest = run()
    except SystemExit as exc:
        _log_system_exit(what, exc)
        raise
    except Exception:  # noqa: BLE001 - report and exit non-zero; the session is already released
        log.exception("%s failed", what)
        return 1
    log.info("done; manifest at %s", manifest)
    return 0


def _main_from_solved(args: argparse.Namespace, cleanup: ExitStack) -> int:
    given = [flag for flag, attr in _FROM_SOLVED_EXCLUDES if getattr(args, attr) is not None]
    if given:
        print(f"refused: --from-solved takes its configuration from the source job's manifest; drop {', '.join(given)}",
              file=sys.stderr)
        return 2
    if not args.job_id:
        print("refused: --from-solved needs --job-id, the fresh job beside the source (for example tight1_v2)",
              file=sys.stderr)
        return 2
    reference = None
    try:
        source = load_source_job(args.from_solved)
        if args.column_order_reference is not None:
            reference = load_reference(args.column_order_reference)
        check_before_copy(source, reference)
        guard = check_guard_before_job()
        environment = environment_state(REPO_ROOT)
        dirs = reexport_job_dirs(source, args.job_id)
    except (OSError, ValueError, ImportError, licence_guard.GuardError, EnvironmentStateError) as exc:
        print(f"refused before the job started: {exc}", file=sys.stderr)
        return 2
    cleanup.callback(_setup_logging(dirs.logs / "run.log", args.log_level))
    write_json(dirs.root / "effective_config.json", config_to_dict(source.config))
    log.info("job %s -> %s: schema-2 re-export of %s, no solve", args.job_id, dirs.root, source.job_dir)
    log.info(_guard_log_line(guard))
    (log.warning if environment.dirty else log.info)(_environment_log_line(environment))
    if reference is not None:
        log.info("column-order reference from %s (%s, %s)", args.column_order_reference, reference.get("source_job"),
                 reference.get("aedt_version"))
    return _run_logged("re-export", functools.partial(
        run_from_solved, source, dirs, args.job_id, reference, aedt_temp=args.temp_settings,
        launch_timeout_s=args.launch_timeout_s, environment=environment))


# --- plan mode: one task of a frozen job-array plan (HF-057, OR-037, OR-045, HF-058) ---------------------------------

class PlanRefused(ValueError):
    """A plan-mode task refused by its pre-flight, before any job directory or AEDT exists (exit 2)."""


@dataclass(frozen=True)
class TaskContext:
    """What ``plan_preflight`` verified and the run uses: the plan (parsed from the bytes whose sha256 was checked),
    the task's index and record, its run configuration (``array_plan.task_config`` of the captured config bytes), the
    column-order reference parsed from its checked bytes (None for a probe), the job id and directory, the plan's
    predicted exit lattice, the policy, the code checkout's state, the temp-directory settings and launch timeout, the
    guard and reservation (``licence_guard.check_job_guard``) and the warnings to log once the job's log exists."""

    plan: dict
    index: int
    cfg: RunConfig
    reference: dict | None
    job_id: str
    expected_lattice: dict
    policy: str
    state: EnvironmentState
    plan_path: Path
    plan_sha256: str
    task: dict
    aedt_temp: tempdir.TempSettings
    launch_timeout_s: float
    guard: dict
    job_root: Path
    array_job_id: str | None = None
    warnings: tuple[str, ...] = ()

    @property
    def plan_record(self) -> dict:
        """The plan block of the manifest and ``start.json``: {path, sha256, index}."""
        return {"path": str(self.plan_path), "sha256": self.plan_sha256, "index": self.index}


def _plan_entry(plan: Mapping, name: str) -> Mapping:
    entry = plan.get(name)
    if not isinstance(entry, Mapping):
        raise PlanRefused(f"the plan's {name} must be an object, got {entry!r}")
    return entry


def _absolute(value: object, what: str) -> Path:
    if not isinstance(value, str) or not os.path.isabs(value):
        raise PlanRefused(f"the plan's {what} must be an absolute path, as every path a plan records is, got {value!r}")
    return Path(value)


def _frozen_bytes(entry: Mapping, what: str) -> tuple[Path, bytes]:
    """A plan input copy's bytes, read once, which must have the sha256 the plan records."""
    copy = _absolute(entry.get("copy"), f"{what}.copy")
    try:
        data = copy.read_bytes()
    except OSError as exc:
        raise PlanRefused(f"the plan's {what} copy {copy} cannot be read ({exc})") from None
    digest = hashlib.sha256(data).hexdigest()
    if digest != entry.get("sha256"):
        raise PlanRefused(f"the plan's {what} copy {copy} has sha256 {digest}, the plan records "
                          f"{entry.get('sha256')!r}: it changed after the plan was made")
    return copy, data


def _resolve_index(index: int | None, env: Mapping[str, str]) -> int:
    value = env.get("SLURM_ARRAY_TASK_ID")
    if value:
        if not _DIGITS.fullmatch(value):
            raise PlanRefused(f"SLURM_ARRAY_TASK_ID={value!r} is not a task index")
        if index is not None and index != int(value):
            raise PlanRefused(f"--plan-index {index} disagrees with SLURM_ARRAY_TASK_ID={value}")
        return int(value)
    if index is None:
        raise PlanRefused("no task index: give --plan-index N, or run as a task of the plan's array "
                          "(SLURM_ARRAY_TASK_ID)")
    return index


def _check_code(plan: Mapping) -> EnvironmentState:
    """The plan's code: the checkout this module was imported from, at the plan's commit with no change to a tracked
    file, and the plan's PyAEDT (``manifest.environment_state``, whose git error or timeout refuses)."""
    code = _plan_entry(plan, "code")
    repo = _absolute(code.get("repo"), "code.repo")
    try:
        same = os.path.samefile(repo, REPO_ROOT)
    except OSError:
        same = False
    if not same:
        raise PlanRefused(f"the plan's code {repo} is not the checkout this runner was imported from ({REPO_ROOT}): "
                          "the task must run the plan's own worktree")
    state = environment_state(repo)
    if state.commit != code.get("commit"):
        raise PlanRefused(f"the code {repo} is at {state.commit}, not the plan's commit {code.get('commit')!r}: its "
                          "HEAD moved after the plan was made")
    if state.dirty:
        raise PlanRefused(f"the code {repo} has uncommitted changes to tracked files; the task runs the plan's commit")
    if state.pyaedt_version != code.get("pyaedt_version"):
        raise PlanRefused(f"PyAEDT {state.pyaedt_version} is installed, and the plan was made with "
                          f"{code.get('pyaedt_version')!r}")
    return state


def _check_record(plan: Mapping, cfg: RunConfig) -> None:
    """The mapping record frozen beside the project copy: the plan's bytes, the record ``stage_project`` reads, and it
    confirms the plan's dataset for the config's design and object. A single run only warns without one (HF-057)."""
    from bbsim.mapping import recorded_dataset_id

    project = _plan_entry(plan, "project")
    entry = project.get("record")
    if not isinstance(entry, Mapping):
        raise PlanRefused(f"the plan's project records no frozen mapping record {{source, copy, sha256}}, got "
                          f"{entry!r}")
    copy, _ = _frozen_bytes(entry, "project.record")
    try:
        recorded, found = recorded_dataset_id(cfg.project.path, cfg.project.design, cfg.geometry.object)
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
        raise PlanRefused(f"the record beside the project copy {cfg.project.path} cannot be read "
                          f"({type(exc).__name__}: {exc})") from None
    if found is None or Path(found).resolve() != copy.resolve():
        raise PlanRefused(f"the record beside the project copy {cfg.project.path} is {found}, not the plan's frozen "
                          f"record {copy}")
    if recorded != plan.get("dataset_id"):
        raise PlanRefused(f"the frozen record {copy} gives {recorded!r} for {cfg.project.design}/{cfg.geometry.object}, "
                          f"not the plan's dataset {plan.get('dataset_id')!r}")


def _check_policy(plan: Mapping, policy: str, cfg: RunConfig, reference: Mapping | None) -> None:
    if cfg.acceptance.kind != "array_run":
        raise PlanRefused(f"the task config's acceptance kind is {cfg.acceptance.kind!r}, not 'array_run'")
    if cfg.project.dataset_id != plan.get("dataset_id"):
        raise PlanRefused(f"the task config's dataset_id {cfg.project.dataset_id!r} is not the plan's "
                          f"{plan.get('dataset_id')!r}")
    pols = sorted(cfg.output.polarizations)
    if policy == "production":
        if cfg.output.schema != 2 or pols != [0, 1] or reference is None:
            raise PlanRefused(f"a production plan runs schema 2, both polarizations, with the HF-041 column-order "
                              f"reference; this task's config writes schema {cfg.output.schema}, polarizations {pols}, "
                              f"and the plan records {'a' if reference is not None else 'no'} reference")
        return
    exc = cfg.excitation
    keys = (len(incident_angle_values(exc.theta_lower_deg, exc.theta_upper_deg, exc.theta_step_deg))
            * len(incident_angle_values(exc.phi_lower_deg, exc.phi_upper_deg, exc.phi_step_deg)))
    if cfg.output.schema != 1 or pols != [0] or keys != 1:
        raise PlanRefused(f"a probe plan runs the schema-1 single-angle diagnostic, one key and Ephi=0 exported "
                          f"(HF-058); this task's config writes schema {cfg.output.schema}, polarizations {pols}, "
                          f"{keys} keys")


def _check_lattice(task: Mapping, cfg: RunConfig, index: int) -> dict:
    """The task's effective step along l, at most c/(4f) (OR-048), and its recorded lattice, the frozen rectangle's
    HF-045 lattice at that step (``array_plan.predict_exit_lattice``, the planner's own). Returns the recorded one."""
    frequency = task.get("frequency_ghz")
    effective = task.get("effective_l_resolution_mm")
    try:
        limit = frequency_grid.max_l_step_mm(frequency)
        if isinstance(effective, bool) or not isinstance(effective, (int, float)) or not effective > 0:
            raise ValueError(f"effective_l_resolution_mm {effective!r} is not a positive step")
        if effective > limit:
            raise ValueError(f"its effective step along l {effective!r} mm exceeds c/(4f) = {limit:.7g} mm at "
                             f"{frequency!r} GHz (OR-048)")
        predicted = json.loads(json.dumps(array_plan.predict_exit_lattice(cfg, effective)))
    except (frequency_grid.GridError, ValueError) as exc:
        raise PlanRefused(f"task {index}: {exc}") from None
    recorded = task.get("lattice")
    if recorded != predicted:
        try:
            detail = array_plan.lattice_problems("the frozen rectangle's HF-045 lattice", predicted["lo_mm"],
                                                 predicted["hi_mm"], predicted["step_mm"], predicted["counts"],
                                                 *lattice_record(recorded))
        except ValueError as exc:
            detail = [str(exc)]
        raise PlanRefused(f"task {index}'s recorded lattice {recorded!r} is not the frozen rectangle's HF-045 lattice "
                          f"{predicted!r}" + (f": {'; '.join(detail)}" if detail else ""))
    if task.get("l_step_mm") != predicted["step_mm"][1] or predicted["step_mm"][1] > limit:
        raise PlanRefused(f"task {index}'s l_step_mm {task.get('l_step_mm')!r} is not its lattice's fitted step along "
                          f"l {predicted['step_mm'][1]!r} mm, at most c/(4f) = {limit:.7g} mm (OR-048)")
    return recorded


def _check_guard(plan: Mapping, env: Mapping[str, str]) -> tuple[dict, str | None]:
    """The plan's guard is the shared one every launch lease takes (OR-045), named by ``BBSIM_LICENCE_GUARD``; an
    array task runs under its own ``array`` reservation, bound to its array job (``licence_guard.check_job_guard``)."""
    guard_dir = _absolute(plan.get("guard_dir"), "guard_dir")
    try:
        here, shared = guard_dir.resolve(strict=True), Path(licence_guard.DEFAULT_ROOT).resolve(strict=True)
    except (OSError, RuntimeError) as exc:
        raise PlanRefused(f"the plan's guard {guard_dir} or the shared guard {licence_guard.DEFAULT_ROOT} cannot be "
                          f"reached ({exc})") from None
    if here != shared:
        raise PlanRefused(f"the plan's guard_dir {guard_dir} is not the shared guard {shared}, whose launch lease every "
                          "AEDT session takes (OR-045)")
    array_job = env.get("SLURM_ARRAY_JOB_ID") or None
    if array_job is not None and not _DIGITS.fullmatch(array_job):
        raise PlanRefused(f"SLURM_ARRAY_JOB_ID={array_job!r} is not a job id")
    found = licence_guard.check_job_guard(kinds=("array",) if array_job else None, env=env, root=guard_dir)
    return found, array_job


def plan_preflight(plan_path: Path, plan_sha256: str, index: int | None, env: Mapping[str, str]) -> TaskContext:
    """Every check of a job-array task before a job directory or AEDT exists (the module docstring's list), reading
    the plan's bytes and each frozen input's once and the job's variables from ``env`` only. Raises ``PlanRefused``,
    ``licence_guard.GuardError``, ``EnvironmentStateError``, ``ImportError`` (no PyAEDT) or another ``ValueError``
    (a config the overrides cannot build), each a refusal (exit 2)."""
    plan_path = Path(os.path.abspath(plan_path))
    try:
        data = plan_path.read_bytes()
    except OSError as exc:
        raise PlanRefused(f"the plan {plan_path} cannot be read ({exc})") from None
    digest = hashlib.sha256(data).hexdigest()
    if digest != plan_sha256:
        raise PlanRefused(f"the plan {plan_path} has sha256 {digest}, not the {plan_sha256} this task was given")
    try:
        plan = json.loads(data.decode("utf-8"))
    except (UnicodeDecodeError, ValueError) as exc:
        raise PlanRefused(f"the plan {plan_path} is not JSON ({exc})") from None
    if not isinstance(plan, dict) or plan.get("schema") != array_plan.PLAN_SCHEMA:
        raise PlanRefused(f"the plan {plan_path} is not a {array_plan.PLAN_SCHEMA} plan")
    index = _resolve_index(index, env)
    tasks = plan.get("tasks") if isinstance(plan.get("tasks"), list) else []
    found = [t for t in tasks if isinstance(t, dict) and type(t.get("index")) is int and t["index"] == index]
    if len(found) != 1:
        raise PlanRefused(f"the plan holds {len(found)} records of task {index}, not 1")
    task = found[0]
    policy = plan.get("policy")
    if policy not in ("production", "probe") or plan.get("probe") is not (policy == "probe"):
        raise PlanRefused(f"the plan's policy {policy!r} and probe flag {plan.get('probe')!r} are not production or "
                          "probe")
    if plan.get("acceptance_kind") != "array_run":
        raise PlanRefused(f"the plan's acceptance_kind is {plan.get('acceptance_kind')!r}, not 'array_run'")

    state = _check_code(plan)
    _, config_bytes = _frozen_bytes(_plan_entry(plan, "config"), "config")
    project = _plan_entry(plan, "project")
    project_copy = _absolute(project.get("copy"), "project.copy")
    try:
        project_sha = sha256_file(project_copy)
    except OSError as exc:
        raise PlanRefused(f"the plan's project copy {project_copy} cannot be read ({exc})") from None
    if project_sha != project.get("sha256"):
        raise PlanRefused(f"the plan's project copy {project_copy} has sha256 {project_sha}, the plan records "
                          f"{project.get('sha256')!r}: it changed after the plan was made")
    reference = None
    entry = plan.get("column_order_reference")
    if policy == "probe" and entry is not None:
        raise PlanRefused("a probe plan records no column-order reference (its schema-1 output neither needs nor "
                          f"applies one), and this one records {entry!r}")
    if policy == "production":
        if not isinstance(entry, Mapping):
            raise PlanRefused("a production plan records the HF-041 column-order reference its schema-2 runs embed, "
                              "and this one records none")
        _, reference_bytes = _frozen_bytes(entry, "column_order_reference")
        try:
            reference = json.loads(reference_bytes.decode("utf-8"))
        except (UnicodeDecodeError, ValueError) as exc:
            raise PlanRefused(f"the plan's column-order reference copy is not JSON ({exc})") from None
        if not isinstance(reference, dict):
            raise PlanRefused("the plan's column-order reference copy is not a JSON object")
        dataset_v2.check_reference_before_aedt(reference, state.pyaedt_version, REPO_ROOT)
    try:
        raw = tomllib.loads(config_bytes.decode("utf-8"))
    except (UnicodeDecodeError, ValueError) as exc:
        raise PlanRefused(f"the plan's config copy is not TOML ({exc})") from None
    cfg = array_plan.task_config(plan, index, raw_config=raw)
    _check_record(plan, cfg)
    _check_policy(plan, policy, cfg, reference)
    check_before_aedt(cfg, reference)
    expected_lattice = _check_lattice(task, cfg, index)
    try:
        temp = tempdir.TempSettings(**task["aedt_temp"])
    except (KeyError, TypeError, ValueError) as exc:
        raise PlanRefused(f"task {index}'s aedt_temp {task.get('aedt_temp')!r} cannot be used: {exc}") from None
    aedt = _plan_entry(plan, "aedt")
    try:
        launch_timeout_setting(aedt.get("launch_timeout_s"))
    except ValueError as exc:
        raise PlanRefused(f"the plan's aedt.launch_timeout_s: {exc}") from None

    warnings = []
    cpus = env.get("SLURM_CPUS_PER_TASK")
    if cpus:
        if not _DIGITS.fullmatch(cpus):
            raise PlanRefused(f"SLURM_CPUS_PER_TASK={cpus!r} is not a number of CPUs")
        if int(cpus) < cfg.solver.cores:
            raise PlanRefused(f"SLURM_CPUS_PER_TASK = {cpus}, fewer than the plan's {cfg.solver.cores} cores")
        if int(cpus) > cfg.solver.cores:
            warnings.append(f"SLURM_CPUS_PER_TASK = {cpus}, more than the plan's {cfg.solver.cores} cores; the solve "
                            f"uses {cfg.solver.cores}")
    attach = [name for name in ATTACH_VARIABLES if name in env]
    if attach:
        raise PlanRefused(f"{', '.join(attach)} set: PyAEDT 1.7.0 would attach to that desktop instead of starting "
                          "its own; the array script unsets every PYAEDT_* variable")
    guard, array_job = _check_guard(plan, env)
    job_id = f"a{array_job}_{index}" if array_job else f"{default_job_id(env)}_{index}"
    root = job_root(cfg.output.root, cfg.project.dataset_id, cfg.solver.frequency_ghz, job_id)
    if os.path.lexists(root):
        raise PlanRefused(f"the job directory {root} already exists")
    return TaskContext(plan=plan, index=index, cfg=cfg, reference=reference, job_id=job_id,
                       expected_lattice=expected_lattice, policy=policy, state=state, plan_path=plan_path.resolve(),
                       plan_sha256=digest, task=task, aedt_temp=temp,
                       launch_timeout_s=float(aedt.get("launch_timeout_s")), guard=guard, job_root=root,
                       array_job_id=array_job, warnings=tuple(warnings))


def slurm_time_limit(job_id: str | None) -> dict:
    """The job's ``TimeLimit`` as ``scontrol show job`` reports it, for ``start.json``: {value, error}. Never raises:
    a record must not stop the run."""
    if not job_id:
        return {"value": None, "error": "not in a Slurm job (no SLURM_JOB_ID)"}
    if not _DIGITS.fullmatch(job_id):
        return {"value": None, "error": f"SLURM_JOB_ID={job_id!r} is not a job id"}
    try:
        run = subprocess.run(["scontrol", "show", "job", job_id], capture_output=True, text=True,
                             timeout=SCONTROL_TIMEOUT_S, check=False)
    except (OSError, subprocess.SubprocessError) as exc:
        return {"value": None, "error": f"scontrol could not run ({exc})"}
    match = re.search(r"\bTimeLimit=(\S+)", run.stdout)
    if run.returncode != 0 or match is None:
        return {"value": None, "error": f"scontrol show job {job_id} exited {run.returncode} without a TimeLimit "
                                        f"({run.stderr.strip()[:200]})"}
    return {"value": match.group(1), "error": None}


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _write_json_atomic(path: Path, data: Any) -> None:
    """``data`` as JSON at ``path`` through a temporary file in the same directory and ``os.replace``, so a reader never
    sees part of it."""
    tmp = path.with_name(f".{path.name}.tmp")
    try:
        write_json(tmp, data)
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()


def start_record(task: TaskContext, dirs: JobDirs, argv: Sequence[str], env: Mapping[str, str]) -> dict:
    """``start.json``, written in the job before AEDT starts: the plan and task, the job, the Slurm context and time
    limit, the host, process and temp directory, the code, PyAEDT and its settings files, AEDT's temp mode, the guard."""
    home_settings = Path.home() / "pyaedt_settings.yaml"
    local = env.get("PYAEDT_LOCAL_SETTINGS_PATH") or None
    return {
        "plan": task.plan_record, "frequency_ghz": task.task.get("frequency_ghz"), "label": task.task.get("label"),
        "class": task.task.get("class"), "job_id": task.job_id, "job_root": str(dirs.root),
        "slurm": {name: env.get(name) for name in START_SLURM_VARIABLES},
        "time_limit": slurm_time_limit(env.get("SLURM_JOB_ID")),
        "host": socket.gethostname(), "pid": os.getpid(), "tmpdir": env.get("TMPDIR"), "argv": list(argv),
        "executable": sys.executable, "code": {"repo": task.plan["code"]["repo"], "commit": task.state.commit},
        "pyaedt_version": task.state.pyaedt_version,
        "pyaedt_settings": {"home_settings": str(home_settings), "home_settings_exists": home_settings.is_file(),
                            "local_settings_path": local,
                            "local_settings_exists": bool(local) and Path(local).exists()},
        "aedt_temp": {"mode": task.aedt_temp.mode,
                      "requested": tempdir.requested_directory(task.aedt_temp, dirs.root)},
        "guard": {"root": task.guard.get("root"), "reservation": task.guard.get("reservation")},
        "created_utc": _utc_now(),
    }


def done_record(task: TaskContext, manifest_path: Path, report: Path, passed: bool) -> dict:
    """``done.json``: the run completed and its report was written, PASS or a deterministic FAIL (HF-057). The digests
    of the manifest, sidecar and report bind the completion to those bytes."""
    manifest = json.loads(Path(manifest_path).read_text())
    aedt = manifest.get("aedt") if isinstance(manifest.get("aedt"), Mapping) else {}
    sidecar = Path(manifest_path).parent / dataset_record_name(task.cfg.project.dataset_id,
                                                               task.cfg.solver.frequency_ghz)
    return {"solve_completed": True, "acceptance_passed": bool(passed),
            "manifest_sha256": sha256_file(manifest_path), "sidecar_sha256": sha256_file(sidecar),
            "report_sha256": sha256_file(report), "report": str(report), "plan_sha256": task.plan_sha256,
            "index": task.index, "aedt": {"pid": aedt.get("pid"), "port": aedt.get("port")},
            "finished_utc": _utc_now()}


def guard_failure_line(task: TaskContext, exc: licence_guard.GuardError) -> str:
    """The last log line of a run whose licence guard stopped it inside ``run_job``, before AEDT started (OR-045).
    Only ``licence_guard.LeaseRefused``, ``acquire_lease``'s refusal for want of capacity, is called a refusal; in an
    array, a refusal because every slot of the reservation is held says that an earlier task's held lease, its
    cleanup unverified, makes it expected (HF-058(2)). Any other guard failure (the lock, a record, the launch scope)
    is a fault, named neutrally."""
    if not isinstance(exc, licence_guard.LeaseRefused):
        return f"job failed: the licence guard failed before AEDT started: {exc}; nothing was solved"
    root = task.guard.get("root")
    line = f"job failed: launch lease refused by the licence guard (OR-045), so nothing was solved: {exc}."
    if exc.reason == "slots" and task.array_job_id:
        return (f"{line} Under the array's reservation {task.guard.get('reservation')}, a launch lease an earlier task "
                f"still holds, its AEDT cleanup not yet verified, keeps its slot (HF-058(2)), so this refusal is "
                f"expected while it is held: list the job's leases with python -m bbsim.licence_guard --root {root} "
                f"status --job {task.array_job_id}, release a held one only after scope-check on its host (release "
                f"LEASE --confirm-cleanup REASON), then rerun task {task.index}")
    return (f"{line} python -m bbsim.licence_guard --root {root} status lists the launch leases holding the "
            f"capacity; rerun task {task.index} once one is released after its cleanup is verified")


def _main_plan(args: argparse.Namespace, cleanup: ExitStack, argv: Sequence[str]) -> int:
    try:
        task = plan_preflight(args.plan, args.plan_sha256, args.plan_index, os.environ)
        dirs = create_job_dirs(task.cfg.output.root, task.cfg.project.dataset_id, task.cfg.solver.frequency_ghz,
                               task.job_id)
    except (OSError, ValueError, ImportError, licence_guard.GuardError, EnvironmentStateError) as exc:
        print(f"refused before the job started: {exc}", file=sys.stderr)
        return 2
    cleanup.callback(_setup_logging(dirs.logs / "run.log", args.log_level))
    _write_json_atomic(dirs.root / START_NAME, start_record(task, dirs, argv, os.environ))
    write_json(dirs.root / "effective_config.json", config_to_dict(task.cfg))
    log.info("job %s -> %s: task %d (%s, class %s) of plan %s (sha256 %s), %s policy, schema %d", task.job_id,
             dirs.root, task.index, task.task.get("label"), task.task.get("class"), task.plan_path, task.plan_sha256,
             task.policy, task.cfg.output.schema)
    log.info(_guard_log_line(task.guard))
    log.info(_environment_log_line(task.state))
    if task.reference is not None:
        log.info("column-order reference from the plan's copy (%s, %s)", task.reference.get("source_job"),
                 task.reference.get("aedt_version"))
    for warning in task.warnings:
        log.warning(warning)
    try:
        manifest = Path(run_job(task.cfg, dirs, task.job_id, reference=task.reference,
                                aedt_temp=task.aedt_temp, launch_timeout_s=task.launch_timeout_s,
                                environment=task.state, expected_lattice=task.expected_lattice,
                                plan_record=task.plan_record))
    except SystemExit as exc:
        _log_system_exit("job", exc)
        raise
    except licence_guard.GuardError as exc:
        log.error(guard_failure_line(task, exc))
        return 1
    except Exception:  # noqa: BLE001 - report and exit non-zero; the session is already released
        log.exception("job failed")
        return 1
    return _judge(task, dirs, manifest)


def _judge(task: TaskContext, dirs: JobDirs, manifest: Path) -> int:
    """The array_run check of a completed run against its plan (``bbsim.run_acceptance``), its rows logged and its
    report written, then ``done.json`` and the run's last log line. A check that does not finish, or a report that
    cannot be written, leaves no ``done.json``: the task is incomplete, never a verdict (HF-057)."""
    from bbsim import run_acceptance

    try:
        result = run_acceptance.check(dirs.root, task.plan_path, task.plan_sha256, task.index, write=False)
        for line in result.format_rows().splitlines():
            log.info("%s", line)
        run_acceptance.write_report(result)
        for line in result.format_footer().splitlines():
            if line:
                log.info("%s", line)
    except SystemExit as exc:
        _log_system_exit("the array_run check", exc)
        raise
    except Exception:  # noqa: BLE001 - an incomplete task, rerun; the solved project stays in the job
        log.exception("job failed: the solve and export completed but the array_run check did not finish, so there "
                      "is no report and no %s; the solved project %s is kept", DONE_NAME, dirs.project)
        return 1
    if result.report is None:
        log.error("job failed: the solve and export completed but the array_run report was not written (%s), so "
                  "there is no %s", result.report_error, DONE_NAME)
        return 1
    try:
        _write_json_atomic(dirs.root / DONE_NAME, done_record(task, manifest, Path(result.report), result.passed))
    except (OSError, ValueError) as exc:
        log.error("job failed: %s could not be written (%s); the report %s stands", DONE_NAME, exc, result.report)
        return 1
    if result.passed:
        log.info("done; manifest at %s", manifest)
        return 0
    log.error("completed; acceptance FAIL (%s); manifest at %s", ", ".join(c.name for c in result.failures) or
              "the report", manifest)
    return 1


def _log_system_exit(what: str, exc: SystemExit) -> None:
    log.error("%s stopped by SystemExit(%s)%s; what the cleanup on the way out did (the session's release, the launch "
              "lease, AEDT's temp directory) is logged above", what, exc.code,
              " from SIGTERM" if exc.code == EXIT_SIGTERM else "")


def main(argv: list[str] | None = None) -> int:
    """The CLI. For the length of the run SIGTERM raises ``SystemExit(143)``, and the run's log handlers are removed
    and closed on the way out, so ``main`` can run twice in one process (HF-057)."""
    with ExitStack() as cleanup:
        cleanup.callback(_handle_sigterm())
        args = parse_args(argv)
        if args.plan is not None:
            return _main_plan(args, cleanup, list(sys.argv[1:] if argv is None else argv))
        if args.from_solved is not None:
            return _main_from_solved(args, cleanup)
        return _main_config(args, cleanup)


def _main_config(args: argparse.Namespace, cleanup: ExitStack) -> int:
    if args.config is None:
        print("refused: one of --config or --from-solved is required", file=sys.stderr)
        return 2
    cfg = load_config(args.config, overrides_from_args(args))
    if not cfg.project.path.is_file():
        print(f"prepared project not found: {cfg.project.path} (run prepare_hfss_project.py first, or pass --project)",
              file=sys.stderr)
        return 2
    reference = None
    try:
        if args.column_order_reference is not None:
            if cfg.output.schema != 2:
                raise ValueError(f"--column-order-reference applies to schema 2 only; {args.config} writes schema "
                                 f"{cfg.output.schema}")
            reference = load_reference(args.column_order_reference)
        check_before_aedt(cfg, reference)
        guard = check_guard_before_job()
        environment = environment_state(REPO_ROOT)
    except (OSError, ValueError, ImportError, licence_guard.GuardError, EnvironmentStateError) as exc:
        print(f"refused before the job started: {exc}", file=sys.stderr)
        return 2
    job_id = args.job_id or default_job_id()
    dirs = create_job_dirs(cfg.output.root, cfg.project.dataset_id, cfg.solver.frequency_ghz, job_id)
    cleanup.callback(_setup_logging(dirs.logs / "run.log", args.log_level))
    write_json(dirs.root / "effective_config.json", config_to_dict(cfg))
    log.info("job %s -> %s (schema %d)", job_id, dirs.root, cfg.output.schema)
    log.info(_guard_log_line(guard))
    (log.warning if environment.dirty else log.info)(_environment_log_line(environment))
    if reference is not None:
        log.info("column-order reference from %s (%s, %s)", args.column_order_reference, reference.get("source_job"),
                 reference.get("aedt_version"))
    warning = allocation_warning(cfg.solver.cores)
    if warning:
        log.warning(warning)
    return _run_logged("job", functools.partial(
        run_job, cfg, dirs, job_id, reference=reference, aedt_temp=args.temp_settings,
        launch_timeout_s=args.launch_timeout_s, environment=environment))


if __name__ == "__main__":
    raise SystemExit(main())
