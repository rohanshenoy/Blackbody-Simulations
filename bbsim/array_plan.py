"""The job-array plan (HF-057, OR-037, OR-044, OR-045, OR-048, S65): the planner, a task's predicted exit lattice and
one task's frozen configuration.

``plan_hfss_array.py`` (``main``) runs on the login node and opens no AEDT. ``build_plan`` checks everything first and
creates nothing until every check has passed; each refusal is one ``PlanError`` line (exit 2):

- ``--run-dir`` does not exist and its name is a plain word (it names the Slurm jobs); ``--tree`` likewise, required
  unless ``--probe`` and refused with it;
- the code checkout is clean at the commit the plan records (``manifest.environment_state``, one seam, injected by the
  tests); the bridge is the checkout's sibling clone ``../bbridge`` and no other repository (``--scripts-dir`` may only
  name its ``scripts/hpc/hfss``; OR-044), with no uncommitted change to a tracked file and every script the plan names
  in its ``HEAD``;
- the grid file parses with at most ``--max-array-size`` entries, and each task's effective step along l is at most
  c/(4f) (``frequency_grid``, OR-048);
- the config loads with no ``solver.frequency_ghz`` and declares ``[acceptance] kind = "array_run"``, so a
  reference-based config is refused (HF-057); it meets BBR025, the rules BBRsim's sidecar loader enforces (symmetry
  ``mirror_l``, ``mirror_g`` and ``end_to_end``, pose ``canonical-z``, ``exit_cs_x`` = z and ``exit_cs_y`` = y), and its
  geometry is the canonical rectangle ``predict_exit_lattice`` requires;
- a production plan takes a schema-2 config of both polarizations whose incidence grid is ``--expect-incidence``
  (OR-036's 35 keys by default) and the HF-041 ``--column-order-reference``; a probe plan the schema-1 single-angle
  diagnostic of one key and Ephi=0 (HF-057, HF-058), with none of ``--tree``, ``--include-job``, a reference or an
  expectation;
- every task's config, built by ``task_config`` from the plan being written, passes the runner's own pre-AEDT check;
- the project exists, and the record beside it (``.mapping.json`` from prepare_project, or ``.build.json``) gives the
  config's dataset_id for its design and object (``mapping.recorded_dataset_id``, the runner's check, which a single
  run only warns about); cores, the throttle L, the classes (which partition the task indices, their throttles summing
  to L), the resources, the stagger and the AEDT temp settings are well formed;
- every ``--sbatch-arg`` passes ``licence_guard.check_user_sbatch_args`` and sets none of ``PLAN_RESOURCE_OPTIONS``,
  which the plan sets itself; ``--licence-cap`` is 1 to 10; ``--guard-dir``
  and an inherited ``BBSIM_LICENCE_GUARD`` are the guard's root; every pending or running array the guard manages
  (comment ``bbsim:``) runs at the plan's ``--cores``, read from ``scontrol show job``, since arrays at two core counts
  share the per-user ``pyaedt_config`` registry entry; and the guard's reconciled usage plus L fits the smaller of the
  ceiling and the cap. That count is advisory: admission is repeated at every submission (OR-045);
- each ``--include-job`` (HF-057, S65) is a schema-2 job of another dataset ID on the expected incidence grid, qualified
  for both polarizations (``acceptance_gate.convergence_problems``), whose files are its sidecar's and whose recorded
  acceptance passes the gate (``required_acceptance``, then ``report_problems`` on its bound report); an ``array_run``
  job's manifest names its plan, whose bytes must have the recorded sha256 and whose policy must be production.

Then it creates ``<RUN>/{inputs,jobs,slurm}`` (and any missing parent), the detached worktree ``<RUN>/code`` at the
commit, where the reference is checked against the exporter files (``dataset_v2.check_reference_before_aedt``), the
byte copies of the config, the reference, the project and its record (beside the project copy under its own name, so
each task's ``run_frequency.stage_project`` confirms the dataset again) in ``<RUN>/inputs/``, and ``plan.json`` and
``plan.sha256`` last. A failure removes everything it created, ``git worktree remove --force`` included, and is a
refusal; only a cleanup that cannot finish is ``PlanLeftovers`` (exit 1), naming what is left.

``plan.json`` (``PLAN_SCHEMA``): schema, created_utc, name, run_dir, tree (null for a probe), dataset_id, probe, policy
(``production`` or ``probe``), acceptance_kind (``array_run``), code {repo (the worktree), source_repo, commit,
pyaedt_version}, bridge {repo, commit, scripts {repo-relative path: sha256 of the HEAD blob}}, config, project and
column_order_reference (null for a probe) {source, copy, sha256}, the project also its record {source, copy, sha256},
grid {source, sha256, entries}, order, tasks [{index, frequency_ghz, label, class, l_resolution_mm (the grid entry's
requested step, may be null), effective_l_resolution_mm (the applied one), l_step_mm (the fitted one, the lattice's
step along l), lattice, aedt_temp {mode, template, restore, strict}}], classes [{name, indices, mem, time, qos,
throttle}], slurm {cores, throttle, stagger {step_s, modulus}, sbatch_args, tmpdir, job_name, comment, max_array_size,
assembly {mem, time, qos}, monitor {enabled, interval_s, time}}, incidence and expect_incidence {theta_deg, phi_deg},
aedt {temp, launch_timeout_s}, include_jobs [{job, dataset_id, frequency_ghz, acceptance_kind, files {path: sha256},
source_plan {path, sha256, index} or null}], licence_cap, guard_dir, jobs_root. Every path it records is absolute and
resolved; an include job's file inside the job is keyed by its path relative to the job, a file outside it (a reference
job's, a plan) by its absolute path.

``task_config`` builds a task's run configuration from the config copy, or from the caller's raw TOML mapping already
captured from that copy, with ``config.build_config`` and exactly five overrides through the validated override path:

- ``solver.frequency_ghz``: the task's ``frequency_ghz``;
- ``project.path``: the plan's ``project.copy``;
- ``output.root``: the plan's ``jobs_root`` (``<RUN>/jobs``);
- ``solver.cores``: the plan's ``slurm.cores``;
- ``exit_field.resolution_mm``: the frozen config's, with index 1 (the step along l; the list runs along p, l, g) the
  task's ``effective_l_resolution_mm``. That field is never null: it is the grid entry's requested step when the
  entry gives one (``l_resolution_mm``), else the config's own, after the c/(4f) rule (``frequency_grid``).

After writing the plan the CLI prints its submission block (``submission_block``; OR-044, OR-045, HF-057): from the
frozen code ``<RUN>/code``, one ``python -m bbsim.array_submit submit`` line per class, the licence monitor when the
plan runs one, and the assembly after every array; then the rerun and throttle lines. Every job goes through the guard
from the recorded bridge blob, and nothing in the block calls sbatch or scontrol. A line refused part way through
prints what remains of the block, with the ids already admitted (``python -m bbsim.array_submit remaining``).
``--submit`` does the same at once (``array_submit.submit_plan``); ``plan_hfss_array.py --log-submission RUN KIND
ID`` records a job submitted by hand, with no admission (``array_submit.log_submission``).

The planner, the runner, the array_run self-check (``bbsim.run_acceptance``) and the assembler all build a task's
config here, and the planner and the self-check predict its lattice with ``predict_exit_lattice``, so none can
disagree on either. A caller that hashes the config copy passes the mapping it parsed from those bytes as
``raw_config``, and the copy is then never reopened. The module imports only the standard library and standard-library
``bbsim`` modules at load time; numpy and the extractor are imported where the lattice is predicted.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import tomllib
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError
from pathlib import Path
from typing import Any

from bbsim import acceptance_gate, aedt_temp, frequency_grid, licence_guard
from bbsim.config import ConfigError, ExitFieldConfig, RunConfig, build_config
from bbsim.manifest import (GIT_TIMEOUT_S, REPO_ROOT, EnvironmentState, EnvironmentStateError, environment_state,
                            sha256_file)
from bbsim.naming import frequency_label
from bbsim.session import DEFAULT_LAUNCH_TIMEOUT_S, launch_timeout_setting

log = logging.getLogger(__name__)

PLAN_SCHEMA = "bbsim-array-plan/2"
ACCEPTANCE_KIND = "array_run"
LATTICE_SOURCE = "config.geometry.expected_box_size_mm"
AREA_REL = 1e-6                  # the rectangle's area against expected_face_area_mm2 (geometry.check_face_area's)
DEFAULT_THROTTLE = {"production": 2, "probe": 3}
EXPECT_THETA_DEG = (90.0, 105.0, 120.0, 135.0, 150.0, 165.0, 180.0)    # OR-036, crack-local
EXPECT_PHI_DEG = (0.0, 22.5, 45.0, 67.5, 90.0)
ANGLE_TOL_DEG = 1e-9
ORDERS = ("longest-first", "ascending")
TMPDIRS = ("run", "node", "off")
BBR025_SYMMETRY = ("mirror_l", "mirror_g", "end_to_end")
CANONICAL_EXIT_CS = ((0.0, 0.0, 1.0), (0.0, 1.0, 0.0))                  # exit_cs_x, exit_cs_y under canonical-z
# The bridge scripts a plan names (OR-044); the monitor only when the plan runs one.
ARRAY_SCRIPT = "run_frequency_array.sbatch"
ASSEMBLY_SCRIPT = "assemble_dataset.sbatch"
MONITOR_SCRIPT = "monitor_licences.sbatch"
MONITOR_MARGIN_S = 1800          # the monitor's default time is the plan's duration estimate plus 30 min
# The resources the plan sets on every job itself, and the planner option that sets each: a --sbatch-arg may not set
# them again, since sbatch would take one of the two values silently (submission refuses them too, HF-057).
PLAN_RESOURCE_OPTIONS = {"--cpus-per-task": "--cores", "--mem": "--mem (or --class mem=)",
                         "--time": "--time (or --class time=)", "--qos": "--qos (or --class qos=)"}
JOB_DIRS = ("inputs", "jobs", "slurm")

_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*")
_MEM = re.compile(r"[1-9][0-9]*[KMGT]?")
_QOS = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*")
_TIME = re.compile(r"(?:([0-9]+)-)?([0-9]+)(?::([0-9]+))?(?::([0-9]+))?")
_INDICES = re.compile(r"([0-9]{1,7})(?:-([0-9]{1,7})(?::([0-9]{1,7}))?)?")
_STAGGER = re.compile(r"([0-9]+)(?::([0-9]+))?")
_CONTROL = re.compile(r"[\x00-\x1f\x7f]")
_COMMIT = re.compile(r"[0-9a-f]{40}|[0-9a-f]{64}")
_CLASS_KEYS = ("indices", "mem", "time", "qos", "throttle")
_CORE_STATES = frozenset({"PENDING", "RUNNING"})       # the queue states whose arrays the cores rule compares


class PlanError(ValueError):
    """The plan is refused, and nothing the planner made is left (exit 2)."""


class LatticeError(PlanError):
    """A config's exit lattice cannot be predicted: it is not the canonical rectangle, or a step is unusable."""


class PlanLeftovers(RuntimeError):
    """A plan failed and what the planner had made could not all be removed (exit 1); the message names it."""


# --- one task's frozen configuration ---------------------------------------------------------------------------------

def _absolute(value: object, what: str) -> str:
    if not isinstance(value, str) or not os.path.isabs(value):
        raise ValueError(f"the plan's {what} must be an absolute path, as every path a plan records is, got {value!r}")
    return value


def _positive(value: object, what: str) -> float | int:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not (math.isfinite(value) and value > 0):
        raise ValueError(f"{what} must be a finite positive number, got {value!r}")
    return value


def _copy(plan: Mapping[str, Any], name: str) -> str:
    entry = plan.get(name)
    if not isinstance(entry, Mapping):
        raise ValueError(f"the plan records no {name} {{source, copy, sha256}}, got {entry!r}")
    return _absolute(entry.get("copy"), f"{name}.copy")


def _task(plan: Mapping[str, Any], index: object) -> Mapping[str, Any]:
    if isinstance(index, bool) or not isinstance(index, int) or index < 0:
        raise ValueError(f"a task index is a nonnegative integer, got {index!r}")
    tasks = plan.get("tasks")
    if not isinstance(tasks, list):
        raise ValueError(f"the plan's tasks must be a list of task records, got {type(tasks).__name__}")
    found = [t for t in tasks if isinstance(t, Mapping) and type(t.get("index")) is int and t["index"] == index]
    if len(found) != 1:
        raise ValueError(f"the plan holds {len(found)} records of task {index}, not 1")
    return found[0]


def task_config(plan: Mapping, index: int, *, raw_config: Mapping | None = None) -> RunConfig:
    """Task ``index``'s run configuration: the plan's frozen config copy (or ``raw_config``, the caller's mapping parsed
    from that copy's captured bytes, in which case the copy is not reopened) built with the five overrides of the module
    docstring. Raises ``ValueError`` for a plan or task record it cannot use (``ConfigError``, a ``ValueError``, for a
    config the overrides cannot build); an unreadable copy raises ``OSError``."""
    if not isinstance(plan, Mapping):
        raise ValueError(f"a plan is an object, got {type(plan).__name__}")
    if plan.get("schema") != PLAN_SCHEMA:
        raise ValueError(f"the plan's schema is {plan.get('schema')!r}, not {PLAN_SCHEMA!r}")
    task = _task(plan, index)
    where = f"task {index}'s"
    frequency = _positive(task.get("frequency_ghz"), f"{where} frequency_ghz")
    effective = float(_positive(task.get("effective_l_resolution_mm"), f"{where} effective_l_resolution_mm"))
    project = _copy(plan, "project")
    jobs_root = _absolute(plan.get("jobs_root"), "jobs_root")
    slurm = plan.get("slurm")
    cores = slurm.get("cores") if isinstance(slurm, Mapping) else None
    if isinstance(cores, bool) or not isinstance(cores, int) or cores < 1:
        raise ValueError(f"the plan's slurm.cores must be a positive integer, got {cores!r}")
    if raw_config is None:
        with open(_copy(plan, "config"), "rb") as fh:
            raw_config = tomllib.load(fh)
    elif not isinstance(raw_config, Mapping):
        raise ValueError(f"raw_config must be the config copy's parsed TOML mapping, got {type(raw_config).__name__}")
    overrides: dict[str, object] = {"solver.frequency_ghz": frequency, "project.path": project,
                                    "output.root": jobs_root, "solver.cores": cores}
    resolution = list(build_config(raw_config, overrides).exit_field.resolution_mm)
    resolution[1] = effective
    return build_config(raw_config, {**overrides, "exit_field.resolution_mm": resolution})


# --- the predicted exit lattice (HF-045, HF-057) ----------------------------------------------------------------------

def canonical_extent(cfg: RunConfig) -> tuple[float, float, float]:
    """The frozen box's extents along the exit axes (p, l, g), in mm: ``geometry.expected_box_size_mm`` projected on
    ``exit_cs_x``, ``exit_cs_y`` and their cross product (``dataset_record.exit_axes``). ``LatticeError`` without a
    box, or without a positive extent along l and g."""
    import numpy as np

    from bbsim.dataset_record import exit_axes

    box = cfg.geometry.expected_box_size_mm
    if box is None:
        raise LatticeError("the frozen config gives no geometry.expected_box_size_mm, so the plan's canonical "
                           "rectangle (its lattice and its entrance area A) cannot be derived")
    sizes = np.asarray(box, dtype=float)
    p, l, g = (float(np.sum(np.abs(row) * sizes)) for row in exit_axes(cfg.geometry.exit_cs_x, cfg.geometry.exit_cs_y))
    if not all(math.isfinite(v) and v > 0.0 for v in (l, g)):
        raise LatticeError(f"the frozen box {list(box)} mm has no positive extent along l and g ({l!r}, {g!r} mm)")
    return p, l, g


def _lattice_step(value: object, what: str) -> float:
    """A lattice step in mm: a finite, positive and normal float (a subnormal one would overflow the node count)."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise LatticeError(f"{what} must be a number of mm, got {value!r}")
    step = float(value)
    if not (math.isfinite(step) and step > 0.0):
        raise LatticeError(f"{what} {value!r} mm must be a finite positive number")
    if step < sys.float_info.min:
        raise LatticeError(f"{what} {value!r} mm is subnormal, too small to step a lattice")
    return step


def predict_exit_lattice(cfg: RunConfig, effective_l_resolution_mm: float) -> dict:
    """The exit lattice a task of ``cfg`` exports at the step along l ``effective_l_resolution_mm``, without AEDT.

    The config must describe the canonical rectangle (rectangles only, HF-057): ``geometry.expected_box_size_mm``
    given, no rotational symmetry (a disc), no curved-wall facets (``solver.wall_normal_deviation_deg``), no hand-set
    ``exit_field.boundary_mm``, and any ``geometry.expected_face_area_mm2`` equal to the rectangle's l x g to
    ``AREA_REL``. Its face, the four corners (0, +-l/2, +-g/2) in the exit frame (x_e, y_e, z_e) = (p, l, g), goes
    through the runner's own ``extract.exit_grid`` with the config's resolution and index 1 replaced, so each face axis
    gets HF-045's fitted step. Returns ``{lo_mm, hi_mm, step_mm, counts, fitted, extent_mm (p, l, g), entrance_area_mm2,
    extent_source}``, bounds centred on the exit face. Every refusal is a ``LatticeError``, a non-finite, non-positive
    or subnormal step at once."""
    import numpy as np

    from bbsim.extract import ExitFaceGeometry, classify_cross_section, exit_grid

    geometry = cfg.geometry
    shape = []
    if "rotational" in geometry.symmetry:
        shape.append("geometry.symmetry declares rotational, which a disc has")
    if cfg.solver.wall_normal_deviation_deg is not None:
        shape.append(f"solver.wall_normal_deviation_deg = {cfg.solver.wall_normal_deviation_deg!r} sets the facets "
                     "of a curved wall")
    if cfg.exit_field.boundary_mm is not None:
        shape.append("exit_field.boundary_mm sets the lattice's bounds by hand, not from the face")
    if shape:
        raise LatticeError("the plan predicts the exit lattice of a canonical rectangle only (HF-057): "
                           + "; ".join(shape))
    p, l, g = canonical_extent(cfg)
    area = l * g
    declared = geometry.expected_face_area_mm2
    if declared is not None and not abs(float(declared) - area) <= AREA_REL * area:
        raise LatticeError(f"geometry.expected_face_area_mm2 = {declared!r} mm2 is not the rectangle's l x g = "
                           f"{area!r} mm2 of geometry.expected_box_size_mm {list(geometry.expected_box_size_mm)} mm "
                           f"(to {AREA_REL:g} relative): the face is not that rectangle")
    resolution = [float(v) for v in cfg.exit_field.resolution_mm]
    resolution[1] = _lattice_step(effective_l_resolution_mm, "the effective step along l")
    for axis, extent, name in ((1, l, "l"), (2, g, "g")):
        if resolution[axis] > 0.0 and not math.isfinite(extent / _lattice_step(resolution[axis],
                                                                                f"the step along {name}")):
            raise LatticeError(f"the step along {name} {resolution[axis]!r} mm gives no finite node count over "
                               f"{extent!r} mm")
    points = np.array([[0.0, y, z] for y in (-l / 2.0, l / 2.0) for z in (-g / 2.0, g / 2.0)])
    face = ExitFaceGeometry(points, 4, "vertices", classify_cross_section(points, 4), 0.0)
    if face.cross_section["shape"] != "rectangle":
        raise LatticeError(f"the box's exit face {l!r} x {g!r} mm does not classify as a rectangle "
                           f"({face.cross_section['shape']})")
    grid = exit_grid(face, ExitFieldConfig(resolution_mm=tuple(resolution), edge_samples=cfg.exit_field.edge_samples))
    return {"lo_mm": [float(v) for v in grid.lo_mm], "hi_mm": [float(v) for v in grid.hi_mm],
            "step_mm": [float(v) for v in grid.step_mm], "counts": [int(n) for n in grid.counts],
            "fitted": [bool(f) for f in grid.fitted], "extent_mm": [p, l, g], "entrance_area_mm2": area,
            "extent_source": LATTICE_SOURCE}


# --- the planner's options -------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class PlanOptions:
    """Every option of the planner's CLI, as given (``build_plan`` checks and parses them). ``None`` takes the
    default: the config's cores, the policy's throttle (3 probe, 2 production), the monitor on for a probe plan only,
    its time the plan's duration estimate plus 30 min, the sibling bridge clone's scripts folder. ``classes`` holds one
    tuple of ``KEY=VALUE`` tokens per ``--class``; ``submit`` is acted on by the submission step, never here."""

    config: Path
    grid: Path
    run_dir: Path
    tree: Path | None = None
    project: Path | None = None
    cores: int | None = None
    mem: str = "32G"
    time: str = "08:00:00"
    qos: str = "normal"
    throttle: int | None = None
    classes: tuple[tuple[str, ...], ...] = ()
    order: str = "longest-first"
    max_array_size: int = frequency_grid.MAX_ENTRIES
    stagger: str = "5:8"
    aedt_temp: str = "task"
    aedt_temp_index: tuple[str, ...] = ()
    aedt_temp_restore: str | None = None
    aedt_temp_strict: bool = False
    launch_timeout_s: float = DEFAULT_LAUNCH_TIMEOUT_S
    tmpdir: str = "run"
    expect_incidence: str | None = None
    licence_cap: int = licence_guard.CEILING
    include_jobs: tuple[Path, ...] = ()
    column_order_reference: Path | None = None
    scripts_dir: Path | None = None
    guard_dir: Path = licence_guard.DEFAULT_ROOT
    sbatch_args: tuple[str, ...] = ()
    monitor: bool | None = None
    monitor_interval_s: int = 60
    monitor_time: str | None = None
    assembly_mem: str = "4G"
    assembly_time: str = "02:00:00"
    assembly_qos: str | None = None
    probe: bool = False
    submit: bool = False


# --- small checks -----------------------------------------------------------------------------------------------------

def _one_line(text: object) -> str:
    return " ".join(str(text).split())


def _integer(value: object, what: str, lo: int, hi: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < lo or (hi is not None and value > hi):
        bounds = f"from {lo} to {hi}" if hi is not None else f"of at least {lo}"
        raise PlanError(f"{what} must be an integer {bounds}, got {value!r}")
    return value


def _mem(value: object, what: str) -> str:
    if not (isinstance(value, str) and _MEM.fullmatch(value)):
        raise PlanError(f"{what} takes a Slurm memory size such as 32G or 4000M, got {value!r}")
    return value


def _qos(value: object, what: str) -> str:
    if not (isinstance(value, str) and _QOS.fullmatch(value)):
        raise PlanError(f"{what} takes a QOS name, got {value!r}")
    return value


def _seconds(value: object, what: str) -> int:
    """A Slurm time limit in seconds: MM, MM:SS, HH:MM:SS, D-HH, D-HH:MM or D-HH:MM:SS, above 0."""
    match = _TIME.fullmatch(value) if isinstance(value, str) else None
    if match:
        days, a, b, c = match.groups()
        if days is not None:
            h, m, s = int(a), int(b or 0), int(c or 0)
            ok, total = h < 24 and m < 60 and s < 60, int(days) * 86400 + h * 3600 + m * 60 + s
        elif c is not None:
            h, m, s = int(a), int(b), int(c)
            ok, total = m < 60 and s < 60, h * 3600 + m * 60 + s
        elif b is not None:
            ok, total = int(b) < 60, int(a) * 60 + int(b)
        else:
            ok, total = True, int(a) * 60
        if ok and total > 0:
            return total
    raise PlanError(f"{what} takes a Slurm time limit above 0 (MM, MM:SS, HH:MM:SS, D-HH, D-HH:MM or D-HH:MM:SS), "
                    f"got {value!r}")


def _time_text(seconds: int) -> str:
    days, rest = divmod(int(seconds), 86400)
    hours, rest = divmod(rest, 3600)
    minutes, secs = divmod(rest, 60)
    return f"{days}-{hours:02d}:{minutes:02d}:{secs:02d}" if days else f"{hours:02d}:{minutes:02d}:{secs:02d}"


def _stagger(value: object) -> dict:
    match = _STAGGER.fullmatch(value) if isinstance(value, str) else None
    if match is None or (match.group(2) is not None and int(match.group(2)) < 1):
        raise PlanError(f"--stagger takes STEP[:MODULUS], whole seconds and a modulus of at least 1, got {value!r}")
    return {"step_s": int(match.group(1)), "modulus": int(match.group(2) or 8)}


def _angles(text: str, what: str) -> list[float]:
    values = []
    for token in text.split(","):
        try:
            value = float(token)
        except ValueError:
            value = math.nan
        if not math.isfinite(value):
            raise PlanError(f"--expect-incidence: {what} {token!r} is not a finite number of degrees")
        values.append(value)
    if any(b <= a for a, b in zip(values, values[1:])):
        raise PlanError(f"--expect-incidence: the {what} list {values} must be strictly ascending")
    return values


def _expected_incidence(text: str | None) -> dict:
    if text is None:
        return {"theta_deg": list(EXPECT_THETA_DEG), "phi_deg": list(EXPECT_PHI_DEG)}
    thetas, sep, phis = text.partition(":") if isinstance(text, str) else ("", "", "")
    if not sep or not thetas or not phis:
        raise PlanError(f"--expect-incidence takes THETAS:PHIS, two comma lists of degrees, got {text!r}")
    return {"theta_deg": _angles(thetas, "theta"), "phi_deg": _angles(phis, "phi")}


def _same_angles(a: Sequence[float], b: Sequence[float]) -> bool:
    return len(a) == len(b) and all(abs(float(x) - float(y)) <= ANGLE_TOL_DEG for x, y in zip(a, b))


def _new_path(value: object, option: str) -> Path:
    """An absolute, resolved path that does not exist yet (nor does the name given for it)."""
    if not isinstance(value, (str, os.PathLike)) or not os.fspath(value):
        raise PlanError(f"{option} must name a path, got {value!r}")
    given = Path(os.path.abspath(os.path.expanduser(os.fspath(value))))
    if _CONTROL.search(str(given)):
        raise PlanError(f"{option} {str(given)!r} holds a control character")
    path = given.resolve()
    for candidate in (given, path):
        if os.path.lexists(candidate):
            raise PlanError(f"{option} {candidate} exists; a plan never reuses or overwrites a directory")
    return path


def _existing_file(value: object, what: str) -> Path:
    try:
        path = Path(os.path.expanduser(os.fspath(value))).resolve(strict=True)    # type: ignore[arg-type]
    except (OSError, RuntimeError, TypeError) as exc:
        raise PlanError(f"{what} {value} does not exist or cannot be reached ({exc})") from None
    if not path.is_file():
        raise PlanError(f"{what} {path} is not a file")
    return path


def _git(repo: Path, *args: str, text: bool = True) -> subprocess.CompletedProcess:
    command = ["git", "-C", str(repo), "--no-optional-locks", *args]
    try:
        return subprocess.run(command, capture_output=True, text=text, timeout=GIT_TIMEOUT_S)
    except subprocess.TimeoutExpired:
        raise PlanError(f"{' '.join(command)} timed out after {GIT_TIMEOUT_S} s") from None
    except OSError as exc:
        raise PlanError(f"git could not run on {repo}: {exc}") from None


def _git_out(repo: Path, *args: str) -> str:
    run = _git(repo, *args)
    if run.returncode != 0:
        raise PlanError(f"git -C {repo} {' '.join(args)} failed with exit {run.returncode}: {_one_line(run.stderr)}")
    return run.stdout


def _worktree_paths(repo: Path) -> set[Path]:
    run = _git(repo, "worktree", "list", "--porcelain")
    if run.returncode != 0:
        raise PlanError(f"git worktree list in {repo} failed: {_one_line(run.stderr)}")
    return {Path(line[len("worktree "):]).resolve() for line in run.stdout.splitlines() if line.startswith("worktree ")}


def _temp_record(settings: aedt_temp.TempSettings) -> dict:
    return {"mode": settings.mode, "template": settings.template, "restore": settings.restore,
            "strict": settings.strict}


# --- the planner's checks --------------------------------------------------------------------------------------------

def _sibling_bridge(repo: Path) -> Path:
    """The bridge clone beside the code checkout ``repo``, ``../bbridge`` (OR-044); for this checkout it is
    ``licence_guard.default_bridge_repo()``, the one clone the guard's path mode takes."""
    return Path(repo).resolve().parent / licence_guard.default_bridge_repo().name


def _bridge(scripts_dir: object, names: Sequence[str], repo: Path) -> dict:
    """The code checkout ``repo``'s sibling bridge clone, whose ``scripts/hpc/hfss`` ``scripts_dir`` names (default
    that folder): clean, with every named script in its ``HEAD``, recorded as {repo, commit, scripts {repo-relative
    path: sha256 of the HEAD blob}}. No other repository is taken, so no other repository's committed scripts can
    reach sbatch through a plan (OR-044)."""
    sibling = _sibling_bridge(repo)
    folder = Path(os.fspath(scripts_dir)) if scripts_dir is not None else sibling / licence_guard.BRIDGE_SCRIPTS_DIR
    try:
        folder = folder.resolve(strict=True)
    except (OSError, RuntimeError) as exc:
        raise PlanError(f"the bridge's scripts folder {folder} cannot be reached ({exc}); --scripts-dir names it, by "
                        f"default the sibling bridge clone's {licence_guard.BRIDGE_SCRIPTS_DIR} (OR-044)") from None
    if not folder.is_dir():
        raise PlanError(f"--scripts-dir {folder} is not a directory")
    top = Path(_git_out(folder, "rev-parse", "--show-toplevel").strip()).resolve()
    try:
        expected = sibling.resolve(strict=True)
    except (OSError, RuntimeError) as exc:
        raise PlanError(f"the code checkout's sibling bridge clone {sibling} cannot be reached ({exc}) "
                        "(OR-044)") from None
    if top != expected:
        raise PlanError(f"--scripts-dir {folder} lies in the git repository {top}, not in the code checkout's sibling "
                        f"bridge clone {expected}: a plan takes its scripts from that clone only (OR-044)")
    if folder != top / licence_guard.BRIDGE_SCRIPTS_DIR:
        raise PlanError(f"--scripts-dir {folder} is not {licence_guard.BRIDGE_SCRIPTS_DIR} of its bridge clone {top}")
    if _git_out(top, "status", "--porcelain", "--untracked-files=no").strip():
        raise PlanError(f"the bridge clone {top} has uncommitted changes to tracked files; commit or restore them, "
                        "since the plan records the scripts of its HEAD (OR-044)")
    commit = _git_out(top, "rev-parse", "--verify", "HEAD").strip()
    scripts = {}
    for name in names:
        rel = (licence_guard.BRIDGE_SCRIPTS_DIR / name).as_posix()
        blob = _git(top, "cat-file", "blob", f"{commit}:{rel}", text=False)
        if blob.returncode != 0:
            raise PlanError(f"the bridge clone {top} has no {rel} at its HEAD {commit}; every script the plan names "
                            "is committed there (OR-044)")
        scripts[rel] = hashlib.sha256(blob.stdout).hexdigest()
    return {"repo": str(top), "commit": commit, "scripts": scripts}


def _parse_indices(text: str, n: int, where: str) -> list[int]:
    indices: list[int] = []
    for part in text.split(","):
        match = _INDICES.fullmatch(part)
        if match is None:
            raise PlanError(f"{where}: {part!r} is not an index, a range a-b or a stepped range a-b:s")
        first, last = int(match.group(1)), int(match.group(2) or match.group(1))
        step = int(match.group(3) or 1)
        if last < first or step < 1:
            raise PlanError(f"{where}: the range {part!r} runs backwards or has no step")
        if last >= n:
            raise PlanError(f"{where}: the plan has tasks 0 to {n - 1}, not {last}")
        indices += range(first, last + 1, step)
    if len(set(indices)) != len(indices):
        raise PlanError(f"{where}: an index repeats")
    return sorted(indices)


def _classes(o: PlanOptions, n: int, default_throttle: int) -> tuple[list[dict], int]:
    """The resource classes and the plan's throttle L: one class of every task unless ``--class`` gives several,
    which must partition the task indices with throttles summing to L (HF-057)."""
    mem, time, qos = _mem(o.mem, "--mem"), o.time, _qos(o.qos, "--qos")
    _seconds(time, "--time")
    throttle = None if o.throttle is None else _integer(o.throttle, "--throttle", 1)
    specs: list[dict[str, Any]] = []
    for number, tokens in enumerate(o.classes or (("indices=" + (f"0-{n - 1}" if n > 1 else "0"),),), start=1):
        where = f"--class {' '.join(tokens)}"
        spec: dict[str, Any] = {}
        for token in tokens:
            key, sep, value = token.partition("=")
            if not sep or key not in _CLASS_KEYS or key in spec or not value:
                raise PlanError(f"{where}: {token!r} is not one of {', '.join(k + '=' for k in _CLASS_KEYS)} given "
                                "once with a value")
            spec[key] = value
        if "indices" not in spec:
            raise PlanError(f"{where}: indices= is required")
        entry = {"name": f"A{number}", "indices": _parse_indices(spec["indices"], n, where),
                 "mem": _mem(spec.get("mem", mem), f"{where}: mem"), "time": spec.get("time", time),
                 "qos": _qos(spec.get("qos", qos), f"{where}: qos"), "throttle": None}
        _seconds(entry["time"], f"{where}: time")
        if "throttle" in spec:
            text = spec["throttle"]
            entry["throttle"] = _integer(int(text) if re.fullmatch(r"[0-9]{1,7}", text) else text,
                                         f"{where}: throttle", 1)
        specs.append(entry)
    owner: dict[int, str] = {}
    for entry in specs:
        for i in entry["indices"]:
            if i in owner:
                raise PlanError(f"--class: task index {i} is in two classes ({owner[i]} and {entry['name']}); the "
                                "classes partition the task indices")
            owner[i] = entry["name"]
    missing = sorted(set(range(n)) - set(owner))
    if missing:
        raise PlanError(f"--class: task indices {missing} are in no class; the classes partition the task indices 0 "
                        f"to {n - 1}")
    if len(specs) == 1:
        (entry,) = specs
        if throttle is not None and entry["throttle"] is not None and entry["throttle"] != throttle:
            raise PlanError(f"--class throttle={entry['throttle']} is not --throttle {throttle}; one class's throttle "
                            "is the plan's")
        total = throttle or entry["throttle"] or default_throttle
        if total > n:
            log.info("the throttle %d is clamped to the plan's %d tasks", total, n)
            total = n
        entry["throttle"] = total
        return specs, total
    if any(entry["throttle"] is None for entry in specs):
        raise PlanError("--class: throttle= is required in every class when several are given")
    total = sum(entry["throttle"] for entry in specs)
    if throttle is not None and throttle != total:
        raise PlanError(f"--class: the class throttles {[e['throttle'] for e in specs]} sum to {total}, not "
                        f"--throttle {throttle}")
    for entry in specs:
        if entry["throttle"] > len(entry["indices"]):
            raise PlanError(f"--class {entry['name']}: throttle {entry['throttle']} exceeds its "
                            f"{len(entry['indices'])} tasks")
    return specs, total


def _temp_settings(o: PlanOptions, n: int) -> tuple[aedt_temp.TempSettings, dict[int, aedt_temp.TempSettings]]:
    try:
        mode, template = aedt_temp.parse_mode(o.aedt_temp)
        default = aedt_temp.TempSettings(mode, restore=o.aedt_temp_restore, strict=o.aedt_temp_strict,
                                         template=template)
    except (TypeError, ValueError) as exc:
        raise PlanError(str(exc)) from None
    overrides: dict[int, aedt_temp.TempSettings] = {}
    for item in o.aedt_temp_index:
        text, sep, mode_text = item.partition("=") if isinstance(item, str) else ("", "", "")
        if not sep or not re.fullmatch(r"[0-9]{1,7}", text):
            raise PlanError(f"--aedt-temp-index takes INDEX=MODE, got {item!r}")
        index = int(text)
        if index >= n:
            raise PlanError(f"--aedt-temp-index {item}: the plan has tasks 0 to {n - 1}, no task {index}")
        if index in overrides:
            raise PlanError(f"--aedt-temp-index gives task {index} twice")
        try:
            mode, template = aedt_temp.parse_mode(mode_text)
            overrides[index] = aedt_temp.TempSettings(mode, restore=o.aedt_temp_restore, strict=o.aedt_temp_strict,
                                                      template=template)
        except ValueError as exc:
            raise PlanError(f"--aedt-temp-index {item}: {exc}") from None
    return default, overrides


def _json_object(inputs: acceptance_gate.InputCapture, path: Path, what: str) -> dict:
    try:
        data = inputs.parse(path, json.load)
    except acceptance_gate.InputChanged as exc:
        raise PlanError(f"{what}: {exc}") from None
    except (OSError, ValueError) as exc:
        raise PlanError(f"{what}: {path} is missing or not JSON ({exc})") from None
    if not isinstance(data, dict):
        raise PlanError(f"{what}: {path} is not a JSON object")
    return data


def _source_plan(inputs: acceptance_gate.InputCapture, manifest: Mapping, dataset_id: str,
                 where: str) -> tuple[acceptance_gate.PlanContext, dict]:
    """An included array_run job's own plan: the one its manifest names, with the recorded sha256 and the production
    policy (a probe's output never enters a tree)."""
    block = manifest.get("plan")
    if not (isinstance(block, Mapping) and isinstance(block.get("path"), str) and os.path.isabs(block["path"])
            and isinstance(block.get("sha256"), str) and type(block.get("index")) is int and block["index"] >= 0):
        raise PlanError(f"{where}: an array_run job's manifest names its plan {{path, sha256, index}}, got {block!r}")
    path = Path(block["path"]).resolve()
    plan = _json_object(inputs, path, f"{where}: its plan")
    digest = inputs.recorded(path)
    if digest != block["sha256"]:
        raise PlanError(f"{where}: its manifest names the plan {path} with sha256 {block['sha256']}, but that file's "
                        f"sha256 is {digest}")
    if plan.get("schema") != PLAN_SCHEMA or plan.get("policy") != "production" or plan.get("dataset_id") != dataset_id:
        raise PlanError(f"{where}: its plan {path} is schema {plan.get('schema')!r}, policy {plan.get('policy')!r}, "
                        f"dataset {plan.get('dataset_id')!r}; an included array_run job comes from a production "
                        f"{PLAN_SCHEMA} plan of its own dataset (HF-057)")
    context = acceptance_gate.PlanContext(str(path), digest, block["index"], "production")
    return context, {"path": str(path), "sha256": digest, "index": block["index"]}


def _include_job(value: object, dataset_id: str, expect: Mapping[str, list[float]]) -> dict:
    """One ``--include-job`` checked as the module docstring says, with the digest of every file its acceptance read."""
    try:
        job = Path(os.fspath(value)).resolve(strict=True)                # type: ignore[arg-type]
    except (OSError, RuntimeError, TypeError) as exc:
        raise PlanError(f"--include-job {value} cannot be reached ({exc})") from None
    where = f"--include-job {job}"
    if not job.is_dir():
        raise PlanError(f"{where} is not a job directory")
    inputs = acceptance_gate.InputCapture()
    manifest = _json_object(inputs, job / "manifest.json", where)
    sidecars = sorted(job.glob("*.dataset.json"))
    if len(sidecars) != 1:
        raise PlanError(f"{where} holds {len(sidecars)} dataset sidecars (*.dataset.json), not 1")
    sidecar = _json_object(inputs, sidecars[0], where)
    version = sidecar.get("schema_version")
    if not (isinstance(version, str) and version.split(".")[0] == "2"):
        raise PlanError(f"{where}: its sidecar declares schema_version {version!r}; an included job must be schema 2 "
                        "(S65, one schema per tree), and a probe's schema-1 output never enters a tree")
    job_dataset, frequency = sidecar.get("dataset_id"), sidecar.get("frequency_ghz")
    if not isinstance(job_dataset, str) or not job_dataset:
        raise PlanError(f"{where}: its sidecar names no dataset ID")
    if job_dataset == dataset_id:
        raise PlanError(f"{where}: its dataset ID {job_dataset} is the plan's own; an included job brings another "
                        "dataset ID into the tree (S65)")
    try:
        frequency = float(_positive(frequency, "its sidecar's frequency_ghz"))
    except ValueError as exc:
        raise PlanError(f"{where}: {exc}") from None
    excitation = sidecar.get("excitation") if isinstance(sidecar.get("excitation"), Mapping) else {}
    thetas, phis = excitation.get("incident_theta_deg"), excitation.get("incident_phi_deg")
    if not (isinstance(thetas, list) and isinstance(phis, list) and _same_angles(thetas, expect["theta_deg"])
            and _same_angles(phis, expect["phi_deg"])):
        raise PlanError(f"{where}: its incidence grid theta {thetas} phi {phis} deg is not the expected theta "
                        f"{expect['theta_deg']} phi {expect['phi_deg']} deg (--expect-incidence, OR-036)")
    problems = acceptance_gate.convergence_problems(manifest)
    if problems:
        raise PlanError(f"{where} is not qualified for both polarizations: {'; '.join(problems)}")
    try:
        requirement = acceptance_gate.required_acceptance(job, manifest, sidecar)
    except acceptance_gate.GateError as exc:
        raise PlanError(f"{where}: {exc}") from None
    context, source_plan = None, None
    if requirement.kind == "array_run":
        context, source_plan = _source_plan(inputs, manifest, job_dataset, where)
    problems = acceptance_gate.sidecar_file_problems(job, sidecar, inputs)
    if not problems and requirement.kind in acceptance_gate.REPORT_NAMES:
        name = acceptance_gate.REPORT_NAMES[requirement.kind]
        report = _json_object(inputs, job / name, where)
        try:
            problems = acceptance_gate.report_problems(report, requirement.kind, job, requirement, inputs,
                                                       context=context)
        except acceptance_gate.GateError as exc:
            problems = [str(exc)]
    if problems:
        shown = "; ".join(problems[:3]) + (f"; and {len(problems) - 3} more" if len(problems) > 3 else "")
        raise PlanError(f"{where}: its {requirement.kind} acceptance does not pass the gate: {shown}")
    files = {}
    for key, digest in inputs.hashes.items():
        path = Path(key)
        files[path.relative_to(job).as_posix() if job in path.parents else key] = digest
    return {"job": str(job), "dataset_id": job_dataset, "frequency_ghz": frequency,
            "acceptance_kind": requirement.kind, "files": dict(sorted(files.items())), "source_plan": source_plan}


def _load_config(path: object, probe: bool) -> tuple[Path, bytes, dict, RunConfig]:
    """The config file's path, bytes, parsed TOML and config: a job-array config of the plan's policy that meets
    BBR025 (HF-057)."""
    config_path = _existing_file(path, "--config")
    config_bytes = config_path.read_bytes()
    try:
        raw = tomllib.loads(config_bytes.decode("utf-8"))
        base = build_config(raw)
    except (UnicodeDecodeError, tomllib.TOMLDecodeError, ConfigError) as exc:
        raise PlanError(f"the config {config_path} cannot be loaded: {exc}") from None
    if base.solver.frequency_ghz is not None:
        raise PlanError(f"the config {config_path} gives solver.frequency_ghz = {base.solver.frequency_ghz!r}; a "
                        "job-array config gives none, since the plan supplies each task's frequency (HF-057)")
    if base.acceptance.kind != ACCEPTANCE_KIND:
        raise PlanError(f"the config {config_path} declares [acceptance] kind = {base.acceptance.kind!r}; the job "
                        f"array runs kind = {ACCEPTANCE_KIND!r} only, each task judged against its frozen plan, and "
                        "refuses a reference-based config (HF-057)")
    geometry, problems = base.geometry, []
    missing = [s for s in BBR025_SYMMETRY if s not in geometry.symmetry]
    if missing:
        problems.append(f"geometry.symmetry {list(geometry.symmetry)} lacks {missing}")
    if geometry.pose != "canonical-z":
        problems.append(f"geometry.pose is {geometry.pose!r}, not 'canonical-z'")
    axes = tuple(tuple(float(c) for c in v) for v in (geometry.exit_cs_x, geometry.exit_cs_y))
    if axes != CANONICAL_EXIT_CS:
        problems.append(f"exit_cs_x, exit_cs_y are {[list(v) for v in axes]}, not "
                        f"{[list(v) for v in CANONICAL_EXIT_CS]}")
    if problems:
        raise PlanError(f"the config {config_path} breaks BBR025, the rules BBRsim's sidecar loader enforces: "
                        + "; ".join(problems))
    output = base.output
    if probe:
        if output.schema != 1:
            raise PlanError(f"a probe plan runs the single-angle diagnostic, which writes schema 1 (HF-057); the "
                            f"config {config_path} writes schema {output.schema}")
        if tuple(output.polarizations) != acceptance_gate.ARRAY_POLICIES["probe"]:
            raise PlanError(f"a probe plan exports Ephi=0 only; the config's output.polarizations is "
                            f"{list(output.polarizations)}")
    else:
        if output.schema != 2:
            raise PlanError(f"a production plan writes schema 2 (HF-057); the config {config_path} writes schema "
                            f"{output.schema}, which only a --probe plan runs")
        if sorted(output.polarizations) != [0, 1]:
            raise PlanError(f"a production plan exports both polarizations; the config's output.polarizations is "
                            f"{list(output.polarizations)}")
    return config_path, config_bytes, raw, base


def _incidence(base: RunConfig, expect_text: str | None, probe: bool) -> tuple[dict, dict]:
    """The config's incidence grid and the expected one: OR-036's (or --expect-incidence's) for production, its one
    key for a probe."""
    from bbsim.sampling import incident_angle_values

    exc = base.excitation
    incidence = {"theta_deg": [float(v) for v in incident_angle_values(exc.theta_lower_deg, exc.theta_upper_deg,
                                                                       exc.theta_step_deg)],
                 "phi_deg": [float(v) for v in incident_angle_values(exc.phi_lower_deg, exc.phi_upper_deg,
                                                                     exc.phi_step_deg)]}
    if probe:
        if len(incidence["theta_deg"]) != 1 or len(incidence["phi_deg"]) != 1:
            raise PlanError(f"a probe plan runs one incidence key; the config declares theta {incidence['theta_deg']} "
                            f"phi {incidence['phi_deg']} deg")
        expect = incidence
    else:
        expect = _expected_incidence(expect_text)
        if not (_same_angles(incidence["theta_deg"], expect["theta_deg"])
                and _same_angles(incidence["phi_deg"], expect["phi_deg"])):
            raise PlanError(f"the config's incidence grid theta {incidence['theta_deg']} phi {incidence['phi_deg']} "
                            f"deg is not the expected theta {expect['theta_deg']} phi {expect['phi_deg']} deg "
                            "(--expect-incidence, OR-036)")
    return incidence, expect


def _reference(path: object) -> tuple[Path, bytes, dict]:
    """The HF-041 column-order reference a production plan freezes: its path, bytes and parsed object."""
    if path is None:
        raise PlanError("a production plan writes schema 2 and needs --column-order-reference, the HF-041 "
                        "reference crack1's accepted 500 GHz solve keeps as scratch/column_order_reference.json")
    reference_path = _existing_file(path, "--column-order-reference")
    reference_bytes = reference_path.read_bytes()
    try:
        reference = json.loads(reference_bytes.decode("utf-8"))
    except (UnicodeDecodeError, ValueError) as err:
        raise PlanError(f"the column-order reference {reference_path} is not JSON ({err})") from None
    if not isinstance(reference, dict):
        raise PlanError(f"the column-order reference {reference_path} is not a JSON object")
    return reference_path, reference_bytes, reference


def _project_record(project: Path, cfg: RunConfig) -> tuple[str, Path, bytes]:
    """The record beside the prepared project (``.mapping.json``, which prepare_project writes, or ``.build.json``)
    read by the runner's own ``mapping.recorded_dataset_id``: it must give the config's dataset_id for its design and
    object, else the plan is refused, where a single run only warns (HF-057). Returns the record's name beside the
    project, its resolved path and its bytes, which the plan freezes beside the project copy so that each task's
    ``run_frequency.stage_project`` confirms the dataset again."""
    from bbsim.mapping import recorded_dataset_id

    design, obj, dataset = cfg.project.design, cfg.geometry.object, cfg.project.dataset_id
    try:
        recorded, record = recorded_dataset_id(project, design, obj)
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
        raise PlanError(f"the record beside the project {project} (.mapping.json or .build.json) cannot be read "
                        f"({type(exc).__name__}: {exc})") from None
    if record is None:
        raise PlanError(f"the project {project} has no .mapping.json or .build.json beside it, so the config's "
                        f"dataset_id {dataset!r} cannot be confirmed; a plan runs a prepared project, whose record "
                        "prepare_project writes (HF-057)")
    if recorded is None:
        raise PlanError(f"{record} has no entry for {design}/{obj}, so the config's dataset_id {dataset!r} cannot be "
                        "confirmed (HF-057)")
    if recorded != dataset:
        raise PlanError(f"the config's dataset_id {dataset!r} disagrees with {record}, which gives {recorded!r} for "
                        f"design {design!r} object {obj!r}")
    try:
        return Path(record).name, Path(record).resolve(strict=True), Path(record).read_bytes()
    except (OSError, RuntimeError) as exc:
        raise PlanError(f"the record {record} cannot be read ({exc})") from None


def _check_cores(guard: Any, cores: int) -> None:
    """Every pending or running array the guard manages (comment ``bbsim:``) must run at the plan's cores per task:
    two arrays at different core counts share the per-user ``pyaedt_config`` registry entry, so one's tasks could solve
    at the other's count (HF-057). Its CPUs/Task comes from ``scontrol show job`` (the comment names no plan file to
    read); a value that cannot be read refuses too."""
    try:
        entries = guard.scheduler.queue()
    except licence_guard.GuardError as err:
        raise PlanError(f"the scheduler's queue cannot be read, so the cores of running bbsim arrays are unknown: "
                        f"{err}") from None
    parents = sorted({e.array_parent for e in entries if e.array_parent and e.state in _CORE_STATES
                      and e.comment.startswith(licence_guard.MANAGED_COMMENT_PREFIX)})
    for parent in parents:
        try:
            found = guard.scheduler.cpus_per_task(parent)
        except licence_guard.GuardError as err:
            raise PlanError(f"the bbsim array {parent}'s cores per task cannot be read ({err}), so --cores {cores} "
                            "cannot be checked against it (HF-057)") from None
        if isinstance(found, bool) or not isinstance(found, int):
            raise PlanError(f"the bbsim array {parent}'s cores per task cannot be read: scontrol show job {parent} "
                            f"prints no CPUs/Task, so --cores {cores} cannot be checked against it (HF-057)")
        if found != cores:
            raise PlanError(f"the bbsim array {parent} runs at {found} cores per task and this plan at --cores "
                            f"{cores}: two arrays at different core counts share the per-user pyaedt_config registry "
                            f"entry, so one's tasks could solve at the other's (HF-057); plan with --cores {found}, or "
                            f"once {parent} has left the queue")


def _check_licence(guard: Any, throttle: int, cap: int) -> None:
    """The guard's reconciled usage plus the plan's throttle must fit min(CEILING, cap); advisory at plan time, since
    admission is repeated at every submission (OR-045)."""
    try:
        usage = guard.usage()
    except licence_guard.GuardError as err:
        raise PlanError(f"the licence guard's usage cannot be read, so the plan cannot be checked against the cap: "
                        f"{err}") from None
    limit = min(licence_guard.CEILING, cap)
    if usage.total + throttle > limit:
        raise PlanError(f"the guard counts {usage.total} in use, and this plan's throttle {throttle} exceeds {limit}, "
                        f"the smaller of the ceiling {licence_guard.CEILING} (OR-045) and --licence-cap {cap}; the "
                        "count is advisory at plan time and admission is repeated at every submission, but a plan "
                        "that does not fit now is refused")
    log.info("licence guard %s: %d in use + this plan's throttle %d fits %d (advisory at plan time: admission is "
             "repeated at every submission, OR-045)", guard.root, usage.total, throttle, limit)
    if usage.uncertain:
        log.warning("licence guard: uncertain reservations %s; no slot is admitted at submission until each is "
                    "resolved", ", ".join(usage.uncertain))


@dataclass
class _Draft:
    """A checked plan, nothing created yet, and what ``_create`` needs to make it."""

    plan: dict
    run_dir: Path
    source_repo: Path
    commit: str
    pyaedt_version: str
    config_bytes: bytes
    project_source: Path
    reference: dict | None
    reference_bytes: bytes | None
    record_bytes: bytes
    design_object: tuple[str, str]


def _draft(o: PlanOptions, env: EnvironmentState, guard: Any, now: datetime, repo: Path) -> _Draft:
    if not isinstance(o, PlanOptions):
        raise PlanError(f"build_plan takes PlanOptions, got {type(o).__name__}")
    if not isinstance(now, datetime) or now.utcoffset() is None:
        raise PlanError(f"the plan's creation time must carry a time zone, got {now!r}")
    if not isinstance(o.probe, bool):
        raise PlanError(f"probe must be true or false, got {o.probe!r}")
    probe = o.probe
    policy = "probe" if probe else "production"
    run_dir = _new_path(o.run_dir, "--run-dir")
    if not _NAME.fullmatch(run_dir.name):
        raise PlanError(f"the run directory's name {run_dir.name!r} names the Slurm jobs (bbsim_<name>), so it holds "
                        "letters, digits, '.', '_' and '-' only and starts with a letter or digit")
    tree = None
    if probe:
        for given, option, why in ((o.tree, "--tree", "its schema-1 diagnostic never enters a tree"),
                                   (o.include_jobs, "--include-job", "it assembles no tree"),
                                   (o.column_order_reference, "--column-order-reference",
                                    "the schema-1 probe takes no column-order reference"),
                                   (o.expect_incidence, "--expect-incidence", "its one key is its config's")):
            if given:
                raise PlanError(f"a probe plan takes no {option}: {why} (HF-057)")
    elif o.tree is None:
        raise PlanError("a production plan needs --tree, the data tree its runs are assembled into")
    else:
        tree = _new_path(o.tree, "--tree")
        if tree == run_dir:
            raise PlanError("--tree and --run-dir name the same directory")

    if not isinstance(env, EnvironmentState) or not isinstance(env.commit, str) or not _COMMIT.fullmatch(env.commit):
        raise PlanError(f"the code checkout's state must be an EnvironmentState with a commit, got {env!r}")
    if env.dirty:
        raise PlanError(f"the code checkout {repo} has uncommitted changes to tracked files; commit or restore them, "
                        "so every task runs the code the plan records")
    if not isinstance(env.pyaedt_version, str) or not env.pyaedt_version.strip():
        raise PlanError(f"the installed PyAEDT version must be known, got {env.pyaedt_version!r}")
    cap = _integer(o.licence_cap, "--licence-cap (at most 10 HFSS solves at a time, OR-045)", 1, licence_guard.CEILING)
    sbatch_args = list(o.sbatch_args)
    try:
        licence_guard.check_user_sbatch_args(sbatch_args)
    except licence_guard.GuardRefused as exc:
        raise PlanError(f"--sbatch-arg: {exc}") from None
    named = {token.split("=", 1)[0] for token in sbatch_args if token.startswith("--")}
    twice = [name for name in PLAN_RESOURCE_OPTIONS if name in named]
    if twice:
        uses = ", ".join(f"{n}, which the plan sets itself: use {PLAN_RESOURCE_OPTIONS[n]}" for n in twice)
        raise PlanError(f"--sbatch-arg sets {uses}; sbatch would take one of the two values silently (HF-057)")
    try:
        guard_dir = Path(os.fspath(o.guard_dir)).resolve(strict=True)
    except (OSError, RuntimeError, TypeError) as exc:
        raise PlanError(f"--guard-dir {o.guard_dir} cannot be reached ({exc})") from None
    if guard_dir != guard.root:
        raise PlanError(f"--guard-dir {guard_dir} is not the guard's root {guard.root}; one guard holds the shared "
                        "cap (OR-045)")
    inherited = os.environ.get(licence_guard.GUARD_ENV)
    if inherited and Path(inherited).resolve() != guard.root:
        raise PlanError(f"the inherited {licence_guard.GUARD_ENV}={inherited} names another guard than {guard.root}; "
                        "one guard holds the shared cap (OR-045)")
    if o.order not in ORDERS:
        raise PlanError(f"--order must be one of {', '.join(ORDERS)}, got {o.order!r}")
    if o.tmpdir not in TMPDIRS:
        raise PlanError(f"--tmpdir must be one of {', '.join(TMPDIRS)}, got {o.tmpdir!r}")
    stagger = _stagger(o.stagger)
    try:
        launch_timeout_setting(o.launch_timeout_s)
    except ValueError as exc:
        raise PlanError(f"--launch-timeout-s: {exc}") from None
    max_array_size = _integer(o.max_array_size, "--max-array-size", 1)
    monitor_interval = _integer(o.monitor_interval_s, "--monitor-interval-s", 1)
    assembly = {"mem": _mem(o.assembly_mem, "--assembly-mem"), "time": o.assembly_time,
                "qos": None if o.assembly_qos is None else _qos(o.assembly_qos, "--assembly-qos")}
    _seconds(o.assembly_time, "--assembly-time")
    if o.monitor is not None and not isinstance(o.monitor, bool):
        raise PlanError(f"monitor must be true, false or None, got {o.monitor!r}")
    monitor = probe if o.monitor is None else o.monitor

    grid_path = _existing_file(o.grid, "--grid")
    try:
        entries, grid_sha256 = frequency_grid.read_grid(grid_path, max_array_size)
    except (OSError, frequency_grid.GridError) as exc:
        raise PlanError(f"--grid: {exc}") from None

    config_path, config_bytes, raw, base = _load_config(o.config, probe)
    cores = base.solver.cores if o.cores is None else _integer(o.cores, "--cores", 1)

    incidence, expect = _incidence(base, o.expect_incidence, probe)
    project_source = _existing_file(o.project if o.project is not None else base.project.path, "the project")
    record_name, record_source, record_bytes = _project_record(project_source, base)
    reference_path, reference_bytes, reference = (None, None, None) if probe else \
        _reference(o.column_order_reference)
    names = [p.name for p in (config_path, project_source, reference_path) if p is not None] + [record_name]
    if len(set(names)) != len(names):
        raise PlanError(f"the inputs copied into <RUN>/inputs/ under their own names share a name: {names}")

    bridge = _bridge(o.scripts_dir, [ARRAY_SCRIPT, ASSEMBLY_SCRIPT] + ([MONITOR_SCRIPT] if monitor else []), repo)
    try:
        source_repo = Path(repo).resolve(strict=True)
    except (OSError, RuntimeError) as err:
        raise PlanError(f"the code checkout {repo} cannot be reached ({err})") from None
    if _git(source_repo, "cat-file", "-e", f"{env.commit}^{{commit}}").returncode != 0:
        raise PlanError(f"the code checkout {source_repo} holds no commit {env.commit}")

    ordered = sorted(entries, key=lambda e: e.frequency_ghz, reverse=o.order == "longest-first")
    n = len(ordered)
    tasks = []
    for index, entry in enumerate(ordered):
        try:
            effective = frequency_grid.effective_l_resolution_mm(entry, base.exit_field.resolution_mm[1])
        except frequency_grid.GridError as err:
            raise PlanError(f"--grid {grid_path}: {err}") from None
        tasks.append({"index": index, "frequency_ghz": entry.frequency_ghz,
                      "label": frequency_label(entry.frequency_ghz), "class": None,
                      "l_resolution_mm": entry.l_resolution_mm, "effective_l_resolution_mm": effective,
                      "l_step_mm": None, "lattice": None, "aedt_temp": None})
    classes, throttle = _classes(o, n, DEFAULT_THROTTLE[policy])
    for entry in classes:
        for i in entry["indices"]:
            tasks[i]["class"] = entry["name"]
    default_temp, temp_overrides = _temp_settings(o, n)
    for task in tasks:
        task["aedt_temp"] = _temp_record(temp_overrides.get(task["index"], default_temp))
    if o.monitor_time is not None:
        monitor_time = o.monitor_time
        _seconds(monitor_time, "--monitor-time")
    else:
        duration = max(math.ceil(len(c["indices"]) / c["throttle"]) * _seconds(c["time"], "time") for c in classes)
        monitor_time = _time_text(duration + MONITOR_MARGIN_S)

    inputs = run_dir / "inputs"
    plan: dict[str, Any] = {
        "schema": PLAN_SCHEMA, "created_utc": now.astimezone(timezone.utc).isoformat(), "name": run_dir.name,
        "run_dir": str(run_dir), "tree": None if tree is None else str(tree), "dataset_id": base.project.dataset_id,
        "probe": probe, "policy": policy, "acceptance_kind": ACCEPTANCE_KIND,
        "code": {"repo": str(run_dir / "code"), "source_repo": str(source_repo), "commit": env.commit,
                 "pyaedt_version": env.pyaedt_version},
        "bridge": bridge,
        "config": {"source": str(config_path), "copy": str(inputs / config_path.name),
                   "sha256": hashlib.sha256(config_bytes).hexdigest()},
        "project": {"source": str(project_source), "copy": str(inputs / project_source.name),
                    "sha256": sha256_file(project_source),
                    "record": {"source": str(record_source), "copy": str(inputs / record_name),
                               "sha256": hashlib.sha256(record_bytes).hexdigest()}},
        "column_order_reference": None if reference_path is None else {
            "source": str(reference_path), "copy": str(inputs / reference_path.name),
            "sha256": hashlib.sha256(reference_bytes).hexdigest()},
        "grid": {"source": str(grid_path), "sha256": grid_sha256,
                 "entries": [{"frequency_ghz": e.frequency_ghz, "l_resolution_mm": e.l_resolution_mm}
                             for e in entries]},
        "order": o.order, "tasks": tasks, "classes": classes,
        "slurm": {"cores": cores, "throttle": throttle, "stagger": stagger, "sbatch_args": sbatch_args,
                  "tmpdir": o.tmpdir, "job_name": f"bbsim_{run_dir.name}", "comment": f"bbsim:{run_dir / 'plan.json'}",
                  "max_array_size": max_array_size, "assembly": assembly,
                  "monitor": {"enabled": monitor, "interval_s": monitor_interval, "time": monitor_time}},
        "incidence": incidence, "expect_incidence": expect,
        "aedt": {"temp": _temp_record(default_temp), "launch_timeout_s": float(o.launch_timeout_s)},
        "include_jobs": [], "licence_cap": cap, "guard_dir": str(guard.root), "jobs_root": str(run_dir / "jobs"),
    }

    from bbsim.run_frequency import check_before_aedt

    for task in tasks:
        what = f"task {task['index']} ({task['label']})"
        try:
            cfg = task_config(plan, task["index"], raw_config=raw)
        except ValueError as err:
            raise PlanError(f"{what}: its config cannot be built from the frozen config: {err}") from None
        lattice = predict_exit_lattice(cfg, task["effective_l_resolution_mm"])
        limit = frequency_grid.max_l_step_mm(task["frequency_ghz"])
        if lattice["step_mm"][1] > limit:
            raise PlanError(f"{what}: the fitted step along l {lattice['step_mm'][1]!r} mm exceeds c/(4f) = "
                            f"{limit:.7g} mm (OR-048)")
        task["lattice"], task["l_step_mm"] = lattice, lattice["step_mm"][1]
        try:
            check_before_aedt(cfg, reference)
        except ValueError as err:
            raise PlanError(f"{what}: the runner's pre-AEDT check refuses it: {err}") from None

    seen_jobs: set[str] = set()
    seen_keys: dict[tuple[str, float], str] = {}
    for value in o.include_jobs:
        record = _include_job(value, base.project.dataset_id, expect)
        key = (record["dataset_id"], record["frequency_ghz"])
        if record["job"] in seen_jobs:
            raise PlanError(f"--include-job {record['job']} is given twice")
        if key in seen_keys:
            raise PlanError(f"--include-job {record['job']} and {seen_keys[key]} are both {key[0]} at {key[1]:g} GHz")
        seen_jobs.add(record["job"])
        seen_keys[key] = record["job"]
        plan["include_jobs"].append(record)

    _check_cores(guard, cores)
    _check_licence(guard, throttle, cap)
    return _Draft(plan=plan, run_dir=run_dir, source_repo=source_repo, commit=env.commit,
                  pyaedt_version=env.pyaedt_version, config_bytes=config_bytes, project_source=project_source,
                  reference=reference, reference_bytes=reference_bytes, record_bytes=record_bytes,
                  design_object=(base.project.design, base.geometry.object))


# --- making the plan --------------------------------------------------------------------------------------------------

def _write_new(path: Path, data: bytes) -> None:
    with open(path, "xb") as fh:
        fh.write(data)


def _cleanup(d: _Draft, top: Path | None, worktree: bool) -> list[str]:
    """Remove the worktree and ``top``, the highest directory this plan made (None: it made none); what is left."""
    code = d.run_dir / "code"
    left: list[str] = []
    try:
        if worktree and os.path.lexists(code):
            _git(d.source_repo, "worktree", "remove", "--force", str(code))
        if top is not None:
            shutil.rmtree(top, ignore_errors=True)
            if os.path.lexists(top):
                left.append(str(top))
        if worktree and code.resolve() in _worktree_paths(d.source_repo):
            left.append(f"the worktree entry of {code} in {d.source_repo} (git -C {d.source_repo} worktree remove "
                        f"--force {code})")
    except PlanError as exc:
        left.append(f"the cleanup could not be checked ({exc})")
    return left


def _create(d: _Draft) -> Path:
    """Make the checked plan: the run directory, the worktree, the input copies, plan.json and plan.sha256 last. Any
    failure, an interruption included, removes everything made first."""
    plan, run_dir = d.plan, d.run_dir
    missing = []
    directory = run_dir
    while not os.path.lexists(directory):
        missing.append(directory)
        directory = directory.parent
    if not missing:
        raise PlanError(f"--run-dir {run_dir} exists now, though it did not when the plan was checked; nothing is "
                        "written into it")
    code = run_dir / "code"
    top: Path | None = None
    worktree = False
    step = "creating the run directory"
    try:
        for directory in reversed(missing):     # one at a time: a directory another process made is never ours
            directory.mkdir()
            top = top or directory
        for name in JOB_DIRS:
            (run_dir / name).mkdir()
        step, worktree = "adding the worktree", True
        _git_out(d.source_repo, "worktree", "add", "--detach", str(code), d.commit)
        head = _git_out(code, "rev-parse", "--verify", "HEAD").strip()
        if head != d.commit or _git_out(code, "status", "--porcelain", "--untracked-files=no").strip():
            raise PlanError(f"the worktree {code} is not clean at the commit {d.commit} (its HEAD is {head})")
        if d.reference is not None:
            from bbsim import dataset_v2

            step = f"checking the column-order reference against the worktree {code}"
            dataset_v2.check_reference_before_aedt(d.reference, d.pyaedt_version, code)
            _write_new(Path(plan["column_order_reference"]["copy"]), d.reference_bytes)
        step = "copying the inputs"
        _write_new(Path(plan["config"]["copy"]), d.config_bytes)
        copy = Path(plan["project"]["copy"])
        shutil.copyfile(d.project_source, copy)
        if sha256_file(copy) != plan["project"]["sha256"]:
            raise PlanError(f"the project {d.project_source} changed while it was copied")
        _write_new(Path(plan["project"]["record"]["copy"]), d.record_bytes)
        from bbsim.mapping import recorded_dataset_id

        confirmed, record = recorded_dataset_id(copy, *d.design_object)
        if confirmed != plan["dataset_id"]:
            raise PlanError(f"the frozen record {record} gives {confirmed!r}, not the plan's dataset "
                            f"{plan['dataset_id']!r}, for {'/'.join(d.design_object)}: it changed while it was read")
        step = "writing plan.json"
        data = (json.dumps(plan, indent=2, sort_keys=True) + "\n").encode()
        tmp = run_dir / ".plan.json.tmp"
        _write_new(tmp, data)
        os.replace(tmp, run_dir / "plan.json")
        _write_new(run_dir / "plan.sha256", f"{hashlib.sha256(data).hexdigest()}  plan.json\n".encode())
    except BaseException as exc:
        left = _cleanup(d, top, worktree)
        if left:
            raise PlanLeftovers(f"{step} failed ({_one_line(exc)}), and the cleanup left {'; '.join(left)}; remove "
                                "them by hand") from exc
        if isinstance(exc, PlanError):
            raise
        if isinstance(exc, (OSError, ValueError)):
            raise PlanError(f"{step}: {_one_line(exc)}; everything the planner made is removed") from exc
        raise
    return run_dir / "plan.json"


def build_plan(options: PlanOptions, *, env: EnvironmentState, guard: Any, now: datetime,
               repo: Path = REPO_ROOT) -> Path:
    """Check ``options`` (the module docstring's list) and, only if every check passes, make the plan: the run
    directory, the worktree of ``repo`` (the code checkout, this one by default) at ``env.commit``, the input copies,
    and ``plan.json`` and ``plan.sha256`` last. ``env`` is ``manifest.environment_state(repo)``; ``guard`` a
    ``licence_guard.Guard`` (its ``root`` and ``usage()``); ``now`` an aware time. Returns ``<RUN>/plan.json``.
    ``PlanError`` refuses with nothing left; ``PlanLeftovers`` names what a failed cleanup left."""
    return _create(_draft(options, env, guard, now, Path(repo)))


# --- the submission block --------------------------------------------------------------------------------------------

def submission_block(plan_path: Path, plan_sha256: str, plan: Mapping[str, Any], *,
                     admitted: Mapping[str, str] | None = None, monitor_submitted: bool = False,
                     assembly_submitted: bool = False) -> str:
    """The block the planner prints (OR-044, OR-045, HF-057), or what remains of it once some of its jobs are in
    ``submissions.log`` (``python -m bbsim.array_submit remaining``): ``admitted`` maps a class already submitted to
    its job id, written as a literal id, and the monitor and the assembly already submitted are left out.

    From the frozen code every job goes through ``python -m bbsim.array_submit``, which captures the recorded bridge
    blob and repeats admission and the cores rule when its line runs. The lines run in one subshell, so a STOP leaves
    the operator's shell open. A STOP after a refusal (exit 2: nothing was reserved or submitted) prints what remains,
    ready to paste, with the ids admitted so far; a STOP after any other failure (an uncertain or unrecorded result)
    names the command that prints it, to run once the failure's own message is dealt with."""
    admitted = dict(admitted or {})
    q = shlex.quote
    code, plan_arg = q(str(plan["code"]["repo"])), q(str(plan_path))
    submit = f"python -m bbsim.array_submit submit {plan_arg} {plan_sha256}"
    remaining = f"python -m bbsim.array_submit remaining {plan_arg} {plan_sha256}"
    names = [entry["name"] for entry in plan["classes"]]
    ids = {name: admitted.get(name, f"${name}") for name in names}
    monitor = plan["slurm"]["monitor"]["enabled"] is True and not monitor_submitted
    lines = ["# Submit through the licence guard (OR-045) from the frozen code, in a shell that has sourced bb_env.sh;",
             "# each line repeats admission when it runs, so a delayed paste is checked again.",
             f"cd {code} && (",
             'bbsim_stop() { echo "STOP: $1 was not submitted (exit $2); the arrays submitted before it keep '
             'running."; if [ "$2" = 2 ]; then echo "Nothing was reserved or submitted for it. Once the guard '
             f'admits it, paste what remains:"; {remaining}; else echo "Its result is uncertain or unrecorded: do '
             f'what the message above says first, then print what remains: cd {code} && {remaining}"; fi; exit 1; }}']
    for name in names:
        if name not in admitted:
            lines.append(f'{name}=$({submit} array --class {name}) || bbsim_stop "class {name}" $?')
    if monitor:
        lines.append(f'M=$({submit} monitor) || echo "the licence monitor was not submitted (exit $?); the arrays '
                     f'keep running; {remaining} prints its line"')
    if not assembly_submitted:
        dependency = "afterany:" + ":".join(ids[name] for name in names)
        lines.append(f'R=$({submit} assembly --dependency "{dependency}") || bbsim_stop "the assembly" $?')
    jobs = " ".join(f"{name}={ids[name]}" for name in names) + (" monitor=$M" if monitor else "") + \
        ("" if assembly_submitted else " assembly=$R")
    lines += [f'echo "submitted: {jobs}"', ")",
              "# Rerun failed indices of one class (the assembly report prints these lines): scancel the pending",
              "# assembly that submissions.log names, then submit the rerun and a new assembly after it:",
              f"#   cd {code} && python -m bbsim.array_submit rerun {plan_arg} {plan_sha256} --indices I,J "
              "[--mem MEM] [--time TIME] [--qos QOS]",
              f"#   cd {code} && {submit} assembly --dependency afterany:RERUN_ID",
              "# Live throttle changes go through the guard from the frozen code, never through scontrol;",
              f"# submissions.log records one only through python -m bbsim.array_submit throttle {plan_arg} "
              f"{plan_sha256} ARRAY_ID THROTTLE:",
              f"#   cd {code} && python -m bbsim.licence_guard throttle ARRAY_ID THROTTLE --cap {plan['licence_cap']}"]
    return "\n".join(lines)


# --- the CLI ----------------------------------------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="plan_hfss_array.py", allow_abbrev=False,
                                description="Freeze a job-array run plan (HF-057), with no AEDT, and print its "
                                            "guarded submission block (OR-045). plan_hfss_array.py --log-submission "
                                            "RUN KIND ID records a job submitted by hand, with no admission.")
    p.add_argument("--config", type=Path, required=True,
                   help="the job array's config, with no frequency and [acceptance] kind = \"array_run\" "
                        "(configs/crack1_grid.toml; crack1_grid_single_angle.toml with --probe)")
    p.add_argument("--grid", type=Path, required=True,
                   help="the grid file (configs/frequency_grids/): '<GHz> [<requested step along l, mm>]' per line")
    p.add_argument("--run-dir", type=Path, required=True, help="the run directory; it must not exist, and its name "
                                                               "names the jobs (bbsim_<name>)")
    p.add_argument("--tree", type=Path, help="the data tree the runs are assembled into (must not exist); required "
                                             "unless --probe, refused with it")
    p.add_argument("--project", type=Path, help="the prepared .aedt project (default the config's project.path)")
    p.add_argument("--cores", type=int, help="cores per task (default the config's solver.cores)")
    p.add_argument("--mem", default="32G", help="memory per task (default %(default)s)")
    p.add_argument("--time", default="08:00:00", help="time limit per task (default %(default)s)")
    p.add_argument("--qos", default="normal", help="QOS of the tasks (default %(default)s)")
    p.add_argument("--throttle", type=int,
                   help="the plan's throttle L, the most tasks solving at once (default 3 with --probe, 2 otherwise); "
                        "clamped to the number of tasks")
    p.add_argument("--class", dest="classes", action="append", nargs="+", default=[], metavar="KEY=VALUE",
                   help="a resource class: indices=0-11,13 [mem=32G] [time=08:00:00] [qos=normal] [throttle=2]; "
                        "classes partition the task indices, and with several each gives a throttle, summing to L")
    p.add_argument("--order", choices=ORDERS, default="longest-first",
                   help="task order: longest-first puts the highest frequency at index 0 (default)")
    p.add_argument("--max-array-size", type=int, default=frequency_grid.MAX_ENTRIES,
                   help="the cluster's MaxArraySize, the most grid entries a plan takes (default %(default)s, Slurm's "
                        "default; the cluster probe records the actual one)")
    p.add_argument("--stagger", default="5:8", metavar="STEP[:MODULUS]",
                   help="each task sleeps (index mod MODULUS) x STEP seconds before it starts (default %(default)s)")
    p.add_argument("--aedt-temp", default="task", metavar="{off,task,env,TEMPLATE}",
                   help="AEDT's temp directory in each task (HF-057): task (default) points it at <job>/aedt_temp and "
                        "sets it to --aedt-temp-restore at the end; a path holding {job} is task mode there, {job} "
                        "expanding to the absolute job directory, so a /scratch/{job} template mirrors the whole "
                        "/resnick/... path of the job under /scratch; env exports TMPDIR and TMP = <job>/aedt_temp "
                        "before AEDT starts; off leaves it")
    p.add_argument("--aedt-temp-index", action="append", default=[], metavar="INDEX=MODE",
                   help="one task's own mode, as --aedt-temp takes it, such as 3=env (repeatable)")
    p.add_argument("--aedt-temp-restore", metavar="PATH",
                   help="task mode: the TempDirectory AEDT is set back to at the end, the host's TempDirectory before "
                        "any plan-mode run, read from AEDT's .cfg by the cluster probe")
    p.add_argument("--aedt-temp-strict", action="store_true",
                   help="task mode: a task fails instead of warning when AEDT does not report the requested directory")
    p.add_argument("--launch-timeout-s", type=float, default=DEFAULT_LAUNCH_TIMEOUT_S,
                   help="seconds PyAEDT waits for AEDT to start (default %(default)s)")
    p.add_argument("--tmpdir", choices=TMPDIRS, default="run",
                   help="each task's TMPDIR: run (<RUN>/slurm/tmp_<A>_<i>, default), node (Slurm's) or off")
    p.add_argument("--expect-incidence", metavar="THETAS:PHIS",
                   help="the incidence grid the config must declare, two comma lists of degrees (default OR-036's "
                        "theta 90 to 180 by 15 and phi 0 to 90 by 22.5); not with --probe")
    p.add_argument("--licence-cap", type=int, default=licence_guard.CEILING,
                   help="this plan's own cap on HFSS solves at once, 1 to 10 (OR-045); it never raises the ceiling")
    p.add_argument("--include-job", dest="include_jobs", action="append", default=[], type=Path, metavar="JOBDIR",
                   help="an accepted schema-2 job of another dataset ID to assemble beside the runs (S65); repeatable")
    p.add_argument("--column-order-reference", type=Path, metavar="JSON",
                   help="schema 2: the HF-041 reference crack1's accepted 500 GHz solve keeps as "
                        "scratch/column_order_reference.json; copied and hashed into the plan")
    p.add_argument("--scripts-dir", type=Path,
                   help="scripts/hpc/hfss of the code checkout's sibling bridge clone ../bbridge, the only bridge a "
                        "plan takes (OR-044); the default")
    p.add_argument("--guard-dir", type=Path, default=licence_guard.DEFAULT_ROOT,
                   help="the shared licence guard (default %(default)s)")
    p.add_argument("--sbatch-arg", dest="sbatch_args", action="append", default=[], metavar="OPT",
                   help="an extra sbatch option from the guard's allowlist, given as --sbatch-arg=--nodelist=NODE "
                        "(repeatable)")
    monitor = p.add_mutually_exclusive_group()
    monitor.add_argument("--monitor-licences", dest="monitor", action="store_true",
                         help="run the licence monitor (the default with --probe)")
    monitor.add_argument("--no-monitor-licences", dest="monitor", action="store_false",
                         help="run no licence monitor (the default without --probe)")
    p.set_defaults(monitor=None)
    p.add_argument("--monitor-interval-s", type=int, default=60, help="seconds between lmstat samples (default 60)")
    p.add_argument("--monitor-time", metavar="TIME",
                   help="the monitor's time limit (default the plan's duration estimate plus 30 min)")
    p.add_argument("--assembly-mem", default="4G", help="the assembly job's memory (default %(default)s)")
    p.add_argument("--assembly-time", default="02:00:00", help="the assembly job's time limit (default %(default)s)")
    p.add_argument("--assembly-qos", help="the assembly job's QOS (default the cluster's)")
    p.add_argument("--probe", action="store_true",
                   help="the schema-1 single-angle diagnostic (HF-057, HF-058): no tree, never assembled")
    p.add_argument("--submit", action="store_true",
                   help="after writing the plan, submit its jobs through the guard as the printed block does")
    return p


def options_from_args(args: argparse.Namespace) -> PlanOptions:
    return PlanOptions(
        config=args.config, grid=args.grid, run_dir=args.run_dir, tree=args.tree, project=args.project,
        cores=args.cores, mem=args.mem, time=args.time, qos=args.qos, throttle=args.throttle,
        classes=tuple(tuple(tokens) for tokens in args.classes), order=args.order,
        max_array_size=args.max_array_size, stagger=args.stagger, aedt_temp=args.aedt_temp,
        aedt_temp_index=tuple(args.aedt_temp_index), aedt_temp_restore=args.aedt_temp_restore,
        aedt_temp_strict=args.aedt_temp_strict, launch_timeout_s=args.launch_timeout_s, tmpdir=args.tmpdir,
        expect_incidence=args.expect_incidence, licence_cap=args.licence_cap, include_jobs=tuple(args.include_jobs),
        column_order_reference=args.column_order_reference, scripts_dir=args.scripts_dir, guard_dir=args.guard_dir,
        sbatch_args=tuple(args.sbatch_args), monitor=args.monitor, monitor_interval_s=args.monitor_interval_s,
        monitor_time=args.monitor_time, assembly_mem=args.assembly_mem, assembly_time=args.assembly_time,
        assembly_qos=args.assembly_qos, probe=args.probe, submit=args.submit)


def _summary(path: Path) -> str:
    plan = json.loads(path.read_text())
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    lines = [f"plan {path}", f"sha256 {digest}",
             f"{plan['policy']} plan of {len(plan['tasks'])} tasks of {plan['dataset_id']}, throttle "
             f"{plan['slurm']['throttle']}, licence cap {plan['licence_cap']}; code {plan['code']['commit']} in "
             f"{plan['code']['repo']}, bridge {plan['bridge']['commit']}"]
    for task in plan["tasks"]:
        lattice = task["lattice"]
        requested = "none" if task["l_resolution_mm"] is None else task["l_resolution_mm"]
        lines.append(f"  task {task['index']}: {task['frequency_ghz']} GHz, class {task['class']}, step along l "
                     f"requested {requested}, effective {task['effective_l_resolution_mm']}, fitted "
                     f"{task['l_step_mm']:.7g} mm, {lattice['counts'][1]} x {lattice['counts'][2]} exit nodes, AEDT "
                     f"temp {task['aedt_temp']['mode']}")
    for job in plan["include_jobs"]:
        lines.append(f"  include {job['dataset_id']} at {job['frequency_ghz']:g} GHz ({job['acceptance_kind']}): "
                     f"{job['job']}")
    return "\n".join(lines)


def main(argv: list[str] | None = None, *, env: EnvironmentState | None = None, guard: Any = None,
         repo: Path = REPO_ROOT, now: datetime | None = None) -> int:
    """The CLI: 0 with the plan written, 2 refused with nothing left, 1 when a failed plan's cleanup left something
    or a ``--submit`` submission failed (the plan stays). ``--log-submission RUN KIND ID`` first records a
    hand-submitted job (``array_submit.log_submission_main``). ``env``, ``guard``, ``repo`` and ``now`` replace the
    checkout's state, the shared guard, this checkout and the clock (tests)."""
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["--log-submission"]:
        from bbsim import array_submit

        return array_submit.log_submission_main(argv[1:])
    try:
        args = build_parser().parse_args(argv)
    except SystemExit as exc:                    # argparse has printed its message
        return 2 if exc.code else 0
    options = options_from_args(args)
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter("%(message)s"))
    level = log.level
    log.addHandler(handler)
    log.setLevel(logging.INFO)
    try:
        try:
            state = env if env is not None else environment_state(Path(repo))
        except EnvironmentStateError as exc:
            raise PlanError(f"the code checkout's state cannot be read: {exc}") from None
        except PackageNotFoundError as exc:
            raise PlanError(f"the installed PyAEDT version, which the plan records and every task checks, cannot be "
                            f"read: {exc}") from None
        shared = guard if guard is not None else licence_guard.Guard(options.guard_dir, licence_guard.SlurmScheduler())
        path = build_plan(options, env=state, guard=shared, now=now or datetime.now(timezone.utc), repo=Path(repo))
    except (PlanError, licence_guard.GuardError) as exc:
        print(f"plan_hfss_array: refused: {_one_line(exc)}", file=sys.stderr)
        return 2
    except PlanLeftovers as exc:
        print(f"plan_hfss_array: {_one_line(exc)}", file=sys.stderr)
        return 1
    finally:
        log.removeHandler(handler)
        log.setLevel(level)
    print(_summary(path))
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    print(submission_block(path, digest, json.loads(path.read_text())))
    if not options.submit:
        return 0
    from bbsim import array_submit

    return array_submit.submit_plan(path, digest, guard=shared, env=env)
