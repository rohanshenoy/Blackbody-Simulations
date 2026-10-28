"""The job-array plan (HF-057, OR-037, OR-048): one task's frozen configuration.

A plan (``plan.json``, schema ``PLAN_SCHEMA``) freezes a byte copy of the config, a copy of the prepared project and
one record per array task. ``task_config`` builds a task's run configuration from the config copy, or from the caller's
raw TOML mapping already captured from that copy, with ``config.build_config`` and exactly five overrides through the
validated override path:

- ``solver.frequency_ghz``: the task's ``frequency_ghz``;
- ``project.path``: the plan's ``project.copy``;
- ``output.root``: the plan's ``jobs_root`` (``<RUN>/jobs``);
- ``solver.cores``: the plan's ``slurm.cores``;
- ``exit_field.resolution_mm``: the frozen config's, with index 1 (the step along l; the list runs along p, l, g) the
  task's ``effective_l_resolution_mm``. That field is never null: it is the grid entry's requested step when the
  entry gives one (``l_resolution_mm``), else the config's own, after the c/(4f) rule (``frequency_grid``).

The runner, the array_run self-check (``bbsim.run_acceptance``) and the assembler all build a task's config here, so
they cannot disagree on it. Every path a plan records is absolute and resolved. A caller that hashes the config copy
passes the mapping it parsed from those bytes as ``raw_config``, and the copy is then never reopened.
"""
from __future__ import annotations

import math
import os
import tomllib
from collections.abc import Mapping
from typing import Any

from bbsim.config import RunConfig, build_config

PLAN_SCHEMA = "bbsim-array-plan/2"


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
