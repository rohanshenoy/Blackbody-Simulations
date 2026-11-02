"""The job-array plan's own assembler: one accepted run per task, the report and the check (HF-057, OR-027, CDX-038,
OR-045).

Usage::

    assemble_hfss_dataset.py --plan PLAN --plan-sha256 HEX [--plan PLAN2 --plan-sha256 HEX2]...
                             [--tree PATH] [--check | --report] [--report-path PATH]
                             [--choose [PLAN_SHA256:]INDEX=JOB]... [--expect-incidence THETAS:PHIS]
                             [--allow-unqualified [PLAN_SHA256:]INDEX]...

``load_plan`` reads a plan once, checks the sha256 the caller recorded for it against those bytes and parses the same
bytes. ``select`` finds every task's candidates and gives each a verdict:

- a directory under ``<jobs_root>/<stem>/`` whose ``start.json`` is readable and names this plan (path, sha256) and
  this task is a candidate; anything else under the plan's ``jobs_root`` is listed as foreign (no readable
  ``start.json``, another plan, another task, a stray file, a symlink);
- ``incomplete``: something completion requires is not there: ``done.json`` (an ``array_run.json`` without it is a task
  stopped between its report and its completion record), the manifest, the sidecar, the report, an Ephi directory or
  file the policy requires, or the ``run.log`` line its verdict requires (``done; manifest at <job>/manifest.json`` for
  a PASS, ``completed; acceptance FAIL (<rows>); manifest at <job>/manifest.json`` for a FAIL), whose absence beside a
  ``done.json`` is a task stopped after that record and before its last log line;
- ``invalid``: what is there disagrees: ``done.json`` malformed or naming another plan, task or report, or with digests
  other than the manifest's, the sidecar's and the report's; the manifest's ``plan`` block not this plan and task; Ephi
  directories, files or sidecar ``files`` other than the policy's (production: both polarizations; probe: Ephi=0
  alone), or a listed file whose sha256, size or rows differ; the report's verdict not ``done.json``'s; the start
  record's or the manifest's commit not the plan's; the sidecar's incidence not the expected one; an acceptance other
  than ``array_run``; or the gate refusing the report: ``acceptance_gate.report_problems`` with the caller's
  ``PlanContext`` (the plan, its sha256, the task index and the policy) and ``acceptance_gate.convergence_problems``
  at the policy's polarizations;
- ``failed``: a deterministic FAIL, a well-formed bound report and ``done.json`` recording it, listed with its failed
  rows; its binding and context are judged as the gate judges a passed report, with its own verdict set aside, so a
  FAIL is told from a report that no longer binds these bytes; a FAIL is never accepted, except:
- ``valid``: complete, and the gate accepts the report. For a task named by ``--allow-unqualified`` a FAIL whose only
  failed row is the waivable ``convergence`` row, beside passed ``convergence.records`` and otherwise usable
  convergence records, is valid with ``waived = "convergence"`` (``report_problems(..., waive={"convergence"})`` and
  ``convergence_problems(..., allow_unqualified=True)``); the report, ``done.json`` and the Slurm state are never
  rewritten, and the override is recorded.

A task with exactly one valid candidate is selected (a FAIL beside it is no ambiguity); several valid ones are
ambiguous until ``--choose`` names one; with none it is ``failed`` when a deterministic FAIL exists, else ``missing``.

``report`` prints, and writes as ``<RUN>/assembly_<utc>.json`` (the first plan's run directory, or ``--report-path``),
every candidate with its verdict and reasons and every foreign entry; for a task with no valid run and no FAIL its
``sacct -j <A>_<index>`` lines (``submissions.log``'s arrays holding it) and the guarded rerun block from the plan's
frozen code (``python -m bbsim.array_submit rerun`` with the class's resources, then a new assembly after the reruns,
any pending assembly of the plan's cancelled first), and the same block with twice the class's time (TIMEOUT) and
twice its memory (OUT_OF_MEMORY), never more cores; ``--choose`` lines (as ``--check`` lines) for an ambiguous task;
an ``--allow-unqualified`` line for a FAIL the override could accept, which it never applies unasked; for every array
of the plan ``sacct``, the ``du -sh`` of each job and its ``.aedtresults``, and the hosts' ``TempDirectory`` lines.
Every line holds real values (paths, digests, indices, job ids), so it is pasted as printed. Its last line is ``REPORT
PASS``/``REPORT FAIL: ...`` or ``CHECK PASS``/``CHECK FAIL: ...``.

The CLI: ``--report`` reports any plan, a probe's too; ``--check`` is the report under the assembly's refusals: a probe
plan, a plan's worktree ``<RUN>/code`` off the plan's commit or with uncommitted changes to tracked files
(``code_state``), an existing tree or ``<TREE>.partial``. The assembly proper (neither option) takes the same refusals
and does not write a tree in this version. Exit codes: 2 refused before anything was judged or written, 1 a task without
exactly one valid run (or a report that could not be written), 0 every task with one.

Imports: the standard library, ``bbsim.naming`` and ``bbsim.acceptance_gate`` at module level (the binding, the
required acceptance and the convergence gate, one implementation shared with the bridge's ``assemble_tree.sbatch``;
HF-057); ``acceptance_gate.report_problems`` imports the checkers (numpy, pandas), so the selection runs in the bb
venv. ``bbsim.array_plan`` (standard library at import) and ``bbsim.manifest`` are imported
where they are used.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shlex
import subprocess
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from bbsim import acceptance_gate as gate
from bbsim.naming import dataset_dir_name, dataset_record_name, frequency_label

VALID, FAILED, INCOMPLETE, INVALID, FOREIGN = "valid", "failed", "incomplete", "invalid", "foreign"
SELECTED, MISSING, AMBIGUOUS, UNRESOLVED = "selected", "missing", "ambiguous", "unresolved"
MODES = ("report", "check", "assembly")
REPORT_SCHEMA = "bbsim-assembly-report/1"
START_NAME, DONE_NAME, LOG_NAME = "start.json", "done.json", "run.log"
KIND = "array_run"
REPORT_NAME = gate.REPORT_NAMES[KIND]
WAIVE = frozenset({"convergence"})
CSV_NAMES = ("waveguide.csv", "far_field.csv")
ANGLE_TOL_DEG = 1e-9
SACCT_TASK = "-n -o JobID%18,State,ExitCode,Elapsed,MaxRSS"
SACCT_ARRAY = "-o JobID%18,State,ExitCode,Elapsed,MaxRSS,NodeList"
TEMP_GREP = "grep -i tempdirectory ~/Ansoft/ElectronicsDesktop2025.2/config/*.cfg"
# The runner's last log line of a completed task (run_frequency's plan mode, LOG_FORMAT "... bbsim.run: <message>")
PASS_LINE = re.compile(r" bbsim\.run: done; manifest at (?P<manifest>.+)$")
FAIL_LINE = re.compile(r" bbsim\.run: completed; acceptance FAIL \((?P<rows>.*)\); manifest at (?P<manifest>.+)$")
_HEX64 = re.compile(r"[0-9a-f]{64}")
_COMMIT = re.compile(r"[0-9a-f]{40}|[0-9a-f]{64}")
_MEM = re.compile(r"([1-9][0-9]*)([KMGT]?)")
_KEY = re.compile(r"(?:(?P<sha>[0-9a-f]{64}):)?(?P<index>[0-9]+)")

Key = tuple[str, int]          # (plan sha256, task index)


class AssembleRefused(ValueError):
    """Refused before anything was judged or written (exit 2); the message names what and why."""


@dataclass(frozen=True)
class LoadedPlan:
    """A plan as ``load_plan`` read it: its resolved path, the sha256 of its bytes and what those bytes parse to."""
    path: Path
    sha256: str
    data: dict

    @property
    def policy(self) -> str:
        return self.data["policy"]

    @property
    def run_dir(self) -> Path:
        return Path(self.data["run_dir"])

    @property
    def code_repo(self) -> Path:
        return Path(self.data["code"]["repo"])

    def tasks(self) -> list[dict]:
        return sorted(self.data["tasks"], key=lambda t: t["index"])

    def task(self, index: int) -> dict:
        return next(t for t in self.data["tasks"] if t["index"] == index)

    def stem(self, task: Mapping) -> str:
        return f"{self.data['dataset_id']}_{task['label']}"


@dataclass
class Candidate:
    """One job directory of a task and its verdict (``valid``, ``failed``, ``incomplete``, ``invalid``), with the
    reasons and, for a run accepted by the convergence override, ``waived = "convergence"``. ``waivable`` marks a FAIL
    the override could accept (the report prints its line; it is never applied unasked)."""
    job: Path
    index: int
    verdict: str
    reasons: list[str]
    waived: str | None = None
    waivable: bool = False


@dataclass
class Selection:
    """``select``'s result. Keys are (plan sha256, task index). ``states`` is each task's ``selected``, ``missing``,
    ``ambiguous``, ``failed`` or ``unresolved`` (a ``--choose`` naming a run that is not valid), ``unresolved`` the
    reason a task has no chosen run; ``foreign_reasons`` why each foreign entry is one; ``notes`` warnings the CLI
    adds (a plan's worktree off its commit under ``--report``)."""
    plans: list[LoadedPlan]
    candidates: dict[Key, list[Candidate]]
    chosen: dict[Key, Candidate]
    foreign: list[Path]
    foreign_reasons: dict[Path, str] = field(default_factory=dict)
    states: dict[Key, str] = field(default_factory=dict)
    unresolved: dict[Key, str] = field(default_factory=dict)
    choose: dict[Key, Path] = field(default_factory=dict)
    allow_unqualified: frozenset[Key] = frozenset()
    expect_incidence: dict | None = None
    notes: list[str] = field(default_factory=list)

    def plan(self, sha256: str) -> LoadedPlan:
        return next(p for p in self.plans if p.sha256 == sha256)


def _one_line(text: object) -> str:
    return " ".join(str(text).split())


def _show(value: object) -> str:
    try:
        return json.dumps(value)
    except (TypeError, ValueError):
        return repr(value)


def _now() -> datetime:
    return datetime.now(timezone.utc)


# --- the plan --------------------------------------------------------------------------------------------------------

def _absolute(value: object, what: str) -> Path:
    if not isinstance(value, str) or not os.path.isabs(value):
        raise AssembleRefused(f"the plan's {what} must be an absolute path, as every path a plan records is, got "
                              f"{value!r}")
    return Path(value)


def load_plan(path: Path | str, sha256: str) -> LoadedPlan:
    """The plan at ``path``, whose bytes must have the sha256 ``sha256`` the caller recorded, parsed from those same
    bytes: a ``bbsim-array-plan/2`` plan of the ``array_run`` kind, its run directory's ``plan.json``, with a policy,
    a dataset, tasks of unique indices whose labels are their frequencies', classes, a worktree and a jobs root.
    Raises ``AssembleRefused``."""
    from bbsim.array_plan import PLAN_SCHEMA

    if not isinstance(sha256, str) or not _HEX64.fullmatch(sha256):
        raise AssembleRefused(f"a plan's sha256 must be 64 lowercase hex digits, got {sha256!r}")
    try:
        resolved = Path(os.fspath(path)).resolve(strict=True)
        data = resolved.read_bytes()
    except (OSError, RuntimeError, TypeError) as exc:
        raise AssembleRefused(f"the plan {path} cannot be read ({exc})") from None
    digest = hashlib.sha256(data).hexdigest()
    if digest != sha256:
        raise AssembleRefused(f"the plan {resolved} has sha256 {digest}, not the {sha256} given: it is not the plan "
                              "that was frozen")
    try:
        plan = json.loads(data.decode("utf-8"))
    except (UnicodeDecodeError, ValueError) as exc:
        raise AssembleRefused(f"the plan {resolved} is not JSON ({exc})") from None
    if not isinstance(plan, dict) or plan.get("schema") != PLAN_SCHEMA:
        raise AssembleRefused(f"the plan {resolved} is not a {PLAN_SCHEMA} plan")
    run_dir = _absolute(plan.get("run_dir"), "run_dir")
    if run_dir.resolve() / "plan.json" != resolved:
        raise AssembleRefused(f"the plan {resolved} is not its run directory's plan.json (run_dir {run_dir})")
    policy = plan.get("policy")
    if policy not in gate.ARRAY_POLICIES or plan.get("probe") is not (policy == "probe"):
        raise AssembleRefused(f"the plan's policy {policy!r} and probe flag {plan.get('probe')!r} are not production "
                              "or probe")
    if plan.get("acceptance_kind") != KIND:
        raise AssembleRefused(f"the plan's acceptance_kind is {plan.get('acceptance_kind')!r}, not {KIND!r}")
    if not isinstance(plan.get("dataset_id"), str) or not plan["dataset_id"]:
        raise AssembleRefused(f"the plan records no dataset_id ({plan.get('dataset_id')!r})")
    code = plan.get("code")
    if not isinstance(code, Mapping):
        raise AssembleRefused("the plan records no code {repo, commit}")
    _absolute(code.get("repo"), "code.repo")
    if not isinstance(code.get("commit"), str) or not _COMMIT.fullmatch(code["commit"]):
        raise AssembleRefused(f"the plan's code.commit {code.get('commit')!r} is not a commit")
    _absolute(plan.get("jobs_root"), "jobs_root")
    tasks = plan.get("tasks")
    if not isinstance(tasks, list) or not tasks or not all(isinstance(t, dict) for t in tasks):
        raise AssembleRefused("the plan's tasks are not a nonempty list of records")
    seen: set[int] = set()
    for task in tasks:
        index, frequency, label = task.get("index"), task.get("frequency_ghz"), task.get("label")
        if type(index) is not int or index < 0 or index in seen:
            raise AssembleRefused(f"the plan's task index {index!r} is not a unique nonnegative integer")
        seen.add(index)
        if isinstance(frequency, bool) or not isinstance(frequency, (int, float)) or not (
                math.isfinite(frequency) and frequency > 0) or label != frequency_label(frequency):
            raise AssembleRefused(f"the plan's task {index} records frequency {frequency!r} and label {label!r}, not "
                                  "a positive frequency and its label")
    classes = plan.get("classes")
    if not isinstance(classes, list) or not all(isinstance(c, Mapping) and isinstance(c.get("name"), str)
                                                 and isinstance(c.get("indices"), list) for c in classes):
        raise AssembleRefused("the plan's classes are not a list of records {name, indices, ...}")
    held = sorted(i for c in classes for i in c["indices"])
    if held != sorted(seen) or len({c["name"] for c in classes}) != len(classes):
        raise AssembleRefused(f"the plan's classes {[c['name'] for c in classes]} do not partition its task indices "
                              f"{sorted(seen)} (they hold {held})")
    return LoadedPlan(resolved, digest, plan)


def code_state(repo: Path) -> tuple[str, bool]:
    """``HEAD`` of the checkout ``repo`` and whether a tracked file differs from it (the git commands of
    ``manifest.environment_state``; untracked files never count), with no ``GIT_*`` variable of the caller's.
    Raises ``AssembleRefused`` when git cannot tell (a worktree whose source checkout moved is no repository)."""
    from bbsim.manifest import GIT_TIMEOUT_S

    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}

    def git(*args: str) -> str:
        command = ["git", "--no-optional-locks", "-C", str(repo), *args]
        try:
            run = subprocess.run(command, capture_output=True, text=True, timeout=GIT_TIMEOUT_S, env=env)
        except subprocess.TimeoutExpired:
            raise AssembleRefused(f"{' '.join(command)} timed out after {GIT_TIMEOUT_S} s") from None
        except OSError as exc:
            raise AssembleRefused(f"git could not run on {repo}: {exc}") from None
        if run.returncode != 0:
            raise AssembleRefused(f"git cannot read the code {repo}: {' '.join(command)} exited {run.returncode}: "
                                  f"{_one_line(run.stderr)}")
        return run.stdout

    commit = git("rev-parse", "--verify", "HEAD").strip()
    if not _COMMIT.fullmatch(commit):
        raise AssembleRefused(f"git rev-parse HEAD in {repo} printed {commit!r}, not a commit")
    return commit, bool(git("status", "--porcelain", "--untracked-files=no").strip())


def worktree_problem(plan: LoadedPlan) -> str | None:
    """Why the plan's worktree ``<RUN>/code`` is not the code the plan froze (another ``HEAD``, a changed tracked file,
    or git cannot read it), or None. Every task refuses such a worktree, and so does the assembly."""
    repo, commit = plan.code_repo, plan.data["code"]["commit"]
    try:
        head, dirty = code_state(repo)
    except AssembleRefused as exc:
        return f"the plan {plan.path}'s worktree cannot be read: {exc}"
    if head != commit:
        return (f"the plan {plan.path}'s worktree {repo} is at {head}, not the plan's commit {commit}: its HEAD moved "
                "after the plan was made, and every task and the assembly run the plan's own code")
    if dirty:
        return (f"the plan {plan.path}'s worktree {repo} has uncommitted changes to tracked files; restore it to the "
                f"plan's commit {commit}")
    return None


# --- one candidate's verdict ------------------------------------------------------------------------------------------

class _Stop(Exception):
    def __init__(self, verdict: str, reasons: Sequence[str]):
        super().__init__(verdict)
        self.verdict, self.reasons = verdict, list(reasons)


def _object(inputs: gate.InputCapture, path: Path, what: str) -> dict:
    """A JSON object read through ``inputs``: missing or unreadable is incomplete, changed while read invalid."""
    try:
        data = inputs.parse(path, json.load)
    except gate.InputChanged as exc:
        raise _Stop(INVALID, [str(exc)]) from None
    except FileNotFoundError:
        raise _Stop(INCOMPLETE, [f"no {what} {path.name}: the run did not complete"]) from None
    except (OSError, ValueError) as exc:
        raise _Stop(INCOMPLETE, [f"its {what} {path.name} is unreadable ({_one_line(exc)})"]) from None
    if not isinstance(data, dict):
        raise _Stop(INVALID, [f"its {what} {path.name} is not a JSON object"])
    return data


def _same_path(value: object, path: Path) -> bool:
    return isinstance(value, str) and os.path.isabs(value) and Path(value).resolve() == path.resolve()


def _plan_block_problems(name: str, block: object, plan: LoadedPlan, index: int) -> list[str]:
    if not isinstance(block, Mapping):
        return [f"its {name} records no plan block {{path, sha256, index}} (got {_show(block)})"]
    problems = []
    if not _same_path(block.get("path"), plan.path):
        problems.append(f"its {name} names the plan {block.get('path')!r}, not {plan.path}")
    if block.get("sha256") != plan.sha256:
        problems.append(f"its {name} names the plan sha256 {block.get('sha256')!r}, not {plan.sha256}")
    if type(block.get("index")) is not int or block["index"] != index:
        problems.append(f"its {name} names task {block.get('index')!r}, not task {index}")
    return problems


def _rows(path: Path, digest: str | None) -> tuple[int | None, str | None]:
    """The data rows of a CSV, the non-blank lines after its header, counted by streaming the file; with the problem
    when those bytes are not the ones ``inputs`` hashed (``digest``)."""
    hashed, rows = hashlib.sha256(), 0
    try:
        with open(path, "rb") as fh:
            header = True
            for line in fh:
                hashed.update(line)
                if header:
                    header = False
                elif line.strip():
                    rows += 1
    except OSError as exc:
        return None, f"{path}: its rows cannot be counted ({_one_line(exc)})"
    if digest is not None and hashed.hexdigest() != digest:
        return None, f"{path} changed while it was being checked: its rows were counted on other bytes"
    return rows, None


def _completion(plan: LoadedPlan, task: Mapping, job: Path, start: Mapping,
                inputs: gate.InputCapture) -> tuple[dict, dict, dict, dict]:
    """What completion requires (``select``'s rules); returns done.json, the manifest, the sidecar and the report, or
    raises ``_Stop`` with every reason of the stage that found one."""
    index, pols = task["index"], gate.ARRAY_POLICIES[plan.policy]
    dataset_id, frequency = plan.data["dataset_id"], task["frequency_ghz"]
    done_path, report_path = job / DONE_NAME, job / REPORT_NAME
    if not done_path.exists():
        if report_path.exists():
            raise _Stop(INCOMPLETE, [f"{REPORT_NAME} exists without {DONE_NAME}: the task was stopped between writing "
                                     "its report and its completion record, so it has no verdict"])
        raise _Stop(INCOMPLETE, [f"no {DONE_NAME} and no {REPORT_NAME}: the run did not complete (it failed, was "
                                 "cancelled or is still running)"])
    done = _object(inputs, done_path, "completion record")
    bad: list[str] = []
    if done.get("solve_completed") is not True:
        bad.append(f"its {DONE_NAME} records solve_completed {_show(done.get('solve_completed'))}, not true")
    if not isinstance(done.get("acceptance_passed"), bool):
        bad.append(f"its {DONE_NAME} records acceptance_passed {_show(done.get('acceptance_passed'))}, not a verdict")
    for name in ("manifest_sha256", "sidecar_sha256", "report_sha256"):
        if not (isinstance(done.get(name), str) and _HEX64.fullmatch(done[name])):
            bad.append(f"its {DONE_NAME} records {name} {_show(done.get(name))}, not a sha256")
    if done.get("plan_sha256") != plan.sha256 or type(done.get("index")) is not int or done["index"] != index:
        bad.append(f"its {DONE_NAME} names plan {done.get('plan_sha256')!r} task {done.get('index')!r}, not "
                   f"{plan.sha256} task {index}")
    if not _same_path(done.get("report"), report_path):
        bad.append(f"its {DONE_NAME} names the report {done.get('report')!r}, not {report_path}")
    if not _same_path(start.get("job_root"), job):
        bad.append(f"its {START_NAME} names the job {start.get('job_root')!r}, not {job}")
    if bad:
        raise _Stop(INVALID, bad)

    sidecar_path = job / dataset_record_name(dataset_id, frequency)
    files = {"manifest": job / "manifest.json", "sidecar": sidecar_path, "report": report_path}
    missing: list[str] = []
    for what, path in files.items():
        if not path.is_file():
            missing.append(f"no {what} {path.name}, which its {DONE_NAME} records")
    expected_dirs = {dataset_dir_name(dataset_id, frequency, e) for e in pols}
    found_dirs = {p.name for p in job.glob("*_Ephi=*")}
    for name in sorted(expected_dirs - found_dirs):
        missing.append(f"no {name}/, which the {plan.policy} policy exports")
    for name in sorted(expected_dirs & found_dirs):
        for file in gate.EPHI_FILES:
            if not (job / name / file).is_file():
                missing.append(f"no {name}/{file}")
    log_path = job / "logs" / LOG_NAME
    if not log_path.is_file():
        missing.append(f"no logs/{LOG_NAME}")
    if missing:
        raise _Stop(INCOMPLETE, missing)
    manifest = _object(inputs, files["manifest"], "manifest")
    sidecar = _object(inputs, sidecar_path, "sidecar")
    report = _object(inputs, report_path, "report")

    for what, path in files.items():
        if inputs.recorded(path) != done[f"{what}_sha256"]:
            bad.append(f"its {what} {path.name} has sha256 {inputs.recorded(path)}, not the {done[f'{what}_sha256']} "
                       f"its {DONE_NAME} recorded: it changed after the run completed")
    bad += _plan_block_problems("manifest", manifest.get("plan"), plan, index)
    for name in sorted(found_dirs - expected_dirs):
        bad.append(f"{name}/ is not among the {plan.policy} policy's Ephi directories {sorted(expected_dirs)}")
    for name in sorted(expected_dirs):
        directory = job / name
        if directory.is_symlink() or not directory.is_dir():
            bad.append(f"{name} is not a directory of the job")
            continue
        extra = sorted(p.name for p in directory.iterdir() if p.name not in gate.EPHI_FILES)
        if extra:
            bad.append(f"{name}/ holds {extra} beside its {list(gate.EPHI_FILES)}, and nothing else belongs there")
        bad += [f"{name}/{file} is a symlink" for file in gate.EPHI_FILES if (directory / file).is_symlink()]
    others = sorted(p.name for p in job.glob("*.dataset.json") if p != sidecar_path)
    if others:
        bad.append(f"the job holds the sidecars {others} beside its own {sidecar_path.name}")
    listed = sidecar.get("files")
    expected_files = {f"{d}/{n}" for d in expected_dirs for n in CSV_NAMES}
    if not isinstance(listed, Mapping) or set(listed) != expected_files:
        bad.append(f"its sidecar lists the files {sorted(listed) if isinstance(listed, Mapping) else listed!r}, not "
                   f"the {plan.policy} policy's {sorted(expected_files)}")
    else:
        bad += gate.sidecar_file_problems(job, sidecar, inputs)
        for rel in sorted(expected_files):
            entry, path = listed[rel], job / rel
            rows, problem = _rows(path, inputs.recorded(path))
            if problem is not None:
                bad.append(problem)
            elif not isinstance(entry, Mapping) or type(entry.get("rows")) is not int or entry["rows"] != rows:
                bad.append(f"{path}: {rows} rows, the sidecar records "
                           f"{entry.get('rows') if isinstance(entry, Mapping) else entry!r}")
    if bad:
        raise _Stop(INVALID, bad)

    try:
        log_text = inputs.read_text(log_path)
    except (OSError, ValueError) as exc:
        raise _Stop(INCOMPLETE, [f"its logs/{LOG_NAME} is unreadable ({_one_line(exc)})"]) from None
    manifest_path = files["manifest"]
    passes = [m for m in map(PASS_LINE.search, log_text.splitlines()) if m and _same_path(m["manifest"],
                                                                                          manifest_path)]
    fails = [m for m in map(FAIL_LINE.search, log_text.splitlines()) if m and _same_path(m["manifest"],
                                                                                        manifest_path)]
    passed = done["acceptance_passed"]
    wanted, other = (passes, fails) if passed else (fails, passes)
    if other:
        raise _Stop(INVALID, [f"its logs/{LOG_NAME} records the {'FAIL' if passed else 'PASS'} completion line, and "
                              f"its {DONE_NAME} the {'PASS' if passed else 'FAIL'} verdict"])
    if not wanted:
        line = (f"done; manifest at {manifest_path}" if passed else
                f"completed; acceptance FAIL (...); manifest at {manifest_path}")
        raise _Stop(INCOMPLETE, [f"its logs/{LOG_NAME} holds no '{line}' line, which its {DONE_NAME} verdict requires: "
                                 f"the task was stopped after {DONE_NAME} and before its last log line"])
    if report.get("passed") is not passed:
        raise _Stop(INVALID, [f"its {REPORT_NAME} records passed {_show(report.get('passed'))}, and its {DONE_NAME} "
                              f"acceptance_passed {_show(passed)}"])
    return done, manifest, sidecar, report


def _run_problems(plan: LoadedPlan, job: Path, start: Mapping, manifest: Mapping, sidecar: Mapping,
                  expect: Mapping) -> list[str]:
    """The commit the start record and the manifest name is the plan's, the manifest's code was clean, and the
    sidecar's incidence is the expected one."""
    problems = []
    commit = plan.data["code"]["commit"]
    code = start.get("code") if isinstance(start.get("code"), Mapping) else {}
    if code.get("commit") != commit:
        problems.append(f"its {START_NAME} records the commit {code.get('commit')!r}, not the plan's {commit}")
    if manifest.get("git_commit") != commit:
        problems.append(f"its manifest records git_commit {manifest.get('git_commit')!r}, not the plan's {commit}")
    state = manifest.get("environment_state")
    if isinstance(state, Mapping) and state.get("dirty") is not False:
        problems.append(f"its manifest records its code dirty ({_show(state.get('dirty'))})")
    excitation = sidecar.get("excitation") if isinstance(sidecar.get("excitation"), Mapping) else {}
    for axis in ("theta", "phi"):
        found, wanted = excitation.get(f"incident_{axis}_deg"), expect.get(f"{axis}_deg")
        if not (isinstance(found, list) and isinstance(wanted, list) and len(found) == len(wanted) and all(
                isinstance(a, (int, float)) and not isinstance(a, bool) and abs(float(a) - float(b)) <= ANGLE_TOL_DEG
                for a, b in zip(found, wanted))):
            problems.append(f"its sidecar's incident {axis} {found!r} deg is not the expected {wanted!r}")
    return problems


def _failed_rows(report: Mapping) -> list[str]:
    checks = report.get("checks") if isinstance(report.get("checks"), list) else []
    return [f"{c.get('name')}: {_one_line(c.get('detail'))}" for c in checks
            if isinstance(c, Mapping) and c.get("severity") == "fail" and c.get("passed") is not True]


def _verdict_set_aside(report: Mapping) -> dict:
    """The report with its verdict and its rows' verdicts set to passed, nothing else changed: what the gate then still
    finds wrong (binding, context, tolerances, required rows) is wrong whatever the verdict. Never accepted."""
    checks = report.get("checks")
    if isinstance(checks, list):
        checks = [dict(c, passed=True) if isinstance(c, Mapping) else c for c in checks]
    return dict(report, passed=True, checks=checks)


def judge(plan: LoadedPlan, task: Mapping, job: Path, start: Mapping, inputs: gate.InputCapture, *,
          waiver: bool = False, expect: Mapping | None = None) -> Candidate:
    """The verdict of one candidate ``job`` of ``task`` (the module docstring's rules), every file read through
    ``inputs``; ``waiver`` when the caller named this task with ``--allow-unqualified``, ``expect`` the incidence
    {theta_deg, phi_deg} its sidecar must declare (the plan's ``expect_incidence`` by default)."""
    index = task["index"]
    pols = gate.ARRAY_POLICIES[plan.policy]
    expect = expect if expect is not None else plan.data.get("expect_incidence") or {}
    try:
        _, manifest, sidecar, report = _completion(plan, task, job, start, inputs)
        problems = _run_problems(plan, job, start, manifest, sidecar, expect)
        if problems:
            raise _Stop(INVALID, problems)
        try:
            requirement = gate.required_acceptance(job, manifest, sidecar)
        except gate.GateError as exc:
            raise _Stop(INVALID, [str(exc)]) from None
        if requirement.kind != KIND:
            raise _Stop(INVALID, [f"its manifest records the acceptance {requirement.kind!r} ({requirement.source}), "
                                  f"not {KIND!r}"])
        context = gate.PlanContext(str(plan.path), plan.sha256, index, plan.policy)

        def problems_of(subject: Mapping, waive: frozenset = frozenset(), allow: bool = False) -> list[str]:
            return (gate.report_problems(subject, KIND, job, requirement, inputs, context=context, waive=waive)
                    + gate.convergence_problems(manifest, allow_unqualified=allow, expected_polarizations=pols))

        if report["passed"] is True:
            problems = problems_of(report)
            if problems:
                raise _Stop(INVALID, problems)
            return Candidate(job, index, VALID, [])
        failed = _failed_rows(report)
        waived = problems_of(report, WAIVE, True)
        if waiver and not waived:
            return Candidate(job, index, VALID, [f"accepted with the convergence override (--allow-unqualified, a "
                                                 f"diagnostic scope, OR-027): its {REPORT_NAME} failed {failed}; the "
                                                 f"report and {DONE_NAME} are unchanged"], waived="convergence")
        note = ([f"the convergence override does not apply: {'; '.join(waived)}"] if waiver else [])
        problems = gate.report_problems(_verdict_set_aside(report), KIND, job, requirement, inputs, context=context)
        if problems:
            raise _Stop(INVALID, problems + note)
        reasons = [f"its {REPORT_NAME} is a deterministic FAIL: {row}" for row in failed] or [
            f"its {REPORT_NAME} records passed false"]
        return Candidate(job, index, FAILED, reasons + note, waivable=not waived)
    except _Stop as stop:
        return Candidate(job, index, stop.verdict, stop.reasons)


# --- the selection ----------------------------------------------------------------------------------------------------

def _start(inputs: gate.InputCapture, job: Path) -> tuple[dict | None, str]:
    path = job / START_NAME
    try:
        start = inputs.parse(path, json.load)
    except (OSError, ValueError) as exc:
        return None, f"no readable {START_NAME} ({_one_line(exc)})"
    if not isinstance(start, dict):
        return None, f"its {START_NAME} is not a JSON object"
    return start, ""


def _scan(plan: LoadedPlan, inputs: gate.InputCapture) -> tuple[dict[int, list[tuple[Path, dict]]], dict[Path, str]]:
    """Every entry under the plan's ``jobs_root``: the candidates of each task index (job, start record), and the
    foreign entries with the reason each is one."""
    root = Path(plan.data["jobs_root"])
    stems = {plan.stem(t): t["index"] for t in plan.tasks()}
    found: dict[int, list[tuple[Path, dict]]] = {i: [] for i in stems.values()}
    foreign: dict[Path, str] = {}
    if not root.is_dir():
        return found, foreign
    root = root.resolve()
    for entry in sorted(root.iterdir()):
        if entry.name not in stems or entry.is_symlink() or not entry.is_dir():
            foreign[entry] = f"not a task directory of the plan {plan.path} (its tasks' stems are {sorted(stems)})"
            continue
        index = stems[entry.name]
        for job in sorted(entry.iterdir()):
            if job.is_symlink() or not job.is_dir():
                foreign[job] = "not a job directory (a stray file or a symlink)"
                continue
            start, why = _start(inputs, job)
            if start is None:
                foreign[job] = why
                continue
            block = start.get("plan") if isinstance(start.get("plan"), Mapping) else {}
            if block.get("sha256") != plan.sha256:
                foreign[job] = f"its {START_NAME} names another plan (sha256 {block.get('sha256')!r})"
            elif not _same_path(block.get("path"), plan.path) or type(block.get("index")) is not int \
                    or block["index"] != index:
                foreign[job] = (f"its {START_NAME} names {block.get('path')!r} task {block.get('index')!r}, not "
                                f"{plan.path} task {index}")
            else:
                found[index].append((job, start))
    return found, foreign


def select(plans: list[LoadedPlan], *, choose: Mapping[Key, Path] | None = None,
           allow_unqualified: set[Key] | frozenset[Key] = frozenset(), inputs: gate.InputCapture | None = None,
           expect_incidence: Mapping | None = None) -> Selection:
    """Every task's candidates and verdicts, and the one run selected per task where there is one (the module
    docstring's rules), every file read through ``inputs`` (a fresh ``InputCapture`` by default). ``choose`` maps a
    task to the job ``--choose`` named; ``allow_unqualified`` the tasks the convergence override may accept;
    ``expect_incidence`` replaces every plan's ``expect_incidence``. Raises ``AssembleRefused`` for a plan given twice,
    or a key or a ``--choose`` job that names no task or no job of that task."""
    inputs = inputs if inputs is not None else gate.InputCapture()
    choose = {key: Path(job).resolve() for key, job in (choose or {}).items()}
    allow = frozenset(allow_unqualified)
    digests = [p.sha256 for p in plans]
    if not plans or len(set(digests)) != len(digests):
        raise AssembleRefused(f"the plans must be given once each, and at least one; got {digests}")
    tasks = {(p.sha256, t["index"]): (p, t) for p in plans for t in p.tasks()}
    for what, keys in (("--choose", choose), ("--allow-unqualified", allow)):
        for key in keys:
            if key not in tasks:
                raise AssembleRefused(f"{what} names task {key[1]} of the plan {key[0]}, which no plan given holds")
    sel = Selection(plans=list(plans), candidates={}, chosen={}, foreign=[], choose=choose, allow_unqualified=allow,
                    expect_incidence=None if expect_incidence is None else dict(expect_incidence))
    for plan in plans:
        found, foreign = _scan(plan, inputs)
        sel.foreign += sorted(foreign)
        sel.foreign_reasons.update(foreign)
        for task in plan.tasks():
            key = (plan.sha256, task["index"])
            sel.candidates[key] = [judge(plan, task, job, start, inputs, waiver=key in allow, expect=expect_incidence)
                                   for job, start in found[task["index"]]]
    for key, candidates in sel.candidates.items():
        valid = [c for c in candidates if c.verdict == VALID]
        if key in choose:
            match = [c for c in candidates if c.job == choose[key]]
            if not match:
                raise AssembleRefused(f"--choose {key[1]}={choose[key]} names no job of task {key[1]} of the plan "
                                      f"{key[0]} (its jobs: {[str(c.job) for c in candidates]})")
            if match[0].verdict == VALID:
                sel.chosen[key], sel.states[key] = match[0], SELECTED
            else:
                sel.states[key] = UNRESOLVED
                sel.unresolved[key] = f"--choose names {choose[key]}, which is {match[0].verdict}"
        elif len(valid) == 1:
            sel.chosen[key], sel.states[key] = valid[0], SELECTED
        elif valid:
            sel.states[key] = AMBIGUOUS
            sel.unresolved[key] = f"{len(valid)} valid runs; --choose names the one to assemble"
        elif any(c.verdict == FAILED for c in candidates):
            sel.states[key] = FAILED
            sel.unresolved[key] = "no valid run, and a deterministic FAIL (no rerun line: a rerun repeats it)"
        else:
            sel.states[key] = MISSING
            sel.unresolved[key] = "no valid run" + ("" if candidates else " and no job")
    return sel


# --- the report -------------------------------------------------------------------------------------------------------

def _submissions(plan: LoadedPlan) -> tuple[list[dict], list[str]]:
    """The plan's lines of ``<RUN>/submissions.log`` (those naming its sha256), and the lines that are not records."""
    path = plan.run_dir / "submissions.log"
    try:
        text = path.read_text(encoding="utf-8")
    except FileNotFoundError:
        return [], [f"{path} does not exist: no job of the plan was submitted through the guard"]
    except (OSError, ValueError) as exc:
        return [], [f"{path} is unreadable ({_one_line(exc)})"]
    records, problems = [], []
    for number, line in enumerate(text.splitlines(), start=1):
        try:
            record = json.loads(line)
        except ValueError:
            record = None
        if not isinstance(record, dict) or not isinstance(record.get("kind"), str) or not isinstance(
                record.get("id"), str) or not record["id"].isdigit():
            problems.append(f"{path} line {number} is not a submission record")
        elif record.get("plan_sha256") == plan.sha256:
            records.append(record)
    return records, problems


def _arrays(records: Sequence[Mapping], candidates: Sequence[Candidate], index: int | None = None) -> list[str]:
    """The array job ids of the plan (those holding ``index`` when given): ``submissions.log``'s, then any a job's
    start record names that the log lacks."""
    ids = [r["id"] for r in records if r.get("kind") == "array" and (
        index is None or (isinstance(r.get("indices"), list) and index in r["indices"]))]
    for cand in candidates:
        if index is None or cand.index == index:
            match = re.fullmatch(r"job_a([0-9]+)_[0-9]+", cand.job.name)
            if match:
                ids.append(match.group(1))
    return list(dict.fromkeys(ids))


def _double_mem(mem: str) -> str:
    match = _MEM.fullmatch(mem)
    return f"{2 * int(match.group(1))}{match.group(2)}" if match else mem


def _double_time(text: str) -> str:
    from bbsim import array_plan

    try:
        return array_plan._time_text(2 * array_plan._seconds(text, "time"))
    except array_plan.PlanError:
        return text


def _rerun_block(plan: LoadedPlan, by_class: Mapping[str, list[int]], assemblies: Sequence[str], variant: str,
                 current_job: str | None) -> list[str]:
    q = shlex.quote
    plan_arg, sha = q(str(plan.path)), plan.sha256
    classes = {c["name"]: c for c in plan.data["classes"]}
    lines = [f"cd {q(str(plan.code_repo))} && (",
             'bbsim_stop() { echo "STOP: $1 was not submitted (exit $2); what was submitted before it keeps running, '
             'and nothing after it was submitted."; exit 1; }']
    for job in assemblies:
        if job != current_job:
            lines.append(f"squeue -h -j {job} -t PENDING -o %i 2>/dev/null | grep -qx {job} && scancel {job}")
    names = []
    for number, (name, indices) in enumerate(sorted(by_class.items()), start=1):
        cls = classes[name]
        mem, time, qos = str(cls.get("mem")), str(cls.get("time")), cls.get("qos")
        if variant == "TIMEOUT":
            time = _double_time(time)
        elif variant == "OUT_OF_MEMORY":
            mem = _double_mem(mem)
        resources = f"--mem {q(mem)} --time {q(time)}" + (f" --qos {q(str(qos))}" if qos is not None else "")
        names.append(f"RERUN{number}")
        indices_arg = ",".join(map(str, sorted(indices)))
        lines.append(f"RERUN{number}=$(python -m bbsim.array_submit rerun {plan_arg} {sha} --indices {indices_arg} "
                     f'{resources}) || bbsim_stop "the rerun of class {name}" $?')
    dependency = "afterany:" + ":".join(f"${n}" for n in names)
    lines.append(f"ASSEMBLY=$(python -m bbsim.array_submit submit {plan_arg} {sha} assembly "
                 f'--dependency "{dependency}") || bbsim_stop "the new assembly" $?')
    lines += ['echo "submitted: ' + " ".join(f"{n}=${n}" for n in names) + ' assembly=$ASSEMBLY"', ")"]
    return lines


def _check_line(sel: Selection, *extra: str) -> str:
    """A ``--check`` line from the first plan's code with the selection's plans, choices and overrides, plus
    ``extra``."""
    q = shlex.quote
    several = len(sel.plans) > 1
    argv = []
    for plan in sel.plans:
        argv += ["--plan", q(str(plan.path)), "--plan-sha256", plan.sha256]
    argv.append("--check")
    for (sha, index), job in sorted(sel.choose.items()):
        argv += ["--choose", f"{sha + ':' if several else ''}{index}={q(str(job))}"]
    for sha, index in sorted(sel.allow_unqualified):
        argv += ["--allow-unqualified", f"{sha}:{index}"]
    if sel.expect_incidence is not None:
        argv += ["--expect-incidence", ",".join(repr(float(v)) for v in sel.expect_incidence["theta_deg"]) + ":" +
                 ",".join(repr(float(v)) for v in sel.expect_incidence["phi_deg"])]
    return f"cd {q(str(sel.plans[0].code_repo))} && python assemble_hfss_dataset.py {' '.join(argv + list(extra))}"


def _task_name(sel: Selection, key: Key) -> str:
    task = sel.plan(key[0]).task(key[1])
    where = f" of plan {key[0][:12]}" if len(sel.plans) > 1 else ""
    return f"task {key[1]} ({task['label']}, class {_class_of(sel.plan(key[0]), key[1])}{where})"


def _class_of(plan: LoadedPlan, index: int) -> str:
    names = [c["name"] for c in plan.data["classes"] if index in c["indices"]]
    return names[0] if len(names) == 1 else "?"


def _lines(sel: Selection, mode: str) -> tuple[list[str], bool, str]:
    """The printed lines (without the verdict line), whether every task has one valid run, and the summary."""
    q = shlex.quote
    current_job = os.environ.get("SLURM_JOB_ID") or None
    keys = [(p.sha256, t["index"]) for p in sel.plans for t in p.tasks()]
    lines = [f"== bbsim.assemble --{mode}: {len(sel.plans)} plan(s), {len(keys)} task(s)"]
    for plan in sel.plans:
        lines.append(f"plan {plan.path} sha256 {plan.sha256}: {plan.policy}, {plan.data['dataset_id']}, "
                     f"{len(plan.tasks())} task(s), code {plan.code_repo} at {plan.data['code']['commit']}")
    lines += [f"WARN: {note}" for note in sel.notes]
    for key in keys:
        state = sel.states[key]
        head = f"{_task_name(sel, key)}: {state}"
        if key in sel.chosen:
            chosen = sel.chosen[key]
            head += f" {chosen.job}" + (" (accepted with the convergence override)" if chosen.waived else "")
        else:
            head += f": {sel.unresolved[key]}"
        lines.append(head)
        for cand in sel.candidates[key]:
            lines.append(f"    {cand.verdict:<10} {cand.job}")
            lines += [f"               {reason}" for reason in cand.reasons]
    for path in sel.foreign:
        lines.append(f"foreign    {path}: {sel.foreign_reasons[path]}")

    unresolved = [k for k in keys if k not in sel.chosen]
    for plan in sel.plans:
        records, problems = _submissions(plan)
        lines += [f"WARN: {p}" for p in problems]
        cands = [c for t in plan.tasks() for c in sel.candidates[(plan.sha256, t["index"])]]
        rerun = [k for k in unresolved if k[0] == plan.sha256 and sel.states[k] == MISSING]
        if rerun:
            lines.append(f"# The plan {plan.path}'s tasks with no valid run and no deterministic FAIL, and their Slurm "
                         "records (a State of TIMEOUT or OUT_OF_MEMORY takes the matching block below):")
            for key in rerun:
                arrays = _arrays(records, cands, key[1])
                lines.append(f"# {_task_name(sel, key)}")
                lines += [f"sacct -j {a}_{key[1]} {SACCT_TASK}" for a in arrays] or [
                    f"# task {key[1]} is in no array submissions.log records"]
            by_class: dict[str, list[int]] = {}
            for key in rerun:
                by_class.setdefault(_class_of(plan, key[1]), []).append(key[1])
            assemblies = [r["id"] for r in records if r.get("kind") == "assembly"]
            lines.append("# Paste exactly ONE of the three blocks below, in a shell that has sourced bb_env.sh: "
                         "each reruns these tasks through the licence guard (OR-045) from the frozen code, cancels a "
                         "pending assembly of the plan first, and submits a new one that waits for the reruns.")
            for variant, title in (("", "Block 1, with each class's resources:"),
                                   ("TIMEOUT", "Block 2, instead of block 1 if sacct shows TIMEOUT: twice each class's "
                                               "time (raise --time, never --cpus-per-task):"),
                                   ("OUT_OF_MEMORY", "Block 3, instead of block 1 if sacct shows OUT_OF_MEMORY: twice "
                                                     "each class's memory (raise --mem, never --cpus-per-task):")):
                lines.append(f"# {title}")
                lines += _rerun_block(plan, by_class, assemblies, variant, current_job)
        elif any(k[0] == plan.sha256 for k in unresolved):
            assemblies = [r["id"] for r in records if r.get("kind") == "assembly" and r["id"] != current_job]
            lines.append(f"# A new assembly of the plan {plan.path} through the guard, once the tasks above are "
                         "settled (a pending one is cancelled first):")
            lines.append(f"cd {q(str(plan.code_repo))} && (" + "".join(
                f"squeue -h -j {a} -t PENDING -o %i 2>/dev/null | grep -qx {a} && scancel {a}; " for a in assemblies)
                + f"python -m bbsim.array_submit submit {q(str(plan.path))} {plan.sha256} assembly)")
    several = len(sel.plans) > 1
    for key in unresolved:
        if sel.states[key] == AMBIGUOUS:
            lines.append(f"# {_task_name(sel, key)} has several valid runs: check with the one to assemble (in a shell "
                         "that has sourced bb_env.sh), which the assembly takes with the same --choose:")
            for cand in sel.candidates[key]:
                if cand.verdict == VALID:
                    lines.append(_check_line(sel, "--choose", f"{key[0] + ':' if several else ''}{key[1]}="
                                                              f"{q(str(cand.job))}"))
        waivable = [c for c in sel.candidates[key] if c.verdict == FAILED and c.waivable]
        if waivable:
            lines.append(f"# {_task_name(sel, key)} failed only its convergence row, on otherwise usable evidence: the "
                         "diagnostic override accepts it only when asked (OR-027), recording it and rewriting nothing "
                         "(in a shell that has sourced bb_env.sh):")
            lines.append(_check_line(sel, "--allow-unqualified", f"{key[0]}:{key[1]}"))
    lines.append("# Every array of the plan(s), each job's disk use, and the hosts' AEDT temp-directory settings:")
    for plan in sel.plans:
        records, _ = _submissions(plan)
        cands = [c for t in plan.tasks() for c in sel.candidates[(plan.sha256, t["index"])]]
        lines += [f"sacct -j {a} {SACCT_ARRAY}" for a in _arrays(records, cands)]
        for cand in cands:
            results = sorted((cand.job / "project").glob("*.aedtresults"))
            lines.append("du -sh " + " ".join(q(str(p)) for p in [cand.job, *results]))
    lines.append(TEMP_GREP)
    counts: dict[str, int] = {}
    for key in unresolved:
        counts[sel.states[key]] = counts.get(sel.states[key], 0) + 1
    summary = ", ".join(f"{n} {state}" for state, n in sorted(counts.items()))
    return lines, not unresolved, summary


def _record(sel: Selection, mode: str, lines: Sequence[str], passed: bool, verdict: str) -> dict:
    from bbsim.manifest import REPO_ROOT

    try:
        commit, dirty = code_state(REPO_ROOT)
        assembler: dict[str, Any] = {"repo": str(REPO_ROOT), "commit": commit, "dirty": dirty}
    except AssembleRefused as exc:
        assembler = {"repo": str(REPO_ROOT), "error": str(exc)}
    tasks = []
    for plan in sel.plans:
        for task in plan.tasks():
            key = (plan.sha256, task["index"])
            chosen = sel.chosen.get(key)
            tasks.append({"plan_sha256": plan.sha256, "index": task["index"], "frequency_ghz": task["frequency_ghz"],
                          "label": task["label"], "class": _class_of(plan, task["index"]), "stem": plan.stem(task),
                          "state": sel.states[key], "chosen": None if chosen is None else str(chosen.job),
                          "waived": None if chosen is None else chosen.waived, "reason": sel.unresolved.get(key),
                          "candidates": [{"job": str(c.job), "verdict": c.verdict, "reasons": c.reasons,
                                          "waived": c.waived, "waivable": c.waivable} for c in sel.candidates[key]]})
    return {"schema": REPORT_SCHEMA, "mode": mode, "created_utc": _now().isoformat(), "passed": passed,
            "verdict": verdict,
            "plans": [{"path": str(p.path), "sha256": p.sha256, "policy": p.policy, "dataset_id": p.data["dataset_id"],
                       "run_dir": str(p.run_dir), "code": {"repo": str(p.code_repo),
                                                           "commit": p.data["code"]["commit"]}} for p in sel.plans],
            "assembler": assembler,
            "choose": [{"plan_sha256": s, "index": i, "job": str(j)} for (s, i), j in sorted(sel.choose.items())],
            "allow_unqualified": [{"plan_sha256": s, "index": i} for s, i in sorted(sel.allow_unqualified)],
            "expect_incidence": sel.expect_incidence, "tasks": tasks,
            "foreign": [{"path": str(p), "reason": sel.foreign_reasons[p]} for p in sel.foreign],
            "notes": list(sel.notes), "lines": list(lines)}


def report_path(sel: Selection, now: datetime | None = None) -> Path:
    """``<RUN>/assembly_<utc>.json`` in the first plan's run directory."""
    stamp = (now or _now()).strftime("%Y%m%dT%H%M%S.%fZ")
    return sel.plans[0].run_dir / f"assembly_{stamp}.json"


def report(selection: Selection, *, path: Path | None, mode: str = "report") -> int:
    """Print the report (the module docstring's lines) and write it as JSON to ``path`` (``report_path`` when None),
    through a temporary file and ``os.replace``; the last line printed is the verdict, ``<MODE> PASS`` or ``<MODE>
    FAIL: ...`` (none for the assembly, which prints its own). Returns 0 when every task has exactly one valid run, 1
    otherwise or when the report cannot be written."""
    if mode not in MODES:
        raise ValueError(f"mode must be one of {', '.join(MODES)}, got {mode!r}")
    lines, passed, summary = _lines(selection, mode)
    word = mode.upper()
    n = sum(len(p.tasks()) for p in selection.plans)
    verdict = (f"{word} PASS: {n} task(s) of {len(selection.plans)} plan(s), one valid run each" if passed else
               f"{word} FAIL: {n - len(selection.chosen)} of {n} task(s) without exactly one valid run ({summary})")
    for line in lines:
        print(line)
    target = Path(path) if path is not None else report_path(selection)
    tmp = target.with_name(f".{target.name}.tmp")
    try:
        tmp.write_text(json.dumps(_record(selection, mode, lines, passed, verdict), indent=2, sort_keys=True) + "\n")
        os.replace(tmp, target)
    except OSError as exc:
        if tmp.exists():
            tmp.unlink()
        print(f"the report could not be written to {target} ({_one_line(exc)})")
        if mode != "assembly":
            print(f"{word} FAIL: the report could not be written")
        return 1
    print(f"report: {target}")
    if mode != "assembly":
        print(verdict)
    return 0 if passed else 1


# --- the CLI ----------------------------------------------------------------------------------------------------------

def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="assemble_hfss_dataset.py", allow_abbrev=False,
                                description="Select one accepted run per task of a job-array plan through the gate, "
                                            "and report (HF-057, OR-027, CDX-038).")
    p.add_argument("--plan", type=Path, action="append", default=[], help="a frozen plan.json (repeatable)")
    p.add_argument("--plan-sha256", action="append", default=[], help="the plan's recorded sha256, one per --plan")
    p.add_argument("--tree", type=Path, help="the tree (default: the plans' recorded tree)")
    mode = p.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true", help="the report under the assembly's refusals")
    mode.add_argument("--report", action="store_true", help="the report only; probe plans are accepted")
    p.add_argument("--report-path", type=Path, help="write the JSON report here (default <RUN>/assembly_<utc>.json)")
    p.add_argument("--choose", action="append", default=[], metavar="[PLAN_SHA256:]INDEX=JOB",
                   help="the valid run to assemble for an ambiguous task")
    p.add_argument("--include-job", action="append", default=[], metavar="JOB",
                   help="a finished job of another dataset to assemble with the plan's runs")
    p.add_argument("--expect-incidence", metavar="THETAS:PHIS",
                   help="the incidence every sidecar declares, two comma lists of degrees (default the plan's)")
    p.add_argument("--allow-unqualified", action="append", default=[], metavar="[PLAN_SHA256:]INDEX",
                   help="accept a task's run whose only failed row is a missed convergence target (OR-027; diagnostic)")
    p.add_argument("--verify", type=Path, metavar="TREE", help="re-check a tree with the standard library only")
    return p


def _key(text: str, plans: Sequence[LoadedPlan], option: str, job: Path | None = None) -> Key:
    """``[PLAN_SHA256:]INDEX`` as a task key. Without a digest the plan is the only one, or for ``--choose`` the one
    whose jobs root holds ``job``; ``--allow-unqualified`` over several plans needs the digest."""
    match = _KEY.fullmatch(text)
    if match is None:
        raise AssembleRefused(f"{option} {text!r} is not [PLAN_SHA256:]INDEX")
    index = int(match["index"])
    if match["sha"] is not None:
        matching = [p for p in plans if p.sha256 == match["sha"]]
        if not matching:
            raise AssembleRefused(f"{option} {text}: no plan given has the sha256 {match['sha']}")
        plan = matching[0]
    elif len(plans) == 1:
        plan = plans[0]
    elif job is not None:
        holders = [p for p in plans if Path(p.data["jobs_root"]).resolve() in job.parents]
        if len(holders) != 1:
            raise AssembleRefused(f"{option} {text}: {len(holders)} of the plans given hold {job}; name the plan, "
                                  "PLAN_SHA256:INDEX")
        plan = holders[0]
    else:
        raise AssembleRefused(f"{option} {text} is ambiguous over {len(plans)} plans: name the plan, "
                              "PLAN_SHA256:INDEX")
    if index not in {t["index"] for t in plan.tasks()}:
        raise AssembleRefused(f"{option} {text}: the plan {plan.path} has no task {index}")
    return plan.sha256, index


def _tree(args: argparse.Namespace, plans: Sequence[LoadedPlan]) -> Path:
    if args.tree is not None:
        return Path(os.path.abspath(args.tree))
    trees = sorted({p.data.get("tree") for p in plans}, key=str)
    if len(trees) != 1 or not isinstance(trees[0], str) or not os.path.isabs(trees[0]):
        raise AssembleRefused(f"the plans record the trees {trees}, not one absolute tree: give --tree")
    return Path(trees[0])


def main(argv: list[str] | None = None) -> int:
    """The CLI (the module docstring's). Refusals print one ``refused:`` line on stderr and exit 2."""
    try:
        args = _parser().parse_args(argv)
    except SystemExit as exc:                     # argparse has printed its message
        return 2 if exc.code else 0
    mode = "check" if args.check else "report" if args.report else "assembly"
    try:
        if args.verify is not None:
            raise AssembleRefused("--verify is not available in this version of bbsim.assemble")
        if not args.plan or len(args.plan) != len(args.plan_sha256):
            raise AssembleRefused(f"give each --plan with its --plan-sha256 ({len(args.plan)} plans, "
                                  f"{len(args.plan_sha256)} digests)")
        plans = [load_plan(path, digest) for path, digest in zip(args.plan, args.plan_sha256)]
        if len({p.sha256 for p in plans}) != len(plans):
            raise AssembleRefused("a plan is given twice")
        choose: dict[Key, Path] = {}
        for text in args.choose:
            index, sep, job = text.partition("=")
            if not sep or not job:
                raise AssembleRefused(f"--choose takes [PLAN_SHA256:]INDEX=JOB, got {text!r}")
            path = Path(os.path.abspath(job)).resolve()
            key = _key(index, plans, "--choose", path)
            if key in choose:
                raise AssembleRefused(f"--choose names task {key[1]} twice")
            choose[key] = path
        allow = {_key(text, plans, "--allow-unqualified") for text in args.allow_unqualified}
        expect = None
        if args.expect_incidence is not None:
            from bbsim import array_plan

            try:
                expect = array_plan._expected_incidence(args.expect_incidence)
            except array_plan.PlanError as exc:
                raise AssembleRefused(str(exc)) from None
        if args.include_job:
            raise AssembleRefused("--include-job is not available in this version of bbsim.assemble")
        notes = []
        if mode == "report":
            notes = [problem for problem in map(worktree_problem, plans) if problem is not None]
        else:
            probes = [str(p.path) for p in plans if p.policy == "probe"]
            if probes:
                raise AssembleRefused(f"the probe plan(s) {probes} are diagnostics, never assembled: --report "
                                      "reports them")
            for plan in plans:
                problem = worktree_problem(plan)
                if problem is not None:
                    raise AssembleRefused(problem)
            included = [str(p.path) for p in plans if p.data.get("include_jobs")]
            if included:
                raise AssembleRefused(f"the plan(s) {included} record include jobs, which this version of "
                                      "bbsim.assemble does not check")
            tree = _tree(args, plans)
            for path in (tree, tree.with_name(tree.name + ".partial")):
                if os.path.lexists(path):
                    raise AssembleRefused(f"{path} exists; the assembly writes a new tree only"
                                          + (" (remove the .partial only after confirming no assembly job is running)"
                                             if path != tree else ""))
            if mode == "assembly":
                raise AssembleRefused("writing the tree is not available in this version of bbsim.assemble; --check "
                                      "runs its selection under the same refusals")
        if args.report_path is not None:
            target = Path(os.path.abspath(args.report_path))
            if os.path.lexists(target) or not target.parent.is_dir():
                raise AssembleRefused(f"--report-path {target} exists or has no directory; a report is a new file")
        else:
            target = None
        selection = select(plans, choose=choose, allow_unqualified=allow, expect_incidence=expect)
        selection.notes += notes
    except AssembleRefused as exc:
        print(f"assemble: refused: {_one_line(exc)}", file=sys.stderr)
        return 2
    return report(selection, path=target, mode=mode)


if __name__ == "__main__":
    sys.exit(main())
