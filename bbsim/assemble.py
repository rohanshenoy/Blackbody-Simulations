"""The job-array plan's own assembler: one accepted run per task, the report, the check and the tree (HF-057, OR-027,
CDX-038, OR-045, S65, G4-012).

Usage::

    assemble_hfss_dataset.py --plan PLAN --plan-sha256 HEX [--plan PLAN2 --plan-sha256 HEX2]...
                             [--tree PATH] [--check | --report] [--report-path PATH]
                             [--choose [PLAN_SHA256:]INDEX=JOB]... [--expect-incidence THETAS:PHIS]
                             [--allow-unqualified [PLAN_SHA256:]INDEX]... [--include-job JOB]...

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
  record's or the manifest's commit not the plan's, or the manifest's ``environment_state`` absent, dirty or at another
  commit; the sidecar's incidence not the expected one; an acceptance other
  than ``array_run``; or the gate refusing the report: ``acceptance_gate.report_problems`` with the caller's
  ``PlanContext`` (the plan, its sha256, the task index and the policy) and ``acceptance_gate.convergence_problems``
  at the policy's polarizations. The sidecar's ``dataset_id``, ``frequency_label`` and ``frequency_ghz`` and each
  CSV's first ``Freq`` token are not compared with the task here: they stand on the run's own ``array_run`` identity
  and table rows, which the gate binds to these bytes, and the sidecar is found by the name the plan's dataset and the
  task's frequency give, so a misnamed one is missing;
- ``failed``: a deterministic FAIL, a well-formed bound report and ``done.json`` recording it, listed with its failed
  rows; its binding and context are judged as the gate judges a passed report, with its own verdict set aside, so a
  FAIL is told from a report that no longer binds these bytes; a FAIL is never accepted, except:
- ``valid``: complete, and the gate accepts the report. For a task named by ``--allow-unqualified`` a FAIL whose only
  failed row is the waivable ``convergence`` row, beside passed ``convergence.records`` and otherwise usable
  convergence records, is valid with ``waived = "convergence"`` (``report_problems(..., waive={"convergence"})`` and
  ``convergence_problems(..., allow_unqualified=True)``); the report, ``done.json`` and the Slurm state are never
  rewritten, and the override is recorded. Where a qualified run of the task passed, the override applies to nothing:
  the FAIL counts as a FAIL beside that PASS, and the report says so.

A task with exactly one valid candidate is selected (a FAIL beside it is no ambiguity); several valid ones are
ambiguous until ``--choose`` names one; with none it is ``failed`` when a deterministic FAIL exists, else ``missing``.

``report`` prints, and writes as ``<RUN>/assembly_<utc>.json`` (the first plan's run directory, or ``--report-path``),
every candidate with its verdict and reasons and every foreign entry; for a task with no valid run and no FAIL its
``sacct -j <A>_<index>`` lines (``submissions.log``'s arrays holding it) and the guarded rerun block from the plan's
frozen code (``python -m bbsim.array_submit rerun`` with the class's resources, then a new assembly after the reruns,
any pending assembly of the plan's cancelled first), and the same block with twice the class's time (TIMEOUT) and twice
its memory (OUT_OF_MEMORY), never more cores. One block of the three is pasted, once every listed task has ended (the
``sacct -X`` line printed above them, when an array holds them): each makes ``<RUN>/rerun_<report stamp>`` first with
``mkdir``, so a second paste of any block of the report stops before submitting (one marker per plan; each plan's blocks
rerun only its own tasks). A plan with tasks still unresolved and nothing to rerun, a plan whose tasks this report
resolved with a ``--choose`` or the convergence override, and a production plan whose every task has one valid run (its
assembly job may have stopped without writing the tree) get the guarded re-assembly instead, except inside a logged
assembly job of the plan (the job those lines would submit); it makes ``<RUN>/reassemble_<report stamp>`` first in the
same way. Every assembly these lines submit carries the report's ``--choose`` and applied ``--allow-unqualified`` of
that plan and its ``--tree`` (``array_submit``'s assembly arguments), so the guarded job assembles what the report
selected. No guarded assembly is printed (``_guard_blocker``; the rerun blocks then rerun and submit none, and the
re-assembly line is left out with the reason) for a report of several plans, whose tree is assembled by hand from the
printed ``--check`` and assembly lines of the report's plans and ``--tree`` (no assembly line where a check would refuse
or the tree checks failed); while the tree or its ``.partial`` is in the way; or when ``--check``'s tree checks failed.
``--choose`` lines for an ambiguous task and an ``--allow-unqualified`` line for a FAIL the override could accept, which
it never applies unasked (nor where a qualified run of the task passed), are ``--check`` lines with the report's
``--tree``, or ``--report`` lines whenever a ``--check`` of the plans would refuse (a probe plan, plans at different
commits or worktrees off them, a ``--tree`` the report was not given, an existing tree or ``<TREE>.partial``); they
carry the overrides that accepted a run of their task. For every array of the plan ``sacct``, the ``du -sh`` of each job
and its ``.aedtresults``, and the hosts' ``TempDirectory`` lines. Every line holds real values (paths, digests, indices,
job ids), so it is pasted as printed. Its last line is ``REPORT PASS``/``REPORT FAIL: ...`` or ``CHECK PASS``/``CHECK
FAIL: ...``.

``assemble`` writes the tree once every task has one valid run (HF-057, S65, G4-012). The include jobs are those the
plans recorded at plan time (``include_records``; a job two plans record must carry one record, and ``--include-job``
must name exactly them). Before anything is created it checks, and any failure is an ``ASSEMBLY FAIL`` with no tree:

- one schema major across every sidecar going into the tree (S65; a probe's schema-1 output never enters one);
- every include job re-validated exactly as at plan time (``array_plan._include_job``: schema 2, another dataset ID, the
  plan's incidence, qualified, its acceptance through the gate with its bound evidence and, for an ``array_run`` job,
  its own production plan), and the record that gives equal to the plan's, every digest included; its seven files
  (sidecar, two CSVs and the manifest copy per polarization) are the record's;
- across one plan's runs, each manifest's ``config`` equal, field by field, to ``array_plan.task_config`` of the plan's
  frozen config copy for that task, so the runs differ only in frequency and the resolution the plan prescribed, and
  one AEDT and one PyAEDT version;
- per dataset ID, every sidecar's ``tree_invariant.invariant`` agreeing with the others' on ``INVARIANT_BLOCKS``
  (BBRsim's reader's projection and comparison, booleans told from numbers), and no frequency twice; a WARN when one ID
  mixes AEDT versions.

Then it makes ``<TREE>.partial`` with an exclusive ``mkdir`` (its parents as needed), copies per dataset into
``waveguides/`` the sidecar and, per polarization, ``waveguide.csv``, ``far_field.csv`` and ``manifest.json``, byte for
byte (never ``project/``, ``scratch/`` or logs), re-hashes every copy and compares it with the digest it was verified at
(the CSVs' bytes and rows with the sidecar's too), writes ``tree_manifest.json`` (``TREE_SCHEMA``) and ``SHA256SUMS``
(``sha256sum -c`` format, the paths under ``waveguides/`` and ``tree_manifest.json`` in C order), checks ``<TREE>`` is
still absent and renames ``.partial`` to it. Any failure after the ``mkdir``, a SIGTERM included, removes ``.partial``
and the parents made for it; what cannot be removed is named, ``.partial`` with its paste-ready ``rm -rf`` line. Its
last line is ``ASSEMBLY PASS: <n> datasets + <m> included -> <TREE>; <id>: <f>, ... GHz x <keys> keys`` per ID, or
``ASSEMBLY FAIL: ...``.

The CLI: ``--report`` reports any plan, a probe's too; ``--check`` is the go/no-go before an assembly: the report,
then, once every task has one valid run, every check of ``assemble`` that writes nothing (``_checked_inputs``), whose
problems make it a ``CHECK FAIL`` (exit 1), all under the assembly's refusals: a probe
plan, plans at different commits or an assembler not imported from one plan's worktree (``assembler_problem``, as the
runner requires of a task), a plan's worktree ``<RUN>/code`` off the plan's commit or with uncommitted changes to
tracked files (``code_state``), an existing tree, or a ``<TREE>.partial`` (refused with the paste-ready ``squeue``
line that shows no assembly job of the plans runs, and the ``rm -rf`` line); ``--report`` prints those as ``WARN:``
lines. The assembly proper (neither option) takes the same refusals, prints the report, and writes the tree when every
task has one valid run. ``--verify`` is not available in this version. Exit codes: 2 refused before anything was
judged or created, 1 a task without exactly one valid run, a failed check or a failure while writing (the
``.partial`` removed), or a report that could not be written, 0 success.

Imports: the standard library, ``bbsim.naming``, ``bbsim.manifest`` (standard library; the checkout's ``REPO_ROOT``
and git timeout), ``bbsim.tree_invariant`` (standard library; its names are this module's too) and
``bbsim.acceptance_gate`` at module level (the binding, the required acceptance and the convergence gate, one
implementation shared with the bridge's ``assemble_tree.sbatch``; HF-057); ``acceptance_gate.report_problems`` imports
the checkers (numpy, pandas), so the selection runs in the bb venv. ``bbsim.array_plan`` and ``bbsim.config``
(standard library at import) are imported where they are used.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
import threading
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from bbsim import acceptance_gate as gate
from bbsim.manifest import GIT_TIMEOUT_S, REPO_ROOT
from bbsim.naming import dataset_dir_name, dataset_record_name, frequency_label
from bbsim.tree_invariant import (INVARIANT_BLOCKS, InvariantError, close, invariant,  # noqa: F401 (this module's too)
                                  invariant_diff, schema_major)

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
TREE_SCHEMA = "bbsim-tree/2"
WAVEGUIDES, TREE_MANIFEST, SUMS = "waveguides", "tree_manifest.json", "SHA256SUMS"
EXIT_SIGTERM = 143

Key = tuple[str, int]          # (plan sha256, task index)


class AssembleRefused(ValueError):
    """Refused before anything was judged or created (exit 2); the message names what and why, and ``lines`` are
    paste-ready commands printed after it (for a leftover ``<TREE>.partial``)."""

    def __init__(self, message: str, lines: Sequence[str] = ()):
        super().__init__(message)
        self.lines = list(lines)


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
    the override could accept (the report prints its line; it is never applied unasked); ``failed_rows`` are a FAIL's
    failed rows."""
    job: Path
    index: int
    verdict: str
    reasons: list[str]
    waived: str | None = None
    waivable: bool = False
    failed_rows: list[str] = field(default_factory=list)
    records: dict = field(default_factory=dict, repr=False)   # a valid run's start, done, manifest, sidecar, report


@dataclass
class Selection:
    """``select``'s result. Keys are (plan sha256, task index). ``states`` is each task's ``selected``, ``missing``,
    ``ambiguous``, ``failed`` or ``unresolved`` (a ``--choose`` naming a run that is not valid), ``unresolved`` the
    reason a task has no chosen run; ``foreign_reasons`` why each foreign entry is one; ``notes`` warnings the CLI
    adds (a plan's worktree, or the assembler's checkout, not the plans' code under ``--report``), ``note_lines``
    the paste-ready commands printed after them; ``tree`` the tree the printed ``--check`` lines name (``--tree``, or
    the one the plans record); ``override_notes`` why a task named by ``--allow-unqualified`` took no override;
    ``check_refusals`` the refusals a ``--check`` of these plans would give (the CLI's), so printed lines are
    ``--report`` lines then; ``inputs`` the one capture every file was read through, which the assembly reads on;
    ``tree_blocked`` why the tree cannot be written (it or its ``.partial`` exists), so no guarded assembly is printed;
    ``tree_checks`` ``--check``'s tree checks {problems, warnings}, None where they did not run."""
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
    tree: Path | None = None
    override_notes: dict[Key, str] = field(default_factory=dict)
    check_refusals: list[str] = field(default_factory=list)
    note_lines: list[str] = field(default_factory=list)
    inputs: gate.InputCapture | None = None
    tree_blocked: str | None = None
    tree_checks: dict | None = None

    def plan(self, sha256: str) -> LoadedPlan:
        return next(p for p in self.plans if p.sha256 == sha256)

    @property
    def applied_overrides(self) -> frozenset[Key]:
        """The tasks with a run the convergence override accepted (OR-027), selected or not: a printed line of an
        ambiguous task of such runs carries the override, or its --choose would name a FAIL."""
        return frozenset(key for key, cands in self.candidates.items() if any(c.waived for c in cands))


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


def _samefile(a: Path, b: Path) -> bool:
    try:
        return os.path.samefile(a, b)
    except OSError:
        return False


def assembler_problem(plans: Sequence[LoadedPlan]) -> str | None:
    """Why this assembler may not judge ``plans`` as their assembly would, or None: the plans are at one commit, and
    the checkout this module was imported from (``REPO_ROOT``) is one plan's worktree ``<RUN>/code``, as the runner and
    ``array_submit job-field`` require of every task (HF-057); ``worktree_problem`` then holds that worktree to the
    commit, clean. Another checkout's gate, rows and tolerances are not the ones the runs were judged by."""
    commits = sorted({p.data["code"]["commit"] for p in plans})
    if len(commits) > 1:
        listing = "; ".join(f"{p.path} at {p.data['code']['commit']}" for p in plans)
        return (f"the plans are at different commits ({listing}): one assembler judges the plans of one commit, from "
                "that commit's worktree")
    if not any(_samefile(p.code_repo, REPO_ROOT) for p in plans):
        return (f"the assembler runs from {REPO_ROOT}, which is no plan's worktree "
                f"({', '.join(sorted(str(p.code_repo) for p in plans))}): run it from the plan's <RUN>/code, the code "
                "every task ran")
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
    """The commit the start record and the manifest name is the plan's, the manifest records its code clean at it (an
    absent record fails), and the sidecar's incidence is the expected one."""
    problems = []
    commit = plan.data["code"]["commit"]
    code = start.get("code") if isinstance(start.get("code"), Mapping) else {}
    if code.get("commit") != commit:
        problems.append(f"its {START_NAME} records the commit {code.get('commit')!r}, not the plan's {commit}")
    if manifest.get("git_commit") != commit:
        problems.append(f"its manifest records git_commit {manifest.get('git_commit')!r}, not the plan's {commit}")
    state = manifest.get("environment_state")
    if not isinstance(state, Mapping):
        problems.append(f"its manifest records no environment_state {{commit, dirty, pyaedt_version}} (got "
                        f"{_show(state)}), so the state of the code it ran is unknown")
    else:
        if state.get("dirty") is not False:
            problems.append(f"its manifest records its code dirty ({_show(state.get('dirty'))}), not clean")
        if state.get("commit") != commit:
            problems.append(f"its manifest's environment_state names the commit {state.get('commit')!r}, not the "
                            f"plan's {commit}")
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
        done, manifest, sidecar, report = _completion(plan, task, job, start, inputs)
        records = {"start": start, "done": done, "manifest": manifest, "sidecar": sidecar, "report": report}
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
            return Candidate(job, index, VALID, [], records=records)
        failed = _failed_rows(report)
        waived = problems_of(report, WAIVE, True)
        if waiver and not waived:
            return Candidate(job, index, VALID, [f"accepted with the convergence override (--allow-unqualified, a "
                                                 f"diagnostic scope, OR-027): its {REPORT_NAME} failed {failed}; the "
                                                 f"report and {DONE_NAME} are unchanged"], waived="convergence",
                             waivable=True, failed_rows=failed, records=records)
        note = ([f"the convergence override does not apply: {'; '.join(waived)}"] if waiver else [])
        problems = gate.report_problems(_verdict_set_aside(report), KIND, job, requirement, inputs, context=context)
        if problems:
            raise _Stop(INVALID, problems + note)
        reasons = [f"its {REPORT_NAME} is a deterministic FAIL: {row}" for row in failed] or [
            f"its {REPORT_NAME} records passed false"]
        return Candidate(job, index, FAILED, reasons + note, waivable=not waived, failed_rows=failed)
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
    task to the job ``--choose`` named; ``allow_unqualified`` the tasks the convergence override may accept, where no
    qualified run passed (a qualified PASS is preferred, and the FAIL beside it counts as a FAIL);
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
                    expect_incidence=None if expect_incidence is None else dict(expect_incidence), inputs=inputs)
    for plan in plans:
        found, foreign = _scan(plan, inputs)
        sel.foreign += sorted(foreign)
        sel.foreign_reasons.update(foreign)
        for task in plan.tasks():
            key = (plan.sha256, task["index"])
            sel.candidates[key] = [judge(plan, task, job, start, inputs, waiver=key in allow, expect=expect_incidence)
                                   for job, start in found[task["index"]]]
    for key, candidates in sel.candidates.items():
        if key in allow and any(c.verdict == VALID and c.waived is None for c in candidates):
            # A qualified run of the task passed, so the override applies to nothing: a FAIL it would have accepted
            # counts as the FAIL it is, beside the PASS (OR-027, HF-057)
            for cand in candidates:
                if cand.waived:
                    cand.verdict, cand.waived, cand.waivable = FAILED, None, False
                    cand.reasons = [f"its {REPORT_NAME} is a deterministic FAIL: {row}" for row in cand.failed_rows] + [
                        "the convergence override was not applied to it: a qualified run of this task passed"]
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
    for key in sorted(allow - sel.applied_overrides):
        sel.override_notes[key] = (
            "a qualified run of the task passed, and is preferred to any FAIL the override would accept"
            if any(c.verdict == VALID for c in sel.candidates[key]) else
            "no run of the task failed its convergence row alone, on otherwise usable evidence")
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


def rerun_marker(plan: LoadedPlan, stamp: str) -> Path:
    """The directory a rerun block of the report ``stamp`` makes first, ``<RUN>/rerun_<stamp>``: an exclusive
    ``mkdir``, so one block of a report's three for a plan is submitted, once."""
    return plan.run_dir / f"rerun_{stamp}"


def reassembly_marker(plan: LoadedPlan, stamp: str) -> Path:
    """The directory the re-assembly line of the report ``stamp`` makes first, ``<RUN>/reassemble_<stamp>``, so that
    line submits one assembly, once."""
    return plan.run_dir / f"reassemble_{stamp}"


def assembly_args(sel: Selection, plan: LoadedPlan) -> list[str]:
    """What a guarded assembly of ``plan`` takes from this report (``array_submit``'s assembly arguments): the
    ``--choose`` that selected a run of the plan's tasks, the ``--allow-unqualified`` that accepted one (OR-027) and,
    for a one-plan report of a production plan, its ``--tree``. A plan's assembly job assembles that plan alone, so a
    tree of several plans is not handed to it."""
    args: list[str] = []
    for (sha, index), job in sorted(sel.choose.items()):
        if sha == plan.sha256 and sel.states.get((sha, index)) == SELECTED:
            args += ["--choose", f"{index}={job}"]
    for sha, index in sorted(sel.applied_overrides):
        if sha == plan.sha256:
            args += ["--allow-unqualified", f"{sha}:{index}"]
    if len(sel.plans) == 1 and plan.policy == "production" and sel.tree is not None:
        args += ["--tree", str(sel.tree)]
    return args


def _quoted(args: Sequence[str]) -> str:
    return "".join(f" {shlex.quote(a)}" for a in args)


def _rerun_block(plan: LoadedPlan, by_class: Mapping[str, list[int]], assemblies: Sequence[str], variant: str,
                 current_job: str | None, marker: Path, extra: Sequence[str] | None) -> list[str]:
    """One rerun block; with ``extra`` None (``_guard_blocker``) it reruns and submits no assembly."""
    q = shlex.quote
    plan_arg, sha = q(str(plan.path)), plan.sha256
    classes = {c["name"]: c for c in plan.data["classes"]}
    stop = (f"STOP: a block of this report was already pasted ({marker} exists, or cannot be made), so nothing was "
            "submitted now: make a new report for what remains")
    lines = [f"cd {q(str(plan.code_repo))} && (",
             'bbsim_stop() { echo "STOP: $1 was not submitted (exit $2); what was submitted before it keeps running, '
             'and nothing after it was submitted."; exit 1; }',
             f"mkdir {q(str(marker))} || {{ echo {q(stop)}; exit 1; }}"]
    for job in assemblies if extra is not None else ():
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
    if extra is None:
        lines += ['echo "submitted: ' + " ".join(f"{n}=${n}" for n in names) + '"', ")"]
        return lines
    dependency = "afterany:" + ":".join(f"${n}" for n in names)
    lines.append(f"ASSEMBLY=$(python -m bbsim.array_submit submit {plan_arg} {sha} assembly "
                 f'--dependency "{dependency}"{_quoted(extra)}) || bbsim_stop "the new assembly" $?')
    lines += ['echo "submitted: ' + " ".join(f"{n}=${n}" for n in names) + ' assembly=$ASSEMBLY"', ")"]
    return lines


def _guard_blocker(sel: Selection) -> str | None:
    """Why this report prints no guarded assembly (a plan's assembly job), or None: a report of several plans (each
    such job assembles its plan alone into that plan's tree, so their tree is assembled by hand), a tree or
    ``<TREE>.partial`` in the way (the job would refuse), or tree checks that failed (it would fail them)."""
    if len(sel.plans) > 1:
        return ("this report holds several plans, and a plan's guarded assembly job assembles that plan alone, so "
                "their tree is assembled by hand (the lines below)")
    if sel.tree_blocked is not None:
        return f"{sel.tree_blocked}; clear that first (the lines above), then make a new report"
    if sel.tree_checks is not None and sel.tree_checks["problems"]:
        return "the tree checks failed (above), and an assembly would fail them too"
    return None


def _line_mode(sel: Selection) -> str:
    """The mode of a printed assembler line: ``--check`` (the assembly's refusals) unless a ``--check`` of these plans
    would refuse (``check_refusals``: a probe plan, the plans' commits or worktrees, no one tree and no ``--tree``, an
    existing tree or ``<TREE>.partial``); then ``--report``."""
    if sel.check_refusals or any(p.policy == "probe" for p in sel.plans) or sel.tree is None:
        return "--report"
    return "--check"


def _check_line(sel: Selection, *extra: str, mode: str | None = None) -> str:
    """An assembler line (``_line_mode``, or with ``mode="assembly"`` the assembly itself) from the first plan's code
    with the selection's plans, tree, choices and the overrides that applied, plus ``extra``."""
    q = shlex.quote
    several = len(sel.plans) > 1
    argv = []
    for plan in sel.plans:
        argv += ["--plan", q(str(plan.path)), "--plan-sha256", plan.sha256]
    mode = mode or _line_mode(sel)
    if mode != "assembly":
        argv.append(mode)
    if mode != "--report":
        argv += ["--tree", q(str(sel.tree))]
    for (sha, index), job in sorted(sel.choose.items()):
        argv += ["--choose", f"{sha + ':' if several else ''}{index}={q(str(job))}"]
    for sha, index in sorted(sel.applied_overrides):
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


def _lines(sel: Selection, mode: str, stamp: str) -> tuple[list[str], bool, str]:
    """The printed lines (without the verdict line) of the report ``stamp``, whether every task has one valid run, and
    the summary."""
    q = shlex.quote
    current_job = os.environ.get("SLURM_JOB_ID") or None
    keys = [(p.sha256, t["index"]) for p in sel.plans for t in p.tasks()]
    lines = [f"== bbsim.assemble --{mode}: {len(sel.plans)} plan(s), {len(keys)} task(s)"]
    for plan in sel.plans:
        lines.append(f"plan {plan.path} sha256 {plan.sha256}: {plan.policy}, {plan.data['dataset_id']}, "
                     f"{len(plan.tasks())} task(s), code {plan.code_repo} at {plan.data['code']['commit']}")
    lines += [f"WARN: {note}" for note in sel.notes]
    lines += sel.note_lines
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
        if key in sel.override_notes:
            lines.append(f"    NOTE: the convergence override (--allow-unqualified) of {_task_name(sel, key)} applied "
                         f"to nothing: {sel.override_notes[key]}")
    for path in sel.foreign:
        lines.append(f"foreign    {path}: {sel.foreign_reasons[path]}")

    if sel.tree_checks is not None:
        checks = sel.tree_checks
        lines += [f"WARN: {w}" for w in checks["warnings"]]
        lines += [f"FAIL: {_one_line(problem)}" for problem in checks["problems"]]
        verdict = f"FAIL, {len(checks['problems'])} problem(s)" if checks["problems"] else "PASS"
        lines.append(f"tree checks: {verdict} (one schema, the include jobs re-validated, each run's config its "
                     "task's, one AEDT and PyAEDT, BBRsim's cross-frequency blocks; nothing written)")
    unresolved = [k for k in keys if k not in sel.chosen]
    blocker = _guard_blocker(sel)
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
            marker = rerun_marker(plan, stamp)
            listed = [f"{a}_{key[1]}" for key in rerun for a in _arrays(records, cands, key[1])]
            when = ("only once sacct shows every listed task has ended (no PENDING, RUNNING, REQUEUED or SUSPENDED "
                    "State in the line below)" if listed else
                    "only once squeue shows none of these tasks queued (no array that submissions.log or the jobs "
                    "record holds them, so no sacct line can show it)")
            then = ("These blocks submit no assembly: " + blocker + "." if blocker else
                    "Each also cancels a pending assembly of the plan first, and submits a new one that waits for the "
                    "reruns.")
            lines.append(f"# Paste exactly ONE of the three blocks below, and {when}, in a shell that has sourced "
                         f"bb_env.sh. Each block first makes {marker}, so this report admits one block of this plan, "
                         "pasted once; each reruns these tasks through the licence guard (OR-045) from the frozen "
                         f"code. {then}")
            if listed:
                lines.append(f"sacct -X -n -j {','.join(listed)} -o JobID%18,State%14")
            for variant, title in (("", "Block 1, with each class's resources:"),
                                   ("TIMEOUT", "Block 2, instead of block 1 if sacct shows TIMEOUT: twice each class's "
                                               "time (raise --time, never --cpus-per-task):"),
                                   ("OUT_OF_MEMORY", "Block 3, instead of block 1 if sacct shows OUT_OF_MEMORY: twice "
                                                     "each class's memory (raise --mem, never --cpus-per-task):")):
                lines.append(f"# {title}")
                lines += _rerun_block(plan, by_class, assemblies, variant, current_job, marker,
                                      None if blocker else assembly_args(sel, plan))
            continue
        extra = assembly_args(sel, plan)
        open_tasks = any(k[0] == plan.sha256 for k in unresolved)
        # Inside a guarded assembly of the plan (its own --report too) the choices are that job's: no second one
        inside = current_job is not None and current_job in [r["id"] for r in records if r.get("kind") == "assembly"]
        decided = any(a in ("--choose", "--allow-unqualified") for a in extra)
        # A production plan whose tasks all have a run gets the line too: its assembly job may have stopped (TIMEOUT,
        # OUT_OF_MEMORY, a SIGKILL during its cleanup) without writing the tree
        ready = mode != "assembly" and not inside and (decided or plan.policy == "production")
        if (open_tasks or ready) and blocker is not None:
            if len(sel.plans) == 1:
                lines.append(f"# No guarded re-assembly is printed: {blocker}.")
        elif open_tasks or ready:
            assemblies = [r["id"] for r in records if r.get("kind") == "assembly" and r["id"] != current_job]
            marker = reassembly_marker(plan, stamp)
            when = ("once the tasks above are settled" if open_tasks else
                    "with this report's choices, which the assembly submitted before does not take" if decided else
                    "should the assembly submitted before not have written the tree (it stopped, timed out or ran out "
                    "of memory), pasted only once squeue shows no assembly of the plan running")
            lines.append(f"# A new assembly of the plan {plan.path} through the guard, {when}; it makes {marker} "
                         "first, so this report submits it once, and a pending assembly of the plan is cancelled first "
                         "(in a shell that has sourced bb_env.sh):")
            stop = (f"STOP: this report's re-assembly was already pasted ({marker} exists, or cannot be made), so "
                    "nothing was submitted now: make a new report")
            lines.append(f"cd {q(str(plan.code_repo))} && (mkdir {q(str(marker))} || {{ echo {q(stop)}; exit 1; }}; "
                         + "".join(f"squeue -h -j {a} -t PENDING -o %i 2>/dev/null | grep -qx {a} && scancel {a}; "
                                   for a in assemblies)
                         + f"python -m bbsim.array_submit submit {q(str(plan.path))} {plan.sha256} assembly"
                         + f"{_quoted(extra)})")
    several = len(sel.plans) > 1
    if several and mode != "assembly":
        lines.append(f"# The tree of these {len(sel.plans)} plans is assembled by hand, from the code of the first, "
                     "once every task has one valid run: check, then assemble (in a shell that has sourced bb_env.sh); "
                     "a plan's guarded assembly job assembles that plan alone:")
        lines.append(_check_line(sel))
        if sel.check_refusals:
            lines.append(f"# No assembly line: the assembly would refuse ({_one_line(sel.check_refusals[0])}).")
        elif sel.tree_checks is not None and sel.tree_checks["problems"]:
            lines.append("# No assembly line: the tree checks failed (above), and the assembly would fail them too.")
        else:
            lines.append(_check_line(sel, mode="assembly"))
    for key in unresolved:
        if sel.states[key] == AMBIGUOUS:
            lines.append(f"# {_task_name(sel, key)} has several valid runs: check with the one to assemble (in a shell "
                         "that has sourced bb_env.sh); that check's report prints the guarded assembly with the same "
                         "--choose:")
            for cand in sel.candidates[key]:
                if cand.verdict == VALID:
                    lines.append(_check_line(sel, "--choose", f"{key[0] + ':' if several else ''}{key[1]}="
                                                              f"{q(str(cand.job))}"))
        waivable = [c for c in sel.candidates[key] if c.verdict == FAILED and c.waivable]
        if waivable and not any(c.verdict == VALID for c in sel.candidates[key]):
            lines.append(f"# {_task_name(sel, key)} failed only its convergence row, on otherwise usable evidence: the "
                         "diagnostic override accepts it only when asked (OR-027), recording it and rewriting nothing "
                         "(in a shell that has sourced bb_env.sh); that check's report prints the guarded assembly "
                         "with the override:")
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
            "allow_unqualified_applied": [{"plan_sha256": s, "index": i} for s, i in sorted(sel.applied_overrides)],
            "override_notes": [{"plan_sha256": s, "index": i, "note": n} for (s, i), n in
                               sorted(sel.override_notes.items())],
            "tree": None if sel.tree is None else str(sel.tree), "tree_checks": sel.tree_checks,
            "expect_incidence": sel.expect_incidence, "tasks": tasks,
            "foreign": [{"path": str(p), "reason": sel.foreign_reasons[p]} for p in sel.foreign],
            "notes": list(sel.notes), "lines": list(lines)}


def report_stamp(now: datetime | None = None) -> str:
    """A report's UTC stamp, which names its JSON file and its rerun markers."""
    return (now or _now()).strftime("%Y%m%dT%H%M%S.%fZ")


def report_path(sel: Selection, now: datetime | None = None, *, stamp: str | None = None) -> Path:
    """``<RUN>/assembly_<stamp>.json`` in the first plan's run directory."""
    return sel.plans[0].run_dir / f"assembly_{stamp or report_stamp(now)}.json"


def report(selection: Selection, *, path: Path | None, mode: str = "report") -> int:
    """Print the report (the module docstring's lines) and write it as JSON to ``path`` (``report_path`` when None),
    through a temporary file and ``os.replace``; the last line printed is the verdict, ``<MODE> PASS`` or ``<MODE>
    FAIL: ...`` (for the assembly only a FAIL: a passing assembly goes on to print its own). Returns 0 when every task
    has exactly one valid run and the tree checks, where they ran (``Selection.tree_checks``, under ``--check``),
    found nothing; 1 otherwise or when the report cannot be written."""
    if mode not in MODES:
        raise ValueError(f"mode must be one of {', '.join(MODES)}, got {mode!r}")
    stamp = report_stamp()
    lines, passed, summary = _lines(selection, mode, stamp)
    word = mode.upper()
    n = sum(len(p.tasks()) for p in selection.plans)
    problems = selection.tree_checks["problems"] if selection.tree_checks is not None else []
    if not passed:
        verdict = f"{word} FAIL: {n - len(selection.chosen)} of {n} task(s) without exactly one valid run ({summary})"
    elif problems:
        passed = False
        verdict = (f"{word} FAIL: every task has one valid run, and the tree checks found {len(problems)} problem(s): "
                   f"{_one_line(problems[0])}")
    else:
        verdict = f"{word} PASS: {n} task(s) of {len(selection.plans)} plan(s), one valid run each" + (
            ", and the tree checks pass" if selection.tree_checks is not None else "")
    for line in lines:
        print(line)
    target = Path(path) if path is not None else report_path(selection, stamp=stamp)
    tmp = target.with_name(f".{target.name}.tmp")
    try:
        tmp.write_text(json.dumps(_record(selection, mode, lines, passed, verdict), indent=2, sort_keys=True) + "\n")
        os.replace(tmp, target)
    except OSError as exc:
        if tmp.exists():
            tmp.unlink()
        print(f"the report could not be written to {target} ({_one_line(exc)})")
        print(f"{word} FAIL: the report could not be written")
        return 1
    print(f"report: {target}")
    if mode != "assembly" or not passed:
        print(verdict)
    return 0 if passed else 1


# --- the tree ---------------------------------------------------------------------------------------------------------

INCLUDE_KEYS = ("job", "dataset_id", "frequency_ghz", "acceptance_kind", "files", "source_plan")


def partial_path(tree: Path) -> Path:
    """``<TREE>.partial``, where the tree is written before its rename."""
    return tree.with_name(tree.name + ".partial")


def tree_problem(tree: Path, plans: Sequence[LoadedPlan]) -> tuple[str, list[str]] | None:
    """Why the assembly cannot write ``tree``, or None: the tree exists, or ``<TREE>.partial`` does (left by an
    assembly that stopped part way, or one writing it now), with the paste-ready lines that settle it: the ``squeue``
    line that shows any assembly job of the plans, and the ``rm -rf`` line. Nothing is removed here."""
    partial = partial_path(tree)
    if os.path.lexists(tree):
        return f"{tree} exists; the assembly writes a new tree only, never into or over one", []
    if not os.path.lexists(partial):
        return None
    names = sorted({p.data["slurm"]["job_name"] + "_asm" for p in plans if isinstance(p.data.get("slurm"), Mapping)
                    and isinstance(p.data["slurm"].get("job_name"), str)})
    squeue = (f"squeue -h -n {shlex.quote(','.join(names))} -o '%i %T %j'" if names else
              "squeue -h -u \"$USER\" -o '%i %T %j'")
    return (f"{partial} exists: an assembly that stopped part way left it, or one is writing it now, and nothing is "
            "removed here. Remove it only once the squeue line below prints no assembly job of the plan(s) and no "
            "assembly run by hand is writing it, with the rm line below",
            [squeue, f"rm -rf {shlex.quote(str(partial))}"])


def include_records(plans: Sequence[LoadedPlan], include_jobs: Sequence[Path | str] | None = None) -> list[dict]:
    """The include jobs the plans recorded at plan time (HF-057, S65), once each, with the digests of the plans that
    record each (``recorded_by``); a job two plans record must carry one record. ``include_jobs`` (``--include-job``),
    when given, must name exactly those jobs: a job the plans did not record has no frozen record to be held to (it is
    planned with ``plan_hfss_array.py --include-job``), and a recorded one is not dropped here. Raises
    ``AssembleRefused``."""
    found: dict[str, dict] = {}
    by: dict[str, list[str]] = {}
    for plan in plans:
        records = plan.data.get("include_jobs", [])
        if not isinstance(records, list):
            raise AssembleRefused(f"the plan {plan.path}'s include_jobs is not a list ({_show(records)})")
        for record in records:
            if not (isinstance(record, Mapping) and set(INCLUDE_KEYS) <= set(record)
                    and isinstance(record["job"], str) and os.path.isabs(record["job"])
                    and isinstance(record["dataset_id"], str) and isinstance(record["files"], Mapping)
                    and _positive(record["frequency_ghz"])
                    and (record["source_plan"] is None or isinstance(record["source_plan"], Mapping))):
                raise AssembleRefused(f"the plan {plan.path} records an include job that is not a record "
                                      f"{{{', '.join(INCLUDE_KEYS)}}}: {_show(record)}")
            key = str(Path(record["job"]).resolve())
            mine = {k: record[k] for k in INCLUDE_KEYS}
            if key in found and found[key] != mine:
                raise AssembleRefused(f"the include job {key} is recorded by the plans {by[key]} and {plan.sha256} "
                                      "with different records; one job has one frozen record")
            found.setdefault(key, mine)
            by.setdefault(key, []).append(plan.sha256)
    if include_jobs is not None:
        given = sorted({str(Path(os.path.abspath(job)).resolve()) for job in include_jobs})
        if given != sorted(found):
            raise AssembleRefused(f"--include-job names {given}, and the plans recorded {sorted(found)}: the assembly "
                                  "takes exactly the include jobs the plans recorded, each held to its frozen "
                                  "record; a job not recorded has none (plan it with plan_hfss_array.py --include-job)")
    return [dict(found[key], recorded_by=by[key]) for key in found]


def _positive(value: object) -> bool:
    return not isinstance(value, bool) and isinstance(value, (int, float)) and math.isfinite(value) and value > 0


@dataclass
class TreeInput:
    """One dataset going into the tree, a selected run of a plan task (``key``) or an include job: every file copied,
    by its path relative to the job, which is its path under ``waveguides/`` too, with the sha256 it was verified at;
    the bytes and rows its sidecar records for each CSV; and its ``tree_manifest.json`` entry."""
    kind: str
    job: Path
    dataset_id: str
    frequency_ghz: float
    label: str
    sidecar: dict
    files: dict[str, str]
    csv: dict[str, tuple[int, int]]
    entry: dict
    key: Key | None = None
    manifest: dict | None = None

    @property
    def stem(self) -> str:
        return f"{self.dataset_id}_{self.label}"


class _Failed(Exception):
    def __init__(self, problems: Sequence[str]):
        super().__init__(problems[0])
        self.problems = list(problems)


def _sidecar_csv(sidecar: Mapping, rel: str) -> tuple[str, int, int]:
    files = sidecar.get("files") if isinstance(sidecar.get("files"), Mapping) else {}
    entry = files.get(rel)
    if not (isinstance(entry, Mapping) and isinstance(entry.get("sha256"), str) and type(entry.get("bytes")) is int
            and type(entry.get("rows")) is int):
        raise ValueError(f"its sidecar records no sha256, bytes and rows for {rel} (got {_show(entry)}), which BBRsim "
                         "checks at load (G4-011)")
    return entry["sha256"], entry["bytes"], entry["rows"]


def _run_input(plan: LoadedPlan, task: Mapping, cand: Candidate,
               inputs: gate.InputCapture) -> tuple[TreeInput, list[str]]:
    """A selected run's files, at the digests its completion record and sidecar give, and its manifest entry."""
    records = cand.records
    done, manifest, sidecar, start = records["done"], records["manifest"], records["sidecar"], records["start"]
    dataset_id, frequency, index = plan.data["dataset_id"], task["frequency_ghz"], task["index"]
    pols = gate.ARRAY_POLICIES[plan.policy]
    problems: list[str] = []
    files = {dataset_record_name(dataset_id, frequency): done["sidecar_sha256"]}
    csv: dict[str, tuple[int, int]] = {}
    for e in pols:
        directory = dataset_dir_name(dataset_id, frequency, e)
        files[f"{directory}/manifest.json"] = done["manifest_sha256"]
        for name in CSV_NAMES:
            rel = f"{directory}/{name}"
            try:
                digest, size, rows = _sidecar_csv(sidecar, rel)
            except ValueError as exc:
                problems.append(f"{cand.job}: {exc}")
                continue
            files[rel], csv[rel] = digest, (size, rows)
    convergence = {}
    for e in pols:
        path = cand.job / "logs" / f"convergence_Ephi{e}.txt"
        try:
            text = inputs.read_text(path)
        except (OSError, ValueError) as exc:
            problems.append(f"{path}: its convergence table cannot be read as it was judged ({_one_line(exc)})")
            continue
        lines = [line for line in text.splitlines() if line.strip()]
        convergence[str(e)] = {"last_line": lines[-1] if lines else None, "sha256": inputs.recorded(path)}
    slurm = start.get("slurm") if isinstance(start.get("slurm"), Mapping) else {}
    versions = manifest.get("versions") if isinstance(manifest.get("versions"), Mapping) else {}
    if not all(isinstance(versions.get(k), str) and versions[k] for k in ("aedt", "pyaedt")):
        problems.append(f"{cand.job}: its manifest records no AEDT and PyAEDT versions (versions "
                        f"{_show(manifest.get('versions'))}), which the runs of one plan must share and show")
    project = manifest.get("project") if isinstance(manifest.get("project"), Mapping) else {}
    exit_field = sidecar.get("exit_field") if isinstance(sidecar.get("exit_field"), Mapping) else {}
    entry = {"plan_sha256": plan.sha256, "index": index, "frequency_ghz": frequency, "label": task["label"],
             "class": _class_of(plan, index), "stem": plan.stem(task), "job": str(cand.job),
             "job_id": start.get("job_id"), "array_job_id": slurm.get("SLURM_ARRAY_JOB_ID"),
             "slurm_job_id": slurm.get("SLURM_JOB_ID"), "aedt_version": versions.get("aedt"),
             "pyaedt_version": versions.get("pyaedt"), "project_sha256": project.get("job_copy_sha256_before_run"),
             "config_sha256": plan.data["config"]["sha256"] if isinstance(plan.data.get("config"), Mapping) else None,
             "sidecar_sha256": done["sidecar_sha256"], "manifest_sha256": done["manifest_sha256"],
             "report_sha256": done["report_sha256"], "l_resolution_mm": task.get("l_resolution_mm"),
             "effective_l_resolution_mm": task.get("effective_l_resolution_mm"), "l_step_mm": task.get("l_step_mm"),
             "exit_grid": exit_field.get("grid"), "exit_points_per_key": exit_field.get("points_per_key_retained"),
             "convergence": convergence, "accepted_with_override": cand.waived is not None}
    item = TreeInput("run", cand.job, dataset_id, float(frequency), task["label"], sidecar, files, csv, entry,
                     key=(plan.sha256, index), manifest=manifest)
    return item, problems


def _include_input(record: Mapping, inputs: gate.InputCapture) -> tuple[TreeInput | None, list[str]]:
    """An include job's seven files at the digests the plan recorded, read from the sidecar its record names."""
    job = Path(record["job"])
    where = f"--include-job {job}"
    sidecars = sorted(job.glob("*.dataset.json"))
    if len(sidecars) != 1:
        return None, [f"{where} holds {len(sidecars)} dataset sidecars (*.dataset.json), not 1"]
    try:
        sidecar = inputs.parse(sidecars[0], json.load)
    except (OSError, ValueError) as exc:
        return None, [f"{where}: its sidecar {sidecars[0].name} cannot be read ({_one_line(exc)})"]
    if not isinstance(sidecar, dict):
        return None, [f"{where}: its sidecar {sidecars[0].name} is not a JSON object"]
    dataset_id, frequency = sidecar.get("dataset_id"), sidecar.get("frequency_ghz")
    if not (isinstance(dataset_id, str) and dataset_id and _positive(frequency)):
        return None, [f"{where}: its sidecar names dataset {dataset_id!r} at {frequency!r} GHz"]
    label, recorded, problems = frequency_label(frequency), record["files"], []
    name = dataset_record_name(dataset_id, frequency)
    if sidecars[0].name != name or sidecar.get("frequency_label") != label:
        problems.append(f"{where}: its sidecar is {sidecars[0].name} with label {sidecar.get('frequency_label')!r}, "
                        f"not {name} with label {label!r}")
    if dataset_id != record["dataset_id"] or float(frequency) != float(record["frequency_ghz"]):
        problems.append(f"{where}: its sidecar is {dataset_id} at {frequency:g} GHz, and the plan recorded "
                        f"{record['dataset_id']} at {record['frequency_ghz']:g} GHz")
    if recorded.get(sidecars[0].name) != inputs.recorded(sidecars[0]):
        problems.append(f"{where}: its sidecar has sha256 {inputs.recorded(sidecars[0])}, not the "
                        f"{recorded.get(sidecars[0].name)} the plan recorded: it changed after planning")
    files, csv = {name: recorded.get(sidecars[0].name)}, {}
    for e in gate.ARRAY_POLICIES["production"]:
        directory = dataset_dir_name(dataset_id, frequency, e)
        for file in gate.EPHI_FILES:
            rel = f"{directory}/{file}"
            if not isinstance(recorded.get(rel), str):
                problems.append(f"{where}: the plan's record of it holds no {rel}, which the tree takes")
                continue
            files[rel] = recorded[rel]
            if file in CSV_NAMES:
                try:
                    digest, size, rows = _sidecar_csv(sidecar, rel)
                except ValueError as exc:
                    problems.append(f"{where}: {exc}")
                    continue
                if digest != recorded[rel]:
                    problems.append(f"{where}: its sidecar records {rel} at sha256 {digest}, the plan {recorded[rel]}")
                csv[rel] = (size, rows)
    provenance = sidecar.get("provenance") if isinstance(sidecar.get("provenance"), Mapping) else {}
    kind = record["acceptance_kind"]
    report_name = gate.REPORT_NAMES.get(kind)
    entry = {"job": str(job), "dataset_id": dataset_id, "frequency_ghz": float(frequency), "label": label,
             "stem": f"{dataset_id}_{label}", "acceptance_kind": kind,
             "report_sha256": recorded.get(report_name) if report_name else None, "evidence": dict(recorded),
             "source_plan": record["source_plan"], "aedt_version": provenance.get("aedt_version"),
             "recorded_by": list(record.get("recorded_by", []))}
    return TreeInput("include", job, dataset_id, float(frequency), label, sidecar, files, csv, entry), problems


def _record_differences(frozen: Mapping, fresh: Mapping) -> list[str]:
    out = []
    for key in INCLUDE_KEYS:
        if key == "files":
            a, b = frozen.get("files", {}), fresh.get("files", {})
            for path in sorted(set(a) | set(b)):
                if a.get(path) != b.get(path):
                    out.append(f"{path} has sha256 {b.get(path)}, the plan recorded {a.get(path)}")
        elif frozen.get(key) != fresh.get(key):
            out.append(f"{key} is {_show(fresh.get(key))}, the plan recorded {_show(frozen.get(key))}")
    return out


def _revalidate(record: Mapping, plans: Mapping[str, LoadedPlan]) -> list[str]:
    """An include job re-validated exactly as at plan time, by each plan that records it (``array_plan._include_job``,
    the planner's own check, with the plan's dataset and incidence), and the record that gives equal to the frozen
    one."""
    from bbsim import array_plan

    problems = []
    for sha in record["recorded_by"]:
        plan = plans[sha]
        try:
            fresh = array_plan._include_job(record["job"], plan.data["dataset_id"], plan.data.get("expect_incidence"))
        except (OSError, ValueError, TypeError, KeyError) as exc:
            problems.append(f"{_one_line(exc)} (re-validated as at plan time, HF-057)")
            continue
        differences = _record_differences(record, fresh)
        if differences:
            shown = "; ".join(differences[:5]) + (f"; and {len(differences) - 5} more" if len(differences) > 5 else "")
            problems.append(f"--include-job {record['job']} changed after planning (the plan {plan.path} froze it): "
                            f"{shown}")
    return problems


def _json_differences(a: object, b: object, path: str) -> list[str]:
    """Where the JSON value ``a`` differs from ``b``, by dotted path; a boolean equals only a boolean."""
    if isinstance(a, Mapping) and isinstance(b, Mapping):
        out: list[str] = []
        for key in sorted(set(a) | set(b), key=str):
            if key not in a:
                out.append(f"{path}.{key} is absent, the task config gives {_show(b[key])}")
            elif key not in b:
                out.append(f"{path}.{key} = {_show(a[key])} is not in the task config")
            else:
                out += _json_differences(a[key], b[key], f"{path}.{key}")
        return out
    if isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        return [d for i, (x, y) in enumerate(zip(a, b)) for d in _json_differences(x, y, f"{path}[{i}]")]
    if isinstance(a, bool) or isinstance(b, bool):
        same = isinstance(a, bool) and isinstance(b, bool) and a == b
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)):
        same = a == b
    else:
        same = type(a) is type(b) and a == b
    return [] if same else [f"{path} = {_show(a)}, the task config gives {_show(b)}"]


def _plan_problems(plan: LoadedPlan, runs: Sequence[TreeInput], inputs: gate.InputCapture) -> list[str]:
    """Across one plan's runs: each manifest's config is its task's config of the plan's frozen copy
    (``array_plan.task_config``), so the runs differ only in frequency and the resolution the plan prescribes; one AEDT
    and one PyAEDT version."""
    import tomllib

    from bbsim import array_plan
    from bbsim.config import config_to_dict

    entry = plan.data.get("config")
    if not (isinstance(entry, Mapping) and isinstance(entry.get("copy"), str)):
        return [f"the plan {plan.path} records no config copy {{source, copy, sha256}}"]
    copy = Path(entry["copy"])
    try:
        raw = inputs.parse(copy, tomllib.load)
    except (OSError, ValueError) as exc:
        return [f"the plan {plan.path}'s config copy {copy} cannot be read ({_one_line(exc)})"]
    problems = []
    if inputs.recorded(copy) != entry.get("sha256"):
        problems.append(f"the plan's config copy {copy} has sha256 {inputs.recorded(copy)}, not the plan's "
                        f"{entry.get('sha256')}")
    for run in runs:
        index = run.key[1]
        try:
            expected = json.loads(json.dumps(config_to_dict(array_plan.task_config(plan.data, index, raw_config=raw))))
        except (OSError, ValueError) as exc:
            problems.append(f"task {index}'s config cannot be built from the plan's copy ({_one_line(exc)})")
            continue
        differences = _json_differences(run.manifest.get("config"), expected, "config")
        if differences:
            problems.append(f"{run.job}: its manifest's config is not task {index}'s config of the plan's frozen copy, "
                            "which differs between tasks only in frequency and the resolution the plan prescribes: "
                            + "; ".join(differences[:5]) + (f"; and {len(differences) - 5} more"
                                                            if len(differences) > 5 else ""))
    versions: dict[tuple, list[str]] = {}
    for run in runs:
        recorded = run.manifest.get("versions") if isinstance(run.manifest.get("versions"), Mapping) else {}
        versions.setdefault((recorded.get("aedt"), recorded.get("pyaedt")), []).append(run.label)
    if len(versions) > 1:
        listing = "; ".join(f"AEDT {a!r} and PyAEDT {b!r} at {', '.join(labels)}"
                            for (a, b), labels in versions.items())
        problems.append(f"the plan {plan.path}'s runs record more than one AEDT and PyAEDT version ({listing}); one "
                        "plan's runs are solved by one AEDT and one PyAEDT")
    return problems


def _aedt_version(sidecar: Mapping) -> str:
    provenance = sidecar.get("provenance")
    return str(provenance.get("aedt_version", "unknown")) if isinstance(provenance, Mapping) else "unknown"


def _dataset_problems(items: Sequence[TreeInput]) -> tuple[list[str], list[str]]:
    """Per dataset ID, BBRsim's cross-frequency invariant (``tree_invariant``, G4-012) of every sidecar agreeing with
    the first's, and no frequency twice; and the warnings: one ID mixing AEDT versions."""
    problems, warnings = [], []
    by_id: dict[str, list[TreeInput]] = {}
    for item in items:
        by_id.setdefault(item.dataset_id, []).append(item)
    for dataset_id, group in by_id.items():
        projected = []
        for item in group:
            try:
                projected.append((item, invariant(item.sidecar)))
            except InvariantError as exc:
                problems.append(f"{item.job}: its sidecar {item.stem}.dataset.json lacks what BBRsim's reader compares "
                                f"across frequencies (G4-012): {exc}")
        if len(projected) == len(group):
            first, first_invariant = projected[0]
            for item, other in projected[1:]:
                differ = invariant_diff(first_invariant, other)
                if differ:
                    problems.append(f"{dataset_id}: the sidecars of {first.stem} ({first.job}) and {item.stem} "
                                    f"({item.job}) differ in {differ}, which BBRsim's reader requires every frequency "
                                    "of one dataset to share (G4-012)")
        seen: dict[float, TreeInput] = {}
        for item in group:
            if item.frequency_ghz in seen:
                problems.append(f"{dataset_id}: {seen[item.frequency_ghz].job} and {item.job} are both "
                                f"{item.frequency_ghz:g} GHz; a dataset holds each frequency once")
            seen.setdefault(item.frequency_ghz, item)
        versions = sorted({_aedt_version(item.sidecar) for item in group})
        if len(versions) > 1:
            warnings.append(f"{dataset_id} mixes the AEDT versions {versions} in one lookup")
    return problems, warnings


def _checked_inputs(sel: Selection, records: Sequence[Mapping],
                    inputs: gate.InputCapture) -> tuple[list[TreeInput], list[str]]:
    """Every dataset going into the tree, checked (the module docstring's list, in three stages, each stopping the
    assembly with every problem it found); the warnings. Raises ``_Failed``."""
    items, problems = [], []
    for plan in sel.plans:
        for task in plan.tasks():
            item, found = _run_input(plan, task, sel.chosen[(plan.sha256, task["index"])], inputs)
            items.append(item)
            problems += found
    for record in records:
        item, found = _include_input(record, inputs)
        problems += found
        if item is not None:
            items.append(item)
    majors: dict[int, list[str]] = {}
    for item in items:
        try:
            majors.setdefault(schema_major(item.sidecar), []).append(str(item.job))
        except InvariantError as exc:
            problems.append(f"{item.job}: {exc}")
    if len(majors) > 1:
        listing = "; ".join(f"schema {major}.x from {jobs}" for major, jobs in sorted(majors.items()))
        problems.insert(0, f"the tree would mix schema majors ({listing}); a tree holds one schema (S65), and a "
                           "probe's schema-1 output never enters one")
    own = {p.data["dataset_id"] for p in sel.plans}
    for item in items:
        if item.kind == "include" and item.dataset_id in own:
            problems.append(f"--include-job {item.job}: its dataset ID {item.dataset_id} is a plan's own; an included "
                            "job brings another dataset ID into the tree")
    if problems:
        raise _Failed(problems)
    plans = {p.sha256: p for p in sel.plans}
    for record in records:
        problems += _revalidate(record, plans)
    if problems:
        raise _Failed(problems)
    for plan in sel.plans:
        problems += _plan_problems(plan, [i for i in items if i.key is not None and i.key[0] == plan.sha256], inputs)
    found, warnings = _dataset_problems(items)
    problems += found
    if problems:
        raise _Failed(problems)
    return items, warnings


def _copy_file(src: Path, dest: Path) -> None:
    shutil.copyfile(src, dest)


def _digest(path: Path, count_rows: bool) -> tuple[str, int, int | None]:
    """The sha256 and size of a file, and for a CSV its rows (non-blank lines after the header), in one pass."""
    hashed, size, rows = hashlib.sha256(), 0, 0
    with open(path, "rb") as fh:
        if count_rows:
            header = True
            for line in fh:
                hashed.update(line)
                size += len(line)
                if header:
                    header = False
                elif line.strip():
                    rows += 1
        else:
            while chunk := fh.read(gate.CHUNK):
                hashed.update(chunk)
                size += len(chunk)
    return hashed.hexdigest(), size, rows if count_rows else None


class _Sigterm:
    """SIGTERM while the tree is written: ``SystemExit(143)``, except while ``.partial`` is being removed, when it is
    held until the removal ends."""

    def __init__(self) -> None:
        self.cleaning = self.pending = False

    def __call__(self, signum: int, frame: object) -> None:
        if self.cleaning:
            self.pending = True
            return
        raise SystemExit(EXIT_SIGTERM)


def _write_tree(tree: Path, plans: Sequence[LoadedPlan], items: Sequence[TreeInput], record: Mapping,
                left: list[str]) -> int:
    """``<TREE>.partial`` by an exclusive ``mkdir``, the copies re-hashed against what was verified, the tree
    manifest, ``SHA256SUMS``, the rename; any failure removes ``.partial``. A SIGTERM handler is installed (in the main
    thread) before ``.partial`` is made and kept through its removal, which a SIGTERM never interrupts: the assembly
    then ends with exit 143. Returns the files ``SHA256SUMS`` lists. Raises ``AssembleRefused`` when another
    ``.partial`` appeared since the checks, ``_Failed`` for every failure, the tree's parent or ``.partial`` that
    cannot be made included. A failure removes ``.partial`` and the parent directories made for the tree; what cannot
    be removed is put in ``left`` and printed (``.partial`` with its paste-ready ``rm -rf`` line)."""
    partial = partial_path(tree)
    made_parents = []                                      # the tree's missing ancestors, deepest first
    ancestor = tree.parent
    while not os.path.lexists(ancestor) and ancestor != ancestor.parent:
        made_parents.append(ancestor)
        ancestor = ancestor.parent
    sigterm = _Sigterm()
    previous = signal.signal(signal.SIGTERM, sigterm) if threading.current_thread() is threading.main_thread() \
        else None
    made = renamed = False
    try:
        try:
            tree.parent.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            blocker = next((p for p in [tree.parent, *tree.parent.parents] if p.exists() and not p.is_dir()), None)
            raise _Failed([f"the tree's parent cannot be made: {blocker} is not a directory" if blocker is not None
                           else f"the tree's parent {tree.parent} cannot be made ({_one_line(exc)})"]) from None
        try:
            os.mkdir(partial)
        except FileExistsError:
            if os.path.lexists(partial):
                raise AssembleRefused(*(tree_problem(tree, plans) or (f"{partial} exists", []))) from None
            raise _Failed([f"{partial} cannot be made in {tree.parent}"]) from None
        except OSError as exc:
            raise _Failed([f"{partial} cannot be made ({_one_line(exc)}); nothing was written"]) from None
        made = True
        sums: dict[str, str] = {}
        waveguides = partial / WAVEGUIDES
        waveguides.mkdir()
        for item in items:
            for rel in sorted(item.files):
                src, dest = item.job / rel, waveguides / rel
                dest.parent.mkdir(exist_ok=True)
                _copy_file(src, dest)
                digest, size, rows = _digest(dest, rel in item.csv)
                if digest != item.files[rel]:
                    raise _Failed([f"the copy of {src} has sha256 {digest}, not the {item.files[rel]} it was verified "
                                   "at: it changed after it was judged, or the copy is not faithful"])
                if rel in item.csv and (size, rows) != item.csv[rel]:
                    raise _Failed([f"the copy of {src} has {size} bytes and {rows} rows, the sidecar records "
                                   f"{item.csv[rel][0]} and {item.csv[rel][1]}"])
                sums[f"{WAVEGUIDES}/{rel}"] = digest
        data = (json.dumps(record, indent=2, sort_keys=True) + "\n").encode("utf-8")
        with open(partial / TREE_MANIFEST, "xb") as fh:
            fh.write(data)
        sums[TREE_MANIFEST] = hashlib.sha256(data).hexdigest()
        with open(partial / SUMS, "x", encoding="utf-8") as fh:
            for path in sorted(sums, key=lambda p: p.encode("utf-8")):          # the order of LC_ALL=C sort
                fh.write(f"{sums[path]}  {path}\n")
        if os.path.lexists(tree):
            raise _Failed([f"{tree} appeared while the tree was written; nothing is renamed over it"])
        os.rename(partial, tree)
        renamed = True
        return len(sums)
    except OSError as exc:
        raise _Failed([f"writing {partial} failed: {_one_line(exc)}"]) from exc
    finally:
        sigterm.cleaning = True
        if not renamed:
            left += _remove_partial(partial if made else None, made_parents)
        if previous is not None:
            signal.signal(signal.SIGTERM, previous)
        if sigterm.pending and not renamed:            # a whole tree, renamed, stays a success
            raise SystemExit(EXIT_SIGTERM)


def _remove_partial(partial: Path | None, made_parents: Sequence[Path]) -> list[str]:
    """After a failure: remove ``partial`` (when this assembly made it) and the empty parents made for the tree,
    deepest first; print and return what is left (``.partial`` with its ``rm -rf`` line, to paste once no assembly
    runs; a parent someone else wrote into is named, never removed)."""
    left = []
    if partial is not None:
        shutil.rmtree(partial, ignore_errors=True)
        if os.path.lexists(partial):
            print(f"WARN: {partial} could not be removed; remove it once no assembly of the plan runs:")
            print(f"rm -rf {shlex.quote(str(partial))}")
            left.append(str(partial))
    for directory in made_parents:
        if not os.path.lexists(directory):
            continue
        try:
            os.rmdir(directory)
        except OSError:
            print(f"WARN: {directory}, made for the tree, is left: it is not empty")
            left.append(str(directory))
    return left


def _assembly_job() -> dict | None:
    """The guarded assembly job this runs in (its ``BBSIM_JOB`` context: the bridge commit and the script it ran), or
    None outside one."""
    text = os.environ.get("BBSIM_JOB")
    if not text:
        return None
    from bbsim import array_submit

    try:
        context = array_submit.decode_job(text)
    except array_submit.SubmitError as exc:
        return {"slurm_job_id": os.environ.get("SLURM_JOB_ID"), "error": str(exc)}
    return {"slurm_job_id": os.environ.get("SLURM_JOB_ID"),
            **{k: context[k] for k in ("kind", "plan_sha256", "script", "script_sha256", "bridge_commit",
                                       "extra_args")}}


def _tree_record(sel: Selection, tree: Path, items: Sequence[TreeInput]) -> dict:
    """``tree_manifest.json`` (``TREE_SCHEMA``)."""
    try:
        commit, dirty = code_state(REPO_ROOT)
        assembler: dict[str, Any] = {"repo": str(REPO_ROOT), "commit": commit, "dirty": dirty}
    except AssembleRefused as exc:
        assembler = {"repo": str(REPO_ROOT), "error": str(exc)}
    datasets: dict[str, dict] = {}
    for item in items:
        excitation = item.sidecar.get("excitation") or {}
        exit_field = item.sidecar.get("exit_field") or {}
        keys = len(excitation.get("incident_phi_deg", [])) * len(excitation.get("incident_theta_deg", []))
        record = datasets.setdefault(item.dataset_id, {
            "schema_major": schema_major(item.sidecar), "runs": 0, "included": 0, "frequencies_ghz": [],
            "incidence_keys": keys, "exit_points_per_key": {}})
        record["runs" if item.kind == "run" else "included"] += 1
        record["frequencies_ghz"].append(item.entry["frequency_ghz"])
        record["exit_points_per_key"][item.label] = exit_field.get("points_per_key_retained")
    for record in datasets.values():
        record["frequencies_ghz"].sort()
    waived = []
    for key, cand in sorted(sel.chosen.items()):
        if cand.waived:
            waived.append({"plan_sha256": key[0], "index": key[1], "failed_row": cand.waived,
                           "failed_rows": list(cand.failed_rows),
                           "report_sha256": cand.records["done"]["report_sha256"], "scope": "diagnostic"})
    return {
        "schema": TREE_SCHEMA, "created_utc": _now().isoformat(), "tree": str(tree),
        "plans": [{"path": str(p.path), "sha256": p.sha256, "commit": p.data["code"]["commit"], "policy": p.policy,
                   "dataset_id": p.data["dataset_id"], "run_dir": str(p.run_dir)} for p in sel.plans],
        "bridge": [{"plan_sha256": p.sha256, "repo": (p.data.get("bridge") or {}).get("repo"),
                    "commit": (p.data.get("bridge") or {}).get("commit")} for p in sel.plans],
        "assembly_job": _assembly_job(), "assembler": assembler, "datasets": datasets,
        "runs": [i.entry for i in items if i.kind == "run"],
        "include_jobs": [i.entry for i in items if i.kind == "include"],
        "overrides": {"allow_unqualified": waived,
                      "choose": [{"plan_sha256": s, "index": i, "job": str(j)}
                                 for (s, i), j in sorted(sel.choose.items())],
                      "expect_incidence": sel.expect_incidence,
                      "tree": str(tree) if any(p.data.get("tree") != str(tree) for p in sel.plans) else None},
    }


def _fail(problems: Sequence[str], what: str) -> int:
    for problem in problems:
        print(f"FAIL: {_one_line(problem)}")
    print(f"ASSEMBLY FAIL: {len(problems)} problem(s), {what}: {_one_line(problems[0])}")
    return 1


def assemble(selection: Selection, *, tree: Path, include_jobs: Sequence[Path]) -> int:
    """Write ``tree`` from the selection's runs and the plans' include jobs (``include_jobs`` must name exactly those
    the plans recorded), as the module docstring says. Returns 0 (``ASSEMBLY PASS``) or 1 (``ASSEMBLY FAIL``, no
    tree, no ``.partial``); raises ``AssembleRefused`` for a probe plan, an existing tree or ``.partial``, or include
    jobs other than the recorded."""
    tree = Path(os.path.abspath(tree))
    plans = selection.plans
    probes = [str(p.path) for p in plans if p.policy == "probe"]
    if probes:
        raise AssembleRefused(f"the probe plan(s) {probes} are diagnostics, never assembled into a tree (S65)")
    problem = tree_problem(tree, plans)
    if problem is not None:
        raise AssembleRefused(*problem)
    records = include_records(plans, include_jobs)
    keys = [(p.sha256, t["index"]) for p in plans for t in p.tasks()]
    open_tasks = [key for key in keys if key not in selection.chosen]
    if open_tasks:
        print(f"ASSEMBLY FAIL: {len(open_tasks)} of {len(keys)} task(s) without exactly one valid run; no tree written")
        return 1
    inputs = selection.inputs if selection.inputs is not None else gate.InputCapture()
    print(f"== bbsim.assemble: the tree {tree} from {len(keys)} run(s) of {len(plans)} plan(s) and {len(records)} "
          "included job(s)")
    try:
        items, warnings = _checked_inputs(selection, records, inputs)
    except _Failed as failed:
        return _fail(failed.problems, "no tree written")
    for warning in warnings:
        print(f"WARN: {warning}")
    for item in items:
        if item.kind == "include":
            print(f"include {item.job}: {item.dataset_id} at {item.frequency_ghz:g} GHz "
                  f"({item.entry['acceptance_kind']}), re-validated, its {len(item.entry['evidence'])} recorded files "
                  "unchanged since planning")
    left: list[str] = []
    try:
        listed = _write_tree(tree, plans, items, _tree_record(selection, tree, items), left)
    except _Failed as failed:
        return _fail(failed.problems, f"no tree written; left: {', '.join(left)} (the lines above)" if left else
                     "no tree written, and nothing of it left")
    print(f"copied {listed - 1} files into {tree / WAVEGUIDES}, each re-hashed against what was verified; {SUMS} lists "
          f"them and {TREE_MANIFEST}")
    summary = []
    for dataset_id in dict.fromkeys(item.dataset_id for item in items):
        group = [item for item in items if item.dataset_id == dataset_id]
        excitation = group[0].sidecar.get("excitation") or {}
        count = len(excitation.get("incident_phi_deg", [])) * len(excitation.get("incident_theta_deg", []))
        frequencies = ", ".join(f"{f:g}" for f in sorted(item.frequency_ghz for item in group))
        summary.append(f"{dataset_id}: {frequencies} GHz x {count} keys")
    runs = sum(item.kind == "run" for item in items)
    print(f"ASSEMBLY PASS: {runs} datasets + {len(items) - runs} included -> {tree}; {'; '.join(summary)}")
    return 0


# --- the CLI ----------------------------------------------------------------------------------------------------------

def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="assemble_hfss_dataset.py", allow_abbrev=False,
                                description="Select one accepted run per task of a job-array plan through the gate, "
                                            "report, and assemble the tree (HF-057, OR-027, CDX-038, S65).")
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
                   help="an include job the plans recorded (given, they must all be named; default the recorded ones)")
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


def parse_choices(plans: Sequence[LoadedPlan], choose: Sequence[str], allow_unqualified: Sequence[str], *,
                  absolute: bool = False) -> tuple[dict[Key, Path], set[Key]]:
    """``--choose [PLAN_SHA256:]INDEX=JOB`` and ``--allow-unqualified [PLAN_SHA256:]INDEX`` as task keys (``_key``). A
    job is resolved against the working directory, or with ``absolute`` must be an absolute path (a guarded job's
    working directory is not the submitter's). Raises ``AssembleRefused``."""
    chosen: dict[Key, Path] = {}
    for text in choose:
        index, sep, job = text.partition("=")
        if not sep or not job:
            raise AssembleRefused(f"--choose takes [PLAN_SHA256:]INDEX=JOB, got {text!r}")
        if absolute and not os.path.isabs(job):
            raise AssembleRefused(f"--choose {text}: the job must be an absolute path")
        path = Path(os.path.abspath(job)).resolve()
        key = _key(index, plans, "--choose", path)
        if key in chosen:
            raise AssembleRefused(f"--choose names task {key[1]} twice")
        chosen[key] = path
    return chosen, {_key(text, plans, "--allow-unqualified") for text in allow_unqualified}


def _tree(args: argparse.Namespace, plans: Sequence[LoadedPlan]) -> Path | None:
    """``--tree``, else the one absolute tree every plan records, else None."""
    if args.tree is not None:
        return Path(os.path.abspath(args.tree))
    trees = sorted({p.data.get("tree") for p in plans}, key=str)
    if len(trees) != 1 or not isinstance(trees[0], str) or not os.path.isabs(trees[0]):
        return None
    return Path(trees[0])


def _refused(exc: AssembleRefused) -> int:
    print(f"assemble: refused: {_one_line(exc)}", file=sys.stderr)
    for line in exc.lines:
        print(line, file=sys.stderr)
    return 2


def main(argv: list[str] | None = None) -> int:
    """The CLI (the module docstring's). Refusals print one ``refused:`` line on stderr, any paste-ready lines after it,
    and exit 2."""
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
        choose, allow = parse_choices(plans, args.choose, args.allow_unqualified)
        expect = None
        if args.expect_incidence is not None:
            from bbsim import array_plan

            try:
                expect = array_plan._expected_incidence(args.expect_incidence)
            except array_plan.PlanError as exc:
                raise AssembleRefused(str(exc)) from None
        records = include_records(plans, args.include_job or None)
        include = [Path(r["job"]) for r in records]
        tree = _tree(args, plans)
        problems = [assembler_problem(plans), *map(worktree_problem, plans)]
        notes = [problem for problem in problems if problem is not None]
        refusals: list[tuple[str, list[str]]] = []     # what a --check of these plans refuses, in order
        probes = [str(p.path) for p in plans if p.policy == "probe"]
        if probes:
            refusals.append((f"the probe plan(s) {probes} are diagnostics, never assembled: --report reports them", []))
        refusals += [(note, []) for note in notes]
        tree_lines: list[str] = []
        blocked = None
        if tree is None:
            refusals.append((f"the plans record the trees {sorted({str(p.data.get('tree')) for p in plans})}, not one "
                             "absolute tree: give --tree", []))
        else:
            found = tree_problem(tree, plans)
            if found is not None:
                refusals.append(found)
                notes.append(found[0])
                tree_lines = found[1]
                blocked = found[0]
        if mode != "report" and refusals:
            raise AssembleRefused(*refusals[0])
        if args.report_path is not None:
            target = Path(os.path.abspath(args.report_path))
            if os.path.lexists(target) or not target.parent.is_dir():
                raise AssembleRefused(f"--report-path {target} exists or has no directory; a report is a new file")
        else:
            target = None
        selection = select(plans, choose=choose, allow_unqualified=allow, expect_incidence=expect)
        selection.notes += notes
        selection.note_lines += tree_lines
        selection.tree = tree
        selection.tree_blocked = blocked
        selection.check_refusals = [message for message, _ in refusals]
    except AssembleRefused as exc:
        return _refused(exc)
    if mode == "check" and len(selection.chosen) == sum(len(p.tasks()) for p in plans):
        # The go/no-go before an assembly: every check that writes nothing (HF-057)
        try:
            _, warnings = _checked_inputs(selection, records, selection.inputs)
            selection.tree_checks = {"problems": [], "warnings": warnings}
        except _Failed as failed:
            selection.tree_checks = {"problems": failed.problems, "warnings": []}
    rc = report(selection, path=target, mode=mode)
    if mode != "assembly" or rc != 0:
        return rc
    try:
        return assemble(selection, tree=tree, include_jobs=include)
    except AssembleRefused as exc:
        return _refused(exc)


if __name__ == "__main__":
    sys.exit(main())
