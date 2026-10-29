"""Guarded submission of a job-array plan's jobs (HF-057, OR-044, OR-045): every array, rerun, licence monitor and
assembly job of a plan goes through ``licence_guard.Guard.submit``, its script the plan's recorded bridge blob given to
sbatch on stdin. Nothing here calls sbatch or scontrol directly.

``submit(plan, sha256, what, ...)`` checks everything before anything is reserved; each refusal is one ``SubmitError``
line (exit 2):

- the plan file's bytes have the given sha256 and schema ``bbsim-array-plan/2``, and it is its ``run_dir``'s
  ``plan.json``;
- the guard is the plan's: ``Guard.root`` is the plan's ``guard_dir``, and an inherited ``BBSIM_LICENCE_GUARD`` names
  the same directory (the guard checks both again; a second guard cannot bypass the cap);
- the plan's worktree ``<RUN>/code`` is clean at the plan's commit (``manifest.environment_state``, one seam);
- the plan's bridge is its code checkout's sibling clone ``../bbridge`` (OR-044), checked before anything is captured;
- an array: its class and indices are the plan's (a rerun's indices lie in one class, and a class already submitted is
  submitted again only as a rerun), its resources are well formed, the installed PyAEDT is the plan's, every pending or
  running ``bbsim:`` array runs at the plan's cores (the planner's rule, repeated before every array submission,
  HF-057; it runs outside the guard lock, so two submissions started at the same instant from separate shells can
  both pass it, and plans at different cores are submitted one after another), a production plan's frozen HF-041
  column-order reference still has its recorded sha256 and still fits the worktree's exporter files and the installed
  PyAEDT (``dataset_v2.check_reference_before_aedt``), and a probe plan records none;
  a monitor: the plan runs one, and an array of the plan is in ``submissions.log`` for it to watch;
- the plan's own sbatch arguments pass ``licence_guard.check_user_sbatch_args`` and set none of the resources the plan
  sets itself (``PLAN_RESOURCE_OPTIONS``); a dependency is ``after[any|ok|notok]:JOBID[:JOBID...]``;
- the script is the git blob of the recorded bridge commit with the recorded sha256 (``capture_script``): a commit or
  a blob that is no longer there, or other bytes, refuses naming the commit, and the live file is never read.

Then ``Guard.submit`` reserves (an array its class throttle, at most its task count; a monitor or an assembly no slot),
repeats admission under its lock and runs ``sbatch --parsable --export=ALL`` with the captured bytes on stdin and:
``--job-name`` (``bbsim_<name>``, with ``_mon`` or ``_asm`` for the monitor and the assembly), ``--comment`` (the plan's
``slurm.comment``, ``bbsim:<RUN>/plan.json``, which marks every job of the plan as managed),
``--array=INDICES%THROTTLE`` (arrays), ``--chdir=<RUN>/code``, ``--output`` under ``<RUN>/slurm``, the resources,
``--dependency`` with ``--kill-on-invalid-dep=yes`` when one is given, and for arrays the plan's own sbatch arguments.
No script path, plan path or context is in argv: the job gets its context in sbatch's environment, which
``--export=ALL`` passes on, ``BBSIM_CODE_DIR`` (``<RUN>/code``) and ``BBSIM_JOB`` (``encode_job``), beside the
guard's ``BBSIM_LICENCE_GUARD`` and ``BBSIM_GUARD_RESERVATION``. Each admitted submission appends one JSON line to
``<RUN>/submissions.log``: {kind, class, id, utc, script, script_sha256, bridge_commit, reservation, admitted (true),
plan_sha256, indices, rerun, dependency}.
``log_submission`` (the planner's ``--log-submission RUN KIND ID``) records a job submitted by hand with
``admitted: false``; it reserves nothing and grants no admission.

``BBSIM_JOB`` is unpadded base64url of compact UTF-8 JSON holding exactly ``JOB_KEYS``, so it matches
``^[A-Za-z0-9_-]+$`` and holds no ``=``, ``,``, space or quote whatever the values; ``plan`` and ``plan_sha256`` are
never null. A job reads one value with ``python -m bbsim.array_submit job-field NAME`` (``job_field``), run from
``BBSIM_CODE_DIR``: it checks the plan's bytes against the context's sha256, ``BBSIM_CODE_DIR`` (the plan's worktree,
the checkout this module was imported from, clean at the plan's commit), ``BBSIM_LICENCE_GUARD`` (the shared guard,
the context's and the plan's ``guard_dir``) and the context against the plan, and refuses a value holding a newline,
carriage return or NUL; a list prints one element per line, null an empty line.

CLI (``--root`` names the shared guard, by default ``licence_guard.DEFAULT_ROOT``): ``submit PLAN SHA256
{array,monitor,assembly} [--class NAME] [--dependency DEP]``; ``rerun PLAN SHA256 --indices I,J [--mem M] [--time T]
[--qos Q]``, nothing else; ``throttle PLAN SHA256 ARRAY_ID THROTTLE``, a live throttle change of an array of this plan
through ``Guard.change_throttle``, the only throttle change ``submissions.log`` records (``python -m
bbsim.licence_guard throttle`` changes it without a line there); ``remaining PLAN SHA256``, what remains of the plan's
printed block after some of its jobs are in ``submissions.log`` (``remaining_block``), its ids already admitted written
out, ready to paste, which every STOP of the block prints; ``job-field NAME``. ``submit`` and ``rerun`` print the job
id alone on stdout. Exit codes: 2 refused before anything was reserved or submitted, 1 failed, 0 success. I1's probe
and lock test are submitted by path through the guard and use none of this. Standard library at import;
``dataset_v2`` (numpy, pandas) only for a production array's reference check.
"""
from __future__ import annotations

import argparse
import base64
import binascii
import hashlib
import json
import os
import re
import shlex
import subprocess
import sys
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError
from pathlib import Path
from typing import Any

from bbsim import array_plan, licence_guard
from bbsim.manifest import GIT_TIMEOUT_S, REPO_ROOT, EnvironmentState, EnvironmentStateError, environment_state

CODE_DIR_ENV = "BBSIM_CODE_DIR"
JOB_ENV = "BBSIM_JOB"
JOB_KEYS = ("plan", "plan_sha256", "kind", "class", "script", "script_sha256", "bridge_commit", "extra_args",
            "code_commit", "guard_dir")
WHATS = ("array", "monitor", "assembly")
SCRIPTS = {"array": array_plan.ARRAY_SCRIPT, "monitor": array_plan.MONITOR_SCRIPT,
           "assembly": array_plan.ASSEMBLY_SCRIPT}
JOB_NAME_SUFFIX = {"array": "", "monitor": "_mon", "assembly": "_asm"}
SUBMISSIONS_LOG = "submissions.log"
MONITOR_MEM = "1G"
# The resources the plan sets on every job itself; a plan's own sbatch arguments may not set them a second time.
PLAN_RESOURCE_OPTIONS = ("--cpus-per-task", "--mem", "--time", "--qos")

_HEX64 = re.compile(r"[0-9a-f]{64}")
_COMMIT = re.compile(r"[0-9a-f]{40}|[0-9a-f]{64}")
_JOB_TEXT = re.compile(r"[A-Za-z0-9_-]+")
_JOB_ID = re.compile(r"[0-9]+")
_DEPENDENCY = re.compile(r"(?:after|afterany|afterok|afternotok)(?::[0-9]+(?:_[0-9]+)?)+")
_JOB_NAME = re.compile(r"bbsim_[A-Za-z0-9][A-Za-z0-9._-]*")
_CONTROL = re.compile(r"[\x00-\x1f\x7f]")


class SubmitError(RuntimeError):
    """Refused before anything was reserved or submitted (exit 2); the message names what and why."""


class SubmissionNotLogged(RuntimeError):
    """The job was submitted and its reservation bound, but its ``submissions.log`` line could not be written
    (exit 1)."""


def _one_line(text: object) -> str:
    return " ".join(str(text).split())


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _git(repo: Path, *args: str) -> subprocess.CompletedProcess:
    """git in ``repo`` with no ``GIT_*`` variable of the caller's (``GIT_DIR`` would name another repository)."""
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    command = ["git", "--no-optional-locks", "-C", str(repo), *args]
    try:
        return subprocess.run(command, capture_output=True, timeout=GIT_TIMEOUT_S, env=env)
    except subprocess.TimeoutExpired:
        raise SubmitError(f"{' '.join(command)} timed out after {GIT_TIMEOUT_S} s") from None
    except OSError as exc:
        raise SubmitError(f"git could not run on {repo}: {exc}") from None


def _git_error(run: subprocess.CompletedProcess) -> str:
    return _one_line(run.stderr.decode(errors="replace"))


# --- the script: the recorded bridge blob ---------------------------------------------------------------------------

def capture_script(bridge_repo: Path, commit: str, relpath: str, sha256: str) -> bytes:
    """The bytes of ``relpath`` at ``commit`` in the bridge clone, read as a git blob, never from the working tree, and
    required to have the sha256 the plan recorded (OR-044). A commit or blob the clone no longer holds (a
    force-push, a gc, a re-clone) or other bytes raise ``SubmitError`` naming the commit: the plan's scripts are
    submitted from that commit only."""
    if not isinstance(commit, str) or not _COMMIT.fullmatch(commit):
        raise SubmitError(f"the recorded bridge commit {commit!r} is not a full commit id")
    if not isinstance(relpath, str) or not relpath or relpath.startswith(("/", "-")) or "\\" in relpath \
            or _CONTROL.search(relpath) or any(part in ("", ".", "..") for part in relpath.split("/")):
        raise SubmitError(f"the bridge script {relpath!r} is not a path inside the bridge clone")
    if not isinstance(sha256, str) or not _HEX64.fullmatch(sha256):
        raise SubmitError(f"the recorded sha256 of {relpath} is {sha256!r}, not a sha256")
    repo = Path(bridge_repo)
    if _git(repo, "cat-file", "-e", f"{commit}^{{commit}}").returncode != 0:
        raise SubmitError(f"the recorded bridge commit {commit} is not in the bridge clone {repo}: it is no longer "
                          "reachable there (a force-push, a gc or a re-clone). A plan's scripts are submitted from its "
                          "commit only, never from the live files; nothing is submitted. Restore the commit in the "
                          "clone, or plan again")
    run = _git(repo, "cat-file", "blob", f"{commit}:{relpath}")
    if run.returncode != 0:
        raise SubmitError(f"the recorded bridge commit {commit} holds no {relpath} in {repo} ({_git_error(run)}); "
                          "nothing is submitted")
    digest = _sha256(run.stdout)
    if digest != sha256:
        raise SubmitError(f"{relpath} at the recorded bridge commit {commit} has sha256 {digest}, not the plan's "
                          f"{sha256}; nothing is submitted")
    return run.stdout


# --- the job's context ----------------------------------------------------------------------------------------------

def _check_context(context: Any) -> dict:
    if not isinstance(context, Mapping):
        raise SubmitError(f"a job context is an object, got {type(context).__name__}")
    keys = set(context)
    if keys != set(JOB_KEYS):
        raise SubmitError(f"a job context holds exactly {', '.join(JOB_KEYS)}; missing "
                          f"{sorted(set(JOB_KEYS) - keys)}, unexpected {sorted(map(str, keys - set(JOB_KEYS)))}")
    for name in ("plan", "guard_dir"):
        value = context[name]
        if not isinstance(value, str) or not os.path.isabs(value):
            raise SubmitError(f"the job context's {name} must be an absolute path, never null, got {value!r}")
    for name in ("plan_sha256", "script_sha256"):
        if not isinstance(context[name], str) or not _HEX64.fullmatch(context[name]):
            raise SubmitError(f"the job context's {name} must be a sha256, never null, got {context[name]!r}")
    for name in ("bridge_commit", "code_commit"):
        if not isinstance(context[name], str) or not _COMMIT.fullmatch(context[name]):
            raise SubmitError(f"the job context's {name} must be a full commit id, got {context[name]!r}")
    if context["kind"] not in WHATS:
        raise SubmitError(f"the job context's kind must be one of {', '.join(WHATS)}, got {context['kind']!r}")
    if context["class"] is not None and not isinstance(context["class"], str):
        raise SubmitError(f"the job context's class must be a name or null, got {context['class']!r}")
    if not isinstance(context["script"], str) or not context["script"]:
        raise SubmitError(f"the job context's script must name the bridge script, got {context['script']!r}")
    extra = context["extra_args"]
    if not isinstance(extra, list) or not all(isinstance(a, str) for a in extra):
        raise SubmitError(f"the job context's extra_args must be a list of strings, got {extra!r}")
    return {name: (list(context[name]) if name == "extra_args" else context[name]) for name in JOB_KEYS}


def encode_job(context: Mapping) -> str:
    """``context`` (exactly ``JOB_KEYS``) as unpadded base64url of compact UTF-8 JSON: only ``[A-Za-z0-9_-]``, so the
    value travels in the environment unchanged and no shell, ``--export`` list or comma can split it."""
    data = json.dumps(_check_context(context), separators=(",", ":"), ensure_ascii=False, sort_keys=True)
    return base64.urlsafe_b64encode(data.encode("utf-8")).rstrip(b"=").decode("ascii")


def decode_job(text: str) -> dict:
    """The context ``encode_job`` made: only the unpadded base64url alphabet, the canonical encoding of UTF-8 JSON with
    exactly ``JOB_KEYS``; anything else raises ``SubmitError``."""
    if not isinstance(text, str) or not _JOB_TEXT.fullmatch(text):
        shown = text[:40] + "..." if isinstance(text, str) and len(text) > 40 else text
        raise SubmitError(f"{JOB_ENV} must be unpadded base64url, [A-Za-z0-9_-] only, got {shown!r}")
    if len(text) % 4 == 1:
        raise SubmitError(f"{JOB_ENV} has a length no base64 text has ({len(text)})")
    try:
        data = base64.urlsafe_b64decode(text + "=" * (-len(text) % 4))
    except (binascii.Error, ValueError) as exc:
        raise SubmitError(f"{JOB_ENV} is not base64url: {exc}") from None
    if base64.urlsafe_b64encode(data).rstrip(b"=").decode("ascii") != text:
        raise SubmitError(f"{JOB_ENV} is not the canonical encoding of its bytes")
    try:
        context = json.loads(data.decode("utf-8"))
    except (UnicodeDecodeError, ValueError) as exc:
        raise SubmitError(f"{JOB_ENV} does not hold UTF-8 JSON: {exc}") from None
    return _check_context(context)


# --- the plan --------------------------------------------------------------------------------------------------------

def read_plan(plan_path: Path | str, plan_sha256: str) -> tuple[Path, dict]:
    """The plan at ``plan_path`` (resolved), whose bytes must have ``plan_sha256``, parsed from those same bytes."""
    if not isinstance(plan_sha256, str) or not _HEX64.fullmatch(plan_sha256):
        raise SubmitError(f"the plan's sha256 must be 64 lowercase hex digits, got {plan_sha256!r}")
    try:
        path = Path(os.fspath(plan_path)).resolve(strict=True)
        data = path.read_bytes()
    except (OSError, RuntimeError, TypeError) as exc:
        raise SubmitError(f"the plan {plan_path} cannot be read ({exc})") from None
    digest = _sha256(data)
    if digest != plan_sha256:
        raise SubmitError(f"the plan {path} has sha256 {digest}, not the given {plan_sha256}: it is not the plan that "
                          "was frozen")
    try:
        plan = json.loads(data.decode("utf-8"))
    except (UnicodeDecodeError, ValueError) as exc:
        raise SubmitError(f"the plan {path} is not JSON ({exc})") from None
    if not isinstance(plan, dict) or plan.get("schema") != array_plan.PLAN_SCHEMA:
        raise SubmitError(f"the plan {path} is not a {array_plan.PLAN_SCHEMA} plan")
    run_dir = plan.get("run_dir")
    if not isinstance(run_dir, str) or not os.path.isabs(run_dir) or Path(run_dir).resolve() / "plan.json" != path:
        raise SubmitError(f"the plan {path} is not its run directory's plan.json (run_dir {run_dir!r})")
    return path, plan


def _get(plan: Mapping, *keys: str) -> Any:
    value: Any = plan
    for key in keys:
        if not isinstance(value, Mapping) or key not in value:
            raise SubmitError(f"the plan records no {'.'.join(keys)}")
        value = value[key]
    return value


def _check_root(plan: Mapping, guard: Any) -> Path:
    """The guard's root, which must be the plan's ``guard_dir`` and any inherited ``BBSIM_LICENCE_GUARD``'s."""
    root = getattr(guard, "root", None)
    if not isinstance(root, Path):
        raise SubmitError(f"submit takes a licence_guard.Guard, got {type(guard).__name__}")
    recorded = _get(plan, "guard_dir")
    if not isinstance(recorded, str) or Path(recorded).resolve() != root:
        raise SubmitError(f"the plan's guard_dir {recorded!r} is not this guard {root}; one guard holds the shared cap "
                          "(OR-045)")
    inherited = os.environ.get(licence_guard.GUARD_ENV)
    if inherited and Path(inherited).resolve() != root:
        raise SubmitError(f"the inherited {licence_guard.GUARD_ENV}={inherited} names another guard than {root}; one "
                          "guard holds the shared cap (OR-045)")
    return root


def _check_code(plan: Mapping, env: EnvironmentState | None) -> EnvironmentState:
    """The plan's worktree, clean at the plan's commit; ``env`` is its ``manifest.environment_state`` (tests inject
    it, with the PyAEDT version a Mac does not have)."""
    repo, commit = Path(_get(plan, "code", "repo")), _get(plan, "code", "commit")
    if env is None:
        try:
            env = environment_state(repo)
        except EnvironmentStateError as exc:
            raise SubmitError(f"the plan's code {repo} cannot be read: {exc}") from None
        except PackageNotFoundError as exc:
            raise SubmitError(f"the installed PyAEDT version cannot be read: {exc}") from None
    if not isinstance(env, EnvironmentState):
        raise SubmitError(f"the code's state must be an EnvironmentState, got {env!r}")
    if env.commit != commit:
        raise SubmitError(f"the plan's code {repo} is at {env.commit}, not the plan's commit {commit}: every job runs "
                          "the code the plan froze")
    if env.dirty:
        raise SubmitError(f"the plan's code {repo} has uncommitted changes to tracked files; restore it to the plan's "
                          f"commit {commit}")
    return env


def _check_bridge(plan: Mapping) -> Path:
    """The plan's bridge, which must be its code checkout's sibling clone ``../bbridge`` (OR-044), as the planner
    admitted it; checked before any script is captured."""
    source, recorded = _get(plan, "code", "source_repo"), _get(plan, "bridge", "repo")
    try:
        expected = array_plan._sibling_bridge(Path(source)).resolve(strict=True)
        actual = Path(recorded).resolve(strict=True)
    except (OSError, RuntimeError, TypeError) as exc:
        raise SubmitError(f"the plan's bridge {recorded!r} or its code checkout's sibling clone cannot be reached "
                          f"({exc}) (OR-044)") from None
    if actual != expected:
        raise SubmitError(f"the plan's bridge {actual} is not its code checkout's sibling clone {expected} (OR-044): "
                          "a plan's scripts come from that clone only")
    return actual


def _plan_check(check: Any, *args: Any) -> Any:
    try:
        return check(*args)
    except array_plan.PlanError as exc:
        raise SubmitError(str(exc)) from None


def _array_class(plan: Mapping, class_name: str | None, indices: str | None) -> tuple[dict, list[int], bool]:
    """The class an array submission is of, its task indices and whether it is a rerun (``indices`` given)."""
    classes = _get(plan, "classes")
    tasks = _get(plan, "tasks")
    if not isinstance(classes, list) or not classes or not isinstance(tasks, list) or not all(
            isinstance(c, Mapping) and isinstance(c.get("name"), str) and isinstance(c.get("indices"), list)
            for c in classes):
        raise SubmitError("the plan's classes or tasks are not recorded as lists of records")
    by_name = {c["name"]: dict(c) for c in classes}
    names = ", ".join(by_name)
    if class_name is not None and class_name not in by_name:
        raise SubmitError(f"the plan has no class {class_name!r}; its classes are {names}")
    if indices is None:
        if class_name is None:
            if len(classes) != 1:
                raise SubmitError(f"the plan has the classes {names}: name the class of the array (--class)")
            class_name = classes[0]["name"]
        return by_name[class_name], sorted(by_name[class_name]["indices"]), False
    if not isinstance(indices, str):
        raise SubmitError(f"--indices takes the task indices as text such as 2,5 or 0-3, got {indices!r}")
    chosen = _plan_check(array_plan._parse_indices, indices, len(tasks), "--indices")
    if class_name is None:
        owners = sorted({c["name"] for c in classes for i in chosen if i in c["indices"]})
        if len(owners) != 1:
            raise SubmitError(f"--indices {indices} lie in the classes {owners}: a rerun is one class's array, so give "
                              "one rerun line per class")
        class_name = owners[0]
    outside = [i for i in chosen if i not in by_name[class_name]["indices"]]
    if outside:
        raise SubmitError(f"--indices {outside} are not in class {class_name} ({by_name[class_name]['indices']})")
    return by_name[class_name], chosen, True


def _format_indices(indices: Sequence[int]) -> str:
    """Slurm's array index list: 0-3,5,7-8."""
    parts: list[str] = []
    values = sorted(set(indices))
    start = prev = values[0]
    for value in values[1:] + [None]:
        if value is not None and value == prev + 1:
            prev = value
            continue
        parts.append(str(start) if start == prev else f"{start}-{prev}")
        if value is not None:
            start = prev = value
    return ",".join(parts)


def _check_reference(plan: Mapping, env: EnvironmentState, policy: str) -> None:
    """A production array's HF-041 column-order reference, re-hashed and re-checked against the worktree's exporter
    files and the installed PyAEDT; a probe plan records none and runs no schema-2 check."""
    reference = _get(plan, "column_order_reference")
    if policy == "probe":
        if reference is not None:
            raise SubmitError("a probe plan takes no column-order reference (its schema-1 diagnostic writes no schema "
                              "2), but this plan records one")
        return
    if not isinstance(reference, Mapping):
        raise SubmitError(f"a production plan records the HF-041 column-order reference, but this one records "
                          f"{reference!r}")
    copy, recorded = reference.get("copy"), reference.get("sha256")
    try:
        data = Path(copy).read_bytes()
    except (OSError, TypeError) as exc:
        raise SubmitError(f"the plan's column-order reference copy {copy!r} cannot be read ({exc})") from None
    if _sha256(data) != recorded:
        raise SubmitError(f"the plan's column-order reference copy {copy} has sha256 {_sha256(data)}, not the plan's "
                          f"{recorded} (HF-041)")
    try:
        parsed = json.loads(data.decode("utf-8"))
    except (UnicodeDecodeError, ValueError) as exc:
        raise SubmitError(f"the plan's column-order reference copy {copy} is not JSON ({exc})") from None
    from bbsim import dataset_v2

    try:
        dataset_v2.check_reference_before_aedt(parsed, env.pyaedt_version, Path(_get(plan, "code", "repo")))
    except (ValueError, OSError) as exc:
        raise SubmitError(f"the plan's column-order reference no longer fits its code and PyAEDT: {_one_line(exc)} "
                          "(HF-041)") from None


def _read_log(run_dir: Path) -> list[dict]:
    path = run_dir / SUBMISSIONS_LOG
    if not path.exists():
        return []
    lines = []
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        try:
            record = json.loads(line)
        except ValueError:
            record = None
        if not isinstance(record, dict) or not isinstance(record.get("kind"), str) \
                or not isinstance(record.get("id"), str):
            raise SubmitError(f"{path} line {number} is not a JSON object with a kind and an id; repair the log "
                              "before submitting more")
        lines.append(record)
    return lines


def _append_log(run_dir: Path, record: dict) -> None:
    with open(run_dir / SUBMISSIONS_LOG, "a", encoding="utf-8") as fh:
        fh.write(json.dumps(record, sort_keys=True) + "\n")
        fh.flush()
        os.fsync(fh.fileno())


def _plan_arrays(lines: Sequence[Mapping], plan_sha256: str) -> list[dict]:
    return [dict(line) for line in lines if line.get("kind") == "array" and line.get("plan_sha256") == plan_sha256]


def _resource(check: Any, value: Any, what: str) -> Any:
    try:
        check(value, what)
    except array_plan.PlanError as exc:
        raise SubmitError(str(exc)) from None
    return value


def _check_sbatch_args(plan: Mapping) -> list[str]:
    args = _get(plan, "slurm", "sbatch_args")
    if not isinstance(args, list) or not all(isinstance(a, str) for a in args):
        raise SubmitError(f"the plan's slurm.sbatch_args must be a list of strings, got {args!r}")
    try:
        licence_guard.check_user_sbatch_args(args)
    except licence_guard.GuardRefused as exc:
        raise SubmitError(f"the plan's sbatch_args: {exc}") from None
    named = {token.split("=", 1)[0] for token in args if token.startswith("--")}
    twice = sorted(named & set(PLAN_RESOURCE_OPTIONS))
    if twice:
        raise SubmitError(f"the plan's sbatch_args set {', '.join(twice)}, which the plan sets itself from --cores and "
                          "its classes (a rerun takes --mem, --time and --qos); sbatch would take one silently")
    return list(args)


# --- submission -------------------------------------------------------------------------------------------------------

def submit(plan_path: Path, plan_sha256: str, what: str, *, class_name: str | None = None, indices: str | None = None,
           mem: str | None = None, time: str | None = None, qos: str | None = None, guard: Any,
           dependency: str | None = None, env: EnvironmentState | None = None) -> str:
    """Submit one job of the plan through ``guard`` (the module docstring's checks, then ``Guard.submit``) and log it.

    ``what`` is ``array`` (``class_name`` and the class's indices; with ``indices`` a rerun of those, its class found
    from them), ``monitor`` or ``assembly``; ``mem``, ``time`` and ``qos`` replace the class's (or the monitor's or the
    assembly's) for this job. ``env`` is the worktree's ``manifest.environment_state`` (read when None). Returns the
    job id. ``SubmitError`` and ``licence_guard.GuardRefused`` refuse with nothing reserved; ``GuardError`` is a
    failure with the reservation kept as the guard says; ``SubmissionNotLogged`` names a submitted job whose log line
    is missing."""
    if what not in WHATS:
        raise SubmitError(f"what must be one of {', '.join(WHATS)}, got {what!r}")
    path, plan = read_plan(plan_path, plan_sha256)
    _check_root(plan, guard)
    state = _check_code(plan, env)
    bridge = _check_bridge(plan)
    run_dir = Path(plan["run_dir"])
    code = _get(plan, "code", "repo")
    policy = _get(plan, "policy")
    if policy not in ("probe", "production") or _get(plan, "probe") is not (policy == "probe"):
        raise SubmitError(f"the plan's policy {policy!r} and probe {plan.get('probe')!r} disagree")
    job_name, comment = _get(plan, "slurm", "job_name"), _get(plan, "slurm", "comment")
    if not isinstance(job_name, str) or not _JOB_NAME.fullmatch(job_name):
        raise SubmitError(f"the plan's slurm.job_name {job_name!r} is not bbsim_<name>")
    if not isinstance(comment, str) or not comment.startswith(licence_guard.MANAGED_COMMENT_PREFIX) \
            or _CONTROL.search(comment):
        raise SubmitError(f"the plan's slurm.comment {comment!r} does not start with "
                          f"{licence_guard.MANAGED_COMMENT_PREFIX}, which marks the plan's jobs for the guard")
    user_args = _check_sbatch_args(plan)
    if dependency is not None and (not isinstance(dependency, str) or not _DEPENDENCY.fullmatch(dependency)):
        raise SubmitError(f"--dependency {dependency!r} is not after[any|ok|notok]:JOBID[:JOBID...]")
    lines = _read_log(run_dir)
    resources: list[str]
    cls: dict | None = None
    chosen: list[int] | None = None
    rerun = False
    extra: list[str] = []
    slots = 0
    if what == "array":
        cls, chosen, rerun = _array_class(plan, class_name, indices)
        if not rerun:
            earlier = [line["id"] for line in _plan_arrays(lines, plan_sha256)
                       if line.get("class") == cls["name"] and line.get("rerun") is not True]
            if earlier:
                raise SubmitError(f"class {cls['name']} was submitted already as job {', '.join(map(str, earlier))} "
                                  f"({run_dir / SUBMISSIONS_LOG}); submit its failed indices with rerun --indices, and "
                                  f"print what remains of the plan's block with python -m bbsim.array_submit "
                                  f"remaining {shlex.quote(str(path))} {plan_sha256}")
        throttle = cls.get("throttle")
        if isinstance(throttle, bool) or not isinstance(throttle, int) or throttle < 1:
            raise SubmitError(f"class {cls['name']}'s throttle {throttle!r} is not a positive integer")
        slots = min(throttle, len(chosen))
        cores = _get(plan, "slurm", "cores")
        if isinstance(cores, bool) or not isinstance(cores, int) or cores < 1:
            raise SubmitError(f"the plan's slurm.cores {cores!r} is not a positive integer")
        resources = [f"--cpus-per-task={cores}",
                     f"--mem={_resource(array_plan._mem, mem if mem is not None else cls.get('mem'), 'mem')}",
                     f"--time={_resource(array_plan._seconds, time if time is not None else cls.get('time'), 'time')}",
                     f"--qos={_resource(array_plan._qos, qos if qos is not None else cls.get('qos'), 'qos')}"]
        installed, planned = state.pyaedt_version, _get(plan, "code", "pyaedt_version")
        if installed != planned:
            raise SubmitError(f"the installed PyAEDT {installed!r} is not the plan's {planned!r}; every task refuses "
                              "it before AEDT")
        _plan_check(array_plan._check_cores, guard, cores)
        _check_reference(plan, state, policy)
    else:
        if class_name is not None or indices is not None:
            raise SubmitError(f"a {what} job takes no class or indices")
        if what == "monitor":
            if _get(plan, "slurm", "monitor", "enabled") is not True:
                raise SubmitError("the plan runs no licence monitor (it was planned without --monitor-licences)")
            extra = [line["id"] for line in _plan_arrays(lines, plan_sha256)]
            if not extra:
                raise SubmitError(f"the licence monitor watches the plan's arrays, and no array of it is in "
                                  f"{run_dir / SUBMISSIONS_LOG} yet: submit the arrays first")
            planned_time, planned_mem, planned_qos = _get(plan, "slurm", "monitor", "time"), MONITOR_MEM, None
        else:
            assembly = _get(plan, "slurm", "assembly")
            planned_time, planned_mem, planned_qos = (_get(plan, "slurm", "assembly", "time"),
                                                      _get(plan, "slurm", "assembly", "mem"), assembly.get("qos"))
        resources = ["--cpus-per-task=1",
                     f"--mem={_resource(array_plan._mem, mem if mem is not None else planned_mem, 'mem')}",
                     f"--time={_resource(array_plan._seconds, time if time is not None else planned_time, 'time')}"]
        chosen_qos = qos if qos is not None else planned_qos
        if chosen_qos is not None:
            resources.append(f"--qos={_resource(array_plan._qos, chosen_qos, 'qos')}")

    relpath = (licence_guard.BRIDGE_SCRIPTS_DIR / SCRIPTS[what]).as_posix()
    scripts = _get(plan, "bridge", "scripts")
    commit = _get(plan, "bridge", "commit")
    if not isinstance(scripts, Mapping) or relpath not in scripts:
        raise SubmitError(f"the plan records no bridge script {relpath}")
    data = capture_script(bridge, commit, relpath, scripts[relpath])
    context = {"plan": str(path), "plan_sha256": plan_sha256, "kind": what, "class": cls["name"] if cls else None,
               "script": relpath, "script_sha256": scripts[relpath], "bridge_commit": commit, "extra_args": extra,
               "code_commit": _get(plan, "code", "commit"), "guard_dir": _get(plan, "guard_dir")}
    encoded = encode_job(context)

    argv = [f"--job-name={job_name}{JOB_NAME_SUFFIX[what]}", f"--comment={comment}"]
    if what == "array":
        argv.append(f"--array={_format_indices(chosen)}%{slots}")
    argv += [f"--chdir={code}",
             f"--output={run_dir / 'slurm'}/{'%x_%A_%a.out' if what == 'array' else '%x_%j.out'}", *resources]
    if dependency is not None:
        argv += [f"--dependency={dependency}", "--kill-on-invalid-dep=yes"]
    if what == "array":
        argv += user_args
    job = guard.submit(kind=what, slots=slots, cap=_get(plan, "licence_cap"), argv=argv, stdin=data,
                       env={CODE_DIR_ENV: code, JOB_ENV: encoded},
                       context={**context, "source": "bbsim.array_submit", "indices": chosen, "rerun": rerun,
                                "dependency": dependency})
    record = {"kind": what, "class": context["class"], "id": job, "utc": _now(), "script": relpath,
              "script_sha256": scripts[relpath], "bridge_commit": commit, "reservation": guard.last_reservation_id,
              "admitted": True, "plan_sha256": plan_sha256, "indices": chosen,
              "rerun": rerun if what == "array" else None, "dependency": dependency}
    try:
        _append_log(run_dir, record)
    except OSError as exc:
        raise SubmissionNotLogged(f"job {job} was submitted (reservation {guard.last_reservation_id}), but its line "
                                  f"could not be written to {run_dir / SUBMISSIONS_LOG} ({exc}); record it with "
                                  f"plan_hfss_array.py --log-submission {run_dir} "
                                  f"{what + ':' + context['class'] if cls else what} {job}") from exc
    return job


def _after_stop(exc: BaseException, plan_path: Path, plan_sha256: str) -> str:
    """What a stopped ``--submit`` prints next: after a refusal (nothing reserved or submitted) what remains of the
    block, ready to paste; after any other failure the command that prints it, to run once that failure is dealt
    with."""
    if isinstance(exc, (SubmitError, licence_guard.GuardRefused)):
        try:
            return "Once the guard admits it, paste what remains:\n" + remaining_block(plan_path, plan_sha256)
        except (SubmitError, OSError) as err:
            return f"what remains cannot be printed ({_one_line(err)})"
    return (f"its result is uncertain or unrecorded: do what the message above says first, then print what remains "
            f"with python -m bbsim.array_submit remaining {shlex.quote(str(plan_path))} {plan_sha256}")


def submit_plan(plan_path: Path, plan_sha256: str, *, guard: Any, env: EnvironmentState | None = None) -> int:
    """The planner's ``--submit``: every class's array, the monitor when the plan runs one, then the assembly after
    every array (``afterany``), each through ``submit``. A refused or failed array stops there (exit 1): the arrays
    before it keep running, and what remains is printed as the block's STOP prints it. A monitor or an assembly that
    is not submitted is reported the same way (exit 1), and the arrays keep running."""
    failures = (SubmitError, SubmissionNotLogged, licence_guard.GuardError, OSError)
    try:
        _, plan = read_plan(plan_path, plan_sha256)
        names = [c["name"] for c in _get(plan, "classes")]
        monitor = _get(plan, "slurm", "monitor", "enabled") is True
    except (SubmitError, KeyError, TypeError) as exc:
        print(f"plan_hfss_array: --submit: refused: {_one_line(exc)}", file=sys.stderr)
        return 1
    arrays = []
    for name in names:
        try:
            job = submit(plan_path, plan_sha256, "array", class_name=name, guard=guard, env=env)
        except failures as exc:
            print(f"plan_hfss_array: --submit: class {name} was not submitted: {_one_line(exc)}. The arrays "
                  f"submitted before it ({', '.join(arrays) or 'none'}) keep running; no monitor or assembly was "
                  "submitted", file=sys.stderr)
            print(_after_stop(exc, plan_path, plan_sha256))
            return 1
        arrays.append(job)
        print(f"array {name}: job {job}")
    stopped: BaseException | None = None
    jobs = [("monitor", {})] if monitor else []
    jobs.append(("assembly", {"dependency": "afterany:" + ":".join(arrays)}))
    for what, kw in jobs:
        try:
            print(f"{what}: job {submit(plan_path, plan_sha256, what, guard=guard, env=env, **kw)}")
        except failures as exc:
            print(f"plan_hfss_array: --submit: the {what} was not submitted: {_one_line(exc)}; the arrays keep "
                  "running", file=sys.stderr)
            stopped = stopped or exc
    if stopped is None:
        return 0
    print(_after_stop(stopped, plan_path, plan_sha256))
    return 1


def change_throttle(plan_path: Path, plan_sha256: str, array_id: str, throttle: int, *, guard: Any) -> None:
    """A live throttle change of one of this plan's arrays, through ``Guard.change_throttle`` with the plan's cap
    (admission repeated for an increase), logged in ``submissions.log``."""
    path, plan = read_plan(plan_path, plan_sha256)
    _check_root(plan, guard)
    run_dir = Path(plan["run_dir"])
    if not isinstance(array_id, str) or not _JOB_ID.fullmatch(array_id):
        raise SubmitError(f"the array's job id must be a plain number, got {array_id!r}")
    found = [line for line in _plan_arrays(_read_log(run_dir), plan_sha256) if line["id"] == array_id]
    if not found:
        raise SubmitError(f"job {array_id} is not an array of the plan {path} in {run_dir / SUBMISSIONS_LOG}")
    guard.change_throttle(array_id, throttle, cap=_get(plan, "licence_cap"))
    try:
        _append_log(run_dir, {"kind": "throttle", "class": found[-1].get("class"), "id": array_id,
                              "throttle": throttle, "utc": _now(), "plan_sha256": plan_sha256, "admitted": True})
    except OSError as exc:
        raise SubmissionNotLogged(f"array {array_id}'s throttle was changed to {throttle} through the guard (its "
                                  f"events.log in {guard.root} records it), but {run_dir / SUBMISSIONS_LOG} could not "
                                  f"record the change ({exc})") from exc


def remaining_block(plan_path: Path, plan_sha256: str) -> str:
    """What remains of the plan's printed block (``array_plan.submission_block``) after the jobs in
    ``submissions.log``: the classes not yet submitted, the monitor when the plan runs one and it is not in the log,
    and the assembly unless one follows the last class's array in the log; the classes already submitted appear as
    their literal ids, so the lines hold no variable they do not set. Nothing is submitted here."""
    path, plan = read_plan(plan_path, plan_sha256)
    lines = _read_log(Path(plan["run_dir"]))
    names = [c.get("name") for c in _get(plan, "classes") if isinstance(c, Mapping)]
    admitted: dict[str, str] = {}
    last_array = -1
    for number, line in enumerate(lines):
        if line["kind"] == "array" and line.get("plan_sha256") == plan_sha256 and line.get("rerun") is not True \
                and line.get("class") in names:
            if not _JOB_ID.fullmatch(line["id"]):
                raise SubmitError(f"{SUBMISSIONS_LOG} names class {line['class']}'s job {line['id']!r}, not a job id")
            admitted[line["class"]] = line["id"]
            last_array = number
    ours = [(number, line["kind"]) for number, line in enumerate(lines) if line.get("plan_sha256") == plan_sha256]
    monitor_submitted = any(kind == "monitor" for _, kind in ours)
    assembly_submitted = len(admitted) == len(names) and any(kind == "assembly" and number > last_array
                                                             for number, kind in ours)
    monitor_needed = _get(plan, "slurm", "monitor", "enabled") is True and not monitor_submitted
    if len(admitted) == len(names) and assembly_submitted and not monitor_needed:
        return (f"nothing remains: every class, the assembly and any monitor of the plan {path} are in "
                f"{Path(plan['run_dir']) / SUBMISSIONS_LOG}")
    try:
        return array_plan.submission_block(path, plan_sha256, plan, admitted=admitted,
                                           monitor_submitted=monitor_submitted, assembly_submitted=assembly_submitted)
    except (KeyError, TypeError) as exc:
        raise SubmitError(f"the plan {path} cannot give its block: {exc!r} is missing or malformed") from None


def log_submission(run_dir: Path | str, kind: str, job_id: str) -> dict:
    """Record in ``<RUN>/submissions.log`` a job submitted by hand: ``kind`` is ``array:<class>``, ``monitor`` or
    ``assembly``, ``job_id`` the plain number sbatch printed. The line says ``admitted: false``: nothing is reserved
    and no admission is granted."""
    try:
        run = Path(os.fspath(run_dir)).resolve(strict=True)
        data = (run / "plan.json").read_bytes()
    except (OSError, RuntimeError, TypeError) as exc:
        raise SubmitError(f"{run_dir} holds no readable plan.json ({exc})") from None
    _, plan = read_plan(run / "plan.json", _sha256(data))
    what, sep, name = kind.partition(":") if isinstance(kind, str) else ("", "", "")
    classes = [c.get("name") for c in _get(plan, "classes") if isinstance(c, Mapping)]
    if what == "array" and sep and name in classes:
        cls = name
    elif what in ("monitor", "assembly") and not sep:
        cls = None
    else:
        raise SubmitError(f"KIND is array:<class> (the plan's classes: {', '.join(map(str, classes))}), monitor or "
                          f"assembly, got {kind!r}")
    if not isinstance(job_id, str) or not _JOB_ID.fullmatch(job_id):
        raise SubmitError(f"ID is the plain job number sbatch --parsable printed (without ;cluster), got {job_id!r}")
    record = {"kind": what, "class": cls, "id": job_id, "utc": _now(), "script": None, "script_sha256": None,
              "bridge_commit": None, "reservation": None, "admitted": False, "plan_sha256": _sha256(data),
              "indices": None, "rerun": None, "dependency": None}
    _append_log(run, record)
    return record


def log_submission_main(argv: Sequence[str]) -> int:
    """``plan_hfss_array.py --log-submission RUN KIND ID``: 0 recorded, 2 refused, 1 the log could not be written."""
    parser = argparse.ArgumentParser(prog="plan_hfss_array.py --log-submission", allow_abbrev=False,
                                     description="Record a job submitted by hand (no admission, admitted: false).")
    parser.add_argument("run", metavar="RUN")
    parser.add_argument("kind", metavar="KIND", help="array:<class>, monitor or assembly")
    parser.add_argument("job_id", metavar="ID")
    try:
        args = parser.parse_args(list(argv))
    except SystemExit as exc:
        return 2 if exc.code else 0
    try:
        record = log_submission(args.run, args.kind, args.job_id)
    except SubmitError as exc:
        print(f"plan_hfss_array: refused: {_one_line(exc)}", file=sys.stderr)
        return 2
    except OSError as exc:
        print(f"plan_hfss_array: {_one_line(exc)}", file=sys.stderr)
        return 1
    print(f"recorded {args.kind} job {record['id']} in {Path(args.run).resolve() / SUBMISSIONS_LOG} with admitted: "
          "false: the guard reserved nothing for it (OR-045)")
    return 0


# --- the job's side: job-field ---------------------------------------------------------------------------------------

def _checkout_state(repo: Path) -> tuple[str, bool]:
    """``HEAD`` and whether tracked files differ from it, by the git commands of ``manifest.environment_state``."""
    head = _git(repo, "rev-parse", "--verify", "HEAD")
    status = _git(repo, "status", "--porcelain", "--untracked-files=no")
    if head.returncode != 0 or status.returncode != 0:
        raise SubmitError(f"git cannot read the code {repo}: {_git_error(head if head.returncode else status)}")
    return head.stdout.decode().strip(), bool(status.stdout.strip())


def job_field(name: str, *, environ: Mapping[str, str] | None = None, root: Path | None = None,
              module_root: Path | None = None) -> str:
    """One value of the job's context, after checking the context against its plan, the code and the guard.

    ``environ`` is the job's environment (``os.environ``), ``root`` the shared guard (``licence_guard.DEFAULT_ROOT``),
    ``module_root`` the checkout this module was imported from (tests name the plan's worktree)."""
    env = os.environ if environ is None else environ
    if name not in JOB_KEYS:
        raise SubmitError(f"job-field takes one of {', '.join(JOB_KEYS)}, not {name!r}")
    text = env.get(JOB_ENV)
    if not text:
        raise SubmitError(f"{JOB_ENV} is not set: this job was not submitted by python -m bbsim.array_submit")
    context = decode_job(text)
    path, plan = read_plan(context["plan"], context["plan_sha256"])

    given = env.get(CODE_DIR_ENV)
    if not given:
        raise SubmitError(f"{CODE_DIR_ENV} is not set: this job was not submitted by python -m bbsim.array_submit")
    try:
        code = Path(given).resolve(strict=True)
        planned = Path(_get(plan, "code", "repo")).resolve(strict=True)
    except (OSError, RuntimeError, TypeError) as exc:
        raise SubmitError(f"{CODE_DIR_ENV}={given} or the plan's code cannot be reached ({exc})") from None
    if code != planned:
        raise SubmitError(f"{CODE_DIR_ENV}={given} is not the plan's code {planned}")
    here = Path(module_root if module_root is not None else REPO_ROOT).resolve()
    if here != code:
        raise SubmitError(f"this module was imported from {here}, not the plan's code {code}: run python from "
                          f"{CODE_DIR_ENV}")
    commit = _get(plan, "code", "commit")
    if context["code_commit"] != commit:
        raise SubmitError(f"the job context's code commit {context['code_commit']} is not the plan's {commit}")
    head, dirty = _checkout_state(code)
    if head != commit:
        raise SubmitError(f"the plan's code {code} is at {head}, not the plan's commit {commit}")
    if dirty:
        raise SubmitError(f"the plan's code {code} has uncommitted changes to tracked files")

    value = env.get(licence_guard.GUARD_ENV)
    if not value:
        raise SubmitError(f"{licence_guard.GUARD_ENV} is not set: the guard sets it on every job it submits")
    shared = Path(root if root is not None else licence_guard.DEFAULT_ROOT)
    try:
        here_root, shared_root = Path(value).resolve(strict=True), shared.resolve(strict=True)
    except (OSError, RuntimeError) as exc:
        raise SubmitError(f"{licence_guard.GUARD_ENV}={value} or the shared guard {shared} cannot be reached "
                          f"({exc})") from None
    claims = {licence_guard.GUARD_ENV: here_root, "the job context's guard_dir": Path(context["guard_dir"]).resolve(),
              "the plan's guard_dir": Path(str(_get(plan, "guard_dir"))).resolve()}
    for what, claimed in claims.items():
        if claimed != shared_root:
            raise SubmitError(f"{what} {claimed} is not the shared guard {shared_root}: a second guard cannot bypass "
                              "the cap (OR-045)")

    if context["bridge_commit"] != _get(plan, "bridge", "commit"):
        raise SubmitError(f"the job context's bridge commit {context['bridge_commit']} is not the plan's "
                          f"{_get(plan, 'bridge', 'commit')}")
    scripts = _get(plan, "bridge", "scripts")
    expected = (licence_guard.BRIDGE_SCRIPTS_DIR / SCRIPTS[context["kind"]]).as_posix()
    recorded = scripts.get(expected) if isinstance(scripts, Mapping) else None
    if context["script"] != expected or recorded != context["script_sha256"]:
        raise SubmitError(f"the job context's script {context['script']} ({context['script_sha256']}) is not the "
                          f"plan's {context['kind']} script {expected} ({recorded})")
    names = [c.get("name") for c in _get(plan, "classes") if isinstance(c, Mapping)]
    if context["kind"] == "array" and context["class"] not in names:
        raise SubmitError(f"the job context's class {context['class']!r} is not one of the plan's classes "
                          f"{', '.join(map(str, names))}")
    if context["kind"] != "array" and context["class"] is not None:
        raise SubmitError(f"a {context['kind']} job's context names no class, but this one names {context['class']!r}")

    raw = context[name]
    values = raw if isinstance(raw, list) else ["" if raw is None else str(raw)]
    for item in values:
        if "\n" in item or "\r" in item or "\0" in item:
            raise SubmitError(f"the job context's {name} holds a newline, carriage return or NUL, which a shell cannot "
                              "read back as one value")
    return "\n".join(values)


# --- the CLI ----------------------------------------------------------------------------------------------------------

def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="python -m bbsim.array_submit", allow_abbrev=False,
                                description="Submit a job-array plan's jobs through the licence guard (OR-045).")
    p.add_argument("--root", type=Path, default=None,
                   help=f"the shared guard directory (default {licence_guard.DEFAULT_ROOT})")
    sub = p.add_subparsers(dest="command", required=True)
    s = sub.add_parser("submit", allow_abbrev=False, help="submit one array (a class), the monitor or the assembly")
    s.add_argument("plan", type=Path)
    s.add_argument("sha256")
    s.add_argument("what", choices=WHATS)
    s.add_argument("--class", dest="class_name", metavar="NAME", help="the array's class (A1, A2, ...)")
    s.add_argument("--dependency", metavar="DEP", help="after[any|ok|notok]:JOBID[:JOBID...], with "
                                                       "--kill-on-invalid-dep=yes")
    r = sub.add_parser("rerun", allow_abbrev=False, help="submit indices of one class again, as a new array")
    r.add_argument("plan", type=Path)
    r.add_argument("sha256")
    r.add_argument("--indices", required=True, metavar="I,J", help="task indices of one class: 2,5 or 0-3")
    r.add_argument("--mem", help="replaces the class's memory (OUT_OF_MEMORY)")
    r.add_argument("--time", help="replaces the class's time limit (TIMEOUT)")
    r.add_argument("--qos", help="replaces the class's QOS")
    t = sub.add_parser("throttle", allow_abbrev=False, help="change the throttle of one of the plan's arrays")
    t.add_argument("plan", type=Path)
    t.add_argument("sha256")
    t.add_argument("array_id", metavar="ARRAY_ID")
    t.add_argument("throttle", type=int)
    rest = sub.add_parser("remaining", allow_abbrev=False, help="print what remains of the plan's block, after the "
                                                                "jobs in submissions.log")
    rest.add_argument("plan", type=Path)
    rest.add_argument("sha256")
    j = sub.add_parser("job-field", allow_abbrev=False, help="one value of this job's context, verified")
    j.add_argument("name", choices=JOB_KEYS)
    return p


def main(argv: list[str] | None = None, *, scheduler: Any = None, env: EnvironmentState | None = None,
         environ: Mapping[str, str] | None = None, module_root: Path | None = None) -> int:
    """The CLI. ``scheduler``, ``env`` (the worktree's state), ``environ`` and ``module_root`` replace Slurm, the
    worktree's ``environment_state``, a job's environment and this module's checkout (tests)."""
    try:
        args = _parser().parse_args(argv)
    except SystemExit as exc:                    # argparse has printed its message
        return 2 if exc.code else 0
    try:
        if args.command == "job-field":
            print(job_field(args.name, environ=environ, root=args.root, module_root=module_root))
            return 0
        if args.command == "remaining":
            print(remaining_block(args.plan, args.sha256))
            return 0
        root = args.root if args.root is not None else licence_guard.DEFAULT_ROOT
        guard = licence_guard.Guard(root, scheduler if scheduler is not None else licence_guard.SlurmScheduler())
        if args.command == "throttle":
            change_throttle(args.plan, args.sha256, args.array_id, args.throttle, guard=guard)
            print(f"THROTTLE {args.array_id} {args.throttle}")
            return 0
        if args.command == "rerun":
            what = "array"
            job = submit(args.plan, args.sha256, what, indices=args.indices, mem=args.mem, time=args.time,
                         qos=args.qos, guard=guard, env=env)
        else:
            what = args.what
            job = submit(args.plan, args.sha256, what, class_name=args.class_name, dependency=args.dependency,
                         guard=guard, env=env)
        print(job)
        print(f"array_submit: {args.command} {what}: job {job}, reservation {guard.last_reservation_id} in "
              f"{guard.root}", file=sys.stderr)
        return 0
    except (SubmitError, licence_guard.GuardRefused) as exc:
        print(f"array_submit: refused: {_one_line(exc)}", file=sys.stderr)
        return 2
    except (licence_guard.GuardError, SubmissionNotLogged, OSError) as exc:
        print(f"array_submit: {_one_line(exc)}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
