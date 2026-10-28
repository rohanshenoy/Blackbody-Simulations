"""One shared admission for HFSS solves on the cluster: at most 10 at a time (OR-045, HF-057, HF-058).

Every licensed submission (job arrays, single solves, re-exports) and every zero-slot job the job array submits
(monitor, assembly, the cluster probe) goes through ``Guard.submit``, under one lock in one shared directory, by
default ``DEFAULT_ROOT``. Arbitrary scripts, other users and hand-started AEDT sessions are outside it; only a
scheduler-enforced limit would cover those.

The lock is an atomic ``os.mkdir(<root>/lock)`` holding ``owner.json`` {nonce, host, boot_id, pid, process_start,
slurm_job, utc}; ``flock`` on a network filesystem is not relied on. A waiting client retries every ``LOCK_POLL_S``
seconds for ``lock_timeout_s`` and then refuses with ``lock_info``'s observation and the recovery command. A lock is
never broken because it is old: a crashed owner's lock is removed only by one operator's ``recover-lock --expected-lock
TOKEN --confirm-quiescent REASON``, after every guard client is stopped and the owner is verified dead (its pid and
start time on its own host), or, when the record is absent or cannot be verified from here, after the whole workflow
is verified quiescent. Recovery serializes on ``<root>/recover.lock``, rechecks the exact lock (``TOKEN`` is the lock
directory's inode, its ctime and the owner's nonce, so a lock written or replaced since refuses), refuses a live
owner, logs the assertion and the observation, and only then moves the lock to ``<root>/lock-evidence/`` with a
``recovery.json``. It never touches a reservation or a lease. Normal release removes only its own lock instance.

Reservations: ``<root>/reservations/<id>.json`` {id, kind, slots, cap, state, job_id, context, created_utc, events,
sbatch_argv, throttle_target}, state ``pending`` (sbatch in flight, or a direct session's slot), ``bound`` (to a Slurm
job id; an array's parent id), ``uncertain`` (sbatch's answer was ambiguous) or ``released``. Every record is written
atomically under the lock, and every change appends one JSON line to ``<root>/events.log``.

Usage is reconciled under the lock from the scheduler's queue (pending, running, completing, suspended), the
reservations and the leases:

- array tasks are grouped by parent; an array whose comment starts ``bbsim:``, or that a reservation is bound to,
  counts its ``ArrayTaskThrottle``, or its remaining tasks when the throttle is 0; a ``bbsim_hfss`` or
  ``bbsim_reexport`` job counts 1, and so does any queued job a slot-bearing reservation is bound to; a held job
  still counts;
- a ``pending`` or ``uncertain`` reservation counts its slots; a ``bound`` one counts, once per job id, the largest of
  its slots, the scheduler's count for that id and its held leases while the job is queued; once the job has left the
  queue it is released only when accounting (``sacct``) calls it terminal and no held lease names it, and counts its
  slots until then;
- a direct session's reservation (kind ``direct``, ``reserve_direct``) is released only when the leases that name it
  are all released, i.e. after its own lease verified cleanup;
- a throttle decrease keeps the reservation at the largest of the new throttle, the array's running tasks and its
  held leases, and drains down to the new throttle as those stop.

Admission: the reconciled total plus the new slots must not exceed ``min(CEILING, cap)``, ``cap`` an integer 1 to 10;
no slot is admitted while an uncertain reservation holds slots. Zero-slot kinds (``monitor``, ``assembly``) are
logged and never counted. A submission reserves before calling ``sbatch --parsable --export=ALL`` and binds the
returned id before the lock is released; a non-zero exit, an id that is not a plain number (``;cluster`` is dropped)
or an error leaves the reservation ``uncertain``, counted, until ``reconcile --resolve RESERVATION=released|job:<id>
--reason TEXT`` after checking ``squeue``. ``sbatch --wait`` is never run under the lock.

Leases (``<root>/leases/<id>.json``) belong to the launch-lease layer; this module reads only their
``reservation_id`` and ``state``. Any state but ``released`` holds, unknown cleanup included (HF-058(2)); a held lease
naming no live reservation counts one slot of its own.

The guard owns the script, the job's context and identity, its array indices and throttle, its working directory and
the submission protocol. User sbatch arguments (``check_user_sbatch_args``) are complete canonical long names from
``ALLOWED_USER_OPTIONS`` only, as ``--opt=value`` or ``--opt value``, passed on byte for byte; abbreviations (sbatch's
getopt_long would take ``--pars`` as ``--parsable``), short options, unknown options, ``--`` and ``:`` are refused, and
so is every option in ``OWNED_SBATCH_OPTIONS``. Internal callers pass identity options as ``--name=value`` and may not
touch the protocol (``--export*``, ``--parsable``, ``--quiet``, ``--wait``, ``--wrap``, ...); in stdin mode no script
pathname is allowed, and the script's own ``#SBATCH`` header may not set the protocol, the array or the working
directory either. sbatch's environment is the caller's without ``SBATCH_*`` and the scheduler-option aliases
(``SCHEDULER_ENV_ALIASES``), plus the caller's validated context and, set by the guard alone, ``GUARD_ENV`` (the
resolved root) and ``RESERVATION_ENV``. An inherited or explicit ``GUARD_ENV`` or a plan's ``guard_dir`` that resolves
to another directory refuses before anything is reserved, so a second guard cannot bypass the cap.

Path mode (the CLI's ``submit``; OR-044, HF-058(3)): one script, one of ``PATH_MODE_SCRIPTS``, whose real path lies in
the bridge clone's ``scripts/hpc/hfss/`` (the sibling ``../bbridge`` of this checkout, or ``--bridge-repo``), tracked
there and byte-equal to its ``HEAD`` blob, then literal script arguments. Its kind and slots follow the script.

The two-node lock test (``lock-test-worker``, ``lock-test-check``): each worker refuses unless this checkout is clean
and at the recorded commit, waits at a barrier until the expected number of distinct hosts has announced itself, then
visits the production lock ``--iterations`` times; inside each visit it creates the test's ``active`` sentinel
directory (an existing one is an overlap), advances the shared sequence and logs an enter and a leave record. The
check judges serialization by the sentinel and the sequence in the shared log, never by clocks.

Exit codes: 2 refused before anything was reserved or submitted, 1 failed, 0 success. Standard library only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import math
import os
import re
import secrets
import socket
import subprocess
import sys
import time
from collections.abc import Callable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Protocol

log = logging.getLogger(__name__)

CEILING = 10                                   # OR-045: at most 10 HFSS solves at a time, every array and job
GUARD_ENV = "BBSIM_LICENCE_GUARD"
RESERVATION_ENV = "BBSIM_GUARD_RESERVATION"
DEFAULT_ROOT = Path("/resnick/groups/golwala/rshenoy/bbsim/licence-guard")
HFSS_CHECKOUT = Path(__file__).resolve().parents[1]
BRIDGE_SCRIPTS_DIR = Path("scripts/hpc/hfss")

LOCK_POLL_S = 0.2
SCHEDULER_TIMEOUT_S = 30.0

KINDS = ("array", "single", "reexport", "monitor", "assembly")
ZERO_SLOT_KINDS = ("monitor", "assembly")
DIRECT_KIND = "direct"
STATES = ("pending", "bound", "uncertain", "released")
SOLVE_JOB_NAMES = frozenset({"bbsim_hfss", "bbsim_reexport"})
MANAGED_COMMENT_PREFIX = "bbsim:"
LEASE_RELEASED = "released"
PENDING_QUEUE_STATES = frozenset({"PENDING"})
TERMINAL_STATES = frozenset({"BOOT_FAIL", "CANCELLED", "COMPLETED", "DEADLINE", "FAILED", "NODE_FAIL",
                             "OUT_OF_MEMORY", "PREEMPTED", "REVOKED", "TIMEOUT"})

# Path mode: script name -> (kind, slots). The two existing wrappers and I1's and I2b's zero-slot scripts.
PATH_MODE_SCRIPTS = {
    "run_frequency.sbatch": ("single", 1),
    "reexport.sbatch": ("reexport", 1),
    "probe_cluster.sbatch": ("monitor", 0),
    "guard_lock_contention.sbatch": ("monitor", 0),
    "scope_check.sbatch": ("monitor", 0),
}

USER_VALUE_OPTIONS = ("--nodelist", "--constraint", "--exclude", "--begin", "--account", "--partition",
                      "--cpus-per-task", "--mem", "--time", "--qos", "--nodes", "--ntasks", "--ntasks-per-node")
USER_FLAG_OPTIONS = ("--exclusive",)           # also --exclusive=<value>; a separate token is a pathname to sbatch
ALLOWED_USER_OPTIONS = frozenset(USER_VALUE_OPTIONS + USER_FLAG_OPTIONS)
OWNED_SBATCH_OPTIONS = frozenset({
    "--export", "--export-file", "--parsable", "--quiet", "-Q", "--job-name", "-J", "--comment", "--array", "-a",
    "--chdir", "-D", "--wrap", "--wait", "-W", "--requeue", "--no-requeue", "--test-only", "--clusters", "-M",
    "--get-user-env",
})
# Internal callers may set identity, array and dependency options, never these (nor an abbreviation of one).
_PROTOCOL_NAMES = ("export", "export-file", "parsable", "quiet", "wait", "wrap", "test-only", "clusters",
                   "get-user-env", "help", "usage", "version")
# A script's own #SBATCH header may set none of these either.
_DIRECTIVE_NAMES = ("array", "wait", "quiet", "parsable", "export", "export-file", "wrap", "test-only", "clusters",
                    "chdir", "requeue", "get-user-env")
_DIRECTIVE_SHORT = frozenset("aWQMD")
SCHEDULER_ENV_ALIASES = ("SLURM_CLUSTERS", "SLURM_HINT", "SLURM_EXIT_ERROR", "SLURM_MEM_PER_CPU",
                         "SLURM_MEM_PER_GPU", "SLURM_MEM_PER_NODE")

LOCK_TEST_LOG = "log.jsonl"
LOCK_TEST_SEQUENCE = "sequence"
LOCK_TEST_ACTIVE = "active"
LOCK_TEST_READY = "ready"
LOCK_TEST_FAILURES = "failures"

_OWNER = "owner.json"
_LONG = re.compile(r"--([A-Za-z0-9][A-Za-z0-9-]*)(?:(=)(.*))?", re.S)
_CONTROL = re.compile(r"[\x00-\x1f\x7f]")
_ENV_NAME = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
_JOB_ID = re.compile(r"[0-9]+")
_INDEX_RANGE = re.compile(r"([0-9]+)(?:-([0-9]+)(?::([0-9]+))?)?")
_ARRAY_TASK = re.compile(r"([0-9]+)_(?:([0-9]+)|\[([^\]]*)\])")
_COMMIT = re.compile(r"[0-9a-f]{40}|[0-9a-f]{64}")


class GuardError(RuntimeError):
    """The guard cannot admit, reconcile or record something; after a reservation exists, it stays recorded."""


class GuardRefused(GuardError):
    """Refused before anything was reserved or submitted (exit 2)."""


@dataclass(frozen=True)
class QueueEntry:
    """One ``squeue`` line: a job, an array task or a range of pending array tasks (``123_[4-9%3]``)."""

    job_id: str
    array_parent: str | None
    name: str
    comment: str
    state: str


class Scheduler(Protocol):
    def queue(self) -> list[QueueEntry]: ...

    def array_throttle(self, parent_id: str) -> int | None: ...

    def terminal_state(self, job_id: str) -> str | None: ...

    def submit(self, argv: list[str], stdin: bytes | None, env: Mapping[str, str]) -> tuple[int, str]: ...

    def set_throttle(self, parent_id: str, throttle: int) -> None: ...


@dataclass(frozen=True)
class Usage:
    """The reconciled count: ``sources`` holds one record per counted array, job, reservation or lease."""

    total: int
    sources: list[dict]
    uncertain: list[str]
    unmanaged: list[str]


def _query_env() -> dict[str, str]:
    return {k: v for k, v in os.environ.items()
            if not (k.startswith(("SBATCH_", "SQUEUE_", "SACCT_", "SCONTROL_")) or k in SCHEDULER_ENV_ALIASES)}


class SlurmScheduler:
    """The scheduler through Slurm's commands, each with a timeout; any failure is a ``GuardError``.

    ``sacct`` runs with ``-P`` so a state such as ``CANCELLED by 1234`` is not truncated to ``CANCELLED+``.
    """

    def __init__(self, user: str | None = None, timeout_s: float = SCHEDULER_TIMEOUT_S):
        self.user = user or _user()
        self.timeout_s = timeout_s

    def _run(self, cmd: list[str], *, stdin: bytes | None = None, env: Mapping[str, str] | None = None
             ) -> subprocess.CompletedProcess:
        try:
            if stdin is None:
                return subprocess.run(cmd, stdin=subprocess.DEVNULL, capture_output=True, timeout=self.timeout_s,
                                      env=dict(env) if env is not None else _query_env())
            return subprocess.run(cmd, input=stdin, capture_output=True, timeout=self.timeout_s,
                                  env=dict(env) if env is not None else _query_env())
        except subprocess.TimeoutExpired as exc:
            raise GuardError(f"{' '.join(cmd[:2])} timed out after {self.timeout_s:g} s") from exc
        except OSError as exc:
            raise GuardError(f"{cmd[0]} could not run: {exc}") from exc

    def _checked(self, cmd: list[str]) -> str:
        run = self._run(cmd)
        if run.returncode != 0:
            err = run.stderr.decode(errors="replace").strip()
            raise GuardError(f"{' '.join(cmd)} exited {run.returncode}: {err}")
        return run.stdout.decode(errors="replace")

    def queue(self) -> list[QueueEntry]:
        out = self._checked(["squeue", "-h", "-u", self.user, "-t", "PD,R,CG,S", "-o", "%i|%F|%j|%k|%T"])
        entries = []
        for line in out.splitlines():
            if not line.strip():
                continue
            fields = line.strip().split("|")
            if len(fields) < 5:
                raise GuardError(f"squeue line not understood: {line!r}")
            job_id, parent, state = fields[0], fields[1], fields[-1]
            name, comment = fields[2], "|".join(fields[3:-1])     # a '|' inside a name or comment lands in comment
            entries.append(QueueEntry(job_id, parent if "_" in job_id else None, name, comment, state))
        return entries

    def array_throttle(self, parent_id: str) -> int | None:
        out = self._checked(["scontrol", "show", "job", parent_id])
        match = re.search(r"\bArrayTaskThrottle=([0-9]+)", out)
        return int(match.group(1)) if match else None

    def terminal_state(self, job_id: str) -> str | None:
        out = self._checked(["sacct", "-n", "-X", "-P", "-o", "State", "-j", job_id])
        states = [line.split()[0].rstrip("+") for line in out.splitlines() if line.strip()]
        if states and all(s in TERMINAL_STATES for s in states):
            return ",".join(sorted(set(states)))
        return None

    def submit(self, argv: list[str], stdin: bytes | None, env: Mapping[str, str]) -> tuple[int, str]:
        run = self._run(["sbatch", *argv], stdin=stdin, env=env)
        if run.returncode != 0 and run.stderr:
            log.warning("sbatch: %s", run.stderr.decode(errors="replace").strip())
        return run.returncode, run.stdout.decode(errors="replace")

    def set_throttle(self, parent_id: str, throttle: int) -> None:
        self._checked(["scontrol", "update", f"JobId={parent_id}", f"ArrayTaskThrottle={throttle}"])


class _NoScheduler:
    """For a guard opened for its lock only (the lock test)."""

    def _refuse(self, *args: Any, **kwargs: Any) -> Any:
        raise GuardError("this guard was opened for its lock only and has no scheduler")

    queue = array_throttle = terminal_state = submit = set_throttle = _refuse


# --- identity of this process and of a lock owner -------------------------------------------------------------

_BOOT_ID: list[str | None] = []
_OWN_START: dict[int, str | None] = {}


def _user() -> str:
    return os.environ.get("USER") or os.environ.get("LOGNAME") or str(os.getuid())


def _boot_id() -> str | None:
    """This host's boot identity: Linux's boot_id, else the BSD boot time."""
    if not _BOOT_ID:
        value = None
        try:
            value = Path("/proc/sys/kernel/random/boot_id").read_text().strip() or None
        except OSError:
            try:
                run = subprocess.run(["sysctl", "-n", "kern.boottime"], capture_output=True, text=True, timeout=10)
                value = (run.stdout.strip() or None) if run.returncode == 0 else None
            except (OSError, subprocess.TimeoutExpired):
                value = None
        _BOOT_ID.append(value)
    return _BOOT_ID[0]


def _process_start(pid: int) -> str | None:
    """A process's start identity on this host (with its pid, a process incarnation); None if gone or unknown."""
    try:
        text = Path(f"/proc/{pid}/stat").read_text()
    except FileNotFoundError:
        if Path("/proc/self/stat").exists():
            return None                         # a /proc host: the process is gone
    except OSError:
        return None
    else:
        fields = text[text.rindex(")") + 2:].split()
        return f"linux-starttime:{fields[19]}" if len(fields) > 19 else None
    try:
        run = subprocess.run(["ps", "-o", "lstart=", "-p", str(pid)], capture_output=True, text=True, timeout=10,
                             env={**os.environ, "LC_ALL": "C"})
    except (OSError, subprocess.TimeoutExpired):
        return None
    value = run.stdout.strip()
    return f"ps-lstart:{value}" if run.returncode == 0 and value else None


def _identity() -> dict:
    return {"host": socket.gethostname(), "pid": os.getpid(), "user": _user()}


def _owner_identity() -> dict:
    pid = os.getpid()
    if pid not in _OWN_START:
        _OWN_START.clear()
        _OWN_START[pid] = _process_start(pid)
    return {"host": socket.gethostname(), "boot_id": _boot_id(), "pid": pid, "process_start": _OWN_START[pid],
            "slurm_job": os.environ.get("SLURM_JOB_ID"), "user": _user()}


def _owner_status(owner: Any) -> tuple[str, str]:
    """``alive``, ``dead`` or ``unknown`` for a lock owner record, with the reason."""
    if not isinstance(owner, dict):
        return "unknown", ("no owner record: its owner died between mkdir and writing the record, or is writing it "
                           "now; only whole-workflow quiescence can establish that it is abandoned")
    host, pid, here = owner.get("host"), owner.get("pid"), socket.gethostname()
    if host != here:
        return "unknown", f"recorded on host {host!r}; this is {here!r}, so its process cannot be checked here"
    recorded_boot, boot = owner.get("boot_id"), _boot_id()
    if recorded_boot and boot and recorded_boot != boot:
        return "dead", "this host has rebooted since the record was written"
    if isinstance(pid, bool) or not isinstance(pid, int) or pid <= 0:
        return "unknown", f"the record holds no valid pid ({pid!r})"
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return "dead", f"no process {pid} on this host"
    except PermissionError:
        pass
    except OSError as exc:
        return "unknown", f"process {pid} cannot be checked: {exc}"
    start, recorded_start = _process_start(pid), owner.get("process_start")
    if start is None or not recorded_start:
        return "unknown", f"process {pid} exists and its start cannot be compared with the record"
    if start == recorded_start:
        return "alive", f"process {pid} is the recorded owner (start {start})"
    return "dead", f"pid {pid} now belongs to another process (start {start}, recorded {recorded_start})"


# --- small file helpers ---------------------------------------------------------------------------------------

def _utc(t: float) -> str:
    return datetime.fromtimestamp(t, timezone.utc).isoformat(timespec="milliseconds")


def _compact(t: float) -> str:
    return datetime.fromtimestamp(t, timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def _write_json_atomic(path: Path, data: Any) -> None:
    tmp = path.with_name(f".{path.name}.tmp-{secrets.token_hex(6)}")
    try:
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(data, fh, indent=2, sort_keys=True)
            fh.write("\n")
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
    except BaseException:
        try:
            tmp.unlink()
        except OSError:
            pass
        raise


def _append_line(path: Path, record: dict) -> None:
    with open(path, "a", encoding="utf-8") as fh:
        fh.write(json.dumps(record, sort_keys=True) + "\n")
        fh.flush()
        os.fsync(fh.fileno())


def _read_json(path: Path) -> tuple[Any, str | None]:
    try:
        return json.loads(path.read_text(encoding="utf-8")), None
    except FileNotFoundError:
        return None, None
    except (OSError, ValueError) as exc:
        return None, f"{path.name} unreadable: {exc}"


def _git(repo: Path, *args: str) -> subprocess.CompletedProcess:
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    try:
        return subprocess.run(["git", "--no-optional-locks", "-C", str(repo), *args], capture_output=True,
                              timeout=SCHEDULER_TIMEOUT_S, env=env)
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise GuardRefused(f"git {' '.join(args)} in {repo} could not run: {exc}") from exc


# --- sbatch arguments -----------------------------------------------------------------------------------------

def _refuse_token(token: Any) -> str:
    if not isinstance(token, str):
        raise GuardRefused(f"sbatch argument {token!r} is not a string")
    if _CONTROL.search(token):
        raise GuardRefused(f"sbatch argument {token!r} holds a control character")
    if token == "--":
        raise GuardRefused("'--' among the sbatch options: an option terminator would change where sbatch's options "
                           "end")
    if token == ":":
        raise GuardRefused("':' among the sbatch options: heterogeneous-job separators are refused")
    return token


def split_user_sbatch_args(args: Sequence[str]) -> tuple[list[str], str | None, list[str]]:
    """User sbatch options, then the first pathname and its literal arguments; refuses any option off the allowlist.

    Only the complete canonical long names of ``ALLOWED_USER_OPTIONS`` pass, as ``--opt=value`` or ``--opt value``,
    unchanged. Everything after the first pathname is the script's own arguments and is not read here.
    """
    options: list[str] = []
    i = 0
    args = list(args)
    while i < len(args):
        token = _refuse_token(args[i])
        if token.startswith("-") and token != "-":
            match = _LONG.fullmatch(token)
            if match is None:
                raise GuardRefused(f"sbatch option {token!r}: short options and bundled or attached short forms are "
                                   "refused; give the complete canonical long name")
            name, has_value, value = f"--{match.group(1)}", match.group(2) is not None, match.group(3)
            if name in OWNED_SBATCH_OPTIONS:
                raise GuardRefused(f"sbatch option {name} is set by the guard, which owns the script, the job's "
                                   "context and identity, its array, its working directory and the submission "
                                   "protocol (exactly one --export=ALL)")
            if name not in ALLOWED_USER_OPTIONS:
                like = sorted(n for n in ALLOWED_USER_OPTIONS | OWNED_SBATCH_OPTIONS
                              if n.startswith("--") and n.startswith(name))
                what = f"an abbreviation of {', '.join(like)}" if like else "not on the guard's allowlist"
                raise GuardRefused(f"sbatch option {name} is {what}; allowed: {', '.join(sorted(ALLOWED_USER_OPTIONS))}"
                                   ", each as its complete name")
            if has_value and value == "":
                raise GuardRefused(f"sbatch option {token!r} has an empty value")
            if name in USER_VALUE_OPTIONS and not has_value:
                if i + 1 >= len(args):
                    raise GuardRefused(f"sbatch option {name} needs a value")
                value_token = _refuse_token(args[i + 1])
                if value_token == "" or value_token.startswith("-"):
                    raise GuardRefused(f"sbatch option {name} is followed by {value_token!r}, not a value; give "
                                       f"{name}=VALUE")
                options += [token, value_token]
                i += 2
                continue
            options.append(token)
            i += 1
            continue
        return options, token, args[i + 1:]
    return options, None, []


def check_user_sbatch_args(args: Sequence[str]) -> None:
    """Refuse user sbatch arguments that are not allowlisted options (stdin mode: no script pathname either)."""
    _, script, _ = split_user_sbatch_args(args)
    if script is not None:
        raise GuardRefused(f"{script!r} among the sbatch options would be read as a script pathname; the guard "
                           "submits the captured script itself")


@dataclass(frozen=True)
class _Option:
    name: str
    value: str | None


def _split_internal(argv: Sequence[str], *, stdin_mode: bool) -> tuple[list[_Option], str | None]:
    """Internal sbatch arguments: long options before the script; in stdin mode, no pathname at all."""
    options: list[_Option] = []
    argv = list(argv)
    i = 0
    while i < len(argv):
        token = _refuse_token(argv[i])
        if token.startswith("-") and token != "-":
            match = _LONG.fullmatch(token)
            if match is None:
                raise GuardRefused(f"internal sbatch argument {token!r}: internal options use canonical long names")
            name = match.group(1)
            owned = [p for p in _PROTOCOL_NAMES if p.startswith(name)]
            if owned:
                raise GuardRefused(f"internal sbatch option --{name} would change the submission protocol "
                                   f"(--{', --'.join(owned)}), which the guard sets itself")
            if match.group(2) is None and f"--{name}" in USER_VALUE_OPTIONS:
                if i + 1 >= len(argv):
                    raise GuardRefused(f"sbatch option --{name} needs a value")
                options.append(_Option(name, _refuse_token(argv[i + 1])))
                i += 2
                continue
            options.append(_Option(name, match.group(3)))
            i += 1
            continue
        if stdin_mode:
            raise GuardRefused(f"{token!r} among the sbatch arguments: sbatch would read it as the script, but the "
                               "guard submits the captured script on stdin (internal options use --name=value)")
        return options, token
    if not stdin_mode:
        raise GuardRefused("a submission by path names its script after the sbatch options")
    return options, None


def _count_indices(spec: str) -> int:
    total = 0
    for part in spec.split(","):
        match = _INDEX_RANGE.fullmatch(part)
        if match is None:
            raise ValueError(f"index range {part!r} not understood")
        lo = int(match.group(1))
        hi = int(match.group(2)) if match.group(2) is not None else lo
        step = int(match.group(3)) if match.group(3) is not None else 1
        if hi < lo or step < 1 or hi - lo > 1_000_000:
            raise ValueError(f"index range {part!r} not understood")
        total += (hi - lo) // step + 1
    return total


def _task_count(job_id: str) -> int:
    match = _ARRAY_TASK.fullmatch(job_id)
    if match is None:
        raise GuardError(f"array task id {job_id!r} from the queue is not understood")
    if match.group(2) is not None:
        return 1
    try:
        return _count_indices(match.group(3).split("%", 1)[0])
    except ValueError as exc:
        raise GuardError(f"array task id {job_id!r} from the queue: {exc}") from exc


def _check_array(kind: str, slots: int, options: list[_Option]) -> None:
    arrays = [o for o in options if o.name == "array"]
    if kind != "array":
        if arrays:
            raise GuardRefused(f"a {kind} submission takes no --array")
        return
    if len(arrays) != 1 or arrays[0].value is None:
        raise GuardRefused("an array submission takes exactly one --array=INDICES[%THROTTLE]")
    indices, sep, throttle = arrays[0].value.partition("%")
    try:
        count = _count_indices(indices)
    except ValueError as exc:
        raise GuardRefused(f"--array={arrays[0].value}: {exc}") from exc
    if sep:
        if not _JOB_ID.fullmatch(throttle) or int(throttle) != slots:
            raise GuardRefused(f"--array={arrays[0].value}: an array reserves exactly its throttle, {slots} slot(s)")
    elif count > slots:
        raise GuardRefused(f"--array={arrays[0].value} has no throttle, so its {count} tasks may all run; it "
                           f"reserves {slots} slot(s)")


def script_directive_problems(script: bytes) -> list[str]:
    """``#SBATCH`` header lines that would set the protocol, the array or the working directory."""
    problems = []
    for number, line in enumerate(script.decode("utf-8", errors="replace").splitlines(), start=1):
        stripped = line.strip()
        if not stripped:
            continue
        if not stripped.startswith("#"):
            break                               # sbatch reads directives only before the first command
        if not stripped.startswith("#SBATCH"):
            continue
        for token in stripped[len("#SBATCH"):].split():
            if token.startswith("#"):
                break
            if token.startswith("--"):
                name = token[2:].split("=", 1)[0]
                if name and any(p.startswith(name) for p in _DIRECTIVE_NAMES):
                    problems.append(f"line {number}: directive {token} sets what the guard owns")
            elif token.startswith("-") and len(token) > 1:
                if token[1] in _DIRECTIVE_SHORT or set(token[2:]) & {"W", "Q"}:
                    problems.append(f"line {number}: directive {token} sets what the guard owns")
    return problems


# --- path mode ------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class PathScript:
    """A path-mode script the guard verified in the bridge clone."""

    path: Path
    name: str
    rel: str
    kind: str
    slots: int
    sha256: str
    bridge_repo: Path
    bridge_commit: str


def default_bridge_repo() -> Path:
    """The bridge clone beside this checkout, ``../bbridge`` (OR-044)."""
    return HFSS_CHECKOUT.parent / "bbridge"


def check_path_mode_script(script: str | os.PathLike, bridge_repo: str | os.PathLike | None = None) -> PathScript:
    """The script, if its real path is a listed script of the bridge clone's folder, tracked and equal to HEAD."""
    bridge = Path(bridge_repo) if bridge_repo is not None else default_bridge_repo()
    try:
        bridge_real = bridge.resolve(strict=True)
    except (OSError, RuntimeError) as exc:
        raise GuardRefused(f"the bridge clone {bridge} cannot be reached: {exc}") from exc
    folder = bridge_real / BRIDGE_SCRIPTS_DIR
    try:
        real = Path(script).resolve(strict=True)
    except (OSError, RuntimeError) as exc:
        raise GuardRefused(f"script {script} does not exist: {exc}") from exc
    if real.parent != folder:
        raise GuardRefused(f"script {script} resolves to {real}, outside the bridge clone's {folder}")
    if real.name not in PATH_MODE_SCRIPTS:
        raise GuardRefused(f"script {real.name} is not one of {', '.join(PATH_MODE_SCRIPTS)}, the scripts the guard "
                           "submits by path")
    if not real.is_file():
        raise GuardRefused(f"script {real} is not a file")
    top = _git(bridge_real, "rev-parse", "--show-toplevel")
    if top.returncode != 0 or Path(top.stdout.decode().strip()).resolve() != bridge_real:
        raise GuardRefused(f"{bridge_real} is not the root of a git repository, so it cannot be the bridge clone")
    rel = real.relative_to(bridge_real).as_posix()
    if _git(bridge_real, "ls-files", "--error-unmatch", "--", rel).returncode != 0:
        raise GuardRefused(f"script {rel} is not tracked in the bridge clone {bridge_real}")
    blob = _git(bridge_real, "cat-file", "blob", f"HEAD:{rel}")
    if blob.returncode != 0:
        raise GuardRefused(f"script {rel} has no blob at the bridge clone's HEAD: "
                           f"{blob.stderr.decode(errors='replace').strip()}")
    data = real.read_bytes()
    if data != blob.stdout:
        raise GuardRefused(f"script {rel} differs from its HEAD blob in the bridge clone; commit or restore it")
    head = _git(bridge_real, "rev-parse", "HEAD")
    if head.returncode != 0:
        raise GuardRefused(f"the bridge clone {bridge_real} has no HEAD commit")
    kind, slots = PATH_MODE_SCRIPTS[real.name]
    return PathScript(path=real, name=real.name, rel=rel, kind=kind, slots=slots,
                      sha256=hashlib.sha256(data).hexdigest(), bridge_repo=bridge_real,
                      bridge_commit=head.stdout.decode().strip())


# --- the guard ------------------------------------------------------------------------------------------------

def _check_cap(cap: Any) -> int:
    if isinstance(cap, bool) or not isinstance(cap, int) or not 1 <= cap <= CEILING:
        raise GuardRefused(f"the cap must be an integer from 1 to {CEILING}, got {cap!r}; a lower cap is an extra "
                           f"restriction and never raises the ceiling of {CEILING}")
    return cap


def _check_kind_slots(kind: Any, slots: Any) -> None:
    if kind not in KINDS:
        raise GuardRefused(f"kind {kind!r} is not one of {', '.join(KINDS)}")
    if isinstance(slots, bool) or not isinstance(slots, int):
        raise GuardRefused(f"slots must be an integer, got {slots!r}")
    if kind in ZERO_SLOT_KINDS and slots != 0:
        raise GuardRefused(f"a {kind} job takes no solve slot (--slots 0), got {slots}")
    if kind in ("single", "reexport") and slots != 1:
        raise GuardRefused(f"a {kind} job is one solve (--slots 1), got {slots}")
    if kind == "array" and not 1 <= slots <= CEILING:
        raise GuardRefused(f"an array reserves 1 to {CEILING} slots, its throttle, got {slots}")


def _check_context(context: Any) -> dict:
    if not isinstance(context, dict):
        raise GuardRefused(f"the submission context must be a mapping, got {type(context).__name__}")
    try:
        return json.loads(json.dumps(context))
    except (TypeError, ValueError) as exc:
        raise GuardRefused(f"the submission context is not JSON: {exc}") from exc


def _scheduler_override(name: str) -> bool:
    return name.startswith("SBATCH_") or name in SCHEDULER_ENV_ALIASES


def _check_env_extras(env: Mapping[str, str] | None) -> dict[str, str]:
    if env is None:
        return {}
    extras = {}
    for name, value in env.items():
        if not isinstance(name, str) or not _ENV_NAME.fullmatch(name):
            raise GuardRefused(f"environment name {name!r} is not a variable name")
        if not isinstance(value, str) or "\0" in value:
            raise GuardRefused(f"environment value of {name} must be a string without NUL")
        if name == RESERVATION_ENV:
            raise GuardRefused(f"{RESERVATION_ENV} is set by the guard to the reservation it creates")
        if _scheduler_override(name):
            raise GuardRefused(f"{name} would change sbatch's options; the guard removes such variables")
        extras[name] = value
    return extras


def _parse_job_id(stdout: str) -> str | None:
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    if len(lines) != 1:
        return None
    job = lines[0].split(";", 1)[0]
    return job if _JOB_ID.fullmatch(job) else None


@dataclass
class _Group:
    key: str
    is_array: bool
    names: set = field(default_factory=set)
    comments: set = field(default_factory=set)
    tasks: int = 0
    active: int = 0


def _group_entries(entries: Sequence[QueueEntry]) -> dict[str, _Group]:
    groups: dict[str, _Group] = {}
    for e in entries:
        key = e.array_parent or e.job_id
        group = groups.setdefault(key, _Group(key, e.array_parent is not None))
        n = _task_count(e.job_id) if e.array_parent else 1
        group.tasks += n
        if e.state not in PENDING_QUEUE_STATES:
            group.active += n
        group.names.add(e.name)
        group.comments.add(e.comment)
    return groups


@dataclass
class _Reconciled:
    usage: Usage
    groups: dict[str, _Group]
    counts: dict[str, int]                  # reservation id -> the slots counted for it (its job's, if bound)
    blocking: list[str]                     # uncertain reservations holding slots
    held: dict[str, int]                    # reservation id -> held leases


class Guard:
    """The shared admission in ``root`` (resolved once), with its scheduler."""

    def __init__(self, root: str | os.PathLike, scheduler: Scheduler, *, lock_timeout_s: float = 60.0,
                 clock: Callable[[], float] = time.time):
        try:
            resolved = Path(root).resolve(strict=True)
        except (OSError, RuntimeError) as exc:
            raise GuardRefused(f"the guard directory {root} cannot be reached ({exc}); it is created once, by hand "
                               "(mkdir -p), never by a client") from exc
        if not resolved.is_dir():
            raise GuardRefused(f"the guard directory {resolved} is not a directory")
        if isinstance(lock_timeout_s, bool) or not isinstance(lock_timeout_s, (int, float)) \
                or not (math.isfinite(lock_timeout_s) and lock_timeout_s >= 0):
            raise GuardRefused(f"lock_timeout_s must be a finite nonnegative number, got {lock_timeout_s!r}")
        self.root = resolved
        self.scheduler = scheduler
        self.lock_timeout_s = float(lock_timeout_s)
        self.clock = clock
        self.last_reservation_id: str | None = None
        self._depth = 0
        self._nonce: str | None = None

    # paths
    @property
    def _lock(self) -> Path:
        return self.root / "lock"

    @property
    def _recovery(self) -> Path:
        return self.root / "recover.lock"

    @property
    def _reservations(self) -> Path:
        return self.root / "reservations"

    @property
    def _leases(self) -> Path:
        return self.root / "leases"

    @property
    def _events(self) -> Path:
        return self.root / "events.log"

    def _now(self) -> str:
        return _utc(self.clock())

    # the lock
    @contextmanager
    def locked(self) -> Iterator[None]:
        """Hold the guard lock (re-entrant within this guard); released in ``finally``, only if still ours."""
        if self._depth:
            self._depth += 1
            try:
                yield
            finally:
                self._depth -= 1
            return
        nonce = self._acquire()
        self._nonce, self._depth = nonce, 1
        try:
            yield
        finally:
            self._depth, self._nonce = 0, None
            self._release(nonce)

    def _acquire(self) -> str:
        deadline = time.monotonic() + self.lock_timeout_s
        while True:
            try:
                os.mkdir(self._lock)
                break
            except FileExistsError:
                pass
            except OSError as exc:
                raise GuardRefused(f"cannot create the guard lock {self._lock}: {exc}") from exc
            if time.monotonic() >= deadline:
                raise GuardRefused(self._lock_timeout_message(self.lock_info()))
            time.sleep(LOCK_POLL_S)
        nonce = secrets.token_hex(16)
        try:
            _write_json_atomic(self._lock / _OWNER, {"nonce": nonce, "utc": self._now(), **_owner_identity()})
        except BaseException:
            # The directory is ours and holds no complete record: remove it, or leave it for recover-lock.
            try:
                children = list(self._lock.iterdir())
            except OSError:
                children = []
            for child in children:
                try:
                    child.unlink()
                except OSError:
                    pass
            try:
                os.rmdir(self._lock)
            except OSError:
                pass
            raise
        return nonce

    def _release(self, nonce: str) -> None:
        owner, _ = _read_json(self._lock / _OWNER)
        if not isinstance(owner, dict) or owner.get("nonce") != nonce:
            log.warning("guard lock %s is no longer this process's (nonce %s); left in place", self._lock, nonce)
            return
        trash = self.root / f".lock-released-{nonce}"
        try:
            os.rename(self._lock, trash)
        except OSError as exc:
            log.warning("guard lock %s could not be released: %s", self._lock, exc)
            return
        for child in trash.iterdir():
            try:
                child.unlink()
            except OSError:
                pass
        try:
            os.rmdir(trash)
        except OSError:
            pass

    def _lock_timeout_message(self, info: dict) -> str:
        token = info.get("token")
        return (f"the guard lock {self._lock} stayed held for {self.lock_timeout_s:g} s and is never broken by age; "
                f"observation (no lock taken): {json.dumps(info, sort_keys=True, default=str)}. If its owner is "
                "alive, wait and retry. If it is dead, or no owner record exists, stop every guard client, verify "
                "the owner is dead (or the whole workflow quiescent), then one operator runs: python -m "
                f"bbsim.licence_guard --root {self.root} recover-lock --expected-lock {token} --confirm-quiescent "
                "'<what was verified>'")

    def lock_info(self) -> dict:
        """The lock as it is now, read without taking it: its identity token and owner record, or their absence."""
        info: dict[str, Any] = {"root": str(self.root), "lock": str(self._lock), "held": False, "token": None,
                                "owner": None, "owner_error": None, "owner_status": None, "owner_status_reason": None}
        recovery, recovery_error = _read_json(self._recovery / _OWNER)
        info["recovery"] = ({"held": True, "owner": recovery, "error": recovery_error}
                            if self._recovery.exists() else {"held": False})
        try:
            st = os.stat(self._lock)
            entries = sorted(os.listdir(self._lock))
        except FileNotFoundError:
            return info
        owner, owner_error = _read_json(self._lock / _OWNER)
        nonce = owner.get("nonce") if isinstance(owner, dict) and isinstance(owner.get("nonce"), str) else None
        status, reason = _owner_status(owner)
        info.update(held=True, inode=st.st_ino, ctime_ns=st.st_ctime_ns, entries=entries,
                    token=f"{st.st_ino}-{st.st_ctime_ns}-{nonce or 'no-owner'}", owner=owner,
                    owner_error=owner_error, owner_status=status, owner_status_reason=reason)
        return info

    def recover_lock(self, expected_lock: str, *, confirm_quiescent: str) -> None:
        """Move a crashed owner's lock to ``lock-evidence/`` after rechecking it; reservations and leases untouched."""
        reason = confirm_quiescent.strip() if isinstance(confirm_quiescent, str) else ""
        if not reason:
            raise GuardRefused("recover-lock needs --confirm-quiescent REASON: what was verified (every guard client "
                               "stopped; the owner dead, or the whole workflow quiescent)")
        if not isinstance(expected_lock, str) or not expected_lock:
            raise GuardRefused("recover-lock needs --expected-lock TOKEN from status --lock-info")
        try:
            os.mkdir(self._recovery)
        except FileExistsError:
            other, _ = _read_json(self._recovery / _OWNER)
            raise GuardRefused(f"another recovery holds {self._recovery} (record {other!r}); recoveries never compete. "
                               "If that operator's process is verified gone, its directory is removed by hand after "
                               "checking, then recovery is run once") from None
        except OSError as exc:
            raise GuardRefused(f"cannot create {self._recovery}: {exc}") from exc
        try:
            _write_json_atomic(self._recovery / _OWNER, {"expected_lock": expected_lock, "utc": self._now(),
                                                         **_identity()})
            info = self.lock_info()
            if not info["held"]:
                raise GuardRefused(f"there is no lock at {self._lock} to recover")
            if info["token"] != expected_lock:
                raise GuardRefused(f"the lock changed: expected {expected_lock}, found {info['token']}; nothing was "
                                   "changed. Read status --lock-info again after stopping every client")
            if info["owner_status"] == "alive":
                raise GuardRefused(f"the lock's owner is alive: {info['owner_status_reason']}; nothing was changed")
            evidence_root = self.root / "lock-evidence"
            evidence_root.mkdir(exist_ok=True)
            dest = evidence_root / f"{_compact(self.clock())}-{re.sub(r'[^A-Za-z0-9_.-]', '_', expected_lock)}"
            if dest.exists():
                raise GuardRefused(f"{dest} exists already")
            observation = {k: info.get(k) for k in ("inode", "ctime_ns", "entries", "owner", "owner_error",
                                                    "owner_status", "owner_status_reason")}
            event = {"utc": self._now(), "event": "lock_recovered", "token": expected_lock, "assertion": reason,
                     "observation": observation, "evidence": str(dest), "by": _identity()}
            self._append_event(event)               # logged before normal acquisition is possible again
            again = self.lock_info()
            if again["token"] != expected_lock:
                self._append_event({"utc": self._now(), "event": "lock_recovery_aborted", "token": expected_lock,
                                    "found": again["token"], "by": _identity()})
                raise GuardRefused(f"the lock changed during recovery (now {again['token']}); it was left in place")
            os.rename(self._lock, dest)
            _write_json_atomic(dest / "recovery.json", event)
        finally:
            try:
                (self._recovery / _OWNER).unlink()
            except OSError:
                pass
            try:
                os.rmdir(self._recovery)
            except OSError:
                pass

    # records
    def _append_event(self, record: dict) -> None:
        _append_line(self._events, record)

    def _new_id(self) -> str:
        while True:
            rid = f"r{_compact(self.clock())}-{secrets.token_hex(4)}"
            if not (self._reservations / f"{rid}.json").exists():
                return rid

    def _read_reservations(self) -> dict[str, dict]:
        records: dict[str, dict] = {}
        if not self._reservations.is_dir():
            return records
        for path in sorted(self._reservations.iterdir()):
            if path.name.startswith(".") or path.suffix != ".json":
                continue
            rec, error = _read_json(path)
            ok = (error is None and isinstance(rec, dict) and rec.get("id") == path.stem
                  and rec.get("state") in STATES and rec.get("kind") in (*KINDS, DIRECT_KIND)
                  and isinstance(rec.get("slots"), int) and not isinstance(rec.get("slots"), bool)
                  and rec["slots"] >= 0 and isinstance(rec.get("events"), list)
                  and (rec["state"] != "bound" or (isinstance(rec.get("job_id"), str)
                                                   and _JOB_ID.fullmatch(rec["job_id"]))))
            if not ok:
                raise GuardError(f"reservation record {path} is unreadable or malformed ({error or 'fields'}); "
                                 "admission is refused until it is repaired")
            records[rec["id"]] = rec
        return records

    def _read_leases(self) -> list[dict]:
        leases = []
        if not self._leases.is_dir():
            return leases
        for path in sorted(self._leases.iterdir()):
            if path.name.startswith(".") or path.suffix != ".json":
                continue
            rec, error = _read_json(path)
            if error is not None or not isinstance(rec, dict) or not isinstance(rec.get("state"), str) \
                    or not (rec.get("reservation_id") is None or isinstance(rec.get("reservation_id"), str)):
                raise GuardError(f"lease record {path} is unreadable or malformed ({error or 'fields'}); admission "
                                 "is refused until it is repaired")
            leases.append({"id": rec.get("id", path.stem), "reservation_id": rec.get("reservation_id"),
                           "state": rec["state"]})
        return leases

    def _change(self, record: dict, state: str, *, event: str, **details: Any) -> None:
        record["state"] = state
        record["events"].append({"utc": self._now(), "event": event, **details})
        _write_json_atomic(self._reservations / f"{record['id']}.json", record)
        self._append_event({"utc": self._now(), "event": event, "reservation": record["id"], "state": state,
                            "kind": record["kind"], "slots": record["slots"], **details, "by": _identity()})

    def _create(self, *, kind: str, slots: int, cap: int, context: dict, sbatch_argv: list[str] | None) -> dict:
        self._reservations.mkdir(exist_ok=True)
        rid = self._new_id()
        record = {"id": rid, "kind": kind, "slots": slots, "cap": cap, "state": "pending", "job_id": None,
                  "context": context, "created_utc": self._now(), "events": [], "sbatch_argv": sbatch_argv,
                  "throttle_target": None}
        self._change(record, "pending", event="reserved", cap=cap, context=context)
        self.last_reservation_id = rid
        return record

    # reconciliation
    def _drain(self, record: dict, active: int, held: int) -> None:
        target = record.get("throttle_target")
        if target is None:
            return
        new = min(record["slots"], max(target, active, held))
        if new == record["slots"] and new != target:
            return
        record["slots"] = new
        if new == target:
            record["throttle_target"] = None
        self._change(record, record["state"], event="throttle_drained" if new == target else "throttle_draining",
                     target=target, active_tasks=active, held_leases=held)

    def _reconcile(self) -> _Reconciled:
        groups = _group_entries(self.scheduler.queue())
        records = self._read_reservations()
        leases = self._read_leases()
        held: dict[str, int] = {}
        lease_states: dict[str | None, list[str]] = {}
        for lease in leases:
            lease_states.setdefault(lease["reservation_id"], []).append(lease["state"])
            if lease["state"] != LEASE_RELEASED and lease["reservation_id"] is not None:
                held[lease["reservation_id"]] = held.get(lease["reservation_id"], 0) + 1
        throttles: dict[str, int | None] = {}

        def scheduler_count(group: _Group, managed_slots: bool) -> int:
            if group.is_array:
                if managed_slots or any(c.startswith(MANAGED_COMMENT_PREFIX) for c in group.comments):
                    if group.key not in throttles:
                        throttles[group.key] = self.scheduler.array_throttle(group.key)
                    throttle = throttles[group.key]
                    return throttle if throttle else group.tasks
                return 0
            return 1 if (group.names & SOLVE_JOB_NAMES or managed_slots) else 0

        total, sources, uncertain, unmanaged, blocking = 0, [], [], [], []
        counts: dict[str, int] = {}
        bound: dict[str, list[dict]] = {}
        for rid, rec in sorted(records.items(), key=lambda kv: (kv[1].get("created_utc", ""), kv[0])):
            state = rec["state"]
            if state == "released":
                continue
            if state == "bound":
                bound.setdefault(rec["job_id"], []).append(rec)
                continue
            if rec["kind"] == DIRECT_KIND and state == "pending":
                states = lease_states.get(rid, [])
                if states and all(s == LEASE_RELEASED for s in states):
                    self._change(rec, "released", event="released_with_lease", leases=len(states))
                    continue
            n = max(rec["slots"], held.get(rid, 0))
            counts[rid] = n
            total += n
            if n:
                sources.append({"source": "reservation", "id": rid, "kind": rec["kind"], "state": state,
                                "job_id": None, "count": n,
                                "detail": f"{state}, held leases {held.get(rid, 0)}"})
            if state == "uncertain":
                uncertain.append(rid)
                if rec["slots"]:
                    blocking.append(rid)
        consumed = set()
        for job, recs in bound.items():
            group = groups.get(job)
            if group is not None:
                consumed.add(job)
                for rec in recs:
                    self._drain(rec, group.active, held.get(rec["id"], 0))
                slots = max(rec["slots"] for rec in recs)
                leases_held = sum(held.get(rec["id"], 0) for rec in recs)
                n = max(slots, scheduler_count(group, slots > 0), leases_held)
                total += n
                for rec in recs:
                    counts[rec["id"]] = n
                if n:
                    sources.append({"source": "reservation", "id": recs[0]["id"], "kind": recs[0]["kind"],
                                    "state": "bound", "job_id": job, "count": n,
                                    "detail": f"queued, slots {slots}, held leases {leases_held}",
                                    "reservations": [rec["id"] for rec in recs]})
                continue
            terminal = self.scheduler.terminal_state(job)
            remaining = []
            for rec in recs:
                h = held.get(rec["id"], 0)
                if terminal is not None and h == 0:
                    self._change(rec, "released", event="released_terminal", job_state=terminal)
                    continue
                self._drain(rec, 0, h)
                remaining.append(rec)
            if not remaining:
                continue
            leases_held = sum(held.get(rec["id"], 0) for rec in remaining)
            n = max(max(rec["slots"] for rec in remaining), leases_held)
            total += n
            for rec in remaining:
                counts[rec["id"]] = n
            if n:
                sources.append({"source": "reservation", "id": remaining[0]["id"], "kind": remaining[0]["kind"],
                                "state": "bound", "job_id": job, "count": n,
                                "detail": f"not queued; accounting {terminal or 'not terminal or unknown'}; "
                                          f"held leases {leases_held}",
                                "reservations": [rec["id"] for rec in remaining]})
        for key, group in groups.items():
            if key in consumed:
                continue
            n = scheduler_count(group, False)
            if n:
                total += n
                unmanaged.append(key)
                sources.append({"source": "array" if group.is_array else "job", "id": key, "count": n,
                                "names": sorted(group.names), "comments": sorted(group.comments),
                                "detail": "submitted outside the guard"})
        live = {rid for rid, rec in records.items() if rec["state"] != "released"}
        for lease in leases:
            if lease["state"] != LEASE_RELEASED and lease["reservation_id"] not in live:
                total += 1
                sources.append({"source": "lease", "id": lease["id"], "count": 1,
                                "detail": f"held ({lease['state']}) and naming no live reservation "
                                          f"({lease['reservation_id']!r})"})
        usage = Usage(total=total, sources=sources, uncertain=uncertain, unmanaged=unmanaged)
        return _Reconciled(usage=usage, groups=groups, counts=counts, blocking=blocking, held=held)

    def _reconcile_refusing(self) -> _Reconciled:
        try:
            return self._reconcile()
        except GuardRefused:
            raise
        except Exception as exc:                # a scheduler or record error, before any reservation
            raise GuardRefused(f"reconciliation failed, so nothing is admitted: {exc}") from exc

    def _admit(self, state: _Reconciled, slots: int, cap: int, what: str) -> None:
        if slots == 0:
            return
        if state.blocking:
            raise GuardRefused(f"{what}: reservation(s) {', '.join(state.blocking)} are uncertain and hold slots, so "
                               "nothing more is admitted until each is resolved: check squeue, then reconcile "
                               "--resolve ID=released|job:JOBID --reason TEXT")
        limit = min(CEILING, cap)
        if state.usage.total + slots > limit:
            raise GuardRefused(f"{what}: {state.usage.total} in use + {slots} requested exceeds {limit}, the smaller "
                               f"of the ceiling {CEILING} (OR-045) and this request's cap {cap}")

    def usage(self) -> Usage:
        """The reconciled usage, under the lock."""
        with self.locked():
            return self._reconcile().usage

    def _check_root_claims(self, env: Mapping[str, str] | None, context: dict) -> None:
        claims = []
        inherited = os.environ.get(GUARD_ENV)
        if inherited:
            claims.append((f"the inherited {GUARD_ENV}", inherited))
        if env is not None and env.get(GUARD_ENV):
            claims.append((f"the explicit {GUARD_ENV}", env[GUARD_ENV]))
        if context.get("guard_dir") is not None:
            claims.append(("the plan's guard_dir", context["guard_dir"]))
        for what, value in claims:
            if not isinstance(value, str) or Path(value).resolve() != self.root:
                raise GuardRefused(f"{what} {value!r} is not this guard ({self.root}); a second guard cannot "
                                   "bypass the shared cap")

    # admission
    def submit(self, *, kind: str, slots: int, cap: int, argv: list[str], stdin: bytes | None = None,
               env: Mapping[str, str] | None = None, context: dict) -> str:
        """Reserve, run ``sbatch --parsable --export=ALL <argv>`` and bind the returned job id, under one lock."""
        _check_kind_slots(kind, slots)
        cap = _check_cap(cap)
        context = _check_context(context)
        self._check_root_claims(env, context)
        extras = _check_env_extras(env)
        if stdin is not None and not isinstance(stdin, (bytes, bytearray)):
            raise GuardRefused("stdin must be the captured script's bytes")
        options, script = _split_internal(argv, stdin_mode=stdin is not None)
        _check_array(kind, slots, options)
        if stdin is not None:
            data = bytes(stdin)
        else:
            try:
                data = Path(script).read_bytes()
            except OSError as exc:
                raise GuardRefused(f"script {script} cannot be read: {exc}") from exc
        problems = script_directive_problems(data)
        if problems:
            raise GuardRefused("the script's #SBATCH header sets what the guard owns: " + "; ".join(problems))
        full = ["--parsable", "--export=ALL", *argv]
        with self.locked():
            state = self._reconcile_refusing()
            self._admit(state, slots, cap, f"a {kind} submission of {slots} slot(s)")
            record = self._create(kind=kind, slots=slots, cap=cap, context=context, sbatch_argv=full)
            rid = record["id"]
            sub_env = {k: v for k, v in os.environ.items() if not _scheduler_override(k) and k != RESERVATION_ENV}
            sub_env.update(extras)
            sub_env[GUARD_ENV] = str(self.root)
            sub_env[RESERVATION_ENV] = rid
            try:
                returncode, stdout = self.scheduler.submit(list(full), data if stdin is not None else None, sub_env)
            except Exception as exc:
                self._change(record, "uncertain", event="submit_error", error=str(exc))
                raise GuardError(self._uncertain_message(record, f"sbatch failed: {exc}")) from exc
            job = _parse_job_id(stdout) if returncode == 0 else None
            if job is None:
                self._change(record, "uncertain", event="submit_ambiguous", returncode=returncode,
                             stdout=stdout[-2000:])
                raise GuardError(self._uncertain_message(record, f"sbatch exited {returncode} with stdout "
                                                                 f"{stdout.strip()[-200:]!r}"))
            record["job_id"] = job
            self._change(record, "bound", event="bound", job_id=job)
            return job

    def _uncertain_message(self, record: dict, why: str) -> str:
        return (f"the submission's result is ambiguous ({why}): reservation {record['id']} stays uncertain and keeps "
                f"its {record['slots']} slot(s) counted, and no further slot is admitted until it is resolved. Check "
                f"squeue -u $USER (and sacct) for the job, then run python -m bbsim.licence_guard --root {self.root} "
                f"reconcile --resolve {record['id']}=job:JOBID --reason TEXT, or {record['id']}=released if no job "
                "exists")

    def reserve_direct(self, *, cap: int, context: dict) -> str:
        """One slot for an AEDT session started outside a guarded job; released only with its cleaned lease."""
        cap = _check_cap(cap)
        context = _check_context(context)
        self._check_root_claims(None, context)
        with self.locked():
            state = self._reconcile_refusing()
            self._admit(state, 1, cap, "a direct session's slot")
            return self._create(kind=DIRECT_KIND, slots=1, cap=cap, context=context, sbatch_argv=None)["id"]

    def change_throttle(self, parent_id: str, throttle: int, *, cap: int) -> None:
        """Change a guarded array's throttle: an increase is admitted like a submission, a decrease drains."""
        cap = _check_cap(cap)
        if not isinstance(parent_id, str) or not _JOB_ID.fullmatch(parent_id):
            raise GuardRefused(f"the array's parent job id must be a plain number, got {parent_id!r}")
        if isinstance(throttle, bool) or not isinstance(throttle, int) or not 1 <= throttle <= CEILING:
            raise GuardRefused(f"the throttle must be an integer from 1 to {CEILING}, got {throttle!r}")
        with self.locked():
            state = self._reconcile_refusing()
            records = [r for r in self._read_reservations().values()
                       if r["state"] == "bound" and r["job_id"] == parent_id and r["kind"] == "array"]
            if len(records) != 1:
                raise GuardRefused(f"no array reservation is bound to job {parent_id}" if not records else
                                   f"several array reservations are bound to job {parent_id}")
            record = records[0]
            group = state.groups.get(parent_id)
            if group is None:
                raise GuardRefused(f"array {parent_id} is not in the queue")
            previous = record["slots"]
            if throttle > previous:
                counted = state.counts.get(record["id"], previous)
                self._admit(state, max(0, throttle - counted), cap, f"raising array {parent_id}'s throttle to "
                                                                     f"{throttle}")
                record["slots"], record["throttle_target"] = throttle, None
                self._change(record, "bound", event="throttle_increase", previous=previous, throttle=throttle)
            elif throttle < previous:
                held = state.held.get(record["id"], 0)
                record["slots"] = max(throttle, group.active, held)
                record["throttle_target"] = throttle if record["slots"] > throttle else None
                self._change(record, "bound", event="throttle_decrease", previous=previous, throttle=throttle,
                             active_tasks=group.active, held_leases=held)
            else:
                return
            try:
                self.scheduler.set_throttle(parent_id, throttle)
            except Exception as exc:
                raise GuardError(f"scontrol could not set array {parent_id}'s throttle to {throttle}: {exc}; the "
                                 f"reservation records {record['slots']} slot(s), the larger count; check scontrol "
                                 "show job before retrying") from exc

    def resolve(self, reservation_id: str, outcome: str, reason: str) -> None:
        """The logged manual resolve of an uncertain (or stuck) reservation: ``released`` or ``job:<id>``."""
        if not isinstance(reason, str) or not reason.strip():
            raise GuardRefused("a resolve needs --reason: what was checked (squeue, sacct)")
        if outcome == "released":
            job = None
        elif isinstance(outcome, str) and outcome.startswith("job:") and _JOB_ID.fullmatch(outcome[4:]):
            job = outcome[4:]
        else:
            raise GuardRefused(f"the outcome must be 'released' or 'job:<numeric id>', got {outcome!r}")
        with self.locked():
            state = self._reconcile_refusing()
            record = self._read_reservations().get(reservation_id)
            if record is None:
                raise GuardRefused(f"there is no reservation {reservation_id!r}")
            if record["state"] == "released":
                raise GuardRefused(f"reservation {reservation_id} is released already")
            held = state.held.get(reservation_id, 0)
            if job is None:
                if held:
                    raise GuardRefused(f"{held} held launch lease(s) name reservation {reservation_id}; a lease is "
                                       "released only after its cleanup is verified, and the reservation with it")
                if record["state"] == "bound" and record["job_id"] in state.groups:
                    raise GuardRefused(f"reservation {reservation_id}'s job {record['job_id']} is still queued")
                self._change(record, "released", event="resolved", outcome=outcome, reason=reason.strip(),
                             previous=record["state"], job_queued=False)
                return
            if record["state"] not in ("pending", "uncertain"):
                raise GuardRefused(f"only a pending or uncertain reservation is bound by a resolve; "
                                   f"{reservation_id} is {record['state']}")
            others = [r["id"] for r in self._read_reservations().values()
                      if r["state"] == "bound" and r["job_id"] == job and r["id"] != reservation_id]
            if others:
                raise GuardRefused(f"job {job} is bound to reservation {', '.join(others)} already")
            previous = record["state"]
            record["job_id"] = job
            self._change(record, "bound", event="resolved", outcome=outcome, reason=reason.strip(), previous=previous,
                         job_queued=job in state.groups)


# --- the two-node lock test -----------------------------------------------------------------------------------

def check_checkout(checkout: Path, expect_commit_file: Path) -> str:
    """This checkout's HEAD, if it is clean and equal to the commit recorded in ``expect_commit_file``."""
    try:
        expected = Path(expect_commit_file).read_text().strip()
    except OSError as exc:
        raise GuardRefused(f"the recorded commit {expect_commit_file} cannot be read: {exc}") from exc
    if not _COMMIT.fullmatch(expected):
        raise GuardRefused(f"{expect_commit_file} holds {expected!r}, not a commit")
    status = _git(Path(checkout), "status", "--porcelain", "--untracked-files=no")
    if status.returncode != 0:
        raise GuardRefused(f"{checkout} is not a git checkout: {status.stderr.decode(errors='replace').strip()}")
    if status.stdout.strip():
        raise GuardRefused(f"the checkout {checkout} has uncommitted changes: "
                           f"{status.stdout.decode(errors='replace').strip()}")
    head = _git(Path(checkout), "rev-parse", "HEAD").stdout.decode().strip()
    if head != expected:
        raise GuardRefused(f"the checkout {checkout} is at {head}, not the recorded {expected}")
    return head


def _read_sequence(test_dir: Path) -> int:
    path = test_dir / LOCK_TEST_SEQUENCE
    if not path.exists():
        return 0
    text = path.read_text().strip()
    if not _JOB_ID.fullmatch(text):
        raise GuardError(f"{path} holds {text!r}, not a count")
    return int(text)


def run_lock_test_worker(root: str | os.PathLike, test_dir: str | os.PathLike, *,
                         expect_commit_file: str | os.PathLike, iterations: int = 200, hosts: int = 2,
                         barrier_timeout_s: float = 300.0, hold_s: float = 0.01, pause_s: float = 0.05,
                         lock_timeout_s: float = 120.0, host: str | None = None, checkout: Path = HFSS_CHECKOUT,
                         clock: Callable[[], float] = time.time) -> int:
    """One worker of the two-node lock test on the production lock in ``root``; artifacts in ``test_dir``.

    Exit 2 when it refuses (the checkout, stale artifacts, the barrier's timeout), 1 on an overlap or a lock error,
    0 after ``iterations`` visits. ``host`` and ``clock`` are injectable for local tests only.
    """
    host = host or socket.gethostname()
    pid = os.getpid()
    job = os.environ.get("SLURM_JOB_ID") or "none"
    test_dir = Path(test_dir)
    tag = f"{re.sub(r'[^A-Za-z0-9_.]', '_', host)}-{pid}"

    def fail(kind: str, message: str, code: int, **extra: Any) -> int:
        try:
            _write_json_atomic(test_dir / LOCK_TEST_FAILURES / f"{kind}-{tag}-{secrets.token_hex(3)}.json",
                               {"kind": kind, "host": host, "pid": pid, "slurm_job": job, "message": message,
                                "utc": _utc(clock()), **extra})
        except OSError as exc:
            print(f"lock-test-worker: cannot record the failure: {exc}", file=sys.stderr)
        print(f"lock-test-worker: {kind}: {message}", file=sys.stderr)
        return code

    try:
        for sub in (LOCK_TEST_READY, LOCK_TEST_FAILURES):
            (test_dir / sub).mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        print(f"lock-test-worker: refused: the test directory {test_dir} cannot be created: {exc}", file=sys.stderr)
        return 2
    try:
        for value, what in ((iterations, "iterations"), (hosts, "hosts")):
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise GuardRefused(f"{what} must be a positive integer, got {value!r}")
        check_checkout(checkout, Path(expect_commit_file))
        guard = Guard(root, _NoScheduler(), lock_timeout_s=lock_timeout_s, clock=clock)
        stale = [n for n in (LOCK_TEST_LOG, LOCK_TEST_SEQUENCE, LOCK_TEST_ACTIVE) if (test_dir / n).exists()]
        if stale:
            raise GuardRefused(f"{test_dir} holds {', '.join(stale)} from an earlier run; use a fresh directory")
    except GuardRefused as exc:
        return fail("refused", str(exc), 2)
    _write_json_atomic(test_dir / LOCK_TEST_READY / f"{tag}.json",
                       {"host": host, "pid": pid, "slurm_job": job, "utc": _utc(clock())})
    deadline = time.monotonic() + barrier_timeout_s
    while True:
        seen = set()
        for path in (test_dir / LOCK_TEST_READY).glob("*.json"):
            rec, _ = _read_json(path)
            if isinstance(rec, dict) and rec.get("slurm_job") == job:
                seen.add(rec.get("host"))
        if len(seen) >= hosts:
            break
        if time.monotonic() >= deadline:
            return fail("refused", f"the start barrier saw {len(seen)} of {hosts} host(s) ({sorted(map(str, seen))}) "
                                   f"after {barrier_timeout_s:g} s", 2)
        time.sleep(0.1)
    log_path = test_dir / LOCK_TEST_LOG
    for visit in range(1, iterations + 1):
        start = time.monotonic()
        try:
            with guard.locked():
                waited = time.monotonic() - start
                try:
                    os.mkdir(test_dir / LOCK_TEST_ACTIVE)
                except FileExistsError:
                    return fail("overlap", f"visit {visit}: {test_dir / LOCK_TEST_ACTIVE} exists while this worker "
                                           "holds the guard lock", 1, visit=visit)
                seq = _read_sequence(test_dir) + 1
                (test_dir / LOCK_TEST_SEQUENCE).write_text(f"{seq}\n")
                record = {"host": host, "pid": pid, "visit": visit, "seq": seq}
                _append_line(log_path, {"event": "enter", **record, "utc": _utc(clock()),
                                        "lock_wait_s": round(waited, 4)})
                if hold_s:
                    time.sleep(hold_s)
                _append_line(log_path, {"event": "leave", **record, "utc": _utc(clock())})
                try:
                    os.rmdir(test_dir / LOCK_TEST_ACTIVE)
                except OSError as exc:
                    return fail("sentinel", f"visit {visit}: the sentinel could not be removed: {exc}", 1,
                                visit=visit)
        except (GuardError, OSError) as exc:
            return fail("lock", f"visit {visit}: {exc}", 1, visit=visit)
        if pause_s:
            time.sleep(pause_s)
    print(f"lock-test-worker: {iterations} visits on {host} (pid {pid})")
    return 0


def check_lock_test(test_dir: str | os.PathLike, *, hosts: int = 2, iterations: int = 200) -> tuple[bool, list[str]]:
    """Judge a lock test's artifacts: distinct hosts, complete visits, a gapless sequence, no overlap, no sentinel."""
    test_dir = Path(test_dir)
    if not test_dir.is_dir():
        raise GuardRefused(f"the lock test directory {test_dir} does not exist")
    problems: list[str] = []
    for path in sorted((test_dir / LOCK_TEST_FAILURES).glob("*.json")) if (test_dir / LOCK_TEST_FAILURES).is_dir() \
            else []:
        rec, error = _read_json(path)
        message = rec.get("message") if isinstance(rec, dict) else error
        problems.append(f"failure record {path.name}: {message}")
    if (test_dir / LOCK_TEST_ACTIVE).exists():
        problems.append(f"the sentinel {test_dir / LOCK_TEST_ACTIVE} is not empty at the end")
    records = []
    log_path = test_dir / LOCK_TEST_LOG
    if not log_path.exists():
        problems.append(f"no log {log_path}")
    else:
        for number, line in enumerate(log_path.read_text().splitlines(), start=1):
            try:
                rec = json.loads(line)
                if not isinstance(rec, dict):
                    raise ValueError("not an object")
                records.append((number, rec))
            except ValueError as exc:
                problems.append(f"log line {number} unreadable: {exc}")
    expected = 1
    open_enter: tuple[int, dict] | None = None
    visits: dict[str, list[tuple[Any, Any]]] = {}
    waits: dict[str, int] = {}

    def who(rec: dict) -> str:
        return f"{rec.get('host')} pid {rec.get('pid')} visit {rec.get('visit')} seq {rec.get('seq')}"

    for number, rec in records:
        event = rec.get("event")
        if event == "enter":
            if open_enter is not None:
                problems.append(f"log line {number}: enter of {who(rec)} before the leave of {who(open_enter[1])}")
            seq = rec.get("seq")
            if seq != expected:
                problems.append(f"log line {number}: sequence {seq!r} where {expected} was expected (a gap or a "
                                "repeat)")
            expected = (seq + 1) if isinstance(seq, int) and not isinstance(seq, bool) else expected + 1
            open_enter = (number, rec)
            if isinstance(rec.get("lock_wait_s"), (int, float)) and rec["lock_wait_s"] > 0.1:
                waits[str(rec.get("host"))] = waits.get(str(rec.get("host")), 0) + 1
        elif event == "leave":
            keys = ("host", "pid", "visit", "seq")
            if open_enter is None or any(rec.get(k) != open_enter[1].get(k) for k in keys):
                problems.append(f"log line {number}: leave of {who(rec)} without its own enter just before")
            else:
                visits.setdefault(str(rec.get("host")), []).append((rec.get("pid"), rec.get("visit")))
            open_enter = None
        else:
            problems.append(f"log line {number}: unknown event {event!r}")
    if open_enter is not None:
        problems.append(f"log line {open_enter[0]}: enter of {who(open_enter[1])} without its leave")
    if len(visits) != hosts:
        problems.append(f"{len(visits)} distinct host(s) completed visits ({', '.join(sorted(visits)) or 'none'}), "
                        f"{hosts} expected")
    for name, done in sorted(visits.items()):
        pids = {pid for pid, _ in done}
        if len(pids) != 1:
            problems.append(f"host {name}: {len(pids)} workers, one expected")
        if [v for _, v in done] != list(range(1, iterations + 1)):
            problems.append(f"host {name}: {len(done)} complete visits, {iterations} expected in order")
    pairs = sum(len(done) for done in visits.values())
    try:
        final = _read_sequence(test_dir)
        if final != pairs:
            problems.append(f"the sequence file holds {final}, but {pairs} complete visits were logged")
    except GuardError as exc:
        problems.append(str(exc))
    lines = [f"host {name}: {len(done)} visits, {waits.get(name, 0)} waited over 0.1 s for the lock"
             for name, done in sorted(visits.items())]
    lines += [f"PROBLEM {p}" for p in problems]
    lines.append(f"LOCK TEST: PASS ({len(visits)} hosts, {pairs} visits, sequence 1 to {pairs})" if not problems
                 else f"LOCK TEST: FAIL ({len(problems)} problem(s))")
    return not problems, lines


# --- the CLI --------------------------------------------------------------------------------------------------

def format_usage(root: Path, usage: Usage) -> str:
    lines = [f"guard {root}", f"usage {usage.total} of {CEILING}"]
    for src in usage.sources:
        job = f" job {src['job_id']}" if src.get("job_id") else ""
        lines.append(f"  {src['source']} {src['id']}{job}: {src['count']} ({src['detail']})")
    if usage.uncertain:
        lines.append(f"uncertain: {' '.join(usage.uncertain)}; no slot is admitted while one holds slots. Check "
                     "squeue, then reconcile --resolve ID=released|job:JOBID --reason TEXT")
    if usage.unmanaged:
        lines.append(f"unmanaged: {' '.join(usage.unmanaged)}; counted, but submitted outside the guard: drain them or "
                     "account for them before licensed work")
    return "\n".join(lines)


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="python -m bbsim.licence_guard", description=__doc__.splitlines()[0])
    p.add_argument("--root", type=Path, default=None,
                   help=f"the shared guard directory (default {DEFAULT_ROOT}); an inherited {GUARD_ENV} naming "
                        "another directory refuses")
    sub = p.add_subparsers(dest="command", required=True)
    status = sub.add_parser("status", help="the reconciled usage (takes the lock)")
    status.add_argument("--lock-info", action="store_true", help="the lock's identity and owner, without the lock")
    submit = sub.add_parser("submit", help="submit one bridge script by path: submit --kind K --slots N --cap N -- "
                                           "<sbatch options> <script> [script args]")
    submit.add_argument("--kind", choices=KINDS, required=True)
    submit.add_argument("--slots", type=int, required=True)
    submit.add_argument("--cap", type=int, required=True)
    submit.add_argument("--bridge-repo", type=Path, default=None,
                        help="the bridge clone (default: ../bbridge beside this checkout)")
    throttle = sub.add_parser("throttle", help="change a guarded array's throttle")
    throttle.add_argument("parent")
    throttle.add_argument("throttle", type=int)
    throttle.add_argument("--cap", type=int, required=True)
    reconcile = sub.add_parser("reconcile", help="reconcile and print; or resolve one reservation")
    reconcile.add_argument("--resolve", metavar="RESERVATION=released|job:ID")
    reconcile.add_argument("--reason")
    recover = sub.add_parser("recover-lock", help="the operator's recovery of a crashed owner's lock")
    recover.add_argument("--expected-lock", required=True, metavar="TOKEN")
    recover.add_argument("--confirm-quiescent", required=True, metavar="REASON")
    worker = sub.add_parser("lock-test-worker", help="one worker of I1's two-node lock test")
    worker.add_argument("--test-dir", type=Path, required=True)
    worker.add_argument("--expect-commit-file", type=Path, required=True)
    worker.add_argument("--iterations", type=int, default=200)
    worker.add_argument("--hosts", type=int, default=2)
    worker.add_argument("--barrier-timeout-s", type=float, default=300.0)
    worker.add_argument("--hold-s", type=float, default=0.01)
    worker.add_argument("--pause-s", type=float, default=0.05)
    worker.add_argument("--lock-timeout-s", type=float, default=120.0)
    check = sub.add_parser("lock-test-check", help="judge a lock test's artifacts")
    check.add_argument("--test-dir", type=Path, required=True)
    check.add_argument("--hosts", type=int, default=2)
    check.add_argument("--iterations", type=int, default=200)
    return p


def _check_inherited_root(root: Path) -> None:
    inherited = os.environ.get(GUARD_ENV)
    if inherited and Path(inherited).resolve() != root.resolve():
        raise GuardRefused(f"the inherited {GUARD_ENV}={inherited} names another guard than {root}; one guard holds "
                           "the shared cap")


def _cmd_submit(args: argparse.Namespace, guard: Guard, sbatch_args: list[str]) -> int:
    options, script, script_args = split_user_sbatch_args(sbatch_args)
    check_user_sbatch_args(options)
    if script is None:
        raise GuardRefused("submit names one bridge script after the sbatch options")
    verified = check_path_mode_script(script, args.bridge_repo)
    if (args.kind, args.slots) != (verified.kind, verified.slots):
        raise GuardRefused(f"{verified.name} is submitted as --kind {verified.kind} --slots {verified.slots}, not "
                           f"--kind {args.kind} --slots {args.slots}")
    context = {"source": "submit by path", "script": str(verified.path), "script_rel": verified.rel,
               "script_sha256": verified.sha256, "bridge_repo": str(verified.bridge_repo),
               "bridge_commit": verified.bridge_commit, "sbatch_options": options, "script_args": script_args,
               "cwd": os.getcwd()}
    job = guard.submit(kind=args.kind, slots=args.slots, cap=args.cap,
                       argv=[*options, str(verified.path), *script_args], context=context)
    print(job)
    print(f"licence_guard: submitted {verified.name} as job {job}, reservation {guard.last_reservation_id} "
          f"({args.kind}, {args.slots} slot(s), cap {args.cap}); bridge {verified.bridge_commit}, script sha256 "
          f"{verified.sha256}", file=sys.stderr)
    return 0


def main(argv: list[str] | None = None, *, scheduler: Scheduler | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    sbatch_args = None
    if "--" in argv:
        cut = argv.index("--")
        argv, sbatch_args = argv[:cut], argv[cut + 1:]
    parser = _parser()
    try:
        args = parser.parse_args(argv)
    except SystemExit as exc:             # argparse has printed its message
        return 2 if exc.code else 0
    try:
        if (args.command == "submit") != (sbatch_args is not None):
            raise GuardRefused("submit, and only submit, takes '--' followed by the sbatch options and the script")
        if args.command == "lock-test-check":
            ok, lines = check_lock_test(args.test_dir, hosts=args.hosts, iterations=args.iterations)
            print("\n".join(lines))
            return 0 if ok else 1
        root = args.root if args.root is not None else DEFAULT_ROOT
        _check_inherited_root(root)
        if args.command == "lock-test-worker":
            return run_lock_test_worker(root, args.test_dir, expect_commit_file=args.expect_commit_file,
                                        iterations=args.iterations, hosts=args.hosts,
                                        barrier_timeout_s=args.barrier_timeout_s, hold_s=args.hold_s,
                                        pause_s=args.pause_s, lock_timeout_s=args.lock_timeout_s)
        guard = Guard(root, scheduler if scheduler is not None else SlurmScheduler())
        if args.command == "status":
            if args.lock_info:
                info = guard.lock_info()
                print(json.dumps(info, indent=2, sort_keys=True, default=str))
                print(f"LOCK HELD {info['token']}" if info["held"] else "LOCK FREE")
                return 0
            print(format_usage(guard.root, guard.usage()))
            return 0
        if args.command == "submit":
            return _cmd_submit(args, guard, sbatch_args)
        if args.command == "throttle":
            guard.change_throttle(args.parent, args.throttle, cap=args.cap)
            print(f"THROTTLE {args.parent} {args.throttle}")
            print(format_usage(guard.root, guard.usage()))
            return 0
        if args.command == "reconcile":
            if args.resolve is not None:
                reservation, sep, outcome = args.resolve.partition("=")
                if not sep:
                    raise GuardRefused("--resolve takes RESERVATION=released or RESERVATION=job:ID")
                if args.reason is None:
                    raise GuardRefused("--resolve needs --reason TEXT")
                guard.resolve(reservation, outcome, args.reason)
                print(f"RESOLVED {reservation} {outcome}")
            elif args.reason is not None:
                raise GuardRefused("--reason goes with --resolve")
            print(format_usage(guard.root, guard.usage()))
            return 0
        if args.command == "recover-lock":
            guard.recover_lock(args.expected_lock, confirm_quiescent=args.confirm_quiescent)
            print(f"LOCK RECOVERED {args.expected_lock}; the stale lock is kept under {guard.root / 'lock-evidence'}; "
                  "reservations and leases are unchanged: reconcile them separately")
            return 0
        raise GuardRefused(f"unknown command {args.command}")
    except GuardRefused as exc:
        print(f"licence_guard: refused: {exc}", file=sys.stderr)
        return 2
    except (GuardError, OSError) as exc:
        print(f"licence_guard: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
