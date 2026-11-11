"""One shared admission for HFSS solves on the cluster: at most 10 at a time (OR-045, HF-057, HF-058).

Every licensed submission (job arrays, single solves, re-exports) and every zero-slot job the job array submits
(monitor, assembly, the cluster probe) goes through ``Guard.submit``, under one lock in one shared directory, by
default ``DEFAULT_ROOT``. Arbitrary scripts, other users and hand-started AEDT sessions are outside it; only a
scheduler-enforced limit would cover those.

The lock is an atomic ``os.mkdir(<root>/lock)`` holding ``owner.json`` {nonce, host, boot_id, pid, process_start,
slurm_job, utc}; ``flock`` on a network filesystem is not relied on. A waiting client retries every ``LOCK_POLL_S``
seconds for ``lock_timeout_s`` (``LOCK_TIMEOUT_S`` for the operator's commands, ``IN_JOB_LOCK_TIMEOUT_S`` for a
job's launch lease) and then refuses with ``lock_info``'s observation and the recovery command. A lock is
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

Usage is reconciled under the lock from the scheduler's queue, the reservations and the leases:

- the queue holds pending, running, completing, suspended, configuring, stopped, requeued and staging-out jobs
  (``QUEUE_STATES``), read with ``SLURM_BITSTR_LEN=0`` so no array index string is cut short;
- array tasks are grouped by parent; an array whose comment starts ``bbsim:``, whose tasks are ``bbsim_hfss`` or
  ``bbsim_reexport`` jobs, or that a reservation is bound to, counts the larger of its ``ArrayTaskThrottle`` and its
  running tasks, or all its remaining tasks when the throttle is 0 or unknown; other arrays count 0; a single
  ``bbsim_hfss`` or ``bbsim_reexport`` job counts 1, and so does any queued job a slot-bearing reservation is bound
  to; held and pending tasks still count;
- a ``pending`` or ``uncertain`` reservation counts its slots, and a pending submission is uncertain (only a crash
  leaves one, since ``submit`` holds the lock from reservation to bind); a ``bound`` one counts, once per job id, the largest of
  its slots, the scheduler's count for that id and its held leases while the job is queued; once the job has left the
  queue it is released only when accounting (``sacct``) calls it terminal and no held lease names it, and counts its
  slots until then;
- a direct session's reservation (kind ``direct``, ``reserve_direct``) is released only when the leases that name it
  are all released, i.e. after its own lease verified cleanup;
- a throttle change sets the reservation to the largest of the requested throttle, the array's running tasks and
  its held leases, admits what that adds like a submission, and drains down to the request as those stop.

Admission: the reconciled total plus the new slots must not exceed ``min(CEILING, cap)``, ``cap`` an integer 1 to 10;
no slot is admitted while an uncertain reservation holds slots. Zero-slot kinds (``monitor``, ``assembly``) are
logged and never counted. The ceiling counts only our own solves, and the licence pool is shared with other users, so
a submission of one slot or more first reads the licence server (``LmutilServer``, ``lmutil lmstat -f`` per feature)
outside the lock, never under it: each of ``HEADROOM_FEATURES``, the five features one solve takes one of, must have
as many free (issued minus in use) as the slots requested, or it refuses before anything is reserved. Under the lock
the same reading must also hold our own counted slots that have not started (no launch lease yet, so not in lmstat's
in-use count; ``queued_ours``). A reading that cannot be read is retried every ``HEADROOM_RETRY_INTERVAL_S`` for up to
``HEADROOM_RETRY_WINDOW_S`` and then refuses, with no override; a full pool, or a feature the server issues none of,
refuses at once. The reading is advisory, recorded as the reservation context's ``headroom`` and printed as one
``licence headroom:`` line on stderr. A throttle raise is judged the same way against a reading taken before the lock;
a lowered throttle never reads the server (OR-045). A submission reserves before calling
``sbatch --parsable --export=ALL`` and binds the returned id before the lock is released; a non-zero exit, an id that
is not a plain number (``;cluster`` is dropped) or any interruption (``KeyboardInterrupt`` included) leaves the
reservation ``uncertain``, counted, until ``reconcile --resolve RESERVATION=released|job:<id> --reason TEXT`` after
checking ``squeue``; a bind that cannot be written after sbatch accepted the job names the job and that command.
``sbatch --wait`` is never run under the lock.

Launch leases (``<root>/leases/<id>.json``): every owned AEDT start takes one (``bbsim.session.aedt_session`` through
``launch_lease_from_env``). Admission reads only their ``reservation_id`` and ``state``: any state but ``released``
holds, unknown cleanup included (HF-058(2)), and a held lease naming no live reservation counts one slot of its own.
A job's lease consumes its reservation's slots (at most ``slots`` held per reservation); a session with no reservation
first takes a direct one-slot reservation; at most ``CEILING`` leases are held in total. A lease records {id,
reservation_id, direct, state (acquired, launching, retained, released), cleanup (unknown, live, stopped), outcome,
cleanup_waited_s, identity, scope, processes, cleanup_checks, cleanup_evidence, events}. Before ``Hfss()`` its launch
scope is recorded (host, boot, the launcher and its ancestors, the numeric ``SLURM_JOB_ID``, the step incarnation, the
launcher's cgroup and the scope root: inside a Slurm job the whole job's cgroup ``job_<SLURM_JOB_ID>``, outside one the
launcher's own; ``ProcScopeObserver``), after it the desktop's process and its descendants. After
``release_desktop`` the session polls ``cleanup_status`` every ``CLEANUP_POLL_S`` for up to ``CLEANUP_WAIT_S`` and
releases the lease only on ``stopped``: the scope root's subtree read in full on its host with no process started since
the launch left in it, so a solver that outlives its desktop, or runs in another step of the job, keeps it held. A
vanished pid, an empty process list, another host, an unreadable cgroup, terminal accounting in ``sacct`` and Slurm's
cgroup setting are never proof. Otherwise the
lease stays held until the operator, after ``scope-check LEASE`` on its host, runs the logged ``release LEASE
--confirm-cleanup REASON``. ``status --job ID`` lists a job's leases with their state, ``cleanup_waited_s`` and
outcome. A job checks, before AEDT, that it runs under the shared guard and its own reservation (``check_job_guard``,
CLI ``check-job``); a second guard named in ``GUARD_ENV`` refuses.

The guard owns the script, the job's context and identity, its array indices and throttle, its working directory and
the submission protocol. User sbatch arguments (``check_user_sbatch_args``) are complete canonical long names from
``ALLOWED_USER_OPTIONS`` only, as ``--opt=value`` or ``--opt value``, passed on byte for byte; abbreviations (sbatch's
getopt_long would take ``--pars`` as ``--parsable``), short options, unknown options, ``--`` and ``:`` are refused, and
so is every option in ``OWNED_SBATCH_OPTIONS``. Internal callers add only the complete names of
``INTERNAL_VALUE_OPTIONS`` (array, job name, comment, dependency and kill-on-invalid-dep, output, error, working
directory) as ``--name=value``, never the protocol (``--export*``, ``--parsable``, ``--quiet``, ``--wait``, ``--wrap``,
...); in stdin
mode no script pathname is allowed, and the script's own ``#SBATCH`` header may not set the protocol, the array or the
working directory, cluster short options or separate heterogeneous-job components either. sbatch's environment is the caller's without ``SBATCH_*`` and the scheduler-option aliases
(``SCHEDULER_ENV_ALIASES``), plus the caller's validated context and, set by the guard alone, ``GUARD_ENV`` (the
resolved root) and ``RESERVATION_ENV``. An inherited or explicit ``GUARD_ENV`` or a plan's ``guard_dir`` that resolves
to another directory refuses before anything is reserved, so a second guard cannot bypass the cap.

Path mode (the CLI's ``submit``, and ``Guard.submit`` without stdin; OR-044, HF-058(3)): one script, one of
``PATH_MODE_SCRIPTS``, whose real path lies in the bridge clone's ``scripts/hpc/hfss/``, tracked there and byte-equal
to its ``HEAD`` blob, then literal script arguments. The bridge clone is only the sibling ``../bbridge`` of the checkout
this module runs from: no caller can name another repository, whose own committed edits would pass every check. Its
kind and slots follow the script.

The two-node lock test (``lock-test-worker``, ``lock-test-check``): each worker refuses unless this checkout is clean
and at the recorded commit, waits at a barrier until the expected number of distinct hosts has announced itself, then
visits the production lock ``--iterations`` times; inside each visit it creates the test's ``active`` sentinel
directory (an existing one is an overlap), advances the shared sequence, lists ``records/`` (every earlier visit's
record must be there: admission discovers reservations by listing a directory, which an NFS client may cache) and
adds its own, and logs an enter and a leave record. The check judges serialization by the sentinel, the sequence and
the listings in the shared log, never by clocks, and needs at least ``LOCK_TEST_MIN_SWITCHES`` hand-overs of the lock
between hosts, so a PASS shows exclusion under contention.

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
# The operator's commands (submit, status, ...) wait LOCK_TIMEOUT_S for the guard lock and then refuse with the
# lock's observation. The lease calls a running job makes (acquire before AEDT, the launch scope and processes, the
# settle after release_desktop) wait IN_JOB_LOCK_TIMEOUT_S: a concurrent submission holds the lock for an squeue, one
# scontrol per counted array and one sacct per bound job that has left the queue, each up to SCHEDULER_TIMEOUT_S, and
# a lease that times out ends a task whose queue time is spent, or leaves its lease held (HF-057, OR-045).
LOCK_TIMEOUT_S = 60.0
IN_JOB_LOCK_TIMEOUT_S = 600.0

KINDS = ("array", "single", "reexport", "monitor", "assembly")
ZERO_SLOT_KINDS = ("monitor", "assembly")
DIRECT_KIND = "direct"
STATES = ("pending", "bound", "uncertain", "released")
SOLVE_JOB_NAMES = frozenset({"bbsim_hfss", "bbsim_reexport"})
MANAGED_COMMENT_PREFIX = "bbsim:"
LEASE_RELEASED = "released"
QUEUE_STATES = "PD,R,CG,S,CF,ST,RQ,SO"
PENDING_QUEUE_STATES = frozenset({"PENDING"})             # every other listed state is an active task
TERMINAL_STATES = frozenset({"BOOT_FAIL", "CANCELLED", "COMPLETED", "DEADLINE", "FAILED", "NODE_FAIL",
                             "OUT_OF_MEMORY", "PREEMPTED", "REVOKED", "TIMEOUT"})

# Path mode: script name -> (kind, slots). The single-solve, modes and re-export wrappers, and I1's and I2b's
# zero-slot scripts.
PATH_MODE_SCRIPTS = {
    "run_frequency.sbatch": ("single", 1),
    "run_modes.sbatch": ("single", 1),          # the fibre bore's ports-only mode run (OR-053)
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
# Internal callers (the job array's submission helper) may add these complete names to the user allowlist, each as
# --name=value; nothing else, and no abbreviation (sbatch's getopt_long would take --ar=0-7 as --array=0-7). The
# dependent assembly sets --kill-on-invalid-dep=yes beside its --dependency, so a dependency Slurm cannot satisfy
# ends it instead of leaving it pending; a user --sbatch-arg cannot.
INTERNAL_VALUE_OPTIONS = ("--array", "--job-name", "--comment", "--dependency", "--kill-on-invalid-dep", "--output",
                          "--error", "--chdir")
ALLOWED_INTERNAL_OPTIONS = frozenset(INTERNAL_VALUE_OPTIONS) | ALLOWED_USER_OPTIONS
# The protocol the guard sets itself, named in the refusal of an internal option.
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
LOCK_TEST_RECORDS = "records"
# A PASS needs the lock handed between hosts at least this often in log order; two back-to-back runs give one.
LOCK_TEST_MIN_SWITCHES = 5

_OWNER = "owner.json"
_LONG = re.compile(r"--([A-Za-z0-9][A-Za-z0-9-]*)(?:(=)(.*))?", re.S)
_CONTROL = re.compile(r"[\x00-\x1f\x7f]")
_ENV_NAME = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
_JOB_ID = re.compile(r"[0-9]+")
_INDEX_RANGE = re.compile(r"([0-9]+)(?:-([0-9]+)(?::([0-9]+))?)?")
_ARRAY_TASK = re.compile(r"([0-9]+)_(?:([0-9]+)|\[([^\]]*)\])")
_COMMIT = re.compile(r"[0-9a-f]{40}|[0-9a-f]{64}")
_RESERVATION_ID = re.compile(r"r[0-9]{8}T[0-9]{6}Z-[0-9a-f]{8}")
_LEASE_ID = re.compile(r"l[0-9]{8}T[0-9]{6}Z-[0-9a-f]{8}")

# Launch leases (HF-058(2)): the in-job wait for an owned AEDT scope to stop after release_desktop, and its poll.
CLEANUP_WAIT_S = 120.0
CLEANUP_POLL_S = 2.0
LEASE_ACQUIRED, LEASE_LAUNCHING, LEASE_RETAINED = "acquired", "launching", "retained"
LEASE_STATES = (LEASE_ACQUIRED, LEASE_LAUNCHING, LEASE_RETAINED, LEASE_RELEASED)
EXITED_STATES = frozenset({"Z", "X"})          # exited, not yet reaped: holds no licence
# cgroup v1 hierarchies in the order tried; Slurm's proctrack/cgroup tracks a job's processes in freezer.
CGROUP_V1_PREFERENCE = ("freezer", "pids", "memory", "cpuset", "cpuacct", "cpu")
AEDT_PROCESS = re.compile(r"ansysedt|AnsysEM|hf3d|hfss", re.IGNORECASE)

# The licence features one 4-core HFSS solve checks out, one of each (measured 2026-11-10; OR-045). The pool is
# shared with other users, so a licensed submission first checks each has room for its slots.
HEADROOM_FEATURES = ("electronics_desktop", "electronics3d_gui", "elec_solve_hfss", "elec_solve_level1",
                     "elec_solve_level2")
LMSTAT_TIMEOUT_S = 60.0
# A reading that fails is retried (OR-045): lmstat could not be read for about a minute on 2026-11-10. Attempts are
# HEADROOM_RETRY_INTERVAL_S apart, and none starts more than HEADROOM_RETRY_WINDOW_S after the first. A full pool is a
# reading, not a failure, and is never retried.
HEADROOM_RETRY_WINDOW_S = 90.0
HEADROOM_RETRY_INTERVAL_S = 15.0
_retry_clock = time.monotonic                   # the retry loop's clock and wait (tests replace them)
_retry_sleep = time.sleep
LDLINUX_ENV = "BBSIM_LDLINUX"
DEFAULT_LDLINUX = "/lib64/ld-linux-x86-64.so.2"
LOADER_EXITS = (126, 127)                       # bash's "cannot execute" and "not found": rerun through the loader
_USERS_OF = re.compile(r"Users of (\S+):\s+\(Total of ([0-9]+) licenses? issued;\s+Total of ([0-9]+) licenses? "
                       r"in use\)")
_LMUTIL_BANNER = re.compile(r"lmutil - Copyright")


class GuardError(RuntimeError):
    """The guard cannot admit, reconcile or record something; after a reservation exists, it stays recorded."""


class GuardRefused(GuardError):
    """Refused before anything was reserved or submitted (exit 2)."""


class LeaseRefused(GuardRefused):
    """``Guard.acquire_lease`` refused a launch lease for want of capacity: ``reason`` is ``"ceiling"`` (``CEILING``
    leases held in all) or ``"slots"`` (every slot of the job's reservation held by a lease, an earlier task's whose
    cleanup is unverified included; HF-058(2)). A lock, record or reservation fault is another ``GuardError``."""

    def __init__(self, message: str, reason: str) -> None:
        super().__init__(message)
        self.reason = reason


class LicenceQueryError(GuardError):
    """The licence server cannot be read for a feature: a variable or lmutil missing, a failed or timed-out lmstat, or
    no ``Users of`` line for it. ``Guard.submit`` refuses a licensed submission on it, with no override (OR-045)."""


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

    def cpus_per_task(self, parent_id: str) -> int | None: ...

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
    env = {k: v for k, v in os.environ.items()
           if not (k.startswith(("SBATCH_", "SQUEUE_", "SACCT_", "SCONTROL_")) or k in SCHEDULER_ENV_ALIASES)}
    # squeue cuts an array's index string at 64 characters unless SLURM_BITSTR_LEN is 0; a cut string cannot be
    # counted, and an uncountable queue refuses every admission.
    env["SLURM_BITSTR_LEN"] = "0"
    return env


class SlurmScheduler:
    """The scheduler through Slurm's commands, each with a timeout; any failure is a ``GuardError``.

    ``squeue`` lists ``QUEUE_STATES``: pending, running, completing, suspended, configuring (about to run), stopped
    (it keeps its processes and licences), requeued and staging out. ``sacct`` runs with ``-P`` so a state such as
    ``CANCELLED by 1234`` is not truncated to ``CANCELLED+``.
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
        out = self._checked(["squeue", "-h", "-u", self.user, "-t", QUEUE_STATES, "-o", "%i|%F|%j|%k|%T"])
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

    def cpus_per_task(self, parent_id: str) -> int | None:
        """The job's (or array's) CPUs/Task, None when ``scontrol show job`` prints none; an array's records that
        disagree are a ``GuardError``. Read by the planner, never by admission (HF-057)."""
        out = self._checked(["scontrol", "show", "job", parent_id])
        values = {int(v) for v in re.findall(r"\bCPUs/Task=([0-9]+)", out)}
        if len(values) > 1:
            raise GuardError(f"scontrol show job {parent_id} prints CPUs/Task {sorted(values)}, not one value")
        return values.pop() if values else None

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

    queue = array_throttle = cpus_per_task = terminal_state = submit = set_throttle = _refuse


# --- the licence server: headroom before a licensed submission (OR-045) ------------------------------------------

class LicenceServer(Protocol):
    def counts(self, features: Sequence[str]) -> dict[str, tuple[int, int]]:
        """(issued, in use) per feature; ``LicenceQueryError`` when the server cannot be read."""
        ...


def parse_lmstat(feature: str, text: str) -> tuple[int, int]:
    """(issued, in use) from ``lmutil lmstat -f <feature>``'s one line ``Users of <feature>:  (Total of N license(s)
    issued;  Total of M license(s) in use)``; none, several or one not in that form is a ``LicenceQueryError``."""
    lines = [line.strip() for line in text.splitlines() if line.strip().startswith(f"Users of {feature}:")]
    if not lines:
        raise LicenceQueryError(f"lmstat printed no 'Users of {feature}:' line ({_lmstat_reason(text)})")
    if len(lines) > 1:
        raise LicenceQueryError(f"lmstat printed {len(lines)} 'Users of {feature}:' lines, not one: {lines[:3]!r}")
    match = _USERS_OF.fullmatch(lines[0])
    if match is None or match.group(1) != feature:
        raise LicenceQueryError(f"lmstat's line for {feature} is not 'Total of N licenses issued;  Total of M "
                                f"licenses in use': {lines[0][:200]!r}")
    return int(match.group(2)), int(match.group(3))


def _lmstat_reason(text: str) -> str:
    """What lmutil's output says went wrong: its first line containing ``Error``, else its last line that is not
    lmutil's copyright banner."""
    shown = [line.strip() for line in text.splitlines() if line.strip()]
    error = next((line for line in shown if "Error" in line), None)
    if error is not None:
        return f"its error: {error[:200]!r}"
    rest = [line for line in shown if not _LMUTIL_BANNER.match(line)]
    if rest:
        return f"its last line: {rest[-1][:200]!r}"
    return "only lmutil's banner" if shown else "no output"


def _run_lmutil(argv: list[str], timeout_s: float) -> tuple[int, str]:
    """Run ``argv`` with stdout and stderr together, as bash would: an executable that cannot be run is exit 126, one
    not found (or whose ELF interpreter is missing) exit 127; a timeout is a ``LicenceQueryError``."""
    try:
        run = subprocess.run(argv, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             timeout=timeout_s)
    except subprocess.TimeoutExpired as exc:
        raise LicenceQueryError(f"lmutil timed out after {timeout_s:g} s") from exc
    except FileNotFoundError as exc:
        return 127, str(exc)
    except OSError as exc:                          # permission denied, exec format error
        return 126, str(exc)
    return run.returncode, run.stdout.decode(errors="replace")


class LmutilServer:
    """The licence server through ``$ANSYS_ROOT/licensingclient/linx64/lmutil lmstat -c $ANSYSLMD_LICENSE_FILE -f
    <feature>``, one query per feature with a ``timeout_s`` timeout, as ``monitor_licences.sbatch`` runs it: an exit of
    126 or 127 reruns the query through the loader ``${BBSIM_LDLINUX:-/lib64/ld-linux-x86-64.so.2}``. The variables
    are read at each ``counts``, so building one touches nothing; a missing variable or lmutil is reported once, before
    any feature is queried, since it is no feature's fault. ``environ`` and ``runner`` (``(argv, timeout_s) ->
    (returncode, output)``) replace ``os.environ`` and the subprocess (tests)."""

    def __init__(self, *, environ: Mapping[str, str] | None = None,
                 runner: Callable[[list[str], float], tuple[int, str]] | None = None,
                 timeout_s: float = LMSTAT_TIMEOUT_S):
        self.environ = environ
        self.runner = runner if runner is not None else _run_lmutil
        self.timeout_s = timeout_s

    def counts(self, features: Sequence[str]) -> dict[str, tuple[int, int]]:
        env = os.environ if self.environ is None else self.environ
        for name in ("ANSYS_ROOT", "ANSYSLMD_LICENSE_FILE"):
            if not env.get(name):
                raise LicenceQueryError(f"{name} is not set: bb_env.sh sets it")
        lmutil = Path(env["ANSYS_ROOT"]) / "licensingclient" / "linx64" / "lmutil"
        if not lmutil.is_file():
            raise LicenceQueryError(f"lmutil is not at {lmutil}")
        found: dict[str, tuple[int, int]] = {}
        for feature in features:
            try:
                found[feature] = self._query(env, lmutil, feature)
            except LicenceQueryError as exc:
                raise LicenceQueryError(f"{feature}: {exc}") from exc
        return found

    def _query(self, env: Mapping[str, str], lmutil: Path, feature: str) -> tuple[int, int]:
        query = [str(lmutil), "lmstat", "-c", env["ANSYSLMD_LICENSE_FILE"], "-f", feature]
        how = "lmutil"
        code, output = self.runner(query, self.timeout_s)
        if code in LOADER_EXITS:
            loader = env.get(LDLINUX_ENV) or DEFAULT_LDLINUX
            how += f" (direct exit {code}) through the loader {loader}"
            code, output = self.runner([loader, *query], self.timeout_s)
        if code != 0:
            raise LicenceQueryError(f"{how} exited {code}: {_lmstat_reason(output)}")
        return parse_lmstat(feature, output)


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


def _split_internal(argv: Sequence[str], *, stdin_mode: bool) -> tuple[list[_Option], int | None]:
    """Internal sbatch arguments: listed long options, then (by path) the script's index; in stdin mode, no pathname.

    Only the complete names of ``ALLOWED_INTERNAL_OPTIONS`` pass. ``INTERNAL_VALUE_OPTIONS`` take ``--name=value``
    with a value; the user allowlist's value options may also take a separate token, as the CLI passes them on.
    """
    options: list[_Option] = []
    argv = list(argv)
    i = 0
    while i < len(argv):
        token = _refuse_token(argv[i])
        if token.startswith("-") and token != "-":
            match = _LONG.fullmatch(token)
            if match is None:
                raise GuardRefused(f"internal sbatch argument {token!r}: internal options use canonical long names")
            name, has_value, value = match.group(1), match.group(2) is not None, match.group(3)
            owned = [p for p in _PROTOCOL_NAMES if p.startswith(name)]
            if owned:
                raise GuardRefused(f"internal sbatch option --{name} would change the submission protocol "
                                   f"(--{', --'.join(owned)}), which the guard sets itself")
            if f"--{name}" not in ALLOWED_INTERNAL_OPTIONS:
                like = sorted(n for n in ALLOWED_INTERNAL_OPTIONS if n.startswith(f"--{name}"))
                what = f"an abbreviation of {', '.join(like)}" if like else "not on the guard's internal list"
                raise GuardRefused(f"internal sbatch option --{name} is {what}; internal options are the complete "
                                   f"names {', '.join(sorted(ALLOWED_INTERNAL_OPTIONS))}")
            if has_value and value == "":
                raise GuardRefused(f"internal sbatch option {token!r} has an empty value")
            if f"--{name}" in INTERNAL_VALUE_OPTIONS and not has_value:
                raise GuardRefused(f"internal sbatch option --{name} takes --{name}=VALUE")
            if not has_value and f"--{name}" in USER_VALUE_OPTIONS:
                if i + 1 >= len(argv):
                    raise GuardRefused(f"sbatch option --{name} needs a value")
                value_token = _refuse_token(argv[i + 1])
                if value_token == "" or value_token.startswith("-"):
                    raise GuardRefused(f"sbatch option --{name} is followed by {value_token!r}, not a value")
                options.append(_Option(name, value_token))
                i += 2
                continue
            options.append(_Option(name, value))
            i += 1
            continue
        if stdin_mode:
            raise GuardRefused(f"{token!r} among the sbatch arguments: sbatch would read it as the script, but the "
                               "guard submits the captured script on stdin (internal options use --name=value)")
        return options, i
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
    """``#SBATCH`` header lines that would set the protocol, the array or the working directory.

    A header gives long options, or single short options with a separate value (``-p expansion``): a cluster or an
    attached value (``-sa0-3`` is ``-s`` and ``-a 0-3``) is refused, as are heterogeneous-job separators.
    """
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
            if token in ("hetjob", "packjob", ":", "--"):
                problems.append(f"line {number}: directive {token} separates heterogeneous-job components or ends "
                                "the options")
            elif token.startswith("--"):
                name = token[2:].split("=", 1)[0]
                if name and any(p.startswith(name) for p in _DIRECTIVE_NAMES):
                    problems.append(f"line {number}: directive {token} sets what the guard owns")
            elif token.startswith("-") and len(token) > 1:
                if len(token) > 2:
                    problems.append(f"line {number}: directive {token} is a short-option cluster or an attached "
                                    "value, which can hide an option the guard owns; give long options")
                elif token[1] in _DIRECTIVE_SHORT:
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
    data: bytes = field(repr=False, default=b"")

    def facts(self) -> dict:
        return {"path": str(self.path), "rel": self.rel, "sha256": self.sha256, "bridge_repo": str(self.bridge_repo),
                "bridge_commit": self.bridge_commit}


def default_bridge_repo() -> Path:
    """The bridge clone beside this checkout, ``../bbridge`` (OR-044): the only one path mode takes.

    No caller can name another: a repository passed in would carry its own committed scripts past every check.
    """
    return HFSS_CHECKOUT.parent / "bbridge"


def check_path_mode_script(script: str | os.PathLike) -> PathScript:
    """The script, if its real path is a listed script of the bridge clone's folder, tracked and equal to HEAD."""
    bridge = default_bridge_repo()
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
                      bridge_commit=head.stdout.decode().strip(), data=data)


# --- the launch scope on this host: /proc and the cgroup filesystem -------------------------------------------

class ProcScopeObserver:
    """Process and cgroup evidence on this host, read from ``/proc`` and the cgroup filesystem (standard library).

    ``launch_scope`` (before ``Hfss()``) records the launcher's identity and its ancestors', the boot-relative clock
    tick at which the launch began, and the launcher's cgroup: the hierarchy of ``/proc/self/cgroup`` whose
    ``cgroup.procs`` (with its descendant cgroups') lists the launcher, cgroup v1's ``freezer`` first, where Slurm's
    proctrack/cgroup tracks a job. It also records the scope root, the cgroup subtree the cleanup check lists: inside
    a Slurm job (``SLURM_JOB_ID``) the whole job's cgroup, the ancestor whose path component is ``job_<SLURM_JOB_ID>``
    as a whole component, in the same hierarchy, so a solver started in another step of the job (through ``srun``,
    say) is seen; outside Slurm, the launcher's own cgroup. Inside a job whose cgroup path has no such component the
    root is an error, never the step's subtree. ``scope`` then counts every process listed in the root's subtree
    (every descendant's ``cgroup.procs``: under cgroup v2 only leaves hold processes) that started at or after that
    tick, other than the launcher and its recorded ancestors: AEDT's processes inherit the launcher's cgroup, and a
    solver reparented away from its exiting desktop stays in it, which a walk of the desktop's process tree would miss.
    Processes of the job that started before the launch (its extern step, slurmstepd) are outside the scope; one
    adopted into the job after it (an ssh session) is inside and holds the lease until it ends (HF-058(2)). An exited
    process not yet reaped (state Z) holds nothing. A cgroup that cannot be read anywhere in the subtree, or that no
    longer lists the living launcher, is incomplete evidence, never proof of cleanup. The roots, host, uid and tick
    rate are injectable for tests.
    """

    def __init__(self, proc_root: str | os.PathLike = "/proc", cgroup_root: str | os.PathLike = "/sys/fs/cgroup", *,
                 host: str | None = None, clock_ticks: int | None = None, uid: int | None = None):
        self.proc_root = Path(proc_root)
        self.cgroup_root = Path(cgroup_root)
        self._host = host
        self._hz = clock_ticks
        self.uid = os.getuid() if uid is None else uid

    def host(self) -> str:
        return self._host or socket.gethostname()

    def _ticks_per_s(self) -> int:
        if self._hz is None:
            try:
                self._hz = int(os.sysconf("SC_CLK_TCK"))
            except (AttributeError, ValueError, OSError):
                self._hz = 100
        return self._hz

    def boot_id(self) -> str | None:
        try:
            return (self.proc_root / "sys" / "kernel" / "random" / "boot_id").read_text().strip() or None
        except OSError:
            return None

    def stat(self, pid: int | str) -> dict | None:
        """pid, comm, state, ppid and start identity from ``<proc>/<pid>/stat``; None if the process is gone."""
        path = self.proc_root / str(pid) / "stat"
        try:
            text = path.read_text()
        except (FileNotFoundError, ProcessLookupError, NotADirectoryError):
            return None
        head, sep, tail = text.rpartition(")")
        try:
            own, comm = head.split(" (", 1)
            fields = tail.split()
            start = int(fields[19])
            return {"pid": int(own), "comm": comm, "state": fields[0], "ppid": int(fields[1]), "starttime": start,
                    "start": f"linux-starttime:{start}"}
        except (IndexError, ValueError) as exc:
            raise OSError(f"{path} is not understood ({exc})") from exc

    def process(self, pid: int) -> dict | None:
        return self.stat(pid)

    def _pids(self) -> list[int]:
        return sorted(int(name) for name in os.listdir(self.proc_root) if name.isdigit())

    def descendants(self, pid: int) -> list[dict]:
        """The processes below ``pid`` in the process tree now (a reparented orphan is not among them)."""
        try:
            pids = self._pids()
        except OSError:
            return []
        stats, children = {}, {}
        for p in pids:
            try:
                st = self.stat(p)
            except OSError:
                continue
            if st is not None:
                stats[p] = st
                children.setdefault(st["ppid"], []).append(p)
        found, todo, seen = [], list(children.get(pid, [])), set()
        while todo:
            p = todo.pop(0)
            if p in seen:
                continue
            seen.add(p)
            found.append({k: stats[p][k] for k in ("pid", "start", "ppid", "comm")})
            todo.extend(children.get(p, []))
        return found

    def _cgroup_lines(self, pid: int | str = "self") -> list[tuple[str, str, str]]:
        lines = []
        for line in (self.proc_root / str(pid) / "cgroup").read_text().splitlines():
            hierarchy, sep, rest = line.partition(":")
            controllers, sep2, path = rest.partition(":")
            if sep and sep2:
                lines.append((hierarchy, controllers, path))
        return lines

    def _candidates(self, lines: list[tuple[str, str, str]]) -> list[tuple[str, str, Path]]:
        v1 = [(c, p) for h, c, p in lines if h != "0" and c and not c.startswith("name=")]
        found: list[tuple[str, str, Path]] = []
        for want in CGROUP_V1_PREFERENCE:
            for controllers, path in v1:
                if want in controllers.split(","):
                    found.append((f"v1 {controllers}", path, self.cgroup_root / controllers / path.lstrip("/")))
        for h, c, path in lines:
            if h == "0" and c == "":
                found.append(("v2", path, self.cgroup_root / path.lstrip("/")))
                found.append(("v2 unified", path, self.cgroup_root / "unified" / path.lstrip("/")))
        unique, seen = [], set()
        for kind, path, directory in found:
            if directory not in seen:
                seen.add(directory)
                unique.append((kind, path, directory))
        return unique

    @staticmethod
    def _listing(directory: Path) -> set[int]:
        """Every pid in the ``cgroup.procs`` of ``directory`` and of its descendant cgroups; OSError if one is
        unreadable."""
        def refuse(exc: OSError) -> None:
            raise exc

        pids: set[int] = set()
        read = 0
        for folder, _, _ in os.walk(directory, onerror=refuse):
            text = (Path(folder) / "cgroup.procs").read_text()
            read += 1
            try:
                pids.update(int(x) for x in text.split())
            except ValueError as exc:
                raise OSError(f"{folder}/cgroup.procs is not understood ({exc})") from exc
        if not read:
            raise OSError(f"{directory} holds no cgroup.procs")
        return pids

    def launch_scope(self, slurm_job_id: str | None = None) -> dict:
        """The launch scope of this process, before AEDT starts; ``cgroup.readable`` is False when it cannot be read.

        ``slurm_job_id`` is this process's ``SLURM_JOB_ID`` (None outside Slurm); ``cgroup.scope_root`` records the
        subtree the cleanup check will list, {kind ``job`` or ``launcher``, dir, error}.
        """
        job = slurm_job_id or None
        root: dict[str, Any] = {"kind": "job" if job is not None else "launcher", "dir": None, "error": None}
        cgroup: dict[str, Any] = {"readable": False, "kind": None, "path": None, "dir": None, "error": None,
                                  "lines": None, "scope_root": root}
        if job is not None and not _JOB_ID.fullmatch(str(job)):
            root["error"] = f"SLURM_JOB_ID={job!r} is not a numeric job id, so the job's cgroup cannot be found"
        scope: dict[str, Any] = {"observer": "proc", "host": self.host(), "boot_id": self.boot_id(),
                                 "launcher": {"pid": os.getpid(), "start": None}, "ancestors": [], "cgroup": cgroup,
                                 "begin_ticks": None, "clock_ticks_per_s": self._ticks_per_s()}
        try:
            me = self.stat("self")
        except OSError as exc:
            me, why = None, str(exc)
        else:
            why = f"{self.proc_root / 'self' / 'stat'} does not exist"
        if me is None:
            cgroup["error"] = f"this host's processes cannot be observed: {why}"
            return scope
        scope["launcher"] = {"pid": me["pid"], "start": me["start"]}
        chain, seen, ppid = [], {me["pid"]}, me["ppid"]
        while ppid > 0 and ppid not in seen and len(chain) < 64:
            seen.add(ppid)
            try:
                st = self.stat(ppid)
            except OSError:
                st = None
            if st is None:
                break
            chain.append({"pid": st["pid"], "start": st["start"]})
            ppid = st["ppid"]
        scope["ancestors"] = chain
        try:
            uptime_s = float((self.proc_root / "uptime").read_text().split()[0])
        except (OSError, ValueError, IndexError) as exc:
            cgroup["error"] = (f"{self.proc_root / 'uptime'} cannot be read ({exc}), so processes started after the "
                               "launch cannot be told from earlier ones")
            return scope
        # One tick of margin toward inclusion: a process started in the tick the launch began is in scope.
        scope["begin_ticks"] = int(uptime_s * self._ticks_per_s()) - 1
        try:
            lines = self._cgroup_lines("self")
        except OSError as exc:
            cgroup["error"] = f"{self.proc_root / 'self' / 'cgroup'} cannot be read: {exc}"
            return scope
        cgroup["lines"] = [":".join(line) for line in lines]
        candidates = self._candidates(lines)
        if job is not None:                     # a hierarchy whose path names the job first, in the usual order
            token = f"job_{job}"
            candidates = sorted(candidates, key=lambda c: token not in c[1].split("/"))
        tried = []
        for kind, path, directory in candidates:
            try:
                listed = self._listing(directory)
            except OSError as exc:
                tried.append(f"{directory}: {exc}")
                continue
            if me["pid"] in listed:
                cgroup.update(readable=True, kind=kind, path=path, dir=str(directory))
                if root["error"] is None:
                    if job is None:
                        root["dir"] = str(directory)
                    else:
                        job_dir = _job_cgroup_dir(path, directory, job)
                        if job_dir is None:
                            root["error"] = (f"SLURM_JOB_ID={job}, but the launcher's cgroup {path} has no path "
                                             f"component job_{job}, so the job's cgroup cannot be found")
                        else:
                            root["dir"] = str(job_dir)
                return scope
            tried.append(f"{directory} does not list the launcher {me['pid']}")
        cgroup["error"] = "no cgroup of the launcher lists it: " + ("; ".join(tried) or "no cgroup line was understood")
        return scope

    def scope(self, scope: Mapping[str, Any]) -> dict:
        """The processes of a recorded launch scope now: ``complete`` only when its cgroup was read in full."""
        cg = scope.get("cgroup") if isinstance(scope.get("cgroup"), Mapping) else {}
        root = cg.get("scope_root") if isinstance(cg.get("scope_root"), Mapping) else {}
        result: dict[str, Any] = {"complete": False, "dir": root.get("dir"), "root": root.get("kind"), "live": [],
                                  "exited": [], "reason": None}
        if not cg.get("readable") or not cg.get("dir"):
            result["reason"] = f"the launch recorded no readable cgroup ({cg.get('error') or 'none recorded'})"
            return result
        if root.get("error") or not root.get("dir"):
            result["reason"] = f"the launch recorded no scope root ({root.get('error') or 'none recorded'})"
            return result
        begin = scope.get("begin_ticks")
        if isinstance(begin, bool) or not isinstance(begin, int):
            result["reason"] = "the launch recorded no start tick, so its processes cannot be told from earlier ones"
            return result
        try:
            listed = self._listing(Path(root["dir"]))
        except OSError as exc:
            what = "the job's cgroup" if root.get("kind") == "job" else "the launcher's cgroup"
            result["reason"] = f"{what} {root['dir']} cannot be read in full now: {exc}"
            return result
        launcher = scope.get("launcher") if isinstance(scope.get("launcher"), Mapping) else {}
        own = {(launcher.get("pid"), launcher.get("start"))}
        own |= {(a.get("pid"), a.get("start")) for a in scope.get("ancestors") or [] if isinstance(a, Mapping)}
        pid = launcher.get("pid")
        try:
            me = self.stat(pid) if isinstance(pid, int) and not isinstance(pid, bool) else None
        except OSError:
            me = None
        if me is not None and me["start"] == launcher.get("start") and me["state"] not in EXITED_STATES \
                and pid not in listed:
            result["reason"] = (f"the launcher {pid} is alive but the scope root {root['dir']} does not list it, so "
                                "that listing is not its scope")
            return result
        for p in sorted(listed):
            try:
                st = self.stat(p)
            except OSError as exc:
                result["reason"] = f"process {p} of the cgroup cannot be read: {exc}"
                return result
            if st is None or (p, st["start"]) in own or st["starttime"] < begin:
                continue
            entry = {k: st[k] for k in ("pid", "start", "ppid", "comm", "state")}
            (result["exited"] if st["state"] in EXITED_STATES else result["live"]).append(entry)
        result["complete"] = True
        return result

    def job_processes(self, job_id: str) -> dict:
        """Every process on this host whose ``/proc/<pid>/cgroup`` has the path component ``job_<job_id>``."""
        token = f"job_{job_id}"
        result: dict[str, Any] = {"complete": False, "processes": [], "reason": None}
        try:
            pids = self._pids()
        except OSError as exc:
            result["reason"] = f"{self.proc_root} cannot be listed: {exc}"
            return result
        for pid in pids:
            try:
                lines = self._cgroup_lines(pid)
            except (FileNotFoundError, ProcessLookupError, NotADirectoryError):
                continue
            except OSError as exc:
                result["reason"] = f"{self.proc_root / str(pid) / 'cgroup'} cannot be read: {exc}"
                return result
            paths = [path for _, _, path in lines if token in path.split("/")]
            if not paths:
                continue
            try:
                st = self.stat(pid)
            except OSError:
                st = None
            result["processes"].append({"pid": pid, "comm": st["comm"] if st else None,
                                        "start": st["start"] if st else None, "state": st["state"] if st else None,
                                        "cgroup": paths[0]})
        result["complete"] = True
        return result

    def aedt_processes(self) -> list[dict]:
        """This user's AEDT processes on this host, each with the Slurm job its cgroup names (None outside one)."""
        found = []
        try:
            pids = self._pids()
        except OSError:
            return found
        for pid in pids:
            folder = self.proc_root / str(pid)
            try:
                if folder.stat().st_uid != self.uid:
                    continue
                command = (folder / "cmdline").read_bytes().replace(b"\0", b" ").decode(errors="replace").strip()
                st = self.stat(pid)
            except OSError:
                continue
            if st is None or not (AEDT_PROCESS.search(command) or AEDT_PROCESS.search(st["comm"])):
                continue
            try:
                lines = self._cgroup_lines(pid)
            except OSError:
                lines = []
            job = next((part[4:] for _, _, path in lines for part in path.split("/")
                        if part.startswith("job_") and part[4:].isdigit()), None)
            found.append({"pid": pid, "comm": st["comm"], "job": job, "command": command[:200]})
        return found


def _job_cgroup_dir(path: str, directory: Path, job_id: str) -> Path | None:
    """The ancestor of ``directory`` (the launcher's cgroup at cgroup ``path``) whose path component is
    ``job_<job_id>`` as a whole component (``job_771`` is not ``job_77``), in the same hierarchy; None if none."""
    token = f"job_{job_id}"
    parts = [part for part in str(path).split("/") if part]
    if token not in parts:
        return None
    job_dir = Path(directory)
    for _ in range(len(parts) - parts.index(token) - 1):
        job_dir = job_dir.parent
    return job_dir if job_dir.name == token else None


def default_scope_observer() -> ProcScopeObserver:
    """The production observer: this host's ``/proc`` and ``/sys/fs/cgroup``."""
    return ProcScopeObserver()


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


def _reservation_problem(rec: Any, error: str | None, stem: str) -> str | None:
    """Why a reservation record read from ``<stem>.json`` is unusable, or None."""
    if error is not None:
        return error
    ok = (isinstance(rec, dict) and rec.get("id") == stem and rec.get("state") in STATES
          and rec.get("kind") in (*KINDS, DIRECT_KIND)
          and isinstance(rec.get("slots"), int) and not isinstance(rec.get("slots"), bool)
          and rec["slots"] >= 0 and isinstance(rec.get("events"), list)
          and (rec["state"] != "bound" or (isinstance(rec.get("job_id"), str) and _JOB_ID.fullmatch(rec["job_id"]))))
    return None if ok else "fields"


def format_headroom(headroom: Mapping[str, Any], slots: int) -> str:
    """The admitted submission's one line: ``licence headroom: <feature> <free> free of <issued>, ... (each solve takes
    one of each; <slots> requested; OR-045)``, the features in ``HEADROOM_FEATURES``'s order."""
    parts = [f"{f} {headroom[f][0] - headroom[f][1]} free of {headroom[f][0]}" for f in HEADROOM_FEATURES]
    return f"licence headroom: {', '.join(parts)} (each solve takes one of each; {slots} requested; OR-045)"


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
    throttles: dict[str, int | None]        # array parent -> ArrayTaskThrottle, for the arrays asked
    queued: int                             # our counted slots that may still start and hold no lease yet (OR-045)


class Guard:
    """The shared admission in ``root`` (resolved once), with its scheduler and this host's launch-scope observer."""

    def __init__(self, root: str | os.PathLike, scheduler: Scheduler, *, lock_timeout_s: float = LOCK_TIMEOUT_S,
                 clock: Callable[[], float] = time.time, observer: Any = None, licences: LicenceServer | None = None):
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
        # The launch-scope observer of this host (ProcScopeObserver); tests pass one reading a fake /proc.
        self.observer = observer if observer is not None else default_scope_observer()
        # The licence server read before a licensed submission (OR-045); None is lmutil's, built at the first one,
        # so the usage, a direct session's slot and the lease calls inside jobs never touch lmutil.
        self.licences = licences
        self.lock_timeout_s = float(lock_timeout_s)
        self.clock = clock
        self.last_reservation_id: str | None = None
        self.last_path_script: PathScript | None = None
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
            problem = _reservation_problem(rec, error, path.stem)
            if problem is not None:
                raise GuardError(f"reservation record {path} is unreadable or malformed ({problem}); "
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
            # A counted array (bound to a slot-bearing reservation, commented bbsim:, or of bbsim_hfss or
            # bbsim_reexport tasks) counts its throttle, or its running tasks if more run (a throttle lowered outside
            # the guard stops no running task); an unthrottled one, or one whose throttle is unknown, counts every
            # remaining task. Pending and held tasks free nothing (OR-045). Arrays of other jobs count 0.
            if group.is_array:
                if (managed_slots or any(c.startswith(MANAGED_COMMENT_PREFIX) for c in group.comments)
                        or group.names & SOLVE_JOB_NAMES):
                    if group.key not in throttles:
                        throttles[group.key] = self.scheduler.array_throttle(group.key)
                    throttle = throttles[group.key]
                    return max(throttle, group.active) if throttle else group.tasks
                return 0
            return 1 if (group.names & SOLVE_JOB_NAMES or managed_slots) else 0

        total, sources, uncertain, unmanaged, blocking = 0, [], [], [], []
        # queued: the slots of ours that may still start a solve and hold no launch lease yet, so lmstat's in-use
        # count does not hold them (OR-045). A held lease is a started task; a job gone from the queue starts nothing;
        # an array queues no more than its remaining tasks; jobs submitted outside the guard and leases naming no live
        # reservation queue nothing here.
        queued = 0
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
            queued += max(0, rec["slots"] - held.get(rid, 0))
            if n:
                what = state if state == "uncertain" or rec["kind"] == DIRECT_KIND else \
                    "pending with no bind: a crashed submission, uncertain"
                sources.append({"source": "reservation", "id": rid, "kind": rec["kind"], "state": state,
                                "job_id": None, "count": n, "detail": f"{what}, held leases {held.get(rid, 0)}"})
            # submit holds the lock from reservation to bind, so only a crash leaves a submission pending
            if state == "uncertain" or rec["kind"] != DIRECT_KIND:
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
                queued += max(0, min(n, group.tasks) - leases_held)
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
        return _Reconciled(usage=usage, groups=groups, counts=counts, blocking=blocking, held=held,
                           throttles=throttles, queued=queued)

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

    def _licence_answer(self) -> tuple[dict[str, tuple[int, int]] | None, str]:
        """One attempt at the licence server: all of ``HEADROOM_FEATURES`` read afresh, as ``(counts, "")``, or
        ``(None, why)`` when it cannot be read (a ``LicenceQueryError``, any other failure, or an answer without two
        nonnegative integers for every feature)."""
        try:
            if self.licences is None:
                self.licences = LmutilServer()
            answer = self.licences.counts(HEADROOM_FEATURES)
        except Exception as exc:                # LicenceQueryError, or a server that failed otherwise
            return None, str(exc) if isinstance(exc, LicenceQueryError) else f"{type(exc).__name__}: {exc}"
        counts: dict[str, tuple[int, int]] = {}
        bad = []
        for feature in HEADROOM_FEATURES:
            value = answer.get(feature) if isinstance(answer, Mapping) else None
            if isinstance(value, (tuple, list)) and len(value) == 2 \
                    and all(isinstance(n, int) and not isinstance(n, bool) and n >= 0 for n in value):
                counts[feature] = (value[0], value[1])
            else:
                bad.append(f"{feature} {value!r}")
        if bad:
            return None, f"its answer holds no issued and in-use counts for {', '.join(bad)}"
        return counts, ""

    def _read_headroom(self, what: str, outcome: str = "nothing is submitted") -> dict:
        """The licence reading for a licensed submission or a throttle raise, taken outside the lock: an lmstat can
        take seconds, and jobs' lease calls wait on the lock (OR-045). An attempt that cannot be read is retried
        ``HEADROOM_RETRY_INTERVAL_S`` later, each attempt reading every feature afresh, and none starting more than
        ``HEADROOM_RETRY_WINDOW_S`` after the first; then it refuses, with no override. A feature the server issues
        none of refuses at once, since nothing will free up. Returns one attempt's complete reading as the reservation
        records it, ``{feature: [issued, in_use], ..., "utc": ...}``; ``_judge_headroom`` judges it."""
        first = _retry_clock()
        attempts = 0
        while True:
            attempts += 1
            counts, why = self._licence_answer()
            if counts is not None:
                break
            now = _retry_clock()
            if now + HEADROOM_RETRY_INTERVAL_S - first > HEADROOM_RETRY_WINDOW_S:
                raise GuardRefused(f"{what}: the licence server cannot be read ({why}) after {attempts} "
                                   f"attempt{'' if attempts == 1 else 's'} over {now - first:g} s, so {outcome}; "
                                   "there is no override (OR-045)")
            _retry_sleep(HEADROOM_RETRY_INTERVAL_S)
        utc = self._now()
        none = [f for f, (issued, _) in counts.items() if issued == 0]
        if none:
            raise GuardRefused(f"{what}: the licence server issues none of {', '.join(none)} (0 issued), so nothing "
                               f"will free up and waiting will not help; each solve takes one of each of "
                               f"{', '.join(HEADROOM_FEATURES)}, so {outcome} (OR-045)")
        return {**{f: [issued, in_use] for f, (issued, in_use) in counts.items()}, "utc": utc}

    @staticmethod
    def _judge_headroom(reading: Mapping[str, Any], slots: int, queued: int, what: str,
                        outcome: str = "nothing was reserved or submitted") -> None:
        """Each of ``HEADROOM_FEATURES`` must have ``slots`` free (issued minus in use) on top of ``queued``, our own
        slots that have not started and so are not in lmstat's in-use count. The pool is shared with other users,
        whom the ceiling does not count, and a solve that cannot check a feature out runs with no solution. The
        reading is advisory, since another user can take a licence at any moment (OR-045)."""
        short = []
        for f in HEADROOM_FEATURES:
            issued, in_use = reading[f]
            if issued - in_use < slots + queued:
                short.append(f"{f} {issued - in_use} free of {issued} ({in_use} in use)"
                             + (f" with {queued} of ours queued, for {slots} requested" if queued else ""))
        if short:
            ours = f" on top of {queued} of ours that have not started" if queued else ""
            raise GuardRefused(f"{what}: the licence pool, shared with other users, has too few free for {slots} "
                               f"solve(s){ours}: {'; '.join(short)}. Each solve takes one of each of "
                               f"{', '.join(HEADROOM_FEATURES)}; {outcome}, so retry once they are free (OR-045)")

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
        """Reserve, run ``sbatch --parsable --export=ALL <argv>`` and bind the returned job id, under one lock.

        By path (``stdin`` None) the script is verified as ``check_path_mode_script`` verifies it, whoever calls:
        ``kind`` and ``slots`` must be the script's, sbatch gets its resolved path, and the context records the
        verified facts as ``path_script`` (also ``last_path_script``). Any interruption of sbatch, an exception or
        not, leaves the reservation uncertain.
        """
        _check_kind_slots(kind, slots)
        cap = _check_cap(cap)
        context = _check_context(context)
        self._check_root_claims(env, context)
        extras = _check_env_extras(env)
        if stdin is not None and not isinstance(stdin, (bytes, bytearray)):
            raise GuardRefused("stdin must be the captured script's bytes")
        argv = list(argv)
        options, index = _split_internal(argv, stdin_mode=stdin is not None)
        _check_array(kind, slots, options)
        self.last_path_script = None
        if stdin is not None:
            data = bytes(stdin)
        else:
            verified = check_path_mode_script(argv[index])
            if (kind, slots) != (verified.kind, verified.slots):
                raise GuardRefused(f"{verified.name} is submitted as --kind {verified.kind} --slots {verified.slots}, "
                                   f"not --kind {kind} --slots {slots}")
            argv[index] = str(verified.path)
            data = verified.data
            context["path_script"] = verified.facts()
            self.last_path_script = verified
        problems = script_directive_problems(data)
        if problems:
            raise GuardRefused("the script's #SBATCH header sets what the guard owns: " + "; ".join(problems))
        what = f"a {kind} submission of {slots} slot(s)"
        headroom = None
        if slots >= 1:                          # the licence server, read before the lock (OR-045)
            headroom = self._read_headroom(what)
            self._judge_headroom(headroom, slots, 0, what)
            context["headroom"] = headroom
        full = ["--parsable", "--export=ALL", *argv]
        with self.locked():
            state = self._reconcile_refusing()
            self._admit(state, slots, cap, what)
            if headroom is not None:            # the same reading, with our own slots that have not started
                self._judge_headroom(headroom, slots, state.queued, what)
                headroom["queued_ours"] = state.queued
                # stderr: a caller's stdout is the job id, which pastes capture
                print(format_headroom(headroom, slots), file=sys.stderr)
            record = self._create(kind=kind, slots=slots, cap=cap, context=context, sbatch_argv=full)
            rid = record["id"]
            sub_env = {k: v for k, v in os.environ.items() if not _scheduler_override(k) and k != RESERVATION_ENV}
            sub_env.update(extras)
            sub_env[GUARD_ENV] = str(self.root)
            sub_env[RESERVATION_ENV] = rid
            try:
                returncode, stdout = self.scheduler.submit(list(full), data if stdin is not None else None, sub_env)
            except BaseException as exc:          # KeyboardInterrupt and SystemExit included
                self._mark_uncertain(record, event="submit_error", error=f"{type(exc).__name__}: {exc}")
                if isinstance(exc, Exception):
                    raise GuardError(self._uncertain_message(record, f"sbatch failed: {exc}")) from exc
                raise
            job = _parse_job_id(stdout) if returncode == 0 else None
            if job is None:
                self._mark_uncertain(record, event="submit_ambiguous", returncode=returncode, stdout=stdout[-2000:])
                raise GuardError(self._uncertain_message(record, f"sbatch exited {returncode} with stdout "
                                                                 f"{stdout.strip()[-200:]!r}"))
            record["job_id"] = job
            try:
                self._change(record, "bound", event="bound", job_id=job)
            except BaseException as exc:
                message = (f"sbatch accepted job {job}, but binding reservation {rid} to it could not be recorded "
                           f"({type(exc).__name__}: {exc}). The reservation stays counted and blocks admission until "
                           f"it is resolved: python -m bbsim.licence_guard --root {self.root} reconcile --resolve "
                           f"{rid}=job:{job} --reason 'sbatch accepted job {job}; the bind write failed' (nothing is "
                           "needed if status already lists it bound)")
                if isinstance(exc, Exception):
                    raise GuardError(message) from exc
                print(f"licence_guard: {message}", file=sys.stderr)
                raise
            return job

    def _mark_uncertain(self, record: dict, *, event: str, **details: Any) -> None:
        try:
            self._change(record, "uncertain", event=event, **details)
        except BaseException as exc:          # the record stays pending, which reconciliation treats as uncertain
            log.warning("reservation %s could not be marked uncertain (%s); it stays pending, which counts the "
                        "same", record["id"], exc)

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
        """Set a guarded array's throttle; the reservation never counts fewer than the tasks that may run.

        The reservation becomes the largest of the requested throttle, the array's running tasks and its held leases,
        with a drain target when that exceeds the request; whatever this adds to what was counted is admitted like a
        submission, and needs the licence headroom a submission needs, from a reading taken before the lock (OR-045).
        Only a request equal to the scheduler's throttle, with no drain pending, changes nothing. A lowered throttle
        is a safety action: it never reads the licence server and never waits on it.
        """
        cap = _check_cap(cap)
        if not isinstance(parent_id, str) or not _JOB_ID.fullmatch(parent_id):
            raise GuardRefused(f"the array's parent job id must be a plain number, got {parent_id!r}")
        if isinstance(throttle, bool) or not isinstance(throttle, int) or not 1 <= throttle <= CEILING:
            raise GuardRefused(f"the throttle must be an integer from 1 to {CEILING}, got {throttle!r}")
        what = f"setting array {parent_id}'s throttle to {throttle}"
        reading = self._read_headroom(what, "the throttle is unchanged") \
            if self._throttle_is_raise(parent_id, throttle) else None
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
            if parent_id in state.throttles:
                in_effect = state.throttles[parent_id]
            else:
                try:
                    in_effect = self.scheduler.array_throttle(parent_id)
                except Exception as exc:        # nothing has changed yet
                    raise GuardRefused(f"array {parent_id}'s throttle cannot be read: {exc}") from exc
            if in_effect == throttle and record.get("throttle_target") is None:
                return
            held = state.held.get(record["id"], 0)
            counted = state.counts.get(record["id"], record["slots"])
            new = max(throttle, group.active, held)
            added = max(0, new - counted)
            self._admit(state, added, cap, what)
            details: dict[str, Any] = {}
            if added:
                if reading is None:             # never an lmstat under the lock
                    raise GuardRefused(f"{what}: it adds {added} slot(s), which need a licence reading taken before "
                                       "the lock, and none was taken because the reservation changed meanwhile; the "
                                       "throttle is unchanged, so retry the command (OR-045)")
                self._judge_headroom(reading, added, state.queued, what, "the throttle is unchanged")
                details["headroom"] = {**reading, "queued_ours": state.queued}
            previous, previous_target = record["slots"], record.get("throttle_target")
            record["slots"] = new
            record["throttle_target"] = throttle if new > throttle else None
            self._change(record, "bound", event="throttle_change", previous=previous, previous_target=previous_target,
                         throttle=throttle, scheduler_throttle=in_effect, active_tasks=group.active,
                         held_leases=held, **details)
            try:
                self.scheduler.set_throttle(parent_id, throttle)
            except Exception as exc:
                raise GuardError(f"scontrol could not set array {parent_id}'s throttle to {throttle}: {exc}; the "
                                 f"reservation records {record['slots']} slot(s), the larger count; check scontrol "
                                 "show job before retrying") from exc

    def _throttle_is_raise(self, parent_id: str, throttle: int) -> bool:
        """Whether a throttle request may add slots, decided before the lock (no lmstat runs under it) from the array's
        reservation as recorded, read without the lock and changing nothing: a request above its recorded slots. A
        record that cannot be read or found is no raise here; the same lookup under the lock refuses it (OR-045)."""
        try:
            records = [r for r in self._read_reservations().values()
                       if r["state"] == "bound" and r["job_id"] == parent_id and r["kind"] == "array"]
        except (GuardError, OSError):
            return False
        return len(records) == 1 and throttle > records[0]["slots"]

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

    # launch leases
    def _lease_file(self, lease_id: Any) -> Path:
        if not isinstance(lease_id, str) or not _LEASE_ID.fullmatch(lease_id):
            raise GuardRefused(f"{lease_id!r} is not a lease id (l<UTC time>-<8 hex digits>, from status --job)")
        return self._leases / f"{lease_id}.json"

    def lease(self, lease_id: str) -> dict:
        """A lease record as stored (read without the lock: records are replaced atomically)."""
        path = self._lease_file(lease_id)
        rec, error = _read_json(path)
        if rec is None and error is None:
            raise GuardRefused(f"there is no lease {lease_id} in {self._leases}")
        if error is not None or not isinstance(rec, dict) or rec.get("id") != lease_id \
                or not isinstance(rec.get("state"), str) or not isinstance(rec.get("processes"), list) \
                or not isinstance(rec.get("events"), list):
            raise GuardError(f"lease record {path} is unreadable or malformed ({error or 'fields'})")
        return rec

    def _new_lease_id(self) -> str:
        while True:
            lid = f"l{_compact(self.clock())}-{secrets.token_hex(4)}"
            if not (self._leases / f"{lid}.json").exists():
                return lid

    def _save_lease(self, rec: dict) -> None:
        _write_json_atomic(self._lease_file(rec["id"]), rec)

    def _lease_event(self, rec: dict, event: str, **details: Any) -> None:
        rec["events"].append({"utc": self._now(), "event": event, **details})
        self._save_lease(rec)
        self._append_event({"utc": self._now(), "event": event, "lease": rec["id"],
                            "reservation": rec.get("reservation_id"), "lease_state": rec["state"], **details,
                            "by": _identity()})

    def acquire_lease(self, *, reservation_id: str | None, identity: dict) -> str:
        """A launch lease before an owned AEDT start; its id.

        With ``reservation_id`` (the job's ``RESERVATION_ENV``) the reservation must hold a solve slot, not be
        released, and, once bound, be bound to this Slurm job; its held leases stay within its slots. Without one, a
        direct one-slot reservation is admitted first (``reserve_direct``), released only with this lease. At most
        ``CEILING`` leases are held in total.
        """
        identity = _check_context(identity)
        if reservation_id is not None and (not isinstance(reservation_id, str)
                                           or not _RESERVATION_ID.fullmatch(reservation_id)):
            raise GuardRefused(f"{reservation_id!r} is not a reservation id")
        with self.locked():
            held = [lease for lease in self._read_leases() if lease["state"] != LEASE_RELEASED]
            if len(held) >= CEILING:
                raise LeaseRefused(f"{len(held)} launch leases are held, the ceiling of {CEILING} (OR-045): no AEDT "
                                   "start is admitted until one is released after its cleanup is verified", "ceiling")
            direct = reservation_id is None
            if direct:
                rid = self.reserve_direct(cap=CEILING, context={"source": "direct session", "identity": identity})
            else:
                rec = self._read_reservations().get(reservation_id)
                if rec is None:
                    raise GuardRefused(f"there is no reservation {reservation_id} in {self.root}")
                problem = _reservation_use_problem(rec, os.environ)
                if problem:
                    raise GuardRefused(problem)
                on_it = sum(1 for lease in held if lease["reservation_id"] == reservation_id)
                if on_it + 1 > rec["slots"]:
                    raise LeaseRefused(f"reservation {reservation_id} holds {rec['slots']} slot(s) and {on_it} held "
                                       "launch lease(s) name it already; a lease is released only after its cleanup "
                                       "is verified", "slots")
                rid = reservation_id
            lid = self._new_lease_id()
            record = {"id": lid, "reservation_id": rid, "direct": direct, "state": LEASE_ACQUIRED,
                      "cleanup": "unknown", "outcome": None, "cleanup_waited_s": None, "identity": identity,
                      "acquired_by": _owner_identity(), "created_utc": self._now(), "scope": None, "processes": [],
                      "cleanup_checks": [], "cleanup_evidence": None, "events": []}
            try:
                self._leases.mkdir(exist_ok=True)
                self._lease_event(record, "lease_acquired", direct=direct, slurm_job=os.environ.get("SLURM_JOB_ID"))
            except BaseException:
                if direct:                      # no lease names it and nothing was launched
                    res = self._read_reservations().get(rid)
                    if res is not None and res["state"] != "released":
                        self._change(res, "released", event="released_unused", reason="its lease could not be written")
                raise
            return lid

    def begin_launch_scope(self, lease_id: str, identity: dict) -> None:
        """Record the owned launch scope before ``Hfss()``: host, boot, launcher and ancestors, the numeric
        ``SLURM_JOB_ID`` (an array task's own), the step incarnation and the readable cgroup; state ``launching``,
        cleanup ``unknown``. Recorded once."""
        identity = _check_context(identity)
        self._lease_file(lease_id)
        try:
            observed = dict(self.observer.launch_scope(slurm_job_id=os.environ.get("SLURM_JOB_ID") or None))
        except Exception as exc:                # no evidence: recorded as such, and cleanup stays unknown
            observed = {"observer": "failed", "host": socket.gethostname(), "boot_id": _boot_id(),
                        "cgroup": {"readable": False, "error": f"the observer failed: {type(exc).__name__}: {exc}"}}
        scope = {**observed, **_slurm_identity(os.environ), "began_utc": self._now(), "identity": identity}
        with self.locked():
            rec = self.lease(lease_id)
            if rec["state"] != LEASE_ACQUIRED or rec.get("scope") is not None:
                raise GuardRefused(f"lease {lease_id} is {rec['state']}: its launch scope is recorded once, before "
                                   "Hfss() is called")
            rec["state"], rec["cleanup"], rec["scope"] = LEASE_LAUNCHING, "unknown", scope
            cg = scope.get("cgroup") if isinstance(scope.get("cgroup"), dict) else {}
            self._lease_event(rec, "lease_launch_scope", host=scope.get("host"), slurm_job_id=scope.get("slurm_job_id"),
                              step_incarnation=scope.get("step_incarnation"), cgroup=cg.get("path"),
                              cgroup_readable=bool(cg.get("readable")), cgroup_error=cg.get("error"),
                              scope_root=cg.get("scope_root"))

    def set_lease_processes(self, lease_id: str, processes: Sequence[dict]) -> None:
        """Append observed desktop and solver processes (pid and start identity) to the lease's scope; earlier
        evidence is never replaced, and an empty observation adds nothing."""
        if not isinstance(processes, (list, tuple)):
            raise GuardRefused(f"processes must be a list of records, got {type(processes).__name__}")
        clean = []
        for p in processes:
            pid = p.get("pid") if isinstance(p, dict) else None
            if isinstance(pid, bool) or not isinstance(pid, int) or pid <= 0:
                raise GuardRefused(f"a process record needs a positive integer pid, got {p!r}")
            start = p.get("start")
            if start is not None and not isinstance(start, str):
                raise GuardRefused(f"process {pid}'s start identity must be a string or null, got {start!r}")
            entry = {"pid": pid, "start": start, "role": str(p.get("role") or "descendant")}
            entry.update({k: p[k] for k in ("comm", "ppid") if k in p})
            clean.append(entry)
        with self.locked():
            rec = self.lease(lease_id)
            if rec["state"] == LEASE_RELEASED:
                raise GuardRefused(f"lease {lease_id} is released; no process is recorded on it")
            known = {(q.get("pid"), q.get("start")) for q in rec["processes"]}
            added = []
            for entry in clean:
                if (entry["pid"], entry["start"]) not in known:
                    known.add((entry["pid"], entry["start"]))
                    added.append({**entry, "observed_utc": self._now()})
            rec["processes"].extend(added)
            self._lease_event(rec, "lease_processes", observed=len(clean),
                              added=[{"pid": e["pid"], "role": e["role"]} for e in added])

    def note_lease(self, lease_id: str, event: str, **details: Any) -> None:
        """Log one event of a lease's launch (a constructor that raised, a desktop release that failed, ...)."""
        with self.locked():
            rec = self.lease(lease_id)
            self._lease_event(rec, event, **json.loads(json.dumps(details, default=str)))

    def _observe_cleanup(self, rec: dict) -> tuple[str, dict]:
        """``stopped``, ``live`` or ``unknown`` for a lease's owned launch scope, with the evidence used.

        ``stopped`` needs independent complete evidence: no launch scope was ever recorded (``Hfss()`` was never
        called), the host has rebooted since, or the scope root's subtree (the whole job's cgroup inside a Slurm job,
        the launcher's own outside one) was read in full on its host and holds none of the launch's processes. A
        recorded desktop or solver alive with its start identity, or any process of that subtree started since the
        launch, is ``live``. Anything else, an empty process list, a vanished or reused pid, another
        host, an unreadable cgroup included, is ``unknown``. Accounting and Slurm's settings are never read here.
        """
        try:
            here = self.observer.host()
        except Exception:                       # noqa: BLE001 - an observer that cannot name its host proves nothing
            here = None
        evidence: dict[str, Any] = {"utc": self._now(), "checked_on": here, "basis": None, "detail": None}
        if rec["state"] == LEASE_RELEASED:
            evidence.update(basis="released", detail="the lease is released")
            return "stopped", evidence
        scope = rec.get("scope")
        if not isinstance(scope, dict):
            if rec["state"] == LEASE_ACQUIRED and not rec.get("processes"):
                evidence.update(basis="never_launched", detail="no launch scope was recorded, so Hfss() was never "
                                                               "called under this lease")
                return "stopped", evidence
            evidence.update(basis="no_scope", detail="the lease records processes but no launch scope")
            return "unknown", evidence
        host = scope.get("host")
        evidence["lease_host"] = host
        if here is None or host != here:
            evidence.update(basis="other_host", detail=f"recorded on host {host!r}; checked on {here!r}, where its "
                                                        "processes cannot be seen: check on the lease's host")
            return "unknown", evidence
        recorded_boot = scope.get("boot_id")
        try:
            boot = self.observer.boot_id()
        except Exception:                       # noqa: BLE001
            boot = None
        if recorded_boot and boot and recorded_boot != boot:
            evidence.update(basis="rebooted", detail=f"host {host} has rebooted since the launch (boot "
                                                     f"{recorded_boot}, now {boot}), which ended every process of it")
            return "stopped", evidence
        recorded = []
        for p in rec.get("processes") or []:
            entry = {"pid": p.get("pid"), "start": p.get("start"), "role": p.get("role")}
            try:
                st = self.observer.process(p.get("pid"))
            except Exception as exc:            # noqa: BLE001
                entry["now"], entry["error"] = "unreadable", str(exc)
            else:
                if st is None:
                    entry["now"] = "gone"
                elif not p.get("start") or not st.get("start"):
                    entry["now"] = "unverifiable"
                elif st["start"] != p.get("start"):
                    entry["now"] = "reused"
                elif st.get("state") in EXITED_STATES:
                    entry["now"] = "exited"
                else:
                    entry["now"], entry["comm"] = "alive", st.get("comm")
            recorded.append(entry)
        evidence["recorded"] = recorded
        alive = [r for r in recorded if r["now"] == "alive"]
        if alive:
            evidence.update(basis="recorded_process_alive", live=alive,
                            detail=f"{len(alive)} recorded process(es) alive")
            return "live", evidence
        try:
            cg = dict(self.observer.scope(scope))
        except Exception as exc:                # noqa: BLE001
            cg = {"complete": False, "reason": f"the observer failed: {type(exc).__name__}: {exc}"}
        evidence["cgroup"] = {k: cg.get(k) for k in ("complete", "root", "dir", "reason", "exited")}
        recorded_root = (scope.get("cgroup") or {}).get("scope_root") if isinstance(scope.get("cgroup"), dict) else None
        evidence["scope_root"] = dict(recorded_root) if isinstance(recorded_root, dict) else None
        where = f"the {'job' if cg.get('root') == 'job' else 'launcher'}'s cgroup {cg.get('dir')}"
        if cg.get("complete"):
            if cg.get("live"):
                evidence.update(basis="cgroup", live=cg["live"],
                                detail=f"{len(cg['live'])} process(es) in {where} started since the launch")
                return "live", evidence
            evidence.update(basis="cgroup_empty", detail=f"{where} holds no process started since the launch, other "
                                                         "than the launcher and its ancestors")
            return "stopped", evidence
        seen = ", ".join(f"{r['pid']} {r['now']}" for r in recorded) or "none"
        evidence.update(basis="incomplete", detail=f"{cg.get('reason')}; recorded processes: {seen}; a vanished pid "
                                                   "alone is not proof of cleanup (HF-058(2))")
        return "unknown", evidence

    def _record_check(self, rec: dict, status: str, evidence: dict) -> None:
        rec["cleanup"] = status
        rec["cleanup_evidence"] = evidence
        checks = rec.setdefault("cleanup_checks", [])
        checks.append({"utc": evidence.get("utc"), "status": status, "basis": evidence.get("basis"),
                       "live": [p.get("pid") for p in evidence.get("live") or []]})
        del checks[:-200]
        self._save_lease(rec)

    def cleanup_status(self, lease_id: str) -> str:
        """Observe the lease's owned launch scope now (``_observe_cleanup``) and record what was seen."""
        with self.locked():
            rec = self.lease(lease_id)
            status, evidence = self._observe_cleanup(rec)
            if rec["state"] != LEASE_RELEASED:
                self._record_check(rec, status, evidence)
            return status

    def release_lease(self, lease_id: str, *, confirm_cleanup: str | None = None,
                      cleanup_waited_s: float | None = None) -> bool:
        """Release a lease when its scope is observed ``stopped``, or on an operator's ``confirm_cleanup`` reason
        attesting verified cleanup (refused while a process of it is seen alive here); otherwise keep it held and
        return False. The release event records who, why, the evidence, the wait and the outcome; a direct reservation
        is released with its cleaned lease. Accounting never releases a lease."""
        reason = None
        if confirm_cleanup is not None:
            reason = confirm_cleanup.strip() if isinstance(confirm_cleanup, str) else ""
            if not reason:
                raise GuardRefused("--confirm-cleanup needs a reason: what verified the cleanup (scope-check on the "
                                   "lease's host showing SCOPE EMPTY, and what else was checked)")
        if cleanup_waited_s is not None and (isinstance(cleanup_waited_s, bool)
                                             or not isinstance(cleanup_waited_s, (int, float))
                                             or not math.isfinite(cleanup_waited_s) or cleanup_waited_s < 0):
            raise GuardRefused(f"cleanup_waited_s must be a finite nonnegative number, got {cleanup_waited_s!r}")
        with self.locked():
            rec = self.lease(lease_id)
            if rec["state"] == LEASE_RELEASED:
                return True
            status, evidence = self._observe_cleanup(rec)
            self._record_check(rec, status, evidence)
            if reason is None and status != "stopped":
                return False
            if reason is not None and status == "live":
                raise GuardRefused(f"lease {lease_id}: its launch scope is seen live on this host "
                                   f"({describe_evidence(evidence)}); a confirmation cannot override processes seen "
                                   "alive")
            outcome = "verified_cleanup" if status == "stopped" else "operator_confirmed"
            rec["state"], rec["outcome"], rec["released_utc"] = LEASE_RELEASED, outcome, self._now()
            if cleanup_waited_s is not None:
                rec["cleanup_waited_s"] = float(cleanup_waited_s)
            rec["release"] = {"by": _identity(), "reason": reason, "cleanup_status": status, "evidence": evidence}
            self._lease_event(rec, "lease_released", outcome=outcome, reason=reason, cleanup_status=status,
                              cleanup_waited_s=rec.get("cleanup_waited_s"), evidence=evidence)
            if rec.get("direct"):
                res = self._read_reservations().get(rec["reservation_id"])
                if res is not None and res["kind"] == DIRECT_KIND and res["state"] != "released":
                    states = [lease["state"] for lease in self._read_leases()
                              if lease["reservation_id"] == res["id"]]
                    if states and all(s == LEASE_RELEASED for s in states):
                        self._change(res, "released", event="released_with_lease", lease=rec["id"], outcome=outcome)
            return True

    def retain_lease(self, lease_id: str, *, cleanup_waited_s: float) -> None:
        """Keep a lease held after the cleanup wait ended without verified cleanup, recording the wait, the outcome
        and the last evidence."""
        with self.locked():
            rec = self.lease(lease_id)
            if rec["state"] == LEASE_RELEASED:
                return
            rec["state"], rec["outcome"], rec["cleanup_waited_s"] = LEASE_RETAINED, "retained", float(cleanup_waited_s)
            self._lease_event(rec, "lease_retained", cleanup_status=rec.get("cleanup"),
                              cleanup_waited_s=rec["cleanup_waited_s"], evidence=rec.get("cleanup_evidence"))

    def leases_of_job(self, job_id: str) -> list[dict]:
        """The leases of a Slurm job or array: its id is the lease's recorded job or array id, or the job its
        reservation is bound to. Unreadable records are listed as such."""
        if not isinstance(job_id, str) or not _JOB_ID.fullmatch(job_id):
            raise GuardRefused(f"the job id must be a plain number, got {job_id!r}")
        try:
            bound = {rid for rid, r in self._read_reservations().items() if r.get("job_id") == job_id}
        except GuardError:
            bound = set()
        found = []
        if self._leases.is_dir():
            for path in sorted(self._leases.glob("*.json")):
                if path.name.startswith("."):
                    continue
                rec, error = _read_json(path)
                if error is not None or not isinstance(rec, dict):
                    found.append({"id": path.stem, "unreadable": error or "not an object"})
                    continue
                scope = rec.get("scope") if isinstance(rec.get("scope"), dict) else {}
                if job_id in (scope.get("slurm_job_id"), scope.get("slurm_array_job_id")) \
                        or rec.get("reservation_id") in bound:
                    found.append(rec)
        return sorted(found, key=lambda r: (str(r.get("created_utc", "")), r["id"]))


# --- launch leases: the guard a job or session uses, and the cleanup wait -------------------------------------

def _own_job_ids(env: Mapping[str, str]) -> list[str]:
    return [v for v in (env.get("SLURM_JOB_ID"), env.get("SLURM_ARRAY_JOB_ID")) if v and _JOB_ID.fullmatch(v)]


def _slurm_identity(env: Mapping[str, str]) -> dict:
    """The Slurm identity of this process: its numeric job id (an array task's own), the array, step and restart."""
    def numeric(name: str) -> str | None:
        value = env.get(name)
        return value if value and _JOB_ID.fullmatch(value) else None

    job, step = numeric("SLURM_JOB_ID"), env.get("SLURM_STEP_ID")
    restart = env.get("SLURM_RESTART_COUNT") or "0"
    return {"slurm_job_id": job, "slurm_array_job_id": numeric("SLURM_ARRAY_JOB_ID"),
            "slurm_array_task_id": numeric("SLURM_ARRAY_TASK_ID"), "slurm_step_id": step,
            "slurm_restart_count": restart, "slurm_node": env.get("SLURMD_NODENAME"),
            "step_incarnation": f"{job}.{step or 'batch'}.{restart}" if job else None}


def _reservation_use_problem(rec: dict, env: Mapping[str, str], kinds: Sequence[str] | None = None) -> str | None:
    """Why this process may not take a launch lease under ``rec``, or None."""
    rid = rec["id"]
    if rec["state"] == "released":
        return f"reservation {rid} is released: it holds no capacity"
    if rec["kind"] in ZERO_SLOT_KINDS or rec["slots"] < 1:
        return f"reservation {rid} ({rec['kind']}, {rec['slots']} slot(s)) holds no solve slot"
    if kinds is not None and rec["kind"] not in kinds:
        return f"reservation {rid} is a {rec['kind']} reservation; this job takes {' or '.join(kinds)}"
    if rec["state"] == "bound":
        own = _own_job_ids(env)
        if not own:
            return (f"reservation {rid} is bound to job {rec['job_id']}, but this process is not in a Slurm job (no "
                    "SLURM_JOB_ID): a reservation admits only the job it was submitted as")
        if rec["job_id"] not in own:
            return f"reservation {rid} is bound to job {rec['job_id']}, not to this job ({', '.join(own)})"
    return None


def check_job_guard(*, kinds: Sequence[str] | None = None, env: Mapping[str, str] | None = None,
                    root: str | os.PathLike | None = None) -> dict:
    """The guard a job or an owned AEDT session uses, verified before anything starts (OR-045, HF-057).

    ``GUARD_ENV`` is required, never defaulted here, and must resolve to the shared ``DEFAULT_ROOT`` (or ``root``):
    a second guard cannot bypass the shared cap. ``RESERVATION_ENV``, required when ``kinds`` is given, must name a
    reservation in that guard that is not released, holds a solve slot, is of one of ``kinds`` and, once bound, is
    bound to this Slurm job (``SLURM_JOB_ID``, or an array task's ``SLURM_ARRAY_JOB_ID``). Returns {root, reservation,
    kind, state, job_id, slurm_job_id}; raises ``GuardRefused``.
    """
    env = os.environ if env is None else env
    value = env.get(GUARD_ENV)
    if not value:
        raise GuardRefused(f"{GUARD_ENV} is not set: every owned AEDT start takes a launch lease from the shared guard "
                           f"(OR-045); a job gets it from python -m bbsim.licence_guard submit, and a session started "
                           f"by hand sets {GUARD_ENV}={DEFAULT_ROOT if root is None else root}")
    shared_root = DEFAULT_ROOT if root is None else Path(root)
    try:
        here = Path(value).resolve(strict=True)
    except (OSError, RuntimeError) as exc:
        raise GuardRefused(f"{GUARD_ENV}={value} cannot be reached ({exc})") from exc
    try:
        shared = Path(shared_root).resolve(strict=True)
    except (OSError, RuntimeError) as exc:
        raise GuardRefused(f"the shared guard {shared_root} cannot be reached ({exc})") from exc
    if here != shared:
        raise GuardRefused(f"{GUARD_ENV}={value} resolves to {here}, not the shared guard {shared}: a second guard "
                           "cannot bypass the cap")
    if not here.is_dir():
        raise GuardRefused(f"the guard {here} is not a directory")
    found = {"root": str(here), "reservation": None, "kind": None, "state": None, "job_id": None,
             "slurm_job_id": env.get("SLURM_JOB_ID") or None}
    rid = env.get(RESERVATION_ENV) or None
    if rid is None:
        if kinds is not None:
            raise GuardRefused(f"{RESERVATION_ENV} is not set: this job runs only under the reservation python -m "
                               "bbsim.licence_guard submit admitted it with")
        return found
    if not _RESERVATION_ID.fullmatch(rid):
        raise GuardRefused(f"{RESERVATION_ENV}={rid!r} is not a reservation id")
    path = here / "reservations" / f"{rid}.json"
    rec, error = _read_json(path)
    if rec is None and error is None:
        raise GuardRefused(f"there is no reservation {rid} in {here}")
    problem = _reservation_problem(rec, error, rid)
    if problem is not None:
        raise GuardRefused(f"reservation record {path} is unreadable or malformed ({problem})")
    problem = _reservation_use_problem(rec, env, kinds)
    if problem is not None:
        raise GuardRefused(problem)
    found.update(reservation=rid, kind=rec["kind"], state=rec["state"], job_id=rec.get("job_id"))
    return found


class LaunchLease(str):
    """A launch lease's id (a ``str``) with the ``Guard`` that holds it."""

    guard: Guard

    def __new__(cls, lease_id: str, guard: Guard) -> LaunchLease:
        obj = super().__new__(cls, lease_id)
        obj.guard = guard
        return obj


def _check_wait(value: Any, what: str = "cleanup_wait_s") -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
        raise GuardRefused(f"{what} must be a finite nonnegative number of seconds, got {value!r}")
    return float(value)


@contextmanager
def launch_lease_from_env(identity: dict, *, cleanup_wait_s: float = CLEANUP_WAIT_S, observer: Any = None,
                          scheduler: Scheduler | None = None) -> Iterator[LaunchLease]:
    """A launch lease from the job's guard (``GUARD_ENV``, required; ``RESERVATION_ENV``, optional): ``GuardError``
    before anything when the guard is unset, unreachable or not the shared one, or the reservation is not this job's.

    On exit, whatever happened inside, ``settle_lease`` polls the lease's scope every ``CLEANUP_POLL_S`` for up to
    ``cleanup_wait_s`` and releases it on ``stopped``, or keeps it held and logs what it saw (HF-058(2)). A failure to
    settle leaves the lease held and is logged, never raised over the body's own exception. Every lock this lease's
    guard takes waits up to ``IN_JOB_LOCK_TIMEOUT_S``, since a concurrent submission can hold it longer than the
    operator's ``LOCK_TIMEOUT_S``.
    """
    wait = _check_wait(cleanup_wait_s)
    found = check_job_guard()
    guard = Guard(found["root"], scheduler if scheduler is not None else SlurmScheduler(),
                  lock_timeout_s=IN_JOB_LOCK_TIMEOUT_S, observer=observer)
    lease = LaunchLease(guard.acquire_lease(reservation_id=found["reservation"], identity=identity), guard)
    try:
        yield lease
    finally:
        try:
            settle_lease(guard, lease, wait_s=wait)
        except Exception as exc:                # noqa: BLE001 - the lease stays held; never mask the body's error
            log.error("launch lease %s could not be settled (%s: %s); it stays held until python -m "
                      "bbsim.licence_guard --root %s release %s --confirm-cleanup REASON after scope-check on its host",
                      lease, type(exc).__name__, exc, guard.root, lease)


def describe_evidence(evidence: Mapping[str, Any] | None) -> str:
    """One line of a cleanup observation: its basis, the processes seen alive and the detail."""
    if not isinstance(evidence, Mapping):
        return "no observation recorded"
    live = []
    for p in evidence.get("live") or []:
        if isinstance(p, Mapping):
            live.append(f"pid {p.get('pid')} ({p.get('comm') or p.get('role') or '?'}, start {p.get('start')})")
    text = f"basis {evidence.get('basis')}"
    root = evidence.get("scope_root")
    if isinstance(root, Mapping):
        text += f"; scope root {root.get('kind')} {root.get('dir') or root.get('error')}"
    if live:
        text += f"; live: {', '.join(live)}"
    if evidence.get("detail"):
        text += f"; {evidence['detail']}"
    return text


def settle_lease(guard: Guard, lease_id: str, *, wait_s: float, poll_s: float | None = None) -> bool:
    """After ``release_desktop``: poll ``cleanup_status`` every ``poll_s`` (``CLEANUP_POLL_S``) for up to ``wait_s``
    and release the lease on ``stopped``, recording ``cleanup_waited_s``; otherwise keep it held, record the wait and
    log the scope seen (R-j, HF-058(2)). True if released."""
    wait = _check_wait(wait_s, "wait_s")
    poll = CLEANUP_POLL_S if poll_s is None else _check_wait(poll_s, "poll_s")
    start = time.monotonic()
    status = "unknown"
    while True:
        status = guard.cleanup_status(lease_id)
        if status == "stopped" and guard.release_lease(lease_id,
                                                       cleanup_waited_s=round(time.monotonic() - start, 3)):
            log.info("launch lease %s released: its launch scope stopped", lease_id)
            return True
        elapsed = time.monotonic() - start
        if elapsed >= wait:
            break
        time.sleep(max(0.0, min(poll, wait - elapsed)))
    waited = round(time.monotonic() - start, 3)
    guard.retain_lease(lease_id, cleanup_waited_s=waited)
    rec = guard.lease(lease_id)
    host = (rec.get("scope") or {}).get("host") if isinstance(rec.get("scope"), dict) else None
    log.warning("launch lease %s retained after %.1f s: cleanup %s (%s). It holds its solve slot until its scope is "
                "verified empty: python -m bbsim.licence_guard --root %s scope-check %s on host %s, then release %s "
                "--confirm-cleanup REASON", lease_id, waited, status, describe_evidence(rec.get("cleanup_evidence")),
                guard.root, lease_id, host, lease_id)
    return False


def scope_check(guard: Guard, lease_id: str) -> tuple[bool, list[str]]:
    """The check behind ``scope-check LEASE``, run on the lease's host: empty only when no recorded process of the
    lease (pid and start identity) is alive there and no process on the host is in a cgroup of the lease's job
    (``job_<SLURM_JOB_ID>`` as a whole path component). It also lists this user's AEDT processes with their jobs."""
    rec = guard.lease(lease_id)
    obs = guard.observer
    scope = rec.get("scope") if isinstance(rec.get("scope"), dict) else {}
    here = obs.host()
    job = scope.get("slurm_job_id")
    lines = [f"lease {lease_id}: state {rec.get('state')}, cleanup {rec.get('cleanup')}, recorded on host "
             f"{scope.get('host')!r}, slurm job {job or 'none'}; checked on {here!r}"]
    blocking = []
    if not scope:
        blocking.append("the lease recorded no launch scope, so neither its host nor its job is known")
    elif scope.get("host") != here:
        blocking.append(f"the lease was recorded on host {scope.get('host')!r}; run scope-check there (this is "
                        f"{here!r})")
    else:
        for p in rec.get("processes") or []:
            pid, role = p.get("pid"), p.get("role")
            try:
                st = obs.process(pid)
            except Exception as exc:            # noqa: BLE001
                blocking.append(f"recorded process {pid} ({role}) cannot be read: {exc}")
                continue
            if st is None:
                lines.append(f"recorded process {pid} ({role}): gone")
            elif not p.get("start"):
                blocking.append(f"recorded process {pid} ({role}) has no recorded start identity, and pid {pid} "
                                f"exists ({st.get('comm')}, start {st.get('start')})")
            elif st.get("start") != p.get("start"):
                lines.append(f"recorded process {pid} ({role}): gone (pid {pid} is now {st.get('comm')}, start "
                             f"{st.get('start')})")
            elif st.get("state") in EXITED_STATES:
                lines.append(f"recorded process {pid} ({role}): exited, not yet reaped")
            else:
                blocking.append(f"recorded process {pid} ({role}, {st.get('comm')}) is alive, start {st.get('start')}")
        if not job:
            blocking.append("the lease recorded no numeric SLURM_JOB_ID, so its job's cgroups cannot be checked")
        else:
            found = obs.job_processes(job)
            if not found.get("complete"):
                blocking.append(f"the processes in cgroups of job_{job} cannot be listed: {found.get('reason')}")
            for p in found.get("processes") or []:
                if p.get("state") in EXITED_STATES:
                    lines.append(f"process {p['pid']} ({p.get('comm')}) in {p.get('cgroup')}: exited, not yet reaped")
                else:
                    blocking.append(f"process {p['pid']} ({p.get('comm')}) is in {p.get('cgroup')}, a cgroup of "
                                    f"job_{job}")
    for p in obs.aedt_processes():
        lines.append(f"aedt process {p['pid']} ({p['comm']}) slurm job {p['job'] or 'none'}: {p.get('command')}")
    lines += [f"NOT EMPTY: {b}" for b in blocking]
    lines.append("SCOPE NOT EMPTY" if blocking else "SCOPE EMPTY")
    return not blocking, lines


def format_job_leases(root: Path, job_id: str, leases: Sequence[dict]) -> str:
    lines = [f"guard {root}", f"leases of job {job_id}:"]
    held, waits = 0, []
    for rec in leases:
        if "unreadable" in rec:
            lines.append(f"  PROBLEM lease record {rec['id']} is unreadable ({rec['unreadable']})")
            held += 1
            continue
        scope = rec.get("scope") if isinstance(rec.get("scope"), dict) else {}
        array = (f" (array {scope['slurm_array_job_id']}_{scope.get('slurm_array_task_id')})"
                 if scope.get("slurm_array_job_id") else "")
        wait = rec.get("cleanup_waited_s")
        if isinstance(wait, (int, float)) and not isinstance(wait, bool):
            waits.append(wait)
        cg = scope.get("cgroup") if isinstance(scope.get("cgroup"), dict) else {}
        root = cg.get("scope_root") if isinstance(cg.get("scope_root"), dict) else {}
        root_text = (f"{root.get('kind')}" + (f" (error: {root['error']})" if root.get("error") else "")
                     if root else "none")
        if rec.get("state") != LEASE_RELEASED:
            held += 1
        lines.append(f"  lease {rec['id']} job {scope.get('slurm_job_id') or 'none'}{array} host "
                     f"{scope.get('host') or 'none'} reservation {rec.get('reservation_id')}: state "
                     f"{rec.get('state')}, cleanup {rec.get('cleanup')}, cleanup_waited_s "
                     f"{wait if wait is not None else 'none'}, outcome {rec.get('outcome') or 'none'}, scope_root "
                     f"{root_text}")
    largest = max(waits) if waits else "none"
    lines.append(f"job {job_id}: {len(leases)} lease(s), {held} held, largest cleanup_waited_s {largest}")
    if held:
        lines.append("a held lease keeps its solve slot until its scope is verified empty: scope-check LEASE on its "
                     "host, then release LEASE --confirm-cleanup REASON")
    return "\n".join(lines)


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
                         clock: Callable[[], float] = time.time,
                         listdir: Callable[[Path], Sequence[str]] = os.listdir) -> int:
    """One worker of the two-node lock test on the production lock in ``root``; artifacts in ``test_dir``.

    Each visit, under the lock, also lists ``records/`` and requires every earlier visit's record there before it
    adds its own: admission finds reservations by listing a directory, and an NFS client can serve a cached listing
    that misses another node's new entry. Exit 2 when it refuses (the checkout, stale artifacts, the barrier's
    timeout), 1 on an overlap, a lock error or a record another visit could not see, 0 after ``iterations`` clean
    visits. ``host``, ``clock`` and ``listdir`` are injectable for local tests only.
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
        for sub in (LOCK_TEST_READY, LOCK_TEST_FAILURES, LOCK_TEST_RECORDS):
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
        if any((test_dir / LOCK_TEST_RECORDS).iterdir()):
            stale.append(LOCK_TEST_RECORDS)
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
    unseen = 0
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
                earlier = {str(i) for i in range(1, seq)}
                names = set(listdir(test_dir / LOCK_TEST_RECORDS))
                if names != earlier:
                    unseen += 1
                    fail("discovery", f"visit {visit} (seq {seq}): records/ listed {len(names & earlier)} of the "
                                      f"{seq - 1} earlier visits' records", 1, visit=visit, seq=seq,
                         missing=sorted(earlier - names, key=int)[:20], extra=sorted(names - earlier)[:20])
                try:
                    os.mkdir(test_dir / LOCK_TEST_RECORDS / str(seq))
                except FileExistsError:
                    unseen += 1
                    fail("discovery", f"visit {visit}: records/{seq} exists already", 1, visit=visit, seq=seq)
                record = {"host": host, "pid": pid, "visit": visit, "seq": seq}
                _append_line(log_path, {"event": "enter", **record, "utc": _utc(clock()),
                                        "lock_wait_s": round(waited, 4), "listed": len(names & earlier)})
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
    if unseen:
        print(f"lock-test-worker: {unseen} visit(s) on {host} missed records written before them", file=sys.stderr)
        return 1
    print(f"lock-test-worker: {iterations} visits on {host} (pid {pid})")
    return 0


def check_lock_test(test_dir: str | os.PathLike, *, hosts: int = 2, iterations: int = 200,
                    min_switches: int = LOCK_TEST_MIN_SWITCHES) -> tuple[bool, list[str]]:
    """Judge a lock test's artifacts, never by clocks.

    It requires: the expected distinct hosts with ``iterations`` complete visits each, enter and leave in pairs in
    log order, a sequence without gaps or repeats, every visit listing all earlier visits' records, ``records/``
    holding exactly one record per visit, no failure record, no sentinel at the end, and, with two or more hosts,
    the lock handed from one host to another at least ``min_switches`` times in log order.
    """
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
    order: list[str] = []

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
            if isinstance(seq, int) and not isinstance(seq, bool) and rec.get("listed") != seq - 1:
                problems.append(f"log line {number}: visit {who(rec)} listed {rec.get('listed')!r} of the {seq - 1} "
                                "earlier visits' records")
            if isinstance(rec.get("lock_wait_s"), (int, float)) and rec["lock_wait_s"] > 0.1:
                waits[str(rec.get("host"))] = waits.get(str(rec.get("host")), 0) + 1
        elif event == "leave":
            keys = ("host", "pid", "visit", "seq")
            if open_enter is None or any(rec.get(k) != open_enter[1].get(k) for k in keys):
                problems.append(f"log line {number}: leave of {who(rec)} without its own enter just before")
            else:
                visits.setdefault(str(rec.get("host")), []).append((rec.get("pid"), rec.get("visit")))
                order.append(str(rec.get("host")))
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
    records_dir = test_dir / LOCK_TEST_RECORDS
    found = sorted(p.name for p in records_dir.iterdir()) if records_dir.is_dir() else None
    if found is None:
        problems.append(f"no records directory {records_dir}")
    elif sorted(found) != sorted(str(i) for i in range(1, pairs + 1)):
        problems.append(f"records/ holds {len(found)} entries, not exactly one per visit 1 to {pairs}")
    switches = sum(1 for a, b in zip(order, order[1:]) if a != b)
    if hosts > 1 and switches < min_switches:
        problems.append(f"the lock passed between hosts {switches} time(s) in log order, fewer than {min_switches}: "
                        "the hosts did not alternate, so the run shows no exclusion under contention")
    lines = [f"host {name}: {len(done)} visits, {waits.get(name, 0)} waited over 0.1 s for the lock"
             for name, done in sorted(visits.items())]
    lines.append(f"hand-overs between hosts: {switches}")
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
    status = sub.add_parser("status", help="the reconciled usage (takes the lock); --lock-info and --job take none")
    which = status.add_mutually_exclusive_group()
    which.add_argument("--lock-info", action="store_true", help="the lock's identity and owner, without the lock")
    which.add_argument("--job", metavar="JOBID",
                       help="the launch leases of a Slurm job, or of an array by its parent id: state, "
                            "cleanup_waited_s and outcome of each")
    submit = sub.add_parser("submit", help="submit one script of the bridge clone ../bbridge by path: submit "
                                           "--kind K --slots N --cap N -- <sbatch options> <script> [script args]")
    submit.add_argument("--kind", choices=KINDS, required=True)
    submit.add_argument("--slots", type=int, required=True)
    submit.add_argument("--cap", type=int, required=True)
    throttle = sub.add_parser("throttle", help="change a guarded array's throttle")
    throttle.add_argument("parent")
    throttle.add_argument("throttle", type=int)
    throttle.add_argument("--cap", type=int, required=True)
    reconcile = sub.add_parser("reconcile", help="reconcile and print; or resolve one reservation")
    reconcile.add_argument("--resolve", metavar="RESERVATION=released|job:ID")
    reconcile.add_argument("--reason")
    release = sub.add_parser("release", help="the operator's logged release of a launch lease whose cleanup was "
                                             "verified (scope-check on its host first)")
    release.add_argument("lease")
    release.add_argument("--confirm-cleanup", required=True, metavar="REASON",
                         help="what verified the cleanup; logged with who ran it and what the guard saw")
    scope = sub.add_parser("scope-check", help="on the lease's host: SCOPE EMPTY only when none of its recorded "
                                               "processes is alive and no process is in a cgroup of its job")
    scope.add_argument("lease")
    job_check = sub.add_parser("check-job", help="a job's check, before AEDT, that it runs under the shared guard "
                                                 "and its own reservation (exit 2 refused)")
    job_check.add_argument("--kind", choices=("single", "reexport", "array"), required=True)
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
    lock_check = sub.add_parser("lock-test-check", help="judge a lock test's artifacts")
    lock_check.add_argument("--test-dir", type=Path, required=True)
    lock_check.add_argument("--hosts", type=int, default=2)
    lock_check.add_argument("--iterations", type=int, default=200)
    lock_check.add_argument("--min-switches", type=int, default=LOCK_TEST_MIN_SWITCHES,
                            help="hand-overs of the lock between hosts a PASS needs (default %(default)s)")
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
    context = {"source": "submit by path", "sbatch_options": options, "script_args": script_args,
               "cwd": os.getcwd()}
    job = guard.submit(kind=args.kind, slots=args.slots, cap=args.cap, argv=[*options, script, *script_args],
                       context=context)
    verified = guard.last_path_script
    print(job)
    print(f"licence_guard: submitted {verified.name} as job {job}, reservation {guard.last_reservation_id} "
          f"({args.kind}, {args.slots} slot(s), cap {args.cap}); bridge {verified.bridge_commit}, script sha256 "
          f"{verified.sha256}", file=sys.stderr)
    return 0


def main(argv: list[str] | None = None, *, scheduler: Scheduler | None = None, observer: Any = None,
         licences: LicenceServer | None = None) -> int:
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
            ok, lines = check_lock_test(args.test_dir, hosts=args.hosts, iterations=args.iterations,
                                        min_switches=args.min_switches)
            print("\n".join(lines))
            return 0 if ok else 1
        if args.command == "check-job":
            # The job's own variables against the shared guard, never defaulted: checked before the root below.
            found = check_job_guard(kinds=(args.kind,), root=args.root)
            bound = f", bound to job {found['job_id']}" if found["job_id"] else ""
            print(f"GUARD: OK {found['root']} reservation {found['reservation']} ({found['kind']}, "
                  f"{found['state']}{bound}); this job {found['slurm_job_id'] or 'none'}")
            return 0
        root = args.root if args.root is not None else DEFAULT_ROOT
        _check_inherited_root(root)
        if args.command == "lock-test-worker":
            return run_lock_test_worker(root, args.test_dir, expect_commit_file=args.expect_commit_file,
                                        iterations=args.iterations, hosts=args.hosts,
                                        barrier_timeout_s=args.barrier_timeout_s, hold_s=args.hold_s,
                                        pause_s=args.pause_s, lock_timeout_s=args.lock_timeout_s)
        guard = Guard(root, scheduler if scheduler is not None else SlurmScheduler(), observer=observer,
                      licences=licences)
        if args.command == "status":
            if args.lock_info:
                info = guard.lock_info()
                print(json.dumps(info, indent=2, sort_keys=True, default=str))
                print(f"LOCK HELD {info['token']}" if info["held"] else "LOCK FREE")
                return 0
            if args.job is not None:
                print(format_job_leases(guard.root, args.job, guard.leases_of_job(args.job)))
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
        if args.command == "release":
            guard.release_lease(args.lease, confirm_cleanup=args.confirm_cleanup)
            print(f"RELEASED {args.lease} ({guard.lease(args.lease).get('outcome')}); the reason, who ran it and what "
                  f"the guard saw are in {guard.root / 'events.log'}")
            return 0
        if args.command == "scope-check":
            ok, lines = scope_check(guard, args.lease)
            print("\n".join(lines))
            return 0 if ok else 1
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
