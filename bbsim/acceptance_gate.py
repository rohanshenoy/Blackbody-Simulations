"""The assembly gate's questions about a finished job (CDX-039, CDX-038, OR-027; ledger OR-040, OR-042, HF-057).

Which acceptance does the job require? ``required_acceptance`` reads the requirement the job's manifest recorded at
solve time, its config's ``[acceptance]`` table: ``kind`` is ``none``, ``reexport_check`` (a ``--from-solved``
re-export, judged against its source), ``grid_acceptance`` (a 35-key solve, judged against the 15-key solve of the
same geometry) or ``array_run`` (a job-array task, judged by its own self-check against its frozen plan, with no
reference job; HF-057), with the ``reference_job`` that judges it for the two reference-based kinds and, for
``grid_acceptance``, an optional ``t_rel_tol``. A job made before CDX-039 recorded none, so three explicit rules
decide for it:

1. its sidecar records ``provenance.reexport``: ``reexport_check`` against that source job;
2. otherwise a schema-1 job: ``none``;
3. otherwise a schema-2 solve: ``grid_acceptance``, since every such solve was a grid35 solve; its reference is the one
   its report names.

The answer never depends on which report files exist, so a 35-key job whose checker never ran is refused, not taken
as needing nothing (CDX-039). A record that is present but malformed is refused; it never falls back to these rules,
and neither does an array task (a manifest with a ``plan`` block), which must declare its acceptance.

Is the job's report still about these bytes? (CDX-038.) A checker reads every input through one ``InputCapture``,
which records the sha256 of the very bytes each reader parsed, and writes them into its report as the ``binding``
block (``binding_block``, version ``BINDING_VERSION``). ``report_problems`` refuses a report whose binding is absent
(a report made before CDX-038), names another job or reference, used other tolerances, leaves an evidence file
(``evidence_files``) uncovered or lists a file whose bytes have changed since; ``sidecar_file_problems`` checks the
job's CSVs and manifests against its own sidecar. Assembly checks both before copying anything.

A grid report must also hold the pointwise pair rows at the severities of the current policy (``expected_rows``,
CDX-026, S42), and its tolerances include that policy, so a report made under another policy is refused.

An ``array_run`` report binds no reference job (``binding.reference_job`` is null, never the subject as a pretend
reference). It is judged only against the caller's ``PlanContext``: the plan file, the sha256 the caller recorded for
it, the task index and the plan's policy (``production`` or ``probe``), never values the report asserts alone. The
plan file is hashed through the same capture and must have that sha256; the report's ``plan`` block must name the
same four; its evidence is the job's core files, each configured polarization's convergence table, the plan, its
frozen config copy and, when the plan records one, its column-order reference copy. ``CHECKERS``, ``REPORT_NAMES``,
``evidence_files``, ``binding_block``, ``expected_tolerances`` and ``expected_rows`` dispatch on every kind
explicitly; an unknown kind is a ``GateError``, never a fall-through to another kind's rules. One override exists,
scoped to the plan's own assembler: ``report_problems(..., waive={"convergence"})`` accepts a failed ``array_run``
report whose only failed row is ``convergence``, recorded as waivable, beside a passed ``convergence.records``.

Did the run converge? (OR-027, CDX-025.) ``convergence_problems`` is the assembly's convergence gate, the one
implementation the bridge's ``assemble_tree.sbatch`` and ``bbsim.assemble`` share: every expected polarization's
record in ``solver.convergence_qualification`` complete and well typed (the runner's projection of
``bbsim.convergence.qualify``, ``CONVERGENCE_RECORD_FIELDS``), finite nonnegative achieved deltas, the configured
target, verdicts consistent with their own fields and with the aggregate ``solver.qualified``. Only an ordinary
failure to qualify (a missed target, or AEDT's own ``Converged : No``) on otherwise usable evidence can be allowed,
and only by the caller's explicit ``allow_unqualified``.

Stdlib only: the bridge's ``assemble_tree.sbatch`` imports this module from the checkout, and
``bbsim.config`` validates its ``[acceptance]`` table with ``validate_requirement_fields``, so this module must not
import the config. ``bbsim.run_acceptance`` and ``bbsim.assemble`` import this module too; it never imports them, or
the other checkers, at module level, only inside ``expected_tolerances`` and ``expected_rows``, since they import it.
"""
from __future__ import annotations

import hashlib
import io
import json
import math
import os
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, TypeVar

KINDS = ("none", "reexport_check", "grid_acceptance", "array_run")
REPORT_NAMES = {"reexport_check": "reexport_check.json", "grid_acceptance": "grid_acceptance.json",
                "array_run": "array_run.json"}
FIELDS = ("kind", "reference_job", "t_rel_tol")
CHECKERS = {"reexport_check": "bbsim.reexport_check", "grid_acceptance": "bbsim.grid_acceptance",
            "array_run": "bbsim.run_acceptance"}
REFERENCED = ("reexport_check", "grid_acceptance")   # the kinds judged against a reference job
BINDING_VERSION = "bbsim-acceptance-binding/1"   # the bridge's assemble_tree.sbatch expects exactly this
CHUNK = 1 << 20                                  # files read only for their hash are streamed in 1 MiB chunks
EPHI_FILES = ("waveguide.csv", "far_field.csv", "manifest.json")
# An array plan's policy and the polarizations its config exports (HF-057): production solves and exports both, the
# schema-1 probe exports Ephi=0 only (both are solved, HF-058). The planner refuses any other combination.
ARRAY_POLICIES = {"production": (0, 1), "probe": (0,)}
WAIVABLE_ROWS = frozenset({"convergence"})       # the only row an override can waive (OR-027, HF-057)
# The runner's record of one polarization's convergence, solver.convergence_qualification["<ephi>"]: its projection of
# bbsim.convergence.qualify. qualify's criterion is not stored; a wrong one shows only in reasons.
CONVERGENCE_RECORD_FIELDS = ("achieved", "target", "passes_completed", "max_passes", "converged", "qualified",
                             "reasons")
UNQUALIFIED = "solver.qualified is false, not true; assembly takes qualified runs only (OR-027)"
T = TypeVar("T")


class GateError(ValueError):
    """A job's acceptance cannot be determined or contradicts its own record: the job is refused."""


class InputChanged(GateError):
    """A file read again during one check or one assembly holds other bytes than at its first read."""


@dataclass(frozen=True)
class Requirement:
    kind: str
    reference_job: str | None      # resolved, or None (kinds none and array_run, a schema-2 solve before CDX-039)
    t_rel_tol: float | None
    source: str                     # how the requirement was determined, printed by the gate


def _is_sha256(value: object) -> bool:
    return isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


@dataclass(frozen=True)
class PlanContext:
    """What the caller expects an ``array_run`` report to have judged (HF-057): the frozen plan file (an
    absolute path, as every path a plan records), the sha256 the caller recorded for it, the task index and the plan's
    policy, ``production`` or ``probe``. The task's checker and the plan's assembler supply it; a report's own claims
    never replace it. Raises ``GateError`` for a malformed field."""
    plan_path: str
    plan_sha256: str
    index: int
    policy: str

    def __post_init__(self) -> None:
        if not isinstance(self.plan_path, (str, os.PathLike)) or not os.path.isabs(os.fspath(self.plan_path)):
            raise GateError(f"a plan context's plan_path must be an absolute path, got {self.plan_path!r}")
        if not _is_sha256(self.plan_sha256):
            raise GateError(f"a plan context's plan_sha256 must be 64 lowercase hex digits, got {self.plan_sha256!r}")
        if isinstance(self.index, bool) or not isinstance(self.index, int) or self.index < 0:
            raise GateError(f"a plan context's task index must be a nonnegative integer, got {self.index!r}")
        if self.policy not in ARRAY_POLICIES:
            raise GateError(f"a plan context's policy must be one of {', '.join(ARRAY_POLICIES)}, got {self.policy!r}")


def validate_requirement_fields(kind: object, reference_job: object, t_rel_tol: object) -> None:
    """The rules of one ``[acceptance]`` record, shared by the config loader and the assembly gate: ``kind`` one of
    ``KINDS``; ``reference_job`` a nonblank path for a reference-based kind and absent for ``none`` and ``array_run``
    (HF-057); ``t_rel_tol`` absent, or for ``grid_acceptance`` a finite positive number (never a bool). Raises
    ``GateError`` naming the field."""
    if not isinstance(kind, str) or kind not in KINDS:
        raise GateError(f"acceptance kind must be one of {', '.join(KINDS)}, got {kind!r}")
    if kind not in REFERENCED:
        if reference_job is not None:
            raise GateError(f"acceptance kind {kind} takes no reference_job, got {reference_job!r}")
    elif reference_job is None:
        raise GateError(f"acceptance kind {kind} needs a reference_job, the job that judges it")
    elif not isinstance(reference_job, (str, os.PathLike)) or not str(reference_job).strip():
        raise GateError(f"acceptance reference_job must be a nonblank path, got {reference_job!r}")
    if t_rel_tol is None:
        return
    if kind != "grid_acceptance":
        raise GateError(f"acceptance t_rel_tol applies to grid_acceptance only, not to kind {kind}")
    if isinstance(t_rel_tol, bool) or not isinstance(t_rel_tol, (int, float)) or not (
            math.isfinite(t_rel_tol) and t_rel_tol > 0):
        raise GateError(f"acceptance t_rel_tol must be a finite positive number, got {t_rel_tol!r}")


def _resolved(path: str | os.PathLike) -> str:
    return str(Path(path).resolve())


def _schema_major(sidecar: Mapping, job: Path) -> str:
    version = sidecar.get("schema_version")
    if not isinstance(version, str) or not version:
        raise GateError(f"{job}: its sidecar's schema_version is {version!r}, so its acceptance cannot be decided")
    return version.split(".")[0]


def _reexport_source(sidecar: Mapping, job: Path) -> str | None:
    provenance = sidecar.get("provenance")
    reexport = provenance.get("reexport") if isinstance(provenance, Mapping) else None
    if reexport is None:
        return None
    source = reexport.get("source_job") if isinstance(reexport, Mapping) else None
    if not isinstance(source, str) or not source.strip():
        raise GateError(f"{job}: its sidecar's provenance.reexport has no usable source_job ({reexport!r})")
    return _resolved(source)


def required_acceptance(job: Path, manifest: Mapping, sidecar: Mapping) -> Requirement:
    """The acceptance ``job`` requires, from its manifest's ``config.acceptance`` (the module docstring's rules for a
    job made before CDX-039, never for an array task). Raises ``GateError`` naming the job and the field for a
    malformed record, one that contradicts the sidecar's ``provenance.reexport``, or an array task declaring none."""
    if not isinstance(manifest, Mapping):
        raise GateError(f"{job}: its manifest must be an object, got {type(manifest).__name__}")
    config = manifest.get("config")
    if config is not None and not isinstance(config, Mapping):
        raise GateError(f"{job}: its manifest's config must be an object, got {type(config).__name__}")
    if "plan" in manifest and not (isinstance(config, Mapping) and "acceptance" in config):
        raise GateError(f"{job}: its manifest records a plan block, so it is a job-array task, but no "
                        "config.acceptance; an array task declares its acceptance, and the rules for jobs made before "
                        "CDX-039 never apply to it (HF-057)")
    source = _reexport_source(sidecar, job)
    if isinstance(config, Mapping) and "acceptance" in config:
        record = config["acceptance"]
        if not isinstance(record, Mapping):
            raise GateError(f"{job}: its manifest's config.acceptance must be an object, got {record!r}")
        unknown = sorted(set(record) - set(FIELDS))
        if unknown:
            raise GateError(f"{job}: its manifest's config.acceptance holds unknown fields {unknown}")
        if "kind" not in record:
            raise GateError(f"{job}: its manifest's config.acceptance has no kind")
        kind, reference, tol = record["kind"], record.get("reference_job"), record.get("t_rel_tol")
        try:
            validate_requirement_fields(kind, reference, tol)
        except GateError as exc:
            raise GateError(f"{job}: its manifest's config.acceptance: {exc}") from None
        if reference is not None and not os.path.isabs(str(reference)):
            raise GateError(f"{job}: its manifest's config.acceptance reference_job {reference!r} is not absolute")
        if kind == "reexport_check":
            if source is None:
                raise GateError(f"{job}: it records kind reexport_check, but its sidecar has no provenance.reexport")
            if _resolved(reference) != source:
                raise GateError(f"{job}: it records reexport_check against {reference}, but its sidecar's "
                                f"provenance.reexport source_job is {source}")
        elif source is not None:
            raise GateError(f"{job}: it records kind {kind}, but its sidecar's provenance.reexport makes it a "
                            "re-export, which requires reexport_check")
        return Requirement(kind, None if reference is None else _resolved(reference),
                           None if tol is None else float(tol), "recorded at solve time")
    if source is not None:
        return Requirement("reexport_check", source, None, "a re-export made before CDX-039 (provenance.reexport)")
    if _schema_major(sidecar, job) == "1":
        return Requirement("none", None, None, "a schema-1 job made before CDX-039")
    return Requirement("grid_acceptance", None, None, "a schema-2 solve made before CDX-039; every such solve is a "
                       "grid35 solve, its reference is the one its report names")


# --- CDX-038: what a report judged, and whether those bytes are still there ------------------------------------------------

def _key(path: str | os.PathLike) -> str:
    return str(Path(path).resolve())


class InputCapture:
    """The bytes one checker invocation (or one assembly) judged, by resolved path: the sha256 and size of the first
    read of each file. A later read of the same file must find the same bytes, else ``InputChanged``. One file's bytes
    are held only while its reader parses them; files read only for their hash are streamed."""

    def __init__(self) -> None:
        self.hashes: dict[str, str] = {}
        self.sizes: dict[str, int] = {}
        self.problems: list[str] = []     # every InputChanged raised, kept for the checker's binding.inputs row
        self.retained_bytes = 0            # nonzero only while a reader runs

    def _record(self, key: str, digest: str, size: int) -> None:
        first = self.hashes.get(key)
        if first is None:
            self.hashes[key], self.sizes[key] = digest, size
        elif first != digest:
            problem = (f"{key} changed while it was being checked: sha256 {digest} on a later read, {first} on the "
                       "first")
            self.problems.append(problem)
            raise InputChanged(problem)

    def parse(self, path: str | os.PathLike, reader: Callable[[BinaryIO], T]) -> T:
        """``reader`` applied to a fresh stream over the file's bytes, which are hashed first and dropped after."""
        key = _key(path)
        data = Path(key).read_bytes()
        self._record(key, hashlib.sha256(data).hexdigest(), len(data))
        self.retained_bytes = len(data)
        try:
            return reader(io.BytesIO(data))
        finally:
            self.retained_bytes = 0
            del data

    def read_text(self, path: str | os.PathLike) -> str:
        return self.parse(path, lambda f: f.read().decode("utf-8"))

    def digest(self, path: str | os.PathLike) -> str:
        key = _key(path)
        digest, size = hashlib.sha256(), 0
        with open(key, "rb") as fh:
            while chunk := fh.read(CHUNK):
                digest.update(chunk)
                size += len(chunk)
        self._record(key, digest.hexdigest(), size)
        return self.hashes[key]

    def size(self, path: str | os.PathLike) -> int:
        key = _key(path)
        if key not in self.sizes:
            self.digest(key)
        return self.sizes[key]

    def recorded(self, path: str | os.PathLike) -> str | None:
        """The digest of the first read of ``path``, or None if it was not read."""
        return self.hashes.get(_key(path))


def job_files(job: Path) -> list[Path]:
    """The core files a binding covers, those that exist: the job's manifest.json, its sidecars and, in every Ephi
    directory, the two CSVs and the manifest copy."""
    job = Path(job)
    files = [job / "manifest.json", *job.glob("*.dataset.json")]
    for d in job.glob("*_Ephi=*"):
        if d.is_dir():
            files += [d / name for name in EPHI_FILES]
    return sorted(p for p in files if p.is_file())


def _plan_copy(plan: object, name: str, *, required: bool) -> Path | None:
    """The absolute ``copy`` path of the plan's ``name`` input record ``{source, copy, sha256}``, or None for an
    optional one the plan records as null. Raises ``GateError`` for anything else."""
    if not isinstance(plan, Mapping):
        raise GateError(f"the plan must be an object to name its {name} copy, got {type(plan).__name__}")
    entry = plan.get(name)
    if entry is None:
        if required:
            raise GateError(f"the plan records no {name} {{source, copy, sha256}}")
        return None
    copy = entry.get("copy") if isinstance(entry, Mapping) else None
    if not isinstance(copy, str) or not os.path.isabs(copy):
        raise GateError(f"the plan's {name}.copy must be an absolute path, as every path a plan records is, got "
                        f"{copy!r}")
    return Path(copy)


def evidence_files(kind: str, subject: Path, reference: Path | None, subject_sidecar: Mapping,
                   poynting: Path | None = None, *, context: PlanContext | None = None,
                   plan: Mapping | None = None) -> list[Path]:
    """Every file a ``kind`` report must bind: both jobs' ``job_files``; for ``reexport_check`` also both jobs'
    convergence tables, the subject's kept raw column-order tables when its sidecar's method is ``reference``, and the
    --poynting record when one was used. For ``array_run`` (no reference job; ``context`` and the parsed ``plan``
    required): the subject's ``job_files``, the convergence table of each polarization the plan's policy configures,
    the plan file, its config copy and, when the plan records one, its column-order reference copy. Required evidence
    is listed even when it is missing, so that a missing file refuses rather than dropping out of the binding. Raises
    ``GateError`` for an unknown kind or arguments the kind does not take."""
    subject = Path(subject)
    if kind == "array_run":
        if reference is not None:
            raise GateError(f"an array_run report binds no reference job, got {reference}")
        if context is None or plan is None:
            raise GateError("an array_run report's evidence is listed from the caller's plan context and the plan")
        files = job_files(subject)
        files += [subject / "logs" / f"convergence_Ephi{e}.txt" for e in ARRAY_POLICIES[context.policy]]
        files += [Path(context.plan_path), _plan_copy(plan, "config", required=True)]
        order = _plan_copy(plan, "column_order_reference", required=False)
        if order is not None:
            files.append(order)
        return sorted(set(files))
    if kind not in REFERENCED:
        raise GateError(f"acceptance kind {kind!r} has no report evidence")
    if reference is None:
        raise GateError(f"a {kind} report binds its reference job, and none was given")
    if context is not None or plan is not None:
        raise GateError(f"a plan context applies to array_run reports only, not to {kind}")
    reference = Path(reference)
    files = job_files(subject) + job_files(reference)
    if kind == "reexport_check":
        files += [job / "logs" / f"convergence_Ephi{e}.txt" for job in (subject, reference) for e in (0, 1)]
        exit_field = subject_sidecar.get("exit_field") if isinstance(subject_sidecar, Mapping) else None
        order = exit_field.get("column_order_check") if isinstance(exit_field, Mapping) else None
        if isinstance(order, Mapping) and order.get("method") == "reference":
            files += [subject / "scratch" / f"column_order_reference_Ephi{e}.csv" for e in (0, 1)]
        if poynting is not None:
            files.append(Path(poynting))
    return sorted(set(files))


def binding_block(kind: str, subject: Path, reference: Path | None, files: Mapping[str, str]) -> dict:
    """The ``binding`` block a ``kind`` report records: the resolved subject and reference jobs and the digest of every
    file read. ``reference_job`` is null for ``array_run`` and only for it; ``GateError`` otherwise."""
    if kind == "array_run":
        if reference is not None:
            raise GateError(f"an array_run binding names no reference job (HF-057), got {reference}")
        reference_job = None
    elif kind in REFERENCED:
        if reference is None:
            raise GateError(f"a {kind} binding names the reference job that judged it; only array_run binds none")
        reference_job = _key(reference)
    else:
        raise GateError(f"acceptance kind {kind!r} has no report binding")
    return {"version": BINDING_VERSION, "kind": kind, "subject_job": _key(subject), "reference_job": reference_job,
            "files": dict(sorted(files.items()))}


def _array_policy(kind: str, policy: str | None) -> None:
    if kind == "array_run":
        if policy not in ARRAY_POLICIES:
            raise GateError(f"array_run is judged under its plan's policy, one of {', '.join(ARRAY_POLICIES)}; got "
                            f"{policy!r}")
    elif kind in REFERENCED:
        if policy is not None:
            raise GateError(f"a policy applies to array_run reports only, not to {kind}")
    else:
        raise GateError(f"acceptance kind {kind!r} has no report to judge")


def expected_tolerances(kind: str, t_rel_tol: float | None, policy: str | None = None) -> dict[str, object]:
    """The tolerances a ``kind`` report must have used: the checker's defaults, with the job's recorded t_rel_tol for
    grid acceptance (whose tolerances include the pointwise comparison's policy), and the plan's policy for an array
    run (``bbsim.run_acceptance``'s versioned policy). An unknown kind raises ``GateError``."""
    _array_policy(kind, policy)
    if kind == "reexport_check":
        from bbsim import reexport_check
        return reexport_check.tolerances()
    if kind == "grid_acceptance":
        from bbsim import grid_acceptance
        return grid_acceptance.tolerances(grid_acceptance.T_REL_TOL if t_rel_tol is None else t_rel_tol)
    if t_rel_tol is not None:
        raise GateError(f"t_rel_tol applies to grid_acceptance only, not to array_run (got {t_rel_tol!r})")
    from bbsim import run_acceptance
    return run_acceptance.tolerances(policy)


def expected_rows(kind: str, policy: str | None = None) -> dict[str, str]:
    """The rows a ``kind`` report must hold, by name, at their severity: the pointwise pair rows for grid acceptance
    (CDX-026, S42), none for a re-export check, the policy's rows for an array run (HF-057). An unknown kind raises
    ``GateError``."""
    _array_policy(kind, policy)
    if kind == "reexport_check":
        return {}
    if kind == "grid_acceptance":
        from bbsim import grid_acceptance
        return grid_acceptance.required_rows()
    from bbsim import run_acceptance
    return run_acceptance.required_rows(policy)


def _read_object(inputs: InputCapture, path: Path) -> dict:
    data = inputs.parse(path, json.load)
    if not isinstance(data, dict):
        raise ValueError("not a JSON object")
    return data


def _whose(path: str, subject: str | None, reference: str | None) -> str:
    if subject and path.startswith(subject + os.sep):
        return "the subject job's"
    if reference and path.startswith(reference + os.sep):
        return "the reference job's"
    return "the evidence file"


def _show(value: object) -> str:
    try:
        return json.dumps(value)
    except (TypeError, ValueError):
        return repr(value)


def _only(checks: list, name: str) -> Mapping | None:
    rows = [c for c in checks if c.get("name") == name]
    return rows[0] if len(rows) == 1 else None


def _waiver_refusal(report: Mapping, checks: list, failed: list[str]) -> str | None:
    """Why the convergence waiver does not apply to this ``array_run`` report, or None when it does: a ``passed:
    false`` verdict whose only failed fail-severity row is ``convergence``, recorded waivable, beside a passed
    ``convergence.records`` (OR-027, HF-057)."""
    if report.get("passed") is not False:
        return f"its verdict records passed {_show(report.get('passed'))}; only a failed verdict is waived"
    if failed != ["convergence"]:
        return f"the rows {failed} failed; only a failure of the convergence row alone is waived"
    row, records = _only(checks, "convergence"), _only(checks, "convergence.records")
    if row is None or row.get("severity") != "fail":
        return "it does not hold exactly one convergence row at severity fail"
    values = row.get("values")
    if not (isinstance(values, Mapping) and values.get("waivable") is True):
        return ("its convergence row is not recorded as waivable (a missed target or AEDT's own nonconvergence, on a "
                "finite achieved delta and convergence records that passed)")
    if records is None or records.get("severity") != "fail" or records.get("passed") is not True:
        return "it does not hold exactly one passed convergence.records row at severity fail"
    return None


def _plan_block_problems(name: str, block: object, context: PlanContext) -> list[str]:
    """Where the report's ``plan`` block ``{path, sha256, index, policy}`` differs from the caller's context."""
    if not isinstance(block, Mapping):
        return [f"{name} records no plan block {{path, sha256, index, policy}}, so it names no plan, task or policy "
                f"(got {_show(block)})"]
    problems = []
    path = block.get("path")
    if not (isinstance(path, str) and os.path.isabs(path) and _key(path) == _key(context.plan_path)):
        problems.append(f"{name} judged against the plan {path!r}, not the caller's {_key(context.plan_path)}")
    if block.get("sha256") != context.plan_sha256:
        problems.append(f"{name} judged the plan with sha256 {block.get('sha256')!r}, not the {context.plan_sha256} "
                        "the caller expects")
    index = block.get("index")
    if isinstance(index, bool) or not isinstance(index, int) or index != context.index:
        problems.append(f"{name} judged task {index!r} of its plan, not the caller's task {context.index}")
    if block.get("policy") != context.policy:
        problems.append(f"{name} judged under the policy {block.get('policy')!r}, not the plan's {context.policy!r}")
    return problems


def report_problems(report: object, kind: str, job: Path, requirement: Requirement, inputs: InputCapture, *,
                    context: PlanContext | None = None, waive: frozenset[str] = frozenset()) -> list[str]:
    """Why ``report``, the ``kind`` report found in ``job``, does not accept the bytes now in place; empty when it
    does. A malformed report or binding is a problem, never an exception.

    An ``array_run`` report is judged only against the caller's ``context`` (HF-057): without one it is refused. Its
    ``plan`` block must equal the context, and the plan file, hashed through ``inputs``, must have the context's
    sha256. ``waive={"convergence"}`` (array_run only; any other row raises ``GateError``) accepts a ``passed: false``
    report whose only failed row is a waivable ``convergence``; every other row, the verdict's consistency, the
    context and the binding are judged as without it (OR-027, HF-057)."""
    if kind not in REPORT_NAMES:
        raise GateError(f"acceptance kind {kind!r} has no report to judge")
    waive = frozenset(waive)
    if waive - WAIVABLE_ROWS:
        raise GateError(f"only the {', '.join(sorted(WAIVABLE_ROWS))} row can be waived (OR-027), not "
                        f"{', '.join(sorted(waive - WAIVABLE_ROWS))}")
    array = kind == "array_run"
    if not array and (waive or context is not None):
        raise GateError(f"a plan context and the convergence waiver apply to array_run reports only, not to {kind}")
    name = REPORT_NAMES[kind]
    if not isinstance(report, Mapping):
        return [f"{name} is not a JSON object"]
    if array and context is None:
        return [f"{name} is judged against the caller's plan context (the plan, its sha256, the task index and the "
                "policy), and none was given: an array run is assembled by its plan's own assembler (HF-057)"]
    policy = context.policy if array else None
    problems: list[str] = []
    checks = report.get("checks")
    if not (isinstance(checks, list) and checks and all(isinstance(c, Mapping) for c in checks)):
        if report.get("passed") is not True:
            problems.append(f"{name} records passed {_show(report.get('passed'))}, not true")
        problems.append(f"{name}'s checks must be a non-empty list of rows, got {type(checks).__name__} "
                        f"{checks!r:.40}: its verdict stands on nothing")
    else:
        failed = [str(c.get("name")) for c in checks if c.get("severity") == "fail" and c.get("passed") is not True]
        refusal = _waiver_refusal(report, checks, failed) if waive else "no waiver was asked for"
        if refusal is not None:
            if report.get("passed") is not True:
                problems.append(f"{name} records passed {_show(report.get('passed'))}, not true")
            if failed:
                problems.append(f"{name}'s own rows {failed} failed, which contradicts a passed verdict")
            if waive and problems:
                problems.append(f"the convergence waiver does not apply to {name}: {refusal}")
        try:
            required = expected_rows(kind, policy)
        except Exception as exc:  # noqa: BLE001 - a checker the checkout cannot import is a refusal, not a traceback
            required = {}
            problems.append(f"the rows {kind} must hold cannot be read from the checkout ({exc})")
        for row, severity in required.items():
            if not any(c.get("name") == row and c.get("severity") == severity for c in checks):
                why = (f"the {policy} policy of {CHECKERS[kind]} requires it" if array else
                       "CDX-026: made under another pointwise policy, or before it")
                problems.append(f"{name} lacks the row {row} at severity {severity} ({why}; copy it aside and rerun "
                                "its checker, no AEDT)")
    binding = report.get("binding")
    if binding is None:
        made = "every array_run report has one" if array else "made before CDX-038"
        rerun = " with the plan, its sha256 and the task index" if array else ""
        return problems + [f"{name} carries no file binding ({made}); copy it aside, then rerun python -m "
                           f"{CHECKERS[kind]}{rerun} on this job, no AEDT, before assembling it"]
    if not isinstance(binding, Mapping):
        return problems + [f"{name}'s binding is not an object ({binding!r})"]
    if binding.get("version") != BINDING_VERSION:
        return problems + [f"{name}'s binding version {binding.get('version')!r} is not {BINDING_VERSION!r}"]
    if binding.get("kind") != kind:
        problems.append(f"{name}'s binding kind {binding.get('kind')!r} is not {kind!r}")
    subject, reference, files = binding.get("subject_job"), binding.get("reference_job"), binding.get("files")
    if not isinstance(subject, str) or subject != _key(job):
        problems.append(f"{name} judged the subject job {subject!r}, not this job {_key(job)}")
    if array:
        if reference is not None:
            problems.append(f"{name}'s binding names the reference job {reference!r}; an array run is judged against "
                            "its plan, never a reference (HF-057)")
            reference = None
    elif not isinstance(reference, str) or not reference:
        return problems + [f"{name}'s binding names no reference job ({reference!r})"]
    if not (isinstance(files, Mapping) and all(isinstance(k, str) and isinstance(v, str) for k, v in files.items())):
        return problems + [f"{name}'s binding files must map paths to sha256 digests"]
    gated_tol = None
    if requirement.kind == kind:
        gated_tol = requirement.t_rel_tol
        required_reference = requirement.reference_job
        if required_reference is not None and reference is not None and _key(reference) != required_reference:
            problems.append(f"{name} judged against the reference job {reference}, but the job requires {kind} "
                            f"against {requirement.reference_job}")
    try:
        expected = expected_tolerances(kind, gated_tol, policy)
    except Exception as exc:  # noqa: BLE001 - a checker the checkout cannot import is a refusal, not a traceback
        expected = None
        problems.append(f"the tolerances {kind} must use cannot be read from the checkout ({exc})")
    if expected is not None and report.get("tolerances") != expected:
        problems.append(f"{name} used tolerances {report.get('tolerances')!r}, not {expected!r}")
    plan = None
    if array:
        problems += _plan_block_problems(name, report.get("plan"), context)
        plan_path = Path(context.plan_path)
        try:
            plan = _read_object(inputs, plan_path)
        except InputChanged as exc:
            problems.append(str(exc))
        except (OSError, ValueError) as exc:
            problems.append(f"the plan {plan_path}, which {name} must be judged against, is unreadable ({exc})")
        else:
            if inputs.recorded(plan_path) != context.plan_sha256:
                problems.append(f"the plan {_key(plan_path)} has sha256 {inputs.recorded(plan_path)}, not the "
                                f"{context.plan_sha256} the caller expects")
    elif not Path(reference).is_dir():
        return problems + [f"{name}'s reference job {reference} is not a directory"]
    sidecars = sorted(Path(job).glob("*.dataset.json"))
    try:
        sidecar = _read_object(inputs, sidecars[0]) if len(sidecars) == 1 else {}
    except (OSError, ValueError) as exc:
        sidecar = {}
        problems.append(f"{sidecars[0]} is unreadable ({exc})")
    evidence: list[Path] = []
    if array:
        if plan is not None:
            try:
                evidence = evidence_files(kind, Path(job), None, sidecar, context=context, plan=plan)
            except GateError as exc:
                problems.append(f"the evidence {name} must bind cannot be listed from the plan ({exc})")
    else:
        poynting = report.get("poynting") if kind == "reexport_check" else None
        evidence = evidence_files(kind, Path(job), Path(reference), sidecar,
                                  Path(poynting) if isinstance(poynting, str) else None)
    for path in evidence:
        if _key(path) not in files:
            problems.append(f"{path} is not covered by the acceptance ({name} binds no digest of it)")
    for path, recorded in files.items():
        whose = _whose(path, subject if isinstance(subject, str) else None, reference)
        try:
            digest = inputs.digest(path)
        except OSError as exc:
            problems.append(f"{whose} {path}, judged by {name}, is missing or unreadable ({exc})")
            continue
        except InputChanged as exc:
            problems.append(str(exc))
            continue
        if digest != recorded:
            problems.append(f"{whose} {path} has changed since {name} judged it (sha256 {digest}, judged {recorded})")
    return problems


# --- OR-027, CDX-025: the convergence gate, shared by assemble_tree.sbatch and bbsim.assemble -------------------------

def _real(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _ordinary(reason: object) -> bool:
    """A reason ``bbsim.convergence.qualify`` gives for an ordinary failure to qualify: AEDT's own ``Converged`` line
    other than Yes, or an achieved delta above the target. Every other reason (no table, a non-finite or negative
    delta, another table target, another criterion) is unusable evidence, which no override accepts."""
    return isinstance(reason, str) and (reason.startswith("AEDT reports Converged : ") or (
        reason.startswith("achieved delta ") and " above the target " in reason))


def _configured_target(manifest: Mapping, solver: Mapping) -> tuple[float | None, list[str]]:
    """The run's configured MaxDeltaE: ``solver.max_delta_e``, which the runner records from its config, equal to
    ``config.solver.max_delta_e`` when the manifest records that too."""
    target = solver.get("max_delta_e")
    if not (_real(target) and math.isfinite(target) and target > 0):
        return None, [f"solver.max_delta_e is {_show(target)}, not a finite positive number, so no recorded target can "
                      "be checked against the configured one (OR-027)"]
    config = manifest.get("config")
    config_solver = config.get("solver") if isinstance(config, Mapping) else None
    if isinstance(config_solver, Mapping) and "max_delta_e" in config_solver:
        recorded = config_solver["max_delta_e"]
        if not (_real(recorded) and recorded == target):
            return None, [f"config.solver.max_delta_e is {_show(recorded)}, but solver.max_delta_e is {_show(target)}: "
                          "the manifest contradicts itself on the configured target (OR-027)"]
    return float(target), []


def _record_problems(where: str, record: object, target: float | None) -> tuple[list[str], bool | None]:
    """What is wrong with one polarization's record, and its ``qualified`` when that is a bool."""
    if not isinstance(record, Mapping):
        return [f"{where} is {_show(record)}, not a record of {list(CONVERGENCE_RECORD_FIELDS)}"], None
    missing = [f for f in CONVERGENCE_RECORD_FIELDS if f not in record]
    extra = sorted(str(k) for k in record if k not in CONVERGENCE_RECORD_FIELDS)
    if missing or extra:
        return [f"{where} is malformed: missing {missing}, unexpected {extra}; a record holds exactly "
                f"{list(CONVERGENCE_RECORD_FIELDS)}, the runner's projection of bbsim.convergence.qualify (CDX-025)"], (
                    record["qualified"] if isinstance(record.get("qualified"), bool) else None)
    problems: list[str] = []
    achieved, recorded_target, qualified = record["achieved"], record["target"], record["qualified"]
    converged, reasons = record["converged"], record["reasons"]
    if not (_real(achieved) and math.isfinite(achieved)):
        problems.append(f"{where}: achieved delta {_show(achieved)} is not a finite number (the table is missing, "
                        "unreadable or not finite)")
        achieved = None
    elif achieved < 0:
        problems.append(f"{where}: achieved delta {achieved!r} is negative")
        achieved = None
    if not (_real(recorded_target) and math.isfinite(recorded_target)):
        problems.append(f"{where}: target {_show(recorded_target)} is not a finite number (the table shows no target)")
    elif target is not None and not math.isclose(recorded_target, target, rel_tol=1e-9):
        problems.append(f"{where}: target {recorded_target!r} is not the configured max_delta_e {target!r}")
    for field in ("passes_completed", "max_passes"):
        value = record[field]
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            problems.append(f"{where}: {field} {_show(value)} is not a nonnegative integer")
    if not (isinstance(converged, str) and converged.strip()):
        problems.append(f"{where}: converged {_show(converged)} is not AEDT's Converged line")
        converged = None
    if not (isinstance(reasons, list) and all(isinstance(r, str) for r in reasons)):
        problems.append(f"{where}: reasons {_show(reasons)} is not a list of strings")
        reasons = None
    if not isinstance(qualified, bool):
        problems.append(f"{where}: qualified {_show(qualified)} is not true or false")
        return problems, None
    if qualified:
        if reasons:
            problems.append(f"{where} records qualified true beside the reasons {reasons}")
        if converged is not None and converged.lower() != "yes":
            problems.append(f"{where} records qualified true, but Converged : {converged}")
        if achieved is not None and target is not None and achieved > target:
            problems.append(f"{where} records qualified true, but its achieved delta {achieved!r} is above the "
                            f"configured target {target!r}")
    elif reasons is not None:
        unusual = [r for r in reasons if not _ordinary(r)]
        if not reasons:
            problems.append(f"{where} records qualified false and gives no reason")
        elif unusual:
            problems.append(f"{where} does not qualify for {unusual}: unusable evidence, not a missed target or AEDT's "
                            "own nonconvergence, which no override accepts (OR-027)")
    return problems, qualified


def convergence_problems(manifest: Mapping, allow_unqualified: bool = False,
                         expected_polarizations: tuple[int, ...] = (0, 1)) -> list[str]:
    """Why the run ``manifest`` records may not be assembled for its convergence (OR-027, CDX-025); empty when it may.

    Every expected polarization has a record in ``solver.convergence_qualification`` (keyed ``"0"`` and ``"1"`` by
    default, as the runner writes them, and no other key) holding exactly ``CONVERGENCE_RECORD_FIELDS``, well typed:
    a finite nonnegative
    achieved delta, the configured target (``solver.max_delta_e``), pass counts, AEDT's Converged line, a verdict
    consistent with those and with its reasons, and an aggregate ``solver.qualified`` equal to all of them. Those
    problems always stand. A run that is consistently unqualified for an ordinary reason (a missed target, or AEDT's
    own ``Converged : No``) gives ``UNQUALIFIED``, the one problem ``allow_unqualified`` removes. Standard library only:
    the bridge's ``assemble_tree.sbatch`` calls it with its strict default. ``expected_polarizations`` other than a
    nonempty subset of (0, 1) raises ``GateError``."""
    pols = tuple(expected_polarizations)
    if not pols or len(set(pols)) != len(pols) or any(type(e) is not int or e not in (0, 1) for e in pols):
        raise GateError(f"expected_polarizations must be a nonempty subset of (0, 1), got {expected_polarizations!r}")
    if not isinstance(manifest, Mapping):
        return [f"the manifest must be an object, got {type(manifest).__name__}, so its convergence is unknown"]
    solver = manifest.get("solver")
    if not isinstance(solver, Mapping):
        return [f"the manifest records no solver block (got {_show(solver)}), so its convergence is unknown (OR-027)"]
    target, problems = _configured_target(manifest, solver)
    expected = [str(e) for e in pols]
    records = solver.get("convergence_qualification")
    verdicts: dict[str, bool] = {}
    if not isinstance(records, Mapping):
        problems.append(f"solver.convergence_qualification is {_show(records)}, not one record per polarization "
                        f"{expected}, so no polarization's convergence is known (CDX-025)")
    else:
        for key in records:
            if key not in expected:
                problems.append(f"solver.convergence_qualification records Ephi={key}, which is not among the expected "
                                f"polarizations {expected}")
        for key in expected:
            if key not in records:
                problems.append(f"solver.convergence_qualification has no record for Ephi={key}")
                continue
            found, verdict = _record_problems(f"solver.convergence_qualification Ephi={key}", records[key], target)
            problems += found
            if verdict is not None:
                verdicts[key] = verdict
    aggregate = solver.get("qualified")
    if not isinstance(aggregate, bool):
        problems.append(f"solver.qualified is {_show(aggregate)}, not true; assembly takes qualified runs only "
                        "(OR-027), and a run records true or false")
    elif len(verdicts) == len(expected) and aggregate != all(verdicts.values()):
        problems.append(f"solver.qualified is {_show(aggregate)}, but the per-polarization records say "
                        f"{_show({k: verdicts[k] for k in expected})}: the manifest contradicts itself (OR-027)")
    elif aggregate is False and not allow_unqualified:
        problems.append(UNQUALIFIED)
    return problems


def sidecar_file_problems(job: Path, sidecar: Mapping, inputs: InputCapture) -> list[str]:
    """What in ``job`` disagrees with its own sidecar: a listed file whose sha256 or size differs, a CSV the sidecar
    does not list, a manifest.json that is not the one the sidecar's provenance.manifest_sha256 names."""
    job = Path(job)
    problems: list[str] = []
    files = sidecar.get("files") if isinstance(sidecar, Mapping) else None
    if not isinstance(files, Mapping):
        return [f"{job}: its sidecar lists no files"]
    for rel, entry in sorted(files.items()):
        path = job / rel
        if not (isinstance(entry, Mapping) and isinstance(entry.get("sha256"), str)):
            problems.append(f"{path}: the sidecar's entry {entry!r} has no sha256")
            continue
        try:
            digest, size = inputs.digest(path), inputs.size(path)
        except OSError as exc:
            problems.append(f"{path}: listed in the sidecar, missing or unreadable ({exc})")
            continue
        except InputChanged as exc:
            problems.append(str(exc))
            continue
        if digest != entry["sha256"]:
            problems.append(f"{path}: sha256 {digest}, the sidecar records {entry['sha256']}")
        if entry.get("bytes") != size:
            problems.append(f"{path}: {size} bytes, the sidecar records {entry.get('bytes')!r}")
    for path in job_files(job):
        if path.suffix == ".csv" and path.relative_to(job).as_posix() not in files:
            problems.append(f"{path}: not listed in the sidecar's files")
    provenance = sidecar.get("provenance")
    recorded = provenance.get("manifest_sha256") if isinstance(provenance, Mapping) else None
    try:
        root = inputs.digest(job / "manifest.json")
    except (OSError, InputChanged) as exc:
        return problems + [f"{job / 'manifest.json'}: unreadable ({exc})"]
    if recorded != root:
        problems.append(f"{job / 'manifest.json'}: sha256 {root}, the sidecar's provenance.manifest_sha256 is "
                        f"{recorded!r}")
    for d in sorted(job.glob("*_Ephi=*")):
        copy = d / "manifest.json"
        try:
            digest = inputs.digest(copy)
        except (OSError, InputChanged) as exc:
            problems.append(f"{copy}: unreadable ({exc})")
            continue
        if digest != root:
            problems.append(f"{copy}: sha256 {digest}, not the job manifest's {root}")
    return problems
