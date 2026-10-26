"""The assembly gate's two questions about a finished job (CDX-039, CDX-038; ledger OR-040, OR-042).

Which acceptance does the job require? ``required_acceptance`` reads the requirement the job's manifest recorded at
solve time, its config's ``[acceptance]`` table: ``kind`` is ``none``, ``reexport_check`` (a ``--from-solved``
re-export, judged against its source) or ``grid_acceptance`` (a 35-key solve, judged against the 15-key solve of the
same geometry), with the ``reference_job`` that judges it and, for ``grid_acceptance``, an optional ``t_rel_tol``. A
job made before CDX-039 recorded none, so three explicit rules decide for it:

1. its sidecar records ``provenance.reexport``: ``reexport_check`` against that source job;
2. otherwise a schema-1 job: ``none``;
3. otherwise a schema-2 solve: ``grid_acceptance``, since every such solve was a grid35 solve; its reference is the one
   its report names.

The answer never depends on which report files exist, so a 35-key job whose checker never ran is refused, not taken
as needing nothing (CDX-039). A record that is present but malformed is refused; it never falls back to these rules.

Is the job's report still about these bytes? (CDX-038.) A checker reads every input through one ``InputCapture``,
which records the sha256 of the very bytes each reader parsed, and writes them into its report as the ``binding``
block (``binding_block``, version ``BINDING_VERSION``). ``report_problems`` refuses a report whose binding is absent
(a report made before CDX-038), names another job or reference, used other tolerances, leaves an evidence file
(``evidence_files``) uncovered or lists a file whose bytes have changed since; ``sidecar_file_problems`` checks the
job's CSVs and manifests against its own sidecar. Assembly checks both before copying anything.

A grid report must also hold the pointwise pair rows at the severities of the current policy (``expected_rows``,
CDX-026, S42), and its tolerances include that policy, so a report made under another policy is refused.

Stdlib only: the bridge's ``assemble_tree.sbatch`` imports this module from the checkout, and
``bbsim.config`` validates its ``[acceptance]`` table with ``validate_requirement_fields``, so this module must not
import the config, and imports the checkers only inside ``expected_tolerances`` and ``expected_rows``, since they
import it.
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

KINDS = ("none", "reexport_check", "grid_acceptance")
REPORT_NAMES = {"reexport_check": "reexport_check.json", "grid_acceptance": "grid_acceptance.json"}
FIELDS = ("kind", "reference_job", "t_rel_tol")
CHECKERS = {"reexport_check": "bbsim.reexport_check", "grid_acceptance": "bbsim.grid_acceptance"}
BINDING_VERSION = "bbsim-acceptance-binding/1"   # the bridge's assemble_tree.sbatch expects exactly this
CHUNK = 1 << 20                                  # files read only for their hash are streamed in 1 MiB chunks
EPHI_FILES = ("waveguide.csv", "far_field.csv", "manifest.json")
T = TypeVar("T")


class GateError(ValueError):
    """A job's acceptance cannot be determined or contradicts its own record: the job is refused."""


class InputChanged(GateError):
    """A file read again during one check or one assembly holds other bytes than at its first read."""


@dataclass(frozen=True)
class Requirement:
    kind: str
    reference_job: str | None      # resolved, or None (kind none, or a schema-2 solve made before CDX-039)
    t_rel_tol: float | None
    source: str                     # how the requirement was determined, printed by the gate


def validate_requirement_fields(kind: object, reference_job: object, t_rel_tol: object) -> None:
    """The rules of one ``[acceptance]`` record, shared by the config loader and the assembly gate: ``kind`` one of
    ``KINDS``; ``reference_job`` a nonblank path for a gated kind and absent for ``none``; ``t_rel_tol`` absent, or for
    ``grid_acceptance`` a finite positive number (never a bool). Raises ``GateError`` naming the field."""
    if not isinstance(kind, str) or kind not in KINDS:
        raise GateError(f"acceptance kind must be one of {', '.join(KINDS)}, got {kind!r}")
    if kind == "none":
        if reference_job is not None:
            raise GateError(f"acceptance kind none takes no reference_job, got {reference_job!r}")
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
    job made before CDX-039). Raises ``GateError`` naming the job and the field for a malformed record or one that
    contradicts the sidecar's ``provenance.reexport``."""
    if not isinstance(manifest, Mapping):
        raise GateError(f"{job}: its manifest must be an object, got {type(manifest).__name__}")
    config = manifest.get("config")
    if config is not None and not isinstance(config, Mapping):
        raise GateError(f"{job}: its manifest's config must be an object, got {type(config).__name__}")
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


def evidence_files(kind: str, subject: Path, reference: Path, subject_sidecar: Mapping,
                   poynting: Path | None = None) -> list[Path]:
    """Every file a ``kind`` report must bind: both jobs' ``job_files``; for ``reexport_check`` also both jobs'
    convergence tables, the subject's kept raw column-order tables when its sidecar's method is ``reference``, and the
    --poynting record when one was used. Required evidence is listed even when it is missing, so that a missing file
    refuses rather than dropping out of the binding."""
    subject, reference = Path(subject), Path(reference)
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


def binding_block(kind: str, subject: Path, reference: Path, files: Mapping[str, str]) -> dict:
    return {"version": BINDING_VERSION, "kind": kind, "subject_job": _key(subject), "reference_job": _key(reference),
            "files": dict(sorted(files.items()))}


def expected_tolerances(kind: str, t_rel_tol: float | None) -> dict[str, object]:
    """The tolerances a ``kind`` report must have used: the checker's defaults, with the job's recorded t_rel_tol for
    grid acceptance (whose tolerances include the pointwise comparison's policy)."""
    if kind == "reexport_check":
        from bbsim import reexport_check
        return reexport_check.tolerances()
    from bbsim import grid_acceptance
    return grid_acceptance.tolerances(grid_acceptance.T_REL_TOL if t_rel_tol is None else t_rel_tol)


def expected_rows(kind: str) -> dict[str, str]:
    """The rows a ``kind`` report must hold, by name, at their severity: the pointwise pair rows for grid acceptance
    (CDX-026, S42), none for a re-export check."""
    if kind == "reexport_check":
        return {}
    from bbsim import grid_acceptance
    return grid_acceptance.required_rows()


def _read_object(inputs: InputCapture, path: Path) -> dict:
    data = inputs.parse(path, json.load)
    if not isinstance(data, dict):
        raise ValueError("not a JSON object")
    return data


def _whose(path: str, subject: str, reference: str) -> str:
    if path.startswith(subject + os.sep):
        return "the subject job's"
    if path.startswith(reference + os.sep):
        return "the reference job's"
    return "the evidence file"


def report_problems(report: object, kind: str, job: Path, requirement: Requirement,
                    inputs: InputCapture) -> list[str]:
    """Why ``report``, the ``kind`` report found in ``job``, does not accept the bytes now in place; empty when it
    does. A malformed report or binding is a problem, never an exception."""
    name = REPORT_NAMES[kind]
    if not isinstance(report, Mapping):
        return [f"{name} is not a JSON object"]
    problems: list[str] = []
    if report.get("passed") is not True:
        problems.append(f"{name} records passed {json.dumps(report.get('passed'))}, not true")
    checks = report.get("checks")
    if not (isinstance(checks, list) and checks and all(isinstance(c, Mapping) for c in checks)):
        problems.append(f"{name}'s checks must be a non-empty list of rows, got {type(checks).__name__} {checks!r:.40}: "
                        "its verdict stands on nothing")
    else:
        failed = [str(c.get("name")) for c in checks if c.get("severity") == "fail" and c.get("passed") is not True]
        if failed:
            problems.append(f"{name}'s own rows {failed} failed, which contradicts a passed verdict")
        try:
            required = expected_rows(kind)
        except Exception as exc:  # noqa: BLE001 - a checker the checkout cannot import is a refusal, not a traceback
            required = {}
            problems.append(f"the rows {kind} must hold cannot be read from the checkout ({exc})")
        for row, severity in required.items():
            if not any(c.get("name") == row and c.get("severity") == severity for c in checks):
                problems.append(f"{name} lacks the row {row} at severity {severity} (CDX-026: made under another "
                                "pointwise policy, or before it; copy it aside and rerun its checker, no AEDT)")
    binding = report.get("binding")
    if binding is None:
        return problems + [f"{name} carries no file binding (made before CDX-038); copy it aside, then rerun python -m "
                           f"{CHECKERS[kind]} on this job, no AEDT, before assembling it"]
    if not isinstance(binding, Mapping):
        return problems + [f"{name}'s binding is not an object ({binding!r})"]
    if binding.get("version") != BINDING_VERSION:
        return problems + [f"{name}'s binding version {binding.get('version')!r} is not {BINDING_VERSION!r}"]
    if binding.get("kind") != kind:
        problems.append(f"{name}'s binding kind {binding.get('kind')!r} is not {kind!r}")
    subject, reference, files = binding.get("subject_job"), binding.get("reference_job"), binding.get("files")
    if not isinstance(subject, str) or subject != _key(job):
        problems.append(f"{name} judged the subject job {subject!r}, not this job {_key(job)}")
    if not isinstance(reference, str) or not reference:
        return problems + [f"{name}'s binding names no reference job ({reference!r})"]
    if not (isinstance(files, Mapping) and all(isinstance(k, str) and isinstance(v, str) for k, v in files.items())):
        return problems + [f"{name}'s binding files must map paths to sha256 digests"]
    gated_tol = None
    if requirement.kind == kind:
        gated_tol = requirement.t_rel_tol
        if requirement.reference_job is not None and _key(reference) != requirement.reference_job:
            problems.append(f"{name} judged against the reference job {reference}, but the job requires {kind} "
                            f"against {requirement.reference_job}")
    try:
        expected = expected_tolerances(kind, gated_tol)
    except Exception as exc:  # noqa: BLE001 - a checker the checkout cannot import is a refusal, not a traceback
        expected = None
        problems.append(f"the tolerances {kind} must use cannot be read from the checkout ({exc})")
    if expected is not None and report.get("tolerances") != expected:
        problems.append(f"{name} used tolerances {report.get('tolerances')!r}, not {expected!r}")
    if not Path(reference).is_dir():
        return problems + [f"{name}'s reference job {reference} is not a directory"]
    sidecars = sorted(Path(job).glob("*.dataset.json"))
    try:
        sidecar = _read_object(inputs, sidecars[0]) if len(sidecars) == 1 else {}
    except (OSError, ValueError) as exc:
        sidecar = {}
        problems.append(f"{sidecars[0]} is unreadable ({exc})")
    poynting = report.get("poynting") if kind == "reexport_check" else None
    for path in evidence_files(kind, Path(job), Path(reference), sidecar,
                               Path(poynting) if isinstance(poynting, str) else None):
        if _key(path) not in files:
            problems.append(f"{path} is not covered by the acceptance ({name} binds no digest of it)")
    for path, recorded in files.items():
        whose = _whose(path, subject if isinstance(subject, str) else "", reference)
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
