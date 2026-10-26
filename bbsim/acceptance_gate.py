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

Stdlib only: the bridge's ``assemble_tree.sbatch`` imports this module from the checkout, and
``bbsim.config`` validates its ``[acceptance]`` table with ``validate_requirement_fields``, so this module must not
import the config or the checkers.
"""
from __future__ import annotations

import math
import os
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

KINDS = ("none", "reexport_check", "grid_acceptance")
REPORT_NAMES = {"reexport_check": "reexport_check.json", "grid_acceptance": "grid_acceptance.json"}
FIELDS = ("kind", "reference_job", "t_rel_tol")


class GateError(ValueError):
    """A job's acceptance cannot be determined or contradicts its own record: the job is refused."""


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
    config = manifest.get("config") if isinstance(manifest, Mapping) else None
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
