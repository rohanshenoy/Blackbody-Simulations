"""The check-row pieces the acceptance CLIs share (``bbsim.reexport_check``, ``bbsim.grid_acceptance``,
``bbsim.run_acceptance``).

A refusal type for input a row cannot accept; readers of a finite number, a nested key and a JSON object that refuse
as such; the relative difference of two numbers; the matcher of a per-key record list to the declared keys; a
fail-severity row whose detail states its rule, then its summary or its first ``MAX_LISTED`` failures; a guard that
turns any exception in a step into failing rows, never a traceback (CDX-030); and the failures / passed pattern of a
result that holds its rows in ``checks``.
"""
from __future__ import annotations

import json
import math
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any

from bbsim import dataset_v2
from bbsim.compare import Check
from bbsim.dataset_v2 import key_text

MAX_LISTED = 5             # a row's detail lists at most this many failures, then "and n more"

_MISSING = object()


class Refused(ValueError):
    """Input a check cannot accept: it fails the row that needed it."""


def at_key(ephi: int, key: tuple[float, float]) -> str:
    """``Ephi=0 at (phi, theta) = (45, 135) deg``: how a row names a polarization at an incidence key."""
    return f"Ephi={ephi} at {key_text(key)}"


def number(value: object, where: str) -> float:
    """``dataset_v2.finite_number``, refused as a row refuses (bool, string, null, NaN, infinity)."""
    try:
        return dataset_v2.finite_number(value, where)
    except dataset_v2.DatasetV2Error as exc:
        raise Refused(str(exc)) from exc


def get(mapping: object, *path: str, where: str, default: object = _MISSING) -> Any:
    """``dataset_v2.require``; a missing step returns ``default`` when one is given, else is refused."""
    try:
        return dataset_v2.require(mapping, *path, where=where)
    except dataset_v2.DatasetV2Error as exc:
        if default is _MISSING:
            raise Refused(str(exc)) from exc
        return default


def read_json(path: Path, what: str, *, read_text: Callable[[Path], str] | None = None) -> dict:
    """A JSON object from ``path``, read with ``read_text`` (a checker passes its ``InputCapture.read_text``, so the
    bytes parsed are the bytes its report binds; CDX-038); an unreadable file or anything but an object is refused."""
    try:
        data = json.loads(read_text(path) if read_text is not None else Path(path).read_text())
    except (OSError, ValueError) as exc:
        raise Refused(f"{what} {path} is unreadable: {exc}") from exc
    if not isinstance(data, dict):
        raise Refused(f"{what} {path} is not a JSON object")
    return data


def relative(a: float, b: float) -> float:
    """|a - b| / |b|, 0 when both are 0, infinity when only b is."""
    if a == b:
        return 0.0
    return abs(a - b) / abs(b) if b != 0.0 else math.inf


def phi_theta(record: object, where: str) -> tuple[float, float]:
    if not isinstance(record, Mapping):
        raise Refused(f"{where} is not an object")
    return number(record.get("phi_deg"), f"{where}.phi_deg"), number(record.get("theta_deg"), f"{where}.theta_deg")


def match_records(records: object, keys: Sequence[tuple[float, float]], where: str,
                  key_of: Callable[[object, str], tuple[float, float]] = phi_theta,
                  polarized: bool = False) -> tuple[dict[int, Any], list[str]]:
    """The records of a per-key list by declared key index, and what is wrong with the list: not a list, a record
    off the declared grid or unreadable, two records at one key, declared keys without a record. ``polarized``: each
    record holds both polarizations, so a repeated key repeats Ephi=0 and Ephi=1 alike."""
    if not isinstance(records, list):
        return {}, [f"{where} must be a list of per-key records, got {type(records).__name__}"]
    found: dict[int, Any] = {}
    problems: list[str] = []
    for i, record in enumerate(records):
        try:
            phi, theta = key_of(record, f"{where}[{i}]")
        except Refused as exc:
            problems.append(str(exc))
            continue
        k = dataset_v2.key_position(keys, phi, theta)
        if k is None:
            problems.append(f"{where} has a record at {key_text((phi, theta))}, outside the declared grid")
        elif k in found:
            repeated = ": Ephi=0 and Ephi=1 each appear twice" if polarized else ""
            problems.append(f"{where} has two records at {key_text(keys[k])}{repeated}")
        else:
            found[k] = record
    missing = [keys[k] for k in range(len(keys)) if k not in found]
    if missing:
        problems.append(f"{where} lacks {len(missing)} declared keys, the first at {key_text(missing[0])}")
    return found, problems


def strict_records(records: object, keys: Sequence[tuple[float, float]], where: str,
                   **kwargs: Any) -> dict[int, Any]:
    """``match_records``' records, refused with every problem it found when there is any."""
    found, problems = match_records(records, keys, where, **kwargs)
    if problems:
        raise Refused("; ".join(problems))
    return found


def row(name: str, rule: str, failures: list[str], summary: str, values: dict | None = None) -> Check:
    """A fail-severity row whose detail always states its rule, then the summary or the first failures."""
    if failures:
        listed = "; ".join(failures[:MAX_LISTED])
        more = f"; and {len(failures) - MAX_LISTED} more" if len(failures) > MAX_LISTED else ""
        detail = f"{rule}: {len(failures)} failure(s): {listed}{more}"
    else:
        detail = f"{rule}: {summary}"
    return Check(name, not failures, "fail", detail, {**(values or {}), "failures": failures})


def guard(names: Sequence[str], step: Callable[[], list[Check]]) -> list[Check]:
    """``step()``'s rows; any exception becomes one failing row per name, never a traceback (CDX-030)."""
    try:
        return step()
    except Exception as exc:  # noqa: BLE001 - every refusal becomes a failing row, never a traceback (CDX-030)
        text = str(exc) if isinstance(exc, Refused) else f"{type(exc).__name__}: {exc}"
        return [Check(name, False, "fail", f"refused: {text}") for name in names]


class RowVerdicts:
    """For a result dataclass holding its rows in ``checks``: the failing fail-severity rows, and whether it passed
    (at least one row, none failing)."""

    checks: list[Check]

    @property
    def failures(self) -> list[Check]:
        return [c for c in self.checks if c.severity == "fail" and not c.passed]

    @property
    def passed(self) -> bool:
        return bool(self.checks) and not self.failures
