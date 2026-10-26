"""The check-row pieces the acceptance CLIs share (``bbsim.reexport_check``, ``bbsim.grid_acceptance``).

A refusal type for input a row cannot accept; readers of a finite number, a nested key and a JSON object that refuse
as such; a fail-severity row whose detail states its rule, then its summary or its first ``MAX_LISTED`` failures; a
guard that turns any exception in a step into failing rows, never a traceback (CDX-030); and the failures / passed
pattern of a result that holds its rows in ``checks``.
"""
from __future__ import annotations

import json
from collections.abc import Callable, Sequence
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
