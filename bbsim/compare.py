"""Compare a candidate export dataset directory against a reference one.

Schema, keys, grids and IngoingPower must match; transmission must agree within
tolerance; field distributions are reported (warning only) because solver
versions differ.

Both dataset schemas are read (schema-2 design, section 2.6). A directory's schema comes from the dataset sidecar
beside it (``<stem>.dataset.json`` next to ``<stem>_Ephi=<n>``, as a job root and an assembled tree hold it), else from
its waveguide.csv column list (``dataset_schema``); anything else fails a check and is never guessed.

- Two 1.0 directories are compared as before.
- Two 2.0 directories are compared natively: OutgoingPower is the exit flux, T = OutgoingPower / (IngoingPower *
  InterceptFactor) per key (OR-034), and the W-scaled fields and InterceptFactor are compared as stored.
- A 2.0 directory against a 1.0 one is compared through the 2.0 side's legacy view (``legacy_view``: the |Re S|
  numerator back in OutgoingPower, the fields divided by W, the 1.0 labels, no InterceptFactor), and the result is
  marked ``mixed_schema``.

Two gates come before any value is compared. Every numeric column of every table must hold finite numbers (CDX-030):
a NaN, an infinity or a non-numeric entry fails the comparison and is named by table, column, row and key, instead of
leaving a perfect-match default in a metric. And every dark record of a 2.0 directory must have t_flux exactly 0 in
the sidecar and OutgoingPower exactly 0 on every row of its key (OR-038); a contradiction fails the comparison, naming
key, polarization and value, before any legacy view or dark masking, and T is never repaired. Exit fields of a 2.0
table are read by position, never by label (HF-038).
"""
from __future__ import annotations

import json
import math
import re
from collections.abc import Mapping
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from bbsim import dataset_v2
from bbsim.extract import ExitGrid
from bbsim.schema import (EXIT_FIELD_SLOTS, FAR_FIELD_COLUMNS, INTERCEPT_FACTOR, KEY_COLUMNS, WAVEGUIDE_COLUMNS,
                          WAVEGUIDE_COLUMNS_V2)

SCHEMA_1, SCHEMA_2 = "1.0", "2.0"
_COLUMNS_BY_SCHEMA = {SCHEMA_1: list(WAVEGUIDE_COLUMNS), SCHEMA_2: list(WAVEGUIDE_COLUMNS_V2)}
SIDECAR_SUFFIX = ".dataset.json"
_EPHI_DIR = re.compile(r"^(?P<stem>.+)_Ephi=(?P<ephi>[0-9]+)$")
_OUTGOING = WAVEGUIDE_COLUMNS_V2.index("OutgoingPower")      # the same position in 1.0
KEY_TOL_DEG = dataset_v2.KEY_TOL_DEG                          # sidecar per-key records against table keys
MAX_LISTED = 5


class ComparisonInputError(ValueError):
    """An input the comparison cannot read under either schema: it fails a check, it is never guessed."""


@dataclass(frozen=True)
class Tolerances:
    coord_m: float = 1e-9
    angle_deg: float = 1e-6
    ingoing_rel: float = 1e-9
    t_atol: float = 1e-3   # reference T at 45-degree incidence is 0.012-0.13; 0.01 would hide 80 % errors there
    t_rtol: float = 0.05
    t_noise_floor: float = 1e-6
    corr_warn: float = 0.95
    rho_atol: float = 0.05         # paired-polarization overlap Re rho; a relative sign flip moves it by about 2
    pair_dark_ratio: float = 1e-3  # a key whose weaker polarization carries less than this fraction of the
                                   # stronger one's T has no meaningful rho (cut-off polarization)
    pair_lit_ratio: float = 1e-2   # above this fraction both polarizations are lit; in between a key is not judged
                                   # for rho, so two solves straddling one threshold near a cutoff do not fail.
                                   # Legacy 500 GHz margins: lit keys >= 0.945, dark keys <= 1.9e-8.
    intercept_rel: float = 1e-12   # 2.0 InterceptFactor, candidate against reference: c_k comes from the declared
                                   # grid and is written at full precision (OR-034)


@dataclass
class Check:
    name: str
    passed: bool
    severity: str  # "fail" | "warn" | "info"
    detail: str
    values: dict = field(default_factory=dict)


@dataclass
class ComparisonResult:
    candidate: str
    reference: str
    checks: list[Check]
    tolerances: dict = field(default_factory=dict)
    mixed_schema: bool = False     # a 2.0 side compared through its legacy view against a 1.0 side
    schemas: dict = field(default_factory=dict)

    @property
    def passed(self) -> bool:
        return all(c.passed for c in self.checks if c.severity == "fail")

    def to_dict(self) -> dict:
        return {"candidate": self.candidate, "reference": self.reference, "passed": self.passed,
                "mixed_schema": self.mixed_schema, "schemas": dict(self.schemas),
                "tolerances": dict(self.tolerances), "checks": [asdict(c) for c in self.checks]}

    def format_table(self) -> str:
        tol = " ".join(f"{k}={v}" for k, v in self.tolerances.items())
        lines = [f"candidate: {self.candidate}", f"reference: {self.reference}"]
        if self.schemas:
            lines.append("schemas: " + ", ".join(f"{k} {v}" for k, v in self.schemas.items())
                         + (" (mixed: the 2.0 side is compared through its legacy view)" if self.mixed_schema else ""))
        lines += [f"tolerances: {tol}", ""]
        for c in self.checks:
            status = "PASS" if c.passed else ("FAIL" if c.severity == "fail" else "WARN")
            lines.append(f"{status:4} {c.name:26} {c.detail}")
        lines += ["", "RESULT: " + ("PASS" if self.passed else "FAIL")]
        return "\n".join(lines)


# --- Which schema a directory holds ----------------------------------------------------------------------------------

def sidecar_path(directory: Path) -> Path | None:
    """The dataset sidecar of a dataset directory, or None: ``<stem>.dataset.json`` beside ``<stem>_Ephi=<n>``, where
    the runner writes it in a job root and the assembly copies it into a tree."""
    d = Path(directory)
    match = _EPHI_DIR.match(d.name)
    if match is None:
        return None
    path = d.parent / f"{match['stem']}{SIDECAR_SUFFIX}"
    return path if path.is_file() else None


def _header(directory: Path) -> list[str]:
    path = Path(directory) / "waveguide.csv"
    try:
        return [str(c) for c in pd.read_csv(path, nrows=0).columns]
    except (OSError, ValueError) as exc:
        raise ComparisonInputError(f"{path} is unreadable: {exc}") from exc


def _read_sidecar(path: Path) -> dict:
    try:
        data = json.loads(Path(path).read_text())
    except (OSError, ValueError) as exc:
        raise ComparisonInputError(f"the sidecar {path} is unreadable: {exc}") from exc
    if not isinstance(data, dict):
        raise ComparisonInputError(f"the sidecar {path} is not a JSON object")
    return data


def _resolve(directory: Path) -> tuple[str, dict | None, Path | None]:
    """(schema, sidecar or None, sidecar path or None) of one dataset directory; see ``dataset_schema``."""
    d = Path(directory)
    header = _header(d)
    by_header = next((v for v, columns in _COLUMNS_BY_SCHEMA.items() if header == columns), None)
    path = sidecar_path(d)
    if path is None:
        if by_header is None:
            raise ComparisonInputError(
                f"{d / 'waveguide.csv'} carries neither the 1.0 nor the 2.0 column list and no dataset sidecar beside "
                f"{d} declares its schema: {header}")
        return by_header, None, None
    sidecar = _read_sidecar(path)
    declared = sidecar.get("schema_version")
    if not isinstance(declared, str) or declared not in _COLUMNS_BY_SCHEMA:
        raise ComparisonInputError(f"the sidecar {path} declares schema_version {declared!r}: expected "
                                   f"{SCHEMA_1!r} or {SCHEMA_2!r}")
    if by_header != declared:
        found = f"the {by_header} column list" if by_header else f"neither column list ({header})"
        raise ComparisonInputError(f"the sidecar {path} declares schema {declared} but {d / 'waveguide.csv'} "
                                   f"carries {found}")
    return declared, sidecar, path


def dataset_schema(directory: Path) -> str:
    """The schema of a dataset directory, "1.0" or "2.0": the ``schema_version`` of the dataset sidecar beside it when
    there is one, else what its waveguide.csv column list is (the 15 1.0 columns, or the 16 2.0 columns ending in
    InterceptFactor).

    Raises ``ComparisonInputError`` for anything else: an unknown version, a header that is neither list, or a
    sidecar that the header contradicts. The comparisons turn that into a failing check; the schema is never guessed.
    """
    return _resolve(directory)[0]


def _schema_check(name: str, problems: list[str], versions: dict[str, str], mixed: bool) -> Check:
    if problems:
        detail = "; ".join(problems)
    else:
        detail = f"candidate {versions['candidate']}, reference {versions['reference']}"
        if mixed:
            detail += ": mixed, the 2.0 side is compared through its legacy view"
        elif versions["candidate"] == SCHEMA_2:
            detail += ": native schema-2 comparison"
    return Check(name, not problems, "fail", detail,
                 {"schemas": dict(versions), "mixed_schema": mixed, "failures": list(problems)})


def _listing(name: str, problems: list[str], ok_detail: str, values: dict | None = None) -> Check:
    """A fail-severity check listing its first problems, or ``ok_detail`` when there are none."""
    if problems:
        more = f"; and {len(problems) - MAX_LISTED} more" if len(problems) > MAX_LISTED else ""
        detail = f"{len(problems)} problem(s): " + "; ".join(problems[:MAX_LISTED]) + more
    else:
        detail = ok_detail
    return Check(name, not problems, "fail", detail, {**(values or {}), "failures": list(problems)})


# --- CDX-030: every numeric value finite, before anything is compared -------------------------------------------------

def _key_text(key: tuple[float, float]) -> str:
    return f"(phi, theta) = ({key[0]:g}, {key[1]:g}) deg"


def _shown(value: object) -> str:
    if isinstance(value, str):
        return repr(value)
    try:
        return repr(float(value))                              # nan, inf, -inf, without numpy's type wrapper
    except (TypeError, ValueError):
        return repr(value)


def _key_at(df: pd.DataFrame, i: int) -> str:
    try:
        return f", at {_key_text((float(df[KEY_COLUMNS[0]].iloc[i]), float(df[KEY_COLUMNS[1]].iloc[i])))}"
    except (KeyError, IndexError, TypeError, ValueError):
        return ""


def _non_finite(where: str, df: pd.DataFrame) -> list[str]:
    """One message per numeric column (every column but Freq) holding a NaN, an infinity or a non-number."""
    problems = []
    for column in df.columns:
        if column == "Freq":
            continue
        numbers = pd.to_numeric(df[column], errors="coerce").to_numpy(dtype=float)
        bad = ~np.isfinite(numbers)
        if bad.any():
            i = int(np.flatnonzero(bad)[0])
            problems.append(f"{where} column {column}: {int(bad.sum())} non-finite or non-numeric value(s), the first "
                            f"at data row {i}{_key_at(df, i)}: {_shown(df[column].iloc[i])}")
    return problems


def _check_finite(name: str, tables: list[tuple[str, pd.DataFrame]]) -> Check:
    problems = [p for where, df in tables for p in _non_finite(where, df)]
    return _listing(name, problems, "every numeric value is finite")


# --- Schema 2: the sidecar's per-key records and the dark rule (OR-034, OR-038) ---------------------------------------

def _get(mapping: object, *path: str, where: str) -> Any:
    node = mapping
    for name in path:
        if not isinstance(node, Mapping) or name not in node:
            raise ComparisonInputError(f"{where} lacks {'.'.join(path)}")
        node = node[name]
    return node


def _number(value: object, where: str, positive: bool = False, nonnegative: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ComparisonInputError(f"{where} must be a number, got {value!r}")
    number = float(value)
    if not math.isfinite(number) or (positive and number <= 0.0) or (nonnegative and number < 0.0):
        kind = "positive " if positive else "non-negative " if nonnegative else ""
        raise ComparisonInputError(f"{where} must be a finite {kind}number, got {value!r}")
    return number


def _find_key(keys: list[tuple[float, float]], key: tuple[float, float]) -> int | None:
    for k, (phi, theta) in enumerate(keys):
        if abs(phi - key[0]) <= KEY_TOL_DEG and abs(theta - key[1]) <= KEY_TOL_DEG:
            return k
    return None


@dataclass(frozen=True)
class _V2Records:
    """A 2.0 sidecar's per-key records, by declared key: c_k, and per polarization t_flux, |Re S|, W and dark."""

    keys: list[tuple[float, float]]
    c_k: list[float]
    t_flux: dict[int, list[float]]
    abs_re_s: dict[int, list[float]]
    scale: dict[int, list[float]]
    dark: dict[int, list[bool]]


def _v2_records(sidecar: Mapping, where: str) -> _V2Records:
    """``transmittance.per_key`` and ``exit_field.scale.per_key``, checked: one record per key in each list, the same
    keys (within 1e-9 deg), every number finite, W positive, |Re S| non-negative, the dark flags booleans."""
    per_key = _get(sidecar, "transmittance", "per_key", where=where)
    scale = _get(sidecar, "exit_field", "scale", "per_key", where=where)
    if not isinstance(per_key, list) or not isinstance(scale, list) or not per_key:
        raise ComparisonInputError(f"{where}: transmittance.per_key and exit_field.scale.per_key must be non-empty "
                                   "lists")
    keys: list[tuple[float, float]] = []
    c_k: list[float] = []
    t_flux: dict[int, list[float]] = {0: [], 1: []}
    abs_re_s: dict[int, list[float]] = {0: [], 1: []}
    for i, rec in enumerate(per_key):
        at = f"{where} transmittance.per_key[{i}]"
        key = (_number(_get(rec, "phi_deg", where=at), f"{at}.phi_deg"),
               _number(_get(rec, "theta_deg", where=at), f"{at}.theta_deg"))
        if _find_key(keys, key) is not None:
            raise ComparisonInputError(f"{at}: a second record at {_key_text(key)}")
        keys.append(key)
        c_k.append(_number(_get(rec, "c_k", where=at), f"{at}.c_k", positive=True))
        pols = _get(rec, "polarizations", where=at)
        if not isinstance(pols, Mapping) or set(pols) != {"0", "1"}:
            raise ComparisonInputError(f"{at}.polarizations must hold exactly '0' and '1'")
        for e in (0, 1):
            t_flux[e].append(_number(_get(pols[str(e)], "t_flux", where=f"{at}.polarizations.{e}"),
                                     f"{at}.polarizations.{e}.t_flux"))
            abs_re_s[e].append(_number(_get(pols[str(e)], "outgoing_power_abs_re_s_w", where=f"{at}.polarizations.{e}"),
                                       f"{at}.polarizations.{e}.outgoing_power_abs_re_s_w", nonnegative=True))
    w: dict[int, list[float | None]] = {e: [None] * len(keys) for e in (0, 1)}
    dark: dict[int, list[bool | None]] = {e: [None] * len(keys) for e in (0, 1)}
    for i, rec in enumerate(scale):
        at = f"{where} exit_field.scale.per_key[{i}]"
        key = (_number(_get(rec, "phi_deg", where=at), f"{at}.phi_deg"),
               _number(_get(rec, "theta_deg", where=at), f"{at}.theta_deg"))
        k = _find_key(keys, key)
        if k is None:
            raise ComparisonInputError(f"{at}: {_key_text(key)} has no transmittance.per_key record")
        if w[0][k] is not None:
            raise ComparisonInputError(f"{at}: a second record at {_key_text(key)}")
        for e in (0, 1):
            w[e][k] = _number(_get(rec, f"w{e}", where=at), f"{at}.w{e}", positive=True)
            flag = _get(rec, f"dark{e}", where=at)
            if not isinstance(flag, bool):
                raise ComparisonInputError(f"{at}.dark{e} must be true or false, got {flag!r}")
            dark[e][k] = flag
    missing = [keys[k] for k in range(len(keys)) if w[0][k] is None]
    if missing:
        raise ComparisonInputError(f"{where}: exit_field.scale.per_key lacks {_key_text(missing[0])} "
                                   f"({len(missing)} key(s) missing)")
    return _V2Records(keys, c_k, t_flux, abs_re_s, {e: [float(v) for v in w[e]] for e in (0, 1)},
                      {e: [bool(v) for v in dark[e]] for e in (0, 1)})


def _row_index(df: pd.DataFrame, records: _V2Records, where: str) -> np.ndarray:
    """The declared key of every row: every row at a declared key, every declared key with rows."""
    phis = df[KEY_COLUMNS[0]].to_numpy(dtype=float)
    thetas = df[KEY_COLUMNS[1]].to_numpy(dtype=float)
    index = np.full(len(df), -1, dtype=int)
    for k, (phi, theta) in enumerate(records.keys):
        mask = (np.abs(phis - phi) <= KEY_TOL_DEG) & (np.abs(thetas - theta) <= KEY_TOL_DEG)
        if not mask.any():
            raise ComparisonInputError(f"{where}: the sidecar declares {_key_text((phi, theta))} but the table has no "
                                       "rows there")
        index[mask] = k
    unmatched = np.flatnonzero(index < 0)
    if unmatched.size:
        i = int(unmatched[0])
        raise ComparisonInputError(f"{where}: {unmatched.size} row(s) at keys the sidecar does not declare, the first "
                                   f"at data row {i}{_key_at(df, i)}")
    return index


def _ephi_of(df: pd.DataFrame, directory: Path | None, where: str) -> int:
    """The one polarization a table holds; a directory named ``..._Ephi=<n>`` must hold that one (BBRsim reads the
    polarization from the name)."""
    values = sorted({float(v) for v in df["Ephi"].to_numpy(dtype=float)})
    if len(values) != 1 or values[0] not in (0.0, 1.0):
        raise ComparisonInputError(f"{where}: the Ephi column holds {values}, not one polarization 0 or 1")
    ephi = int(values[0])
    match = _EPHI_DIR.match(Path(directory).name) if directory is not None else None
    if match is not None and int(match["ephi"]) != ephi:
        raise ComparisonInputError(f"{where}: the directory is named for Ephi={match['ephi']} but its table holds "
                                   f"Ephi={ephi}")
    return ephi


def _dark_contradictions(df: pd.DataFrame, records: _V2Records, index: np.ndarray, ephi: int,
                         where: str) -> list[str]:
    """OR-038: a dark record needs t_flux exactly 0 and OutgoingPower exactly 0 on every row of its key. One message
    per contradiction, naming key, polarization and value; nothing is repaired."""
    power = df.iloc[:, _OUTGOING].to_numpy(dtype=float)
    problems = []
    for k, key in enumerate(records.keys):
        if not records.dark[ephi][k]:
            continue
        at = f"{where} Ephi={ephi} at {_key_text(key)}"
        t = records.t_flux[ephi][k]
        if t != 0.0:
            problems.append(f"{at}: dark{ephi} is true but the sidecar's t_flux is {t!r}, not exactly 0")
        rows = np.flatnonzero(index == k)
        nonzero = np.flatnonzero(power[rows] != 0.0)
        if nonzero.size:
            j = int(nonzero[0])
            problems.append(f"{at}: dark{ephi} is true but OutgoingPower is {float(power[rows[j]])!r} on data row "
                            f"{int(rows[j])} (row {j + 1} of the key's {rows.size}; {nonzero.size} nonzero), not "
                            "exactly 0")
    return problems


def _v2_lattice(sidecar: Mapping, where: str) -> tuple[str, ExitGrid | None]:
    """The exit cross-section's shape and, for a rectangle, the declared lattice as an ``ExitGrid`` (mm) for
    ``dataset_v2.node_weights``; each axis's step must span min to max in count - 1 steps (HF-040)."""
    shape = _get(sidecar, "exit_field", "cross_section", "shape", where=where)
    if shape not in dataset_v2.QUADRATURE_BY_SHAPE:
        raise ComparisonInputError(f"{where}: exit_field.cross_section.shape {shape!r} is not one of "
                                   f"{sorted(dataset_v2.QUADRATURE_BY_SHAPE)}")
    quadrature = _get(sidecar, "exit_field", "quadrature", where=where)
    if quadrature != dataset_v2.QUADRATURE_BY_SHAPE[shape]:
        raise ComparisonInputError(f"{where}: exit_field.quadrature {quadrature!r} is not the {shape}'s "
                                   f"{dataset_v2.QUADRATURE_BY_SHAPE[shape]!r}")
    if shape == "disc":
        return shape, None
    lo, hi, steps, counts = [0.0], [0.0], [0.0], []
    for axis in ("y_e", "z_e"):
        at = f"{where} exit_field.grid.{axis}"
        rec = _get(sidecar, "exit_field", "grid", axis, where=where)
        lo_m = _number(_get(rec, "min", where=at), f"{at}.min")
        hi_m = _number(_get(rec, "max", where=at), f"{at}.max")
        step_m = _number(_get(rec, "step", where=at), f"{at}.step", positive=True)
        count = _get(rec, "count", where=at)
        if isinstance(count, bool) or not isinstance(count, int) or count < 2:
            raise ComparisonInputError(f"{at}.count must be an integer of at least 2, got {count!r}")
        if abs((hi_m - lo_m) / (count - 1) - step_m) > dataset_v2.NODE_TOL * step_m:
            raise ComparisonInputError(f"{at}: the step {step_m!r} m does not span {lo_m!r} to {hi_m!r} m in "
                                       f"{count - 1} steps")
        lo.append(lo_m * 1e3)
        hi.append(hi_m * 1e3)
        steps.append(step_m * 1e3)
        counts.append(count)
    grid = ExitGrid(lo_mm=tuple(lo), hi_mm=tuple(hi), resolution_mm=tuple(steps), step_mm=tuple(steps))
    if list(grid.counts[1:]) != counts:
        raise ComparisonInputError(f"{where}: the declared counts {counts} are not the lattice's "
                                   f"{list(grid.counts[1:])}")
    return shape, grid


@dataclass
class _V2Table:
    """One 2.0 directory, checked: its sidecar's records, the key of every row, its polarization."""

    sidecar: dict
    records: _V2Records
    index: np.ndarray
    ephi: int


def _v2_table(df: pd.DataFrame, sidecar: dict | None, path: Path | None, directory: Path | None,
              where: str) -> _V2Table:
    if sidecar is None:
        raise ComparisonInputError(f"{where}: no dataset sidecar beside {directory}; a 2.0 table is read only with "
                                   "its sidecar (dark flags, W, |Re S|, the declared grid)")
    records = _v2_records(sidecar, f"the sidecar {path}" if path else "the sidecar")
    return _V2Table(sidecar, records, _row_index(df, records, where), _ephi_of(df, directory, where))


def legacy_view(waveguide: pd.DataFrame, sidecar: dict, ephi: int) -> pd.DataFrame:
    """The schema-1.0 table a 2.0 waveguide table was made from, for mixed-version comparisons only (design 2.6).

    OutgoingPower becomes the key's ``transmittance.per_key[].polarizations[ephi].outgoing_power_abs_re_s_w`` (the
    |Re S| numerator of 1.0), IngoingPower stays A I, the six exit-field values (by position) are divided by the
    key's W (``exit_field.scale`` w0 / w1; 1 on a dark record), the 1.0 labels go back over the unchanged positions
    and InterceptFactor is dropped. Raises ``ComparisonInputError`` before building anything on a table that is not
    the 2.0 column list, holds a non-finite value or another polarization, does not match the sidecar's keys, or
    holds a dark record whose t_flux or OutgoingPower is not exactly zero (OR-038).
    """
    if isinstance(ephi, bool) or ephi not in (0, 1):
        raise ComparisonInputError(f"ephi must be 0 or 1, got {ephi!r}")
    if list(waveguide.columns) != list(WAVEGUIDE_COLUMNS_V2):
        raise ComparisonInputError(f"a legacy view needs the 2.0 columns {WAVEGUIDE_COLUMNS_V2}, got "
                                   f"{list(waveguide.columns)}")
    where = "the 2.0 waveguide table"
    problems = _non_finite(where, waveguide)
    if problems:
        raise ComparisonInputError("; ".join(problems))
    table = _v2_table(waveguide, sidecar, None, None, where)
    if table.ephi != ephi:
        raise ComparisonInputError(f"{where} holds Ephi={table.ephi}, not {ephi}")
    problems = _dark_contradictions(waveguide, table.records, table.index, ephi, where)
    if problems:
        raise ComparisonInputError("; ".join(problems))
    w = np.asarray(table.records.scale[ephi], dtype=float)[table.index]
    view = waveguide.iloc[:, :len(WAVEGUIDE_COLUMNS)].copy()
    view.columns = list(WAVEGUIDE_COLUMNS)
    view["OutgoingPower"] = np.asarray(table.records.abs_re_s[ephi], dtype=float)[table.index]
    fields = waveguide.iloc[:, EXIT_FIELD_SLOTS].to_numpy(dtype=float) / w[:, None]
    for j, position in enumerate(range(EXIT_FIELD_SLOTS.start, EXIT_FIELD_SLOTS.stop)):
        view[view.columns[position]] = fields[:, j]
    return view


# --- Per-table checks -------------------------------------------------------------------------------------------------

def _keys(df: pd.DataFrame) -> set[tuple[float, float]]:
    return {(float(a), float(b)) for a, b in df[KEY_COLUMNS].drop_duplicates().itertuples(index=False, name=None)}


def _check_columns(table: str, cand: pd.DataFrame, ref: pd.DataFrame, expected: list[str]) -> Check:
    ok = list(cand.columns) == expected and list(ref.columns) == expected
    detail = "column order matches Geant4 contract" if ok else \
        f"candidate {list(cand.columns)} / reference {list(ref.columns)} / expected {expected}"
    return Check(f"{table}.columns", ok, "fail", detail)


def _check_constant(table: str, col: str, cand: pd.DataFrame, ref: pd.DataFrame) -> Check:
    cv, rv = sorted(cand[col].astype(str).unique()), sorted(ref[col].astype(str).unique())
    ok = cv == rv and len(cv) == 1
    return Check(f"{table}.{col}", ok, "fail", f"{cv} vs {rv}", {"candidate": cv, "reference": rv})


def _check_keys(table: str, cand: pd.DataFrame, ref: pd.DataFrame) -> Check:
    ck, rk = _keys(cand), _keys(ref)
    ok = ck == rk
    detail = f"{len(ck)} incident angles" if ok else f"missing {sorted(rk - ck)} extra {sorted(ck - rk)}"
    return Check(f"{table}.keys", ok, "fail", detail, {"missing": sorted(rk - ck), "extra": sorted(ck - rk)})


def _check_ingoing(cand: pd.DataFrame, ref: pd.DataFrame, tol: Tolerances) -> Check:
    c, r = float(cand["IngoingPower"].iloc[0]), float(ref["IngoingPower"].iloc[0])
    rel = abs(c - r) / abs(r) if r else float("inf")
    return Check("waveguide.ingoing_power", rel <= tol.ingoing_rel, "fail",
                 f"{c:.9e} vs {r:.9e} (rel {rel:.2e})", {"candidate": c, "reference": r, "rel": rel})


def _transmission(df: pd.DataFrame, intercept: bool = False) -> dict[tuple[float, float], float]:
    """T per key from the key's first row: OutgoingPower / IngoingPower (1.0), or with ``intercept`` OutgoingPower /
    (IngoingPower * InterceptFactor) (2.0, OR-034)."""
    first = df.drop_duplicates(KEY_COLUMNS)
    if intercept:
        cols = KEY_COLUMNS + ["OutgoingPower", "IngoingPower", INTERCEPT_FACTOR]
        return {(float(a), float(b)): float(o) / (float(i) * float(c))
                for a, b, o, i, c in first[cols].itertuples(index=False, name=None)}
    return {(float(a), float(b)): float(o) / float(i)
            for a, b, o, i in first[KEY_COLUMNS + ["OutgoingPower", "IngoingPower"]].itertuples(index=False, name=None)}


def _check_transmission(cand: pd.DataFrame, ref: pd.DataFrame, tol: Tolerances, intercept: bool = False) -> Check:
    tc, tr = _transmission(cand, intercept), _transmission(ref, intercept)
    rows, failures, noise = [], [], []
    worst_rel, worst_key = 0.0, None
    for key in sorted(set(tc) & set(tr)):
        c, r = tc[key], tr[key]
        if c < tol.t_noise_floor and r < tol.t_noise_floor:
            noise.append(key)
            rows.append({"key": list(key), "candidate": c, "reference": r, "status": "noise"})
            continue
        rel = abs(c - r) / abs(r) if r else float("inf")
        if rel > worst_rel:
            worst_rel, worst_key = rel, key
        ok = abs(c - r) <= max(tol.t_atol, tol.t_rtol * abs(r))
        rows.append({"key": list(key), "candidate": c, "reference": r, "rel_diff": rel,
                     "status": "pass" if ok else "fail"})
        if not ok:
            failures.append(key)
    n_pass = len(rows) - len(failures) - len(noise)
    detail = f"{n_pass} pass, {len(noise)} below noise floor, {len(failures)} fail"
    if failures:
        detail += f": {failures}"
    if worst_key is not None:
        detail += f"; max rel diff {worst_rel:.3f} at {worst_key}"
    return Check("waveguide.transmission", not failures, "fail", detail,
                 {"rows": rows, "max_rel_diff": worst_rel, "max_rel_diff_key": list(worst_key) if worst_key else None})


def _check_intercept_factor(cand: pd.DataFrame, ref: pd.DataFrame, tol: Tolerances) -> Check:
    """2.0: InterceptFactor constant within every key on each side, and the candidate's equal to the reference's."""
    failures = []
    first = {}
    for label, df in (("candidate", cand), ("reference", ref)):
        values = df[INTERCEPT_FACTOR].to_numpy(dtype=float)
        first[label] = {}
        for key in sorted(_keys(df)):
            v = values[((df[KEY_COLUMNS[0]] == key[0]) & (df[KEY_COLUMNS[1]] == key[1])).to_numpy()]
            if np.any(v != v[0]):
                failures.append(f"{label} {_key_text(key)}: InterceptFactor varies within the key "
                                f"({float(v.min())!r} to {float(v.max())!r})")
            first[label][key] = float(v[0])
    worst = 0.0
    for key in sorted(set(first["candidate"]) & set(first["reference"])):
        c, r = first["candidate"][key], first["reference"][key]
        rel = abs(c - r) / abs(r)
        worst = max(worst, rel)
        if not rel <= tol.intercept_rel:
            failures.append(f"{_key_text(key)}: candidate {c!r} against reference {r!r} (relative {rel:.3g})")
    compared = len(set(first["candidate"]) & set(first["reference"]))
    return _listing("waveguide.intercept_factor", failures,
                    f"{compared} keys equal (max relative diff {worst:.2e}), each constant within its key",
                    {"max_rel_diff": worst})


def _magnitude(df: pd.DataFrame, pairs: list[tuple[str, str]]) -> np.ndarray:
    total = np.zeros(len(df))
    for re_col, im_col in pairs:
        total += df[re_col].to_numpy() ** 2 + df[im_col].to_numpy() ** 2
    return np.sqrt(total)


# |E| from (real, imaginary) column pairs. The 1.0 waveguide pairs are taken by label as before: |E| sums all six
# values whatever their order. A 2.0 table's pairs are its columns at the stored positions (9, 10), (11, 12), (13, 14).
WAVEGUIDE_PAIRS_V1 = [("Ex_real", "Ex_imag"), ("Ey_real", "Ey_imag"), ("Ez_real", "Ez_imag")]
FAR_FIELD_PAIRS = [("rEtheta_real", "rEtheta_imag"), ("rEphi_real", "rEphi_imag")]


def _slot_pairs(df: pd.DataFrame) -> list[tuple[str, str]]:
    names = list(df.columns)
    return [(names[i], names[i + 1]) for i in range(EXIT_FIELD_SLOTS.start, EXIT_FIELD_SLOTS.stop, 2)]


def _check_points_and_fields(table: str, cand: pd.DataFrame, ref: pd.DataFrame, coord_cols: list[str],
                             pairs: list[tuple[str, str]], coord_tol: float, tol: Tolerances) -> list[Check]:
    point_fail: list[tuple[tuple[float, float], str]] = []
    worst_coord, min_corr, max_nrmsd = 0.0, 1.0, 0.0
    compared = 0
    for key in sorted(_keys(cand) & _keys(ref)):
        mc = (cand[KEY_COLUMNS[0]] == key[0]) & (cand[KEY_COLUMNS[1]] == key[1])
        mr = (ref[KEY_COLUMNS[0]] == key[0]) & (ref[KEY_COLUMNS[1]] == key[1])
        c = cand[mc].sort_values(coord_cols).reset_index(drop=True)
        r = ref[mr].sort_values(coord_cols).reset_index(drop=True)
        if len(c) != len(r):
            point_fail.append((key, f"{len(c)} vs {len(r)} points"))
            continue
        dc = float(np.max(np.abs(c[coord_cols].to_numpy(dtype=float) - r[coord_cols].to_numpy(dtype=float))))
        worst_coord = max(worst_coord, dc)
        if dc > coord_tol:
            point_fail.append((key, f"coordinate mismatch {dc:.3e}"))
            continue
        compared += 1
        mag_c, mag_r = _magnitude(c, pairs), _magnitude(r, pairs)
        scale = float(np.max(mag_r)) or 1.0
        max_nrmsd = max(max_nrmsd, float(np.sqrt(np.mean((mag_c - mag_r) ** 2)) / scale))
        if np.std(mag_c) > 0 and np.std(mag_r) > 0:
            min_corr = min(min_corr, float(np.corrcoef(mag_c, mag_r)[0, 1]))
    points = Check(f"{table}.points", not point_fail, "fail",
                   f"grids identical (worst coordinate diff {worst_coord:.2e})" if not point_fail
                   else f"{len(point_fail)} angles differ: {point_fail[:5]}",
                   {"failures": [(list(k), d) for k, d in point_fail]})
    if compared == 0:
        fields = Check(f"{table}.fields", False, "warn", "no angles with matching grids; fields not compared",
                       {"compared_angles": 0})
    else:
        fields = Check(f"{table}.fields", min_corr >= tol.corr_warn, "warn",
                       f"{compared} angles: min |E| correlation {min_corr:.4f}, max normalized RMS diff {max_nrmsd:.4f}",
                       {"compared_angles": compared, "min_corr": min_corr, "max_nrmsd": max_nrmsd})
    return [points, fields]


def _waveguide_values(cand: pd.DataFrame, ref: pd.DataFrame, tol: Tolerances, native_v2: bool) -> list[Check]:
    """The waveguide checks after the gates: 1.0 (or legacy-view) tables as before; native 2.0 tables with T divided by
    InterceptFactor, InterceptFactor compared and the fields taken by position."""
    checks = [_check_constant("waveguide", "Freq", cand, ref), _check_constant("waveguide", "Ephi", cand, ref),
              _check_keys("waveguide", cand, ref), _check_ingoing(cand, ref, tol),
              _check_transmission(cand, ref, tol, intercept=native_v2)]
    if native_v2:
        checks.append(_check_intercept_factor(cand, ref, tol))
    pairs = _slot_pairs(cand) if native_v2 else WAVEGUIDE_PAIRS_V1
    return checks + _check_points_and_fields("waveguide", cand, ref, ["X", "Y", "Z"], pairs, tol.coord_m, tol)


def _compare_waveguide(tables: dict[str, pd.DataFrame], resolved: dict[str, tuple], dirs: dict[str, Path],
                       tol: Tolerances) -> list[Check]:
    """Columns, the finite gate, for 2.0 sides the sidecar records and the dark rule, then the value checks: native
    when both sides are 2.0, through the 2.0 side's legacy view when one is (mixed)."""
    versions = {label: resolved[label][0] for label in tables}
    expected = {label: _COLUMNS_BY_SCHEMA[v] for label, v in versions.items()}
    cand, ref = tables["candidate"], tables["reference"]
    if versions["candidate"] == versions["reference"]:
        columns = _check_columns("waveguide", cand, ref, expected["candidate"])
    else:
        ok = all(list(tables[label].columns) == expected[label] for label in tables)
        columns = Check("waveguide.columns", ok, "fail",
                        f"candidate {versions['candidate']} and reference {versions['reference']} column lists" if ok
                        else f"candidate {list(cand.columns)} / reference {list(ref.columns)}")
    checks = [columns]
    if not columns.passed:
        return checks
    finite = _check_finite("waveguide.finite", [(f"{label} waveguide.csv", df) for label, df in tables.items()])
    checks.append(finite)
    if not finite.passed:
        return checks
    v2_sides = [label for label in tables if versions[label] == SCHEMA_2]
    if not v2_sides:
        return checks + _waveguide_values(cand, ref, tol, native_v2=False)
    v2: dict[str, _V2Table] = {}
    problems = []
    for label in v2_sides:
        _, sidecar, path = resolved[label]
        try:
            v2[label] = _v2_table(tables[label], sidecar, path, dirs[label], f"{label} waveguide.csv")
        except ComparisonInputError as exc:
            problems.append(f"{label}: {exc}")
    checks.append(_listing("dataset.sidecar", problems,
                           "; ".join(f"{label}: {len(v2[label].records.keys)} keys, Ephi={v2[label].ephi}"
                                     for label in v2_sides) if not problems else ""))
    if problems:
        return checks
    dark = [p for label in v2_sides for p in _dark_contradictions(tables[label], v2[label].records, v2[label].index,
                                                                   v2[label].ephi, f"{label} waveguide.csv")]
    n_dark = sum(sum(v2[label].records.dark[v2[label].ephi]) for label in v2_sides)
    checks.append(_listing("waveguide.dark_records", dark,
                           f"{n_dark} dark record(s), each with t_flux and OutgoingPower exactly 0 (OR-038)"))
    if dark:
        return checks
    if len(v2_sides) == 2:
        return checks + _waveguide_values(cand, ref, tol, native_v2=True)
    label = v2_sides[0]
    view = dict(tables)
    view[label] = legacy_view(tables[label], v2[label].sidecar, v2[label].ephi)
    return checks + _waveguide_values(view["candidate"], view["reference"], tol, native_v2=False)


def compare_datasets(candidate_dir: Path, reference_dir: Path, tol: Tolerances = Tolerances()) -> ComparisonResult:
    """Compare one dataset directory with another; either may be schema 1.0 or 2.0 (see the module docstring)."""
    candidate_dir, reference_dir = Path(candidate_dir), Path(reference_dir)
    dirs = {"candidate": candidate_dir, "reference": reference_dir}

    waveguides = {label: pd.read_csv(d / "waveguide.csv") for label, d in dirs.items()}
    resolved, problems = {}, []
    for label, d in dirs.items():
        try:
            resolved[label] = _resolve(d)
        except ComparisonInputError as exc:
            problems.append(f"{label}: {exc}")
    versions = {label: r[0] for label, r in resolved.items()}
    mixed = not problems and versions["candidate"] != versions["reference"]
    checks = [_schema_check("dataset.schema", problems, versions, mixed)]
    if problems:
        expected = _COLUMNS_BY_SCHEMA[next(iter(versions.values()), SCHEMA_1)]
        checks.append(_check_columns("waveguide", waveguides["candidate"], waveguides["reference"], expected))
    else:
        checks += _compare_waveguide(waveguides, resolved, dirs, tol)

    cf, rf = pd.read_csv(candidate_dir / "far_field.csv"), pd.read_csv(reference_dir / "far_field.csv")
    cols = _check_columns("far_field", cf, rf, FAR_FIELD_COLUMNS)
    checks.append(cols)
    if cols.passed:
        finite = _check_finite("far_field.finite", [("candidate far_field.csv", cf), ("reference far_field.csv", rf)])
        checks.append(finite)
        if finite.passed:
            checks += [_check_constant("far_field", "Freq", cf, rf), _check_constant("far_field", "Ephi", cf, rf),
                       _check_keys("far_field", cf, rf)]
            checks += _check_points_and_fields("far_field", cf, rf, ["Theta", "Phi"], FAR_FIELD_PAIRS, tol.angle_deg,
                                               tol)

    return ComparisonResult(str(candidate_dir), str(reference_dir), checks, asdict(tol), mixed_schema=mixed,
                            schemas=versions)


EXIT_VALUES = ["Ex_real", "Ey_real", "Ez_real", "Ex_imag", "Ey_imag", "Ez_imag"]
FAR_VALUES = ["rEtheta_real", "rEtheta_imag", "rEphi_real", "rEphi_imag"]


def _overlap(a: np.ndarray, b: np.ndarray, weight: np.ndarray | None = None) -> float | None:
    """Re of sum(a . conj b) / sqrt(sum|a|^2 sum|b|^2) from stored real/imaginary parts.

    Equal to the sum of products of the stored reals however they pair into complex components, so the
    exit-field column order (HF-038) does not change it."""
    w = np.ones(len(a)) if weight is None else weight
    p0, p1 = float(np.sum(w[:, None] * a * a)), float(np.sum(w[:, None] * b * b))
    if p0 <= 0.0 or p1 <= 0.0:
        return None
    return float(np.sum(w[:, None] * a * b) / np.sqrt(p0 * p1))


def _paired_keys(tables, tol: Tolerances):
    """Per key common to the two polarizations, its (exit, exit, far, far) rows, and the problems: keys present in
    one polarization only, and keys whose rows do not pair by order, as BBRsim pairs them."""
    (w0, f0), (w1, f1) = tables
    paired, problems = {}, []
    keys0, keys1 = _keys(w0), _keys(w1)
    if keys0 != keys1:
        problems.append(f"incidence keys differ between the polarizations: {sorted(keys0 ^ keys1)}")
    for key in sorted(keys0 & keys1):
        def sel(df):
            return df[(df[KEY_COLUMNS[0]] == key[0]) & (df[KEY_COLUMNS[1]] == key[1])].reset_index(drop=True)
        a, b, fa, fb = sel(w0), sel(w1), sel(f0), sel(f1)
        if len(a) != len(b) or not np.allclose(a[["X", "Y", "Z"]].to_numpy(float), b[["X", "Y", "Z"]].to_numpy(float),
                                               rtol=0.0, atol=tol.coord_m):
            problems.append(f"{key}: exit points do not pair by row")
            continue
        if len(fa) != len(fb) or not np.allclose(fa[["Theta", "Phi"]].to_numpy(float), fb[["Theta", "Phi"]].to_numpy(float),
                                                 rtol=0.0, atol=tol.angle_deg):
            problems.append(f"{key}: far-field grids do not pair by row")
            continue
        paired[key] = (a, b, fa, fb)
    return paired, problems


def _state(t: list[float], tol: Tolerances) -> str:
    return ("dark" if max(t) < tol.t_noise_floor or min(t) < tol.pair_dark_ratio * max(t)
            else "lit" if min(t) > tol.pair_lit_ratio * max(t) else "between")


def _pair_stats(tables, tol: Tolerances):
    """Per key of one 1.0 (Ephi=0, Ephi=1) pair: T0, T1, exit and far-field Re rho, and whether the key is dark.

    Rows pair by order within a key, as BBRsim pairs them; a pairing BBRsim would refuse is a problem."""
    paired, problems = _paired_keys(tables, tol)
    stats = {}
    for key, (a, b, fa, fb) in paired.items():
        t = [float(df["OutgoingPower"].iloc[0]) / float(df["IngoingPower"].iloc[0]) for df in (a, b)]
        sin_t = np.sin(np.radians(fa["Theta"].to_numpy(float)))
        stats[key] = {
            "T0": t[0], "T1": t[1],
            "state": _state(t, tol),
            "rho_exit": _overlap(a[EXIT_VALUES].to_numpy(float), b[EXIT_VALUES].to_numpy(float)),
            "rho_far": _overlap(fa[FAR_VALUES].to_numpy(float), fb[FAR_VALUES].to_numpy(float), sin_t),
        }
    return stats, problems


def _pair_stats_v2(tables, v2: list[_V2Table], lattice: tuple[str, ExitGrid | None], tol: Tolerances):
    """Per key of one native 2.0 pair (schema-2 design 2.6, HF-042): T_p = OutgoingPower / (IngoingPower *
    InterceptFactor); exit Re rho on E_t = (Ex, Ey), read by position, with the dimensionless quadrature weights of the
    declared lattice (trapezoid by each point's node for a rectangle, 1 for a disc); dark masking: a key where either
    basis is dark has rho 0 for the exit and the far field, its stored fields unread, and is not judged for rho."""
    paired, problems = _paired_keys(tables, tol)
    shape, grid = lattice
    stats = {}
    for key, (a, b, fa, fb) in paired.items():
        t = [float(df["OutgoingPower"].iloc[0])
             / (float(df["IngoingPower"].iloc[0]) * float(df[INTERCEPT_FACTOR].iloc[0])) for df in (a, b)]
        dark = []
        for table in v2:
            k = _find_key(table.records.keys, key)
            dark.append(bool(k is not None and table.records.dark[table.ephi][k]))
        masked = any(dark)
        try:
            weights = (dataset_v2.node_weights(a["Y"].to_numpy(float), a["Z"].to_numpy(float), grid)
                       if shape == "rectangle" else np.ones(len(a)))
            rho_exit = dataset_v2.rho(a.iloc[:, EXIT_FIELD_SLOTS].to_numpy(float),
                                      b.iloc[:, EXIT_FIELD_SLOTS].to_numpy(float), weights, dark[0], dark[1])
        except dataset_v2.DatasetV2Error as exc:
            problems.append(f"{key}: {exc}")
            continue
        sin_t = np.sin(np.radians(fa["Theta"].to_numpy(float)))
        stats[key] = {
            "T0": t[0], "T1": t[1],
            "state": "dark" if masked else _state(t, tol),
            "rho_exit": rho_exit,
            "rho_far": (0.0 if masked
                        else _overlap(fa[FAR_VALUES].to_numpy(float), fb[FAR_VALUES].to_numpy(float), sin_t)),
            "masked": masked,
        }
    return stats, problems


def _check_pair_ephi(candidate, reference) -> Check:
    wrong = []
    for label, tables in (("candidate", candidate), ("reference", reference)):
        for expected, (wg, ff) in enumerate(tables):
            found = sorted({int(v) for v in wg["Ephi"].unique()} | {int(v) for v in ff["Ephi"].unique()})
            if found != [expected]:
                wrong.append(f"{label} directory {expected + 1} holds Ephi {found}, expected [{expected}]")
    return Check("pair.ephi", not wrong, "fail", "Ephi=0 then Ephi=1 on both sides" if not wrong else "; ".join(wrong))


def _check_rho(name: str, field: str, cs: dict, rs: dict, tol: Tolerances) -> Check:
    rows, failures, dark = [], [], 0
    for key in sorted(set(cs) & set(rs)):
        c, r = cs[key], rs[key]
        if {c["state"], r["state"]} == {"dark", "lit"}:
            # A polarization cut off on one side and lit on the other is a difference, not a key to skip.
            rows.append({"key": list(key), "candidate_state": c["state"], "reference_state": r["state"], "status": "fail"})
            failures.append(key)
            continue
        if c["state"] != "lit" or r["state"] != "lit" or c[field] is None or r[field] is None:
            dark += 1
            continue
        ok = abs(c[field] - r[field]) <= tol.rho_atol
        rows.append({"key": list(key), "candidate": c[field], "reference": r[field], "status": "pass" if ok else "fail"})
        if not ok:
            failures.append(key)
    detail = f"{len(rows) - len(failures)} pass, {dark} dark or near-dark (not judged), {len(failures)} fail"
    if failures:
        detail += f": {failures}"
    values = {"rows": rows}
    masked = [list(k) for k in sorted(set(cs) | set(rs)) if cs.get(k, {}).get("masked") or rs.get(k, {}).get("masked")]
    if masked:
        detail += f"; {len(masked)} with a dark basis (rho = 0, masked)"
        values["masked"] = masked
    return Check(name, not failures, "fail", detail, values)


def _check_mixed(cs: dict, rs: dict, tol: Tolerances) -> Check:
    rows, failures = [], []
    for key in sorted(set(cs) & set(rs)):
        pred = {}
        for label, s in (("candidate", cs[key]), ("reference", rs[key])):
            cross = np.sqrt(s["T0"] * s["T1"]) * (s["rho_exit"] or 0.0)
            pred[label] = (0.5 * (s["T0"] + s["T1"]) + cross, 0.5 * (s["T0"] + s["T1"]) - cross)
        for sign, i in (("plus", 0), ("minus", 1)):
            c, r = pred["candidate"][i], pred["reference"][i]
            if c < tol.t_noise_floor and r < tol.t_noise_floor:
                continue
            ok = abs(c - r) <= max(tol.t_atol, tol.t_rtol * abs(r))
            rows.append({"key": list(key), "state": sign, "candidate": c, "reference": r, "status": "pass" if ok else "fail"})
            if not ok:
                failures.append((key, sign))
    detail = f"{len(rows) - len(failures)} pass, {len(failures)} fail"
    if failures:
        detail += f": {failures}"
    return Check("pair.mixed_transmission", not failures, "fail", detail, {"rows": rows})


def compare_polarization_pair(candidate: tuple[Path, Path], reference: tuple[Path, Path],
                              tol: Tolerances = Tolerances()) -> ComparisonResult:
    """Compare what coherent transport uses and the per-directory comparison cannot see (CDX-005).

    Each argument is (Ephi=0 directory, Ephi=1 directory). The checks are: the exit-field and far-field
    overlap Re rho of the two polarizations (BBRsim's polarization cross term; the far field weighted by
    sin(Theta)), and the transmission of equal-amplitude mixed states, T = (T0+T1)/2 +- sqrt(T0 T1) Re rho. A
    relative sign or phase change between the polarizations fails them; a phase common to both does not.
    Keys where one polarization is cut off are skipped for rho. At keys where rho is near 0 (normal
    incidence, E along one axis) a relative sign is invisible, and harmless, since equal-plus and equal-minus
    transmissions are then equal; the check sees the phase wherever it matters to BBRsim. T is the producer's
    OutgoingPower / IngoingPower, not BBRsim's load-time value divided by the key's top eigenvalue above 1.

    Both directories of a side share one schema. Two 2.0 sides are compared natively (HF-042): T is OutgoingPower /
    (IngoingPower * InterceptFactor), uncapped; exit rho is taken on E_t with the quadrature weights of the sidecar's
    declared lattice; a key with a dark basis on either polarization has rho 0 (dark masking) and is not judged for
    rho. A 2.0 side against a 1.0 side is compared through the 2.0 side's legacy views, and the result is marked
    ``mixed_schema``. Before anything is compared every numeric value must be finite (CDX-030), and every dark record
    of a 2.0 side must carry t_flux and OutgoingPower exactly 0 (OR-038).
    """
    sides = {"candidate": tuple(Path(d) for d in candidate), "reference": tuple(Path(d) for d in reference)}
    names = (f"{candidate[0]} + {candidate[1]}", f"{reference[0]} + {reference[1]}")
    resolved: dict[str, list[tuple]] = {}
    versions: dict[str, str] = {}
    problems: list[str] = []
    for label, dirs in sides.items():
        found = []
        for i, d in enumerate(dirs):
            try:
                found.append(_resolve(d))
            except ComparisonInputError as exc:
                problems.append(f"{label} directory {i + 1}: {exc}")
        if len(found) == 2 and found[0][0] != found[1][0]:
            problems.append(f"{label}: directory 1 is schema {found[0][0]} and directory 2 schema {found[1][0]}; the "
                            "two polarizations of one dataset share one schema")
        elif len(found) == 2:
            resolved[label], versions[label] = found, found[0][0]
    mixed = not problems and versions["candidate"] != versions["reference"]
    checks = [_schema_check("pair.schema", problems, versions, mixed)]

    def result() -> ComparisonResult:
        return ComparisonResult(names[0], names[1], checks, asdict(tol), mixed_schema=mixed, schemas=versions)

    if problems:
        return result()
    tables = {label: [(pd.read_csv(d / "waveguide.csv"), pd.read_csv(d / "far_field.csv")) for d in dirs]
              for label, dirs in sides.items()}
    finite = _check_finite("pair.finite", [(f"{label} directory {i + 1} {name}", df)
                                           for label, pair in tables.items() for i, (wg, ff) in enumerate(pair)
                                           for name, df in (("waveguide.csv", wg), ("far_field.csv", ff))])
    checks.append(finite)
    if not finite.passed:
        return result()

    # 2.0 sides: the sidecar's records and lattice, then the dark rule, before any legacy view or masking (OR-038).
    v2_sides = [label for label in sides if versions[label] == SCHEMA_2]
    v2: dict[str, list[_V2Table]] = {}
    lattice: dict[str, tuple[str, ExitGrid | None]] = {}
    problems = []
    for label in v2_sides:
        (_, sidecar0, path0), (_, sidecar1, path1) = resolved[label]
        try:
            if sidecar0 is not None and sidecar0 != sidecar1:
                raise ComparisonInputError(f"the two directories' sidecars differ ({path0}, {path1}): not one dataset")
            v2[label] = [_v2_table(tables[label][i][0], sidecar0, path0, sides[label][i],
                                   f"{label} directory {i + 1} waveguide.csv") for i in (0, 1)]
            lattice[label] = _v2_lattice(sidecar0, f"the sidecar {path0}")
        except ComparisonInputError as exc:
            problems.append(f"{label}: {exc}")
    if v2_sides:
        checks.append(_listing("pair.sidecar", problems, "; ".join(
            f"{label}: {len(v2[label][0].records.keys)} keys, {lattice[label][0]} lattice" for label in v2_sides)
            if not problems else ""))
        if problems:
            return result()
        dark = [p for label in v2_sides for i, table in enumerate(v2[label])
                for p in _dark_contradictions(tables[label][i][0], table.records, table.index, table.ephi,
                                              f"{label} directory {i + 1} waveguide.csv")]
        checks.append(_listing("pair.dark_records", dark,
                               "every dark record carries t_flux and OutgoingPower exactly 0 (OR-038)"))
        if dark:
            return result()

    if mixed:
        label = v2_sides[0]
        tables[label] = [(legacy_view(wg, v2[label][i].sidecar, v2[label][i].ephi), ff)
                         for i, (wg, ff) in enumerate(tables[label])]
    cand, ref = tables["candidate"], tables["reference"]
    checks.append(_check_pair_ephi(cand, ref))
    if v2_sides and not mixed:
        (cs, cp), (rs, rp) = (_pair_stats_v2(tables[label], v2[label], lattice[label], tol) for label in sides)
    else:
        (cs, cp), (rs, rp) = _pair_stats(cand, tol), _pair_stats(ref, tol)
    problems = [f"candidate {p}" for p in cp] + [f"reference {p}" for p in rp]
    for (cw, _), (rw, _), e in zip(cand, ref, (0, 1)):
        if _keys(cw) != _keys(rw):
            problems.append(f"Ephi={e} incidence keys differ from the reference: missing {sorted(_keys(rw) - _keys(cw))} "
                            f"extra {sorted(_keys(cw) - _keys(rw))}")
    checks.append(Check("pair.points", not problems, "fail",
                        "both pairs pair by row at every key" if not problems else "; ".join(problems[:5])))
    checks += [_check_rho("pair.exit_overlap", "rho_exit", cs, rs, tol),
               _check_rho("pair.far_field_overlap", "rho_far", cs, rs, tol),
               _check_mixed(cs, rs, tol)]
    return result()
