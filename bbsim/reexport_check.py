"""Acceptance of a schema-2 re-export against its schema-1 source (schema-2 design, section 2.5).

Usage::

    python -m bbsim.reexport_check --source JOB --fresh JOB [--poynting JSON]

``--source`` is a finished schema-1 job; ``--fresh`` is the job ``run_frequency --from-solved`` made from it on a copy
of its solved project (HF-050). On the HPC this is the gate before the ``_v2`` trees are assembled. Run it from the
checkout that made the re-export: the embedded column-order reference is checked against this checkout's exporter
files, as the runner checks it.

Only a declared schema-1 source and a declared schema-2 fresh job are compared. The versions come from each job's
``*.dataset.json`` sidecar, never from a config without ``output.schema``; any other pair, a missing sidecar or a
declaration that a manifest or a CSV header contradicts fails before any value is read. Then, in order:

1. Identity: dataset, frequency, keys, lattice and cross-section; the fresh sidecar's re-export record names this
   source's manifest; every value finite; every exit point (X, Y, Z) equal to 1e-12 m in the same order, and the
   far-field angles likewise.
2. Column order (HF-041, HF-052): the fresh job's embedded reference is validated as the runner validates it; for
   the crack1 500 GHz reference job its kept raw tables are re-evaluated and the source's raw fields are judged by
   the same rule. Other jobs are not held to an Ex-dominance rule; step 3 compares their slots by position.
3. Exit field (CDX-026), per key and polarization, by position. A source whose six values are all zero needs a dark
   record (W = 1, exactly zero fresh values) and no W is fitted; a dark record keeps W = 1 and the source values
   (OR-038); otherwise the least-squares W must be positive and equal the sidecar's w0 / w1 to 1e-6, and every value
   lie within eps of max |W source|. Independently of W, sum_i w_i |E_t,i|^2 h_y h_z / (2 eta0) from the fresh Ex and
   Ey (quadrature weights from the declared grid) must equal OutgoingPower to 1e-6 on every non-dark record.
4. Far field, not rescaled: within eps of the source; a zero source key exactly zero.
5. Rim (HF-029, discs only): the source's zero-filled points come back exactly zero and the zero-fill counts match.
6. Powers: the fresh manifest's raw powers complete and finite; the dark flag recomputed (|flux| <= 1e-12 A I);
   stored OutgoingPower zero on every row of a dark record and the raw flux otherwise; t_flux = OutgoingPower / (A I);
   IngoingPower = A I; InterceptFactor constant within each key and equal to c_key, itself recomputed from the
   declared grid, to 1e-12 (OR-034); the fresh |Re S| equal to the source's OutgoingPower to 1e-4; with
   ``--poynting`` (hpc/diagnose_poynting.py's poynting.json), the raw flux equal to step 8b's ``flux_z_w`` to 1e-4,
   absolute 1e-12 A I at zero (CDX-009).
7. Convergence (CDX-025): the source's saved tables re-qualified at its recorded target; both polarizations must
   qualify and agree field by field with the fresh qualification, and with the source manifest's record if it has one.

Missing, malformed or non-finite input fails the row that needs it; nothing is compared with NaN. Every row, the
largest deviations and the verdict go to ``<fresh>/reexport_check.json``; the printed table ends in
``REEXPORT: PASS`` or ``REEXPORT: FAIL <n>`` (n failing rows). Exit code 0 on PASS, 1 on FAIL, 2 when a job
directory or the --poynting file is missing.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from collections.abc import Callable, Mapping, Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from bbsim import dataset_v2, schema
from bbsim.compare import Check
from bbsim.convergence import qualify
from bbsim.extract import ExitGrid
from bbsim.manifest import REPO_ROOT, sha256_file, write_json
from bbsim.naming import dataset_dir_name

REPORT_NAME = "reexport_check.json"
SOURCE_SCHEMA = "1.0"
FRESH_SCHEMA = "2.0"
SUPPORTED = f"a schema {SOURCE_SCHEMA!r} source re-exported as schema {FRESH_SCHEMA!r}"

EPS = 1e-6                 # field match: |fresh - W source| <= EPS * max |W source| (CDX-026)
SCALE_REL = 1e-6           # least-squares W against the sidecar's w0 / w1
NORMALIZATION_REL = 1e-6   # the fields' weighted tangential integral against OutgoingPower
INTERCEPT_REL = 1e-12      # InterceptFactor and the sidecar's c_k against c_key recomputed from the grid (OR-034)
ABS_RE_S_REL = 1e-4        # fresh |Re S| against the source's OutgoingPower, as step 8b compared them
POYNTING_REL = 1e-4        # fresh raw flux against step 8b's flux_z_w
STORED_REL = 1e-12         # a stored OutgoingPower, t_flux, IngoingPower or |Re S| against the value it records
COORD_M = 1e-12            # exit points, source against fresh: BBRsim pairs rows to this
KEY_DEG = dataset_v2.KEY_TOL_DEG
DARK_REL = dataset_v2.DARK_REL
MAX_LISTED = 5

_WG = {name: i for i, name in enumerate(schema.WAVEGUIDE_COLUMNS_V2)}   # positions 0-8 are the same in 1.0
_FF = {name: i for i, name in enumerate(schema.FAR_FIELD_COLUMNS)}
FIELDS = schema.EXIT_FIELD_SLOTS                                        # the six stored exit-field slots
FAR_VALUES = slice(_FF["rEphi_real"], _FF["rEtheta_imag"] + 1)          # the four far-field values
OUTGOING, INGOING, INTERCEPT = _WG["OutgoingPower"], _WG["IngoingPower"], _WG[schema.INTERCEPT_FACTOR]
XYZ = slice(_WG["X"], _WG["Z"] + 1)
QUALIFICATION_FIELDS = ("achieved", "target", "passes_completed", "max_passes", "converged", "qualified", "reasons")
EXIT_ROWS = ("exit_field.dark", "exit_field.scale", "exit_field.values", "exit_field.normalization")

_MISSING = object()


class _Refused(ValueError):
    """Input a check cannot accept: it fails the row that needed it."""


# --- small helpers ---------------------------------------------------------------------------------------------------

def _key_text(key: tuple[float, float]) -> str:
    return f"(phi, theta) = ({key[0]:g}, {key[1]:g}) deg"


def _where(ephi: int, key: tuple[float, float]) -> str:
    return f"Ephi={ephi} at {_key_text(key)}"


def _number(value: object, where: str) -> float:
    """A finite real number; anything else (bool, string, null, NaN, infinity) is refused."""
    if isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)):
        raise _Refused(f"{where} must be a number, got {value!r}")
    number = float(value)
    if not math.isfinite(number):
        raise _Refused(f"{where} is not finite: {value!r}")
    return number


def _integer(value: object, where: str) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise _Refused(f"{where} must be an integer, got {value!r}")
    return int(value)


def _get(mapping: object, *path: str, where: str, default: object = _MISSING) -> Any:
    node = mapping
    for name in path:
        if not isinstance(node, Mapping) or name not in node:
            if default is _MISSING:
                raise _Refused(f"{where} lacks {'.'.join(path)}")
            return default
        node = node[name]
    return node


def _read_json(path: Path, what: str) -> dict:
    try:
        data = json.loads(Path(path).read_text())
    except (OSError, ValueError) as exc:
        raise _Refused(f"{what} {path} is unreadable: {exc}") from exc
    if not isinstance(data, dict):
        raise _Refused(f"{what} {path} is not a JSON object")
    return data


def _relative(a: float, b: float) -> float:
    """|a - b| / |b|, 0 when both are 0, infinity when only b is."""
    if a == b:
        return 0.0
    return abs(a - b) / abs(b) if b != 0.0 else math.inf


def _same_value(a: object, b: object, rel: float = 1e-12) -> bool:
    numeric = (int, float, np.integer, np.floating)
    if isinstance(a, numeric) and isinstance(b, numeric) and not isinstance(a, bool) and not isinstance(b, bool):
        return math.isfinite(float(a)) and math.isfinite(float(b)) and math.isclose(float(a), float(b), rel_tol=rel)
    if isinstance(a, Mapping) and isinstance(b, Mapping):
        return set(a) == set(b) and all(_same_value(a[k], b[k], rel) for k in a)
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        return len(a) == len(b) and all(_same_value(x, y, rel) for x, y in zip(a, b))
    return a == b


def _key_index(keys: Sequence[tuple[float, float]], phi: float, theta: float) -> int | None:
    for k, (p, t) in enumerate(keys):
        if abs(p - phi) <= KEY_DEG and abs(t - theta) <= KEY_DEG:
            return k
    return None


def _phi_theta(record: object, where: str) -> tuple[float, float]:
    if not isinstance(record, Mapping):
        raise _Refused(f"{where} is not an object")
    return _number(record.get("phi_deg"), f"{where}.phi_deg"), _number(record.get("theta_deg"), f"{where}.theta_deg")


def _rim_key(record: object, where: str) -> tuple[float, float]:
    key = record.get("key") if isinstance(record, Mapping) else None
    if not isinstance(key, (list, tuple)) or len(key) != 2:
        raise _Refused(f"{where}.key must be [phi, theta], got {key!r}")
    return _number(key[0], f"{where}.key[0]"), _number(key[1], f"{where}.key[1]")


def _match_records(records: object, keys: Sequence[tuple[float, float]], where: str,
                   key_of: Callable[[object, str], tuple[float, float]] = _phi_theta) -> tuple[dict[int, Any], list[str]]:
    """The records of a per-key list by declared key index, and what is wrong with the list: not a list, a record
    off the declared grid or unreadable, two records at one key, declared keys without a record."""
    if not isinstance(records, list):
        return {}, [f"{where} must be a list of per-key records, got {type(records).__name__}"]
    found: dict[int, Any] = {}
    problems: list[str] = []
    for i, record in enumerate(records):
        try:
            phi, theta = key_of(record, f"{where}[{i}]")
        except _Refused as exc:
            problems.append(str(exc))
            continue
        k = _key_index(keys, phi, theta)
        if k is None:
            problems.append(f"{where} has a record at {_key_text((phi, theta))}, outside the declared grid")
        elif k in found:
            problems.append(f"{where} has two records at {_key_text(keys[k])}")
        else:
            found[k] = record
    missing = [keys[k] for k in range(len(keys)) if k not in found]
    if missing:
        problems.append(f"{where} lacks {len(missing)} declared keys, the first at {_key_text(missing[0])}")
    return found, problems


def _strict(records: object, keys: Sequence[tuple[float, float]], where: str, **kwargs: Any) -> dict[int, Any]:
    found, problems = _match_records(records, keys, where, **kwargs)
    if problems:
        raise _Refused("; ".join(problems))
    return found


def _row_keys(phis: np.ndarray, thetas: np.ndarray, keys: Sequence[tuple[float, float]]) -> np.ndarray:
    """The declared key index of every table row, -1 for a row at no declared key."""
    index = np.full(len(phis), -1, dtype=int)
    for k, (phi, theta) in enumerate(keys):
        index[(np.abs(phis - phi) <= KEY_DEG) & (np.abs(thetas - theta) <= KEY_DEG)] = k
    return index


def _row(name: str, rule: str, failures: list[str], summary: str, values: dict | None = None) -> Check:
    """A fail-severity row whose detail always states its rule, then the summary or the first failures."""
    if failures:
        listed = "; ".join(failures[:MAX_LISTED])
        more = f"; and {len(failures) - MAX_LISTED} more" if len(failures) > MAX_LISTED else ""
        detail = f"{rule}: {len(failures)} failure(s): {listed}{more}"
    else:
        detail = f"{rule}: {summary}"
    return Check(name, not failures, "fail", detail, {**(values or {}), "failures": failures})


def _exit_grid(record: object, where: str) -> tuple[ExitGrid, tuple[float, float]]:
    """The declared exit lattice as an ``ExitGrid`` (mm) for ``dataset_v2.node_weights``, and its steps in metres.

    Each axis's declared step must span min to max in count - 1 steps, to ``dataset_v2.NODE_TOL`` of the step."""
    lo, hi, steps_m, counts = [0.0], [0.0], [], []
    for axis in ("y_e", "z_e"):
        rec = _get(record, axis, where=where)
        lo_m = _number(_get(rec, "min", where=f"{where}.{axis}"), f"{where}.{axis}.min")
        hi_m = _number(_get(rec, "max", where=f"{where}.{axis}"), f"{where}.{axis}.max")
        step_m = _number(_get(rec, "step", where=f"{where}.{axis}"), f"{where}.{axis}.step")
        count = _integer(_get(rec, "count", where=f"{where}.{axis}"), f"{where}.{axis}.count")
        if step_m <= 0.0 or count < 1:
            raise _Refused(f"{where}.{axis}: step {step_m!r} m and count {count} must be positive")
        if count > 1 and abs((hi_m - lo_m) / (count - 1) - step_m) > dataset_v2.NODE_TOL * step_m:
            raise _Refused(f"{where}.{axis}: the declared step {step_m!r} m does not span {lo_m!r} to {hi_m!r} m in "
                           f"{count - 1} steps")
        lo.append(lo_m * 1e3)
        hi.append(hi_m * 1e3)
        steps_m.append(step_m)
        counts.append(count)
    steps_mm = (0.0, steps_m[0] * 1e3, steps_m[1] * 1e3)
    grid = ExitGrid(lo_mm=tuple(lo), hi_mm=tuple(hi), resolution_mm=steps_mm, step_mm=steps_mm)
    if tuple(grid.counts[1:]) != tuple(counts):
        raise _Refused(f"{where}: the declared counts {counts} are not the lattice's {list(grid.counts[1:])}")
    return grid, (steps_m[0], steps_m[1])


# --- the result ------------------------------------------------------------------------------------------------------

@dataclass
class ReexportResult:
    """The rows in the design's order, the largest deviation of each compared quantity (value, limit, where), the
    tolerances used, and the verdict."""

    source: str
    fresh: str
    poynting: str | None
    checks: list[Check]
    largest_deviation: dict[str, dict[str, Any]]
    tolerances: dict[str, float]

    @property
    def failures(self) -> list[Check]:
        return [c for c in self.checks if c.severity == "fail" and not c.passed]

    @property
    def passed(self) -> bool:
        return bool(self.checks) and not self.failures

    @property
    def verdict(self) -> str:
        return "REEXPORT: PASS" if self.passed else f"REEXPORT: FAIL {len(self.failures)}"

    def to_dict(self) -> dict[str, Any]:
        return {"source": self.source, "fresh": self.fresh, "poynting": self.poynting, "passed": self.passed,
                "failed": len(self.failures), "verdict": self.verdict, "tolerances": dict(self.tolerances),
                "checks": [asdict(c) for c in self.checks],
                "largest_deviations": {k: dict(v) for k, v in self.largest_deviation.items()}}

    def format_table(self) -> str:
        lines = [f"source:   {self.source}", f"fresh:    {self.fresh}"]
        if self.poynting:
            lines.append(f"poynting: {self.poynting}")
        lines += ["tolerances: " + " ".join(f"{k}={v:g}" for k, v in self.tolerances.items()), ""]
        for c in self.checks:
            status = ("PASS" if c.passed else "FAIL") if c.severity == "fail" else ("INFO" if c.passed else "SKIP")
            lines.append(f"{status:4} {c.name:25} {c.detail}")
        if self.largest_deviation:
            lines += ["", "largest deviations:"]
            for name, d in self.largest_deviation.items():
                lines.append(f"  {name:16} {d['value']:.3g} (limit {d['limit']:.3g}) at {d['where']}")
        lines += ["", self.verdict]
        return "\n".join(lines)


# --- the acceptance --------------------------------------------------------------------------------------------------

@dataclass
class _Job:
    label: str
    root: Path
    sidecar: dict = field(default_factory=dict)
    manifest: dict = field(default_factory=dict)
    # Per Ephi, every CSV column as float by its position; column 0 (the Freq label) is NaN, its labels kept apart.
    waveguide: dict[int, np.ndarray] = field(default_factory=dict)
    far_field: dict[int, np.ndarray] = field(default_factory=dict)
    freq_labels: dict[tuple[str, int], set[str]] = field(default_factory=dict)

    def csv(self, ephi: int, name: str) -> Path:
        dataset_id = _get(self.sidecar, "dataset_id", where=f"the {self.label} sidecar")
        frequency = _number(_get(self.sidecar, "frequency_ghz", where=f"the {self.label} sidecar"),
                            f"the {self.label} sidecar's frequency_ghz")
        if not isinstance(dataset_id, str) or not dataset_id:
            raise _Refused(f"the {self.label} sidecar's dataset_id must be a non-empty string, got {dataset_id!r}")
        return self.root / dataset_dir_name(dataset_id, frequency, ephi) / name


class _Acceptance:
    def __init__(self, source: Path, fresh: Path, poynting: Path | None, eps: float):
        self.source, self.fresh = _Job("source", source), _Job("fresh", fresh)
        self.poynting, self.eps = poynting, eps
        self.largest: dict[str, dict[str, Any]] = {}
        self.keys: list[tuple[float, float]] | None = None
        self.theta_list: list[float] = []
        self.grid: ExitGrid | None = None
        self.h_m: tuple[float, float] = (0.0, 0.0)
        self.shape: str | None = None
        self.tables_ok = False
        self.aligned = False
        self.rows: dict[int, list[np.ndarray]] = {}
        self.ff_rows: dict[int, list[np.ndarray]] = {}
        self._cache: dict[str, tuple[Any, BaseException | None]] = {}

    # -- plumbing ------------------------------------------------------------------------------------------------
    def _note(self, name: str, value: float, where: str, limit: float) -> None:
        current = self.largest.get(name)
        if current is None or not value <= current["value"]:      # a NaN replaces, so it shows
            self.largest[name] = {"value": float(value), "limit": float(limit), "where": where}

    def _cached(self, name: str, compute: Callable[[], Any]) -> Any:
        if name not in self._cache:
            try:
                self._cache[name] = (compute(), None)
            except Exception as exc:  # noqa: BLE001 - re-raised to every row that needs it
                self._cache[name] = (None, exc)
        value, exc = self._cache[name]
        if exc is not None:
            raise exc
        return value

    @staticmethod
    def _guard(names: Sequence[str], step: Callable[[], list[Check]]) -> list[Check]:
        try:
            return step()
        except Exception as exc:  # noqa: BLE001 - every refusal becomes a failing row, never a traceback (CDX-030)
            text = str(exc) if isinstance(exc, _Refused) else f"{type(exc).__name__}: {exc}"
            return [Check(name, False, "fail", f"refused: {text}") for name in names]

    def _needs(self, names: Sequence[str], step: Callable[[], list[Check]], *, keys: bool = False,
               tables: bool = False, aligned: bool = False) -> list[Check]:
        why = None
        if (keys or aligned) and self.keys is None:
            why = "the declared keys could not be read (identity.dataset)"
        elif (tables or aligned) and not self.tables_ok:
            why = "the tables are not usable (identity.finite)"
        elif aligned and (self.grid is None or self.shape is None):
            why = "the declared lattice could not be read (identity.dataset)"
        elif aligned and not self.aligned:
            why = "the source and fresh rows do not align (identity.points)"
        if why is not None:
            return [Check(name, False, "info", f"not run: {why}") for name in names]
        return self._guard(names, step)

    def run(self) -> list[Check]:
        checks = self._guard(["schema.declared"], lambda: [self._schemas()])
        if not checks[0].passed:
            return checks
        checks += self._guard(["identity.dataset"], lambda: [self._identity_dataset()])
        checks += self._guard(["identity.finite"], lambda: [self._load_tables()])
        checks += self._needs(["identity.points"], lambda: [self._identity_points()], keys=True, tables=True)
        checks += self._guard(["column_order.reference"], lambda: [self._order_reference()])
        checks += self._needs(["column_order.raw_tables"], lambda: [self._order_raw_tables()], tables=True)
        checks += self._needs(EXIT_ROWS, self._exit_field, aligned=True)
        checks += self._needs(["far_field.values"], lambda: [self._far_field()], aligned=True)
        checks += self._needs(["rim.zero_fill"], lambda: [self._rim()], aligned=True)
        checks += self._needs(["powers.raw_records"], lambda: [self._raw_records()], keys=True)
        for name, step in (("powers.dark", self._dark), ("powers.outgoing", self._outgoing),
                           ("powers.t_flux", self._t_flux), ("powers.ingoing", self._ingoing),
                           ("powers.intercept_factor", self._intercept_constant), ("powers.c_key", self._c_key),
                           ("powers.abs_re_s", self._abs_re_s), ("powers.poynting", self._poynting)):
            checks += self._needs([name], lambda step=step: [step()], aligned=True)
        checks += self._guard(["convergence"], lambda: [self._convergence()])
        return checks

    # -- 0. declared schemas -------------------------------------------------------------------------------------
    def _declaration(self, job: _Job, problems: list[str]) -> object:
        sidecars = sorted(job.root.glob("*.dataset.json"))
        if len(sidecars) != 1:
            names = f": {[p.name for p in sidecars]}" if sidecars else ""
            problems.append(f"the {job.label} job holds {len(sidecars)} dataset sidecars (*.dataset.json), not 1{names}")
            return None
        try:
            job.sidecar = _read_json(sidecars[0], f"the {job.label} sidecar")
        except _Refused as exc:
            problems.append(str(exc))
            return None
        if "schema_version" not in job.sidecar:
            problems.append(f"the {job.label} sidecar {sidecars[0].name} declares no schema_version")
            return None
        return job.sidecar["schema_version"]

    def _contradictions(self) -> list[str]:
        problems: list[str] = []
        for job, number, columns in ((self.source, 1, schema.WAVEGUIDE_COLUMNS),
                                     (self.fresh, 2, schema.WAVEGUIDE_COLUMNS_V2)):
            version = job.sidecar["schema_version"]
            try:
                job.manifest = _read_json(job.root / "manifest.json", f"the {job.label} manifest")
            except _Refused as exc:
                problems.append(str(exc))
                continue
            recorded = _get(job.manifest, "config", "output", "schema", where="", default=None)
            declares = type(recorded) is int and recorded == number      # config.py records the integer 1 or 2
            if job is self.fresh and not declares:
                problems.append(f"the fresh manifest's config.output.schema is {recorded!r}, which contradicts its "
                                f"{version!r} sidecar (a schema-2 run records 2)")
            if job is self.source and recorded is not None and not declares:
                problems.append(f"the source manifest's config.output.schema is {recorded!r}, which contradicts its "
                                f"{version!r} sidecar")
            if _get(job.sidecar, "exit_field", "columns", where="", default=None) != list(columns):
                problems.append(f"the {job.label} sidecar's exit_field.columns are not the {version} waveguide columns")
            if _get(job.sidecar, "far_field", "columns", where="", default=None) != list(schema.FAR_FIELD_COLUMNS):
                problems.append(f"the {job.label} sidecar's far_field.columns are not the far-field columns")
            for ephi in (0, 1):
                for name, expected in (("waveguide.csv", columns), ("far_field.csv", schema.FAR_FIELD_COLUMNS)):
                    try:
                        path = job.csv(ephi, name)
                        header = list(pd.read_csv(path, nrows=0).columns)
                    except (_Refused, OSError, ValueError) as exc:
                        problems.append(f"the {job.label} {name} of Ephi={ephi} is unreadable: {exc}")
                        continue
                    if header != list(expected):
                        problems.append(f"the {job.label} {path.parent.name}/{name} header {header} contradicts its "
                                        f"{version!r} sidecar, which declares {list(expected)}")
        return problems

    def _schemas(self) -> Check:
        problems: list[str] = []
        declared = {job.label: self._declaration(job, problems) for job in (self.source, self.fresh)}
        sv, fv = declared["source"], declared["fresh"]
        if (sv, fv) != (SOURCE_SCHEMA, FRESH_SCHEMA):
            problems.insert(0, f"source declares {sv!r}, fresh declares {fv!r}: the supported comparison is "
                               f"{SUPPORTED}")
        else:
            problems += self._contradictions()
        detail = "; ".join(problems) if problems else (f"source {sv!r} -> fresh {fv!r} from their sidecars; "
                                                       "manifests and CSV headers agree")
        return Check("schema.declared", not problems, "fail", detail,
                     {"source": sv, "fresh": fv, "supported": SUPPORTED, "failures": problems})

    # -- 1. identity -----------------------------------------------------------------------------------------------
    @staticmethod
    def _declared_axes(sidecar: dict, where: str) -> tuple[list[float], list[float]]:
        axes = []
        for name in ("incident_phi_deg", "incident_theta_deg"):
            values = _get(sidecar, "excitation", name, where=where)
            if not isinstance(values, list) or not values:
                raise _Refused(f"{where}'s excitation.{name} must be a non-empty list, got {values!r}")
            numbers = [_number(v, f"{where}'s excitation.{name}") for v in values]
            ordered = sorted(numbers)
            if any(b - a <= KEY_DEG for a, b in zip(ordered, ordered[1:])):
                raise _Refused(f"{where}'s excitation.{name} repeats a value: {numbers}")
            axes.append(numbers)
        return axes[0], axes[1]

    def _identity_dataset(self) -> Check:
        s, f = self.source.sidecar, self.fresh.sidecar
        phis, thetas = self._declared_axes(f, "the fresh sidecar")
        self.keys = [(p, t) for p in phis for t in thetas]
        self.theta_list = thetas
        problems: list[str] = []
        try:
            src_phis, src_thetas = self._declared_axes(s, "the source sidecar")
            for name, a, b in (("incident_phi_deg", src_phis, phis), ("incident_theta_deg", src_thetas, thetas)):
                if len(a) != len(b) or any(abs(x - y) > KEY_DEG for x, y in zip(a, b)):
                    problems.append(f"declared {name}: source {a}, fresh {b}")
        except _Refused as exc:
            problems.append(str(exc))
        for name in ("dataset_id", "frequency_label"):
            if s.get(name) != f.get(name) or not isinstance(f.get(name), str):
                problems.append(f"{name}: source {s.get(name)!r}, fresh {f.get(name)!r}")
        if not _same_value(s.get("frequency_ghz"), f.get("frequency_ghz"), rel=0.0):
            problems.append(f"frequency_ghz: source {s.get('frequency_ghz')!r}, fresh {f.get('frequency_ghz')!r}")

        grid_s = _get(s, "exit_field", "grid", where="the source sidecar")
        grid_f = _get(f, "exit_field", "grid", where="the fresh sidecar")
        for axis in ("y_e", "z_e"):
            for name in ("min", "max"):
                a = _number(_get(grid_s, axis, name, where="the source exit_field.grid"), f"source grid {axis}.{name}")
                b = _number(_get(grid_f, axis, name, where="the fresh exit_field.grid"), f"fresh grid {axis}.{name}")
                if abs(a - b) > COORD_M:
                    problems.append(f"exit_field.grid.{axis}.{name}: source {a!r} m, fresh {b!r} m")
            if _get(grid_s, axis, "count", where="the source exit_field.grid") != _get(
                    grid_f, axis, "count", where="the fresh exit_field.grid"):
                problems.append(f"exit_field.grid.{axis}.count: source {grid_s[axis]['count']!r}, "
                                f"fresh {grid_f[axis]['count']!r}")
        self.grid, self.h_m = _exit_grid(grid_f, "the fresh sidecar's exit_field.grid")

        cs_s = _get(s, "exit_field", "cross_section", where="the source sidecar")
        cs_f = _get(f, "exit_field", "cross_section", where="the fresh sidecar")
        if not _same_value(cs_s, cs_f):
            problems.append(f"exit_field.cross_section: source {cs_s!r}, fresh {cs_f!r}")
        shape = cs_f.get("shape") if isinstance(cs_f, Mapping) else None
        if shape not in dataset_v2.QUADRATURE_BY_SHAPE:
            raise _Refused(f"the fresh cross-section shape {shape!r} is not one of "
                           f"{sorted(dataset_v2.QUADRATURE_BY_SHAPE)}")
        self.shape = shape
        quadrature = _get(f, "exit_field", "quadrature", where="the fresh sidecar", default=None)
        if quadrature != dataset_v2.QUADRATURE_BY_SHAPE[self.shape]:
            problems.append(f"the fresh exit_field.quadrature {quadrature!r} is not the {self.shape}'s "
                            f"{dataset_v2.QUADRATURE_BY_SHAPE[self.shape]!r}")
        retained = [_get(j.sidecar, "exit_field", "points_per_key_retained", where=f"the {j.label} sidecar")
                    for j in (self.source, self.fresh)]
        if retained[0] != retained[1]:
            problems.append(f"exit_field.points_per_key_retained: source {retained[0]!r}, fresh {retained[1]!r}")

        source_sha = sha256_file(self.source.root / "manifest.json")
        reexport = _get(f, "provenance", "reexport", where="the fresh sidecar", default=None)
        if not isinstance(reexport, Mapping):
            problems.append("the fresh sidecar has no provenance.reexport: it is not a re-export")
        elif reexport.get("source_manifest_sha256") != source_sha:
            problems.append(f"the fresh sidecar's provenance.reexport.source_manifest_sha256 "
                            f"{reexport.get('source_manifest_sha256')!r} is not the source manifest's {source_sha}: "
                            "the fresh job was not re-exported from this source")
        return _row("identity.dataset", "same dataset, frequency, keys, lattice and cross-section; the fresh job "
                    "re-exports this source", problems,
                    f"{f.get('dataset_id')} at {f.get('frequency_label')}, {len(self.keys)} keys, {self.shape} lattice "
                    f"{self.grid.counts[1]} x {self.grid.counts[2]}, {retained[1]} points per key")

    def _load_tables(self) -> Check:
        problems: list[str] = []
        for job in (self.source, self.fresh):
            for ephi in (0, 1):
                for name, store in (("waveguide.csv", job.waveguide), ("far_field.csv", job.far_field)):
                    path = job.csv(ephi, name)
                    where = f"the {job.label} {path.parent.name}/{name}"
                    frame = pd.read_csv(path)
                    try:
                        numeric = frame.iloc[:, 1:].to_numpy(dtype=float)
                    except (TypeError, ValueError) as exc:
                        problems.append(f"{where} holds non-numeric values: {exc}")
                        continue
                    bad = ~np.isfinite(numeric)
                    if bad.any():
                        i, j = (int(v) for v in np.argwhere(bad)[0])
                        problems.append(f"{where} holds {int(bad.sum())} non-finite values, the first in "
                                        f"{frame.columns[j + 1]} on row {i}, at "
                                        f"{_key_text((numeric[i, 1], numeric[i, 2]))}")
                        continue
                    store[ephi] = np.column_stack([np.full(len(frame), np.nan), numeric])
                    job.freq_labels[(name, ephi)] = set(frame.iloc[:, 0].astype(str))
        self.tables_ok = not problems
        return _row("identity.finite", "every value of the four CSVs of both jobs numeric and finite", problems,
                    "all finite")

    def _identity_points(self) -> Check:
        problems: list[str] = []
        aligned = True
        label = self.fresh.sidecar.get("frequency_label")
        retained = self.fresh.sidecar["exit_field"]["points_per_key_retained"]
        keys = self.keys
        worst, worst_where = 0.0, "every row"
        for ephi in (0, 1):
            for kind, tables in (("waveguide", (self.source.waveguide, self.fresh.waveguide)),
                                 ("far_field", (self.source.far_field, self.fresh.far_field))):
                src, fresh = tables[0][ephi], tables[1][ephi]
                where = f"{kind} Ephi={ephi}"
                for job in (self.source, self.fresh):
                    labels = job.freq_labels[(f"{kind}.csv", ephi)]
                    if labels != {label}:
                        problems.append(f"the {job.label} {where} Freq labels {sorted(labels)}, not {label!r}")
                for job, table in ((self.source, src), (self.fresh, fresh)):
                    if not np.all(table[:, _WG["Ephi"]] == ephi):
                        problems.append(f"the {job.label} {where} Ephi column is not {ephi} on every row")
                if len(src) != len(fresh):
                    problems.append(f"{where}: the source has {len(src)} rows, the fresh {len(fresh)}")
                    aligned = False
                    continue
                index_s = _row_keys(src[:, _WG["IWavePhi"]], src[:, _WG["IWaveTheta"]], keys)
                index_f = _row_keys(fresh[:, _WG["IWavePhi"]], fresh[:, _WG["IWaveTheta"]], keys)
                if (index_f < 0).any() or (index_s < 0).any():
                    problems.append(f"{where}: rows at keys outside the declared grid")
                    aligned = False
                    continue
                if not np.array_equal(index_s, index_f):
                    i = int(np.flatnonzero(index_s != index_f)[0])
                    problems.append(f"{where}: the keys differ from row {i} on (source {_key_text(keys[index_s[i]])}, "
                                    f"fresh {_key_text(keys[index_f[i]])})")
                    aligned = False
                    continue
                rows = [np.flatnonzero(index_f == k) for k in range(len(keys))]
                empty = [keys[k] for k, r in enumerate(rows) if r.size == 0]
                if empty:
                    problems.append(f"{where}: no rows at the declared {_key_text(empty[0])}")
                    aligned = False
                    continue
                if kind == "waveguide":
                    self.rows[ephi] = rows
                    for k, r in enumerate(rows):
                        if r.size != retained:
                            problems.append(f"{_where(ephi, keys[k])}: {r.size} exit points, the sidecar declares "
                                            f"{retained}")
                    diff = np.abs(src[:, XYZ] - fresh[:, XYZ])
                    i, axis = (int(v) for v in np.unravel_index(int(np.argmax(diff)), diff.shape))
                    if diff[i, axis] > worst:
                        worst, worst_where = float(diff[i, axis]), f"{_where(ephi, keys[index_f[i]])}, {'XYZ'[axis]}"
                    if not diff[i, axis] <= COORD_M:
                        start = int(rows[index_f[i]][0])
                        problems.append(f"{_where(ephi, keys[index_f[i]])}: exit point {i - start} of the key moved "
                                        f"{diff[i, axis]:.3g} m in {'XYZ'[axis]} (limit {COORD_M:g} m)")
                else:
                    self.ff_rows[ephi] = rows
                    angles = np.abs(src[:, _FF["Phi"]:_FF["Theta"] + 1] - fresh[:, _FF["Phi"]:_FF["Theta"] + 1])
                    if not np.all(angles <= KEY_DEG):
                        i = int(np.flatnonzero(~np.all(angles <= KEY_DEG, axis=1))[0])
                        problems.append(f"far field {_where(ephi, keys[index_f[i]])}: the (Phi, Theta) of row {i} "
                                        f"differ by {float(angles[i].max()):.3g} deg")
        self.aligned = aligned
        self._note("points_m", worst, worst_where, COORD_M)
        return _row("identity.points", f"exit points X, Y, Z equal to {COORD_M:g} m and far-field angles to "
                    f"{KEY_DEG:g} deg, row by row in the same order", problems,
                    f"every row matches, the largest exit-point difference {worst:.3g} m")

    # -- 2. column order -------------------------------------------------------------------------------------------
    def _is_reference_dataset(self) -> bool:
        f = self.fresh.sidecar
        return (f.get("dataset_id") == dataset_v2.REFERENCE_DATASET
                and _same_value(f.get("frequency_ghz"), dataset_v2.REFERENCE_GHZ, rel=0.0))

    def _order_block(self) -> tuple[str, dict]:
        block = _get(self.fresh.sidecar, "exit_field", "column_order_check", where="the fresh sidecar")
        if not isinstance(block, Mapping) or block.get("rule") != dataset_v2.ORDER_RULE \
                or block.get("method") not in dataset_v2.ORDER_METHODS or not isinstance(block.get("reference"), Mapping):
            raise _Refused(f"the fresh exit_field.column_order_check must hold rule {dataset_v2.ORDER_RULE!r}, a method "
                           f"in {dataset_v2.ORDER_METHODS} and a reference object")
        return block["method"], dict(block["reference"])

    def _order_reference(self) -> Check:
        method, reference = self._order_block()
        provenance = _get(self.fresh.sidecar, "provenance", where="the fresh sidecar")
        aedt, pyaedt = provenance.get("aedt_version"), provenance.get("pyaedt_version")
        problems: list[str] = []
        try:
            dataset_v2.check_reference_compatible(reference, aedt, pyaedt, REPO_ROOT)
        except dataset_v2.DatasetV2Error as exc:
            problems.append(str(exc))
        if method == "reference" and not self._is_reference_dataset():
            problems.append(f"method 'reference' on {self.fresh.sidecar.get('dataset_id')!r} at "
                            f"{self.fresh.sidecar.get('frequency_ghz')!r} GHz: only {dataset_v2.REFERENCE_DATASET} at "
                            f"{dataset_v2.REFERENCE_GHZ:g} GHz measures the reference")
        return _row("column_order.reference", "the embedded reference complete, passed, and from this checkout's "
                    "exporter files and the job's AEDT / PyAEDT (HF-041, HF-052)", problems,
                    f"method {method!r}, AEDT {aedt!r}, PyAEDT {pyaedt!r}")

    @staticmethod
    def _case_values(table: np.ndarray, phi: float, theta: float) -> np.ndarray:
        mask = (np.abs(table[:, _WG["IWavePhi"]] - phi) <= KEY_DEG) & (np.abs(table[:, _WG["IWaveTheta"]] - theta)
                                                                        <= KEY_DEG)
        return table[mask, FIELDS]

    def _order_raw_tables(self) -> Check:
        rule = "the reference job's kept raw tables reproduce its cases; the source's raw fields pass the same rule"
        if not self._is_reference_dataset():
            return Check("column_order.raw_tables", True, "info",
                         f"not the column-order reference ({dataset_v2.REFERENCE_DATASET} at "
                         f"{dataset_v2.REFERENCE_GHZ:g} GHz): exit_field.values compares its slots by position")
        method, reference = self._order_block()
        problems: list[str] = []
        notes: list[str] = []
        if method == "reference":
            raw = {}
            for ephi in (0, 1):
                path = self.fresh.root / "scratch" / f"column_order_reference_Ephi{ephi}.csv"
                expected = _get(reference, "raw_table_sha256", str(ephi), where="the embedded reference")
                if not path.is_file():
                    problems.append(f"the reference job keeps no raw table {path}")
                    continue
                if sha256_file(path) != expected:
                    problems.append(f"{path.name} does not hash to the reference's raw_table_sha256[{ephi!r}]")
                    continue
                frame = pd.read_csv(path)
                if list(frame.columns) != list(schema.WAVEGUIDE_COLUMNS):
                    problems.append(f"{path.name} does not carry the 1.0 columns it was hashed with")
                    continue
                raw[ephi] = np.column_stack([np.full(len(frame), np.nan), frame.iloc[:, 1:].to_numpy(dtype=float)])
            for case in _get(reference, "cases", where="the embedded reference"):
                ephi, phi, theta = case["ephi"], case["phi_deg"], case["theta_deg"]
                if ephi not in raw:
                    continue
                got = dataset_v2.order_case(self._case_values(raw[ephi], phi, theta), ephi, phi, theta)
                for name in ("energy_fraction_12", "energy_fraction_2", "slot2_over_slot4", "discriminating",
                             "passed"):
                    if not _same_value(got[name], case[name]):
                        problems.append(f"{_where(ephi, (phi, theta))}: the kept raw table gives {name} "
                                        f"{got[name]!r}, the reference records {case[name]!r}")
            notes.append("the kept raw tables reproduce both cases")
        else:
            notes.append("method 'shared_exporter' keeps no raw tables")
        cases = []
        for ephi, phi, theta in dataset_v2.REFERENCE_CASES:
            values = self._case_values(self.source.waveguide[ephi], phi, theta)
            if values.shape[0] == 0:
                problems.append(f"the source has no rows at {_where(ephi, (phi, theta))}")
                continue
            case = dataset_v2.order_case(values, ephi, phi, theta)
            cases.append(case)
            if not case["passed"]:
                problems.append(f"the source's raw fields at {_where(ephi, (phi, theta))} fail the rule: (P1+P2)/P = "
                                f"{case['energy_fraction_12']:.6g}, P2/P = {case['energy_fraction_2']:.6g}, P2/P4 = "
                                f"{case['slot2_over_slot4']!r}")
        if cases and not any(c["discriminating"] for c in cases):
            problems.append("no source case discriminates (P2/P below "
                            f"{dataset_v2.DISCRIMINATING_FRACTION_2_MIN:g} in both)")
        notes.append("the source's raw fields pass the rule with "
                     f"{sum(c['discriminating'] for c in cases)} discriminating case(s)")
        return _row("column_order.raw_tables", rule, problems, "; ".join(notes), {"source_cases": cases})

    # -- 3. exit field -----------------------------------------------------------------------------------------------
    def _scale_records(self) -> dict[int, tuple[dict[int, float], dict[int, bool]]]:
        def parse() -> dict[int, tuple[dict[int, float], dict[int, bool]]]:
            where = "the fresh sidecar's exit_field.scale.per_key"
            found = _strict(_get(self.fresh.sidecar, "exit_field", "scale", "per_key", where="the fresh sidecar"),
                            self.keys, where)
            out = {}
            for k, record in found.items():
                at = f"{where} at {_key_text(self.keys[k])}"
                w = {e: _number(record.get(f"w{e}"), f"{at}: w{e}") for e in (0, 1)}
                dark = {}
                for e in (0, 1):
                    flag = record.get(f"dark{e}")
                    if not isinstance(flag, bool):
                        raise _Refused(f"{at}: dark{e} must be true or false, got {flag!r}")
                    dark[e] = flag
                out[k] = (w, dark)
            return out
        return self._cached("scale", parse)

    def _exit_field(self) -> list[Check]:
        scale = self._scale_records()
        dark_fail, scale_fail, value_fail, norm_fail = [], [], [], []
        dark_records, zero_sources, fitted = [], [], []
        n_normalized = 0
        h_y, h_z = self.h_m
        for ephi in (0, 1):
            src, fresh = self.source.waveguide[ephi], self.fresh.waveguide[ephi]
            for k, key in enumerate(self.keys):
                where = _where(ephi, key)
                r = self.rows[ephi][k]
                s, f = src[r, FIELDS], fresh[r, FIELDS]
                w_side, dark = scale[k][0][ephi], scale[k][1][ephi]
                s_max = float(np.max(np.abs(s)))
                if s_max == 0.0:
                    # Nothing to fit W against: the record must be dark, W = 1 and the fresh values exactly zero.
                    zero_sources.append([ephi, *key])
                    if not dark:
                        dark_fail.append(f"{where}: the source's six exit-field columns are all zero, so no W can be "
                                         f"fitted and the record must be dark, but dark{ephi} is false")
                    if w_side != 1.0:
                        dark_fail.append(f"{where}: a zero source needs W = 1, the sidecar has w{ephi} = {w_side!r}")
                    if np.any(f != 0.0):
                        dark_fail.append(f"{where}: a zero source must come back exactly zero; the largest fresh "
                                         f"|value| is {float(np.max(np.abs(f)))!r}")
                elif dark:
                    dark_records.append([ephi, *key])
                    if w_side != 1.0:
                        dark_fail.append(f"{where}: a dark record keeps W = 1, the sidecar has w{ephi} = {w_side!r}")
                    dev = float(np.max(np.abs(f - s))) / s_max
                    self._note("exit_field", dev, where, self.eps)
                    if not dev <= self.eps:
                        dark_fail.append(f"{where}: a dark record keeps the source values, but they differ by "
                                         f"{dev:.3g} of max |source| (limit {self.eps:g})")
                else:
                    fitted.append([ephi, *key])
                    w_fit = float(np.sum(f * s) / np.sum(s * s))
                    if not (math.isfinite(w_fit) and w_fit > 0.0):
                        scale_fail.append(f"{where}: the least-squares W is {w_fit:.9g}, not positive: a sign flip "
                                          "or a negative scale")
                    rel = _relative(w_fit, w_side)
                    self._note("scale_w", rel, where, SCALE_REL)
                    if math.isfinite(w_fit) and w_fit > 0.0 and not (w_side > 0.0 and rel <= SCALE_REL):
                        scale_fail.append(f"{where}: fitted W {w_fit:.12g} against the sidecar's w{ephi} = "
                                          f"{w_side!r} (relative difference {rel:.3g})")
                    ws = w_fit * s
                    ws_max = float(np.max(np.abs(ws)))
                    dev = float(np.max(np.abs(f - ws))) / ws_max if ws_max > 0.0 else math.inf
                    self._note("exit_field", dev, where, self.eps)
                    if not dev <= self.eps:
                        value_fail.append(f"{where}: |fresh - W source| reaches {dev:.3g} of max |W source| "
                                          f"(W = {w_fit:.9g})")
                if not dark:
                    # Independent of W: the declared normalization integral of the stored Ex, Ey (HF-040).
                    n_normalized += 1
                    try:
                        weights = (dataset_v2.node_weights(fresh[r, _WG["Y"]], fresh[r, _WG["Z"]], self.grid)
                                   if self.shape == "rectangle" else np.ones(r.size))
                    except dataset_v2.DatasetV2Error as exc:
                        norm_fail.append(f"{where}: {exc}")
                        continue
                    e_t = dataset_v2.tangential(f)
                    p_fields = float(np.sum(weights * np.sum(np.abs(e_t) ** 2, axis=1))) * h_y * h_z / (2.0 * dataset_v2.ETA0)
                    out = float(fresh[r[0], OUTGOING])
                    if not (math.isfinite(p_fields) and p_fields > 0.0 and out > 0.0):
                        norm_fail.append(f"{where}: a non-dark record needs a positive field integral and "
                                         f"OutgoingPower, got {p_fields!r} W and {out!r} W")
                        continue
                    rel = _relative(p_fields, out)
                    self._note("normalization", rel, where, NORMALIZATION_REL)
                    if not rel <= NORMALIZATION_REL:
                        norm_fail.append(f"{where}: the field integral {p_fields!r} W against OutgoingPower {out!r} W "
                                         f"(relative difference {rel:.3g})")
        n_fitted = len(fitted)
        return [
            _row("exit_field.dark", f"dark records keep W = 1 and the source values within eps = {self.eps:g}; a zero "
                 "source needs a dark record and comes back exactly zero (OR-038)", dark_fail,
                 f"{len(dark_records)} dark records, {len(zero_sources)} zero sources",
                 {"dark_records": dark_records, "zero_source_records": zero_sources}),
            _row("exit_field.scale", f"least-squares W > 0 and equal to the sidecar's w0 / w1 to {SCALE_REL:g}",
                 scale_fail, f"{n_fitted} records", {"fitted": fitted}),
            _row("exit_field.values", f"|fresh - W source| <= {self.eps:g} max |W source| over the six slots, by "
                 "position (CDX-026)", value_fail, f"{n_fitted} records"),
            _row("exit_field.normalization", f"sum w |E_t|^2 h_y h_z / (2 eta0) from the fresh Ex, Ey equals "
                 f"OutgoingPower to {NORMALIZATION_REL:g} on every non-dark record (HF-040)", norm_fail,
                 f"{n_normalized} non-dark records"),
        ]

    # -- 4. far field ------------------------------------------------------------------------------------------------
    def _far_field(self) -> Check:
        problems: list[str] = []
        compared = 0
        for ephi in (0, 1):
            src, fresh = self.source.far_field[ephi], self.fresh.far_field[ephi]
            for k, key in enumerate(self.keys):
                where = _where(ephi, key)
                r = self.ff_rows[ephi][k]
                s, f = src[r, FAR_VALUES], fresh[r, FAR_VALUES]
                compared += 1
                s_max = float(np.max(np.abs(s)))
                if s_max == 0.0:
                    if np.any(f != 0.0):
                        problems.append(f"{where}: a zero source far field must come back exactly zero; the largest "
                                        f"fresh |value| is {float(np.max(np.abs(f)))!r}")
                    continue
                dev = float(np.max(np.abs(f - s))) / s_max
                self._note("far_field", dev, where, self.eps)
                if not dev <= self.eps:
                    problems.append(f"{where}: |fresh - source| reaches {dev:.3g} of max |source|")
        return _row("far_field.values", f"stored far field not rescaled (W = 1): |fresh - source| <= {self.eps:g} "
                    "max |source|, a zero source exactly zero", problems, f"{compared} records")

    # -- 5. rim ------------------------------------------------------------------------------------------------------
    def _rim_counts(self, job: _Job) -> dict[tuple[int, int], int]:
        per_key = _get(job.manifest, "exit_field", "rim", "per_key", where=f"the {job.label} manifest")
        out = {}
        for ephi in (0, 1):
            where = f"the {job.label} manifest's exit_field.rim.per_key[{str(ephi)!r}]"
            found = _strict(per_key.get(str(ephi)) if isinstance(per_key, Mapping) else None, self.keys, where,
                            key_of=_rim_key)
            for k, record in found.items():
                n = _integer(record.get("unsolved_inside_zeroed"), f"{where} at {_key_text(self.keys[k])}")
                if n < 0:
                    raise _Refused(f"{where} at {_key_text(self.keys[k])}: a negative zero-fill count {n}")
                out[(ephi, k)] = n
        return out

    def _rim(self) -> Check:
        if self.shape != "disc":
            return Check("rim.zero_fill", True, "info", f"a {self.shape}: the rim rule (HF-029) applies to discs only")
        source_counts, fresh_counts = self._rim_counts(self.source), self._rim_counts(self.fresh)
        problems: list[str] = []
        totals = {"0": 0, "1": 0}
        for ephi in (0, 1):
            src, fresh = self.source.waveguide[ephi], self.fresh.waveguide[ephi]
            for k, key in enumerate(self.keys):
                where = _where(ephi, key)
                r = self.rows[ephi][k]
                zero = np.all(src[r, FIELDS] == 0.0, axis=1)
                n = int(zero.sum())
                totals[str(ephi)] += n
                if n != source_counts[(ephi, k)]:
                    problems.append(f"{where}: the source holds {n} all-zero exit points, its manifest records "
                                    f"{source_counts[(ephi, k)]} zero-filled points")
                if fresh_counts[(ephi, k)] != source_counts[(ephi, k)]:
                    problems.append(f"{where}: the fresh job zero-filled {fresh_counts[(ephi, k)]} points, the "
                                    f"source {source_counts[(ephi, k)]}")
                back = fresh[r[zero], FIELDS]
                if np.any(back != 0.0):
                    problems.append(f"{where}: {int(np.any(back != 0.0, axis=1).sum())} of the source's {n} "
                                    f"zero-filled points come back nonzero (largest |value| "
                                    f"{float(np.max(np.abs(back)))!r})")
        return _row("rim.zero_fill", "the source's zero-filled points come back exactly zero and the zero-fill counts "
                    "match per key and polarization (HF-029)", problems,
                    f"{totals['0']} (Ephi=0) and {totals['1']} (Ephi=1) zero-filled points", {"zero_filled": totals})

    # -- 6. powers ---------------------------------------------------------------------------------------------------
    def _ai(self) -> float:
        ai = _number(_get(self.fresh.sidecar, "excitation", "incoming_power_w", where="the fresh sidecar"),
                     "the fresh sidecar's excitation.incoming_power_w (A I)")
        if ai <= 0.0:
            raise _Refused(f"A I must be positive, got {ai!r} W")
        return ai

    def _raw_powers(self) -> tuple[dict[tuple[int, int], tuple[float, float]], list[str]]:
        def parse() -> tuple[dict[tuple[int, int], tuple[float, float]], list[str]]:
            where = "the fresh manifest's exit_field.raw_powers_per_key"
            found, problems = _match_records(
                _get(self.fresh.manifest, "exit_field", "raw_powers_per_key", where="the fresh manifest"),
                self.keys, where)
            out = {}
            for k, record in sorted(found.items()):
                polarizations = record.get("polarizations")
                if not isinstance(polarizations, Mapping) or set(polarizations) != {"0", "1"}:
                    problems.append(f"{where} at {_key_text(self.keys[k])}: polarizations must hold exactly '0' and "
                                    "'1'")
                    continue
                for ephi in (0, 1):
                    at = f"the raw powers of {_where(ephi, self.keys[k])}"
                    try:
                        flux = _number(_get(polarizations[str(ephi)], "exit_flux_re_s_dot_n_w", where=at),
                                       f"{at}: exit_flux_re_s_dot_n_w")
                        absolute = _number(_get(polarizations[str(ephi)], "outgoing_power_abs_re_s_w", where=at),
                                           f"{at}: outgoing_power_abs_re_s_w")
                    except _Refused as exc:
                        problems.append(str(exc))
                        continue
                    if absolute < 0.0:
                        problems.append(f"{at}: |Re S| integral {absolute!r} W is negative")
                        continue
                    out[(ephi, k)] = (flux, absolute)
            return out, problems
        return self._cached("raw", parse)

    def _transmittance(self) -> dict[int, dict[str, Any]]:
        def parse() -> dict[int, dict[str, Any]]:
            where = "the fresh sidecar's transmittance.per_key"
            found = _strict(_get(self.fresh.sidecar, "transmittance", "per_key", where="the fresh sidecar"),
                            self.keys, where)
            out = {}
            for k, record in found.items():
                at = f"{where} at {_key_text(self.keys[k])}"
                polarizations = record.get("polarizations")
                if not isinstance(polarizations, Mapping) or set(polarizations) != {"0", "1"}:
                    raise _Refused(f"{at}: polarizations must hold exactly '0' and '1'")
                out[k] = {"c_k": _number(record.get("c_k"), f"{at}: c_k"),
                          "t_flux": {e: _number(_get(polarizations[str(e)], "t_flux", where=at), f"{at}: t_flux[{e}]")
                                     for e in (0, 1)},
                          "abs": {e: _number(_get(polarizations[str(e)], "outgoing_power_abs_re_s_w", where=at),
                                             f"{at}: outgoing_power_abs_re_s_w[{e}]") for e in (0, 1)}}
            return out
        return self._cached("transmittance", parse)

    def _unchecked(self, raw: Mapping) -> str:
        missing = 2 * len(self.keys) - len(raw)
        return f"; {missing} records without a raw power record not checked (powers.raw_records)" if missing else ""

    def _raw_records(self) -> Check:
        raw, problems = self._raw_powers()
        return _row("powers.raw_records", "one finite signed flux and |Re S| per declared key and polarization in the "
                    "fresh manifest's exit_field.raw_powers_per_key", problems, f"{len(raw)} records")

    def _is_dark(self, flux: float) -> bool:
        return abs(flux) <= DARK_REL * self._ai()

    def _dark(self) -> Check:
        raw, _ = self._raw_powers()
        scale = self._scale_records()
        floor = DARK_REL * self._ai()
        problems = []
        for (ephi, k), (flux, _) in sorted(raw.items()):
            computed, declared = abs(flux) <= floor, scale[k][1][ephi]
            if computed != declared:
                problems.append(f"{_where(ephi, self.keys[k])}: the raw flux {flux!r} W is "
                                f"{'' if computed else 'not '}a numerical zero, but dark{ephi} is {declared}")
        n_dark = sum(1 for flux, _ in raw.values() if abs(flux) <= floor)
        return _row("powers.dark", f"dark<p> recomputed as |raw flux| <= {DARK_REL:g} A I = {floor:.6g} W", problems,
                    f"{n_dark} dark of {len(raw)} records{self._unchecked(raw)}")

    def _outgoing(self) -> Check:
        raw, _ = self._raw_powers()
        scale = self._scale_records()
        problems = []
        for (ephi, k), (flux, _) in sorted(raw.items()):
            where = _where(ephi, self.keys[k])
            stored = self.fresh.waveguide[ephi][self.rows[ephi][k], OUTGOING]
            if self._is_dark(flux) or scale[k][1][ephi]:
                if np.any(stored != 0.0):
                    problems.append(f"{where}: a dark record stores OutgoingPower 0 on every row; "
                                    f"{int(np.sum(stored != 0.0))} rows hold up to {float(np.max(np.abs(stored)))!r} W")
            else:
                dev = float(np.max(np.abs(stored - flux)))
                if not dev <= STORED_REL * abs(flux):
                    problems.append(f"{where}: OutgoingPower differs from the raw flux {flux!r} W by up to {dev:.3g} W")
        return _row("powers.outgoing", "stored OutgoingPower 0 on every row of a dark record (OR-038) and the raw flux "
                    f"to {STORED_REL:g} on every row otherwise", problems, f"{len(raw)} records{self._unchecked(raw)}")

    def _t_flux(self) -> Check:
        raw, _ = self._raw_powers()
        scale = self._scale_records()
        per_key = self._transmittance()
        ai = self._ai()
        problems = []
        for (ephi, k), (flux, _) in sorted(raw.items()):
            where = _where(ephi, self.keys[k])
            t = per_key[k]["t_flux"][ephi]
            if self._is_dark(flux) or scale[k][1][ephi]:
                if t != 0.0:
                    problems.append(f"{where}: a dark record's t_flux must be exactly 0, got {t!r}")
                continue
            expected = float(self.fresh.waveguide[ephi][self.rows[ephi][k][0], OUTGOING]) / ai
            if not _relative(t, expected) <= STORED_REL:
                problems.append(f"{where}: t_flux {t!r} against OutgoingPower / (A I) = {expected!r}")
        return _row("powers.t_flux", f"t_flux = stored OutgoingPower / (A I) to {STORED_REL:g}, exactly 0 when dark",
                    problems, f"{len(raw)} records{self._unchecked(raw)}")

    def _ingoing(self) -> Check:
        ai = self._ai()
        problems = []
        for job in (self.source, self.fresh):
            for ephi in (0, 1):
                table = job.waveguide[ephi]
                rel = np.abs(table[:, INGOING] - ai) / ai
                if not np.all(rel <= STORED_REL):
                    i = int(np.argmax(rel))
                    key = (float(table[i, _WG["IWavePhi"]]), float(table[i, _WG["IWaveTheta"]]))
                    problems.append(f"the {job.label} {_where(ephi, key)}: IngoingPower {float(table[i, INGOING])!r} "
                                    f"W against A I = {ai!r} W")
        return _row("powers.ingoing", f"IngoingPower = A I on every row of both jobs, to {STORED_REL:g} (OR-034)",
                    problems, f"A I = {ai!r} W")

    def _intercept_constant(self) -> Check:
        problems = []
        for ephi in (0, 1):
            fresh = self.fresh.waveguide[ephi]
            for k, key in enumerate(self.keys):
                column = fresh[self.rows[ephi][k], INTERCEPT]
                if np.any(column != column[0]):
                    i = int(np.flatnonzero(column != column[0])[0])
                    problems.append(f"{_where(ephi, key)}: InterceptFactor is not constant within the key: row {i} of "
                                    f"the key holds {float(column[i])!r}, its first row {float(column[0])!r}")
        return _row("powers.intercept_factor", "InterceptFactor constant within each key (OR-034)", problems,
                    f"{2 * len(self.keys)} records")

    def _c_key(self) -> Check:
        where = "the fresh sidecar's transmittance.c_key"
        found = _strict(_get(self.fresh.sidecar, "transmittance", "c_key", where="the fresh sidecar"), self.keys, where)
        per_key = self._transmittance()
        problems = []
        for k, key in enumerate(self.keys):
            c = _number(found[k].get("c_k"), f"{where} at {_key_text(key)}: c_k")
            recomputed = dataset_v2.c_key(key[0], key[1], self.theta_list)
            rel = _relative(c, recomputed)
            self._note("intercept_factor", rel, _key_text(key), INTERCEPT_REL)
            if not rel <= INTERCEPT_REL:
                problems.append(f"at {_key_text(key)}: the sidecar's c_key {c!r} against {recomputed!r} recomputed "
                                "from the declared grid")
            if not _relative(per_key[k]["c_k"], c) <= INTERCEPT_REL:
                problems.append(f"at {_key_text(key)}: transmittance.per_key c_k {per_key[k]['c_k']!r} against "
                                f"c_key {c!r}")
            for ephi in (0, 1):
                first = float(self.fresh.waveguide[ephi][self.rows[ephi][k][0], INTERCEPT])
                rel = _relative(first, c)
                self._note("intercept_factor", rel, _where(ephi, key), INTERCEPT_REL)
                if not rel <= INTERCEPT_REL:
                    problems.append(f"{_where(ephi, key)}: InterceptFactor {first!r} against the sidecar's c_key "
                                    f"{c!r}")
        return _row("powers.c_key", f"InterceptFactor and transmittance c_k equal to c_key, itself equal to |cos theta| "
                    f"or cell_90 recomputed from the declared grid, to {INTERCEPT_REL:g} (OR-034)", problems,
                    f"{len(self.keys)} keys, cell_90 = {dataset_v2.cell_90(self.theta_list)!r}")

    def _abs_re_s(self) -> Check:
        raw, _ = self._raw_powers()
        per_key = self._transmittance()
        floor = DARK_REL * self._ai()
        problems = []
        for (ephi, k), (_, absolute) in sorted(raw.items()):
            where = _where(ephi, self.keys[k])
            source = float(self.source.waveguide[ephi][self.rows[ephi][k][0], OUTGOING])
            diff = abs(absolute - source)
            self._note("abs_re_s", _relative(absolute, source) if diff > floor else 0.0, where, ABS_RE_S_REL)
            if not diff <= max(ABS_RE_S_REL * abs(source), floor):
                problems.append(f"{where}: the fresh |Re S| {absolute!r} W against the source's OutgoingPower "
                                f"{source!r} W (relative difference {_relative(absolute, source):.3g})")
            side = per_key[k]["abs"][ephi]
            if not abs(side - absolute) <= STORED_REL * abs(absolute):
                problems.append(f"{where}: the sidecar's outgoing_power_abs_re_s_w {side!r} W is not the manifest's "
                                f"raw |Re S| {absolute!r} W")
        return _row("powers.abs_re_s", f"the fresh |Re S| (raw_powers_per_key) equals the source's OutgoingPower to "
                    f"{ABS_RE_S_REL:g} relative (absolute {DARK_REL:g} A I at zero), as step 8b compared them; a "
                    "changed definition of either shows here", problems, f"{len(raw)} records{self._unchecked(raw)}")

    def _poynting(self) -> Check:
        rule = (f"the fresh raw flux equals step 8b's flux_z_w to {POYNTING_REL:g} relative (absolute {DARK_REL:g} A I "
                "at zero), every key and polarization")
        if self.poynting is None:
            return Check("powers.poynting", True, "info", "not requested (--poynting)")
        data = _read_json(self.poynting, "the --poynting record")
        raw, _ = self._raw_powers()
        floor = DARK_REL * self._ai()
        problems = []
        source_sha = sha256_file(self.source.root / "manifest.json")
        if data.get("manifest_sha256") != source_sha:
            problems.append(f"the step-8b record was made on a job whose manifest_sha256 is "
                            f"{data.get('manifest_sha256')!r}, not this source's {source_sha}")
        records = data.get("keys")
        if not isinstance(records, list):
            raise _Refused(f"the --poynting record's keys must be a list, got {type(records).__name__}")
        seen: dict[tuple[int, int], float] = {}
        for i, record in enumerate(records):
            at = f"the step-8b record keys[{i}]"
            try:
                ephi = _integer(_get(record, "ephi", where=at), f"{at}.ephi")
                phi, theta = _phi_theta(record, at)
                f8b = _number(_get(record, "exit", "flux_z_w", where=at), f"{at} ({_where(ephi, (phi, theta))}) "
                              "exit.flux_z_w")
            except _Refused as exc:
                problems.append(str(exc))
                continue
            k = _key_index(self.keys, phi, theta)
            if ephi not in (0, 1) or k is None:
                problems.append(f"{at}: Ephi={ephi} at {_key_text((phi, theta))} is not a declared record")
            elif (ephi, k) in seen:
                problems.append(f"{at}: a second record of {_where(ephi, self.keys[k])}")
            else:
                seen[(ephi, k)] = f8b
        missing = [(e, k) for e in (0, 1) for k in range(len(self.keys)) if (e, k) not in seen]
        if missing:
            problems.append(f"the step-8b record lacks {len(missing)} records, the first "
                            f"{_where(missing[0][0], self.keys[missing[0][1]])}")
        for (ephi, k), f8b in sorted(seen.items()):
            if (ephi, k) not in raw:
                continue
            flux = raw[(ephi, k)][0]
            diff = abs(flux - f8b)
            where = _where(ephi, self.keys[k])
            self._note("poynting", _relative(flux, f8b) if diff > floor else 0.0, where, POYNTING_REL)
            if not diff <= max(POYNTING_REL * abs(f8b), floor):
                problems.append(f"{where}: the raw flux {flux!r} W against step 8b's flux_z_w {f8b!r} W (relative "
                                f"difference {_relative(flux, f8b):.3g})")
        return _row("powers.poynting", rule, problems, f"{len(seen)} records",
                    {"step_8b_status": data.get("status"), "step_8b_manifest_sha256": data.get("manifest_sha256")})

    # -- 7. convergence ----------------------------------------------------------------------------------------------
    @staticmethod
    def _differences(a: object, b: object) -> list[str]:
        a = a if isinstance(a, Mapping) else {}
        b = b if isinstance(b, Mapping) else {}
        return [f"{name} {a.get(name, 'missing')!r} vs {b.get(name, 'missing')!r}" for name in QUALIFICATION_FIELDS
                if not _same_value(a.get(name, _MISSING), b.get(name, _MISSING))]

    def _convergence(self) -> Check:
        src_m, fresh_m = self.source.manifest, self.fresh.manifest
        target = _number(_get(src_m, "solver", "max_delta_e", where="the source manifest"),
                         "the source manifest's solver.max_delta_e")
        if target <= 0.0:
            raise _Refused(f"the source manifest's solver.max_delta_e must be positive, got {target!r}")
        problems = []
        fresh_target = _get(fresh_m, "solver", "max_delta_e", where="the fresh manifest", default=None)
        if not _same_value(fresh_target, target, rel=0.0):
            problems.append(f"the fresh manifest's solver.max_delta_e {fresh_target!r} is not the source's {target!r}")
        source_recorded = _get(src_m, "solver", "convergence_qualification", where="", default=None)
        fresh_recorded = _get(fresh_m, "solver", "convergence_qualification", where="the fresh manifest")
        tables, source_q, fresh_q = {}, {}, {}
        for ephi in (0, 1):
            e = str(ephi)
            path = self.source.root / "logs" / f"convergence_Ephi{ephi}.txt"
            try:
                text = path.read_text()
            except OSError as exc:
                problems.append(f"Ephi={ephi}: the source convergence table {path} is unreadable ({exc}): it is no "
                                "evidence of convergence")
                continue
            tables[e] = {"path": str(path), "sha256": sha256_file(path)}
            q = qualify(text, target)
            source_q[e] = {name: q[name] for name in QUALIFICATION_FIELDS}
            if not q["qualified"]:
                problems.append(f"Ephi={ephi}: the source table does not qualify at target {target:g}: "
                                + "; ".join(q["reasons"]))
            if source_recorded is not None:
                record = source_recorded.get(e) if isinstance(source_recorded, Mapping) else None
                diff = self._differences(record, source_q[e])
                if diff:
                    problems.append(f"Ephi={ephi}: the source manifest's solver.convergence_qualification disagrees "
                                    f"with its table: {', '.join(diff)}")
            record = fresh_recorded.get(e) if isinstance(fresh_recorded, Mapping) else None
            diff = self._differences(source_q[e], record)
            if diff:
                problems.append(f"Ephi={ephi}: the source qualification recomputed from its table against the fresh "
                                f"manifest's: {', '.join(diff)}")
            fresh_path = self.fresh.root / "logs" / f"convergence_Ephi{ephi}.txt"
            try:
                fq = qualify(fresh_path.read_text(), target)
            except OSError as exc:
                problems.append(f"Ephi={ephi}: the fresh convergence table {fresh_path} is unreadable ({exc})")
                continue
            fresh_q[e] = {name: fq[name] for name in QUALIFICATION_FIELDS}
            diff = self._differences(fresh_q[e], record)
            if diff:
                problems.append(f"Ephi={ephi}: the fresh table recomputed against the fresh manifest's record: "
                                f"{', '.join(diff)}")
        recorded = "the source manifest's record" if source_recorded is not None else "(the source manifest records none)"
        return _row("convergence", f"the source tables re-qualified at its target {target:g} (CDX-025): both "
                    "polarizations qualify and agree field by field with the fresh qualification", problems,
                    f"both qualify; agree with the fresh tables and record, and with {recorded}",
                    {"target": target, "source_tables": tables, "source_recomputed": source_q,
                     "fresh_recomputed": fresh_q, "source_manifest_has_qualification": source_recorded is not None})


# --- entry points ----------------------------------------------------------------------------------------------------

def check(source_job: Path, fresh_job: Path, poynting: Path | None = None, eps: float = EPS) -> ReexportResult:
    """Run the acceptance of ``fresh_job`` (schema 2) against ``source_job`` (schema 1) and write
    ``<fresh_job>/reexport_check.json``. ``poynting`` is step 8b's poynting.json for the source job, or None.

    Raises ``FileNotFoundError`` when a job directory or the poynting file is missing (a usage error); every
    problem with the jobs' contents becomes a failing row instead.
    """
    source_job, fresh_job = Path(source_job), Path(fresh_job)
    for label, path in (("source", source_job), ("fresh", fresh_job)):
        if not path.is_dir():
            raise FileNotFoundError(f"the {label} job directory {path} does not exist")
    if poynting is not None and not Path(poynting).is_file():
        raise FileNotFoundError(f"the --poynting file {poynting} does not exist")
    if not (isinstance(eps, (int, float)) and math.isfinite(eps) and eps > 0.0):
        raise ValueError(f"eps must be a finite positive number, got {eps!r}")
    acceptance = _Acceptance(source_job, fresh_job, None if poynting is None else Path(poynting), float(eps))
    checks = acceptance.run()
    tolerances = {"eps": float(eps), "scale_rel": SCALE_REL, "normalization_rel": NORMALIZATION_REL,
                  "intercept_rel": INTERCEPT_REL, "abs_re_s_rel": ABS_RE_S_REL, "poynting_rel": POYNTING_REL,
                  "stored_rel": STORED_REL, "coord_m": COORD_M, "key_deg": KEY_DEG, "dark_rel": DARK_REL}
    result = ReexportResult(str(source_job), str(fresh_job), None if poynting is None else str(poynting), checks,
                            acceptance.largest, tolerances)
    write_json(fresh_job / REPORT_NAME, result.to_dict())
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m bbsim.reexport_check", description=__doc__.splitlines()[0])
    parser.add_argument("--source", type=Path, required=True, help="the schema-1 source job directory")
    parser.add_argument("--fresh", type=Path, required=True,
                        help="the schema-2 job run_frequency --from-solved made from it")
    parser.add_argument("--poynting", type=Path, default=None,
                        help="step 8b's diagnose_poynting/poynting.json of the source job")
    try:
        args = parser.parse_args(argv)
    except SystemExit as exc:             # argparse has printed its message
        return 2 if exc.code else 0
    try:
        result = check(args.source, args.fresh, args.poynting)
    except FileNotFoundError as exc:
        print(f"reexport_check: {exc}", file=sys.stderr)
        return 2
    print(result.format_table())
    return 0 if result.passed else 1


if __name__ == "__main__":
    sys.exit(main())
