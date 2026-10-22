"""Key-by-key acceptance of a 35-key solve against the 15-key solve of the same geometry (OR-036, HPC step 10b).

Usage::

    python -m bbsim.grid_acceptance --new JOB --reference JOB [--t-rel-tol 0.005]

``--new`` is a finished schema-2 job on OR-036's incidence grid (theta 90 to 180 deg in 15 deg steps, phi 0 to 90 deg
in 22.5 deg steps, both in the crack-local frame: 35 keys); ``--reference`` is the schema-1 job of the same dataset and
frequency on the 15-key grid (theta 0 to 180 deg and phi 0 to 90 deg in 45 deg steps): crack1 and crack2
``job_tight1`` (step 5d), the round gap ``job_reference1`` (step 8). The two grids share nine keys, ``SHARED_KEYS``,
phi {0, 45, 90} x theta {90, 135, 180} deg; this checks that the new solve reproduces the old one there, so that its
26 new keys can be trusted.

The compared quantity is the legacy magnitude-power ratio on both sides, not the schema-2 photon transmission:
T_new = outgoing_power_abs_re_s_w / IngoingPower, the |Re S| integral the new sidecar keeps in
``transmittance.per_key`` (matched by key value, not by position) over the new table's IngoingPower at the key;
T_ref = OutgoingPower / IngoingPower from the reference table's first row at the key, as schema 1 defines T. Neither
the stored schema-2 OutgoingPower (the exit flux) nor its InterceptFactor c_k enters: the flux differs from |Re S| by
0.2 % at normal incidence and 10 to 18 % on oblique keys (step 8b, CDX-009), and c_k at theta = 90 deg follows the
declared grid (cell_90 0.0838 on 15 keys, 0.0366 on 35), so either would hide the mesh and solver differences this
measures behind the decided numerator and divisor changes. The tolerance, 0.5 % of T_ref, sits above the mesh noise
measured at MaxDeltaE 0.005 (0.1 % between separately meshed polarizations in step 5d, 0.3 % for the round gap's
2.5 deg facet check in step 8).

Rows, in order (the first stops the rest when it fails, as nothing else can be read in the right way):

1. schema.declared: each job holds exactly one ``*.dataset.json`` sidecar; the new one declares schema 2.x, the
   reference 1.x.
2. identity: the same dataset_id and frequency_label on both.
3. grid.new: the new sidecar declares exactly OR-036's product grid, 5 phi x 7 theta = 35 keys.
4. grid.reference: the reference sidecar declares every shared key.
5. tables: both polarizations' ``waveguide.csv`` of both jobs readable, their leading columns Freq, Ephi, IWavePhi,
   IWaveTheta, OutgoingPower, IngoingPower (the same in both schemas, read by position), every key angle finite.
6. powers.ingoing: IngoingPower finite, positive and constant within each table, and equal between the new and the
   reference tables, to 1e-9 relative (same geometry, same Ei: the key's constant A I).
7. transmittance.shared_keys: per shared key and polarization, where T_ref >= 1e-3 x the polarization's largest
   T_ref over the shared keys, |T_new - T_ref| <= t_rel_tol x T_ref; keys below that band are reported as noise and
   do not fail (crack1's E-along-slot polarization at normal incidence carries T about 1e-10); a noise row whose
   T_new lies inside the band carries a note saying so.

Every number read must be finite: a missing, malformed or non-finite value fails the row that reads it, naming the key
and polarization where there is one, and never escapes as an exception (CDX-030: nothing is compared with NaN). The
rows, the per-key table and the verdict go to ``<new>/grid_acceptance.json`` (written only into a job whose sidecar
declares schema 2, never into a 15-key source job), after the table is printed; a report that cannot be written is
named on its own ``report:`` line and fails the run, without a traceback. The printed table ends in
``GRID35: PASS`` or ``GRID35: FAIL <n>`` (n failing rows). Exit code 0 on PASS, 1 on FAIL, 2 when a job directory is
missing or the arguments are wrong.
"""
from __future__ import annotations

import argparse
import math
import sys
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from bbsim import dataset_v2
from bbsim.acceptance_rows import RowVerdicts, guard
from bbsim.acceptance_rows import Refused as _Refused
from bbsim.acceptance_rows import at_key as _where
from bbsim.acceptance_rows import get as _get
from bbsim.acceptance_rows import number as _number
from bbsim.acceptance_rows import read_json as _read_json
from bbsim.acceptance_rows import row as _row
from bbsim.compare import Check
from bbsim.dataset_v2 import key_text
from bbsim.manifest import write_json
from bbsim.naming import dataset_dir_name

REPORT_NAME = "grid_acceptance.json"
# The nine keys both grids declare: phi in {0, 45, 90} and theta in {90, 135, 180} deg, all reachable.
SHARED_KEYS = [(p, t) for p in (0.0, 45.0, 90.0) for t in (90.0, 135.0, 180.0)]
T_REL_TOL = 0.005          # |T_new - T_ref| <= T_REL_TOL * T_ref on the compared keys
NOISE_REL = 1e-3           # a key whose T_ref is below NOISE_REL * the polarization's largest T_ref is noise
INGOING_REL = 1e-9         # IngoingPower within a table and between the jobs
GRID_PHI_DEG = (0.0, 22.5, 45.0, 67.5, 90.0)                              # OR-036
GRID_THETA_DEG = (90.0, 105.0, 120.0, 135.0, 150.0, 165.0, 180.0)        # OR-036
KEY_DEG = dataset_v2.KEY_TOL_DEG
NEW_MAJOR, REFERENCE_MAJOR = "2", "1"
# The waveguide table's first six columns, the same labels and positions in schema 1.0 and 2.0.
LEADING_COLUMNS = ["Freq", "Ephi", "IWavePhi", "IWaveTheta", "OutgoingPower", "IngoingPower"]
PHI, THETA, OUTGOING, INGOING = 2, 3, 4, 5


# --- small helpers ---------------------------------------------------------------------------------------------------

def _axis(sidecar: Mapping, name: str, where: str) -> list[float]:
    values = _get(sidecar, "excitation", name, where=where)
    if not isinstance(values, list) or not values:
        raise _Refused(f"{where}'s excitation.{name} must be a non-empty list, got {values!r}")
    numbers = [_number(v, f"{where}'s excitation.{name}") for v in values]
    ordered = sorted(numbers)
    if any(b - a <= KEY_DEG for a, b in zip(ordered, ordered[1:])):
        raise _Refused(f"{where}'s excitation.{name} repeats a value: {numbers}")
    return ordered


def _same_axis(declared: Sequence[float], expected: Sequence[float]) -> bool:
    return len(declared) == len(expected) and all(abs(a - b) <= KEY_DEG for a, b in zip(declared, expected))


def _holds(axis: Sequence[float], value: float) -> bool:
    return any(abs(a - value) <= KEY_DEG for a in axis)


def _finite_or_none(value: float | None) -> float | None:
    return float(value) if value is not None and math.isfinite(value) else None


# --- the result ------------------------------------------------------------------------------------------------------

@dataclass
class GridResult(RowVerdicts):
    """The rows in order, one comparison record per shared key and polarization (``phi_deg``, ``theta_deg``,
    ``ephi``, ``t_new``, ``t_ref``, ``rel_diff`` = (T_new - T_ref) / T_ref, ``status`` PASS, FAIL or noise, ``note``),
    the tolerance used, where the JSON report goes and whether it was written, and the verdict. A report that could
    not be written fails the run."""

    new_job: str
    reference_job: str
    t_rel_tol: float
    checks: list[Check]
    rows: list[dict[str, Any]] = field(default_factory=list)
    report_target: str | None = None   # <new>/grid_acceptance.json for a declared schema-2 new job, else None
    report: str | None = None          # the JSON written (``write_report``)
    report_error: str | None = None    # why it could not be written

    @property
    def passed(self) -> bool:
        return super().passed and self.report_error is None

    @property
    def verdict(self) -> str:
        if self.passed:
            return "GRID35: PASS"
        return f"GRID35: FAIL {len(self.failures) + (1 if self.report_error else 0)}"

    def to_dict(self) -> dict[str, Any]:
        return {"new_job": self.new_job, "reference_job": self.reference_job, "passed": self.passed,
                "failed": len(self.failures), "verdict": self.verdict,
                "tolerances": {"t_rel_tol": self.t_rel_tol, "noise_rel": NOISE_REL, "ingoing_rel": INGOING_REL,
                               "key_deg": KEY_DEG},
                "shared_keys": [list(k) for k in SHARED_KEYS],
                "checks": [asdict(c) for c in self.checks], "rows": [dict(r) for r in self.rows]}

    def format_table(self) -> str:
        return self.format_rows() + "\n" + self.format_footer()

    def format_footer(self) -> str:
        """The report line and the verdict: what ``main`` prints after writing the report."""
        if self.report:
            report = f"report: {self.report}"
        elif self.report_error:
            report = f"report: {self.report_error}"
        else:
            report = "report: not written (the new job is not a declared schema-2 job; a 15-key source is never written)"
        return "\n".join(["", report, "", self.verdict])

    def format_rows(self) -> str:
        """The header, the check rows and the per-key table: what ``main`` prints before writing the report."""
        lines = [f"new:       {self.new_job}", f"reference: {self.reference_job}",
                 f"tolerances: t_rel_tol={self.t_rel_tol:g} noise_rel={NOISE_REL:g} ingoing_rel={INGOING_REL:g}", ""]
        for c in self.checks:
            status = ("PASS" if c.passed else "FAIL") if c.severity == "fail" else ("INFO" if c.passed else "SKIP")
            lines.append(f"{status:4} {c.name:25} {c.detail}")
        if self.rows:
            lines += ["", f"{'phi':>6} {'theta':>6} {'Ephi':>4} {'T_new':>16} {'T_ref':>16} {'rel_diff':>10}  status"]
            for r in self.rows:
                t_new = "-" if r["t_new"] is None else f"{r['t_new']:.10g}"
                t_ref = "-" if r["t_ref"] is None else f"{r['t_ref']:.10g}"
                rel = "-" if r["rel_diff"] is None else f"{r['rel_diff']:+.4%}"
                note = f"  ({r['note']})" if r.get("note") else ""
                lines.append(f"{r['phi_deg']:>6g} {r['theta_deg']:>6g} {r['ephi']:>4d} {t_new:>16} {t_ref:>16} "
                             f"{rel:>10}  {r['status']}{note}")
        return "\n".join(lines)


# --- the acceptance --------------------------------------------------------------------------------------------------

@dataclass
class _Job:
    label: str
    root: Path
    sidecar: dict = field(default_factory=dict)
    # Per Ephi, the columns this acceptance reads, as floats: phi, theta, IngoingPower, and (reference only)
    # OutgoingPower. The new job's stored OutgoingPower, the schema-2 flux, is never read.
    tables: dict[int, dict[str, np.ndarray]] = field(default_factory=dict)

    def waveguide(self, ephi: int) -> Path:
        where = f"the {self.label} sidecar"
        dataset_id = _get(self.sidecar, "dataset_id", where=where)
        if not isinstance(dataset_id, str) or not dataset_id:
            raise _Refused(f"{where}'s dataset_id must be a non-empty string, got {dataset_id!r}")
        frequency = _number(_get(self.sidecar, "frequency_ghz", where=where), f"{where}'s frequency_ghz")
        return self.root / dataset_dir_name(dataset_id, frequency, ephi) / "waveguide.csv"


class _Acceptance:
    def __init__(self, new: Path, reference: Path, t_rel_tol: float):
        self.new, self.reference = _Job("new", new), _Job("reference", reference)
        self.tol = t_rel_tol
        self.rows: list[dict[str, Any]] = []

    @staticmethod
    def _guard(name: str, step: Any) -> Check:
        return guard([name], lambda: [step()])[0]     # a refusal is a failing row, never a traceback (CDX-030)

    def run(self) -> list[Check]:
        checks = [self._guard("schema.declared", self._schemas)]
        if not checks[0].passed:
            return checks
        for name, step in (("identity", self._identity), ("grid.new", self._grid_new),
                           ("grid.reference", self._grid_reference), ("tables", self._tables),
                           ("powers.ingoing", self._ingoing), ("transmittance.shared_keys", self._transmittance)):
            checks.append(self._guard(name, step))
        return checks

    # -- 1. declared schemas -------------------------------------------------------------------------------------
    def _schemas(self) -> Check:
        problems: list[str] = []
        declared: dict[str, object] = {}
        for job, major in ((self.new, NEW_MAJOR), (self.reference, REFERENCE_MAJOR)):
            sidecars = sorted(job.root.glob("*.dataset.json"))
            if len(sidecars) != 1:
                names = f": {[p.name for p in sidecars]}" if sidecars else ""
                problems.append(f"the {job.label} job holds {len(sidecars)} dataset sidecars (*.dataset.json), "
                                f"not 1{names}")
                continue
            try:
                job.sidecar = _read_json(sidecars[0], f"the {job.label} sidecar")
            except _Refused as exc:
                problems.append(str(exc))
                continue
            version = job.sidecar.get("schema_version")
            declared[job.label] = version
            if not (isinstance(version, str) and version.split(".")[0] == major):
                problems.append(f"the {job.label} sidecar {sidecars[0].name} declares schema_version {version!r}; "
                                f"the {job.label} job must be schema {major}.x")
        return _row("schema.declared", f"the new job declares schema {NEW_MAJOR}.x, the reference schema "
                    f"{REFERENCE_MAJOR}.x, each in its one sidecar", problems,
                    f"new {declared.get('new')!r}, reference {declared.get('reference')!r}", {"declared": declared})

    # -- 2. identity ---------------------------------------------------------------------------------------------
    def _identity(self) -> Check:
        problems = []
        for name in ("dataset_id", "frequency_label"):
            a, b = self.new.sidecar.get(name), self.reference.sidecar.get(name)
            if not isinstance(a, str) or not a or a != b:
                problems.append(f"{name}: new {a!r}, reference {b!r}")
        return _row("identity", "the same dataset_id and frequency_label on both", problems,
                    f"{self.new.sidecar.get('dataset_id')} at {self.new.sidecar.get('frequency_label')}")

    # -- 3, 4. declared grids ------------------------------------------------------------------------------------
    def _grid_new(self) -> Check:
        where = "the new sidecar"
        phis, thetas = _axis(self.new.sidecar, "incident_phi_deg", where), _axis(self.new.sidecar,
                                                                                  "incident_theta_deg", where)
        problems = []
        if not (_same_axis(phis, GRID_PHI_DEG) and _same_axis(thetas, GRID_THETA_DEG)):
            problems.append(f"the new job declares {len(phis)} x {len(thetas)} = {len(phis) * len(thetas)} keys, phi "
                            f"{phis} x theta {thetas} deg, not OR-036's 35")
        return _row("grid.new", f"the new job declares OR-036's 35-key grid, phi {list(GRID_PHI_DEG)} x theta "
                    f"{list(GRID_THETA_DEG)} deg", problems, f"{len(phis) * len(thetas)} keys")

    def _grid_reference(self) -> Check:
        where = "the reference sidecar"
        phis = _axis(self.reference.sidecar, "incident_phi_deg", where)
        thetas = _axis(self.reference.sidecar, "incident_theta_deg", where)
        problems = [f"the reference does not declare {key_text(key)}" for key in SHARED_KEYS
                    if not (_holds(phis, key[0]) and _holds(thetas, key[1]))]
        return _row("grid.reference", "the reference job declares all nine shared keys", problems,
                    f"{len(phis) * len(thetas)} keys declared, the nine shared among them")

    # -- 5. tables -----------------------------------------------------------------------------------------------
    def _tables(self) -> Check:
        problems = []
        for job in (self.new, self.reference):
            columns = {"phi": PHI, "theta": THETA, "in": INGOING}
            if job is self.reference:
                columns["out"] = OUTGOING
            for ephi in (0, 1):
                at = f"the {job.label} Ephi={ephi} waveguide.csv"
                try:
                    path = job.waveguide(ephi)
                    header = list(pd.read_csv(path, nrows=0).columns)
                    if header[:len(LEADING_COLUMNS)] != LEADING_COLUMNS:
                        raise _Refused(f"{at} ({path}) begins {header[:len(LEADING_COLUMNS)]}, not "
                                       f"{LEADING_COLUMNS}")
                    raw = pd.read_csv(path, usecols=sorted(columns.values()), header=0)
                except _Refused as exc:
                    problems.append(str(exc))
                    continue
                except (OSError, ValueError) as exc:
                    problems.append(f"{at} is unreadable: {exc}")
                    continue
                values = {name: pd.to_numeric(raw.iloc[:, sorted(columns.values()).index(pos)],
                                              errors="coerce").to_numpy(dtype=float)
                          for name, pos in columns.items()}
                bad = ~(np.isfinite(values["phi"]) & np.isfinite(values["theta"]))
                if bad.any():
                    problems.append(f"{at}: {int(bad.sum())} rows with a non-finite IWavePhi or IWaveTheta, the "
                                    f"first at row {int(np.flatnonzero(bad)[0])}")
                    continue
                if not len(values["phi"]):
                    problems.append(f"{at} has no rows")
                    continue
                job.tables[ephi] = values
        summary = "" if problems else "; ".join(
            f"{job.label} {len(job.tables[0]['phi'])} / {len(job.tables[1]['phi'])} rows (Ephi=0 / 1)"
            for job in (self.new, self.reference))
        return _row("tables", "both polarizations' waveguide.csv of both jobs readable, leading columns "
                    f"{LEADING_COLUMNS}, every key angle finite", problems, summary)

    def _key_rows(self, job: _Job, ephi: int, key: tuple[float, float]) -> np.ndarray:
        table = job.tables[ephi]
        return np.flatnonzero((np.abs(table["phi"] - key[0]) <= KEY_DEG) & (np.abs(table["theta"] - key[1]) <= KEY_DEG))

    # -- 6. IngoingPower ------------------------------------------------------------------------------------------
    def _ingoing(self) -> Check:
        problems = []
        constants: dict[tuple[str, int], float] = {}
        for job in (self.reference, self.new):
            for ephi in (0, 1):
                at = f"the {job.label} Ephi={ephi} table"
                if ephi not in job.tables:
                    problems.append(f"{at} is not usable (tables)")
                    continue
                table = job.tables[ephi]
                values = table["in"]
                bad = ~np.isfinite(values) | ~(values > 0.0)
                if bad.any():
                    i = int(np.flatnonzero(bad)[0])
                    key = (float(table["phi"][i]), float(table["theta"][i]))
                    problems.append(f"{at}: IngoingPower {float(values[i])!r} W at {key_text(key)} is not a finite "
                                    "positive number")
                    continue
                first = float(values[0])
                spread = np.abs(values - first) / first
                if not np.all(spread <= INGOING_REL):
                    i = int(np.argmax(spread))
                    key = (float(table["phi"][i]), float(table["theta"][i]))
                    problems.append(f"{at}: IngoingPower {float(values[i])!r} W at {key_text(key)} differs from its "
                                    f"first row's {first!r} W by {float(spread[i]):.3g} relative, beyond "
                                    f"{INGOING_REL:g}: not constant")
                    continue
                constants[(job.label, ephi)] = first
        if constants:
            base_label, base = next(iter(constants.items()))
            for (label, ephi), value in constants.items():
                rel = abs(value - base) / base
                if not rel <= INGOING_REL:
                    problems.append(f"the {label} Ephi={ephi} IngoingPower {value!r} W against the {base_label[0]} "
                                    f"Ephi={base_label[1]} {base!r} W: {rel:.3g} relative, beyond {INGOING_REL:g}")
        return _row("powers.ingoing", f"IngoingPower finite, positive and constant within each table and equal between "
                    f"the jobs, to {INGOING_REL:g} relative (the key's constant A I)", problems,
                    f"A I = {next(iter(constants.values()))!r} W in all four tables" if constants else "")

    # -- 7. the comparison ---------------------------------------------------------------------------------------
    def _new_records(self) -> tuple[dict[int, Mapping], dict[int, str], list[str]]:
        """The new sidecar's transmittance.per_key records at the shared keys, by shared key index; the shared keys
        whose record is unusable, with why; and what is wrong with the list itself."""
        where = "the new sidecar's transmittance.per_key"
        records = _get(self.new.sidecar, "transmittance", "per_key", where="the new sidecar")
        if not isinstance(records, list):
            raise _Refused(f"{where} must be a list of per-key records, got {type(records).__name__}")
        found: dict[int, Mapping] = {}
        broken: dict[int, str] = {}
        problems: list[str] = []
        for i, record in enumerate(records):
            at = f"{where}[{i}]"
            if not isinstance(record, Mapping):
                problems.append(f"{at} is not an object")
                continue
            try:
                key = (_number(record.get("phi_deg"), f"{at}.phi_deg"),
                       _number(record.get("theta_deg"), f"{at}.theta_deg"))
            except _Refused as exc:
                problems.append(str(exc))
                continue
            k = dataset_v2.key_position(SHARED_KEYS, *key)
            if k is None:
                continue                       # one of the 26 keys the 15-key grid does not have
            if k in found or k in broken:
                broken[k] = f"{where} holds two records at {key_text(SHARED_KEYS[k])}"
                found.pop(k, None)
            else:
                found[k] = record
        for k in range(len(SHARED_KEYS)):
            if k not in found and k not in broken:
                broken[k] = f"{where} has no record at {key_text(SHARED_KEYS[k])}"
        return found, broken, problems

    def _t_ref(self, ephi: int, key: tuple[float, float]) -> float:
        job = self.reference
        if ephi not in job.tables:
            raise _Refused(f"the reference Ephi={ephi} table is not usable (tables)")
        rows = self._key_rows(job, ephi, key)
        if not rows.size:
            raise _Refused(f"the reference Ephi={ephi} table has no rows at {key_text(key)}")
        out, inn = job.tables[ephi]["out"][rows], job.tables[ephi]["in"][rows]
        if not np.all(np.isfinite(out)):
            raise _Refused(f"the reference OutgoingPower is not finite ({float(out[~np.isfinite(out)][0])!r} W)")
        if not (math.isfinite(inn[0]) and inn[0] > 0.0):
            raise _Refused(f"the reference IngoingPower {float(inn[0])!r} W is not a finite positive number")
        if out[0] < 0.0:
            raise _Refused(f"the reference OutgoingPower {float(out[0])!r} W is negative")
        return float(out[0]) / float(inn[0])

    def _t_new(self, ephi: int, key: tuple[float, float], record: Mapping) -> float:
        at = f"the new sidecar's transmittance.per_key at {key_text(key)}"
        polarizations = record.get("polarizations")
        if not isinstance(polarizations, Mapping) or str(ephi) not in polarizations:
            raise _Refused(f"{at} has no polarization {str(ephi)!r}")
        absolute = _number(_get(polarizations[str(ephi)], "outgoing_power_abs_re_s_w", where=at),
                           f"{at}: outgoing_power_abs_re_s_w[{ephi}]")
        if absolute < 0.0:
            raise _Refused(f"{at}: outgoing_power_abs_re_s_w[{ephi}] {absolute!r} W is negative")
        if ephi not in self.new.tables:
            raise _Refused(f"the new Ephi={ephi} table is not usable (tables)")
        rows = self._key_rows(self.new, ephi, key)
        if not rows.size:
            raise _Refused(f"the new Ephi={ephi} table has no rows at {key_text(key)}")
        inn = float(self.new.tables[ephi]["in"][rows[0]])
        if not (math.isfinite(inn) and inn > 0.0):
            raise _Refused(f"the new IngoingPower {inn!r} W is not a finite positive number")
        return absolute / inn

    def _transmittance(self) -> Check:
        records, broken, problems = self._new_records()
        values: dict[tuple[int, int], tuple[float | None, float | None, str | None]] = {}
        for k, key in enumerate(SHARED_KEYS):
            for ephi in (0, 1):
                t_new = t_ref = None
                notes = []
                try:
                    t_ref = self._t_ref(ephi, key)
                except _Refused as exc:
                    notes.append(str(exc))
                if k in broken:
                    notes.append(broken[k])
                else:
                    try:
                        t_new = self._t_new(ephi, key, records[k])
                    except _Refused as exc:
                        notes.append(str(exc))
                values[(k, ephi)] = (t_new, t_ref, "; ".join(notes) or None)
        largest = {}
        for ephi in (0, 1):
            refs = [v[1] for (k, e), v in values.items() if e == ephi and v[1] is not None]
            largest[ephi] = max(refs) if refs and max(refs) > 0.0 else None
            if largest[ephi] is None:
                problems.append(f"Ephi={ephi}: no positive T_ref among the shared keys, so no key can be compared")
        failures: list[str] = list(problems)
        counts = {"PASS": 0, "FAIL": 0, "noise": 0}
        noted = 0
        self.rows = []
        for k, key in enumerate(SHARED_KEYS):
            for ephi in (0, 1):
                t_new, t_ref, note = values[(k, ephi)]
                rel = (t_new - t_ref) / t_ref if t_new is not None and t_ref is not None and t_ref != 0.0 else None
                if note is not None:
                    status = "FAIL"
                    failures.append(f"{_where(ephi, key)}: {note}")
                elif largest[ephi] is None:
                    status, note = "FAIL", "no positive T_ref for this polarization"
                elif t_ref >= NOISE_REL * largest[ephi]:
                    if abs(t_new - t_ref) <= self.tol * t_ref:
                        status = "PASS"
                    else:
                        status = "FAIL"
                        failures.append(f"{_where(ephi, key)}: T_new {t_new!r} against T_ref {t_ref!r}, "
                                        f"{rel:+.4%}, beyond {self.tol:g}")
                else:
                    status = "noise"
                    if t_new >= NOISE_REL * largest[ephi]:
                        noted += 1
                        note = (f"T_new {t_new:.6g} is above the noise band ({NOISE_REL:g} x the largest T_ref = "
                                f"{NOISE_REL * largest[ephi]:.6g}) while T_ref {t_ref:.6g} is below it; not a failure")
                counts[status] += 1
                self.rows.append({"phi_deg": key[0], "theta_deg": key[1], "ephi": ephi,
                                  "t_new": _finite_or_none(t_new), "t_ref": _finite_or_none(t_ref),
                                  "rel_diff": _finite_or_none(rel), "status": status, "note": note})
        return _row("transmittance.shared_keys",
                    f"T_new = |Re S| / IngoingPower (new sidecar) against T_ref = OutgoingPower / IngoingPower "
                    f"(reference), no c_k and no flux: |T_new - T_ref| <= {self.tol:g} T_ref where T_ref >= "
                    f"{NOISE_REL:g} of the polarization's largest", failures,
                    f"{counts['PASS']} compared within {self.tol:g}, {counts['noise']} below the band reported as "
                    "noise" + (f" ({noted} with T_new above it, noted on the row)" if noted else ""), {"largest_t_ref": {str(e): _finite_or_none(v) for e, v in largest.items()},
                              "counts": counts})


# --- entry points ----------------------------------------------------------------------------------------------------

def write_report(result: GridResult) -> None:
    """Write ``result.to_dict()`` to ``result.report_target`` (nothing when it is None). An OSError is recorded in
    ``result.report_error``, which fails the run, and never raised."""
    if result.report_target is None:
        return
    try:
        write_json(Path(result.report_target), result.to_dict())
    except OSError as exc:
        result.report_error = f"{result.report_target} could not be written ({exc}); the run fails"
        return
    result.report = result.report_target


def compare(new_job: Path, reference_job: Path, t_rel_tol: float = T_REL_TOL, write: bool = True) -> GridResult:
    """Accept ``new_job`` (schema 2, 35 keys) against ``reference_job`` (schema 1, 15 keys) at the nine shared keys
    and, with ``write``, write ``<new_job>/grid_acceptance.json`` (``write_report``); only into a job whose sidecar
    declares schema 2, so that a swapped call never writes into a 15-key source job.

    Raises ``FileNotFoundError`` when a job directory is missing and ``ValueError`` for a tolerance that is not a
    finite positive number (usage errors); every problem with the jobs' contents becomes a failing row instead.
    """
    new_job, reference_job = Path(new_job), Path(reference_job)
    for label, path in (("new", new_job), ("reference", reference_job)):
        if not path.is_dir():
            raise FileNotFoundError(f"the {label} job directory {path} does not exist")
    if isinstance(t_rel_tol, bool) or not isinstance(t_rel_tol, (int, float)) or not (
            math.isfinite(t_rel_tol) and t_rel_tol > 0.0):
        raise ValueError(f"t_rel_tol must be a finite positive number, got {t_rel_tol!r}")
    acceptance = _Acceptance(new_job, reference_job, float(t_rel_tol))
    checks = acceptance.run()
    version = acceptance.new.sidecar.get("schema_version")
    target = new_job / REPORT_NAME if isinstance(version, str) and version.split(".")[0] == NEW_MAJOR else None
    result = GridResult(str(new_job), str(reference_job), float(t_rel_tol), checks, acceptance.rows,
                        report_target=None if target is None else str(target))
    if write:
        write_report(result)
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m bbsim.grid_acceptance", description=__doc__.splitlines()[0])
    parser.add_argument("--new", type=Path, required=True, help="the schema-2 job on the 35-key grid")
    parser.add_argument("--reference", type=Path, required=True,
                        help="the schema-1 job on the 15-key grid (job_tight1 for the cracks, job_reference1 for the "
                             "round gap)")
    parser.add_argument("--t-rel-tol", type=float, default=T_REL_TOL,
                        help=f"relative tolerance on T at the compared keys (default {T_REL_TOL:g})")
    try:
        args = parser.parse_args(argv)
    except SystemExit as exc:             # argparse has printed its message
        return 2 if exc.code else 0
    try:
        result = compare(args.new, args.reference, args.t_rel_tol, write=False)
    except (FileNotFoundError, ValueError) as exc:
        print(f"grid_acceptance: {exc}", file=sys.stderr)
        return 2
    print(result.format_rows(), flush=True)     # the table reaches the log before the report is written
    write_report(result)
    print(result.format_footer())
    return 0 if result.passed else 1


if __name__ == "__main__":
    sys.exit(main())
