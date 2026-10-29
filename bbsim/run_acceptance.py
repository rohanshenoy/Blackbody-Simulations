"""The array_run self-check: one job-array run judged against its frozen plan (HF-057, HF-054, CDX-038, CDX-040, OR-027).

Usage::

    python -m bbsim.run_acceptance --plan PLAN --plan-sha256 HEX --index N JOB

A frequency the job array solves has no earlier solve of its geometry to pair with, so the pointwise pair check (S42,
HF-054) cannot judge it. The run is accepted instead by recomputing its own stored quantities, with the functions that
made them, and comparing them with each other and with its frozen plan: a self-consistency check, not an independent
validation of the physics at a new frequency. ``JOB`` is task ``N``'s job; ``PLAN`` is the plan file, which must have
the sha256 ``HEX`` the caller recorded, never a digest the job asserts. The plan's ``policy`` decides the rules:

- ``production``: a schema-2 run on the plan's incidence grid, both polarizations exported, the HF-041 column-order
  reference the plan froze embedded;
- ``probe``: the schema-1 single-angle diagnostic (one key, Ephi=0 exported, both polarizations solved; HF-058), never
  a tree input. It records no raw flux, t_flux, W or rho, so ``INAPPLICABLE_TO_PROBE`` are not run, and its
  OutgoingPower keeps schema 1's |Re S| definition.

Every file is read once through one ``acceptance_gate.InputCapture``, which hashes the very bytes each reader parses:
the plan, its config copy and column-order reference copy, the sidecar, the manifest and its copies, the CSVs and the
convergence tables. Every value is recomputed with the producer's own functions (``dataset_v2``, ``sampling``,
``convergence``, ``array_plan.predict_exit_lattice``, ``field_pairing.schema1_grid``), never a second formula. The
task's config is ``array_plan.task_config`` of the captured config bytes. A is the entrance area and
I = Ei^2 / (2 c mu0).

Rows, in order (every row is always present, at severity fail; a row whose input another row found unusable fails
naming it):

1. schema.declared: the job's one sidecar declares schema 2.x (production) or 1.x (probe).
2. identity: dataset ID, frequency and label of the sidecar, the manifest and the plan task agree; the plan is an
   ``array_run`` plan of this policy; the config copy and the column-order reference copy have the plan's sha256 (a
   probe plan records no reference), and a production sidecar embeds that reference; the manifest's effective config
   equals the task config built from the captured config bytes (the manifest records no config digest, and a digest
   of its effective config is not the frozen TOML's); its pre-run project digests are the plan's project copy's (the
   solved project's later bytes are never compared with them); a recorded ``plan`` block names this plan and task;
   the job's CSVs and manifest copies are those its sidecar lists (``acceptance_gate.sidecar_file_problems``).
3. grid: the sidecar's declared keys equal the frozen config's incidence (``sampling.incident_angle_values``), the
   exported polarizations the config's and the policy's, and ``solver.solved_variations`` holds exactly Ephi 0 and 1
   under either policy: the setup's sweep always solves both (HF-058).
4. tables: only the policy's Ephi directories, each holding its two CSVs and the manifest copy; schema 2's columns
   (``schema.WAVEGUIDE_COLUMNS_V2``) or schema 1's (``schema.WAVEGUIDE_COLUMNS``) and ``schema.FAR_FIELD_COLUMNS``; every
   value finite.
5. tables.keys (CDX-040): every row of each CSV at exactly one declared key (``dataset_v2.row_keys``), every declared
   key with rows, the same rows per key (the manifest's retained exit points and far-field points per key); between
   the polarizations the same keys row by row, exit points equal to 1e-12 m and far-field angles to 1e-9 deg.
6. tables.identity (CDX-040): every row's Freq is the plan task's label and its Ephi the directory's.
7. powers.ingoing: IngoingPower finite, positive, constant in each table, equal between the polarizations, and equal
   to ``expected_incoming_power`` (A I from the frozen amplitude and the verified entrance area), as must be the
   manifest's and the sidecar's recorded incoming power and amplitude: agreement among wrong values is not enough.
8. powers.outgoing. Production: OutgoingPower nonnegative and constant within each key; equal to the manifest's
   signed raw flux when lit and exactly 0 when dark (|raw flux| <= ``DARK_REL`` A I, OR-038); equal to
   ``t_flux * A I``; the sidecar's |Re S| the manifest's; ``flux_over_abs_re_s`` equal to ``dataset_v2.sanity_ratio``.
   No T <= 1 gate: this normalization permits T > 1. Probe: finite, nonnegative and constant within the key only.
9. intercept_factor (production; OR-034): every row's InterceptFactor and the sidecar's ``c_key`` and per-key ``c_k``
   equal ``dataset_v2.c_key`` of the declared thetas, ``cell_90`` and its band their ``dataset_v2`` values.
10. exit_field.normalization (production; HF-040, OR-038): a lit basis has a finite positive W, and its stored scaled
    fields, with ``dataset_v2.node_weights`` on the declared lattice (``quadrature_weights`` for a disc), give
    ``dataset_v2.scale_factor`` W' with |W'^2 - 1| <= 1e-6: W is never undone. Without the unscaled export the value of
    a positive W is provenance, not checked. A dark basis has W = 1, OutgoingPower and t_flux exactly 0, its flag set,
    its stored noise fields kept.
11. sampling.normalization (production; HF-042): per lit basis the position normalizer (the stored weighted
    |E_t|^2 sum) and the direction normalizer (``dataset_v2.hemisphere_power_w`` with the manifest's far-field steps)
    finite and positive; ``ff_over_flux`` = that hemisphere power over the stored OutgoingPower, null when dark;
    ``dataset_v2.rho`` equal to the manifest's ``rho_per_key``. The far field stays a pattern.
12. convergence.records (CDX-025): each configured polarization's ``logs/convergence_Ephi<e>.txt`` read and
    ``convergence.qualify``-ed at the frozen config's max_delta_e; its projection (``CONVERGENCE_RECORD_FIELDS``) equals
    the manifest's record field by field; the table's criterion and target are the required ones (``criterion`` is
    not stored, so it is judged here); malformed, negative or non-finite evidence fails; the records and the
    aggregate ``solver.qualified`` consistent (``acceptance_gate.convergence_problems``). Never waivable.
13. convergence (OR-027): every configured polarization qualifies. ``values["waivable"]`` is true exactly when
    convergence.records passed and every unqualified polarization has a finite achieved delta (a missed target or
    AEDT's own ``Converged : No``): only the plan's own assembler may waive it.
14. exit_grid (HF-045, OR-037, OR-048): the manifest's step along l is the task's ``effective_l_resolution_mm`` (the
    requested ``l_resolution_mm`` when given, else the frozen config's), at most ``frequency_grid.max_l_step_mm`` at
    every frequency, as are the plan's and the declared fitted step; the plan's predicted lattice is the HF-045 lattice
    of the frozen rectangle (``array_plan.predict_exit_lattice``, the planner's own); the sidecar's declared lattice
    (``dataset_v2.declared_grid``, or ``field_pairing.schema1_grid`` for schema 1) and the manifest's ``step_mm`` equal
    it (bounds and steps to ``NODE_TOL`` of the step, counts exactly); every key's exit points lie on it, each node once
    (``dataset_v2.node_weights``). At 1010 GHz on crack1 that is 10 mm / 135 = 0.074074 mm and 136 x 51, never the
    requested 0.0742 mm as a step.
15. binding.inputs (CDX-038): every file read once and unchanged at the end, every evidence file
    (``acceptance_gate.evidence_files``) bound.

The report, ``<JOB>/array_run.json``, holds the rows, the versioned tolerances and required rows of the policy, the plan
block {path, sha256, index, policy}, the inapplicable rows and the binding (``reference_job`` null: an array run is
judged against its plan, never a reference). It is written atomically (a temporary file, then ``os.replace``) after
the table is printed; a report that cannot be written fails the run. The table ends in ``RUN ACCEPTANCE: PASS`` or
``RUN ACCEPTANCE: FAIL <n>``. Exit code 0 on PASS, 1 on FAIL, 2 when the job is missing or the plan is refused
(unreadable, another sha256, another schema or policy, no such task).

``bbsim.acceptance_gate`` imports this module only inside its dispatch functions; this module imports the gate.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
import tomllib
from collections.abc import Callable, Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from bbsim import acceptance_gate, array_plan, dataset_v2, field_pairing, frequency_grid, schema
from bbsim.acceptance_rows import RowVerdicts, guard
from bbsim.acceptance_rows import Refused as _Refused
from bbsim.acceptance_rows import at_key as _where
from bbsim.acceptance_rows import get as _get
from bbsim.acceptance_rows import number as _number
from bbsim.acceptance_rows import read_json as _read_json
from bbsim.acceptance_rows import relative as _relative
from bbsim.acceptance_rows import row as _row
from bbsim.acceptance_rows import strict_records as _by_key
from bbsim.compare import Check
from bbsim.config import RunConfig, build_config, config_to_dict
from bbsim.convergence import CRITERION, qualify
from bbsim.dataset_v2 import DatasetV2Error, key_text
from bbsim.extract import PLANE_TOL_M, ExitGrid
from bbsim.manifest import write_json
from bbsim.naming import dataset_dir_name, dataset_record_name, frequency_label
from bbsim.sampling import incident_angle_values, incoming_power_w

POLICY_VERSION = "bbsim-array-run/1"
POLICIES = ("production", "probe")
CHECKER = acceptance_gate.CHECKERS["array_run"]
REPORT_NAME = acceptance_gate.REPORT_NAMES["array_run"]
ROWS = ("schema.declared", "identity", "grid", "tables", "tables.keys", "tables.identity", "powers.ingoing",
        "powers.outgoing", "intercept_factor", "exit_field.normalization", "sampling.normalization",
        "convergence.records", "convergence", "exit_grid", "binding.inputs")
# A schema-1 run records no raw flux, t_flux, W or rho: these rows are never run on a probe, nor fabricated for it.
INAPPLICABLE_TO_PROBE = ("intercept_factor", "exit_field.normalization", "sampling.normalization")
SCHEMA_MAJOR = {"production": "2", "probe": "1"}

NORMALIZATION_REL = 1e-6   # |W'^2 - 1| of the stored fields against the stored power (reexport_check's, HF-040)
STORED_REL = 1e-12         # a stored power or ratio against the value it records
INTERCEPT_REL = 1e-12      # InterceptFactor, c_k and cell_90 against their recomputation (OR-034)
DARK_REL = dataset_v2.DARK_REL
FF_OVER_FLUX_REL = 1e-9    # ff_over_flux against the hemisphere power over the stored flux
RHO_ABS = 1e-9             # rho against the manifest's rho_per_key
AREA_REL = array_plan.AREA_REL  # the entrance face against the frozen rectangle (geometry.check_face_area's)
NODE_TOL = dataset_v2.NODE_TOL
KEY_DEG = dataset_v2.KEY_TOL_DEG
LATTICE_SOURCE = array_plan.LATTICE_SOURCE

_WG = {name: i for i, name in enumerate(schema.WAVEGUIDE_COLUMNS_V2)}   # positions 0-14 are the same in 1.0
_FF = {name: i for i, name in enumerate(schema.FAR_FIELD_COLUMNS)}
EPHI, PHI, THETA = _WG["Ephi"], _WG["IWavePhi"], _WG["IWaveTheta"]
OUTGOING, INGOING, INTERCEPT = _WG["OutgoingPower"], _WG["IngoingPower"], _WG[schema.INTERCEPT_FACTOR]
X, Y, Z = _WG["X"], _WG["Y"], _WG["Z"]
FIELDS = schema.EXIT_FIELD_SLOTS                                         # the six stored exit-field slots
FAR_POINT = slice(_FF["Phi"], _FF["rEtheta_imag"] + 1)                  # Phi, Theta and the four far-field values
KINDS = ("waveguide", "far_field")


def _policy(policy: object) -> str:
    if policy not in POLICIES:
        raise ValueError(f"an array_run policy is one of {', '.join(POLICIES)}, got {policy!r}")
    return policy


def tolerances(policy: str) -> dict[str, object]:
    """The versioned tolerances a report records under ``policy`` and assembly expects
    (``acceptance_gate.expected_tolerances``)."""
    return {"policy_version": POLICY_VERSION, "policy": _policy(policy), "normalization_rel": NORMALIZATION_REL,
            "stored_rel": STORED_REL, "intercept_rel": INTERCEPT_REL, "dark_rel": DARK_REL,
            "ff_over_flux_rel": FF_OVER_FLUX_REL, "rho_abs": RHO_ABS, "area_rel": AREA_REL, "node_tol": NODE_TOL,
            "key_deg": KEY_DEG}


def required_rows(policy: str) -> dict[str, str]:
    """The rows a report must hold under ``policy``, by name, at their severity (``acceptance_gate.expected_rows``)."""
    probe = _policy(policy) == "probe"
    return {name: "fail" for name in ROWS if not (probe and name in INAPPLICABLE_TO_PROBE)}


class PlanRefused(ValueError):
    """The plan cannot be read, is not the one the caller expects or names no such task: nothing is judged."""


# --- small helpers ---------------------------------------------------------------------------------------------------

def _is_sha256(value: object) -> bool:
    return isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


def _real(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _positive(value: object) -> bool:
    return _real(value) and math.isfinite(value) and value > 0


def _within(a: float, b: float, rel: float) -> bool:
    return _relative(float(a), float(b)) <= rel


def _same(a: object, b: object) -> bool:
    """JSON equality, a bool never equal to a number."""
    if isinstance(a, bool) or isinstance(b, bool):
        return a is b
    if isinstance(a, Mapping) and isinstance(b, Mapping):
        return set(a) == set(b) and all(_same(a[k], b[k]) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(_same(x, y) for x, y in zip(a, b))
    return (type(a) is type(b) or (_real(a) and _real(b))) and a == b


def _differences(a: object, b: object, path: str) -> list[str]:
    """Where the JSON value ``a`` differs from ``b``, by dotted path."""
    if isinstance(a, Mapping) and isinstance(b, Mapping):
        out: list[str] = []
        for key in sorted(set(a) | set(b), key=str):
            if key not in a:
                out.append(f"{path}.{key} is absent (the task config gives {b[key]!r})")
            elif key not in b:
                out.append(f"{path}.{key} = {a[key]!r} is not in the task config")
            else:
                out += _differences(a[key], b[key], f"{path}.{key}")
        return out
    return [] if _same(a, b) else [f"{path} = {a!r}, the task config gives {b!r}"]


def _polarizations(record: Mapping, where: str, pols: Sequence[int]) -> Mapping:
    values = record.get("polarizations")
    expected = {str(e) for e in pols}
    if not isinstance(values, Mapping) or set(values) != expected:
        raise _Refused(f"{where}: polarizations must hold exactly {sorted(expected)}, got {values!r}")
    return values


def _optional_number(value: object, where: str) -> float | None:
    return None if value is None else _number(value, where)


def _read_table(stream: Any, columns: Sequence[str], at: str) -> tuple[np.ndarray, dict[str, int]]:
    """A CSV's values by position, column 0 (the Freq label) NaN, and its distinct Freq labels with the first row of
    each; the header must be exactly ``columns``. Parsed from the one captured copy of the file's bytes."""
    header = list(pd.read_csv(stream, nrows=0).columns)
    if header != list(columns):
        raise _Refused(f"{at} has the columns {header}, not {list(columns)}")
    stream.seek(0)
    frame = pd.read_csv(stream, dtype={"Freq": str}, keep_default_na=False)
    if frame.empty:
        raise _Refused(f"{at} has no rows")
    labels = {str(value): int(i) for i, value in frame["Freq"].drop_duplicates().items()}
    numeric = frame.iloc[:, 1:]
    try:
        values = numeric.to_numpy(dtype=float)
    except (TypeError, ValueError):
        values = numeric.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    out = np.empty((len(frame), len(columns)))
    out[:, 0] = np.nan
    out[:, 1:] = values
    return out, labels


# --- the frozen rectangle and A I ---------------------------------------------------------------------------------------

def _canonical_extent(cfg: RunConfig) -> tuple[float, float, float]:
    """The frozen box's extents along the exit axes (p, l, g), in mm: ``array_plan.canonical_extent``, the planner's."""
    try:
        return array_plan.canonical_extent(cfg)
    except array_plan.LatticeError as exc:
        raise _Refused(str(exc)) from None


@dataclass(frozen=True)
class _Lattice:
    """A plan task's predicted exit lattice, in the exit frame (x_e, y_e, z_e) = (p, l, g), mm."""

    lo: tuple[float, float, float]
    hi: tuple[float, float, float]
    step: tuple[float, float, float]
    counts: tuple[int, int, int]
    fitted: tuple[bool, bool, bool]
    extent: tuple[float, float, float]
    area: float


def _plan_lattice(task: Mapping) -> _Lattice:
    where = "the plan task's lattice"
    record = task.get("lattice")
    if not isinstance(record, Mapping):
        raise _Refused(f"{where} must be an object, got {record!r}")

    def vector(name: str) -> tuple[float, float, float]:
        values = record.get(name)
        if not isinstance(values, list) or len(values) != 3:
            raise _Refused(f"{where}.{name} must be a list of 3 numbers, got {values!r}")
        return tuple(_number(v, f"{where}.{name}[{i}]") for i, v in enumerate(values))   # type: ignore[return-value]

    counts = record.get("counts")
    if not (isinstance(counts, list) and len(counts) == 3 and all(type(c) is int and c >= 1 for c in counts)):
        raise _Refused(f"{where}.counts must be 3 positive integers, got {counts!r}")
    fitted = record.get("fitted")
    if not (isinstance(fitted, list) and len(fitted) == 3 and all(isinstance(f, bool) for f in fitted)):
        raise _Refused(f"{where}.fitted must be 3 booleans, got {fitted!r}")
    if record.get("extent_source") != LATTICE_SOURCE:
        raise _Refused(f"{where}.extent_source is {record.get('extent_source')!r}, not {LATTICE_SOURCE!r}")
    return _Lattice(vector("lo_mm"), vector("hi_mm"), vector("step_mm"), tuple(counts), tuple(fitted),
                    vector("extent_mm"), _number(record.get("entrance_area_mm2"), f"{where}.entrance_area_mm2"))


def expected_incoming_power(cfg: RunConfig, task: Mapping, manifest: Mapping, sidecar: Mapping) -> float:
    """A I for the plan's canonical rectangle, independently of every recorded power.

    The nominal entrance area is the frozen box's extent along l times its extent along g (``_canonical_extent``).
    The plan task's predicted full-face lattice extents and recorded extents and area, ``geometry.expected_face_area_mm2``
    when given, the manifest's entrance and exit face areas, and the sidecar's geometry extents and exit cross-section
    must all agree with it to ``AREA_REL`` (the geometry's own face-area tolerance). Only then is the actual entrance
    area the manifest records taken as verified, and A I = ``sampling.incoming_power_w(Ei, area)`` with the frozen
    config's Ei, its mm^2-to-m^2 conversion and constants: the nominal dimensions constrain the geometry but do not
    replace the actual area, and no cosine enters. Raises ``acceptance_rows.Refused`` naming every disagreement."""
    p, l, g = _canonical_extent(cfg)
    area = l * g
    lattice = _plan_lattice(task)
    failures: list[str] = []

    def compare(what: str, value: float, expected: float) -> None:
        if not _within(value, expected, AREA_REL):
            failures.append(f"{what} is {value!r}, the frozen rectangle gives {expected!r} (relative difference "
                            f"{_relative(value, expected):.3g})")

    compare("the plan's predicted lattice extent along l (hi - lo, mm)", lattice.hi[1] - lattice.lo[1], l)
    compare("the plan's predicted lattice extent along g (hi - lo, mm)", lattice.hi[2] - lattice.lo[2], g)
    for value, expected, name in zip(lattice.extent, (p, l, g), "plg"):
        compare(f"the plan's lattice extent_mm along {name}", value, expected)
    compare("the plan's lattice entrance_area_mm2", lattice.area, area)
    if cfg.geometry.expected_face_area_mm2 is not None:
        compare("the frozen config's geometry.expected_face_area_mm2", float(cfg.geometry.expected_face_area_mm2), area)
    faces = {}
    for face in ("entrance_face", "exit_face"):
        faces[face] = _number(_get(manifest, "geometry", face, "area_mm2", where="the manifest"),
                              f"the manifest's geometry.{face}.area_mm2")
        compare(f"the manifest's geometry.{face}.area_mm2", faces[face], area)
    extent = _get(sidecar, "geometry", "extent_mm", where="the sidecar")
    for expected, name in zip((p, l, g), "plg"):
        compare(f"the sidecar's geometry.extent_mm.{name}",
                _number(_get(extent, name, where="the sidecar's geometry.extent_mm"),
                        f"the sidecar's geometry.extent_mm.{name}"), expected)
    section = _get(sidecar, "exit_field", "cross_section", where="the sidecar")
    if not (isinstance(section, Mapping) and section.get("shape") == "rectangle"):
        failures.append(f"the sidecar's exit cross-section is {section!r}, not a rectangle")
    else:
        for name, expected in (("y_e_half_m", l), ("z_e_half_m", g)):
            half = _number(section.get(name), f"the sidecar's exit_field.cross_section.{name}")
            compare(f"twice the sidecar's exit_field.cross_section.{name} (mm)", 2.0 * half * 1e3, expected)
    if failures:
        raise _Refused(f"the entrance face is not the plan's canonical rectangle to {AREA_REL:g} relative, so A I "
                       "cannot be derived from it: " + "; ".join(failures))
    return incoming_power_w(cfg.excitation.ei_v_per_m, faces["entrance_face"])


def _nominal_lattice(cfg: RunConfig) -> dict:
    """The HF-045 lattice of the frozen rectangle at the task config's step along l, bounds centred on the exit face:
    ``array_plan.predict_exit_lattice``, the planner's own prediction, so the two never disagree on it."""
    try:
        return array_plan.predict_exit_lattice(cfg, cfg.exit_field.resolution_mm[1])
    except array_plan.LatticeError as exc:
        raise _Refused(str(exc)) from None


def _node_tol_mm(step: float) -> float:
    """A bound's tolerance on an axis: ``NODE_TOL`` of its step, the exit-plane tolerance on the stepless x_e axis."""
    return NODE_TOL * step if step > 0.0 else PLANE_TOL_M * 1e3


# --- the result ------------------------------------------------------------------------------------------------------

@dataclass
class RunAcceptanceResult(RowVerdicts):
    """The rows in order, the plan block the run was judged against, the binding, where the report goes and whether
    it was written, and the verdict. A report that could not be written fails the run."""

    job: str
    plan: dict[str, Any]               # {path, sha256, index, policy}
    checks: list[Check]
    binding: dict | None = None        # acceptance_gate.binding_block over the bytes the rows read (CDX-038)
    report_target: str | None = None   # <job>/array_run.json
    report: str | None = None          # the JSON written (``write_report``)
    report_error: str | None = None    # why it could not be written

    @property
    def policy(self) -> str:
        return self.plan["policy"]

    @property
    def passed(self) -> bool:
        return super().passed and self.report_error is None

    @property
    def verdict(self) -> str:
        if self.passed:
            return "RUN ACCEPTANCE: PASS"
        return f"RUN ACCEPTANCE: FAIL {len(self.failures) + (1 if self.report_error else 0)}"

    def to_dict(self) -> dict[str, Any]:
        return {"passed": self.passed, "checker": CHECKER, "failed": len(self.failures), "verdict": self.verdict,
                "job": self.job, "checks": [asdict(c) for c in self.checks], "tolerances": tolerances(self.policy),
                "required_rows": required_rows(self.policy), "plan": dict(self.plan), "binding": self.binding,
                "inapplicable": list(INAPPLICABLE_TO_PROBE) if self.policy == "probe" else []}

    def format_table(self) -> str:
        return self.format_rows() + "\n" + self.format_footer()

    def format_footer(self) -> str:
        """The report line and the verdict: what ``main`` prints after writing the report."""
        if self.report:
            report = f"report: {self.report}"
        elif self.report_error:
            report = f"report: {self.report_error}"
        else:
            report = "report: not written"
        return "\n".join(["", report, "", self.verdict])

    def format_rows(self) -> str:
        """The header and the rows: what ``main`` prints before writing the report."""
        lines = [f"job:    {self.job}", f"plan:   {self.plan['path']} (sha256 {self.plan['sha256']}), task "
                 f"{self.plan['index']}, policy {self.policy}",
                 "tolerances: " + " ".join(f"{k}={v}" for k, v in tolerances(self.policy).items()), ""]
        for c in self.checks:
            lines.append(f"{'PASS' if c.passed else 'FAIL':4} {c.name:25} {c.detail}")
        if self.policy == "probe":
            lines.append(f"     not run under the probe policy: {', '.join(INAPPLICABLE_TO_PROBE)}")
        return "\n".join(lines)


# --- the acceptance --------------------------------------------------------------------------------------------------

class _Acceptance:
    def __init__(self, job: Path, plan: dict, task: Mapping, context: acceptance_gate.PlanContext,
                 inputs: acceptance_gate.InputCapture):
        self.job, self.plan, self.task, self.context, self.inputs = job, plan, task, context, inputs
        self.policy, self.index = context.policy, context.index
        self.probe = self.policy == "probe"
        self.pols = acceptance_gate.ARRAY_POLICIES[self.policy]
        self.wg_columns = schema.WAVEGUIDE_COLUMNS if self.probe else schema.WAVEGUIDE_COLUMNS_V2
        self.dataset_id: str = plan["dataset_id"]
        self.frequency: float = task["frequency_ghz"]
        self.label = frequency_label(self.frequency)
        self.verdicts: dict[str, bool] = {}
        self._cache: dict[Any, tuple[Any, BaseException | None]] = {}

    def _cached(self, name: Any, compute: Callable[[], Any]) -> Any:
        if name not in self._cache:
            try:
                self._cache[name] = (compute(), None)
            except Exception as exc:  # noqa: BLE001 - re-raised to every row that needs it
                self._cache[name] = (None, exc)
        value, exc = self._cache[name]
        if exc is not None:
            raise exc
        return value

    def run(self) -> list[Check]:
        steps = {"schema.declared": self._schema_declared, "identity": self._identity, "grid": self._grid_row,
                 "tables": self._tables_row, "tables.keys": self._table_keys, "tables.identity": self._table_identity,
                 "powers.ingoing": self._ingoing, "powers.outgoing": self._outgoing,
                 "intercept_factor": self._intercept, "exit_field.normalization": self._normalization,
                 "sampling.normalization": self._sampling, "convergence.records": self._convergence_records,
                 "convergence": self._convergence, "exit_grid": self._exit_grid,
                 "binding.inputs": self._binding_inputs}
        checks: list[Check] = []
        for name in required_rows(self.policy):
            (check,) = guard([name], lambda step=steps[name]: [step()])   # a refusal is a failing row (CDX-030)
            if name == "convergence":
                check.values.setdefault("waivable", False)
            self.verdicts[name] = check.passed
            checks.append(check)
        return checks

    # -- what the rows read, each once through the capture ----------------------------------------------------------
    def _copy(self, name: str) -> Path:
        entry = self.plan.get(name)
        copy = entry.get("copy") if isinstance(entry, Mapping) else None
        if not isinstance(copy, str) or not os.path.isabs(copy):
            raise _Refused(f"the plan's {name}.copy must be an absolute path, got {copy!r}")
        return Path(copy)

    def _sidecar(self) -> dict:
        def read() -> dict:
            sidecars = sorted(self.job.glob("*.dataset.json"))
            if len(sidecars) != 1:
                names = f": {[p.name for p in sidecars]}" if sidecars else ""
                raise _Refused(f"the job holds {len(sidecars)} dataset sidecars (*.dataset.json), not 1{names}")
            return _read_json(sidecars[0], "the sidecar", read_text=self.inputs.read_text)
        return self._cached("sidecar", read)

    def _manifest(self) -> dict:
        return self._cached("manifest", lambda: _read_json(self.job / "manifest.json", "the job manifest",
                                                           read_text=self.inputs.read_text))

    def _raw_config(self) -> dict:
        def read() -> dict:
            copy = self._copy("config")
            try:
                return self.inputs.parse(copy, tomllib.load)
            except (OSError, ValueError) as exc:
                raise _Refused(f"the plan's config copy {copy} is unreadable ({exc})") from exc
        return self._cached("raw_config", read)

    def _cfg(self) -> RunConfig:
        def build() -> RunConfig:
            raw = self._raw_config()
            try:
                return array_plan.task_config(self.plan, self.index, raw_config=raw)
            except ValueError as exc:
                raise _Refused(f"task {self.index}'s config cannot be built from the plan's config copy ({exc})") from exc
        return self._cached("cfg", build)

    def _reference(self) -> dict:
        def read() -> dict:
            copy = self._copy("column_order_reference")
            try:
                reference = self.inputs.parse(copy, json.load)
            except (OSError, ValueError) as exc:
                raise _Refused(f"the plan's column-order reference copy {copy} is unreadable ({exc})") from exc
            if not isinstance(reference, dict):
                raise _Refused(f"the plan's column-order reference copy {copy} is not a JSON object")
            return reference
        return self._cached("reference", read)

    def _keys(self) -> list[tuple[float, float]]:
        def read() -> list[tuple[float, float]]:
            try:
                return field_pairing.declared_keys(self._sidecar(), "the sidecar")
            except field_pairing.PairingRefused as exc:
                raise _Refused(str(exc)) from exc
        return self._cached("keys", read)

    def _thetas(self) -> list[float]:
        return sorted({theta for _, theta in self._keys()})

    def _tables(self) -> tuple[dict[int, dict[str, Any]], list[str]]:
        def load() -> tuple[dict[int, dict[str, Any]], list[str]]:
            problems: list[str] = []
            names = {e: dataset_dir_name(self.dataset_id, self.frequency, e) for e in self.pols}
            unexpected = sorted(p.name for p in self.job.glob("*_Ephi=*") if p.name not in names.values())
            if unexpected:
                problems.append(f"the job holds {unexpected} beside the {self.policy} policy's {sorted(names.values())}"
                                f" (Ephi {list(self.pols)} exported)")
            tables: dict[int, dict[str, Any]] = {}
            for e in self.pols:
                d = self.job / names[e]
                if not d.is_dir():
                    problems.append(f"the Ephi={e} directory {d.name} is missing")
                    continue
                held = sorted(p.name for p in d.iterdir())
                if held != sorted(acceptance_gate.EPHI_FILES):
                    problems.append(f"{d.name} holds {held}, not {sorted(acceptance_gate.EPHI_FILES)}")
                entry: dict[str, Any] = {}
                for kind, columns in zip(KINDS, (self.wg_columns, schema.FAR_FIELD_COLUMNS)):
                    path, at = d / f"{kind}.csv", f"the Ephi={e} {kind}.csv"
                    try:
                        values, labels = self.inputs.parse(path, lambda f, c=columns, at=at: _read_table(f, c, at))
                    except _Refused as exc:
                        problems.append(str(exc))
                        continue
                    except (OSError, ValueError) as exc:
                        problems.append(f"{at} ({path}) is unreadable: {exc}")
                        continue
                    bad = ~np.isfinite(values[:, 1:])
                    if bad.any():
                        i = int(np.flatnonzero(bad.any(axis=1))[0])
                        problems.append(f"{at}: {int(bad.any(axis=1).sum())} rows with a non-finite value, the first "
                                        f"at row {i} in {columns[1 + int(np.flatnonzero(bad[i])[0])]}")
                        continue
                    entry[kind], entry[f"{kind}_freq"] = values, labels
                # Only readable, finite tables with the expected columns are usable by the other rows; a layout problem
                # (another directory, a stray file) fails this row alone.
                if all(kind in entry for kind in KINDS):
                    tables[e] = entry
            return tables, problems
        return self._cached("tables", load)

    def _usable(self) -> dict[int, dict[str, Any]]:
        tables, _ = self._tables()
        missing = [e for e in self.pols if e not in tables]
        if missing:
            raise _Refused(f"the Ephi={missing[0]} tables are not usable (tables)")
        return tables

    def _index(self) -> tuple[dict[tuple[str, int], list[np.ndarray]], list[str]]:
        def index() -> tuple[dict[tuple[str, int], list[np.ndarray]], list[str]]:
            tables, keys = self._usable(), self._keys()
            rows: dict[tuple[str, int], list[np.ndarray]] = {}
            problems: list[str] = []
            for e in self.pols:
                for kind in KINDS:
                    values, at = tables[e][kind], f"the Ephi={e} {kind}.csv"
                    frame = pd.DataFrame({"IWavePhi": values[:, PHI], "IWaveTheta": values[:, THETA]})
                    try:
                        key_of_row = dataset_v2.row_keys(frame, keys, at)
                    except DatasetV2Error as exc:
                        problems.append(str(exc))
                        continue
                    counts = np.bincount(key_of_row, minlength=len(keys))
                    uneven = np.flatnonzero(counts != counts[0])
                    if uneven.size:
                        k = int(uneven[0])
                        problems.append(f"{at} has {int(counts[k])} rows at {key_text(keys[k])} but {int(counts[0])} "
                                        f"at {key_text(keys[0])}: every key holds the same points")
                    rows[(kind, e)] = [np.flatnonzero(key_of_row == k) for k in range(len(keys))]
            return rows, problems
        return self._cached("index", index)

    def _rows(self, kind: str, e: int) -> list[np.ndarray]:
        rows, _ = self._index()
        if (kind, e) not in rows:
            raise _Refused(f"the Ephi={e} {kind}.csv rows do not lie at the declared keys (tables.keys)")
        return rows[(kind, e)]

    def _grid(self) -> ExitGrid:
        def parse() -> ExitGrid:
            where = "the sidecar's exit_field.grid"
            record = _get(self._sidecar(), "exit_field", "grid", where="the sidecar")
            try:
                if self.probe:
                    return field_pairing.schema1_grid(record, where)
                return dataset_v2.declared_grid(record, where)
            except (field_pairing.PairingRefused, DatasetV2Error) as exc:
                raise _Refused(str(exc)) from exc
        return self._cached("grid", parse)

    def _shape(self) -> object:
        section = _get(self._sidecar(), "exit_field", "cross_section", where="the sidecar")
        return section.get("shape") if isinstance(section, Mapping) else None

    def _weights(self, e: int, k: int) -> np.ndarray:
        """The quadrature weights of the key's exit points from their own lattice nodes (``dataset_v2.node_weights``),
        refused unless they are the whole declared lattice, each node once; 1 on every point of a disc."""
        def compute() -> np.ndarray:
            rows, values, grid = self._rows("waveguide", e)[k], self._usable()[e]["waveguide"], self._grid()
            shape = self._shape()
            try:
                if shape == "rectangle":
                    return dataset_v2.node_weights(values[rows, Y], values[rows, Z], grid)
                if shape == "disc":
                    return dataset_v2.quadrature_weights("disc", grid.counts[1], grid.counts[2], rows.size)
            except DatasetV2Error as exc:
                raise _Refused(f"{_where(e, self._keys()[k])}: {exc}") from exc
            raise _Refused(f"the sidecar's exit cross-section shape {shape!r} has no quadrature weights")
        return self._cached(("weights", e, k), compute)

    def _ai(self) -> float:
        return self._cached("ai", lambda: expected_incoming_power(self._cfg(), self.task, self._manifest(),
                                                                  self._sidecar()))

    def _raw(self) -> dict[tuple[int, int], tuple[float, float]]:
        def parse() -> dict[tuple[int, int], tuple[float, float]]:
            keys, where = self._keys(), "the manifest's exit_field.raw_powers_per_key"
            found = _by_key(_get(self._manifest(), "exit_field", "raw_powers_per_key", where="the manifest"), keys,
                            where)
            out = {}
            for k, record in found.items():
                values = _polarizations(record, f"{where} at {key_text(keys[k])}", self.pols)
                for e in self.pols:
                    at = f"the raw powers of {_where(e, keys[k])}"
                    flux = _number(_get(values[str(e)], "exit_flux_re_s_dot_n_w", where=at),
                                   f"{at}: exit_flux_re_s_dot_n_w")
                    absolute = _number(_get(values[str(e)], "outgoing_power_abs_re_s_w", where=at),
                                       f"{at}: outgoing_power_abs_re_s_w")
                    if absolute < 0.0:
                        raise _Refused(f"{at}: the |Re S| integral {absolute!r} W is negative")
                    out[(e, k)] = (flux, absolute)
            return out
        return self._cached("raw", parse)

    def _dark(self) -> dict[tuple[int, int], bool]:
        """Per basis, whether its raw flux is a numerical zero: |raw flux| <= DARK_REL A I (OR-038)."""
        return self._cached("dark", lambda: {ek: abs(flux) <= DARK_REL * self._ai()
                                             for ek, (flux, _) in self._raw().items()})

    def _per_key(self) -> dict[int, dict[str, Any]]:
        def parse() -> dict[int, dict[str, Any]]:
            keys, where = self._keys(), "the sidecar's transmittance.per_key"
            out = {}
            for k, record in _by_key(_get(self._sidecar(), "transmittance", "per_key", where="the sidecar"), keys,
                                     where).items():
                at = f"{where} at {key_text(keys[k])}"
                values = _polarizations(record, at, self.pols)
                out[k] = {"c_k": _number(record.get("c_k"), f"{at}: c_k"),
                          "t_flux": {e: _number(_get(values[str(e)], "t_flux", where=at), f"{at}: t_flux[{e}]")
                                     for e in self.pols},
                          "abs": {e: _number(_get(values[str(e)], "outgoing_power_abs_re_s_w", where=at),
                                             f"{at}: outgoing_power_abs_re_s_w[{e}]") for e in self.pols}}
            return out
        return self._cached("per_key", parse)

    def _per_basis(self, records: object, where: str, field: str | None = None) -> dict[tuple[int, int], float | None]:
        """A per-key list whose records hold one number or null per polarization, under ``polarizations`` (or under
        ``field`` within each polarization), by (Ephi, key index)."""
        keys = self._keys()
        out = {}
        for k, record in _by_key(records, keys, where).items():
            values = _polarizations(record, f"{where} at {key_text(keys[k])}", self.pols)
            for e in self.pols:
                value = values[str(e)] if field is None else _get(values[str(e)], field, where=where)
                out[(e, k)] = _optional_number(value, f"{where} at {_where(e, keys[k])}")
        return out

    def _scale(self) -> dict[int, Mapping]:
        def parse() -> dict[int, Mapping]:
            where = "the sidecar's exit_field.scale.per_key"
            return _by_key(_get(self._sidecar(), "exit_field", "scale", "per_key", where="the sidecar"), self._keys(),
                           where)
        return self._cached("scale", parse)

    def _qualifications(self) -> dict[int, dict | _Refused]:
        def compute() -> dict[int, dict | _Refused]:
            target = self._cfg().solver.max_delta_e
            out: dict[int, dict | _Refused] = {}
            for e in self.pols:
                path = self.job / "logs" / f"convergence_Ephi{e}.txt"
                try:
                    text = self.inputs.read_text(path)
                except (OSError, ValueError) as exc:
                    out[e] = _Refused(f"the Ephi={e} convergence table {path} is unreadable ({exc}): it is no evidence "
                                      "of convergence")
                    continue
                out[e] = qualify(text, target)
            return out
        return self._cached("qualifications", compute)

    # -- 1. schema.declared ---------------------------------------------------------------------------------------
    def _schema_declared(self) -> Check:
        major, problems, version = SCHEMA_MAJOR[self.policy], [], None
        try:
            version = self._sidecar().get("schema_version")
        except _Refused as exc:
            problems.append(str(exc))
        else:
            if not (isinstance(version, str) and version.split(".")[0] == major):
                problems.append(f"the sidecar declares schema_version {version!r}; the plan's {self.policy} policy "
                                f"requires schema {major}.x")
        return _row("schema.declared", f"the job's one sidecar declares schema {major}.x, as the plan's {self.policy} "
                    "policy requires", problems, f"schema_version {version!r}",
                    {"declared": version, "policy": self.policy})

    # -- 2. identity ------------------------------------------------------------------------------------------------
    def _identity(self) -> Check:
        problems: list[str] = []
        plan, task = self.plan, self.task
        if plan.get("acceptance_kind") != "array_run":
            problems.append(f"the plan's acceptance_kind is {plan.get('acceptance_kind')!r}, not 'array_run'")
        if plan.get("probe") is not self.probe:
            problems.append(f"the plan's probe flag {plan.get('probe')!r} contradicts its {self.policy} policy")
        if task.get("label") != self.label:
            problems.append(f"the plan task's label {task.get('label')!r} is not {self.label!r}, the label of its "
                            f"frequency {self.frequency!r} GHz")
        # The frozen inputs: the captured bytes must be the ones the plan hashed.
        for name, read in (("config", self._raw_config), ("column_order_reference", self._reference)):
            if name == "column_order_reference":
                entry = plan.get(name)
                if self.probe:
                    if entry is not None:
                        problems.append(f"a probe plan records no column-order reference (schema 1 neither needs nor "
                                        f"applies one), this one records {entry!r}")
                    continue
                if entry is None:
                    problems.append("a production plan records the HF-041 column-order reference its schema-2 runs "
                                    "embed, and this one records none")
                    continue
            try:
                read()
            except _Refused as exc:
                problems.append(str(exc))
                continue
            digest, recorded = self.inputs.recorded(self._copy(name)), plan[name].get("sha256")
            if digest != recorded:
                problems.append(f"the plan's {name} copy {self._copy(name)} has sha256 {digest}, the plan records "
                                f"{recorded!r}")
        cfg = None
        try:
            cfg = self._cfg()
        except _Refused as exc:
            problems.append(str(exc))
        else:
            if cfg.acceptance.kind != "array_run":
                problems.append(f"the task config's acceptance kind is {cfg.acceptance.kind!r}, not 'array_run'")
            if cfg.project.dataset_id != self.dataset_id:
                problems.append(f"the task config's dataset_id {cfg.project.dataset_id!r} is not the plan's "
                                f"{self.dataset_id!r}")
        try:
            problems += self._manifest_identity(cfg)
        except _Refused as exc:
            problems.append(str(exc))
        try:
            problems += self._sidecar_identity()
        except _Refused as exc:
            problems.append(str(exc))
        what = "config, project and column-order reference" if not self.probe else "config and project"
        return _row("identity", "dataset ID, frequency and label agree across the sidecar, the manifest and the plan "
                    f"task; the {what} copies are the plan's bytes; the manifest's effective config is the task "
                    "config of the frozen copy; the job's files are those its sidecar lists", problems,
                    f"{self.dataset_id} at {self.label}, task {self.index} of the plan")

    def _manifest_identity(self, cfg: RunConfig | None) -> list[str]:
        manifest, problems = self._manifest(), []
        if cfg is not None:
            expected = json.loads(json.dumps(config_to_dict(cfg)))
            diffs = _differences(manifest.get("config"), expected, "config")
            if diffs:
                problems.append(f"the manifest's effective config is not the task config built from the plan's copy: "
                                f"{'; '.join(diffs[:5])}" + (f"; and {len(diffs) - 5} more" if len(diffs) > 5 else ""))
        project, plan_project = manifest.get("project"), self.plan.get("project")
        if not (isinstance(project, Mapping) and isinstance(plan_project, Mapping)):
            problems.append(f"the manifest's project {project!r} or the plan's {plan_project!r} is not an object")
        else:
            for field in ("job_copy_sha256_before_run", "source_sha256"):
                if project.get(field) != plan_project.get("sha256"):
                    problems.append(f"the manifest's project.{field} {project.get(field)!r} is not the plan's project "
                                    f"digest {plan_project.get('sha256')!r}")
            source, copy = project.get("source_path"), plan_project.get("copy")
            if not (isinstance(source, str) and isinstance(copy, str) and Path(source).resolve() == Path(copy).resolve()):
                problems.append(f"the manifest's project.source_path {source!r} is not the plan's project copy {copy!r}")
        label = _get(manifest, "solver", "frequency_label", where="the manifest", default=None)
        if label != self.label:
            problems.append(f"the manifest's solver.frequency_label {label!r} is not {self.label!r}")
        recorded = manifest.get("plan")
        if recorded is not None:            # the runner's plan mode records which plan and task it ran
            context = self.context
            if not (isinstance(recorded, Mapping) and isinstance(recorded.get("path"), str)
                    and Path(recorded["path"]).resolve() == Path(context.plan_path).resolve()
                    and recorded.get("sha256") == context.plan_sha256 and type(recorded.get("index")) is int
                    and recorded["index"] == context.index):
                problems.append(f"the manifest's plan block {recorded!r} does not name this plan "
                                f"({context.plan_path}, {context.plan_sha256}) and task {context.index}")
        return problems

    def _sidecar_identity(self) -> list[str]:
        sidecar, problems = self._sidecar(), []
        if sidecar.get("dataset_id") != self.dataset_id:
            problems.append(f"the sidecar's dataset_id {sidecar.get('dataset_id')!r} is not the plan's "
                            f"{self.dataset_id!r}")
        frequency = sidecar.get("frequency_ghz")
        if not (_real(frequency) and frequency == self.frequency):
            problems.append(f"the sidecar's frequency_ghz {frequency!r} is not the plan task's {self.frequency!r}")
        if sidecar.get("frequency_label") != self.label:
            problems.append(f"the sidecar's frequency_label {sidecar.get('frequency_label')!r} is not {self.label!r}")
        name = dataset_record_name(self.dataset_id, self.frequency)
        if not (self.job / name).is_file():
            problems.append(f"the job's sidecar is not named {name}")
        if not self.probe and self.plan.get("column_order_reference") is not None:
            block = _get(sidecar, "exit_field", "column_order_check", where="the sidecar")
            try:
                reference = self._reference()
            except _Refused:
                reference = None                # reported with the copy
            if not (isinstance(block, Mapping) and block.get("method") == "shared_exporter"
                    and reference is not None and _same(block.get("reference"), reference)):
                problems.append("the sidecar's exit_field.column_order_check does not embed the plan's column-order "
                                "reference by the shared_exporter method (HF-041)")
        return problems + acceptance_gate.sidecar_file_problems(self.job, sidecar, self.inputs)

    # -- 3. grid ------------------------------------------------------------------------------------------------------
    def _grid_row(self) -> Check:
        cfg, manifest, declared = self._cfg(), self._manifest(), self._keys()
        exc = cfg.excitation
        expected = [(float(p), float(t))
                    for p in incident_angle_values(exc.phi_lower_deg, exc.phi_upper_deg, exc.phi_step_deg)
                    for t in incident_angle_values(exc.theta_lower_deg, exc.theta_upper_deg, exc.theta_step_deg)]
        problems = []
        if len(declared) != len(expected) or any(abs(a[0] - b[0]) > KEY_DEG or abs(a[1] - b[1]) > KEY_DEG
                                                 for a, b in zip(declared, expected)):
            problems.append(f"the sidecar declares {len(declared)} keys, the frozen config's excitation gives "
                            f"{len(expected)} (phi outer, theta inner, to {KEY_DEG:g} deg)")
        pols = list(cfg.output.polarizations)
        if tuple(pols) != self.pols:
            problems.append(f"the task config exports Ephi {pols}; the plan's {self.policy} policy exports "
                            f"{list(self.pols)}")
        recorded = _get(manifest, "excitation", "polarizations", where="the manifest", default=None)
        if not _same(recorded, pols):
            problems.append(f"the manifest's excitation.polarizations {recorded!r} is not the config's {pols}")
        outputs = manifest.get("outputs")
        if not (isinstance(outputs, Mapping) and set(outputs) == {str(e) for e in pols}):
            problems.append(f"the manifest's outputs {sorted(outputs) if isinstance(outputs, Mapping) else outputs!r} "
                            f"are not Ephi {pols}")
        solved = _get(manifest, "solver", "solved_variations", where="the manifest", default=None)
        found = ([re.findall(r"Ephi='([^']*)'", v) for v in solved]
                 if isinstance(solved, list) and all(isinstance(v, str) for v in solved) else None)
        if found is None or any(len(f) != 1 for f in found) or sorted(f[0] for f in found) != ["0", "1"]:
            problems.append(f"the manifest's solver.solved_variations {solved!r} does not hold exactly Ephi 0 and 1, "
                            "one variation each: the setup's sweep always solves both (HF-058)")
        return _row("grid", "the sidecar's declared keys are the frozen config's incidence, the exported polarizations "
                    "the config's and the policy's, and the solved variations exactly Ephi 0 and 1", problems,
                    f"{len(declared)} keys, Ephi {pols} exported, Ephi 0 and 1 solved")

    # -- 4. tables ------------------------------------------------------------------------------------------------------
    def _tables_row(self) -> Check:
        tables, problems = self._tables()
        summary = "; ".join(f"Ephi={e}: {len(t['waveguide'])} waveguide and {len(t['far_field'])} far-field rows"
                            for e, t in sorted(tables.items()))
        columns = "schema.WAVEGUIDE_COLUMNS" if self.probe else "schema.WAVEGUIDE_COLUMNS_V2"
        return _row("tables", f"only the {self.policy} policy's Ephi {list(self.pols)} directories, each with its two "
                    f"CSVs and the manifest copy; columns {columns} and schema.FAR_FIELD_COLUMNS; every value finite",
                    problems, summary)

    # -- 5. tables.keys -----------------------------------------------------------------------------------------------
    def _table_keys(self) -> Check:
        rows, problems = self._index()
        problems = list(problems)
        manifest, tables, keys = self._manifest(), self._usable(), self._keys()
        expected = {"waveguide": _get(manifest, "exit_field", "retained_points_per_key", where="the manifest"),
                    "far_field": _get(manifest, "far_field", "points_per_angle", where="the manifest")}
        for (kind, e), per_key in sorted(rows.items()):
            counts = sorted({int(r.size) for r in per_key})
            if counts != [expected[kind]]:
                problems.append(f"the Ephi={e} {kind}.csv holds {counts} rows per key, the manifest records "
                                f"{expected[kind]!r}")
        if len(self.pols) == 2:
            for kind, coords, tol, unit in (("waveguide", [X, Y, Z], field_pairing.COORD_M, "m"),
                                            ("far_field", [_FF["Phi"], _FF["Theta"]], field_pairing.ANGLE_DEG, "deg")):
                a, b = tables[0][kind], tables[1][kind]
                if a.shape[0] != b.shape[0]:
                    problems.append(f"the {kind}.csv of Ephi=0 holds {a.shape[0]} rows, of Ephi=1 {b.shape[0]}")
                    continue
                if not np.array_equal(a[:, [PHI, THETA]], b[:, [PHI, THETA]]):
                    i = int(np.flatnonzero(np.any(a[:, [PHI, THETA]] != b[:, [PHI, THETA]], axis=1))[0])
                    problems.append(f"the {kind}.csv of Ephi=0 and Ephi=1 differ in key at row {i}: BBRsim pairs the "
                                    "polarizations' rows by index")
                    continue
                bad = ~(np.abs(a[:, coords] - b[:, coords]) <= tol)
                if bad.any():
                    i = int(np.flatnonzero(bad.any(axis=1))[0])
                    problems.append(f"the {kind}.csv of Ephi=0 and Ephi=1 differ at row {i} "
                                    f"({key_text((float(a[i, PHI]), float(a[i, THETA])))}) by more than {tol:g} {unit}")
        return _row("tables.keys", "every row of each CSV at exactly one declared key, every declared key with the "
                    "manifest's points per key; between the polarizations the same keys, exit points (1e-12 m) and "
                    "far-field angles (1e-9 deg) row by row", problems,
                    f"{len(keys)} keys x {expected['waveguide']} exit and {expected['far_field']} far-field points, "
                    f"Ephi {list(self.pols)}")

    # -- 6. tables.identity -----------------------------------------------------------------------------------------
    def _table_identity(self) -> Check:
        problems = []
        for e, entry in sorted(self._usable().items()):
            for kind in KINDS:
                at = f"the Ephi={e} {kind}.csv"
                wrong = {label: row for label, row in entry[f"{kind}_freq"].items() if label != self.label}
                if wrong:
                    label, row = min(wrong.items(), key=lambda item: item[1])
                    problems.append(f"{at}: Freq {label!r} (first at row {row}), not the plan task's label "
                                    f"{self.label!r}")
                bad = np.flatnonzero(~(entry[kind][:, EPHI] == e))
                if bad.size:
                    i = int(bad[0])
                    problems.append(f"{at}: Ephi {float(entry[kind][i, EPHI])!r} on {bad.size} rows (the first at row "
                                    f"{i}), not the directory's {e}")
        return _row("tables.identity", "every row's Freq the plan task's label and its Ephi the directory's (the "
                    "producer's record; BBRsim reads the directory)", problems,
                    f"Freq {self.label} and Ephi {list(self.pols)} on every row")

    # -- 7. powers.ingoing --------------------------------------------------------------------------------------------
    def _ingoing(self) -> Check:
        problems: list[str] = []
        try:
            ai = self._ai()
        except _Refused as exc:
            problems.append(str(exc))
            ai = None
        constants: dict[int, float] = {}
        for e, entry in sorted(self._usable().items()):
            values, at = entry["waveguide"][:, INGOING], f"the Ephi={e} waveguide.csv"
            bad = ~(values > 0.0)
            if bad.any():
                problems.append(f"{at}: IngoingPower {float(values[np.flatnonzero(bad)[0]])!r} W is not positive")
                continue
            first = float(values[0])
            spread = np.abs(values - first) / first
            if not np.all(spread <= STORED_REL):
                i = int(np.argmax(spread))
                problems.append(f"{at}: IngoingPower {float(values[i])!r} W at row {i} differs from its first row's "
                                f"{first!r} W by {float(spread[i]):.3g} relative: not constant")
                continue
            constants[e] = first
            if ai is not None and not _within(first, ai, STORED_REL):
                problems.append(f"{at}: IngoingPower {first!r} W is not A I = {ai!r} W")
        if len(constants) == 2 and not _within(constants[1], constants[0], STORED_REL):
            problems.append(f"IngoingPower differs between the polarizations: {constants[0]!r} and {constants[1]!r} W")
        cfg = self._cfg()
        for label, read in (("manifest", self._manifest), ("sidecar", self._sidecar)):
            try:
                block = _get(read(), "excitation", where=f"the {label}")
                recorded = _number(_get(block, "incoming_power_w", where=f"the {label}'s excitation"),
                                   f"the {label}'s excitation.incoming_power_w")
                amplitude = _number(_get(block, "ei_v_per_m", where=f"the {label}'s excitation"),
                                    f"the {label}'s excitation.ei_v_per_m")
            except _Refused as exc:
                problems.append(str(exc))
                continue
            if ai is not None and not _within(recorded, ai, STORED_REL):
                problems.append(f"the {label}'s excitation.incoming_power_w {recorded!r} W is not A I = {ai!r} W")
            if amplitude != cfg.excitation.ei_v_per_m:
                problems.append(f"the {label}'s excitation.ei_v_per_m {amplitude!r} is not the frozen config's "
                                f"{cfg.excitation.ei_v_per_m!r}")
        return _row("powers.ingoing", f"IngoingPower finite, positive, constant in each table and equal between the "
                    f"polarizations; it and the recorded incoming powers equal A I = sampling.incoming_power_w(Ei, A) "
                    f"from the frozen amplitude and the verified entrance area, to {STORED_REL:g}", problems,
                    f"A I = {ai!r} W" if ai is not None else "")

    # -- 8. powers.outgoing -------------------------------------------------------------------------------------------
    def _outgoing(self) -> Check:
        if self.probe:
            return self._outgoing_probe()
        tables, keys, ai = self._usable(), self._keys(), self._ai()
        raw, per_key = self._raw(), self._per_key()
        ratios = self._per_basis(_get(self._manifest(), "transmittance_checks", "flux_over_abs_re_s",
                                      where="the manifest"), "the manifest's transmittance_checks.flux_over_abs_re_s")
        dark = self._dark()
        problems: list[str] = []
        for e in self.pols:
            rows = self._rows("waveguide", e)
            for k, key in enumerate(keys):
                where = _where(e, key)
                stored = tables[e]["waveguide"][rows[k], OUTGOING]
                first = float(stored[0])
                if (stored < 0.0).any():
                    problems.append(f"{where}: OutgoingPower {float(stored.min())!r} W is negative")
                if np.any(np.abs(stored - first) > STORED_REL * abs(first)):
                    problems.append(f"{where}: OutgoingPower is not constant within the key (from {first!r} W to "
                                    f"{float(stored[np.argmax(np.abs(stored - first))])!r} W)")
                flux, absolute = raw[(e, k)]
                t = per_key[k]["t_flux"][e]
                if dark[(e, k)]:
                    if np.any(stored != 0.0):
                        problems.append(f"{where}: the raw flux {flux!r} W is a numerical zero, so OutgoingPower is 0, "
                                        f"not {float(stored[np.flatnonzero(stored)[0]])!r} W (OR-038)")
                    if t != 0.0:
                        problems.append(f"{where}: a dark basis's t_flux is 0, not {t!r}")
                else:
                    if not _within(first, flux, STORED_REL):
                        problems.append(f"{where}: OutgoingPower {first!r} W is not the manifest's raw flux {flux!r} W")
                    if not _within(t * ai, first, STORED_REL):
                        problems.append(f"{where}: t_flux * A I = {t * ai!r} W is not OutgoingPower {first!r} W")
                if not _within(per_key[k]["abs"][e], absolute, STORED_REL):
                    problems.append(f"{where}: the sidecar's |Re S| {per_key[k]['abs'][e]!r} W is not the manifest's "
                                    f"{absolute!r} W")
                try:
                    expected = dataset_v2.sanity_ratio(flux, absolute, ai)
                except DatasetV2Error as exc:
                    problems.append(f"{where}: {exc}")
                    continue
                recorded = ratios[(e, k)]
                if (expected is None) != (recorded is None) or (
                        expected is not None and not _within(recorded, expected, STORED_REL)):
                    problems.append(f"{where}: flux_over_abs_re_s {recorded!r} is not flux / |Re S| = {expected!r}")
        n_dark = sum(dark.values())
        return _row("powers.outgoing", f"OutgoingPower nonnegative and constant within each key, the manifest's signed "
                    f"raw flux when lit and exactly 0 when dark (|raw flux| <= {DARK_REL:g} A I, OR-038), t_flux * A I; "
                    "|Re S| and flux / |Re S| as recorded; no T <= 1 gate", problems,
                    f"{len(raw)} bases, {n_dark} dark")

    def _outgoing_probe(self) -> Check:
        tables, keys, problems = self._usable(), self._keys(), []
        for e in self.pols:
            rows = self._rows("waveguide", e)
            for k, key in enumerate(keys):
                stored = tables[e]["waveguide"][rows[k], OUTGOING]
                first = float(stored[0])
                if (stored < 0.0).any():
                    problems.append(f"{_where(e, key)}: OutgoingPower {float(stored.min())!r} W is negative")
                if np.any(np.abs(stored - first) > STORED_REL * abs(first)):
                    problems.append(f"{_where(e, key)}: OutgoingPower is not constant within the key (from {first!r} W "
                                    f"to {float(stored[np.argmax(np.abs(stored - first))])!r} W)")
        return _row("powers.outgoing", "schema 1's |Re S| integral: OutgoingPower nonnegative and constant within the "
                    "key; no T <= 1 gate and no flux record (a schema-1 run records none)", problems,
                    f"{len(keys)} key, Ephi {list(self.pols)}")

    # -- 9. intercept_factor ------------------------------------------------------------------------------------------
    def _intercept(self) -> Check:
        tables, keys, thetas, sidecar = self._usable(), self._keys(), self._thetas(), self._sidecar()
        problems: list[str] = []
        cell = dataset_v2.cell_90(thetas)
        band = dataset_v2.cell_90_band_deg(thetas)
        recorded = _number(_get(sidecar, "transmittance", "cell_90", where="the sidecar"),
                           "the sidecar's transmittance.cell_90")
        if not _within(recorded, cell, INTERCEPT_REL):
            problems.append(f"the sidecar's cell_90 {recorded!r} is not {cell!r}, recomputed from the declared thetas")
        assumes = _get(sidecar, "transmittance", "cell_90_assumes", where="the sidecar")
        recorded_band = assumes.get("band_deg") if isinstance(assumes, Mapping) else None
        if not (isinstance(recorded_band, list) and len(recorded_band) == 2
                and all(_real(v) and _within(v, b, INTERCEPT_REL) for v, b in zip(recorded_band, band))):
            problems.append(f"the sidecar's cell_90 band {recorded_band!r} is not {list(band)} deg")
        if not (isinstance(assumes, Mapping) and assumes.get("emitter") == dataset_v2.CELL_90_EMITTER
                and assumes.get("cos_floor") == dataset_v2.COS_FLOOR):
            problems.append(f"the sidecar's cell_90_assumes {assumes!r} is not emitter "
                            f"{dataset_v2.CELL_90_EMITTER!r} with cos_floor {dataset_v2.COS_FLOOR!r} (HF-049)")
        c_records = _by_key(_get(sidecar, "transmittance", "c_key", where="the sidecar"), keys,
                            "the sidecar's transmittance.c_key")
        per_key = self._per_key()
        for k, key in enumerate(keys):
            c = dataset_v2.c_key(key[0], key[1], thetas)
            for what, value in (("transmittance.c_key", _number(c_records[k].get("c_k"), "c_key c_k")),
                                ("transmittance.per_key c_k", per_key[k]["c_k"])):
                if not _within(value, c, INTERCEPT_REL):
                    problems.append(f"at {key_text(key)}: the sidecar's {what} {value!r} is not c_key = {c!r}")
            for e in self.pols:
                column = tables[e]["waveguide"][self._rows("waveguide", e)[k], INTERCEPT]
                bad = ~(np.abs(column - c) <= INTERCEPT_REL * abs(c))
                if bad.any():
                    problems.append(f"{_where(e, key)}: InterceptFactor {float(column[np.flatnonzero(bad)[0]])!r} on "
                                    f"{int(bad.sum())} rows is not c_key = {c!r}")
        return _row("intercept_factor", f"every row's InterceptFactor and the sidecar's c_key and per-key c_k equal "
                    f"dataset_v2.c_key of the declared thetas, cell_90 and its band theirs, to {INTERCEPT_REL:g} "
                    "(OR-034)", problems, f"{len(keys)} keys, cell_90 = {cell!r} over {list(band)} deg")

    # -- 10. exit_field.normalization ---------------------------------------------------------------------------------
    def _normalization(self) -> Check:
        tables, keys, grid, ai = self._usable(), self._keys(), self._grid(), self._ai()
        scale, dark, per_key = self._scale(), self._dark(), self._per_key()
        h_y, h_z = grid.step_mm[1] * 1e-3, grid.step_mm[2] * 1e-3
        problems: list[str] = []
        lit = 0
        for e in self.pols:
            rows = self._rows("waveguide", e)
            for k, key in enumerate(keys):
                where = _where(e, key)
                values = tables[e]["waveguide"][rows[k]]
                w, flag = scale[k].get(f"w{e}"), scale[k].get(f"dark{e}")
                if not isinstance(flag, bool):
                    problems.append(f"{where}: dark{e} must be true or false, got {flag!r}")
                elif flag != dark[(e, k)]:
                    problems.append(f"{where}: dark{e} is {flag}, but the raw flux is {'' if dark[(e, k)] else 'not '}a "
                                    f"numerical zero (|flux| <= {DARK_REL:g} A I)")
                if dark[(e, k)]:
                    if not (_real(w) and w == 1.0):
                        problems.append(f"{where}: a dark basis keeps W = 1, the sidecar has w{e} = {w!r}")
                    if np.any(values[:, OUTGOING] != 0.0) or per_key[k]["t_flux"][e] != 0.0:
                        problems.append(f"{where}: a dark basis stores OutgoingPower and t_flux exactly 0")
                    continue
                lit += 1
                if not _positive(w):
                    problems.append(f"{where}: W = {w!r} is not a finite positive number")
                try:
                    weights = self._weights(e, k)
                    w_check, is_dark = dataset_v2.scale_factor(values[:, FIELDS], weights, h_y, h_z,
                                                               float(values[0, OUTGOING]), ai)
                except (_Refused, DatasetV2Error) as exc:
                    problems.append(f"{where}: {exc}")
                    continue
                if is_dark:
                    problems.append(f"{where}: OutgoingPower {float(values[0, OUTGOING])!r} W of a lit basis is a "
                                    "numerical zero")
                elif not abs(w_check ** 2 - 1.0) <= NORMALIZATION_REL:
                    problems.append(f"{where}: the stored fields carry {w_check ** 2:.12g} times OutgoingPower, beyond "
                                    f"{NORMALIZATION_REL:g}")
        return _row("exit_field.normalization", f"lit bases: W finite positive and the stored scaled fields' "
                    f"sum w |E_t|^2 h_y h_z / (2 eta0) equal to OutgoingPower to {NORMALIZATION_REL:g} "
                    "(dataset_v2.scale_factor, W never undone); dark bases: W = 1, OutgoingPower and t_flux 0, the flag "
                    "set (HF-040, OR-038)", problems, f"{lit} lit bases, {len(dark) - lit} dark")

    # -- 11. sampling.normalization -----------------------------------------------------------------------------------
    def _sampling(self) -> Check:
        tables, keys, dark, manifest = self._usable(), self._keys(), self._dark(), self._manifest()
        ff_over = self._per_basis(_get(self._sidecar(), "far_field", "ff_over_flux", where="the sidecar"),
                                  "the sidecar's far_field.ff_over_flux")
        rho_records = _by_key(_get(manifest, "transmittance_checks", "rho_per_key", where="the manifest"), keys,
                              "the manifest's transmittance_checks.rho_per_key")
        steps = [_number(_get(manifest, "far_field", name, where="the manifest"), f"the manifest's far_field.{name}")
                 for name in ("theta_step_deg", "phi_step_deg")]
        problems: list[str] = []
        for k, key in enumerate(keys):
            for e in self.pols:
                where, recorded = _where(e, key), ff_over[(e, k)]
                if dark[(e, k)]:
                    if recorded is not None:
                        problems.append(f"{where}: a dark basis's ff_over_flux is null, not {recorded!r}")
                    continue
                values = tables[e]["waveguide"][self._rows("waveguide", e)[k]]
                try:
                    position = dataset_v2._weighted_power(values[:, FIELDS], self._weights(e, k))
                except (_Refused, DatasetV2Error) as exc:
                    problems.append(f"{where}: {exc}")
                    continue
                far = tables[e]["far_field"][self._rows("far_field", e)[k], FAR_POINT]
                hemisphere = dataset_v2.hemisphere_power_w(
                    pd.DataFrame(far, columns=schema.FAR_FIELD_COLUMNS[FAR_POINT]), *steps)
                for name, value in (("position", position), ("direction", hemisphere)):
                    if not (math.isfinite(value) and value > 0.0):
                        problems.append(f"{where}: the {name} normalizer is {value!r}, not finite and positive")
                expected = hemisphere / float(values[0, OUTGOING])
                if recorded is None or not _within(recorded, expected, FF_OVER_FLUX_REL):
                    problems.append(f"{where}: ff_over_flux {recorded!r} is not the hemisphere power over the stored "
                                    f"flux, {expected!r}")
            if len(self.pols) == 2:
                try:
                    stored = [tables[e]["waveguide"][self._rows("waveguide", e)[k], FIELDS] for e in (0, 1)]
                    value = dataset_v2.rho(stored[0], stored[1], self._weights(0, k), dark[(0, k)], dark[(1, k)])
                except (_Refused, DatasetV2Error) as exc:
                    problems.append(f"at {key_text(key)}: rho: {exc}")
                    continue
                recorded_rho = _number(rho_records[k].get("rho"), f"the manifest's rho_per_key at {key_text(key)}")
                if not abs(recorded_rho - value) <= RHO_ABS:
                    problems.append(f"at {key_text(key)}: the manifest's rho {recorded_rho!r} is not {value!r}, "
                                    f"recomputed on the stored fields (beyond {RHO_ABS:g})")
        return _row("sampling.normalization", "per lit basis the position normalizer (sum w |E_t|^2) and the direction "
                    "normalizer (dataset_v2.hemisphere_power_w) finite and positive; ff_over_flux = hemisphere power / "
                    f"stored flux to {FF_OVER_FLUX_REL:g}, null when dark; dataset_v2.rho equal to the manifest's "
                    f"rho_per_key to {RHO_ABS:g} (HF-042)", problems,
                    f"{len(keys)} keys; the far field stays a pattern")

    # -- 12. convergence.records --------------------------------------------------------------------------------------
    def _convergence_records(self) -> Check:
        cfg, manifest = self._cfg(), self._manifest()
        target = cfg.solver.max_delta_e
        fields = acceptance_gate.CONVERGENCE_RECORD_FIELDS
        solver = manifest.get("solver")
        records = solver.get("convergence_qualification") if isinstance(solver, Mapping) else None
        problems: list[str] = []
        recorded_target = solver.get("max_delta_e") if isinstance(solver, Mapping) else None
        if not _same(recorded_target, target):
            problems.append(f"the manifest's solver.max_delta_e {recorded_target!r} is not the frozen config's {target!r}")
        for e, q in sorted(self._qualifications().items()):
            if isinstance(q, Exception):
                problems.append(str(q))
                continue
            at = f"the Ephi={e} convergence table"
            if q["criterion"] != CRITERION:
                problems.append(f"{at} gives the criterion {q['criterion']!r}, not {CRITERION!r}: never waivable")
            if q["target"] is None or not math.isclose(q["target"], target, rel_tol=1e-9):
                problems.append(f"{at} gives the target {q['target']!r}, not the frozen config's max_delta_e "
                                f"{target!r}: never waivable")
            if q["achieved"] is None or q["achieved"] < 0.0:
                problems.append(f"{at} gives no finite nonnegative achieved delta ({'; '.join(q['reasons'])})")
            for name in ("converged", "passes_completed", "max_passes"):
                if q[name] is None:
                    problems.append(f"{at} has no readable {name}")
            projection = {name: q[name] for name in fields}
            record = records.get(str(e)) if isinstance(records, Mapping) else None
            if not isinstance(record, Mapping):
                problems.append(f"the manifest records no solver.convergence_qualification for Ephi={e}")
                continue
            diffs = [f"{name} {record.get(name, 'missing')!r}, the table gives {projection[name]!r}"
                     for name in fields if name not in record or not _same(record[name], projection[name])]
            if diffs:
                problems.append(f"Ephi={e}: the manifest's record disagrees with its table: {'; '.join(diffs)}")
        problems += acceptance_gate.convergence_problems(manifest, allow_unqualified=True,
                                                         expected_polarizations=self.pols)
        return _row("convergence.records", f"each configured polarization's table qualified at the frozen config's "
                    f"max_delta_e {target:g} (CDX-025) equals the manifest's record field by field, with criterion "
                    f"{CRITERION!r} and that target; the records and solver.qualified consistent; never waivable",
                    list(dict.fromkeys(problems)), f"Ephi {list(self.pols)} records agree with their tables")

    # -- 13. convergence ----------------------------------------------------------------------------------------------
    def _convergence(self) -> Check:
        target = self._cfg().solver.max_delta_e
        qualifications = self._qualifications()
        problems, qualified, finite = [], {}, True
        for e, q in sorted(qualifications.items()):
            usable = not isinstance(q, Exception)
            qualified[str(e)] = usable and q["qualified"]
            if not qualified[str(e)]:
                finite = finite and usable and q["achieved"] is not None
                why = str(q) if not usable else "; ".join(q["reasons"])
                problems.append(f"Ephi={e} does not qualify at the target {target:g}: {why}")
        waivable = bool(self.verdicts.get("convergence.records") is True and finite)
        return _row("convergence", f"every configured polarization qualifies at the frozen config's max_delta_e "
                    f"{target:g} (OR-027); waivable only beside passed convergence.records, each unqualified "
                    "polarization with a finite achieved delta", problems, "every configured polarization qualifies",
                    {"waivable": waivable, "qualified": qualified})

    # -- 14. exit_grid ------------------------------------------------------------------------------------------------
    def _exit_grid(self) -> Check:
        cfg, manifest, sidecar, task = self._cfg(), self._manifest(), self._sidecar(), self.task
        problems: list[str] = []
        limit = frequency_grid.max_l_step_mm(self.frequency)
        effective = _number(task.get("effective_l_resolution_mm"), "the plan task's effective_l_resolution_mm")
        requested = _optional_number(task.get("l_resolution_mm"), "the plan task's l_resolution_mm")
        own = build_config(self._raw_config()).exit_field.resolution_mm[1]
        if effective != (own if requested is None else requested):
            problems.append(f"the plan task's effective_l_resolution_mm {effective!r} is not "
                            + (f"its requested l_resolution_mm {requested!r}" if requested is not None else
                               f"the frozen config's own step along l {own!r} (no step requested)"))
        if effective > limit:
            problems.append(f"the effective step along l {effective!r} mm exceeds c/(4f) = {limit:.7g} mm (OR-048)")
        for label, block in (("manifest", manifest.get("exit_field")), ("sidecar", sidecar.get("exit_field"))):
            resolution = block.get("resolution_mm") if isinstance(block, Mapping) else None
            if not (isinstance(resolution, list) and len(resolution) == 3 and _same(resolution[1], effective)):
                problems.append(f"the {label}'s exit_field.resolution_mm {resolution!r} does not apply the effective "
                                f"step along l {effective!r} mm")
        lattice = _plan_lattice(task)
        l_step = _number(task.get("l_step_mm"), "the plan task's l_step_mm")
        if not _within(l_step, lattice.step[1], NODE_TOL):
            problems.append(f"the plan task's l_step_mm {l_step!r} is not its lattice's step along l "
                            f"{lattice.step[1]!r}")
        if lattice.step[1] > limit:
            problems.append(f"the plan's fitted step along l {lattice.step[1]!r} mm exceeds c/(4f) = {limit:.7g} mm")
        nominal = _nominal_lattice(cfg)
        problems += self._lattice_problems("the frozen rectangle's HF-045 lattice", nominal["lo_mm"], nominal["hi_mm"],
                                           nominal["step_mm"], nominal["counts"], lattice, axes=(0, 1, 2))
        if tuple(nominal["fitted"]) != lattice.fitted:
            problems.append(f"the plan's lattice fitted {list(lattice.fitted)}, the HF-045 rule gives "
                            f"{list(nominal['fitted'])}")
        if self._shape() != "rectangle":
            problems.append(f"the exit cross-section is {self._shape()!r}: the plan predicts a rectangle")
        try:
            grid = self._grid()
        except _Refused as exc:
            problems.append(str(exc))
            grid = None
        if grid is not None:
            problems += self._lattice_problems("the sidecar's declared lattice", grid.lo_mm, grid.hi_mm, grid.step_mm,
                                               grid.counts, lattice, axes=(1, 2))
            if grid.step_mm[1] > limit:
                problems.append(f"the declared step along l {grid.step_mm[1]!r} mm exceeds c/(4f) = {limit:.7g} mm")
        exit_field = manifest.get("exit_field") if isinstance(manifest.get("exit_field"), Mapping) else {}
        step_mm = exit_field.get("step_mm")
        if not (isinstance(step_mm, list) and len(step_mm) == 3 and all(_real(v) for v in step_mm) and all(
                abs(v - s) <= _node_tol_mm(s) for v, s in zip(step_mm, lattice.step))):
            problems.append(f"the manifest's exit_field.step_mm {step_mm!r} is not the plan's {list(lattice.step)}")
        points = lattice.counts[1] * lattice.counts[2]
        for name in ("lattice_points", "retained_points_per_key"):
            value = exit_field.get(name)
            if type(value) is not int or value != points:
                problems.append(f"the manifest's exit_field.{name} {value!r} is not the lattice's {points}")
        if grid is not None:
            for e in self.pols:
                for k in range(len(self._keys())):
                    try:
                        self._weights(e, k)
                    except _Refused as exc:
                        problems.append(str(exc))
        problems = list(dict.fromkeys(problems))
        return _row("exit_grid", f"the effective step along l {effective!r} mm applied, at most c/(4f) = {limit:.7g} mm "
                    "(OR-048); the plan's lattice the frozen rectangle's HF-045 lattice; the declared lattice and the "
                    f"manifest's step_mm the plan's (to {NODE_TOL:g} of the step); every key's exit points on it, each "
                    "node once", problems,
                    f"{lattice.counts[1]} x {lattice.counts[2]} nodes, step along l {lattice.step[1]!r} mm")

    @staticmethod
    def _lattice_problems(what: str, lo: Sequence[float], hi: Sequence[float], step: Sequence[float],
                          counts: Sequence[int], lattice: _Lattice, axes: Sequence[int]) -> list[str]:
        problems = []
        for axis in axes:
            name = ("x_e", "y_e", "z_e")[axis]
            tol = _node_tol_mm(lattice.step[axis])
            if not (abs(lo[axis] - lattice.lo[axis]) <= tol and abs(hi[axis] - lattice.hi[axis]) <= tol):
                problems.append(f"{what} spans {lo[axis]!r} to {hi[axis]!r} mm along {name}, the plan's lattice "
                                f"{lattice.lo[axis]!r} to {lattice.hi[axis]!r} mm")
            if not abs(step[axis] - lattice.step[axis]) <= NODE_TOL * lattice.step[axis]:
                problems.append(f"{what} steps {step[axis]!r} mm along {name}, the plan's lattice "
                                f"{lattice.step[axis]!r} mm")
            if counts[axis] != lattice.counts[axis]:
                problems.append(f"{what} has {counts[axis]} nodes along {name}, the plan's lattice "
                                f"{lattice.counts[axis]}")
        return problems

    # -- 15. binding.inputs -------------------------------------------------------------------------------------------
    def _binding_inputs(self) -> Check:
        problems: list[str] = []
        try:
            sidecar = self._sidecar()
        except _Refused:
            sidecar = {}
        files = acceptance_gate.evidence_files("array_run", self.job, None, sidecar, context=self.context,
                                               plan=self.plan)
        for path in files:
            try:
                self.inputs.digest(path)
            except OSError as exc:
                problems.append(f"{path} is required evidence and is missing or unreadable ({exc})")
            except acceptance_gate.InputChanged:
                pass                                     # in self.inputs.problems
        problems = list(dict.fromkeys(self.inputs.problems + problems))
        return _row("binding.inputs", "every file read once and unchanged until the end of the check, every evidence "
                    "file present and hashed", problems, f"{len(self.inputs.hashes)} files bound")


# --- entry points ----------------------------------------------------------------------------------------------------

def _read_plan(inputs: acceptance_gate.InputCapture, plan_path: Path, plan_sha256: object,
               index: object) -> tuple[dict, Mapping, str]:
    """The plan, read once through ``inputs``, task ``index``'s record and the plan's policy; ``PlanRefused`` for a
    plan that is unreadable, has another sha256 than the caller expects, or names no such task."""
    if not _is_sha256(plan_sha256):
        raise PlanRefused(f"the expected plan sha256 must be 64 lowercase hex digits, got {plan_sha256!r}")
    if isinstance(index, bool) or not isinstance(index, int) or index < 0:
        raise PlanRefused(f"a task index is a nonnegative integer, got {index!r}")
    try:
        plan = inputs.parse(plan_path, json.load)
    except (OSError, ValueError) as exc:
        raise PlanRefused(f"the plan {plan_path} is unreadable ({exc})") from None
    digest = inputs.recorded(plan_path)
    if digest != plan_sha256:
        raise PlanRefused(f"the plan {plan_path} has sha256 {digest}, not the {plan_sha256} the caller expects")
    if not isinstance(plan, dict):
        raise PlanRefused(f"the plan {plan_path} is not a JSON object")
    if plan.get("schema") != array_plan.PLAN_SCHEMA:
        raise PlanRefused(f"the plan's schema is {plan.get('schema')!r}, not {array_plan.PLAN_SCHEMA!r}")
    policy = plan.get("policy")
    if policy not in POLICIES:
        raise PlanRefused(f"the plan's policy is {policy!r}, not one of {', '.join(POLICIES)}")
    dataset_id = plan.get("dataset_id")
    if not isinstance(dataset_id, str) or not dataset_id:
        raise PlanRefused(f"the plan's dataset_id must be a non-empty string, got {dataset_id!r}")
    tasks = plan.get("tasks")
    found = ([t for t in tasks if isinstance(t, Mapping) and type(t.get("index")) is int and t["index"] == index]
             if isinstance(tasks, list) else [])
    if len(found) != 1:
        raise PlanRefused(f"the plan holds {len(found)} records of task {index}, not 1")
    if not _positive(found[0].get("frequency_ghz")):
        raise PlanRefused(f"task {index}'s frequency_ghz must be a finite positive number, got "
                          f"{found[0].get('frequency_ghz')!r}")
    return plan, found[0], policy


def _write_atomic(path: Path, data: Mapping) -> None:
    """``data`` as JSON at ``path``, through a temporary file in the same directory and ``os.replace``."""
    tmp = path.with_name(f".{path.name}.tmp")
    try:
        write_json(tmp, data)
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()


def write_report(result: RunAcceptanceResult) -> None:
    """Write ``result.to_dict()`` to ``result.report_target`` atomically. An OSError is recorded in
    ``result.report_error``, which fails the run, and never raised."""
    if result.report_target is None:
        return
    try:
        _write_atomic(Path(result.report_target), result.to_dict())
    except OSError as exc:
        result.report_error = f"{result.report_target} could not be written ({exc}); the run fails"
        return
    result.report = result.report_target


def check(job: Path, plan_path: Path, plan_sha256: str, index: int, *, write: bool = True) -> RunAcceptanceResult:
    """Judge task ``index``'s job ``job`` against the plan at ``plan_path``, which must have the sha256 the caller
    recorded, and with ``write`` write ``<job>/array_run.json`` (``write_report``).

    Raises ``FileNotFoundError`` when the job directory is missing and ``PlanRefused`` when the plan is unreadable,
    is not the expected one or names no such task (usage errors); every problem with the job's contents becomes a
    failing row instead."""
    job = Path(job)
    if not job.is_dir():
        raise FileNotFoundError(f"the job directory {job} does not exist")
    plan_path = Path(os.path.abspath(plan_path))
    inputs = acceptance_gate.InputCapture()
    plan, task, policy = _read_plan(inputs, plan_path, plan_sha256, index)
    context = acceptance_gate.PlanContext(str(plan_path.resolve()), plan_sha256, index, policy)
    acceptance = _Acceptance(job, plan, task, context, inputs)
    checks = acceptance.run()
    result = RunAcceptanceResult(str(job), {"path": context.plan_path, "sha256": plan_sha256, "index": index,
                                            "policy": policy}, checks,
                                 binding=acceptance_gate.binding_block("array_run", job, None, inputs.hashes),
                                 report_target=str(job / REPORT_NAME))
    if write:
        write_report(result)
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m bbsim.run_acceptance", description=__doc__.splitlines()[0])
    parser.add_argument("--plan", type=Path, required=True, help="the job array's frozen plan.json")
    parser.add_argument("--plan-sha256", required=True, help="the plan's sha256 the caller recorded (plan.sha256)")
    parser.add_argument("--index", type=int, required=True, help="the task's index in the plan")
    parser.add_argument("job", type=Path, help="the task's job directory")
    try:
        args = parser.parse_args(argv)
    except SystemExit as exc:             # argparse has printed its message
        return 2 if exc.code else 0
    try:
        result = check(args.job, args.plan, args.plan_sha256, args.index, write=False)
    except (FileNotFoundError, PlanRefused) as exc:
        print(f"run_acceptance: refused: {exc}", file=sys.stderr)
        return 2
    print(result.format_rows(), flush=True)     # the table reaches the log before the report is written
    write_report(result)
    print(result.format_footer())
    return 0 if result.passed else 1


if __name__ == "__main__":
    sys.exit(main())
