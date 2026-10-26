"""A pointwise, sign-sensitive comparison of two solves of one geometry and frequency (CDX-026, S42, HF-054).

``compare.py --pair`` judges a dataset's two polarizations by integrated quantities (the overlap Re rho and the mixed
transmission), so a relative sign flip of Ephi=1 at a key whose integrated overlap is zero passes it, yet reverses
BBRsim's mixed-state exit-position and direction draws. This module compares the local quantities instead: per shared
key, for the exit field and the far field, with V_e,i the field vector of polarization e at point i and w_i its weight,

    q_e,i = w_i |V_e,i|^2,  c_i = w_i V_0,i . conj(V_1,i),  Q = N_0 + N_1 with N_e = sum_i q_e,i,
    g_e,i = q_e,i / Q,      h_i = c_i / Q,

and compares g_0, g_1 and the complex h of a new solve point by point with an earlier solve's. For any mixture
a V_0 + b V_1 the point density over Q is |a|^2 g_0 + |b|^2 g_1 + 2 Re(a conj(b) h), so these three keep the mixture's
distribution; a relative sign negates h, a common phase or a common positive scale leaves all three unchanged.

Conventions (recorded with every result):

- ``legacy_all_components``: two schema-1.x solves; the three complex exit components read by position (HF-038), weight
  1 on every point.
- ``raw_hfss_tangential``: either side schema 2.x; both sides' underlying HFSS (Ex, Ey), a 2.x table's per-key W undone
  (``compare.legacy_view``), on the declared lattice: trapezoid node weights for a rectangle (``dataset_v2.node_weights``
  on the schema-2 declaration, after point identity), 1 on every retained disc point (zero-filled rim points included,
  HF-029). A schema-1 sidecar's ``exit_field.grid`` has min, max and count but no step; its step is (max - min) /
  (count - 1) per axis, and the lattice must equal the schema-2 declaration to 1e-12 m with the same counts.
- ``outward_sin_theta``: the far field, V = (rEphi, rEtheta), weight sin(Theta) where sin(Theta) cos(Phi) >= 1e-6
  (BBRsim's outward guard) and 0 elsewhere; W never applies to it.

Each side is read through one ``acceptance_gate.InputCapture``, one complete table at a time, and validated in its own
schema over every declared key before anything is selected: finite values, the declared keys and nothing else, one
polarization per directory, equal rows per key, and for 2.x the sidecar's records, InterceptFactor and OR-038's dark
rule. Only then are the compared keys' rows kept. Lattices are checked, never interpolated: exit points equal to
1e-12 m and far-field angles to 1e-9 deg in the same row order, within each side and between the sides.
"""
from __future__ import annotations

import json
from collections.abc import Collection, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

from bbsim import compare, dataset_v2
from bbsim.compare import Check, ComparisonInputError, Tolerances
from bbsim.dataset_v2 import DatasetV2Error, key_text
from bbsim.schema import EXIT_FIELD_SLOTS, FAR_FIELD_COLUMNS, WAVEGUIDE_COLUMNS, WAVEGUIDE_COLUMNS_V2

if TYPE_CHECKING:
    from bbsim.acceptance_gate import InputCapture
    from bbsim.extract import ExitGrid

POLICY = "bbsim-pointwise-pair/1"
MODE = "measure"                    # "measure": report only; "enforce": judge against TOL
TOL: float | None = None            # set from the measured October values when enforcement starts
SUPPORT_MIN = dataset_v2.OUTWARD_GUARD    # far-field support: sin(Theta) cos(Phi) >= this
COORD_M = 1e-12                     # exit-point identity, m
ANGLE_DEG = 1e-9                    # far-field angle identity, deg
DARK_FLOOR = Tolerances.t_noise_floor     # the dark/lit band of compare.py's pair check, on m = |Re S| / (A I)
DARK_RATIO = Tolerances.pair_dark_ratio
LIT_RATIO = Tolerances.pair_lit_ratio
EXIT_1_1 = "legacy_all_components"
EXIT_RAW = "raw_hfss_tangential"
FAR = "outward_sin_theta"
SCHEMA_1, SCHEMA_2 = compare.SCHEMA_1, compare.SCHEMA_2
_EXIT_COORDS = ["X", "Y", "Z"]
_FAR_COORDS = ["Theta", "Phi"]


def policy() -> dict[str, object]:
    """The comparison policy a grid acceptance report records and assembly expects, read at call time."""
    return {"fields_pointwise_policy": POLICY, "fields_pointwise_mode": MODE, "fields_pointwise_tol": TOL,
            "fields_pointwise_support_min": SUPPORT_MIN, "fields_pointwise_coord_m": COORD_M,
            "fields_pointwise_angle_deg": ANGLE_DEG, "fields_pointwise_dark_floor": DARK_FLOOR,
            "fields_pointwise_dark_ratio": DARK_RATIO, "fields_pointwise_lit_ratio": LIT_RATIO,
            "fields_pointwise_exit_conventions": [EXIT_1_1, EXIT_RAW], "fields_pointwise_far_convention": FAR}


class PairingRefused(ValueError):
    """An input the pointwise comparison cannot use as given: it is refused, never repaired or guessed."""


@dataclass
class Side:
    """One dataset (its Ephi=0 and Ephi=1 directories), validated in its own schema. ``exit`` holds the 1.0 column
    labels with a 2.x table's W undone; ``exit`` and ``far`` hold only the kept keys' rows, while ``declared``, ``m``
    and ``dark_flag`` cover every declared key. ``m`` is |Re S| / (A I) per polarization from the key's first row
    (schema 1 OutgoingPower / IngoingPower, schema 2 ``outgoing_power_abs_re_s_w`` / IngoingPower), before c_k or the
    flux substitution. ``shape`` and ``grid`` are the declared exit lattice, None when the side declares none or (schema
    1, ``lattice_refusal``) declares one that cannot be used; that is refused only by a comparison that needs it."""

    label: str
    dirs: tuple[Path, Path]
    schema: str
    sidecar: dict | None
    declared: list[tuple[float, float]]
    exit: dict[int, pd.DataFrame]
    far: dict[int, pd.DataFrame]
    m: dict[tuple[float, float], tuple[float, float]]
    dark_flag: dict[tuple[float, float], tuple[bool, bool]]
    shape: str | None = None
    grid: ExitGrid | None = None
    lattice_refusal: str | None = None     # why a schema-1 side's declared lattice is unusable, refused when needed


# --- reading ---------------------------------------------------------------------------------------------------------

def job_dataset_dirs(job: Path) -> tuple[Path, Path]:
    """The (Ephi=0, Ephi=1) directories of a job root's one dataset, named by its one ``<stem>.dataset.json``."""
    job = Path(job)
    sidecars = sorted(job.glob(f"*{compare.SIDECAR_SUFFIX}")) if job.is_dir() else []
    if len(sidecars) != 1:
        raise PairingRefused(f"{job} holds {len(sidecars)} dataset sidecars (*{compare.SIDECAR_SUFFIX}), not exactly "
                             "one: a job root holds one dataset")
    stem = sidecars[0].name[:-len(compare.SIDECAR_SUFFIX)]
    dirs = (job / f"{stem}_Ephi=0", job / f"{stem}_Ephi=1")
    missing = [str(d) for d in dirs if not d.is_dir()]
    if missing:
        raise PairingRefused(f"{job}: the sidecar {sidecars[0].name} has no dataset directory {', '.join(missing)}")
    return dirs


def _refuse(where: str, exc: Exception) -> PairingRefused:
    return PairingRefused(f"{where}: {exc}")


def declared_keys(sidecar: Mapping, where: str) -> list[tuple[float, float]]:
    """The product grid of ``excitation.incident_phi_deg`` x ``incident_theta_deg``, phi outer."""
    axes = []
    for name in ("incident_phi_deg", "incident_theta_deg"):
        try:
            values = dataset_v2.require(sidecar, "excitation", name, where=where)
            if not isinstance(values, list) or not values:
                raise DatasetV2Error(f"{where} excitation.{name} must be a non-empty list, got {values!r}")
            numbers = [dataset_v2.finite_number(v, f"{where} excitation.{name}") for v in values]
        except DatasetV2Error as exc:
            raise PairingRefused(str(exc)) from exc
        ordered = sorted(numbers)
        if any(b - a <= dataset_v2.KEY_TOL_DEG for a, b in zip(ordered, ordered[1:])):
            raise PairingRefused(f"{where} excitation.{name} repeats a value: {numbers}")
        axes.append(numbers)
    return [(p, t) for p in axes[0] for t in axes[1]]


def _schema1_grid(record: object, where: str) -> ExitGrid:
    """A schema-1 ``exit_field.grid`` (min, max, count per axis; a step only if declared) as a ``declared_grid``,
    the step derived as (max - min) / (count - 1)."""
    full: dict[str, dict] = {}
    for axis in ("y_e", "z_e"):
        at = f"{where}.{axis}"
        try:
            rec = dataset_v2.require(record, axis, where=where)
            lo = dataset_v2.finite_number(dataset_v2.require(rec, "min", where=at), f"{at}.min")
            hi = dataset_v2.finite_number(dataset_v2.require(rec, "max", where=at), f"{at}.max")
            count = dataset_v2.require(rec, "count", where=at)
        except DatasetV2Error as exc:
            raise PairingRefused(str(exc)) from exc
        if isinstance(count, (bool, np.bool_)) or not isinstance(count, (int, np.integer)) or count < 2:
            raise PairingRefused(f"{at}.count is {count!r}: a declared lattice needs an integer count of at least 2")
        full[axis] = {"min": lo, "max": hi, "count": int(count),
                      "step": rec["step"] if "step" in rec else (hi - lo) / (int(count) - 1)}
    try:
        return dataset_v2.declared_grid(full, where)
    except DatasetV2Error as exc:
        raise PairingRefused(str(exc)) from exc


def _lattice(sidecar: dict | None, schema: str, where: str) -> tuple[str | None, ExitGrid | None, str | None]:
    """(shape, grid, refusal): a 2.x side's declared lattice, refused now; a 1.x side's, whose refusal is kept for a
    comparison that needs it (two 1.x sides never do)."""
    if sidecar is None:
        return None, None, None
    if schema == SCHEMA_2:
        try:
            return (*compare._v2_lattice(sidecar, where), None)
        except ComparisonInputError as exc:
            raise PairingRefused(str(exc)) from exc
    exit_field = sidecar.get("exit_field")
    if not isinstance(exit_field, Mapping):
        return None, None, None
    section = exit_field.get("cross_section")
    shape = section.get("shape") if isinstance(section, Mapping) else None
    if shape is not None and shape not in dataset_v2.QUADRATURE_BY_SHAPE:
        return None, None, (f"{where} exit_field.cross_section.shape {shape!r} is not one of "
                            f"{sorted(dataset_v2.QUADRATURE_BY_SHAPE)}")
    if "grid" not in exit_field:
        return shape, None, None
    try:
        return shape, _schema1_grid(exit_field["grid"], f"{where} exit_field.grid"), None
    except PairingRefused as exc:
        return shape, None, str(exc)


def _key_rows(table: pd.DataFrame, keys: Sequence[tuple[float, float]], where: str) -> np.ndarray:
    """Every row at exactly one declared key, every declared key with rows, every key the same number of rows."""
    try:
        index = dataset_v2.row_keys(table, keys, where)
    except DatasetV2Error as exc:
        raise PairingRefused(str(exc)) from exc
    counts = np.bincount(index, minlength=len(keys))
    if np.any(counts != counts[0]):
        k = int(np.flatnonzero(counts != counts[0])[0])
        raise PairingRefused(f"{where} has {int(counts[k])} rows at {key_text(keys[k])} but {int(counts[0])} at "
                             f"{key_text(keys[0])}: every key holds the same points")
    return index


def _kept(table: pd.DataFrame, index: np.ndarray, positions: Sequence[int] | None) -> pd.DataFrame:
    if positions is None:
        return table.reset_index(drop=True)
    return table[np.isin(index, positions)].reset_index(drop=True)


def _read_csv(stream, where: str) -> pd.DataFrame:
    try:
        return pd.read_csv(stream)
    except ValueError as exc:                            # pandas' parser errors and undecodable bytes
        raise _refuse(where, exc) from exc


def _schema_of(columns: list[str], where: str) -> str:
    if columns == list(WAVEGUIDE_COLUMNS):
        return SCHEMA_1
    if columns == list(WAVEGUIDE_COLUMNS_V2):
        return SCHEMA_2
    raise PairingRefused(f"{where} carries neither the 1.0 nor the 2.0 column list: {columns}")


def read_side(label: str, dirs: tuple[Path, Path], inputs: InputCapture,
              keep: Collection[tuple[float, float]] | None = None) -> Side:
    """Read and validate one dataset's two polarizations through ``inputs``, keeping only ``keep``'s rows (all rows
    when None). Raises ``PairingRefused`` naming the file and key of the first problem."""
    dirs = (Path(dirs[0]), Path(dirs[1]))
    paths = [compare.sidecar_path(d) for d in dirs]
    if paths[0] != paths[1]:
        raise PairingRefused(f"{label}: {dirs[0]} and {dirs[1]} are not one dataset's two polarizations (sidecars "
                             f"{paths[0]} and {paths[1]})")
    sidecar_path = paths[0]
    sidecar: dict | None = None
    where_sidecar = f"{label} sidecar {sidecar_path}"
    if sidecar_path is not None:
        def read_sidecar(stream):
            try:
                return json.load(stream)
            except ValueError as exc:                    # malformed JSON or bytes
                raise _refuse(where_sidecar, exc) from exc

        try:
            sidecar = inputs.parse(sidecar_path, read_sidecar)
        except OSError as exc:
            raise _refuse(where_sidecar, exc) from exc
        if not isinstance(sidecar, dict):
            raise PairingRefused(f"{where_sidecar} is not a JSON object")
    declared = declared_keys(sidecar, where_sidecar) if sidecar is not None else None
    state: dict[str, Any] = {"schema": None, "declared": declared, "positions": None}

    def positions(where: str) -> list[int] | None:
        if keep is None:
            return None
        found = []
        for key in keep:
            k = dataset_v2.key_position(state["declared"], *key)
            if k is None:
                raise PairingRefused(f"{where}: the compared key {key_text(key)} is not declared")
            found.append(k)
        return found

    exit_tables: dict[int, pd.DataFrame] = {}
    far_tables: dict[int, pd.DataFrame] = {}
    m: dict[int, list[float]] = {}

    for e, d in enumerate(dirs):
        path = d / "waveguide.csv"
        where = f"{label} Ephi={e} {path}"

        def read_waveguide(stream, e=e, d=d, where=where):
            df = _read_csv(stream, where)
            schema = _schema_of([str(c) for c in df.columns], where)
            if state["schema"] is None:
                state["schema"] = schema
            elif schema != state["schema"]:
                raise PairingRefused(f"{where} is schema {schema} but Ephi=0's table is {state['schema']}")
            if sidecar is not None and sidecar.get("schema_version") != schema:
                raise PairingRefused(f"{where} carries the {schema} column list but its sidecar declares schema "
                                     f"{sidecar.get('schema_version')!r}")
            try:
                problems = compare._non_finite(where, df)
                if problems:
                    raise ComparisonInputError("; ".join(problems))
                if compare._ephi_of(df, d, where) != e:
                    raise ComparisonInputError(f"{where} holds another polarization than Ephi={e}")
                view = compare.legacy_view(df, sidecar, e) if schema == SCHEMA_2 else df
            except ComparisonInputError as exc:
                raise _refuse(where, exc) from exc
            if state["declared"] is None:
                state["declared"] = sorted(compare._keys(view))
            index = _key_rows(view, state["declared"], where)
            incoming = view["IngoingPower"].to_numpy(float)
            outgoing = view["OutgoingPower"].to_numpy(float)
            if not (incoming > 0.0).all() or not (outgoing >= 0.0).all():
                raise PairingRefused(f"{where}: IngoingPower must be positive and |Re S| non-negative")
            first = [int(np.flatnonzero(index == k)[0]) for k in range(len(state["declared"]))]
            return _kept(view, index, positions(where)), [float(outgoing[i] / incoming[i]) for i in first]

        try:
            exit_tables[e], m[e] = inputs.parse(path, read_waveguide)
        except OSError as exc:
            raise _refuse(where, exc) from exc

    for e, d in enumerate(dirs):
        path = d / "far_field.csv"
        where = f"{label} Ephi={e} {path}"

        def read_far(stream, e=e, d=d, where=where):
            df = _read_csv(stream, where)
            if [str(c) for c in df.columns] != list(FAR_FIELD_COLUMNS):
                raise PairingRefused(f"{where} does not carry the far-field columns {FAR_FIELD_COLUMNS}: "
                                     f"{list(df.columns)}")
            try:
                problems = compare._non_finite(where, df)
                if problems:
                    raise ComparisonInputError("; ".join(problems))
                if compare._ephi_of(df, d, where) != e:
                    raise ComparisonInputError(f"{where} holds another polarization than Ephi={e}")
            except ComparisonInputError as exc:
                raise _refuse(where, exc) from exc
            index = _key_rows(df, state["declared"], where)
            return _kept(df, index, positions(where))

        try:
            far_tables[e] = inputs.parse(path, read_far)
        except OSError as exc:
            raise _refuse(where, exc) from exc

    schema = state["schema"]
    keys = state["declared"]
    dark = {key: (False, False) for key in keys}
    if schema == SCHEMA_2:
        try:
            records = compare._v2_records(sidecar, where_sidecar)
        except ComparisonInputError as exc:
            raise PairingRefused(str(exc)) from exc
        for k, key in enumerate(records.keys):
            at = dataset_v2.key_position(keys, *key)
            if at is not None:
                dark[keys[at]] = (records.dark[0][k], records.dark[1][k])
    shape, grid, lattice_refusal = _lattice(sidecar, schema, where_sidecar)
    return Side(label, dirs, schema, sidecar, list(keys), exit_tables, far_tables,
                {key: (m[0][k], m[1][k]) for k, key in enumerate(keys)}, dark, shape, grid, lattice_refusal)


# --- alignment -------------------------------------------------------------------------------------------------------

def _rows(table: pd.DataFrame, key: tuple[float, float]) -> pd.DataFrame:
    mask = ((np.abs(table["IWavePhi"].to_numpy(float) - key[0]) <= dataset_v2.KEY_TOL_DEG)
            & (np.abs(table["IWaveTheta"].to_numpy(float) - key[1]) <= dataset_v2.KEY_TOL_DEG))
    return table[mask].reset_index(drop=True)


def _at_key(side: Side, key: tuple[float, float]) -> tuple[float, float]:
    """The side's own declared key within KEY_TOL_DEG of ``key``."""
    k = dataset_v2.key_position(side.declared, *key)
    if k is None:
        raise PairingRefused(f"{side.label} does not declare {key_text(key)}")
    return side.declared[k]


def shared_keys(a: Side, b: Side) -> list[tuple[float, float]]:
    """The keys both sides declare, in ``a``'s order; none is a refusal."""
    keys = [key for key in a.declared if dataset_v2.key_position(b.declared, *key) is not None]
    if not keys:
        raise PairingRefused(f"{a.label} and {b.label} share no incidence key: {a.label} declares "
                             f"{len(a.declared)}, {b.label} {len(b.declared)}")
    return keys


def _point_difference(x: pd.DataFrame, y: pd.DataFrame, cols: list[str], tol: float) -> str | None:
    if len(x) == 0 or len(y) == 0:
        return f"{len(x)} against {len(y)} points (a side holds no rows at the key)"
    if len(x) != len(y):
        return f"{len(x)} against {len(y)} points"
    diff = np.abs(x[cols].to_numpy(float) - y[cols].to_numpy(float))
    bad = np.flatnonzero((diff > tol).any(axis=1))
    if bad.size:
        i = int(bad[0])
        return (f"row {i} of the key differs: {cols} {x[cols].iloc[i].tolist()} against {y[cols].iloc[i].tolist()} "
                f"(largest difference {float(diff[i].max()):.3g}, limit {tol:g}; {bad.size} row(s))")
    return None


def _lattice_problems(a: Side, b: Side) -> list[str]:
    """The declared lattice a mixed or 2.x comparison needs: both sides declare a shape and grid, and they agree."""
    if a.schema == SCHEMA_1 and b.schema == SCHEMA_1:
        return []
    problems = []
    for side in (a, b):
        if side.lattice_refusal is not None:
            problems.append(f"{side.label} (schema {side.schema}): {side.lattice_refusal}; the {EXIT_RAW} comparison "
                            "needs both sides' declared lattice")
        elif side.shape is None or side.grid is None:
            problems.append(f"{side.label} (schema {side.schema}) declares no exit_field.grid and cross-section shape: "
                            f"the {EXIT_RAW} comparison needs both sides' declared lattice and never guesses it")
    if problems:
        return problems
    if a.shape != b.shape:
        problems.append(f"the exit cross-sections differ: {a.label} {a.shape}, {b.label} {b.shape}")
    tol_mm = COORD_M * 1e3
    for i, axis in ((1, "y_e"), (2, "z_e")):
        ga, gb = a.grid, b.grid
        if ga.counts[i] != gb.counts[i]:
            problems.append(f"the declared {axis} lattices differ: {ga.counts[i]} nodes in {a.label}, "
                            f"{gb.counts[i]} in {b.label}")
            continue
        for name, va, vb in (("min", ga.lo_mm[i], gb.lo_mm[i]), ("max", ga.hi_mm[i], gb.hi_mm[i]),
                             ("step", ga.step_mm[i], gb.step_mm[i])):
            if abs(va - vb) > tol_mm:
                problems.append(f"the declared {axis} {name} differs: {va * 1e-3!r} m in {a.label}, {vb * 1e-3!r} m "
                                f"in {b.label} (limit {COORD_M:g} m)")
    return problems


def check_alignment(a: Side, b: Side, keys: Sequence[tuple[float, float]]) -> list[str]:
    """Every problem with comparing ``a`` and ``b`` point by point at ``keys``: the declared lattices (for a
    ``raw_hfss_tangential`` comparison), and per key the exit points (1e-12 m) and far-field angles (1e-9 deg), in row
    order, within each side's two polarizations and between the sides. Empty when aligned."""
    problems = _lattice_problems(a, b)
    for key in keys:
        for side in (a, b):
            if dataset_v2.key_position(side.declared, *key) is None:
                problems.append(f"{side.label} does not declare {key_text(key)}")
        for name, cols, tol in (("exit", _EXIT_COORDS, COORD_M), ("far", _FAR_COORDS, ANGLE_DEG)):
            tables = {side.label: [_rows(getattr(side, name)[e], key) for e in (0, 1)] for side in (a, b)}
            for side in (a, b):
                diff = _point_difference(*tables[side.label], cols, tol)
                if diff:
                    problems.append(f"{side.label} Ephi=0 and Ephi=1 {name} points at {key_text(key)}: {diff}")
            for e in (0, 1):
                diff = _point_difference(tables[a.label][e], tables[b.label][e], cols, tol)
                if diff:
                    problems.append(f"{a.label} and {b.label} Ephi={e} {name} points at {key_text(key)}: {diff}")
    return problems


def exit_weights(a: Side, b: Side, key: tuple[float, float]) -> tuple[str, np.ndarray]:
    """(convention, weights) of the exit-field comparison at ``key``, for aligned sides (``check_alignment``)."""
    problems = _lattice_problems(a, b)
    if problems:
        raise PairingRefused("; ".join(problems))
    rows = _rows(a.exit[0], key)
    if a.schema == SCHEMA_1 and b.schema == SCHEMA_1:
        return EXIT_1_1, np.ones(len(rows))
    if a.shape == "disc":
        return EXIT_RAW, np.ones(len(rows))
    grid = (a if a.schema == SCHEMA_2 else b).grid
    try:
        return EXIT_RAW, dataset_v2.node_weights(rows["Y"].to_numpy(float), rows["Z"].to_numpy(float), grid)
    except DatasetV2Error as exc:
        raise PairingRefused(f"{a.label} exit points at {key_text(key)}: {exc}") from exc


# --- densities and deviations ----------------------------------------------------------------------------------------

def far_weights(theta_deg: np.ndarray, phi_deg: np.ndarray) -> np.ndarray:
    """sin(Theta) where sin(Theta) cos(Phi) >= SUPPORT_MIN (an outward direction, BBRsim's kMinNormalComponent), else 0;
    angles in degrees."""
    sin_t = np.sin(np.radians(np.asarray(theta_deg, dtype=float)))
    outward = sin_t * np.cos(np.radians(np.asarray(phi_deg, dtype=float)))
    return np.where(outward >= SUPPORT_MIN, sin_t, 0.0)


def _vectors(table: pd.DataFrame, field: str, convention: str) -> np.ndarray:
    """Complex field vectors per row: the exit slots by position (HF-038), all three components or (Ex, Ey); the far
    field's (rEphi, rEtheta)."""
    if field == "far":
        v = table[["rEphi_real", "rEphi_imag", "rEtheta_real", "rEtheta_imag"]].to_numpy(float)
    else:
        v = table.iloc[:, EXIT_FIELD_SLOTS].to_numpy(float)
    c = v[:, ::2] + 1j * v[:, 1::2]
    return c[:, :2] if (field == "exit" and convention == EXIT_RAW) else c


def _own_convention(side: Side) -> str:
    return EXIT_1_1 if side.schema == SCHEMA_1 else EXIT_RAW


def _field_inputs(side: Side, key: tuple[float, float], field: str, convention: str,
                  weights: np.ndarray | None) -> tuple[dict[int, np.ndarray], np.ndarray]:
    tables = side.far if field == "far" else side.exit
    rows = {e: _rows(tables[e], key) for e in (0, 1)}
    if weights is None:
        if field == "far":
            weights = far_weights(rows[0]["Theta"].to_numpy(float), rows[0]["Phi"].to_numpy(float))
        else:
            weights = exit_weights(side, side, key)[1]
    return {e: _vectors(rows[e], field, convention) for e in (0, 1)}, np.asarray(weights, dtype=float)


def _densities(vectors: Mapping[int, np.ndarray], weights: np.ndarray, bases: Sequence[int]
               ) -> tuple[dict[int, np.ndarray], np.ndarray | None, float]:
    q = {e: weights * np.sum(np.abs(vectors[e]) ** 2, axis=1) for e in bases}
    total = float(sum(np.sum(q[e]) for e in bases))
    g = {e: q[e] / total for e in bases}
    h = (weights * np.sum(vectors[0] * np.conj(vectors[1]), axis=1) / total) if len(bases) == 2 else None
    return g, h, total


def densities(side: Side, key: tuple[float, float], field: str, convention: str | None = None,
              weights: np.ndarray | None = None, bases: Sequence[int] = (0, 1)
              ) -> tuple[dict[int, np.ndarray], np.ndarray | None]:
    """(g per basis, h) of one side at ``key`` for ``field`` "exit" or "far", normalized by Q over ``bases``; h only
    when both bases are included. ``convention`` and ``weights`` default to the side's own (far: ``far_weights``)."""
    vectors, w = _field_inputs(side, key, field, convention or _own_convention(side), weights)
    g, h, _ = _densities(vectors, w, list(bases))
    return g, h


@dataclass
class KeyMetrics:
    """One key and one field ("exit" or "far"): the deviations judged (``d``: "g0", "g1", "h"), those skipped with
    their reason, the normalizations (Q per side, the max of each g_ref, B_ref), the bases forming Q (the same set on
    both sides), each side's m and dark/lit state, whether the key went from dark to lit, and the exit points whose
    zero fill (all six values exactly 0) differs between the two solves."""

    key: tuple[float, float]
    field: str
    convention: str
    m_new: tuple[float, float]
    m_ref: tuple[float, float]
    states: dict[str, str]
    d: dict[str, float]
    skipped: dict[str, str]
    denominators: dict[str, float]
    transition: bool
    bases: list[int]
    zero_fill_differs: list[dict]


@dataclass
class PairingResult:
    """Per shared key and field the metrics, the conventions and both schemas, every input refusal, the coverage
    counts and the largest deviation of each metric."""

    conventions: dict[str, str]
    schemas: dict[str, str]
    keys: list[KeyMetrics]
    refusals: list[str]
    coverage: dict[str, int]
    largest: dict[str, dict]

    def to_dict(self) -> dict:
        records = []
        for km in self.keys:
            records.append({"phi_deg": km.key[0], "theta_deg": km.key[1], "field": km.field,
                            "convention": km.convention, "m_new": list(km.m_new), "m_ref": list(km.m_ref),
                            "states": dict(km.states), "d": dict(km.d), "skipped": dict(km.skipped),
                            "denominators": dict(km.denominators), "transition": km.transition,
                            "bases": list(km.bases), "zero_fill_differs": list(km.zero_fill_differs)})
        return {"conventions": dict(self.conventions), "schemas": dict(self.schemas), "keys": records,
                "refusals": list(self.refusals), "coverage": dict(self.coverage), "largest": dict(self.largest)}


def _state(m: tuple[float, float], flags: tuple[bool, bool]) -> str:
    """compare.py's band on m (dark below the floor or below DARK_RATIO of the larger, lit above LIT_RATIO), dark
    outright when either basis carries an explicit schema-2 dark flag."""
    if any(flags) or max(m) < DARK_FLOOR or min(m) < DARK_RATIO * max(m):
        return "dark"
    return "lit" if min(m) > LIT_RATIO * max(m) else "between"


def _transmitting(m: float, flag: bool) -> bool:
    return not flag and m >= DARK_FLOOR


def _not_transmitting(label: str, m: float, flag: bool) -> str:
    return f"{label}: dark flag" if flag else f"{label}: m {m:.3g} below {DARK_FLOOR:g}"


def _zero_fill(new_rows: Mapping[int, pd.DataFrame], ref_rows: Mapping[int, pd.DataFrame]) -> list[dict]:
    found = []
    for e in (0, 1):
        zn = ~new_rows[e].iloc[:, EXIT_FIELD_SLOTS].to_numpy(float).any(axis=1)
        zr = ~ref_rows[e].iloc[:, EXIT_FIELD_SLOTS].to_numpy(float).any(axis=1)
        for i in np.flatnonzero(zn != zr):
            found.append({"ephi": e, "point": int(i), "y_m": float(new_rows[e]["Y"].iloc[i]),
                          "z_m": float(new_rows[e]["Z"].iloc[i]), "zero_on": "new" if zn[i] else "ref"})
    return found


def _key_field(new: Side, ref: Side, key: tuple[float, float], field: str, convention: str,
               weights: np.ndarray | None, refusals: list[str]) -> KeyMetrics:
    kn, kr = _at_key(new, key), _at_key(ref, key)
    m_new, m_ref = new.m[kn], ref.m[kr]
    flags_new, flags_ref = new.dark_flag[kn], ref.dark_flag[kr]
    states = {"new": _state(m_new, flags_new), "ref": _state(m_ref, flags_ref)}
    km = KeyMetrics(key, field, convention, m_new, m_ref, states, {}, {}, {}, set(states.values()) == {"dark", "lit"},
                    [], [])
    sides = {"new": (new, m_new, flags_new), "ref": (ref, m_ref, flags_ref)}
    inputs = {label: _field_inputs(side, key, field, convention, weights) for label, (side, _, _) in sides.items()}
    if field == "exit":
        km.zero_fill_differs = _zero_fill({e: _rows(new.exit[e], key) for e in (0, 1)},
                                          {e: _rows(ref.exit[e], key) for e in (0, 1)})
    # BBR027: a basis a side declares transmitting must have a nonzero norm under this field's weights.
    bad = False
    for label, (side, m, flags) in sides.items():
        vectors, w = inputs[label]
        for e in (0, 1):
            if _transmitting(m[e], flags[e]) and not np.sum(w * np.sum(np.abs(vectors[e]) ** 2, axis=1)) > 0.0:
                what = "outward far-field" if field == "far" else "exit-field"
                refusals.append(f"BBR027: {side.label} Ephi={e} at {key_text(key)} is transmitting (m {m[e]:.3g}) but "
                                f"its {what} norm is zero")
                bad = True
    if bad:
        km.skipped = {name: "refused (BBR027)" for name in ("g0", "g1", "h")}
        return km
    for e in (0, 1):
        if _transmitting(m_new[e], flags_new[e]) and _transmitting(m_ref[e], flags_ref[e]):
            km.bases.append(e)
        else:
            reasons = [_not_transmitting(label, m[e], flags[e]) for label, (_, m, flags) in sides.items()
                       if not _transmitting(m[e], flags[e])]
            km.skipped[f"g{e}"] = f"Ephi={e} not compared: " + "; ".join(reasons)
    if not km.bases:
        km.skipped["h"] = "no basis compared"
        return km
    norm = {}
    for label in sides:
        vectors, w = inputs[label]
        norm[label] = _densities(vectors, w, km.bases)
        km.denominators[f"Q_{label}"] = norm[label][2]
    (g_new, h_new, _), (g_ref, h_ref, _) = norm["new"], norm["ref"]
    for e in km.bases:
        top = float(np.max(g_ref[e]))
        km.denominators[f"max_g_ref{e}"] = top
        km.d[f"g{e}"] = float(np.max(np.abs(g_new[e] - g_ref[e]))) / top
    if len(km.bases) < 2:
        km.skipped["h"] = "a basis is not compared"
    elif states != {"new": "lit", "ref": "lit"}:
        km.skipped["h"] = f"cross term not judged: new {states['new']}, reference {states['ref']}"
    else:
        b_ref = float(np.max((g_ref[0] + g_ref[1]) / 2.0))
        km.denominators["B_ref"] = b_ref
        km.d["h"] = float(np.max(np.abs(h_new - h_ref))) / b_ref
    for name, value in list(km.d.items()):
        if not np.isfinite(value):
            refusals.append(f"{field} {name} at {key_text(key)} is not finite ({value!r}): an input or "
                            "classification inconsistency")
            del km.d[name]
            km.skipped[name] = "not finite"
    return km


def compare_sides(new: Side, ref: Side, keys: Sequence[tuple[float, float]]) -> PairingResult:
    """The pointwise comparison of ``new`` with ``ref`` at ``keys``. Alignment problems are refusals, and then no
    metric is formed."""
    keys = [tuple(map(float, k)) for k in keys]
    refusals = [] if keys else [f"{new.label} and {ref.label}: no shared key to compare"]
    refusals += check_alignment(new, ref, keys)
    conventions = {"exit": EXIT_1_1 if new.schema == SCHEMA_1 and ref.schema == SCHEMA_1 else EXIT_RAW, "far": FAR}
    schemas = {"new": new.schema, "ref": ref.schema}
    metrics: list[KeyMetrics] = []
    if not refusals:
        for key in keys:
            try:
                weights = exit_weights(new, ref, key)[1]
            except PairingRefused as exc:
                refusals.append(str(exc))
                continue
            metrics.append(_key_field(new, ref, key, "exit", conventions["exit"], weights, refusals))
            metrics.append(_key_field(new, ref, key, "far", FAR, None, refusals))
    return PairingResult(conventions, schemas, metrics, refusals, _coverage(keys, metrics), _largest(metrics))


def _coverage(keys: Sequence[tuple[float, float]], metrics: Sequence[KeyMetrics]) -> dict[str, int]:
    by_key: dict[tuple[float, float], list[KeyMetrics]] = {}
    for km in metrics:
        by_key.setdefault(km.key, []).append(km)
    return {"shared_keys": len(keys),
            "keys_with_judged_metrics": sum(1 for kms in by_key.values() if any(km.d for km in kms)),
            "metrics_judged": sum(len(km.d) for km in metrics),
            "metrics_skipped": sum(len(km.skipped) for km in metrics),
            "all_dark_keys": sum(1 for kms in by_key.values() if all(not km.bases and "refused (BBR027)"
                                                                    not in km.skipped.values() for km in kms)),
            "transitions": sum(1 for kms in by_key.values() if kms[0].transition),
            "zero_fill_differences": sum(len(km.zero_fill_differs) for km in metrics)}


def _largest(metrics: Sequence[KeyMetrics]) -> dict[str, dict]:
    largest: dict[str, dict] = {}
    for km in metrics:
        for name, value in km.d.items():
            label = f"{km.field}.{name}"
            if label not in largest or value > largest[label]["d"]:
                largest[label] = {"phi_deg": km.key[0], "theta_deg": km.key[1], "d": value}
    return largest


def compare_jobs(new_job: Path, ref_job: Path, inputs: InputCapture,
                 keys: Sequence[tuple[float, float]]) -> PairingResult:
    """``compare_sides`` on two job roots' datasets, read through ``inputs`` keeping only ``keys``; an input that
    cannot be read or validated becomes the result's refusal."""
    try:
        new = read_side("new", job_dataset_dirs(new_job), inputs, keys)
        ref = read_side("reference", job_dataset_dirs(ref_job), inputs, keys)
    except PairingRefused as exc:
        return PairingResult({}, {}, [], [str(exc)], _coverage(keys, []), {})
    return compare_sides(new, ref, keys)


def verdict(result: PairingResult, mode: str | None = None, tol: float | None = None) -> tuple[str, list[str]]:
    """("not judged" | "pass" | "fail", failures) under the policy: ``mode`` and ``tol`` default to the module's
    MODE and TOL at call time. Measure mode judges nothing. Enforce mode (TOL finite and positive) fails a dark/lit
    transition and any judged deviation above TOL, passes when at least one metric was judged, and is "not judged"
    when every metric was legitimately skipped. Refusals are the inputs row's, in either mode. Any other policy is
    refused."""
    mode = MODE if mode is None else mode
    tol = TOL if tol is None and mode == MODE else tol
    if mode == "measure":
        if tol is not None:
            raise PairingRefused(f"policy {POLICY}: measure mode records no tolerance, got {tol!r}")
        return "not judged", []
    if mode != "enforce":
        raise PairingRefused(f"policy {POLICY}: unknown mode {mode!r}, expected 'measure' or 'enforce'")
    if isinstance(tol, bool) or not isinstance(tol, (int, float)) or not np.isfinite(tol) or tol <= 0.0:
        raise PairingRefused(f"policy {POLICY}: enforce mode needs a finite positive tolerance, got {tol!r}")
    failures = []
    seen = set()
    for km in result.keys:
        if km.transition and km.key not in seen:
            seen.add(km.key)
            failures.append(f"{key_text(km.key)}: dark/lit transition (new {km.states['new']}, reference "
                            f"{km.states['ref']})")
        for name, value in km.d.items():
            if value > tol:
                failures.append(f"{km.field} {name} at {key_text(km.key)}: d = {value:.4g} above {tol:g}")
    if failures:
        return "fail", failures
    return ("pass", []) if any(km.d for km in result.keys) else ("not judged", [])


# --- the two rows a checker reports ------------------------------------------------------------------------------------

def _band() -> str:
    return (f"m = |Re S| / (A I) per polarization; dark below {DARK_FLOOR:g} or below {DARK_RATIO:g} of the larger m, "
            f"lit above {LIT_RATIO:g} of it, between otherwise")


def checks(prefix: str, result: PairingResult) -> list[Check]:
    """``<prefix>.pointwise.inputs`` (severity fail: every refusal, an empty key intersection included) and
    ``<prefix>.pointwise`` (severity info while MODE is measure, else fail) for ``result`` under the current policy.
    A metric row whose verdict is "not judged" carries that verdict in its values; it claims no numerical agreement."""
    from bbsim.acceptance_rows import MAX_LISTED, row

    cov = result.coverage
    inputs = row(f"{prefix}.pointwise.inputs",
                 f"both solves validated in their own schema over every declared key; exit points equal to {COORD_M:g} m "
                 f"and far-field angles to {ANGLE_DEG:g} deg at each compared key, in row order; the declared lattice "
                 f"for a {EXIT_RAW} comparison; no transmitting basis with a zero norm (BBR027)",
                 list(result.refusals),
                 f"{cov.get('shared_keys', 0)} shared keys aligned; exit {result.conventions.get('exit')}, far "
                 f"{result.conventions.get('far')}; schemas new {result.schemas.get('new')}, reference "
                 f"{result.schemas.get('ref')}", {"conventions": dict(result.conventions), "schemas": dict(result.schemas)})
    mode, tol = MODE, TOL
    try:
        status, failures = verdict(result)
    except PairingRefused as exc:
        status, failures = "fail", [str(exc)]
    if status == "not judged":
        said = f"not judged ({mode} mode)" if mode == "measure" else "not judged (every metric legitimately skipped)"
    elif status == "pass":
        said = f"pass, every judged d <= {tol:g}"
    else:
        more = f"; and {len(failures) - MAX_LISTED} more" if len(failures) > MAX_LISTED else ""
        said = f"{len(failures)} failure(s): " + "; ".join(failures[:MAX_LISTED]) + more
    largest = ", ".join(f"{name} {v['d']:.3g} at ({v['phi_deg']:g}, {v['theta_deg']:g})"
                        for name, v in sorted(result.largest.items()))
    detail = (f"g0, g1 and h of the two solves point by point ({_band()}): {said}; {cov.get('metrics_judged', 0)} "
              f"metrics at {cov.get('keys_with_judged_metrics', 0)} keys, {cov.get('metrics_skipped', 0)} skipped "
              f"({cov.get('all_dark_keys', 0)} all-dark keys, {cov.get('transitions', 0)} dark/lit transitions)"
              + (f"; largest {largest}" if largest else ""))
    values = {"verdict": status, "mode": mode, "tol": tol, "policy": POLICY, **result.to_dict(), "failures": failures}
    metric = Check(f"{prefix}.pointwise", status != "fail", "info" if mode == "measure" else "fail", detail, values)
    return [inputs, metric]
