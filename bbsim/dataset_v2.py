"""Arithmetic of the 2.0 dataset: the per-key divisor, quadrature weights, the field scale, rho and the
far-field hemisphere power.

Everything here is pure: arrays and numbers in, numbers out, no AEDT. The runner (later tasks) supplies the raw
exit-field table of one key and polarization, the stored exit flux and the declared grid. The only file access is
hashing the three exporter source files for the column-order reference.

Stored exit-field slots (columns 9-14 of the waveguide table, ``schema.EXIT_FIELD_SLOTS``) are interleaved,
whatever the 1.0 headers say: Re Ex, Im Ex, Re Ey, Im Ey, Re Ez, Im Ez (ledger HF-038). Fields are always taken
by position, never by label.

Ledger: HF-040 (|E_t|^2 sampler, quadrature weights), HF-042 (rho from the sampler's density), HF-049 and
OR-036 (cell_90 from the declared grid), OR-034 (the InterceptFactor column), CDX-015 (far-field fraction),
HF-041 and HF-052 (the exporter column-order reference).
"""
from __future__ import annotations

import hashlib
import math
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

from . import schema
from .manifest import sha256_file
from .sampling import SPEED_OF_LIGHT_M_PER_S, VACUUM_PERMEABILITY_H_PER_M

# A flux at or below DARK_REL * (entrance power A*I) is a numerical zero: the key is dark.
DARK_REL = 1e-12
# The 90 degree cell's divisor floors cos(alpha) here (HF-049, BBRsim's uniform-theta emitter).
COS_FLOOR = 0.01
# Accepted range of flux / |Re S| on reachable keys (CDX-009); both ends are inclusive in the check that uses them.
SANITY_BOUNDS = (0.5, 1.01)
# A polarization whose flux is below WEAK_REL of that polarization's largest flux is exempt from the sanity bound.
WEAK_REL = 1e-3
ETA0 = VACUUM_PERMEABILITY_H_PER_M * SPEED_OF_LIGHT_M_PER_S

# The outward guard on far-field rows: sin(Theta) cos(Phi) must exceed this (BBRsim's sampler uses the same).
OUTWARD_GUARD = 1e-6
_THETA_90_TOL_DEG = 1e-9


class DatasetV2Error(ValueError):
    """A 2.0 dataset quantity cannot be formed: the run stops before anything is written."""


def _is_90(theta_deg: float) -> bool:
    return abs(float(theta_deg) - 90.0) <= _THETA_90_TOL_DEG


def cell_90(theta_list: Sequence[float], cos_floor: float = COS_FLOOR) -> float:
    """Divisor of the 90 degree key: 1 / <1 / max(cos alpha, floor)> over the outermost declared band.

    Crossings are uniform in alpha (BBRsim's uniform-theta emitter) from 90 deg - D/2 to 90 deg, where D is the
    gap from 90 deg to the next larger declared theta. The mean has a closed form through the secant integral
    ln((1 + sin a) / cos a); above acos(floor) the integrand is the constant 1/floor. A declared grid with no
    90 deg entry, or none above it, is refused: the band is never guessed.
    """
    if not 0.0 < cos_floor < 1.0:
        raise DatasetV2Error(f"cos_floor must lie in (0, 1), got {cos_floor!r}")
    thetas = sorted({float(t) for t in theta_list})
    if not any(_is_90(t) for t in thetas):
        raise DatasetV2Error(f"cell_90 needs a theta = 90 deg key in the declared grid, got {thetas}")
    above = [t for t in thetas if t > 90.0 + _THETA_90_TOL_DEG]
    if not above:
        raise DatasetV2Error(f"cell_90 needs a declared theta above 90 deg to set the band, got {thetas}")
    delta_deg = above[0] - 90.0
    a0 = math.radians(90.0 - delta_deg / 2.0)
    half_pi = math.pi / 2.0
    a_floor = math.acos(cos_floor)
    if a0 >= a_floor:
        mean_inverse = 1.0 / cos_floor                       # the whole band lies under the floor
    else:
        def sec_integral(a: float) -> float:
            return math.log((1.0 + math.sin(a)) / math.cos(a))
        mean_inverse = (sec_integral(a_floor) - sec_integral(a0) + (half_pi - a_floor) / cos_floor) / (half_pi - a0)
    return 1.0 / mean_inverse


def c_key(phi_deg: float, theta_deg: float, theta_list: Sequence[float]) -> float:
    """The key's intercept factor c_k: |cos theta| off 90 deg, ``cell_90`` at 90 deg.

    ``phi_deg`` is part of the key but does not enter. Keys with theta below 90 deg (unreachable from the
    exit side of the crack) get |cos theta| too; the value is recorded only.
    """
    if _is_90(theta_deg):
        return cell_90(theta_list)
    return abs(math.cos(math.radians(float(theta_deg))))


def quadrature_weights(shape: str, ny: int, nz: int, n_points: int) -> np.ndarray:
    """Dimensionless weights of the retained exit points, in export order (Y outer, Z inner).

    Rectangle: 1 inside, 1/2 on a wall row, 1/4 at a corner, wall rows being the first and last node index
    along each axis (trapezoid rule). Disc: 1 on every retained point, zero-filled rim points included.
    """
    if shape == "rectangle":
        if ny < 2 or nz < 2:
            raise DatasetV2Error(
                f"a rectangle exit lattice needs at least 2 nodes per axis, got {ny} x {nz}: with fewer than 2 "
                "the first and last node coincide, so the wall weight and the step are undefined")
        if n_points != ny * nz:
            raise DatasetV2Error(f"a {ny} x {nz} rectangle lattice has {ny * nz} points, got {n_points}")
        wy = np.ones(ny)
        wy[0] = wy[-1] = 0.5
        wz = np.ones(nz)
        wz[0] = wz[-1] = 0.5
        return np.outer(wy, wz).reshape(-1)
    if shape == "disc":
        if n_points < 1:
            raise DatasetV2Error(f"a disc exit lattice needs at least one point, got {n_points}")
        return np.ones(int(n_points))
    raise DatasetV2Error(f"unknown exit cross-section shape {shape!r}: expected 'rectangle' or 'disc'")


def tangential(values: np.ndarray) -> np.ndarray:
    """Complex (n, 2) = (Ex, Ey) from the six stored slots (Re Ex, Im Ex, Re Ey, Im Ey, Re Ez, Im Ez), by position."""
    v = np.asarray(values, dtype=float)
    if v.ndim != 2 or v.shape[1] != 6:
        raise DatasetV2Error(f"exit-field values must have shape (n, 6), got {v.shape}")
    return np.stack([v[:, 0] + 1j * v[:, 1], v[:, 2] + 1j * v[:, 3]], axis=1)


def _weighted_power(values: np.ndarray, weights: np.ndarray) -> float:
    """Sum_i w_i |E_t,i|^2 (dimensionless sum; callers apply the area and 2 eta0)."""
    et = tangential(values)
    w = np.asarray(weights, dtype=float)
    if w.shape != (et.shape[0],):
        raise DatasetV2Error(f"{w.shape} weights for {et.shape[0]} exit points")
    return float(np.sum(w * np.sum(np.abs(et) ** 2, axis=1)))


def scale_factor(values: np.ndarray, weights: np.ndarray, h_y_m: float, h_z_m: float,
                 flux_w: float, incoming_w: float) -> tuple[float, bool]:
    """(W, dark): the real factor that makes the weighted tangential field carry the stored exit flux.

    Dark (|flux| <= DARK_REL * incoming): W = 1, the scale rule does not apply. A flux below minus that floor
    stops the run. Otherwise S = sum_i w_i |E_t,i|^2 h_y h_z / (2 eta0) must be finite and positive, and
    W = sqrt(flux / S), so that the scaled field reproduces the flux exactly.
    """
    flux = float(flux_w)
    incoming = float(incoming_w)
    if not math.isfinite(flux):
        raise DatasetV2Error(f"exit flux is not finite: {flux_w!r}")
    if not (math.isfinite(incoming) and incoming > 0.0):
        raise DatasetV2Error(f"entrance power must be finite and positive, got {incoming_w!r}")
    floor = DARK_REL * incoming
    if abs(flux) <= floor:
        return 1.0, True
    if flux < -floor:
        raise DatasetV2Error(
            f"negative exit flux {flux!r} W lies below the dark floor -{floor!r} W (|flux| <= {DARK_REL} A I "
            "is a numerical zero): the key is unphysical")
    s = _weighted_power(values, weights) * float(h_y_m) * float(h_z_m) / (2.0 * ETA0)
    if not (math.isfinite(s) and s > 0.0):
        raise DatasetV2Error(
            f"the weighted tangential field integral is {s!r} for a non-dark exit flux {flux!r} W: "
            "it must be finite and positive to scale the field")
    return math.sqrt(flux / s), False


def rho(values0: np.ndarray, values1: np.ndarray, weights: np.ndarray, dark0: bool, dark1: bool) -> float:
    """Re rho = Re sum_i w_i (E_t0,i . conj E_t1,i) / sqrt(sum w |E_t0|^2 * sum w |E_t1|^2); 0 if either is dark."""
    if dark0 or dark1:
        return 0.0
    e0 = tangential(values0)
    e1 = tangential(values1)
    w = np.asarray(weights, dtype=float)
    if not (e0.shape == e1.shape and w.shape == (e0.shape[0],)):
        raise DatasetV2Error(f"rho needs matching tables, got {e0.shape}, {e1.shape} and {w.shape} weights")
    cross = float(np.real(np.sum(w * np.sum(e0 * np.conj(e1), axis=1))))
    p0 = _weighted_power(values0, w)
    p1 = _weighted_power(values1, w)
    if not (math.isfinite(p0) and math.isfinite(p1) and p0 > 0.0 and p1 > 0.0 and math.isfinite(cross)):
        raise DatasetV2Error(
            f"rho is undefined for a non-dark pair: weighted tangential powers {p0!r} and {p1!r}, cross {cross!r}")
    return cross / (math.sqrt(p0) * math.sqrt(p1))


def hemisphere_power_w(ff_rows: pd.DataFrame, dtheta_deg: float, dphi_deg: float) -> float:
    """Far-field power over the outward hemisphere the CSV spans, in W.

    Rows with sin(Theta) cos(Phi) > OUTWARD_GUARD are summed with weight sin(Theta) dTheta dPhi (radians) and
    integrand (|rE_theta|^2 + |rE_phi|^2) / (2 eta0).
    """
    needed = ["Theta", "Phi", "rEtheta_real", "rEtheta_imag", "rEphi_real", "rEphi_imag"]
    missing = [c for c in needed if c not in ff_rows.columns]
    if missing:
        raise DatasetV2Error(f"far-field table lacks columns {missing}")
    theta = np.radians(ff_rows["Theta"].to_numpy(dtype=float))
    phi = np.radians(ff_rows["Phi"].to_numpy(dtype=float))
    sin_theta = np.sin(theta)
    outward = sin_theta * np.cos(phi) > OUTWARD_GUARD
    intensity = (ff_rows["rEtheta_real"].to_numpy(dtype=float) ** 2 + ff_rows["rEtheta_imag"].to_numpy(dtype=float) ** 2
                 + ff_rows["rEphi_real"].to_numpy(dtype=float) ** 2 + ff_rows["rEphi_imag"].to_numpy(dtype=float) ** 2)
    d_omega_weight = sin_theta * math.radians(float(dtheta_deg)) * math.radians(float(dphi_deg))
    return float(np.sum((d_omega_weight * intensity)[outward]) / (2.0 * ETA0))


def sanity_ratio(flux_w: float, abs_w: float, incoming_w: float) -> float | None:
    """flux / |Re S| over the exit face (CDX-009): the ratio, None for a dark key with |Re S| = 0.

    A positive flux with |Re S| = 0 cannot be a real field: it raises.
    """
    flux = float(flux_w)
    absolute = float(abs_w)
    incoming = float(incoming_w)
    if not (math.isfinite(flux) and math.isfinite(absolute) and math.isfinite(incoming)):
        raise DatasetV2Error(f"sanity ratio needs finite powers, got flux {flux_w!r}, |Re S| {abs_w!r}")
    if absolute == 0.0:
        if abs(flux) <= DARK_REL * incoming:
            return None
        raise DatasetV2Error(f"exit flux {flux!r} W is not dark but the exit |Re S| integral is zero")
    return flux / absolute


# --- The exporter column-order reference (HF-041, HF-052) -------------------------------------------------------
#
# The six stored exit-field values are interleaved (Re Ex, Im Ex, Re Ey, Im Ey, Re Ez, Im Ez) whatever the 1.0
# headers say. The order is tested once, on a known reference: crack1 at 500 GHz, whose normal-incidence field at
# (phi, theta) = (0, 180) for Ephi=0 and (90, 180) for Ephi=1 is almost entirely Ex. Every other job embeds that
# reference and is allowed to run only while the exporter source files and the AEDT and PyAEDT versions are the
# ones the reference was measured with. Field composition of other geometries or frequencies is never tested.

REFERENCE_DATASET = "InfParallelPlate_crack1Rohan"
REFERENCE_GHZ = 500.0
# (ephi, phi_deg, theta_deg)
REFERENCE_CASES = ((0, 0.0, 180.0), (1, 90.0, 180.0))
# Source files that decide how the fields are exported and parsed; repository-relative, hashed as they are.
EXPORTER_FILES = ("bbsim/extract.py", "bbsim/readers.py", "bbsim/schema.py")
ORDER_RULE = "reference_exporter_interleaved_v1"
ORDER_METHODS = ("reference", "shared_exporter")

# A case passes when (P1 + P2) / P >= ENERGY_FRACTION_12_MIN; it discriminates when P2 / P >=
# DISCRIMINATING_FRACTION_2_MIN, and then P2 >= SLOT2_OVER_SLOT4_MIN * P4 as well.
ENERGY_FRACTION_12_MIN = 0.99
DISCRIMINATING_FRACTION_2_MIN = 0.01
SLOT2_OVER_SLOT4_MIN = 1e3

_REFERENCE_KEYS = ("dataset_id", "frequency_ghz", "source_job", "git_commit", "aedt_version", "pyaedt_version",
                   "exporter_sha256", "raw_table_sha256", "cases", "passed")
_CASE_KEYS = ("ephi", "phi_deg", "theta_deg", "energy_fraction_12", "energy_fraction_2", "slot2_over_slot4",
              "discriminating", "passed")


def _judge(fraction_12: float, fraction_2: float, ratio: float | None) -> tuple[bool, bool]:
    """(discriminating, passed) from a case's three metrics; the single place the thresholds are applied.

    ``ratio`` is P2 / P4, None when P4 = 0, in which case P2 >= 1e3 P4 holds trivially.
    """
    discriminating = fraction_2 >= DISCRIMINATING_FRACTION_2_MIN
    ratio_ok = (not discriminating) or ratio is None or ratio >= SLOT2_OVER_SLOT4_MIN
    return discriminating, bool(fraction_12 >= ENERGY_FRACTION_12_MIN and ratio_ok)


def order_case(values: np.ndarray, ephi: int, phi: float, theta: float) -> dict[str, Any]:
    """Evaluate one reference case on the (n, 6) stored slots of its key, by position.

    With P_j = sum_i value_ij^2 and P = sum_j P_j: energy_fraction_12 = (P1 + P2) / P, energy_fraction_2 =
    P2 / P, slot2_over_slot4 = P2 / P4 (None when P4 = 0). A table with no energy fails; a table with a
    non-finite value is refused.
    """
    v = np.asarray(values, dtype=float)
    if v.ndim != 2 or v.shape[1] != 6:
        raise DatasetV2Error(f"column-order case needs exit-field values of shape (n, 6), got {v.shape}")
    if v.shape[0] < 1:
        raise DatasetV2Error("column-order case has no rows")
    if not np.all(np.isfinite(v)):
        raise DatasetV2Error("column-order case has non-finite exit-field values")
    p = [float(x) for x in np.sum(v * v, axis=0)]
    total = math.fsum(p)
    if not math.isfinite(total):
        raise DatasetV2Error("column-order case energy is not finite")
    if total > 0.0:
        fraction_12 = (p[0] + p[1]) / total
        fraction_2 = p[1] / total
    else:
        fraction_12 = fraction_2 = 0.0
    ratio = p[1] / p[3] if p[3] > 0.0 else None
    discriminating, passed = _judge(fraction_12, fraction_2, ratio)
    return {
        "ephi": int(ephi), "phi_deg": float(phi), "theta_deg": float(theta),
        "energy_fraction_12": float(fraction_12), "energy_fraction_2": float(fraction_2),
        "slot2_over_slot4": None if ratio is None else float(ratio),
        "discriminating": discriminating, "passed": passed,
    }


def exporter_hashes(repo_root: Path) -> dict[str, str]:
    """SHA-256 of each file in ``EXPORTER_FILES`` under ``repo_root``, keyed by its repository-relative path."""
    out: dict[str, str] = {}
    for name in EXPORTER_FILES:
        path = Path(repo_root) / name
        if not path.is_file():
            raise DatasetV2Error(f"exporter file {name} not found under {repo_root}: the reference cannot name "
                                 "the implementation that produced it")
        out[name] = sha256_file(path)
    return out


def order_check_block(method: str, reference: dict) -> dict[str, Any]:
    """The sidecar's ``exit_field.column_order_check``: the rule, how this job obtained the evidence, the evidence."""
    if method not in ORDER_METHODS:
        raise DatasetV2Error(f"column-order method must be one of {ORDER_METHODS}, got {method!r}")
    return {"rule": ORDER_RULE, "method": method, "reference": reference}


def _require_version(name: str, value: object) -> str:
    if not isinstance(value, str) or not value.strip():
        raise DatasetV2Error(f"column-order reference needs a non-empty {name} string (version evidence), "
                             f"got {value!r}")
    return value


def _reference_case_rows(table: pd.DataFrame, ephi: int, phi: float, theta: float) -> np.ndarray:
    """The six stored slots (columns 9-14, by position) of the rows at one incidence key of a raw table."""
    if list(table.columns) != list(schema.WAVEGUIDE_COLUMNS):
        raise DatasetV2Error(
            f"the raw reference table for Ephi={ephi} must carry the 1.0 waveguide columns "
            f"{schema.WAVEGUIDE_COLUMNS}, got {list(table.columns)}: the evidence is hashed as exported, "
            "before any schema-2 relabelling")
    try:
        phis = table["IWavePhi"].to_numpy(dtype=float)
        thetas = table["IWaveTheta"].to_numpy(dtype=float)
        mask = np.isclose(phis, phi, rtol=0.0, atol=1e-9) & np.isclose(thetas, theta, rtol=0.0, atol=1e-9)
        values = table.iloc[:, schema.EXIT_FIELD_SLOTS].to_numpy(dtype=float)[mask]
    except (TypeError, ValueError) as exc:
        raise DatasetV2Error(f"the raw reference table for Ephi={ephi} has a non-numeric key or field: {exc}") from exc
    if values.shape[0] == 0:
        raise DatasetV2Error(f"the raw reference table for Ephi={ephi} has no rows at (phi, theta) = "
                             f"({phi:g}, {theta:g}) deg")
    return values


def build_reference(raw_tables: dict[int, pd.DataFrame], dataset_id: str, frequency_ghz: float,
                    source_job: str | None, git_commit: str | None, aedt_version: str, pyaedt_version: str,
                    repo_root: Path) -> dict[str, Any]:
    """The column-order reference object, from the fresh unscaled crack1 500 GHz tables of both polarizations.

    Only ``REFERENCE_DATASET`` at ``REFERENCE_GHZ`` is eligible. Both cases must pass and at least one must
    discriminate; otherwise the run stops with ``DatasetV2Error``. ``raw_tables`` maps Ephi (0, 1) to the raw
    1.0-labelled waveguide table; each is hashed from ``to_csv(index=False)`` as given.
    """
    try:
        frequency_ok = float(frequency_ghz) == REFERENCE_GHZ
    except (TypeError, ValueError):
        frequency_ok = False
    if dataset_id != REFERENCE_DATASET or not frequency_ok:
        raise DatasetV2Error(
            f"the column-order reference must be measured on dataset {REFERENCE_DATASET} at {REFERENCE_GHZ:g} GHz, "
            f"got {dataset_id!r} at {frequency_ghz!r} GHz")
    _require_version("aedt_version", aedt_version)
    _require_version("pyaedt_version", pyaedt_version)
    if set(raw_tables) != {0, 1}:
        raise DatasetV2Error(f"the column-order reference needs the raw tables of Ephi 0 and 1, got {sorted(raw_tables)}")

    cases = [order_case(_reference_case_rows(raw_tables[ephi], ephi, phi, theta), ephi, phi, theta)
             for ephi, phi, theta in REFERENCE_CASES]
    for case in cases:
        if not case["passed"]:
            raise DatasetV2Error(
                f"column-order reference case Ephi={case['ephi']} at (phi, theta) = ({case['phi_deg']:g}, "
                f"{case['theta_deg']:g}) deg failed: (P1+P2)/P = {case['energy_fraction_12']!r} (need >= "
                f"{ENERGY_FRACTION_12_MIN}), P2/P = {case['energy_fraction_2']!r}, P2/P4 = "
                f"{case['slot2_over_slot4']!r}; the stored slots are not Re Ex, Im Ex, Re Ey, ...")
    if not any(case["discriminating"] for case in cases):
        raise DatasetV2Error(
            "no column-order reference case discriminates the interleaved order from the label order: both have "
            f"P2/P below {DISCRIMINATING_FRACTION_2_MIN} ({[c['energy_fraction_2'] for c in cases]}); a reference "
            "whose Ex is nearly real cannot tell Re Ex, Im Ex from Re Ex, Re Ey")

    reference = {
        "dataset_id": dataset_id,
        "frequency_ghz": REFERENCE_GHZ,
        "source_job": source_job,
        "git_commit": git_commit,
        "aedt_version": aedt_version,
        "pyaedt_version": pyaedt_version,
        "exporter_sha256": exporter_hashes(repo_root),
        "raw_table_sha256": {str(ephi): hashlib.sha256(raw_tables[ephi].to_csv(index=False).encode()).hexdigest()
                             for ephi in (0, 1)},
        "cases": cases,
        "passed": True,
    }
    check_reference_compatible(reference, aedt_version, pyaedt_version, repo_root)       # self-consistency
    return reference


def _is_number(x: object) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def _is_sha256(x: object) -> bool:
    return isinstance(x, str) and len(x) == 64 and all(c in "0123456789abcdef" for c in x)


def _validate_reference_case(case: object, index: int) -> tuple[int, float, float]:
    """Check one case's shape, metrics and flags; return its (ephi, phi, theta) key."""
    where = f"column-order reference case {index}"
    if not isinstance(case, Mapping):
        raise DatasetV2Error(f"{where} must be an object, got {type(case).__name__}")
    if set(case) != set(_CASE_KEYS):
        raise DatasetV2Error(f"{where} must hold exactly {list(_CASE_KEYS)}; missing "
                             f"{[k for k in _CASE_KEYS if k not in case]}, unexpected {sorted(set(case) - set(_CASE_KEYS))}")
    ephi = case["ephi"]
    if type(ephi) is not int or not _is_number(case["phi_deg"]) or not _is_number(case["theta_deg"]):
        raise DatasetV2Error(f"{where} needs an integer ephi and numeric phi_deg and theta_deg, got "
                             f"({case['ephi']!r}, {case['phi_deg']!r}, {case['theta_deg']!r})")
    key = (ephi, float(case["phi_deg"]), float(case["theta_deg"]))
    f12, f2, ratio = case["energy_fraction_12"], case["energy_fraction_2"], case["slot2_over_slot4"]
    for name, value in (("energy_fraction_12", f12), ("energy_fraction_2", f2)):
        if not _is_number(value) or not 0.0 <= value <= 1.0:
            raise DatasetV2Error(f"{where}: {name} must be a finite number in [0, 1], got {value!r}")
    if f2 > f12:
        raise DatasetV2Error(f"{where}: P2/P = {f2!r} exceeds (P1+P2)/P = {f12!r}, which is impossible")
    if ratio is not None and (not _is_number(ratio) or ratio < 0.0):
        raise DatasetV2Error(f"{where}: slot2_over_slot4 must be null (P4 = 0) or a finite number >= 0, got {ratio!r}")
    for name in ("discriminating", "passed"):
        if not isinstance(case[name], bool):
            raise DatasetV2Error(f"{where}: {name} must be a boolean, got {case[name]!r}")
    discriminating, passed = _judge(float(f12), float(f2), None if ratio is None else float(ratio))
    if case["discriminating"] != discriminating:
        raise DatasetV2Error(f"{where}: declared discriminating={case['discriminating']} but its metrics give "
                             f"{discriminating} (P2/P = {f2!r}, threshold {DISCRIMINATING_FRACTION_2_MIN})")
    if case["passed"] != passed:
        raise DatasetV2Error(f"{where}: declared passed={case['passed']} but its metrics give {passed} "
                             f"((P1+P2)/P = {f12!r}, P2/P4 = {ratio!r})")
    return key


def check_reference_compatible(reference: dict, aedt_version: str, pyaedt_version: str, repo_root: Path) -> None:
    """Refuse a reference that is incomplete, inconsistent, failed, or from other code or versions.

    Validates the embedded evidence only (the declared metrics, recomputed flags and hashes); re-evaluating the
    raw tables belongs to ``build_reference`` and to acceptance of the reference job, not to every
    shared-exporter dataset. A top-level ``passed: true`` is never trusted: both cases must pass and at least one
    must discriminate. Then every exporter hash and both versions must equal the current ones; the first
    mismatch is named.
    """
    if not isinstance(reference, Mapping):
        raise DatasetV2Error(f"the column-order reference must be an object, got {type(reference).__name__}")
    _require_version("current aedt_version", aedt_version)
    _require_version("current pyaedt_version", pyaedt_version)
    if set(reference) != set(_REFERENCE_KEYS):
        raise DatasetV2Error(
            f"the column-order reference must hold exactly {list(_REFERENCE_KEYS)}; missing "
            f"{[k for k in _REFERENCE_KEYS if k not in reference]}, unexpected "
            f"{sorted(set(reference) - set(_REFERENCE_KEYS))}")
    if reference["dataset_id"] != REFERENCE_DATASET:
        raise DatasetV2Error(f"the column-order reference must be dataset {REFERENCE_DATASET}, "
                             f"got {reference['dataset_id']!r}")
    if not _is_number(reference["frequency_ghz"]) or float(reference["frequency_ghz"]) != REFERENCE_GHZ:
        raise DatasetV2Error(f"the column-order reference must be at {REFERENCE_GHZ:g} GHz, "
                             f"got {reference['frequency_ghz']!r}")
    for name in ("source_job", "git_commit"):
        if reference[name] is not None and not isinstance(reference[name], str):
            raise DatasetV2Error(f"column-order reference {name} must be a string or null, got {reference[name]!r}")
    for name in ("aedt_version", "pyaedt_version"):
        _require_version(f"reference {name}", reference[name])

    exporter = reference["exporter_sha256"]
    if not isinstance(exporter, Mapping) or set(exporter) != set(EXPORTER_FILES) or not all(
            _is_sha256(exporter[k]) for k in EXPORTER_FILES):
        raise DatasetV2Error(f"column-order reference exporter_sha256 must map exactly {list(EXPORTER_FILES)} to "
                             f"lowercase SHA-256 hex digests, got {exporter!r}")
    raw = reference["raw_table_sha256"]
    if not isinstance(raw, Mapping) or set(raw) != {"0", "1"} or not all(_is_sha256(raw[k]) for k in ("0", "1")):
        raise DatasetV2Error(f"column-order reference raw_table_sha256 must map exactly '0' and '1' to lowercase "
                             f"SHA-256 hex digests, got {raw!r}")

    cases = reference["cases"]
    if not isinstance(cases, (list, tuple)) or len(cases) != len(REFERENCE_CASES):
        raise DatasetV2Error(f"the column-order reference must hold exactly {len(REFERENCE_CASES)} cases "
                             f"{list(REFERENCE_CASES)}, got {cases!r}")
    keys = [_validate_reference_case(case, i) for i, case in enumerate(cases)]
    if sorted(keys) != sorted(REFERENCE_CASES):
        raise DatasetV2Error(f"the column-order reference cases must be exactly the two distinct keys "
                             f"{list(REFERENCE_CASES)} (ephi, phi_deg, theta_deg), got {keys}")
    all_passed = all(case["passed"] for case in cases)
    any_discriminating = any(case["discriminating"] for case in cases)
    overall = all_passed and any_discriminating
    if not isinstance(reference["passed"], bool) or reference["passed"] != overall:
        raise DatasetV2Error(f"column-order reference declared passed={reference['passed']!r} but its cases give "
                             f"{overall} (all pass: {all_passed}, one discriminates: {any_discriminating})")
    if not overall:
        raise DatasetV2Error("the column-order reference did not pass: both cases must pass and at least one "
                             f"must discriminate (all pass: {all_passed}, one discriminates: {any_discriminating})")

    current = exporter_hashes(repo_root)
    for name in EXPORTER_FILES:
        if exporter[name] != current[name]:
            raise DatasetV2Error(
                f"the column-order reference was measured with a different {name} (reference sha256 "
                f"{exporter[name]}, current {current[name]}): export the reference again")
    for name, now in (("aedt_version", aedt_version), ("pyaedt_version", pyaedt_version)):
        if reference[name] != now:
            raise DatasetV2Error(
                f"the column-order reference was measured with {name} {reference[name]!r}, the current one is "
                f"{now!r}: export the reference again")
