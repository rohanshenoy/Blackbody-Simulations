"""Arithmetic of the 2.0 dataset: the per-key divisor, quadrature weights, the field scale, rho and the
far-field hemisphere power.

Everything here is pure: arrays and numbers in, numbers out, no AEDT, no file access. The runner (later tasks)
supplies the raw exit-field table of one key and polarization, the stored exit flux and the declared grid.

Stored exit-field slots (columns 9-14 of the waveguide table, ``schema.EXIT_FIELD_SLOTS``) are interleaved,
whatever the 1.0 headers say: Re Ex, Im Ex, Re Ey, Im Ey, Re Ez, Im Ez (ledger HF-038). Fields are always taken
by position, never by label.

Ledger: HF-040 (|E_t|^2 sampler, quadrature weights), HF-042 (rho from the sampler's density), HF-049 and
OR-036 (cell_90 from the declared grid), OR-034 (the InterceptFactor column), CDX-015 (far-field fraction).
"""
from __future__ import annotations

import math
from typing import Sequence

import numpy as np
import pandas as pd

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
