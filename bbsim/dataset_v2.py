"""Arithmetic of the 2.0 dataset: the per-key divisor, quadrature weights, the field scale, rho and the
far-field hemisphere power; ``build`` assembles them into the 2.0 tables and sidecar blocks. It also holds the readers
the 2.0 producer and consumers share (``key_text``, ``key_index``, ``row_keys``, ``declared_grid``, ...), so the runner,
``bbsim.compare`` and ``bbsim.reexport_check`` read a key, a row's key and a declared lattice the same way.

Everything here is pure: arrays and numbers in, numbers out, no AEDT. The runner supplies one ``Extracted`` per
polarization (the raw 1.0-labelled tables and the per-key exit flux and |Re S|), the exit lattice and the declared
grid. The only file access is hashing the three exporter source files for the column-order reference.

Stored exit-field slots (columns 9-14 of the waveguide table, ``schema.EXIT_FIELD_SLOTS``) are interleaved,
whatever the 1.0 headers say: Re Ex, Im Ex, Re Ey, Im Ey, Re Ez, Im Ez (ledger HF-038). Fields are always taken
by position, never by label.

Ledger: HF-040 (|E_t|^2 sampler, quadrature weights), HF-042 (rho from the sampler's density), HF-049 and
OR-036 (cell_90 from the declared grid), OR-034 (the InterceptFactor column), CDX-015 (far-field fraction),
HF-041 and HF-052 (the exporter column-order reference), CDX-009 (flux numerator, sanity bound), OR-038 (dark
records).
"""
from __future__ import annotations

import dataclasses
import hashlib
import math
from collections.abc import Mapping
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

import numpy as np
import pandas as pd

from . import schema
from .manifest import sha256_file
from .sampling import SPEED_OF_LIGHT_M_PER_S, VACUUM_PERMEABILITY_H_PER_M

if TYPE_CHECKING:
    from .extract import ExitGrid

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
# A retained exit point must lie within this fraction of the lattice step of its node (BBRsim's C4 node tolerance).
NODE_TOL = 1e-6
_THETA_90_TOL_DEG = 1e-9


class DatasetV2Error(ValueError):
    """A 2.0 dataset quantity cannot be formed: the run stops before anything is written."""


def _is_90(theta_deg: float) -> bool:
    return abs(float(theta_deg) - 90.0) <= _THETA_90_TOL_DEG


def cell_90_band_deg(theta_list: Sequence[float]) -> tuple[float, float]:
    """(90 - D/2, 90): the 90 deg cell's band, D the gap from 90 deg to the next larger declared theta.

    A declared grid with no 90 deg entry, or none above it, is refused: the band is never guessed.
    """
    thetas = sorted({float(t) for t in theta_list})
    if not any(_is_90(t) for t in thetas):
        raise DatasetV2Error(f"cell_90 needs a theta = 90 deg key in the declared grid, got {thetas}")
    above = [t for t in thetas if t > 90.0 + _THETA_90_TOL_DEG]
    if not above:
        raise DatasetV2Error(f"cell_90 needs a declared theta above 90 deg to set the band, got {thetas}")
    delta_deg = above[0] - 90.0
    return 90.0 - delta_deg / 2.0, 90.0


def cell_90(theta_list: Sequence[float], cos_floor: float = COS_FLOOR) -> float:
    """Divisor of the 90 degree key: 1 / <1 / max(cos alpha, floor)> over the outermost declared band.

    Crossings are uniform in alpha (BBRsim's uniform-theta emitter) from 90 deg - D/2 to 90 deg, where D is the
    gap from 90 deg to the next larger declared theta. The mean has a closed form through the secant integral
    ln((1 + sin a) / cos a); above acos(floor) the integrand is the constant 1/floor. A declared grid with no
    90 deg entry, or none above it, is refused: the band is never guessed.
    """
    if not 0.0 < cos_floor < 1.0:
        raise DatasetV2Error(f"cos_floor must lie in (0, 1), got {cos_floor!r}")
    lo_deg, _ = cell_90_band_deg(theta_list)
    a0 = math.radians(lo_deg)
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


def _axis_weights(n: int) -> np.ndarray:
    """Trapezoid weights along one axis of n nodes: 1/2 on the first and last node, 1 between."""
    w = np.ones(n)
    w[0] = w[-1] = 0.5
    return w


def _rectangle_axes(ny: int, nz: int) -> tuple[np.ndarray, np.ndarray]:
    if ny < 2 or nz < 2:
        raise DatasetV2Error(
            f"a rectangle exit lattice needs at least 2 nodes per axis, got {ny} x {nz}: with fewer than 2 "
            "the first and last node coincide, so the wall weight and the step are undefined")
    return _axis_weights(ny), _axis_weights(nz)


def quadrature_weights(shape: str, ny: int, nz: int, n_points: int) -> np.ndarray:
    """Dimensionless weights of the retained exit points, in export order (Y outer, Z inner).

    Rectangle: 1 inside, 1/2 on a wall row, 1/4 at a corner, wall rows being the first and last node index
    along each axis (trapezoid rule). Disc: 1 on every retained point, zero-filled rim points included.

    This is the weight of the node at each row *assuming* export order; ``node_weights`` assigns the same weights
    from the points' own coordinates, as BBRsim does, and is what ``build`` uses for a rectangle.
    """
    if shape == "rectangle":
        wy, wz = _rectangle_axes(ny, nz)
        if n_points != ny * nz:
            raise DatasetV2Error(f"a {ny} x {nz} rectangle lattice has {ny * nz} points, got {n_points}")
        return np.outer(wy, wz).reshape(-1)
    if shape == "disc":
        if n_points < 1:
            raise DatasetV2Error(f"a disc exit lattice needs at least one point, got {n_points}")
        return np.ones(int(n_points))
    raise DatasetV2Error(f"unknown exit cross-section shape {shape!r}: expected 'rectangle' or 'disc'")


def _lattice_nodes(y_m: np.ndarray, z_m: np.ndarray, grid: ExitGrid) -> tuple[np.ndarray, np.ndarray, int, int]:
    """Node indices (iy, iz) of each exit point of a rectangle lattice, plus the lattice's (ny, nz).

    Node = round((coordinate - min) / step) per axis, with min = the grid's y_e / z_e lower bound and step its step
    along that axis, both in metres: BBRsim's AxisNode. Refuses (``DatasetV2Error``) a count other than ny * nz, a
    non-finite coordinate, a point off its node by more than ``NODE_TOL`` of the step, a node outside
    [0, count - 1], and two points at one node (so a passing call covers every node exactly once).
    """
    ny, nz = int(grid.counts[1]), int(grid.counts[2])
    _rectangle_axes(ny, nz)                                     # at least 2 nodes per axis
    y = np.asarray(y_m, dtype=float)
    z = np.asarray(z_m, dtype=float)
    if y.ndim != 1 or y.shape != z.shape:
        raise DatasetV2Error(f"exit-point Y and Z must be equal-length 1-D arrays, got shapes {y.shape} and {z.shape}")
    if y.size != ny * nz:
        raise DatasetV2Error(f"a {ny} x {nz} rectangle lattice has {ny * nz} points, got {y.size}")
    if not (np.all(np.isfinite(y)) and np.all(np.isfinite(z))):
        raise DatasetV2Error("exit-point Y and Z must be finite")
    indices = []
    for axis, coord, count, lo_mm, step_mm in (("y_e", y, ny, grid.lo_mm[1], grid.step_mm[1]),
                                               ("z_e", z, nz, grid.lo_mm[2], grid.step_mm[2])):
        step = _positive_finite(f"the exit lattice step along {axis}", step_mm) * 1e-3
        steps_from_min = (coord - lo_mm * 1e-3) / step
        node = np.rint(steps_from_min)
        outside = (node < 0) | (node > count - 1)
        if outside.any():
            i = int(np.flatnonzero(outside)[0])
            raise DatasetV2Error(f"exit point {i} lies outside the {count}-node {axis} lattice: node {node[i]:g} of "
                                 f"0..{count - 1} (coordinate {coord[i]!r} m)")
        off = np.abs(steps_from_min - node)
        if (off > NODE_TOL).any():
            i = int(np.argmax(off))
            raise DatasetV2Error(f"exit point {i} is off its {axis} node by {off[i]:g} of the step "
                                 f"(limit {NODE_TOL:g}): coordinate {coord[i]!r} m, nearest node {node[i]:g}")
        indices.append(node.astype(int))
    iy, iz = indices
    flat = iy * nz + iz
    seen = np.bincount(flat, minlength=ny * nz)
    if (seen > 1).any():
        node = int(np.flatnonzero(seen > 1)[0])
        raise DatasetV2Error(f"two exit points sit at the same node (iy, iz) = ({node // nz}, {node % nz}): "
                             "a rectangle lattice holds each node exactly once")
    return iy, iz, ny, nz


def _node_weights(iy: np.ndarray, iz: np.ndarray, ny: int, nz: int) -> np.ndarray:
    wy, wz = _rectangle_axes(ny, nz)
    return wy[iy] * wz[iz]


def node_weights(y_m: np.ndarray, z_m: np.ndarray, grid: ExitGrid) -> np.ndarray:
    """Trapezoid weights of a rectangle's retained exit points from each point's own grid node.

    The weight of the point at node (iy, iz) is (1/2 if iy is 0 or ny - 1 else 1) times (1/2 if iz is 0 or nz - 1
    else 1), the nodes found as BBRsim finds them (AxisNode): i = round((coordinate - min) / step) on Y and on Z,
    in metres, against the grid's y_e / z_e bounds and steps. The result follows the input row order whatever
    that is; in export order (Y outer, Z inner) it equals ``quadrature_weights("rectangle", ny, nz, n)``.

    The point set must be the whole ny x nz lattice, every point within ``NODE_TOL`` of the step of its node.
    Anything else raises ``DatasetV2Error``. Rectangles only: a disc keeps weight 1 on every retained point.
    """
    return _node_weights(*_lattice_nodes(y_m, z_m, grid))


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
    """A version string that is evidence: non-empty and not ``session.aedt_versions``' fallback ``"unknown (...)"``,
    which two runs that both failed to read their version would share (Ruling 5)."""
    if not isinstance(value, str) or not value.strip():
        raise DatasetV2Error(f"column-order reference needs a non-empty {name} string (version evidence), "
                             f"got {value!r}")
    if value.strip().lower().startswith("unknown"):
        raise DatasetV2Error(f"column-order reference needs a known {name} (version evidence), got {value!r}: "
                             "the version could not be read")
    return value


def check_versions(aedt_version: object, pyaedt_version: object) -> None:
    """Refuse current AEDT / PyAEDT versions that could not serve as reference evidence, as ``build_reference`` and
    ``check_reference_compatible`` would; a runner calls it before a solve that would otherwise be wasted."""
    _require_version("aedt_version", aedt_version)
    _require_version("pyaedt_version", pyaedt_version)


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


# --- The 2.0 tables and sidecar blocks from one job's raw extraction (CDX-009, HF-040, HF-042, OR-034, OR-038) ----
#
# Strings are verbatim from the field list (docs/superpowers/specs/2026-10-13-schema-2-sidecar-fields.md, section 2).

TRANSMITTANCE_DEFINITION = "OutgoingPower / (IngoingPower * InterceptFactor)"
TRANSMITTANCE_NUMERATOR = "exit_flux_re_s_dot_n"
TRANSMITTANCE_INCOMING = "area_times_intensity_undivided"
TRANSMITTANCE_ABOVE_ONE = "per_photon_cap"
CELL_90_EMITTER = "uniform_theta"
QUADRATURE_BY_SHAPE = {"rectangle": "trapezoid", "disc": "uniform_disc"}
TANGENTIAL_COMPONENTS = ("Ex", "Ey")
SCALE_RULE = "sum_i w_i abs(E_t,i)^2 h_y h_z / (2 eta0) = outgoing_power_w"
SAMPLER = "abs_e_t_squared_times_quadrature"
CROSS_TERM_DENSITY = "e_t_quadrature"
CROSS_TERM_DEFINITION = ("Re rho = sum_i w_i Re(E_t0,i . conj E_t1,i) / sqrt(sum_i w_i abs(E_t0,i)^2 * "
                         "sum_i w_i abs(E_t1,i)^2)")
CROSS_TERM_ABOVE_LIMIT = "warn_once"
FAR_FIELD_NORMALIZATION = "pattern_only"
# Declared keys and table rows match within this, in degrees (field list: per-key lists matched to 1e-9).
KEY_TOL_DEG = 1e-9
_FAR_FIELD_NUMERIC = ("IWavePhi", "IWaveTheta", "Phi", "Theta",
                      "rEphi_real", "rEphi_imag", "rEtheta_real", "rEtheta_imag")


@dataclasses.dataclass
class Extracted:
    """One polarization's raw extraction: the waveguide table with the 1.0 labels (values in stored, interleaved
    order), the far-field table, and per (phi, theta) key the exit flux Re S.n and |Re S| over the exit face, in W."""

    waveguide: pd.DataFrame
    far_field: pd.DataFrame
    flux_w: dict[tuple[float, float], float]
    abs_w: dict[tuple[float, float], float]


@dataclasses.dataclass
class DatasetV2Result:
    """What ``build`` returns.

    ``waveguide``: per Ephi, the 2.0 table (``schema.WAVEGUIDE_COLUMNS_V2``). ``far_field``: per Ephi, the input
    table unchanged. ``blocks``: the sidecar blocks ``transmittance``, ``exit_field`` (``column_order_check``,
    ``tangential_components``, ``quadrature``, ``scale``, ``sampler``), ``cross_term``, ``far_field``
    (``normalization``, ``ff_over_flux``). ``raw_powers_per_key``: the manifest's ``exit_field.raw_powers_per_key``
    (signed flux before the dark clamp, and |Re S|). ``flux_over_abs_re_s``: per key, flux / |Re S| per polarization
    or null (manifest only). ``rho_per_key``: per key, Re rho on the sampler's density, 0 where either basis is dark
    (a diagnostic; rho is not stored per key in the sidecar). Every per-key list has one record per declared key,
    phi outer and theta inner.
    """

    waveguide: dict[int, pd.DataFrame]
    far_field: dict[int, pd.DataFrame]
    blocks: dict[str, Any]
    raw_powers_per_key: list[dict[str, Any]]
    flux_over_abs_re_s: list[dict[str, Any]]
    rho_per_key: list[dict[str, Any]] = dataclasses.field(default_factory=list)


# --- Readers shared by the 2.0 producer and consumers (bbsim.compare, bbsim.reexport_check) -------------------------
#
# One implementation of each, so that a declared key, a row's key and a declared exit lattice are read the same way by
# the runner, the comparison and the re-export acceptance (HF-040, HF-052).

KEY_UNMATCHED = -1     # key_index: a row at no declared key
KEY_AMBIGUOUS = -2     # key_index: a row within KEY_TOL_DEG of two declared keys


def key_text(key: tuple[float, float]) -> str:
    """``(phi, theta) = (45, 135) deg``: how every diagnostic names an incidence key."""
    return f"(phi, theta) = ({key[0]:g}, {key[1]:g}) deg"


def finite_number(value: object, where: str) -> float:
    """A finite real number (int, float or their numpy types); a bool, a string, None, NaN or an infinity is refused."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float, np.integer, np.floating)):
        raise DatasetV2Error(f"{where} must be a number, got {value!r}")
    number = float(value)
    if not math.isfinite(number):
        raise DatasetV2Error(f"{where} is not finite: {value!r}")
    return number


def require(mapping: object, *path: str, where: str) -> Any:
    """The value at ``path`` in nested mappings; a missing step is refused, naming the whole path."""
    node = mapping
    for name in path:
        if not isinstance(node, Mapping) or name not in node:
            raise DatasetV2Error(f"{where} lacks {'.'.join(path)}")
        node = node[name]
    return node


def key_position(keys: Sequence[tuple[float, float]], phi: float, theta: float) -> int | None:
    """The index of the declared key within KEY_TOL_DEG of (phi, theta) on both angles, or None."""
    for k, (p, t) in enumerate(keys):
        if abs(p - phi) <= KEY_TOL_DEG and abs(t - theta) <= KEY_TOL_DEG:
            return k
    return None


def key_index(phis: np.ndarray, thetas: np.ndarray, keys: Sequence[tuple[float, float]]) -> np.ndarray:
    """The declared key index of every row (within KEY_TOL_DEG on both angles): KEY_UNMATCHED for a row at no
    declared key, KEY_AMBIGUOUS for a row within the tolerance of two declared keys."""
    phis = np.asarray(phis, dtype=float)
    thetas = np.asarray(thetas, dtype=float)
    index = np.full(phis.shape, KEY_UNMATCHED, dtype=int)
    for k, (phi, theta) in enumerate(keys):
        mask = (np.abs(phis - phi) <= KEY_TOL_DEG) & (np.abs(thetas - theta) <= KEY_TOL_DEG)
        taken = index != KEY_UNMATCHED
        index[mask & taken] = KEY_AMBIGUOUS
        index[mask & ~taken] = k
    return index


def row_keys(table: pd.DataFrame, keys: Sequence[tuple[float, float]], where: str) -> np.ndarray:
    """``key_index`` of a table's rows, refused unless every row lies at exactly one declared key and every declared
    key has rows."""
    try:
        phis = table["IWavePhi"].to_numpy(dtype=float)
        thetas = table["IWaveTheta"].to_numpy(dtype=float)
    except (KeyError, TypeError, ValueError) as exc:
        raise DatasetV2Error(f"{where}: unreadable IWavePhi / IWaveTheta: {exc}") from exc
    if not (np.all(np.isfinite(phis)) and np.all(np.isfinite(thetas))):
        raise DatasetV2Error(f"{where}: non-finite IWavePhi / IWaveTheta values")
    index = key_index(phis, thetas, keys)
    ambiguous = np.flatnonzero(index == KEY_AMBIGUOUS)
    if ambiguous.size:
        i = int(ambiguous[0])
        raise DatasetV2Error(f"{where}: rows at {key_text((float(phis[i]), float(thetas[i])))} match two declared keys")
    for k, key in enumerate(keys):
        if not np.any(index == k):
            raise DatasetV2Error(f"{where} has no rows at the declared key {key_text(key)}")
    unmatched = np.flatnonzero(index == KEY_UNMATCHED)
    if unmatched.size:
        first = (float(phis[unmatched[0]]), float(thetas[unmatched[0]]))
        raise DatasetV2Error(f"{where} has {unmatched.size} rows at keys outside the declared grid, the first at "
                             f"{key_text(first)}")
    return index


def declared_grid(record: object, where: str = "exit_field.grid") -> ExitGrid:
    """The exit lattice a sidecar declares (``exit_field.grid``: ``y_e`` and ``z_e``, each ``min``, ``max``, ``count``
    and ``step`` in metres) as an ``ExitGrid`` in mm, for ``node_weights``. The one parser of a declared grid.

    Refused: a missing axis or field; a min, max or step that is not a finite number; a step that is not positive; a
    count that is not an integer of at least 2 (with one node the first and last node coincide, so the wall weight
    and the step are undefined, as the producer refuses); a step that does not span min to max in count - 1 steps, to
    ``NODE_TOL`` of the step; counts that the lattice would not reproduce.
    """
    from .extract import ExitGrid                       # local: this module stays free of the extraction imports

    lo, hi, steps_mm, counts = [0.0], [0.0], [0.0], []
    for axis in ("y_e", "z_e"):
        at = f"{where}.{axis}"
        rec = require(record, axis, where=where)
        lo_m = finite_number(require(rec, "min", where=at), f"{at}.min")
        hi_m = finite_number(require(rec, "max", where=at), f"{at}.max")
        step_m = finite_number(require(rec, "step", where=at), f"{at}.step")
        count = require(rec, "count", where=at)
        if isinstance(count, (bool, np.bool_)) or not isinstance(count, (int, np.integer)):
            raise DatasetV2Error(f"{at}.count must be an integer, got {count!r}")
        if step_m <= 0.0:
            raise DatasetV2Error(f"{at}.step must be positive, got {step_m!r}")
        if count < 2:
            raise DatasetV2Error(f"{at}.count is {int(count)}: an exit lattice needs at least 2 nodes per axis "
                                 "(with fewer the first and last node coincide, so the wall weight and the step are "
                                 "undefined)")
        if abs((hi_m - lo_m) / (count - 1) - step_m) > NODE_TOL * step_m:
            raise DatasetV2Error(f"{at}: the step {step_m!r} m does not span {lo_m!r} to {hi_m!r} m in "
                                 f"{int(count) - 1} steps")
        lo.append(lo_m * 1e3)
        hi.append(hi_m * 1e3)
        steps_mm.append(step_m * 1e3)
        counts.append(int(count))
    grid = ExitGrid(lo_mm=tuple(lo), hi_mm=tuple(hi), resolution_mm=tuple(steps_mm), step_mm=tuple(steps_mm))
    if list(grid.counts[1:]) != counts:
        raise DatasetV2Error(f"{where}: the declared counts {counts} are not the lattice's {list(grid.counts[1:])}")
    return grid


def _positive_finite(name: str, value: object) -> float:
    try:
        number = float(value)                                   # type: ignore[arg-type]
    except (TypeError, ValueError):
        number = math.nan
    if isinstance(value, bool) or not (math.isfinite(number) and number > 0.0):
        raise DatasetV2Error(f"{name} must be a finite positive number, got {value!r}")
    return number


def _declared_keys(phi_values: Sequence[float], theta_values: Sequence[float]) -> list[tuple[float, float]]:
    """The declared incidence keys, phi outer and theta inner (the runner's extraction order)."""
    axes = []
    for name, values in (("incident phi", phi_values), ("incident theta", theta_values)):
        try:
            numbers = [float(v) for v in values]
        except (TypeError, ValueError) as exc:
            raise DatasetV2Error(f"{name} values must be numbers, got {values!r}") from exc
        if not numbers:
            raise DatasetV2Error(f"no {name} values declared")
        if not all(math.isfinite(v) for v in numbers):
            raise DatasetV2Error(f"{name} values must be finite, got {numbers}")
        ordered = sorted(numbers)
        for a, b in zip(ordered, ordered[1:]):
            if b - a <= KEY_TOL_DEG:
                raise DatasetV2Error(f"duplicate {name} value {a!r} in the declared grid {numbers}")
        axes.append(numbers)
    return [(phi, theta) for phi in axes[0] for theta in axes[1]]


def _checked_order_check(order_check: object) -> dict[str, Any]:
    """The ``column_order_check`` block as ``order_check_block`` builds it, with a passed reference."""
    if not isinstance(order_check, Mapping) or set(order_check) != {"rule", "method", "reference"}:
        raise DatasetV2Error(f"the column-order check block must hold exactly rule, method and reference, "
                             f"got {order_check!r}")
    if order_check["rule"] != ORDER_RULE or order_check["method"] not in ORDER_METHODS:
        raise DatasetV2Error(f"the column-order check block must have rule {ORDER_RULE!r} and a method in "
                             f"{ORDER_METHODS}, got {order_check['rule']!r} and {order_check['method']!r}")
    reference = order_check["reference"]
    if not isinstance(reference, Mapping) or reference.get("passed") is not True:
        raise DatasetV2Error("the column-order check block must embed a passed reference: a missing or failed "
                             "reference stops the run")
    return dict(order_check)


def _key_powers(powers: object, keys: list[tuple[float, float]], where: str) -> list[float]:
    """The per-key powers of one polarization, in declared-key order: exactly one finite value per declared key."""
    if not isinstance(powers, Mapping):
        raise DatasetV2Error(f"{where} must map (phi, theta) to W, got {type(powers).__name__}")
    out: list[float | None] = [None] * len(keys)
    for raw_key, value in powers.items():
        try:
            phi, theta = (float(x) for x in raw_key)
        except (TypeError, ValueError) as exc:
            raise DatasetV2Error(f"{where}: key {raw_key!r} is not a (phi, theta) pair") from exc
        k = key_position(keys, phi, theta)
        if k is None:
            raise DatasetV2Error(f"{where} has a value at {key_text((phi, theta))}, outside the declared grid")
        if out[k] is not None:
            raise DatasetV2Error(f"{where} has two values at {key_text(keys[k])}")
        try:
            number = float(value)
        except (TypeError, ValueError):
            number = math.nan
        if isinstance(value, bool) or not math.isfinite(number):
            raise DatasetV2Error(f"{where} at {key_text(keys[k])} is not a finite number: {value!r}")
        out[k] = number
    missing = [keys[k] for k, v in enumerate(out) if v is None]
    if missing:
        raise DatasetV2Error(f"{where} lacks the declared key {key_text(missing[0])} "
                             f"({len(missing)} declared keys missing)")
    return [float(v) for v in out]                              # type: ignore[arg-type]


def _check_polarization_column(table: pd.DataFrame, ephi: int, where: str) -> None:
    try:
        ephis = table["Ephi"].to_numpy(dtype=float)
    except (KeyError, TypeError, ValueError) as exc:
        raise DatasetV2Error(f"{where}: unreadable Ephi column: {exc}") from exc
    if not np.all(ephis == ephi):
        raise DatasetV2Error(f"{where}: the Ephi column holds {sorted(set(ephis.tolist()))}, not {ephi}: the "
                             "tables are not this polarization's")


def build(extracted: Mapping[int, Extracted], grid: ExitGrid, cross_section: Mapping[str, Any], incoming_w: float,
          phi_values: Sequence[float], theta_values: Sequence[float], ff_steps_deg: tuple[float, float],
          polarization_filter_limit_ghz: float, order_check: Mapping[str, Any]) -> DatasetV2Result:
    """The 2.0 tables and sidecar blocks from one job's raw extraction of both polarizations.

    Per declared key and polarization: the quadrature weights (trapezoid by each point's grid node, found from its
    Y and Z as BBRsim finds it, for a rectangle; uniform for a disc), the scale W and the dark flag
    (``scale_factor``), the sanity ratio flux / |Re S|, the far-field hemisphere power over the stored flux, and
    Re rho on the scaled fields. ``cell_90`` is computed once over the whole declared theta list, so a grid without
    theta = 90 deg is refused for every run (Ruling 4).

    Row order within a key is free for a rectangle: W and the blocks do not depend on it. Tables keep their rows
    in input order and every stored value in its position: OutgoingPower becomes the key's stored flux (0 when
    dark), IngoingPower ``incoming_w`` on every row, columns 9-14 times the key's W (1 when dark: the stored fields
    are unchanged), the labels become ``schema.WAVEGUIDE_COLUMNS_V2`` and ``InterceptFactor`` (c_k) is appended
    last. Far-field tables are returned unchanged.

    Stops (``DatasetV2Error``, nothing returned) on: a polarization other than exactly 0 and 1; a declared key
    missing from, or an undeclared key present in, any table or power map of either polarization; a non-finite
    flux, |Re S|, field or far-field value, or a negative |Re S|; a negative flux beyond the dark floor; a
    non-dark record whose tangential integral is not positive; different point counts of the two polarizations
    at a key; for a rectangle, a key that is not the whole lattice (a point off its node by more than ``NODE_TOL``
    of the step, outside the lattice, or two points at one node) and polarizations whose rows are not at the same
    nodes row by row (BBRsim pairs them by index); a positive flux with |Re S| = 0 at any key; and, on keys with
    theta >= 90 deg whose stored flux exceeds ``WEAK_REL`` times the polarization's largest stored flux,
    flux / |Re S| outside ``SANITY_BOUNDS``.
    Dark records (stored flux 0) are never inside that band. Ratios of keys below 90 deg are reported only.
    """
    if not isinstance(extracted, Mapping) or set(extracted) != {0, 1}:
        got = sorted(extracted) if isinstance(extracted, Mapping) else type(extracted).__name__
        raise DatasetV2Error(f"schema 2 needs the extraction of both polarizations, Ephi 0 and 1, got {got}")
    incoming = _positive_finite("the entrance power A*I (incoming_w)", incoming_w)
    limit = _positive_finite("polarization_filter_limit_ghz", polarization_filter_limit_ghz)
    try:
        dtheta_ff, dphi_ff = ff_steps_deg
    except (TypeError, ValueError) as exc:
        raise DatasetV2Error(f"the far-field steps must be (theta step, phi step) in deg, got {ff_steps_deg!r}") from exc
    dtheta_ff = _positive_finite("the far-field step in theta", dtheta_ff)
    dphi_ff = _positive_finite("the far-field step in phi", dphi_ff)
    order_block = _checked_order_check(order_check)

    keys = _declared_keys(phi_values, theta_values)
    cell90 = cell_90(theta_values)                              # Ruling 4: refuses a grid without theta = 90 deg
    band_lo, band_hi = cell_90_band_deg(theta_values)
    c_values = [c_key(phi, theta, theta_values) for phi, theta in keys]

    shape = cross_section.get("shape") if isinstance(cross_section, Mapping) else None
    if shape not in QUADRATURE_BY_SHAPE:
        raise DatasetV2Error(f"unknown exit cross-section shape {shape!r}: expected one of {sorted(QUADRATURE_BY_SHAPE)}")
    ny, nz = int(grid.counts[1]), int(grid.counts[2])
    h_y = _positive_finite("the exit lattice step along y_e", grid.step_mm[1]) * 1e-3
    h_z = _positive_finite("the exit lattice step along z_e", grid.step_mm[2]) * 1e-3

    # Per polarization: the tables' key index, the per-key powers and the raw field values, all checked.
    row_key: dict[int, np.ndarray] = {}
    ff_key: dict[int, np.ndarray] = {}
    values_all: dict[int, np.ndarray] = {}
    position_y: dict[int, np.ndarray] = {}
    position_z: dict[int, np.ndarray] = {}
    flux: dict[int, list[float]] = {}
    absolute: dict[int, list[float]] = {}
    for e in (0, 1):
        ex = extracted[e]
        wg_where, ff_where = f"the waveguide table of Ephi={e}", f"the far-field table of Ephi={e}"
        if list(ex.waveguide.columns) != list(schema.WAVEGUIDE_COLUMNS):
            raise DatasetV2Error(f"{wg_where} must carry the 1.0 columns {schema.WAVEGUIDE_COLUMNS}, got "
                                 f"{list(ex.waveguide.columns)}")
        if list(ex.far_field.columns) != list(schema.FAR_FIELD_COLUMNS):
            raise DatasetV2Error(f"{ff_where} must carry the columns {schema.FAR_FIELD_COLUMNS}, got "
                                 f"{list(ex.far_field.columns)}")
        _check_polarization_column(ex.waveguide, e, wg_where)
        _check_polarization_column(ex.far_field, e, ff_where)
        row_key[e] = row_keys(ex.waveguide, keys, wg_where)
        ff_key[e] = row_keys(ex.far_field, keys, ff_where)
        flux[e] = _key_powers(ex.flux_w, keys, f"the exit flux of Ephi={e}")
        absolute[e] = _key_powers(ex.abs_w, keys, f"the |Re S| integral of Ephi={e}")
        for k, a in enumerate(absolute[e]):
            if a < 0.0:
                raise DatasetV2Error(f"the |Re S| integral of Ephi={e} at {key_text(keys[k])} is negative: {a!r}")
        try:
            values_all[e] = ex.waveguide.iloc[:, schema.EXIT_FIELD_SLOTS].to_numpy(dtype=float)
        except (TypeError, ValueError) as exc:
            raise DatasetV2Error(f"{wg_where} has non-numeric exit-field values: {exc}") from exc
        bad = ~np.all(np.isfinite(values_all[e]), axis=1)
        if bad.any():
            raise DatasetV2Error(f"{wg_where} at {key_text(keys[row_key[e][np.flatnonzero(bad)[0]]])} has "
                                 f"{int(bad.sum())} rows with non-finite exit-field values")
        if shape == "rectangle":                                # a rectangle's weights follow each point's node
            try:
                position_y[e] = ex.waveguide["Y"].to_numpy(dtype=float)
                position_z[e] = ex.waveguide["Z"].to_numpy(dtype=float)
            except (KeyError, TypeError, ValueError) as exc:
                raise DatasetV2Error(f"{wg_where}: unreadable exit-point Y / Z: {exc}") from exc
        for column in _FAR_FIELD_NUMERIC:                       # column by column: the far field can be large
            try:
                finite = np.isfinite(ex.far_field[column].to_numpy(dtype=float))
            except (TypeError, ValueError) as exc:
                raise DatasetV2Error(f"{ff_where} has non-numeric values in {column}: {exc}") from exc
            if not finite.all():
                raise DatasetV2Error(f"{ff_where} has non-finite values in {column}, the first at "
                                     f"{key_text(keys[ff_key[e][np.flatnonzero(~finite)[0]]])}")

    # Per key and polarization: W, dark, the sanity ratio, ff_over_flux; then rho on the pair.
    n_keys = len(keys)
    scale = {e: [1.0] * n_keys for e in (0, 1)}
    dark = {e: [False] * n_keys for e in (0, 1)}
    stored = {e: [0.0] * n_keys for e in (0, 1)}
    ratio: dict[int, list[float | None]] = {e: [None] * n_keys for e in (0, 1)}
    ff_over: dict[int, list[float | None]] = {e: [None] * n_keys for e in (0, 1)}
    rho_values = [0.0] * n_keys
    for k, key in enumerate(keys):
        rows = {e: np.flatnonzero(row_key[e] == k) for e in (0, 1)}
        if rows[0].size != rows[1].size:
            raise DatasetV2Error(f"at {key_text(key)} Ephi=0 has {rows[0].size} exit points and Ephi=1 has "
                                 f"{rows[1].size}: BBRsim pairs the polarizations' rows by index")
        weights = {}
        if shape == "rectangle":
            # Each point's weight comes from its own lattice node, as BBRsim assigns it, never from its row.
            nodes = {}
            for e in (0, 1):
                try:
                    nodes[e] = _lattice_nodes(position_y[e][rows[e]], position_z[e][rows[e]], grid)
                except DatasetV2Error as exc:
                    raise DatasetV2Error(f"Ephi={e} at {key_text(key)}: {exc}") from exc
            if not (np.array_equal(nodes[0][0], nodes[1][0]) and np.array_equal(nodes[0][1], nodes[1][1])):
                raise DatasetV2Error(f"at {key_text(key)} Ephi=0 and Ephi=1 hold their exit points at different "
                                     "lattice nodes, row by row: BBRsim pairs the polarizations' rows by index")
            weights = {e: _node_weights(*nodes[e]) for e in (0, 1)}
        else:
            try:
                weights = {e: quadrature_weights(shape, ny, nz, int(rows[e].size)) for e in (0, 1)}
            except DatasetV2Error as exc:
                raise DatasetV2Error(f"at {key_text(key)}: {exc}") from exc
        scaled = {}
        for e in (0, 1):
            where = f"Ephi={e} at {key_text(key)}"
            values = values_all[e][rows[e]]
            try:
                w, is_dark = scale_factor(values, weights[e], h_y, h_z, flux[e][k], incoming)
                ratio[e][k] = sanity_ratio(flux[e][k], absolute[e][k], incoming)
            except DatasetV2Error as exc:
                raise DatasetV2Error(f"{where}: {exc}") from exc
            scale[e][k], dark[e][k] = float(w), bool(is_dark)
            stored[e][k] = 0.0 if is_dark else flux[e][k]
            if not is_dark:
                ff_rows = extracted[e].far_field.iloc[np.flatnonzero(ff_key[e] == k)]
                ff_over[e][k] = hemisphere_power_w(ff_rows, dtheta_ff, dphi_ff) / stored[e][k]
            scaled[e] = values * w
        try:
            rho_values[k] = rho(scaled[0], scaled[1], weights[0], dark[0][k], dark[1][k])
        except DatasetV2Error as exc:
            raise DatasetV2Error(f"at {key_text(key)}: {exc}") from exc

    # Sanity bound (CDX-009) on reachable keys of non-weak flux; dark records are stored as 0 and never inside.
    lo_bound, hi_bound = SANITY_BOUNDS
    failures = []
    for e in (0, 1):
        largest = max(stored[e])
        for k, (phi, theta) in enumerate(keys):
            if theta < 90.0 - _THETA_90_TOL_DEG or not stored[e][k] > WEAK_REL * largest:
                continue
            r = ratio[e][k]
            if r is None or not lo_bound <= r <= hi_bound:
                failures.append(f"Ephi={e} at {key_text((phi, theta))}: flux / |Re S| = {r!r}")
    if failures:
        raise DatasetV2Error(f"exit flux over |Re S| outside [{lo_bound}, {hi_bound}] on reachable keys: "
                             + "; ".join(failures))

    # The 2.0 tables.
    tables: dict[int, pd.DataFrame] = {}
    for e in (0, 1):
        idx = row_key[e]
        out = extracted[e].waveguide.copy()
        out.columns = list(schema.WAVEGUIDE_COLUMNS_V2[:15])
        out["OutgoingPower"] = np.asarray(stored[e], dtype=float)[idx]
        out["IngoingPower"] = incoming
        fields = values_all[e] * np.asarray(scale[e], dtype=float)[idx][:, None]
        for j, name in enumerate(schema.WAVEGUIDE_COLUMNS_V2[9:15]):
            out[name] = fields[:, j]
        out[schema.INTERCEPT_FACTOR] = np.asarray(c_values, dtype=float)[idx]
        tables[e] = out.loc[:, list(schema.WAVEGUIDE_COLUMNS_V2)]

    def head(k: int) -> dict[str, float]:
        return {"phi_deg": float(keys[k][0]), "theta_deg": float(keys[k][1])}

    blocks = {
        "transmittance": {
            "definition": TRANSMITTANCE_DEFINITION,
            "numerator": TRANSMITTANCE_NUMERATOR,
            "incoming": TRANSMITTANCE_INCOMING,
            "c_key": [{**head(k), "c_k": float(c_values[k])} for k in range(n_keys)],
            "cell_90": float(cell90),
            "cell_90_assumes": {"emitter": CELL_90_EMITTER, "cos_floor": COS_FLOOR, "band_deg": [band_lo, band_hi]},
            "above_one": TRANSMITTANCE_ABOVE_ONE,
            "per_key": [{**head(k), "c_k": float(c_values[k]), "polarizations": {
                str(e): {"t_flux": stored[e][k] / incoming, "outgoing_power_abs_re_s_w": absolute[e][k]}
                for e in (0, 1)}} for k in range(n_keys)],
        },
        "exit_field": {
            "column_order_check": order_block,
            "tangential_components": list(TANGENTIAL_COMPONENTS),
            "quadrature": QUADRATURE_BY_SHAPE[shape],
            "scale": {"rule": SCALE_RULE, "per_key": [
                {**head(k), "w0": scale[0][k], "w1": scale[1][k], "dark0": dark[0][k], "dark1": dark[1][k]}
                for k in range(n_keys)]},
            "sampler": SAMPLER,
        },
        "cross_term": {
            "density": CROSS_TERM_DENSITY,
            "definition": CROSS_TERM_DEFINITION,
            "valid_below_ghz": limit,
            "above_limit": CROSS_TERM_ABOVE_LIMIT,
        },
        "far_field": {
            "normalization": FAR_FIELD_NORMALIZATION,
            "ff_over_flux": [{**head(k), "polarizations": {str(e): ff_over[e][k] for e in (0, 1)}}
                             for k in range(n_keys)],
        },
    }
    raw_powers = [{**head(k), "polarizations": {
        str(e): {"exit_flux_re_s_dot_n_w": flux[e][k], "outgoing_power_abs_re_s_w": absolute[e][k]} for e in (0, 1)}}
        for k in range(n_keys)]
    ratios = [{**head(k), "polarizations": {str(e): ratio[e][k] for e in (0, 1)}} for k in range(n_keys)]
    rho_records = [{**head(k), "rho": float(rho_values[k])} for k in range(n_keys)]
    return DatasetV2Result(waveguide=tables, far_field={e: extracted[e].far_field for e in (0, 1)}, blocks=blocks,
                           raw_powers_per_key=raw_powers, flux_over_abs_re_s=ratios, rho_per_key=rho_records)
