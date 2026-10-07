"""Parsers for HFSS field exports and a bounded wait for export files.

Formats are exactly those consumed by legacy/bbsim1freq.py (read_hfss_field,
read_hfss_far_field). The readers never delete files; callers own scratch paths.
"""
from __future__ import annotations

import logging
import math
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from bbsim.schema import EXIT_FIELD_COLUMNS

log = logging.getLogger(__name__)


class FieldFileError(ValueError):
    """An export file is malformed, truncated, or inconsistent with the requested grid."""


class ExportTimeoutError(RuntimeError):
    """An export file did not appear and settle within the timeout."""


def read_exit_field_fld(path: Path, keep_unsolved: bool = False) -> pd.DataFrame:
    """Read an ExportOnGrid/CalculatorWrite .fld: two header lines, then 9 floats per line.

    A row whose tokens after the three coordinates are all ``nan`` is a grid point outside the solved
    region and is skipped, whatever their number: the form AEDT 2025 R2 writes for such a point is
    unverified, and the legacy reader skipped any row ending in ``nan``. Skipped rows without the full
    9 tokens are reported in one warning per file. Any other non-finite value is an error: BBRsim
    rejects non-finite fields at load (BBR013).

    With ``keep_unsolved`` those rows are returned instead, in file order, with nan field values and
    finite coordinates, for a caller that maps every lattice point (``retain_disc_lattice``).
    """
    lines = Path(path).read_text().splitlines()
    rows: list[list[float]] = []
    unsolved: list[int] = []                 # indices into rows of the all-nan rows kept with keep_unsolved
    odd_outside: list[tuple[int, int]] = []  # (line number, token count) of all-nan rows without 9 tokens
    for lineno, line in enumerate(lines[2:], start=3):
        parts = line.split()
        if not parts:
            continue
        field_tokens = parts[3:]
        if field_tokens and all(token.lower() == "nan" for token in field_tokens):
            if len(parts) != len(EXIT_FIELD_COLUMNS):
                odd_outside.append((lineno, len(parts)))
            if keep_unsolved:
                try:
                    point = [float(token) for token in parts[:3]]
                except ValueError:
                    raise FieldFileError(f"{path}: line {lineno} has a non-numeric coordinate: {line.strip()!r}") from None
                if not all(math.isfinite(v) for v in point):
                    raise FieldFileError(f"{path}: line {lineno} has a non-finite coordinate: {line.strip()!r}")
                unsolved.append(len(rows))
                rows.append(point + [math.nan] * (len(EXIT_FIELD_COLUMNS) - 3))
            continue
        if len(parts) != len(EXIT_FIELD_COLUMNS):
            raise FieldFileError(f"{path}: line {lineno} has {len(parts)} fields, expected {len(EXIT_FIELD_COLUMNS)}")
        try:
            values = [float(token) for token in parts]
        except ValueError:
            raise FieldFileError(f"{path}: line {lineno} has a non-numeric field: {line.strip()!r}") from None
        if not all(math.isfinite(v) for v in values):
            raise FieldFileError(f"{path}: line {lineno} has a non-finite value outside an all-nan row: {line.strip()!r}")
        rows.append(values)
    if odd_outside:
        counts = dict(sorted(Counter(n for _, n in odd_outside).items()))
        fate = ("kept as unsolved lattice points (nan field) for the caller's rim rule" if keep_unsolved
                else "skipped as outside points")
        log.warning("%s: %d all-nan rows with a token count other than %d (counts %s, first at line %d) %s; HFSS "
                    "writes unsolved points in another form, or the export was cut short",
                    path, len(odd_outside), len(EXIT_FIELD_COLUMNS), counts, odd_outside[0][0], fate)
    if len(rows) == len(unsolved):
        raise FieldFileError(f"{path}: no numeric data rows")
    return pd.DataFrame(rows, columns=EXIT_FIELD_COLUMNS)


def read_calculator_scalar(path: Path) -> float:
    """Read the value a Fields Calculator ``CalculatorWrite`` wrote for a scalar expression.

    The file holds header lines and then the value; PyAEDT 1.7.0 reads its last line the same way.
    The last token of the last non-empty line is taken, so a ``x y z value`` row is read too.
    A missing or non-finite value is an error: the power enters the dataset BBRsim loads.
    """
    path = Path(path)
    lines = [line for line in path.read_text().splitlines() if line.strip()]
    if not lines:
        raise FieldFileError(f"{path}: no numeric value (empty file)")
    token = lines[-1].split()[-1]
    try:
        value = float(token)
    except ValueError:
        raise FieldFileError(f"{path}: no numeric value on the last line: {lines[-1].strip()!r}") from None
    if not math.isfinite(value):
        raise FieldFileError(f"{path}: non-finite value {token!r} on the last line")
    return value


def read_far_field_ffd(path: Path, theta_step_deg: float, phi_step_deg: float) -> pd.DataFrame:
    """Read an infinite-sphere .ffd export.

    Line 0: theta_start theta_end n_theta. Line 1: phi_start phi_end n_phi.
    Lines 4+: rEtheta_real rEtheta_imag rEphi_real rEphi_imag, theta outer, phi inner.
    """
    lines = Path(path).read_text().splitlines()
    if len(lines) < 4:
        raise FieldFileError(f"{path}: fewer than 4 header lines")
    theta_start, theta_end, n_theta = (float(v) for v in lines[0].split()[:3])
    phi_start, phi_end, n_phi = (float(v) for v in lines[1].split()[:3])
    n_theta, n_phi = int(n_theta), int(n_phi)

    theta_vals = np.arange(theta_start, theta_end + 1e-8, theta_step_deg)
    phi_vals = np.arange(phi_start, phi_end + 1e-8, phi_step_deg)
    if len(theta_vals) != n_theta:
        raise FieldFileError(f"{path}: header has {n_theta} theta points but step {theta_step_deg} gives {len(theta_vals)}")
    if len(phi_vals) != n_phi:
        raise FieldFileError(f"{path}: header has {n_phi} phi points but step {phi_step_deg} gives {len(phi_vals)}")

    data_lines = [line for line in lines[4:] if line.strip()]
    expected = n_theta * n_phi
    if len(data_lines) != expected:
        raise FieldFileError(f"{path}: expected {expected} data lines, found {len(data_lines)}")
    data = np.array([[float(v) for v in line.split()[:4]] for line in data_lines])

    theta_grid, phi_grid = np.meshgrid(theta_vals, phi_vals, indexing="ij")
    return pd.DataFrame({
        "Phi": phi_grid.ravel(),
        "Theta": theta_grid.ravel(),
        "rEphi_real": data[:, 2],
        "rEphi_imag": data[:, 3],
        "rEtheta_real": data[:, 0],
        "rEtheta_imag": data[:, 1],
    })


def wait_for_file(path: Path, timeout_s: float, poll_s: float = 0.05, stable_polls: int = 2) -> int:
    """Block until ``path`` exists, is non-empty, and its size is unchanged for ``stable_polls`` polls.

    Returns the final size. Raises ExportTimeoutError after ``timeout_s`` seconds.
    """
    path = Path(path)
    deadline = time.monotonic() + timeout_s
    last_size = -1
    stable = 0
    while time.monotonic() < deadline:
        if path.exists():
            size = path.stat().st_size
            if size > 0 and size == last_size:
                stable += 1
            else:
                stable = 0
            last_size = size
            if stable >= stable_polls:
                return size
        time.sleep(poll_s)
    raise ExportTimeoutError(f"{path} did not complete within {timeout_s}s (last size {last_size} bytes)")
