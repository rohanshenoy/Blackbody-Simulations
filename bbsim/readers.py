"""Parsers for HFSS field exports and a bounded wait for export files.

Formats are exactly those consumed by legacy/bbsim1freq.py (read_hfss_field,
read_hfss_far_field). The readers never delete files; callers own scratch paths.
"""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd

from bbsim.schema import EXIT_FIELD_COLUMNS


class FieldFileError(ValueError):
    """An export file is malformed, truncated, or inconsistent with the requested grid."""


class ExportTimeoutError(RuntimeError):
    """An export file did not appear and settle within the timeout."""


def read_exit_field_fld(path: Path) -> pd.DataFrame:
    """Read an ExportOnGrid/CalculatorWrite .fld: two header lines, then 9 floats per line.

    Lines whose last token is ``nan`` (points outside the solved region) are skipped.
    """
    lines = Path(path).read_text().splitlines()
    rows: list[list[float]] = []
    for lineno, line in enumerate(lines[2:], start=3):
        parts = line.split()
        if not parts:
            continue
        if parts[-1].lower() == "nan":
            continue
        if len(parts) != len(EXIT_FIELD_COLUMNS):
            raise FieldFileError(f"{path}: line {lineno} has {len(parts)} fields, expected {len(EXIT_FIELD_COLUMNS)}")
        rows.append([float(p) for p in parts])
    if not rows:
        raise FieldFileError(f"{path}: no numeric data rows")
    return pd.DataFrame(rows, columns=EXIT_FIELD_COLUMNS)


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
