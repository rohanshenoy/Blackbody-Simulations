"""Ideal-PEC circular-guide cutoffs for the round bore, counted per orientation (OR-053).

dataset_record.disc_entries lists each TE/TM family once with its degeneracy; a mode count for HFSS's
wave port must weight every family by that degeneracy, so a family with degeneracy 2 counts twice.
"""
from __future__ import annotations

import math

from bbsim.dataset_record import C as _C, disc_entries

_MARGIN = 2.0
_TE11_X = 1.8411837813406593  # first zero of J1'


def ideal_cutoffs_ghz(radius_m: float, limit_ghz: float) -> list[float]:
    """Sorted cutoffs up to limit_ghz, each repeated by its degeneracy."""
    out: list[float] = []
    for cutoff, _family, _n, _p, degeneracy in disc_entries(radius_m, limit_ghz):
        out.extend([cutoff] * degeneracy)
    return sorted(out)


def propagating_count(radius_m: float, frequency_ghz: float) -> int:
    """Number of ideal modes (per orientation) with cutoff at or below the frequency."""
    return len(ideal_cutoffs_ghz(radius_m, frequency_ghz))


def cutoff_distance(radius_m: float, frequency_ghz: float) -> float:
    """Smallest |f/fc - 1| over the ideal cutoffs, including the next one above f."""
    lowest = _TE11_X * _C / (2 * math.pi * radius_m) / 1e9
    limit = max(_MARGIN * frequency_ghz, _MARGIN * lowest)
    return min(abs(frequency_ghz / c - 1.0) for c in ideal_cutoffs_ghz(radius_m, limit))
