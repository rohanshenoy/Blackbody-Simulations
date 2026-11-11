"""The round bore's ideal modes; the ports-only mode study's export, self-checks and convergence comparison (OR-053).

dataset_record.disc_entries lists each TE/TM family once with its degeneracy; a mode count for HFSS's
wave port must weight every family by that degeneracy, so a family with degeneracy 2 counts twice.

The export, ``modes.csv``, is what OR-053's acceptance check reads: ``# key=value`` provenance lines, then one row per
frequency, port and port mode with gamma = alpha + j beta as AEDT reports it. The producer's self-checks are the rows
it refuses on before anything is handed over. The convergence comparison of two exports pairs each export's
propagating modes with the ideal modes by the acceptance check's own rule, a joint assignment on implied cutoff and
alpha, never by AEDT's mode index; modes of equal ideal cutoff are compared by their summed alpha, which a mesh that
mixes a degenerate set leaves unchanged. Run ``python -m bbsim.modes compare A.csv B.csv --bore-um D``.
"""
from __future__ import annotations

import argparse
import csv
import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple

from bbsim.dataset_record import C as _C, disc_entries, mode_label

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


class GammaRow(NamedTuple):
    """One port mode's propagation constant read from AEDT: alpha = Re Gamma (Np/m), beta = Im Gamma (rad/m)."""
    frequency_ghz: float
    port: int
    mode: int
    alpha_np_per_m: float
    beta_rad_per_m: float


# --- The ideal modes and the wall model the convergence comparison pairs by --------------------------------------

CSV_COLUMNS = ("frequency_ghz", "port", "mode", "alpha_np_per_m", "beta_rad_per_m")
PROVENANCE_KEYS = ("dataset_id", "bore_um", "wall_kind", "conductivity_s_per_m", "thickness_um", "modes_per_port",
                   "bbsim_commit", "tree_clean", "aedt_version", "pyaedt_version")
PORTS = (1, 2)
FREQUENCY_RTOL = 1e-9        # a row's frequency against the run's, after the export's %.10g
SIGN_TOL = 1e-9              # a pec wall's alpha and beta may fall below zero by SIGN_TOL k (round-off); a lossy
                             # wall's may not fall below zero at all, since the acceptance check refuses that on read
BETA_K_TOL = 1e-6            # beta <= k (1 + BETA_K_TOL)
PEC_ALPHA_TOL = 1e-6         # a pec wall: alpha <= PEC_ALPHA_TOL beta for every propagating mode
IMPLIED_CUTOFF_TOL = 0.01    # implied cutoffs against the ideal ones, relative
PORTS_AGREE_TOL = 1e-3       # port 2's propagating betas against port 1's, relative (a warning)
GROUP_RTOL = 1e-6            # ideal cutoffs equal to this form one group (degenerate orientations, TE0m with TM1m)
BORE_RTOL = 1e-9             # compare's --bore-um against each export's bore_um
# The acceptance check's pairing scales (OR-053): the implied cutoff's tolerance, and alpha's at f >= NEAR f_c and
# nearer cutoff. The pairing must stay the acceptance check's, so these track it.
FC_TOL = 0.005
ALPHA_TOL = 0.05
ALPHA_TOL_NEAR = 0.25
NEAR = 1.2
_MU0 = 1.25663706212e-6      # H/m, CODATA 2018, the acceptance check's wall model's value
_Z0 = _MU0 * _C
_NAMED = 8                   # offenders named in a check row's detail


class IdealMode(NamedTuple):
    """One orientation of an ideal-PEC circular-guide mode: family "TE" or "TM", azimuthal order n, radial order p."""
    cutoff_ghz: float
    family: str
    n: int
    p: int

    @property
    def label(self) -> str:
        return mode_label(self.family, self.n, self.p)


def ideal_modes(radius_m: float, limit_ghz: float) -> list[IdealMode]:
    """Every orientation of each mode with cutoff up to limit_ghz, by rising cutoff; equal cutoffs TE first, then by
    n and p. Entry i has the cutoff of ideal_cutoffs_ghz(radius_m, limit_ghz)[i]."""
    out: list[IdealMode] = []
    for cutoff, family, n, p, degeneracy in disc_entries(radius_m, limit_ghz):
        out.extend([IdealMode(cutoff, family, n, p)] * degeneracy)
    return out


def _k(frequency_ghz: float) -> float:
    return 2 * math.pi * frequency_ghz * 1e9 / _C


def model_alpha_np_per_m(mode: IdealMode, radius_m: float, frequency_ghz: float, rho_ohm_m: float) -> float:
    """The mode's field attenuation by a good-conductor wall of resistivity rho, Np/m: R_s / (a Z0 sqrt(1 - (fc/f)^2))
    times (fc/f)^2 + n^2 / (x^2 - n^2) for TE and 1 for TM, R_s = sqrt(pi f mu0 rho), x the mode's Bessel zero
    (Jackson section 8.5; the acceptance check's model). NaN at or below cutoff, 0 for rho = 0."""
    x_root = mode.cutoff_ghz * 1e9 * 2 * math.pi * radius_m / _C
    x = 2 * math.pi * frequency_ghz * 1e9 * radius_m / _C
    if not x_root < x:
        return math.nan
    if rho_ohm_m == 0.0:
        return 0.0
    r2 = (x_root / x) ** 2
    factor = r2 + mode.n * mode.n / (x_root * x_root - mode.n * mode.n) if mode.family == "TE" else 1.0
    surface_resistance = math.sqrt(math.pi * frequency_ghz * 1e9 * _MU0 * rho_ohm_m)
    return surface_resistance / (radius_m * _Z0 * math.sqrt(1.0 - r2)) * factor


def implied_cutoff_ghz(row: GammaRow) -> float:
    """c sqrt(k^2 - (beta - alpha)^2) / (2 pi), GHz: a good-conductor wall raises beta by alpha to first order, and
    removing it leaves the perfect-conductor beta (the acceptance check's formula). Zero once beta - alpha exceeds k."""
    k = _k(row.frequency_ghz)
    b = row.beta_rad_per_m - row.alpha_np_per_m
    return _C * math.sqrt(max(k * k - b * b, 0.0)) / (2 * math.pi) / 1e9


def _propagating(rows: Iterable[GammaRow]) -> list[GammaRow]:
    return [r for r in rows if r.beta_rad_per_m > r.alpha_np_per_m]


# --- The export ---------------------------------------------------------------------------------------------------

def write_modes_csv(path: Path, rows: list[GammaRow], provenance: dict[str, str]) -> None:
    """``# key=value`` for each of PROVENANCE_KEYS in that order, the header, then the rows sorted by (frequency,
    port, mode): frequency %.10g, port and mode as integers, alpha and beta %.12g. Non-finite values are written as
    they are; the self-checks judge them. Refused before anything is written: a provenance key missing or unknown,
    or a value holding a line break."""
    missing = [k for k in PROVENANCE_KEYS if k not in provenance]
    unknown = sorted(set(provenance) - set(PROVENANCE_KEYS))
    if missing or unknown:
        raise ValueError(f"provenance keys missing: {', '.join(missing) or 'none'}; unknown: "
                         f"{', '.join(unknown) or 'none'}; expected {', '.join(PROVENANCE_KEYS)}")
    lines = []
    for key in PROVENANCE_KEYS:
        value = str(provenance[key])
        if "\n" in value or "\r" in value:
            raise ValueError(f"provenance {key} holds a line break: {value!r}")
        lines.append(f"# {key}={value}")
    lines.append(",".join(CSV_COLUMNS))
    for r in sorted(rows, key=lambda r: (r.frequency_ghz, r.port, r.mode)):
        lines.append("%.10g,%d,%d,%.12g,%.12g" % (r.frequency_ghz, r.port, r.mode, r.alpha_np_per_m,
                                                   r.beta_rad_per_m))
    Path(path).write_text("\n".join(lines) + "\n", newline="\n")


def read_modes_csv(path: Path) -> list[GammaRow]:
    """The rows of an export, in file order. ``#`` and blank lines are skipped and other columns ignored. ValueError
    names a line the csv module cannot parse, a missing column, a value that is not a number, a repeated (frequency,
    port, mode) or an export with no rows; NaN and inf pass (the self-checks judge them). OSError if the file cannot
    be read."""
    numbered = [(number, line) for number, line in enumerate(Path(path).read_text().splitlines(), start=1)
                if line.strip() and not line.lstrip().startswith("#")]
    reader = csv.DictReader([line for _, line in numbered])
    try:
        fieldnames = reader.fieldnames
    except csv.Error as exc:
        raise ValueError(f"{path}: line {numbered[0][0]}: not readable as CSV: {exc}") from exc
    missing = [c for c in CSV_COLUMNS if c not in (fieldnames or [])]
    if missing:
        raise ValueError(f"{path}: missing column(s) {', '.join(missing)}; expected {', '.join(CSV_COLUMNS)}")
    rows: list[GammaRow] = []
    seen: set[tuple[float, int, int]] = set()
    records = iter(reader)
    for number, _ in numbered[1:]:
        try:
            record = next(records)
        except StopIteration:
            break
        except csv.Error as exc:
            raise ValueError(f"{path}: line {number}: not readable as CSV: {exc}") from exc
        try:
            row = GammaRow(float(record["frequency_ghz"]), int(record["port"]), int(record["mode"]),
                           float(record["alpha_np_per_m"]), float(record["beta_rad_per_m"]))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{path}: line {number}: not a number: {exc}") from exc
        key = (row.frequency_ghz, row.port, row.mode)
        if key in seen:
            raise ValueError(f"{path}: line {number}: frequency {row.frequency_ghz:g} GHz, port {row.port}, "
                             f"mode {row.mode} repeated")
        seen.add(key)
        rows.append(row)
    if not rows:
        raise ValueError(f"{path}: no rows")
    return rows


def read_modes_provenance(path: Path) -> dict[str, str]:
    """The ``# key=value`` lines before the header, as write_modes_csv wrote them."""
    out: dict[str, str] = {}
    for line in Path(path).read_text().splitlines():
        if not line.startswith("#"):
            if line.strip():
                break
            continue
        body = line[1:]
        body = body[1:] if body.startswith(" ") else body
        key, sep, value = body.partition("=")
        if sep and key.strip():
            out[key.strip()] = value
    return out


# --- The producer's self-checks -----------------------------------------------------------------------------------

class CheckRow(NamedTuple):
    """One self-check: severity "fail" stops the producer (exit 1); "warn" is recorded only."""
    name: str
    passed: bool
    severity: str
    detail: str


def _where(row: GammaRow) -> str:
    return f"{row.frequency_ghz:g} GHz port {row.port} mode {row.mode}"


def _row(name: str, problems: list[str], passed_detail: str, severity: str = "fail") -> CheckRow:
    if not problems:
        return CheckRow(name, True, severity, passed_detail)
    more = len(problems) - _NAMED
    return CheckRow(name, False, severity, "; ".join(problems[:_NAMED]) + (f"; and {more} more" if more > 0 else ""))


def self_check(rows: Sequence[GammaRow], radius_m: float, modes: int, frequencies_ghz: Iterable[float],
               wall_kind: str) -> list[CheckRow]:
    """The producer's checks of an export before it is handed over (OR-053), in order: rows.count, values.finite,
    values.sign, pec.alpha_zero (pec walls only), beta.below_k, propagating.count, implied_cutoff, ports.agree (warn).
    "Propagating" is beta > alpha; counts are per orientation, as ideal_cutoffs_ghz. Ports are 1 and 2. values.sign
    refuses any negative alpha or beta unless the wall is pec, as the acceptance check refuses one when it reads the
    export; a pec export, which that check never reads, keeps a round-off tolerance of SIGN_TOL k."""
    rows = list(rows)
    freqs = [float(f) for f in frequencies_ghz]
    at: dict[tuple[float, int], list[GammaRow]] = {(f, p): [] for f in freqs for p in PORTS}
    stray: dict[tuple[float, int], int] = {}
    for r in rows:
        f = next((g for g in freqs if math.isclose(r.frequency_ghz, g, rel_tol=FREQUENCY_RTOL, abs_tol=0.0)), None)
        if f is None or r.port not in PORTS:
            stray[(r.frequency_ghz, r.port)] = stray.get((r.frequency_ghz, r.port), 0) + 1
        else:
            at[(f, r.port)].append(r)
    finite = [r for r in rows if all(math.isfinite(v) for v in (r.frequency_ghz, r.alpha_np_per_m, r.beta_rad_per_m))]
    out = [_rows_count(at, stray, modes), _values_finite(rows), _values_sign(finite, wall_kind == "pec")]
    if wall_kind == "pec":
        out.append(_pec_alpha_zero(finite))
    out += [_beta_below_k(finite), _propagating_count(at, radius_m), _implied_cutoff(at, radius_m), _ports_agree(at)]
    return out


def _rows_count(at, stray, modes: int) -> CheckRow:
    problems = []
    for (f, port), here in at.items():
        numbers = sorted(r.mode for r in here)
        if len(here) != modes:
            problems.append(f"{f:g} GHz port {port}: {len(here)} rows, expected {modes}")
        elif numbers != list(range(1, modes + 1)):
            problems.append(f"{f:g} GHz port {port}: mode numbers are not 1 to {modes} "
                            f"(from {numbers[0]} to {numbers[-1]}, {len(set(numbers))} distinct)")
    def order(item):
        (f, port), _count = item
        return (math.isnan(f), 0.0 if math.isnan(f) else f, port)

    for (f, port), count in sorted(stray.items(), key=order):
        problems.append(f"{f:g} GHz port {port}: {count} rows at a frequency or port the run did not ask for")
    return _row("rows.count", problems, f"{len(at)} (frequency, port) pairs of {modes} rows each, modes 1 to {modes}")


def _values_finite(rows) -> CheckRow:
    problems = [f"{_where(r)}: alpha {r.alpha_np_per_m!r}, beta {r.beta_rad_per_m!r}" for r in rows
                if not all(math.isfinite(v) for v in (r.frequency_ghz, r.alpha_np_per_m, r.beta_rad_per_m))]
    return _row("values.finite", problems, f"{len(rows)} rows, every frequency, alpha and beta finite")


def _values_sign(finite, pec: bool) -> CheckRow:
    problems = []
    for r in finite:
        floor = -SIGN_TOL * _k(r.frequency_ghz) if pec else 0.0
        allowed = (f"round-off allowed to {floor:.3g} on a pec wall" if pec
                   else "tolerance 0: the acceptance check refuses a negative value")
        for name, value in (("alpha", r.alpha_np_per_m), ("beta", r.beta_rad_per_m)):
            if value < floor:
                problems.append(f"{_where(r)}: {name} {value:.6g} below zero ({allowed})")
    return _row("values.sign", problems,
                f"no alpha or beta below -{SIGN_TOL:g} k (round-off on a pec wall)" if pec
                else "no alpha or beta below zero (tolerance 0, as the acceptance check)")


def _pec_alpha_zero(finite) -> CheckRow:
    problems = [f"{_where(r)}: alpha {r.alpha_np_per_m:.6g} Np/m on a pec wall, beta {r.beta_rad_per_m:.6g} rad/m"
                for r in _propagating(finite) if r.alpha_np_per_m > PEC_ALPHA_TOL * r.beta_rad_per_m]
    return _row("pec.alpha_zero", problems, f"every propagating mode has alpha <= {PEC_ALPHA_TOL:g} beta")


def _beta_below_k(finite) -> CheckRow:
    problems = [f"{_where(r)}: beta {r.beta_rad_per_m:.6g} rad/m above k = {_k(r.frequency_ghz):.6g} rad/m"
                for r in finite if r.beta_rad_per_m > _k(r.frequency_ghz) * (1 + BETA_K_TOL)]
    return _row("beta.below_k", problems, f"every beta <= k (1 + {BETA_K_TOL:g})")


def _propagating_count(at, radius_m: float) -> CheckRow:
    problems, counts = [], []
    for (f, port), here in at.items():
        n, ideal = len(_propagating(here)), propagating_count(radius_m, f)
        counts.append(ideal)
        if n != ideal:
            problems.append(f"{f:g} GHz port {port}: {n} modes with beta > alpha, the ideal guide {ideal}")
    span = f" ({min(counts)} to {max(counts)})" if counts else ""
    return _row("propagating.count", problems,
                f"the modes with beta > alpha match the ideal guide's count at every frequency and port{span}")


def _implied_cutoff(at, radius_m: float) -> CheckRow:
    problems = []
    for (f, port), here in at.items():
        implied = sorted(implied_cutoff_ghz(r) for r in _propagating(here))
        ideal = ideal_cutoffs_ghz(radius_m, f)
        if len(implied) != len(ideal):
            problems.append(f"{f:g} GHz port {port}: count mismatch, {len(implied)} implied cutoffs against "
                            f"{len(ideal)} ideal ones")
            continue
        off = [(abs(i / c - 1), i, c) for i, c in zip(implied, ideal) if not abs(i / c - 1) <= IMPLIED_CUTOFF_TOL]
        if off:
            _, i, c = max(off, key=lambda t: t[0] if math.isfinite(t[0]) else math.inf)
            problems.append(f"{f:g} GHz port {port}: {len(off)} of {len(ideal)} implied cutoffs beyond "
                            f"{IMPLIED_CUTOFF_TOL:.0%}, the worst {i:.6g} GHz against the ideal {c:.6g} GHz")
    return _row("implied_cutoff", problems,
                f"sorted implied cutoffs within {IMPLIED_CUTOFF_TOL:.0%} of the ideal ones at every frequency and port")


def _ports_agree(at) -> CheckRow:
    problems = []
    for f in sorted({f for f, _ in at}):
        b1, b2 = (sorted(r.beta_rad_per_m for r in _propagating(at[(f, p)])) for p in PORTS)
        if len(b1) != len(b2):
            problems.append(f"{f:g} GHz: port 1 has {len(b1)} propagating modes, port 2 {len(b2)}")
            continue
        # relative to the larger of the two, so a beta of zero (a propagating row with alpha below it) is compared too
        off = [d for d in (0.0 if x == y else abs(y - x) / max(abs(x), abs(y)) for x, y in zip(b1, b2))
               if not d <= PORTS_AGREE_TOL]
        if off:
            worst = max(off, key=lambda d: math.inf if math.isnan(d) else d)
            problems.append(f"{f:g} GHz: the ports' sorted propagating betas differ by up to {worst:.3g} relative")
    return _row("ports.agree", problems,
                f"port 2's sorted propagating betas within {PORTS_AGREE_TOL:g} of port 1's at every frequency",
                severity="warn")


# --- The convergence comparison -----------------------------------------------------------------------------------

def _pair_cost(row: GammaRow, mode: IdealMode, alpha_model: float, frequency_ghz: float, with_alpha: bool) -> float:
    value = ((implied_cutoff_ghz(row) / mode.cutoff_ghz - 1) / FC_TOL) ** 2
    if with_alpha:
        tol = ALPHA_TOL if frequency_ghz / mode.cutoff_ghz >= NEAR else ALPHA_TOL_NEAR
        value += (math.log(max(row.alpha_np_per_m, 1e-300) / alpha_model) / tol) ** 2
    return value if math.isfinite(value) else 1e30


def assign_to_ideal(rows: Sequence[GammaRow], radius_m: float,
                    rho_ohm_m: float | None) -> list[tuple[GammaRow, IdealMode, float]] | None:
    """The propagating rows of one frequency and port paired with the ideal modes by the acceptance check's rule
    (OR-053): the assignment that best agrees in implied cutoff and in alpha against the wall model at rho, each
    scaled by its tolerance. rho None (a pec wall) pairs by implied cutoff alone. Returns (row, ideal mode, model
    alpha) in the ideal modes' order, or None when the propagating count differs from the ideal guide's."""
    rows = list(rows)
    if len({(r.frequency_ghz, r.port) for r in rows}) > 1:
        raise ValueError("assign_to_ideal takes the rows of one frequency and port")
    if not rows:
        return []
    f = rows[0].frequency_ghz
    hfss = sorted(_propagating(rows), key=lambda r: r.mode)
    model = ideal_modes(radius_m, f)
    if len(hfss) != len(model):
        return None
    if not model:
        return []
    alphas = [math.nan if rho_ohm_m is None else model_alpha_np_per_m(m, radius_m, f, rho_ohm_m) for m in model]
    from scipy.optimize import linear_sum_assignment

    cost = [[_pair_cost(r, m, al, f, rho_ohm_m is not None) for m, al in zip(model, alphas)] for r in hfss]
    row_index, column = linear_sum_assignment(cost)
    return [(hfss[i], model[j], alphas[j]) for j, i in sorted(zip(column, row_index))]


class Change(NamedTuple):
    """A relative change from export A to export B at one frequency and port, naming the mode or group compared and
    AEDT's mode numbers in each export."""
    relative: float
    frequency_ghz: float
    port: int
    what: str
    modes_a: tuple[int, ...]
    modes_b: tuple[int, ...]
    value_a: float
    value_b: float


@dataclass(frozen=True)
class CompareResult:
    """The largest relative changes, each with its entry (None when nothing of that kind was compared), the (frequency,
    port) pairs compared and the reasons others were not. Recorded, not gated."""
    cutoff: Change | None
    group_alpha: Change | None
    single_alpha: Change | None
    compared: tuple[tuple[float, int], ...]
    skipped: tuple[str, ...]

    @property
    def max_cutoff_change(self) -> float | None:
        return None if self.cutoff is None else self.cutoff.relative

    @property
    def max_group_alpha_change(self) -> float | None:
        return None if self.group_alpha is None else self.group_alpha.relative

    @property
    def max_single_alpha_change(self) -> float | None:
        return None if self.single_alpha is None else self.single_alpha.relative


def _relative(a: float, b: float) -> float:
    if a == b:
        return 0.0
    return math.inf if a == 0 else abs(b - a) / abs(a)


def _worse(current: Change | None, new: Change) -> Change:
    def size(c: Change) -> float:
        return math.inf if math.isnan(c.relative) else c.relative
    return new if current is None or size(new) > size(current) else current


def _group_label(group: list[IdealMode]) -> str:
    counts: dict[str, int] = {}
    for m in group:
        counts[m.label] = counts.get(m.label, 0) + 1
    names = " + ".join(label if n == 1 else f"{label} x{n}" for label, n in counts.items())
    return f"{names} at {group[0].cutoff_ghz:.6g} GHz"


def _groups(pairs) -> list[tuple[int, int]]:
    """Index ranges of equal ideal cutoff (GROUP_RTOL) in the ideal modes' order."""
    out, start = [], 0
    for j in range(1, len(pairs) + 1):
        first = pairs[start][1].cutoff_ghz
        if j == len(pairs) or abs(pairs[j][1].cutoff_ghz - first) > GROUP_RTOL * first:
            out.append((start, j))
            start = j
    return out


def compare_exports(a: Sequence[GammaRow], b: Sequence[GammaRow], radius_m: float, *,
                    rho_ohm_m: float | None) -> CompareResult:
    """Two exports of one bore (a convergence repeat), compared at every (frequency, port) both hold. Each export's
    propagating modes are paired with the ideal modes by assign_to_ideal (the acceptance check's rule; rho the wall's
    resistivity, None for a pec wall). A group of equal ideal cutoff is compared by its summed alpha, which a mesh
    mixing the group leaves unchanged; a single mode by its alpha; implied cutoffs mode by mode, sorted within a
    group. Never by AEDT's mode index. A pair either export cannot be paired at is skipped with its reason."""
    if rho_ohm_m is not None and not (math.isfinite(rho_ohm_m) and rho_ohm_m > 0):
        raise ValueError(f"rho_ohm_m must be positive and finite, or None for a pec wall; not {rho_ohm_m!r}")
    by_a: dict[tuple[float, int], list[GammaRow]] = {}
    by_b: dict[tuple[float, int], list[GammaRow]] = {}
    for rows, into in ((a, by_a), (b, by_b)):
        for r in rows:
            into.setdefault((r.frequency_ghz, r.port), []).append(r)
    cutoff = group = single = None
    compared: list[tuple[float, int]] = []
    skipped: list[str] = []
    for key in sorted(set(by_a) | set(by_b)):
        f, port = key
        where = f"{f:g} GHz, port {port}"
        if key not in by_b or key not in by_a:
            skipped.append(f"{where}: only in export {'A' if key in by_a else 'B'}")
            continue
        pa, pb = (assign_to_ideal(by[key], radius_m, rho_ohm_m) for by in (by_a, by_b))
        if pa is None or pb is None:
            ideal = propagating_count(radius_m, f)
            bad = [f"export {name} has {len(_propagating(by[key]))} propagating modes, the ideal guide {ideal}"
                   for name, by, pairs in (("A", by_a, pa), ("B", by_b, pb)) if pairs is None]
            skipped.append(f"{where}: {'; '.join(bad)}; not compared")
            continue
        compared.append(key)
        for start, stop in _groups(pa):
            ga, gb = pa[start:stop], pb[start:stop]
            what = _group_label([m for _, m, _ in ga])
            for (ca, ma), (cb, mb) in zip(sorted((implied_cutoff_ghz(r), r.mode) for r, _, _ in ga),
                                          sorted((implied_cutoff_ghz(r), r.mode) for r, _, _ in gb)):
                cutoff = _worse(cutoff, Change(_relative(ca, cb), f, port, what, (ma,), (mb,), ca, cb))
            sa, sb = (sum(r.alpha_np_per_m for r, _, _ in g) for g in (ga, gb))
            change = Change(_relative(sa, sb), f, port, what, tuple(sorted(r.mode for r, _, _ in ga)),
                            tuple(sorted(r.mode for r, _, _ in gb)), sa, sb)
            if stop - start > 1:
                group = _worse(group, change)
            else:
                single = _worse(single, change)
    return CompareResult(cutoff, group, single, tuple(compared), tuple(skipped))


# --- The command line ---------------------------------------------------------------------------------------------

def _wall_rho(provenance: dict[str, str]) -> tuple[float | None, str]:
    kind = provenance.get("wall_kind")
    if kind == "pec":
        return None, "export A's wall is pec: modes paired by implied cutoff alone"
    if kind == "finite_conductivity":
        text = provenance.get("conductivity_s_per_m", "")
        try:
            sigma = float(text)
        except ValueError:
            sigma = math.nan
        if not (math.isfinite(sigma) and sigma > 0):
            raise ValueError(f"export A's conductivity_s_per_m is {text!r}")
        return 1.0 / sigma, f"1/conductivity_s_per_m of export A ({text} S/m)"
    raise ValueError(f"export A's wall_kind is {kind!r}, not pec or finite_conductivity")


def _describe(name: str, change: Change | None, unit: str) -> str:
    if change is None:
        return f"{name}: none compared"
    modes_a = ", ".join(str(m) for m in change.modes_a)
    modes_b = ", ".join(str(m) for m in change.modes_b)
    return (f"{name}: largest relative change {change.relative:.3g} at {change.frequency_ghz:g} GHz port "
            f"{change.port}, {change.what} (A mode {modes_a}: {change.value_a:.6g} {unit}; "
            f"B mode {modes_b}: {change.value_b:.6g} {unit})")


def _compare_command(args: argparse.Namespace) -> int:
    for name, value in (("--bore-um", args.bore_um), ("--rho", args.rho)):
        if value is not None and not (math.isfinite(value) and value > 0):
            print(f"compare: cannot use {name} {value!r}: it must be positive and finite")
            return 2
    try:
        a, b = read_modes_csv(args.a), read_modes_csv(args.b)
        provenance_a, provenance_b = read_modes_provenance(args.a), read_modes_provenance(args.b)
    except (OSError, ValueError) as exc:
        print(f"compare: cannot read an export: {exc}")
        return 2
    notes = []
    for name, provenance in (("A", provenance_a), ("B", provenance_b)):
        text = provenance.get("bore_um")
        if text is None:
            notes.append(f"note: export {name} carries no bore_um; --bore-um is not checked against it")
            continue
        try:
            bore = float(text)
        except ValueError:
            bore = math.nan
        if not math.isclose(bore, args.bore_um, rel_tol=BORE_RTOL, abs_tol=0.0):
            print(f"compare: cannot use --bore-um {args.bore_um:g}: export {name}'s bore_um is {text}")
            return 2
    if args.rho is not None:
        rho, source = args.rho, "--rho"
    else:
        try:
            rho, source = _wall_rho(provenance_a)
        except ValueError as exc:
            print(f"compare: cannot tell the wall's resistivity: {exc}; pass --rho")
            return 2
    result = compare_exports(a, b, args.bore_um * 1e-6 / 2, rho_ohm_m=rho)
    print(f"compare: A {args.a}")
    print(f"compare: B {args.b}")
    print(f"bore {args.bore_um:g} um; " + (f"rho {rho:.6g} ohm m ({source})" if rho is not None else source))
    for note in notes:
        print(note)
    for key in ("wall_kind", "conductivity_s_per_m"):
        if provenance_a.get(key) != provenance_b.get(key):
            print(f"note: {key} differs: A {provenance_a.get(key)!r}, B {provenance_b.get(key)!r}")
    print(f"compared {len(result.compared)} (frequency, port) pairs; skipped {len(result.skipped)}")
    print(_describe("implied cutoff", result.cutoff, "GHz"))
    print(_describe("group alpha sum", result.group_alpha, "Np/m"))
    print(_describe("single-mode alpha", result.single_alpha, "Np/m"))
    for reason in result.skipped:
        print(f"skipped: {reason}")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m bbsim.modes",
                                     description="The fibre bore's mode study (OR-053).")
    commands = parser.add_subparsers(dest="command", required=True)
    compare = commands.add_parser("compare", help="compare two mode exports (a convergence repeat); recorded, "
                                                  "not gated: exit 0 whatever the changes, 2 on an unusable file")
    compare.add_argument("a", type=Path, help="export A (the reference)")
    compare.add_argument("b", type=Path, help="export B")
    compare.add_argument("--bore-um", type=float, required=True, help="bore diameter, um")
    compare.add_argument("--rho", type=float, default=None,
                         help="wall resistivity for the pairing, ohm m; default 1/conductivity_s_per_m from export "
                              "A's provenance, and the implied cutoff alone for a pec wall")
    args = parser.parse_args(argv)
    return _compare_command(args)


if __name__ == "__main__":
    raise SystemExit(main())
