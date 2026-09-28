"""Compare a candidate export dataset directory against a reference one.

Schema, keys, grids and IngoingPower must match; transmission must agree within
tolerance; field distributions are reported (warning only) because solver
versions differ.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from bbsim.schema import FAR_FIELD_COLUMNS, KEY_COLUMNS, WAVEGUIDE_COLUMNS


@dataclass(frozen=True)
class Tolerances:
    coord_m: float = 1e-9
    angle_deg: float = 1e-6
    ingoing_rel: float = 1e-9
    t_atol: float = 1e-3   # reference T at 45-degree incidence is 0.012-0.13; 0.01 would hide 80 % errors there
    t_rtol: float = 0.05
    t_noise_floor: float = 1e-6
    corr_warn: float = 0.95


@dataclass
class Check:
    name: str
    passed: bool
    severity: str  # "fail" | "warn" | "info"
    detail: str
    values: dict = field(default_factory=dict)


@dataclass
class ComparisonResult:
    candidate: str
    reference: str
    checks: list[Check]
    tolerances: dict = field(default_factory=dict)

    @property
    def passed(self) -> bool:
        return all(c.passed for c in self.checks if c.severity == "fail")

    def to_dict(self) -> dict:
        return {"candidate": self.candidate, "reference": self.reference, "passed": self.passed,
                "tolerances": dict(self.tolerances), "checks": [asdict(c) for c in self.checks]}

    def format_table(self) -> str:
        tol = " ".join(f"{k}={v}" for k, v in self.tolerances.items())
        lines = [f"candidate: {self.candidate}", f"reference: {self.reference}", f"tolerances: {tol}", ""]
        for c in self.checks:
            status = "PASS" if c.passed else ("FAIL" if c.severity == "fail" else "WARN")
            lines.append(f"{status:4} {c.name:26} {c.detail}")
        lines += ["", "RESULT: " + ("PASS" if self.passed else "FAIL")]
        return "\n".join(lines)


def _keys(df: pd.DataFrame) -> set[tuple[float, float]]:
    return {(float(a), float(b)) for a, b in df[KEY_COLUMNS].drop_duplicates().itertuples(index=False, name=None)}


def _check_columns(table: str, cand: pd.DataFrame, ref: pd.DataFrame, expected: list[str]) -> Check:
    ok = list(cand.columns) == expected and list(ref.columns) == expected
    detail = "column order matches Geant4 contract" if ok else \
        f"candidate {list(cand.columns)} / reference {list(ref.columns)} / expected {expected}"
    return Check(f"{table}.columns", ok, "fail", detail)


def _check_constant(table: str, col: str, cand: pd.DataFrame, ref: pd.DataFrame) -> Check:
    cv, rv = sorted(cand[col].astype(str).unique()), sorted(ref[col].astype(str).unique())
    ok = cv == rv and len(cv) == 1
    return Check(f"{table}.{col}", ok, "fail", f"{cv} vs {rv}", {"candidate": cv, "reference": rv})


def _check_keys(table: str, cand: pd.DataFrame, ref: pd.DataFrame) -> Check:
    ck, rk = _keys(cand), _keys(ref)
    ok = ck == rk
    detail = f"{len(ck)} incident angles" if ok else f"missing {sorted(rk - ck)} extra {sorted(ck - rk)}"
    return Check(f"{table}.keys", ok, "fail", detail, {"missing": sorted(rk - ck), "extra": sorted(ck - rk)})


def _check_ingoing(cand: pd.DataFrame, ref: pd.DataFrame, tol: Tolerances) -> Check:
    c, r = float(cand["IngoingPower"].iloc[0]), float(ref["IngoingPower"].iloc[0])
    rel = abs(c - r) / abs(r) if r else float("inf")
    return Check("waveguide.ingoing_power", rel <= tol.ingoing_rel, "fail",
                 f"{c:.9e} vs {r:.9e} (rel {rel:.2e})", {"candidate": c, "reference": r, "rel": rel})


def _transmission(df: pd.DataFrame) -> dict[tuple[float, float], float]:
    first = df.drop_duplicates(KEY_COLUMNS)
    return {(float(a), float(b)): float(o) / float(i)
            for a, b, o, i in first[KEY_COLUMNS + ["OutgoingPower", "IngoingPower"]].itertuples(index=False, name=None)}


def _check_transmission(cand: pd.DataFrame, ref: pd.DataFrame, tol: Tolerances) -> Check:
    tc, tr = _transmission(cand), _transmission(ref)
    rows, failures, noise = [], [], []
    worst_rel, worst_key = 0.0, None
    for key in sorted(set(tc) & set(tr)):
        c, r = tc[key], tr[key]
        if c < tol.t_noise_floor and r < tol.t_noise_floor:
            noise.append(key)
            rows.append({"key": list(key), "candidate": c, "reference": r, "status": "noise"})
            continue
        rel = abs(c - r) / abs(r) if r else float("inf")
        if rel > worst_rel:
            worst_rel, worst_key = rel, key
        ok = abs(c - r) <= max(tol.t_atol, tol.t_rtol * abs(r))
        rows.append({"key": list(key), "candidate": c, "reference": r, "rel_diff": rel,
                     "status": "pass" if ok else "fail"})
        if not ok:
            failures.append(key)
    n_pass = len(rows) - len(failures) - len(noise)
    detail = f"{n_pass} pass, {len(noise)} below noise floor, {len(failures)} fail"
    if failures:
        detail += f": {failures}"
    if worst_key is not None:
        detail += f"; max rel diff {worst_rel:.3f} at {worst_key}"
    return Check("waveguide.transmission", not failures, "fail", detail,
                 {"rows": rows, "max_rel_diff": worst_rel, "max_rel_diff_key": list(worst_key) if worst_key else None})


def _magnitude(df: pd.DataFrame, prefixes: list[str]) -> np.ndarray:
    total = np.zeros(len(df))
    for p in prefixes:
        total += df[f"{p}_real"].to_numpy() ** 2 + df[f"{p}_imag"].to_numpy() ** 2
    return np.sqrt(total)


def _check_points_and_fields(table: str, cand: pd.DataFrame, ref: pd.DataFrame, coord_cols: list[str],
                             prefixes: list[str], coord_tol: float, tol: Tolerances) -> list[Check]:
    point_fail: list[tuple[tuple[float, float], str]] = []
    worst_coord, min_corr, max_nrmsd = 0.0, 1.0, 0.0
    compared = 0
    for key in sorted(_keys(cand) & _keys(ref)):
        mc = (cand[KEY_COLUMNS[0]] == key[0]) & (cand[KEY_COLUMNS[1]] == key[1])
        mr = (ref[KEY_COLUMNS[0]] == key[0]) & (ref[KEY_COLUMNS[1]] == key[1])
        c = cand[mc].sort_values(coord_cols).reset_index(drop=True)
        r = ref[mr].sort_values(coord_cols).reset_index(drop=True)
        if len(c) != len(r):
            point_fail.append((key, f"{len(c)} vs {len(r)} points"))
            continue
        dc = float(np.max(np.abs(c[coord_cols].to_numpy(dtype=float) - r[coord_cols].to_numpy(dtype=float))))
        worst_coord = max(worst_coord, dc)
        if dc > coord_tol:
            point_fail.append((key, f"coordinate mismatch {dc:.3e}"))
            continue
        compared += 1
        mag_c, mag_r = _magnitude(c, prefixes), _magnitude(r, prefixes)
        scale = float(np.max(mag_r)) or 1.0
        max_nrmsd = max(max_nrmsd, float(np.sqrt(np.mean((mag_c - mag_r) ** 2)) / scale))
        if np.std(mag_c) > 0 and np.std(mag_r) > 0:
            min_corr = min(min_corr, float(np.corrcoef(mag_c, mag_r)[0, 1]))
    points = Check(f"{table}.points", not point_fail, "fail",
                   f"grids identical (worst coordinate diff {worst_coord:.2e})" if not point_fail
                   else f"{len(point_fail)} angles differ: {point_fail[:5]}",
                   {"failures": [(list(k), d) for k, d in point_fail]})
    if compared == 0:
        fields = Check(f"{table}.fields", False, "warn", "no angles with matching grids; fields not compared",
                       {"compared_angles": 0})
    else:
        fields = Check(f"{table}.fields", min_corr >= tol.corr_warn, "warn",
                       f"{compared} angles: min |E| correlation {min_corr:.4f}, max normalized RMS diff {max_nrmsd:.4f}",
                       {"compared_angles": compared, "min_corr": min_corr, "max_nrmsd": max_nrmsd})
    return [points, fields]


def compare_datasets(candidate_dir: Path, reference_dir: Path, tol: Tolerances = Tolerances()) -> ComparisonResult:
    candidate_dir, reference_dir = Path(candidate_dir), Path(reference_dir)
    checks: list[Check] = []

    cw, rw = pd.read_csv(candidate_dir / "waveguide.csv"), pd.read_csv(reference_dir / "waveguide.csv")
    cols = _check_columns("waveguide", cw, rw, WAVEGUIDE_COLUMNS)
    checks.append(cols)
    if cols.passed:
        checks += [_check_constant("waveguide", "Freq", cw, rw), _check_constant("waveguide", "Ephi", cw, rw),
                   _check_keys("waveguide", cw, rw), _check_ingoing(cw, rw, tol), _check_transmission(cw, rw, tol)]
        checks += _check_points_and_fields("waveguide", cw, rw, ["X", "Y", "Z"], ["Ex", "Ey", "Ez"], tol.coord_m, tol)

    cf, rf = pd.read_csv(candidate_dir / "far_field.csv"), pd.read_csv(reference_dir / "far_field.csv")
    cols = _check_columns("far_field", cf, rf, FAR_FIELD_COLUMNS)
    checks.append(cols)
    if cols.passed:
        checks += [_check_constant("far_field", "Freq", cf, rf), _check_constant("far_field", "Ephi", cf, rf),
                   _check_keys("far_field", cf, rf)]
        checks += _check_points_and_fields("far_field", cf, rf, ["Theta", "Phi"], ["rEtheta", "rEphi"], tol.angle_deg, tol)

    return ComparisonResult(str(candidate_dir), str(reference_dir), checks, asdict(tol))
