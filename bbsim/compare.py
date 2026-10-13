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
    rho_atol: float = 0.05         # paired-polarization overlap Re rho; a relative sign flip moves it by about 2
    pair_dark_ratio: float = 1e-3  # a key whose weaker polarization carries less than this fraction of the
                                   # stronger one's T has no meaningful rho (cut-off polarization)
    pair_lit_ratio: float = 1e-2   # above this fraction both polarizations are lit; in between a key is not judged
                                   # for rho, so two solves straddling one threshold near a cutoff do not fail.
                                   # Legacy 500 GHz margins: lit keys >= 0.945, dark keys <= 1.9e-8.


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


EXIT_VALUES = ["Ex_real", "Ey_real", "Ez_real", "Ex_imag", "Ey_imag", "Ez_imag"]
FAR_VALUES = ["rEtheta_real", "rEtheta_imag", "rEphi_real", "rEphi_imag"]


def _overlap(a: np.ndarray, b: np.ndarray, weight: np.ndarray | None = None) -> float | None:
    """Re of sum(a . conj b) / sqrt(sum|a|^2 sum|b|^2) from stored real/imaginary parts.

    Equal to the sum of products of the stored reals however they pair into complex components, so the
    exit-field column order (HF-038) does not change it."""
    w = np.ones(len(a)) if weight is None else weight
    p0, p1 = float(np.sum(w[:, None] * a * a)), float(np.sum(w[:, None] * b * b))
    if p0 <= 0.0 or p1 <= 0.0:
        return None
    return float(np.sum(w[:, None] * a * b) / np.sqrt(p0 * p1))


def _pair_stats(tables, tol: Tolerances):
    """Per key of one (Ephi=0, Ephi=1) pair: T0, T1, exit and far-field Re rho, and whether the key is dark.

    Rows pair by order within a key, as BBRsim pairs them; a pairing BBRsim would refuse is a problem."""
    (w0, f0), (w1, f1) = tables
    stats, problems = {}, []
    keys0, keys1 = _keys(w0), _keys(w1)
    if keys0 != keys1:
        problems.append(f"incidence keys differ between the polarizations: {sorted(keys0 ^ keys1)}")
    for key in sorted(keys0 & keys1):
        def sel(df):
            return df[(df[KEY_COLUMNS[0]] == key[0]) & (df[KEY_COLUMNS[1]] == key[1])].reset_index(drop=True)
        a, b, fa, fb = sel(w0), sel(w1), sel(f0), sel(f1)
        if len(a) != len(b) or not np.allclose(a[["X", "Y", "Z"]].to_numpy(float), b[["X", "Y", "Z"]].to_numpy(float),
                                               rtol=0.0, atol=tol.coord_m):
            problems.append(f"{key}: exit points do not pair by row")
            continue
        if len(fa) != len(fb) or not np.allclose(fa[["Theta", "Phi"]].to_numpy(float), fb[["Theta", "Phi"]].to_numpy(float),
                                                 rtol=0.0, atol=tol.angle_deg):
            problems.append(f"{key}: far-field grids do not pair by row")
            continue
        t = [float(df["OutgoingPower"].iloc[0]) / float(df["IngoingPower"].iloc[0]) for df in (a, b)]
        sin_t = np.sin(np.radians(fa["Theta"].to_numpy(float)))
        stats[key] = {
            "T0": t[0], "T1": t[1],
            "state": ("dark" if max(t) < tol.t_noise_floor or min(t) < tol.pair_dark_ratio * max(t)
                      else "lit" if min(t) > tol.pair_lit_ratio * max(t) else "between"),
            "rho_exit": _overlap(a[EXIT_VALUES].to_numpy(float), b[EXIT_VALUES].to_numpy(float)),
            "rho_far": _overlap(fa[FAR_VALUES].to_numpy(float), fb[FAR_VALUES].to_numpy(float), sin_t),
        }
    return stats, problems


def _check_pair_ephi(candidate, reference) -> Check:
    wrong = []
    for label, tables in (("candidate", candidate), ("reference", reference)):
        for expected, (wg, ff) in enumerate(tables):
            found = sorted({int(v) for v in wg["Ephi"].unique()} | {int(v) for v in ff["Ephi"].unique()})
            if found != [expected]:
                wrong.append(f"{label} directory {expected + 1} holds Ephi {found}, expected [{expected}]")
    return Check("pair.ephi", not wrong, "fail", "Ephi=0 then Ephi=1 on both sides" if not wrong else "; ".join(wrong))


def _check_rho(name: str, field: str, cs: dict, rs: dict, tol: Tolerances) -> Check:
    rows, failures, dark = [], [], 0
    for key in sorted(set(cs) & set(rs)):
        c, r = cs[key], rs[key]
        if {c["state"], r["state"]} == {"dark", "lit"}:
            # A polarization cut off on one side and lit on the other is a difference, not a key to skip.
            rows.append({"key": list(key), "candidate_state": c["state"], "reference_state": r["state"], "status": "fail"})
            failures.append(key)
            continue
        if c["state"] != "lit" or r["state"] != "lit" or c[field] is None or r[field] is None:
            dark += 1
            continue
        ok = abs(c[field] - r[field]) <= tol.rho_atol
        rows.append({"key": list(key), "candidate": c[field], "reference": r[field], "status": "pass" if ok else "fail"})
        if not ok:
            failures.append(key)
    detail = f"{len(rows) - len(failures)} pass, {dark} dark or near-dark (not judged), {len(failures)} fail"
    if failures:
        detail += f": {failures}"
    return Check(name, not failures, "fail", detail, {"rows": rows})


def _check_mixed(cs: dict, rs: dict, tol: Tolerances) -> Check:
    rows, failures = [], []
    for key in sorted(set(cs) & set(rs)):
        pred = {}
        for label, s in (("candidate", cs[key]), ("reference", rs[key])):
            cross = np.sqrt(s["T0"] * s["T1"]) * (s["rho_exit"] or 0.0)
            pred[label] = (0.5 * (s["T0"] + s["T1"]) + cross, 0.5 * (s["T0"] + s["T1"]) - cross)
        for sign, i in (("plus", 0), ("minus", 1)):
            c, r = pred["candidate"][i], pred["reference"][i]
            if c < tol.t_noise_floor and r < tol.t_noise_floor:
                continue
            ok = abs(c - r) <= max(tol.t_atol, tol.t_rtol * abs(r))
            rows.append({"key": list(key), "state": sign, "candidate": c, "reference": r, "status": "pass" if ok else "fail"})
            if not ok:
                failures.append((key, sign))
    detail = f"{len(rows) - len(failures)} pass, {len(failures)} fail"
    if failures:
        detail += f": {failures}"
    return Check("pair.mixed_transmission", not failures, "fail", detail, {"rows": rows})


def compare_polarization_pair(candidate: tuple[Path, Path], reference: tuple[Path, Path],
                              tol: Tolerances = Tolerances()) -> ComparisonResult:
    """Compare what coherent transport uses and the per-directory comparison cannot see (CDX-005).

    Each argument is (Ephi=0 directory, Ephi=1 directory). The checks are: the exit-field and far-field
    overlap Re rho of the two polarizations (BBRsim's polarization cross term; the far field weighted by
    sin(Theta)), and the transmission of equal-amplitude mixed states, T = (T0+T1)/2 +- sqrt(T0 T1) Re rho. A
    relative sign or phase change between the polarizations fails them; a phase common to both does not.
    Keys where one polarization is cut off are skipped for rho. At keys where rho is near 0 (normal
    incidence, E along one axis) a relative sign is invisible, and harmless, since equal-plus and equal-minus
    transmissions are then equal; the check sees the phase wherever it matters to BBRsim. T is the producer's
    OutgoingPower / IngoingPower, not BBRsim's load-time value divided by the key's top eigenvalue above 1.
    """
    def load(dirs):
        return [(pd.read_csv(Path(d) / "waveguide.csv"), pd.read_csv(Path(d) / "far_field.csv")) for d in dirs]

    cand, ref = load(candidate), load(reference)
    checks = [_check_pair_ephi(cand, ref)]
    (cs, cp), (rs, rp) = _pair_stats(cand, tol), _pair_stats(ref, tol)
    problems = [f"candidate {p}" for p in cp] + [f"reference {p}" for p in rp]
    for (cw, _), (rw, _), e in zip(cand, ref, (0, 1)):
        if _keys(cw) != _keys(rw):
            problems.append(f"Ephi={e} incidence keys differ from the reference: missing {sorted(_keys(rw) - _keys(cw))} "
                            f"extra {sorted(_keys(cw) - _keys(rw))}")
    checks.append(Check("pair.points", not problems, "fail",
                        "both pairs pair by row at every key" if not problems else "; ".join(problems[:5])))
    checks += [_check_rho("pair.exit_overlap", "rho_exit", cs, rs, tol),
               _check_rho("pair.far_field_overlap", "rho_far", cs, rs, tol),
               _check_mixed(cs, rs, tol)]
    return ComparisonResult(f"{candidate[0]} + {candidate[1]}", f"{reference[0]} + {reference[1]}", checks, asdict(tol))
