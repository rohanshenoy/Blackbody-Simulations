#!/usr/bin/env python3
"""Per-key detail of a job compared with the Windows reference: convergence, transmission, self-consistency, fields.

Written for HPC step 5 (job reference1, 2026-10-06), whose comparison failed on transmission only.
Run on the login node after `source /home/rshenoy/BBRSim/bb_env.sh` (pandas, a few hundred MB read):

    python hpc/compare_detail.py /home/rshenoy/BBRSim/outputs/InfParallelPlate_crack1Rohan_500GHz/job_reference1

Prints, for both polarizations:
  1. convergence: passes completed and the last rows of each convergence table, and the solve time;
  2. transmission T = OutgoingPower / IngoingPower per incidence key, candidate against reference, with
     the relative difference;
  3. self-consistency at normal incidence (theta = 180 deg), where physics fixes the answer for this
     doubly mirror-symmetric gap: (0, 180) Ephi=0 and (90, 180) Ephi=1 are the same incident field (E along
     the gap), (45, 180) must be half of it for each polarization, and (0, 180) Ephi=1 and (90, 180) Ephi=0
     (E along the long side, below the TE01 cutoff) must vanish. Each polarization is a separate parametric
     variation with its own adaptive mesh, so a mismatch between them measures mesh noise;
  4. field agreement per key: |E| correlation and RMS difference (normalized by the reference maximum) of
     the exit field and the far field.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

KEYS = ["IWavePhi", "IWaveTheta"]
WG = ["Ex", "Ey", "Ez"]
FF = ["rEtheta", "rEphi"]


def transmission(df: pd.DataFrame) -> pd.Series:
    first = df.drop_duplicates(KEYS).set_index(KEYS)
    return first["OutgoingPower"] / first["IngoingPower"]


def magnitude(df: pd.DataFrame, prefixes: list[str]) -> np.ndarray:
    return np.sqrt(sum(df[f"{p}_real"].to_numpy() ** 2 + df[f"{p}_imag"].to_numpy() ** 2 for p in prefixes))


def field_agreement(cand: pd.DataFrame, ref: pd.DataFrame, coords: list[str], prefixes: list[str]) -> dict:
    out = {}
    for key, r in ref.groupby(KEYS):
        c = cand[(cand[KEYS[0]] == key[0]) & (cand[KEYS[1]] == key[1])]
        if len(c) != len(r):
            out[key] = (np.nan, np.nan)
            continue
        c = c.sort_values(coords)
        r = r.sort_values(coords)
        mc, mr = magnitude(c, prefixes), magnitude(r, prefixes)
        scale = float(np.max(mr)) or 1.0
        corr = float(np.corrcoef(mc, mr)[0, 1]) if np.std(mc) > 0 and np.std(mr) > 0 else np.nan
        out[key] = (corr, float(np.sqrt(np.mean((mc - mr) ** 2)) / scale))
    return out


def convergence_summary(path: Path) -> str:
    if not path.is_file():
        return f"{path.name}: missing"
    text = path.read_text(errors="replace")
    completed = re.search(r"Completed\s*:\s*(\S+)", text)
    rows = [line.strip() for line in text.splitlines() if re.match(r"^\s*\d+\s*\|", line)]
    tail = " ; ".join(rows[-3:]) if rows else "no pass rows"
    return f"passes completed {completed.group(1) if completed else '?'}; last rows: {tail}"


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("job", type=Path, help="job directory written by run_hfss_frequency.py")
    p.add_argument("--stem", default="InfParallelPlate_crack1Rohan_500GHz")
    p.add_argument("--reference-root", type=Path, default=Path("/home/rshenoy/BBRSim/BBRSimulation/data/waveguides"))
    args = p.parse_args(argv)
    job, ref_root, stem = args.job, args.reference_root, args.stem

    print("== 1. convergence and solve time")
    for e in (0, 1):
        print(f"  Ephi={e}: {convergence_summary(job / 'logs' / f'convergence_Ephi{e}.txt')}")
    log = job / "logs" / "run.log"
    if log.is_file():
        for line in log.read_text(errors="replace").splitlines():
            if "solve finished in" in line or "solved variations" in line:
                print("  " + line.split(" INFO ", 1)[-1])

    data = {}
    for e in (0, 1):
        for side, root in (("cand", job), ("ref", ref_root)):
            d = root / f"{stem}_Ephi={e}"
            data[(side, e, "wg")] = pd.read_csv(d / "waveguide.csv")
            data[(side, e, "ff")] = pd.read_csv(d / "far_field.csv")
    t = {(side, e): transmission(data[(side, e, "wg")]) for side in ("cand", "ref") for e in (0, 1)}

    print("== 2. transmission per key (rel = candidate / reference - 1)")
    print(f"  {'phi':>5} {'theta':>6} | {'ref E0':>10} {'cand E0':>10} {'rel E0':>8} | {'ref E1':>10} {'cand E1':>10} {'rel E1':>8}")
    for key in t[("ref", 0)].index:
        cells = []
        for e in (0, 1):
            r, c = t[("ref", e)].get(key, np.nan), t[("cand", e)].get(key, np.nan)
            rel = (c / r - 1) if r and abs(r) > 1e-12 else np.nan
            cells.append(f"{r:10.3e} {c:10.3e} {rel:+8.3%}" if np.isfinite(rel) else f"{r:10.3e} {c:10.3e} {'-':>8}")
        print(f"  {key[0]:5.0f} {key[1]:6.0f} | {cells[0]} | {cells[1]}")

    print("== 3. self-consistency at normal incidence (theta = 180)")
    for side in ("ref", "cand"):
        g = lambda e, phi: float(t[(side, e)].get((float(phi), 180.0), np.nan))  # noqa: E731
        gap0, gap1 = g(0, 0), g(1, 90)
        print(f"  {side:4}: E along gap: (0,180) E0 = {gap0:.5f}, (90,180) E1 = {gap1:.5f}, ratio {gap0 / gap1:.4f}")
        print(f"        (45,180): E0 = {g(0, 45):.5f} (half of E0 at 0: {g(0, 45) / gap0:.4f} of it), "
              f"E1 = {g(1, 45):.5f} (half of E1 at 90: {g(1, 45) / gap1:.4f} of it)")
        print(f"        E along the long side (cut off): (0,180) E1 = {g(1, 0):.3e}, (90,180) E0 = {g(0, 90):.3e}")

    print("== 4. field agreement per key: |E| correlation, normalized RMS difference")
    for e in (0, 1):
        wg = field_agreement(data[("cand", e, "wg")], data[("ref", e, "wg")], ["X", "Y", "Z"], WG)
        ff = field_agreement(data[("cand", e, "ff")], data[("ref", e, "ff")], ["Theta", "Phi"], FF)
        print(f"  Ephi={e}: {'phi':>5} {'theta':>6} {'T ref':>10} | {'wg corr':>8} {'wg nrms':>8} | {'ff corr':>8} {'ff nrms':>8}")
        for key in t[("ref", e)].index:
            wc, wn = wg.get(key, (np.nan, np.nan))
            fc, fn = ff.get(key, (np.nan, np.nan))
            print(f"          {key[0]:5.0f} {key[1]:6.0f} {t[('ref', e)][key]:10.3e} | {wc:8.4f} {wn:8.4f} | {fc:8.4f} {fn:8.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
