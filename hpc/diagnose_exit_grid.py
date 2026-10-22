#!/usr/bin/env python3
"""Exit-lattice diagnostic for the round gap's curved rim (HPC step 7a, after step 7 failed).

Step 7 (job roundgap1, Slurm 4077742, 2026-10-06) solved, then stopped at its first exit-field export: HFSS
had written field values at 50 lattice points outside the declared 50 um disc, and the runner's containment
check refuses any valued point outside the declared cross-section. The runner deletes each raw export once
it is read, so those points are lost.

This script runs the runner unchanged (bbsim.run_frequency.main, with run_hfss_frequency.py's arguments)
with two instruments: every raw exit-field export is copied to <job>/diagnose_exit_grid/ before it is read,
and the containment check, and since the rim rule (ledger HF-029) also that rule, print their verdicts
instead of raising, so the run continues through both polarizations. Each polarization is a separate parametric variation with its own adaptive mesh. Then, for
every kept export, it maps which lattice points HFSS filled and which it left nan against their distance
from the declared rim:

  1. lattice rows; valued points inside and outside the disc; nan points inside it;
  2. how far beyond the rim the valued outside points lie, and how deep inside it the nan inside points lie;
  3. |E| at the valued outside points against the ring just inside the rim: a continuation of the interior
     field, or junk such as zeros;
  4. the twenty lattice points exactly on the rim;
  5. across exports: whether the valued points inside the disc are the same set (Geant4 pairs the two
     polarizations' exit points by row index), which is what keeping only the points inside would retain.

Run on HPC in a batch job after `source /home/rshenoy/BBRSim/bb_env.sh` (uses the HFSS solver licence):

    python hpc/diagnose_exit_grid.py --config configs/round_gap_r50um_2000GHz_single_angle.toml \\
        --polarizations 0 1 --cores 4 --job-id diag_exitgrid1

The job directory it writes keeps any outside points in its CSV tables and is not a dataset. The report
starts at the line "== runner exit status"; exit 0 when every kept export was analysed, 1 otherwise.
"""
from __future__ import annotations

import os
import shutil
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import bbsim.extract as extract  # noqa: E402
from bbsim import run_frequency  # noqa: E402
from bbsim.readers import FieldFileError  # noqa: E402

KEEP_DIR = "diagnose_exit_grid"
BIN_EDGES_NM = [0, 20, 50, 100, 200, 500, 1000, 2000, np.inf]
RING_M = 2e-6          # the ring just inside the rim whose |E| the outside values are compared with
RIM_TOL_M = 1e-12      # a lattice point this close to the declared radius is on the rim

kept: list[Path] = []
declared: list[dict] = []
_exports: Counter = Counter()
_read = extract.read_exit_field_fld
_check = extract.check_points_within
_retain = extract.retain_disc_lattice


def keeping_read(path, **kwargs):
    """Copy the raw export to <job>/diagnose_exit_grid/ (the runner deletes it after reading), then read it."""
    path = Path(path)
    keep = path.parent.parent / KEEP_DIR
    keep.mkdir(exist_ok=True)
    _exports[path.stem] += 1
    dest = keep / f"{path.stem}_{_exports[path.stem]:02d}.fld"
    shutil.copy2(path, dest)
    kept.append(dest)
    return _read(path, **kwargs)


def reporting_retain(df, cross_section, where, **kwargs):
    """The rim rule (ledger HF-029), reported instead of raised: on a refusal the run goes on with no band."""
    declared.append(dict(cross_section))
    try:
        return _retain(df, cross_section, where, **kwargs)
    except FieldFileError as exc:
        print(f"rim rule FAIL (reported, not raised): {exc}", flush=True)
        return _retain(df, cross_section, where, **{**kwargs, "band_rel": float("inf")})


def reporting_check(df, cross_section, where, **kwargs):
    """The runner's containment check, reported instead of raised."""
    declared.append(dict(cross_section))
    try:
        _check(df, cross_section, where, **kwargs)
        print(f"containment PASS: {where}", flush=True)
    except FieldFileError as exc:
        print(f"containment FAIL (reported, not raised): {exc}", flush=True)


@dataclass
class Lattice:
    name: str
    yz_nm: np.ndarray       # (N, 2) integer coordinates in nm, to compare exports point by point
    valued: np.ndarray      # HFSS wrote field values (not nan)
    outside: np.ndarray     # outside the declared disc, by the runner's own criterion
    rim: np.ndarray
    r_m: np.ndarray
    mag: np.ndarray         # |E| from the six field components; nan where not valued

    def points(self, mask: np.ndarray) -> set[tuple[int, int]]:
        return {tuple(p) for p in self.yz_nm[mask].tolist()}


def read_lattice(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Every row of a raw ExportOnGrid file: coordinates (m), whether HFSS wrote values, the six components."""
    coords, valued, fields = [], [], []
    for line in path.read_text().splitlines()[2:]:
        parts = line.split()
        if not parts:
            continue
        coords.append([float(t) for t in parts[:3]])
        tokens = parts[3:]
        is_nan = bool(tokens) and all(t.lower() == "nan" for t in tokens)
        valued.append(not is_nan)
        fields.append([float(t) for t in tokens] if not is_nan and len(tokens) == 6 else [np.nan] * 6)
    return np.array(coords, dtype=float), np.array(valued, dtype=bool), np.array(fields, dtype=float)


def histogram(values_nm: np.ndarray) -> str:
    counts, _ = np.histogram(values_nm, bins=BIN_EDGES_NM)
    labels = [f"{lo:g}-{hi:g}" if np.isfinite(hi) else f"{lo:g}+" for lo, hi in zip(BIN_EDGES_NM, BIN_EDGES_NM[1:])]
    return ", ".join(f"{label} nm: {n}" for label, n in zip(labels, counts) if n)


def analyse(path: Path, radius_m: float, rel_tol: float) -> Lattice:
    xyz, valued, fields = read_lattice(path)
    y, z = xyz[:, 1], xyz[:, 2]
    r = np.hypot(y, z)
    lat = Lattice(
        name=path.name,
        yz_nm=np.rint(xyz[:, 1:3] * 1e9).astype(np.int64),
        valued=valued,
        outside=y * y + z * z > radius_m ** 2 * (1 + rel_tol),
        rim=np.abs(r - radius_m) <= RIM_TOL_M,
        r_m=r,
        mag=np.where(valued, np.sqrt(np.sum(fields ** 2, axis=1)), np.nan),
    )
    vin, vout = valued & ~lat.outside, valued & lat.outside
    nin, nout = ~valued & ~lat.outside, ~valued & lat.outside
    print(f"== {path.name}: {len(r)} lattice rows, {int(valued.sum())} valued, {int((~valued).sum())} nan; "
          f"inside the disc {int((~lat.outside).sum())} (rim {int(lat.rim.sum())}); max |X| {np.max(np.abs(xyz[:, 0])):.1e} m")
    print(f"   valued inside {int(vin.sum())}, valued outside {int(vout.sum())}, nan inside {int(nin.sum())}, "
          f"nan outside {int(nout.sum())}; rim points valued {int((valued & lat.rim).sum())} of {int(lat.rim.sum())}")
    if vout.any():
        beyond_nm = (r[vout] - radius_m) * 1e9
        print(f"   valued outside, distance beyond the rim: max {beyond_nm.max():.1f} nm; {histogram(beyond_nm)}")
        ring = vin & (radius_m - r < RING_M)
        ring_median = float(np.median(lat.mag[ring])) if ring.any() else float("nan")
        ratio = lat.mag[vout] / ring_median
        print(f"   |E| there over the median |E| within {RING_M * 1e6:g} um inside the rim ({ring_median:.4g} V/m): "
              f"median {np.median(ratio):.3f}, min {ratio.min():.3f}, max {ratio.max():.3f}")
        far = np.flatnonzero(vout)[np.argsort(-beyond_nm)][:8]
        print("   farthest (Y um, Z um, nm beyond): " + "; ".join(
            f"({y[i] * 1e6:.3f}, {z[i] * 1e6:.3f}, {(r[i] - radius_m) * 1e9:.1f})" for i in far))
    if nin.any():
        depth_nm = (radius_m - r[nin]) * 1e9
        print(f"   nan inside, depth inside the rim: max {depth_nm.max():.1f} nm; {histogram(depth_nm)}")
        deep = np.flatnonzero(nin)[np.argsort(-depth_nm)][:8]
        print("   deepest (Y um, Z um, nm inside): " + "; ".join(
            f"({y[i] * 1e6:.3f}, {z[i] * 1e6:.3f}, {(radius_m - r[i]) * 1e9:.1f})" for i in deep))
    return lat


def compare(lattices: list[Lattice]) -> bool:
    """True when every export holds the same valued points inside the disc."""
    same = True
    base = lattices[0]
    for other in lattices[1:]:
        a, b = base.points(base.valued & ~base.outside), other.points(other.valued & ~other.outside)
        c, d = base.points(base.valued & base.outside), other.points(other.valued & other.outside)
        lattice_same = base.points(np.ones_like(base.valued)) == other.points(np.ones_like(other.valued))
        print(f"== {base.name} against {other.name}: lattice {'identical' if lattice_same else 'DIFFERS'}; "
              f"valued inside {'identical' if a == b else f'DIFFER ({len(a - b)} only in the first, {len(b - a)} only in the second)'}; "
              f"valued outside {'identical' if c == d else f'differ ({len(c - d)} only in the first, {len(d - c)} only in the second, {len(c & d)} shared)'}")
        same = same and a == b
    return same


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    stale = sorted(k for k in os.environ if k.startswith("PYAEDT_"))
    for k in stale:  # as hpc/run_frequency.sbatch: a stale port or pid would attach PyAEDT to another desktop
        del os.environ[k]
    if stale:
        print(f"removed from the environment: {stale}", flush=True)
    extract.read_exit_field_fld = keeping_read
    extract.retain_disc_lattice = reporting_retain
    extract.check_points_within = reporting_check
    status = run_frequency.main(argv)
    sys.stderr.flush()
    print(f"\n== runner exit status {status} (the job log is in logs/run.log of the job directory)", flush=True)
    for verdict in sorted(set(f"{d}" for d in declared)):
        print(f"   declared cross-section: {verdict}")
    if not kept:
        print("no exit-field export was kept: the run stopped before the first export")
        return 1
    print(f"   kept exports: {kept[0].parent}")
    disc = next((d for d in declared if d.get("shape") == "disc"), None)
    if disc is None:
        print(f"the declared cross-section is not a disc: {declared[:1]}")
        return 1
    lattices = [analyse(p, float(disc["radius_m"]), extract.CONTAIN_REL_TOL) for p in kept]
    same = compare(lattices) if len(lattices) > 1 else True
    parts = []
    for lat in lattices:
        vout = lat.valued & lat.outside
        beyond = f", max {(lat.r_m[vout].max() - float(disc['radius_m'])) * 1e9:.0f} nm beyond" if vout.any() else ""
        parts.append(f"{lat.name}: valued inside {int((lat.valued & ~lat.outside).sum())} of "
                     f"{int((~lat.outside).sum())}, valued outside {int(vout.sum())}{beyond}")
    print("EXIT GRID: " + "; ".join(parts) + f"; inside sets {'identical' if same else 'DIFFER'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
