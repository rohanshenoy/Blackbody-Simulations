#!/usr/bin/env python3
"""Settle how HFSS maps (IWavePhi, IWaveTheta) to the incident propagation direction, without a new solve.

Run on HPC in a batch job after the crack1 full reference run (hpc/README.md step 5), on that job's own
project copy, which this script never saves:
    python hpc/check_incident_direction.py --project <job>/project/ParallelPlateGaps.aedt \
        --design parallel_plate_gap_50um_500GHz --out <job>/incident_direction.json

At (IWavePhi, IWaveTheta) = (90, 135) deg, Ephi = 0, it exports E for TotalFields and for ScatteredFields
on two lines through the gap near the entrance (along global y and along global z), forms the incident
field E_inc = E_total - E_scattered, and fits the phase slope along each line. The ratio k_y / k_z does
not depend on the e^{+jwt} / e^{-jwt} convention:
    arrival-direction reading, k = -r_hat(theta, phi) = (0, -0.707, +0.707):  k_y / k_z = -1
    alternative, k = (sin t cos p, sin t sin p, -cos t) = (0, +0.707, +0.707): k_y / k_z = +1
The ratio cannot tell k from -k. As a self-check outside the verdict, the JSON also reports whether the z
slope is negative, as HFSS's e^{+jwt} (E_inc ~ e^{-jk.r}) requires for a wave entering the gap (k_z > 0 in
both readings above).
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bbsim.extract import intrinsics  # noqa: E402
from bbsim.readers import read_exit_field_fld, wait_for_file  # noqa: E402

PHI_DEG, THETA_DEG, EPHI = 90.0, 135.0, 0
FREQ_GHZ = 500.0
FREQ_LABEL = "500GHz"
SOLUTION = "500GHz : LastAdaptive"
# Lines inside the crack1 gap (x 0..0.05, y -5..5, z 0..1 mm) through the entrance-face centre column.
LINES_MM = {
    "y": ((0.025, -0.1, 0.01), (0.025, 0.1, 0.01), ("0mm", "0.002mm", "0mm")),
    "z": ((0.025, 0.0, 0.01), (0.025, 0.0, 0.21), ("0mm", "0mm", "0.002mm")),
}
EXPECTED_RATIO = -1.0
TOLERANCE = 0.05
SLOPE_Z_NOTE = ("Self-check, not in the verdict: under HFSS's e^{+jwt}, E_inc ~ e^{-jk.r}, so a wave entering the "
                "gap (k_z > 0) needs a negative z slope. The ratio cannot tell k from -k; false means k points out "
                "of the gap or the phasors use e^{-jwt}.")


def phase_slope(positions_m: np.ndarray, field: np.ndarray) -> float:
    """Least-squares slope (rad/m) of the unwrapped phase of a complex field sampled along a line."""
    return float(np.polyfit(np.asarray(positions_m, dtype=float), np.unwrap(np.angle(field)), 1)[0])


def dominant_component(fields: np.ndarray) -> np.ndarray:
    """The Cartesian component (column of an (N, 3) complex array) with the largest mean magnitude."""
    return fields[:, int(np.argmax(np.abs(fields).mean(axis=0)))]


def verdict(ratio: float) -> str:
    if abs(ratio - EXPECTED_RATIO) <= TOLERANCE:
        return "PASS: arrival direction, k = -r_hat(theta, phi)"
    if abs(ratio + EXPECTED_RATIO) <= TOLERANCE:
        return "ALTERNATIVE: k = (sin t cos p, sin t sin p, -cos t)"
    return "INCONCLUSIVE"


def _complex_fields(df) -> np.ndarray:
    return np.column_stack([df[f"E{c}_real"].to_numpy() + 1j * df[f"E{c}_imag"].to_numpy() for c in "xyz"])


def _export_line(hfss, field_type: str, axis: str, scratch: Path):
    lo, hi, step = LINES_MM[axis]
    hfss.odesign.GetModule("Solutions").EditSources(["FieldType:=", field_type])
    calculator = hfss.odesign.GetModule("FieldsReporter")
    calculator.CalcStack("clear")
    calculator.EnterQty("E")
    path = scratch / f"{field_type}_{axis}.fld"
    calculator.ExportOnGrid(
        str(path), [f"{v}mm" for v in lo], [f"{v}mm" for v in hi], list(step), SOLUTION,
        intrinsics(EPHI, FREQ_LABEL, PHI_DEG, THETA_DEG) + ["Phase:=", "0deg"],
        ["NAME:ExportOption", "IncludePtInOutput:=", True, "RefCSName:=", "Global",
         "PtInSI:=", True, "FieldInRefCS:=", False],
        "Cartesian", ["0mm", "0mm", "0mm"], False)
    wait_for_file(path, 300.0)
    df = read_exit_field_fld(path)
    path.unlink()
    return df


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--project", type=Path, required=True, help="The solved step-5 job's project copy")
    p.add_argument("--design", required=True, help="The solved run design, e.g. parallel_plate_gap_50um_500GHz")
    p.add_argument("--out", type=Path, required=True, help="Where to write the JSON result")
    p.add_argument("--aedt-version", default="2025.2")
    args = p.parse_args(argv)

    from bbsim.session import aedt_session  # lazy PyAEDT import

    slopes: dict[str, float] = {}
    with tempfile.TemporaryDirectory(prefix="incident_check_") as tmp:
        scratch = Path(tmp)
        with aedt_session(args.project, args.design, args.aedt_version) as hfss:
            for axis in ("y", "z"):
                total = _export_line(hfss, "TotalFields", axis, scratch)
                scattered = _export_line(hfss, "ScatteredFields", axis, scratch)
                incident = _complex_fields(total) - _complex_fields(scattered)
                positions = total["Y" if axis == "y" else "Z"].to_numpy()
                slopes[axis] = phase_slope(positions, dominant_component(incident))
            hfss.odesign.GetModule("Solutions").EditSources(["FieldType:=", "TotalFields"])
    k0 = 2 * math.pi * FREQ_GHZ * 1e9 / 299792458.0
    ratio = slopes["y"] / slopes["z"]
    result = {
        "iwave_phi_deg": PHI_DEG, "iwave_theta_deg": THETA_DEG, "ephi": EPHI, "frequency_ghz": FREQ_GHZ,
        "slope_y_rad_per_m": slopes["y"], "slope_z_rad_per_m": slopes["z"], "k0_rad_per_m": k0,
        "ratio_ky_over_kz": ratio, "magnitude_over_k0": math.hypot(slopes["y"], slopes["z"]) / k0,
        "expected_ratio": EXPECTED_RATIO, "verdict": verdict(ratio),
        "slope_z_negative_as_expected": bool(slopes["z"] < 0), "slope_z_note": SLOPE_Z_NOTE,
    }
    text = json.dumps(result, indent=2)
    print(text)
    args.out.write_text(text + "\n")
    return 0 if result["verdict"].startswith("PASS") else 1


if __name__ == "__main__":
    raise SystemExit(main())
