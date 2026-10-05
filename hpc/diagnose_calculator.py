#!/usr/bin/env python3
"""Diagnose a failed field-calculator evaluation (ClcEval) on a solved job, without a new solve.

Run on HPC in a batch job, after source /home/rshenoy/BBRSim/bb_env.sh:
    python hpc/diagnose_calculator.py --job /home/rshenoy/BBRSim/outputs/InfParallelPlate_crack1Rohan_500GHz/job_baseline1

It copies the job's project file and its results folder into <job>/diagnose_<utc>/, so the job's own
copy is never opened. It opens the copy headlessly and prints what AEDT holds: designs, setups, solved
variations, design variables, the plane-wave definition. Then it tries the outgoing-power evaluation the
runner makes, and several variants of its arguments, printing AEDT's own messages after each attempt.
It also writes everything to diagnose.json in the copy folder.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bbsim.extract import intrinsics  # noqa: E402
from bbsim.naming import frequency_design_name, frequency_label, last_adaptive_solution  # noqa: E402


def variants(ephi: int, freq: str, phi: float, theta: float) -> list[tuple[str, list]]:
    return [
        ("as the runner (Ephi int, angles x.0deg)", intrinsics(ephi, freq, phi, theta)),
        ("Ephi as a string", ["Ephi:=", str(ephi), "Freq:=", freq, "IWavePhi:=", f"{phi}deg", "IWaveTheta:=", f"{theta}deg"]),
        ("whole-degree angles", ["Ephi:=", str(ephi), "Freq:=", freq, "IWavePhi:=", f"{int(phi)}deg", "IWaveTheta:=", f"{int(theta)}deg"]),
        ("runner arguments plus Phase", intrinsics(ephi, freq, phi, theta) + ["Phase:=", "0deg"]),
        ("no incident angles", ["Ephi:=", str(ephi), "Freq:=", freq]),
        ("Freq only", ["Freq:=", freq]),
    ]


def messages(hfss, seen: set) -> list[str]:
    """AEDT message-manager entries not printed before (all severities)."""
    try:
        new = [str(m) for m in hfss.odesktop.GetMessages(hfss.project_name, hfss.design_name, 0)]
    except Exception as exc:  # noqa: BLE001
        return [f"(GetMessages failed: {type(exc).__name__}: {exc})"]
    out = [m for m in new if m not in seen]
    seen.update(out)
    return out


def attempt(label: str, fn, report: dict, key: str):
    try:
        value = fn()
        report[key] = {"ok": True, "value": repr(value)}
        print(f"  {label}: OK -> {value!r}")
        return value
    except Exception as exc:  # noqa: BLE001 - this is a diagnostic: record and go on
        report[key] = {"ok": False, "error": f"{type(exc).__name__}: {exc}"}
        print(f"  {label}: FAILED -> {type(exc).__name__}: {exc}")
        return None


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--job", type=Path, required=True, help="Job directory written by run_hfss_frequency.py")
    p.add_argument("--aedt-version", default="2025.2")
    args = p.parse_args(argv)

    job = args.job.resolve()
    cfg = json.loads((job / "effective_config.json").read_text())
    freq = float(cfg["solver"]["frequency_ghz"])
    design = frequency_design_name(cfg["project"]["design"], freq)
    label = frequency_label(freq)
    solution = last_adaptive_solution(freq)
    phi = float(cfg["excitation"]["phi_lower_deg"])
    theta = float(cfg["excitation"]["theta_lower_deg"])

    projects = sorted((job / "project").glob("*.aedt"))
    if len(projects) != 1:
        print(f"expected one .aedt in {job / 'project'}, found {[q.name for q in projects]}")
        return 2
    src = projects[0]
    results = src.with_name(src.name + "results")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    dest = job / f"diagnose_{stamp}"
    dest.mkdir()
    shutil.copy2(src, dest / src.name)
    if results.is_dir():
        shutil.copytree(results, dest / results.name)
    print(f"copied {src.name}{' and ' + results.name if results.is_dir() else ' (no results folder found)'} to {dest}")
    print(f"design {design}, solution {solution}, Freq {label}, IWavePhi {phi}, IWaveTheta {theta}")

    from bbsim.session import aedt_session, aedt_versions  # lazy PyAEDT import

    report: dict = {"job": str(job), "copy": str(dest), "design": design, "solution": solution}
    seen: set = set()
    with aedt_session(dest / src.name, design, args.aedt_version) as hfss:
        report["versions"] = aedt_versions(hfss)
        print("versions:", report["versions"])
        print("== what AEDT holds")
        attempt("designs", lambda: list(hfss.design_list), report, "designs")
        attempt("setups", lambda: list(hfss.odesign.GetModule("AnalysisSetup").GetSetups()), report, "setups")
        attempt("solutions", lambda: list(hfss.existing_analysis_sweeps), report, "solutions")
        attempt("solved variations", lambda: list(hfss.odesign.ListVariations(solution)), report, "variations")
        attempt("design variables", lambda: {str(v): str(hfss.odesign.GetVariableValue(v))
                                             for v in hfss.odesign.GetVariables()}, report, "variables")
        attempt("excitations", lambda: list(hfss.odesign.GetModule("BoundarySetup").GetExcitations()), report, "excitations")
        attempt("plane wave definition", lambda: {k: str(v) for b in hfss.boundaries if "Incident" in str(b.type)
                                                 for k, v in dict(b.props).items()}, report, "plane_wave")
        for m in messages(hfss, seen):
            print("    aedt:", m)

        fields = hfss.odesign.GetModule("FieldsReporter")
        print("== control: a constant on the calculator stack")
        def control():
            fields.CalcStack("clear")
            fields.EnterScalar(1)
            fields.ClcEval(solution, ["Freq:=", label], "Fields")
            value = fields.GetTopEntryValue(solution, ["Freq:=", label])
            fields.CalcStack("pop")
            return value
        attempt("EnterScalar(1) evaluated at Freq only", control, report, "control")
        for m in messages(hfss, seen):
            print("    aedt:", m)

        print("== outgoing power, as the runner evaluates it, then variants")
        report["variants"] = {}
        for name, var in variants(0, label, phi, theta):
            def evaluate(var=var):
                fields.CalcStack("clear")
                fields.CopyNamedExprToStack("outgoing_power")
                fields.ClcEval(solution, var, "Fields")
                value = fields.GetTopEntryValue(solution, var)
                fields.CalcStack("pop")
                return value
            print(f"variant: {name}: {var}")
            attempt(name, evaluate, report["variants"], name)
            for m in messages(hfss, seen):
                print("    aedt:", m)

    (dest / "diagnose.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"report written to {dest / 'diagnose.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
