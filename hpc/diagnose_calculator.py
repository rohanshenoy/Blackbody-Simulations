#!/usr/bin/env python3
"""Extraction diagnostic for a solved job whose run failed after the solve (no new solve).

Written after HPC step 4, job baseline1 (2026-10-05): the solve finished in 66.6 s, then the first
field-calculator read-back (ClcEval over gRPC) failed. Run on HPC in a batch job, after
`source /home/rshenoy/BBRSim/bb_env.sh`:

    python hpc/diagnose_calculator.py --job /home/rshenoy/BBRSim/outputs/InfParallelPlate_crack1Rohan_500GHz/job_baseline1

Parts, each reported on its own so one failure does not hide the rest:
  1. the job's AEDT results folder: files, sizes, modification times, tails of solver logs (no AEDT);
  2. a copy of the project and its results in <job>/diagnose_<utc>/, opened headlessly, and what AEDT
     holds: designs, setups, solutions, solved variations, variables, excitations, the plane wave, the
     setup and parametric-sweep properties, the convergence table, profile and mesh statistics of both
     polarizations;
  3. the runner's own extraction on the copy: exit face, outgoing power through CalculatorWrite (the
     fix for the ClcEval failure), exit-field export, far-field export. The CSV tables are written into
     the copy folder and the step 4 Expect values are printed (rows per table, T = outgoing/ingoing
     power). If the power evaluation fails, raw CalculatorWrite variants and a bare ExportOnGrid are
     tried so the file formats and AEDT's messages are on record;
  4. for the record, the legacy read-back (ClcEval then GetTopEntryValue), expected to fail over gRPC.
Everything is also written to diagnose.json in the copy folder. Exit 0 when the extraction passed,
1 when it did not, 2 when the job directory is unusable.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bbsim.config import build_config  # noqa: E402
from bbsim.extract import intrinsics  # noqa: E402
from bbsim.naming import (  # noqa: E402
    dataset_dir_name,
    frequency_design_name,
    frequency_label,
    last_adaptive_solution,
    radiation_sphere_name,
    setup_name,
)
from bbsim.readers import read_exit_field_fld, wait_for_file  # noqa: E402
from bbsim.sampling import far_field_grid, incident_angle_values, incoming_power_w  # noqa: E402

TEXT_SUFFIXES = {".log", ".txt", ".prof", ".err", ".out"}
MAX_TEXT_BYTES = 200_000
PROBE_WAIT_S = 60.0


def utc_stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def _jsonable(value):
    try:
        json.dumps(value)
        return value
    except TypeError:
        return repr(value)


def list_results(results: Path, limit: int = 80) -> dict:
    """Files, sizes, mtimes and log tails of an .aedtresults folder: evidence of the solve without AEDT."""
    if not results.is_dir():
        print(f"no results folder at {results}")
        return {"present": False}
    files = sorted(p for p in results.rglob("*") if p.is_file())
    entries = []
    for p in files:
        st = p.stat()
        entries.append({"path": str(p.relative_to(results)), "bytes": st.st_size,
                        "mtime_utc": datetime.fromtimestamp(st.st_mtime, timezone.utc).isoformat(timespec="seconds")})
    total = sum(e["bytes"] for e in entries)
    print(f"results folder {results}: {len(entries)} files, {total / 1e6:.1f} MB")
    for e in entries[:limit]:
        print(f"  {e['bytes']:>12}  {e['mtime_utc']}  {e['path']}")
    if len(entries) > limit:
        print(f"  ... {len(entries) - limit} more files")
    tails: dict[str, list[str]] = {}
    for p in files:
        if p.suffix.lower() in TEXT_SUFFIXES and p.stat().st_size <= MAX_TEXT_BYTES:
            tail = p.read_text(errors="replace").splitlines()[-25:]
            tails[str(p.relative_to(results))] = tail
            print(f"-- tail of {p.relative_to(results)}")
            for line in tail:
                print("     " + line)
    return {"present": True, "files": entries, "total_bytes": total, "log_tails": tails}


def messages(hfss, seen: set) -> list[str]:
    """AEDT message-manager entries not printed before (all severities)."""
    try:
        new = [str(m) for m in hfss.odesktop.GetMessages(hfss.project_name, hfss.design_name, 0)]
    except Exception as exc:  # noqa: BLE001
        return [f"(GetMessages failed: {type(exc).__name__}: {exc})"]
    out = [m for m in new if m not in seen]
    seen.update(out)
    return out


def attempt(report: dict, key: str, label: str, fn, hfss=None, seen=None):
    """Run one probe; record and print its result or its exception, then AEDT's new messages."""
    try:
        value = fn()
    except Exception as exc:  # noqa: BLE001 - this is a diagnostic: record and go on
        report[key] = {"ok": False, "error": f"{type(exc).__name__}: {exc}"}
        print(f"  {label}: FAILED -> {type(exc).__name__}: {exc}")
        value = None
    else:
        shape = getattr(value, "shape", None)
        shown = f"<table {shape[0]} rows x {shape[1]} columns>" if shape is not None and len(shape) == 2 else repr(value)
        report[key] = {"ok": True, "value": shown if shape is not None else _jsonable(value)}
        print(f"  {label}: OK -> {shown if len(shown) <= 600 else shown[:600] + ' ...'}")
    if hfss is not None and seen is not None:
        for m in messages(hfss, seen):
            print("    aedt:", m)
    return value


def file_tail(path: Path, n: int = 12) -> list[str]:
    if not Path(path).is_file():
        return ["(no file written)"]
    lines = Path(path).read_text(errors="replace").splitlines()
    return (["..."] if len(lines) > n else []) + lines[-n:]


def raw_calculator_probes(hfss, ctx, ephi: int, dest: Path, report: dict, seen: set) -> None:
    """CalculatorWrite with other variation lists; each file is printed so the format is on record."""
    fields = hfss.odesign.GetModule("FieldsReporter")
    base = intrinsics(ephi, ctx.freq_label, float(ctx.phi_values[0]), float(ctx.theta_values[0]))
    variants = {
        "runner arguments without Phase": base,
        "Ephi as a string, with Phase": ["Ephi:=", str(ephi)] + base[2:] + ["Phase:=", "0deg"],
        "no incident angles": ["Ephi:=", str(ephi), "Freq:=", ctx.freq_label, "Phase:=", "0deg"],
        "Freq and Phase only": ["Freq:=", ctx.freq_label, "Phase:=", "0deg"],
    }
    for n, (name, var) in enumerate(variants.items()):
        path = dest / f"probe_power_Ephi{ephi}_{n}.fld"

        def write(var=var, path=path):
            if path.exists():
                path.unlink()
            fields.CalcStack("clear")
            fields.CopyNamedExprToStack("outgoing_power")
            fields.CalculatorWrite(str(path), ["Solution:=", ctx.solution], var)
            wait_for_file(path, PROBE_WAIT_S)
            return file_tail(path)

        print(f"  probe {name}: {var}")
        attempt(report, f"probe_power_{ephi}_{n}", name, write, hfss, seen)
    try:
        fields.CalcStack("clear")
    except Exception:  # noqa: BLE001
        pass


def export_on_grid_probe(hfss, ctx, ephi: int, lattice, dest: Path, report: dict, seen: set) -> None:
    """The exit-field export alone, without any calculator read-back: do field data exist for this variation?"""
    fields = hfss.odesign.GetModule("FieldsReporter")
    range_min, range_max = lattice.range_strings(hfss.modeler.model_units)
    resolution = [f"{v}mm" for v in ctx.exit_field.resolution_mm]
    path = dest / f"probe_exitfield_Ephi{ephi}.fld"
    args = intrinsics(ephi, ctx.freq_label, float(ctx.phi_values[0]), float(ctx.theta_values[0])) + ["Phase:=", "0deg"]

    def export():
        if path.exists():
            path.unlink()
        fields.CalcStack("clear")
        fields.EnterQty("E")
        fields.CalcOp("Smooth")
        fields.ExportOnGrid(
            str(path), range_min, range_max, resolution, ctx.solution, args,
            ["NAME:ExportOption", "IncludePtInOutput:=", True, "RefCSName:=", ctx.exit_cs_name,
             "PtInSI:=", True, "FieldInRefCS:=", False],
            "Cartesian", ["0mm", "0mm", "0mm"], False,
        )
        wait_for_file(path, ctx.timeout_s)
        df = read_exit_field_fld(path)
        return f"{len(df)} points, {len(df.columns)} columns (kept: {path.name})"

    attempt(report, f"probe_exitfield_{ephi}", f"ExportOnGrid alone, Ephi={ephi}", export, hfss, seen)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--job", type=Path, required=True, help="Job directory written by run_hfss_frequency.py")
    p.add_argument("--aedt-version", default="2025.2")
    args = p.parse_args(argv)

    job = args.job.resolve()
    cfg_path = job / "effective_config.json"
    if not cfg_path.is_file():
        print(f"no effective_config.json in {job}")
        return 2
    cfg = build_config(json.loads(cfg_path.read_text()))
    freq = cfg.solver.frequency_ghz
    label = frequency_label(freq)
    setup = setup_name(freq)
    solution = last_adaptive_solution(freq)
    design = frequency_design_name(cfg.project.design, freq)
    projects = sorted((job / "project").glob("*.aedt"))
    if len(projects) != 1:
        print(f"expected one .aedt in {job / 'project'}, found {[q.name for q in projects]}")
        return 2
    src = projects[0]
    results = src.with_name(src.name + "results")

    report: dict = {"job": str(job), "design": design, "setup": setup, "solution": solution, "frequency_label": label}
    print("== 1. results folder of the job (no AEDT)")
    report["results_folder"] = list_results(results)

    dest = job / f"diagnose_{utc_stamp()}"
    dest.mkdir()
    scratch = dest / "scratch"
    scratch.mkdir()
    shutil.copy2(src, dest / src.name)
    if results.is_dir():
        shutil.copytree(results, dest / results.name)
    report["copy"] = str(dest)
    print(f"copied {src.name}{' and ' + results.name if results.is_dir() else ' (no results folder found)'} to {dest}")
    print(f"design {design}, setup {setup}, solution {solution}")

    from bbsim.extract import (  # lazy with PyAEDT below; these need pandas and numpy only
        ExtractionContext,
        describe_exit_face,
        evaluate_outgoing_power,
        exit_grid,
        extract_far_field,
        extract_waveguide,
    )
    from bbsim.geometry import select_faces
    from bbsim.hfss_setup import EXIT_CS, export_convergence_text
    from bbsim.session import aedt_session, aedt_versions  # lazy PyAEDT import

    seen: set = set()
    solved_variations = None
    extraction_ok = False
    try:
        with aedt_session(dest / src.name, design, args.aedt_version) as hfss:
            def probe(key, label_, fn):
                return attempt(report, key, label_, fn, hfss, seen)

            print("== 2. what AEDT holds")
            report["versions"] = aedt_versions(hfss)
            print("versions:", report["versions"])
            probe("designs", "designs", lambda: list(hfss.design_list))
            probe("setups", "setups", lambda: list(hfss.odesign.GetModule("AnalysisSetup").GetSetups()))
            probe("solutions", "solutions", lambda: list(hfss.existing_analysis_sweeps))
            solved_variations = probe("solved_variations", f"solved variations of {solution} (Solutions.GetAvailableVariations)",
                                      lambda: [str(v) for v in hfss.osolution.GetAvailableVariations(solution)])
            probe("available_variations", "PyAEDT available_variations.variations",
                  lambda: hfss.available_variations.variations(solution))
            probe("simulations_running", "simulations running", lambda: hfss.are_there_simulations_running)
            probe("variables", "design variables",
                  lambda: {str(v): str(hfss.odesign.GetVariableValue(v)) for v in hfss.odesign.GetVariables()})
            probe("excitations", "excitations", lambda: list(hfss.odesign.GetModule("BoundarySetup").GetExcitations()))
            probe("plane_wave", "plane wave definition",
                  lambda: {k: str(v) for b in hfss.boundaries if "Incident" in str(b.type) for k, v in dict(b.props).items()})
            probe("setup_props", f"setup {setup} properties",
                  lambda: {s.name: {k: str(v) for k, v in dict(s.props).items()
                                    if k in ("Frequency", "MaxDeltaE", "MaximumPasses", "MinimumPasses",
                                             "MinimumConvergedPasses", "SaveRadFieldsOnly", "IsEnabled")}
                           for s in hfss.setups})
            probe("parametric_props", "parametric sweep properties (SaveFields lives in ProdOptiSetupDataV2)",
                  lambda: {s.name: {"IsEnabled": str(dict(s.props).get("IsEnabled")),
                                    "ProdOptiSetupDataV2": str(dict(s.props).get("ProdOptiSetupDataV2")),
                                    "Sim. Setups": str(dict(s.props).get("Sim. Setups"))}
                           for s in hfss.parametrics.setups})
            for ephi in (0, 1):
                variation = f"Ephi='{ephi}'"
                text = probe(f"convergence_{ephi}", f"convergence table {variation}",
                             lambda v=variation, e=ephi: export_convergence_text(hfss, setup, v, dest / f"convergence_Ephi{e}.txt"))
                if text:
                    for line in text.splitlines()[-6:]:
                        print("     " + line)
                probe(f"profile_{ephi}", f"solver profile {variation}",
                      lambda v=variation, e=ephi: file_tail(Path(hfss.export_profile(setup, v, str(dest / f"profile_Ephi{e}.prof")))))
                probe(f"mesh_stats_{ephi}", f"mesh statistics {variation}",
                      lambda v=variation, e=ephi: file_tail(Path(hfss.export_mesh_stats(setup, v, str(dest / f"mesh_Ephi{e}.ms")))))

            print("== 3. the runner's extraction on the solved copy")
            exc = cfg.excitation
            ff = cfg.far_field
            phi_values = incident_angle_values(exc.phi_lower_deg, exc.phi_upper_deg, exc.phi_step_deg)
            theta_values = incident_angle_values(exc.theta_lower_deg, exc.theta_upper_deg, exc.theta_step_deg)
            faces = probe("faces", "entrance and exit faces", lambda: select_faces(hfss, cfg.geometry))
            if faces is None:
                print("  extraction not run: the faces could not be selected")
            else:
                entrance, exit_face = faces
                grid = far_field_grid(ff.theta_lower_deg, ff.theta_upper_deg, ff.phi_lower_deg, ff.phi_upper_deg,
                                      ff.a_mm, ff.b_mm, freq, ff.fineness, ff.min_coarseness_deg, ff.max_coarseness_deg)
                ctx = ExtractionContext(
                    frequency_ghz=freq, solution=solution, sphere_name=radiation_sphere_name(freq),
                    phi_values=phi_values, theta_values=theta_values, grid=grid, exit_field=cfg.exit_field,
                    exit_face=exit_face, exit_cs_name=EXIT_CS, exit_cs_x=cfg.geometry.exit_cs_x, exit_cs_y=cfg.geometry.exit_cs_y,
                    scratch_dir=scratch, timeout_s=cfg.output.export_timeout_s,
                    incoming_power_w=incoming_power_w(exc.ei_v_per_m, entrance.area_mm2),
                )
                report["incoming_power_w"] = ctx.incoming_power_w
                face_geometry = probe("exit_face", "exit-face description", lambda: describe_exit_face(hfss, ctx))
                lattice = exit_grid(face_geometry, cfg.exit_field) if face_geometry is not None and cfg.exit_field.manual else None
                n_angles = len(phi_values) * len(theta_values)
                expect_wg = lattice.lattice_points * n_angles if lattice is not None else None
                expect_ff = grid.points_per_angle * n_angles
                print(f"  expected rows: waveguide {expect_wg} at most (lattice x {n_angles} angles), far_field {expect_ff}")
                extraction_ok = face_geometry is not None
                for ephi in cfg.output.polarizations:
                    power = probe(f"power_{ephi}", f"outgoing power Ephi={ephi} (CalculatorWrite)",
                                  lambda e=ephi: evaluate_outgoing_power(hfss, ctx, e))
                    if power is not None:
                        report[f"transmission_{ephi}"] = {}
                        for (phi, theta), value in power.items():
                            t = value / ctx.incoming_power_w
                            report[f"transmission_{ephi}"][f"phi={phi} theta={theta}"] = t
                            print(f"    Ephi={ephi} phi={phi} theta={theta}: P_out={value:.6e} W  T={t:.6f}")
                    else:
                        extraction_ok = False
                        raw_calculator_probes(hfss, ctx, ephi, dest, report, seen)
                        if lattice is not None:
                            export_on_grid_probe(hfss, ctx, ephi, lattice, dest, report, seen)
                    out_dir = dest / dataset_dir_name(cfg.project.dataset_id, freq, ephi)
                    out_dir.mkdir(exist_ok=True)
                    if face_geometry is not None:
                        wg = probe(f"waveguide_{ephi}", f"waveguide table Ephi={ephi} (power + ExportOnGrid)",
                                   lambda e=ephi: extract_waveguide(hfss, ctx, e, face_geometry))
                        if wg is not None:
                            wg.to_csv(out_dir / "waveguide.csv", index=False)
                            print(f"    waveguide.csv: {len(wg)} rows (expected {expect_wg} at most) -> {out_dir}")
                            report[f"waveguide_rows_{ephi}"] = len(wg)
                        else:
                            extraction_ok = False
                    far = probe(f"far_field_{ephi}", f"far-field table Ephi={ephi} (ExportFieldsToFile)",
                                lambda e=ephi: extract_far_field(hfss, ctx, e))
                    if far is not None:
                        far.to_csv(out_dir / "far_field.csv", index=False)
                        print(f"    far_field.csv: {len(far)} rows (expected {expect_ff}) -> {out_dir}")
                        report[f"far_field_rows_{ephi}"] = len(far)
                    else:
                        extraction_ok = False

            print("== 4. legacy read-back, for the record (expected to fail over gRPC)")
            fields = hfss.odesign.GetModule("FieldsReporter")
            legacy_args = intrinsics(0, label, float(phi_values[0]), float(theta_values[0]))

            def legacy():
                fields.CalcStack("clear")
                fields.CopyNamedExprToStack("outgoing_power")
                fields.ClcEval(solution, legacy_args, "Fields")
                return fields.GetTopEntryValue(solution, legacy_args)

            probe("legacy_readback", "ClcEval then GetTopEntryValue", legacy)
            try:
                fields.CalcStack("clear")
            except Exception:  # noqa: BLE001
                pass
    finally:
        report["extraction_ok"] = extraction_ok
        (dest / "diagnose.json").write_text(json.dumps(report, indent=2, default=str) + "\n")
        print(f"report written to {dest / 'diagnose.json'}")

    fields_present = bool(solved_variations) or any(report.get(f"power_{e}", {}).get("ok") for e in (0, 1))
    print(f"FIELDS: {'present' if fields_present else 'not confirmed'} (solved variations: {solved_variations})")
    print(f"EXTRACTION {'PASS' if extraction_ok else 'FAIL'}: tables in {dest}")
    return 0 if extraction_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
