#!/usr/bin/env python3
"""Solve-and-probe diagnostic for a job whose extraction failed after its solve (HPC step 4a).

History: baseline1 (2026-10-05) and baseline2 (2026-10-06) both solved: AEDT reported "Normal
completion of simulation" and listed both polarizations as solved variations. Both then failed at
the first evaluation of the named expression ``outgoing_power``: ClcEval in baseline1, CalculatorWrite
in baseline2, each as a bare "Failed to execute gRPC AEDT command". Their job copies hold no solution,
because the runner then saved the project only before the solve and the failed run closed it unsaved.

Run on HPC in a batch job, after `source /home/rshenoy/BBRSim/bb_env.sh`:

    python hpc/diagnose_calculator.py --job <job directory> [--solve]

Parts, each reported on its own so one failure does not hide the rest:
  1. the job's results folder (no AEDT);
  2. a copy of the job's project and results in <job>/diagnose_<utc>/, opened headlessly; what AEDT holds;
  2b. with --solve: the copy is solved (HFSS solver licence), its solved variations listed, its results
      folder listed before and after a save (does AEDT write solution data only on save?);
  3. the calculator's inputs: the named expressions and the exit face list exist;
  4. a matrix of calculator evaluations, each written with CalculatorWrite and read back: a constant,
     then ``outgoing_power`` with the runner's arguments and with variants (Ephi as text, whole-degree
     angles, without the incident angles, PyAEDT's own evaluate), each with T = value / ingoing power;
  5. the field exports alone: the exit-field ExportOnGrid and the far-field ExportFieldsToFile, with
     the runner's arguments and with variants;
  6. the runner's own extraction on the copy, writing the CSV tables beside it;
  7. the legacy read-back (ClcEval then GetTopEntryValue), for the record;
  8. solver profile and mesh statistics, last, so a failing export cannot pre-empt the rest.
After every AEDT call AEDT's new messages are printed. Everything is written to diagnose.json in the
copy folder, and a summary of all probes closes the output. Exit 0 when the runner's extraction
passed, 1 when it did not, 2 when the job directory is unusable.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bbsim.config import build_config  # noqa: E402
from bbsim.extract import intrinsic_variation_key, intrinsics  # noqa: E402
from bbsim.naming import (  # noqa: E402
    dataset_dir_name,
    frequency_design_name,
    frequency_label,
    last_adaptive_solution,
    radiation_sphere_name,
    setup_name,
)
from bbsim.readers import read_calculator_scalar, read_exit_field_fld, read_far_field_ffd, wait_for_file  # noqa: E402
from bbsim.sampling import far_field_grid, incident_angle_values, incoming_power_w  # noqa: E402

TEXT_SUFFIXES = {".log", ".txt", ".prof", ".err", ".out", ".profile", ".asol"}
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


def list_results(results: Path, limit: int = 60, tails: bool = True) -> dict:
    """Files, sizes, mtimes and log tails of an .aedtresults folder: evidence of a solve without AEDT."""
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
    log_tails: dict[str, list[str]] = {}
    if tails:
        for p in files:
            if p.suffix.lower() in TEXT_SUFFIXES and p.stat().st_size <= MAX_TEXT_BYTES:
                tail = p.read_text(errors="replace").splitlines()[-40:]
                log_tails[str(p.relative_to(results))] = tail
                print(f"-- tail of {p.relative_to(results)}")
                for line in tail:
                    print("     " + line)
    return {"present": True, "files": entries, "total_bytes": total, "log_tails": log_tails}


def file_head(path: Path, n: int = 6) -> list[str]:
    if not Path(path).is_file():
        return ["(no file written)"]
    lines = Path(path).read_text(errors="replace").splitlines()
    return lines[:n] + (["..."] if len(lines) > n else [])


def file_tail(path: Path, n: int = 12) -> list[str]:
    if not Path(path).is_file():
        return ["(no file written)"]
    lines = Path(path).read_text(errors="replace").splitlines()
    return (["..."] if len(lines) > n else []) + lines[-n:]


class Probes:
    """Runs probes, records each result, and prints AEDT's new messages after every one."""

    def __init__(self) -> None:
        self.report: dict = {}
        self.order: list[str] = []
        self.hfss = None
        self.seen: set[str] = set()

    def new_messages(self) -> list[str]:
        from bbsim.hfss_setup import aedt_messages

        if self.hfss is None:
            return []
        out = [m for m in aedt_messages(self.hfss) if m not in self.seen]
        self.seen.update(out)
        return out

    def run(self, key: str, label: str, fn):
        try:
            value = fn()
        except Exception as exc:  # noqa: BLE001 - a diagnostic records and goes on
            self.report[key] = {"label": label, "ok": False, "error": f"{type(exc).__name__}: {exc}"}
            print(f"  {label}: FAILED -> {type(exc).__name__}: {exc}")
            value = None
        else:
            shape = getattr(value, "shape", None)
            shown = f"<table {shape[0]} rows x {shape[1]} columns>" if shape is not None and len(shape) == 2 else repr(value)
            self.report[key] = {"label": label, "ok": True, "value": shown if shape is not None else _jsonable(value)}
            print(f"  {label}: OK -> {shown if len(shown) <= 600 else shown[:600] + ' ...'}")
        self.order.append(key)
        for m in self.new_messages():
            print("    aedt:", m)
        return value

    def summary(self) -> None:
        print("== summary of probes")
        for key in self.order:
            r = self.report[key]
            detail = r.get("value") if r["ok"] else r.get("error")
            text = str(detail)
            print(f"  {'OK    ' if r['ok'] else 'FAILED'}  {r['label']}: {text if len(text) <= 160 else text[:160] + ' ...'}")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--job", type=Path, required=True, help="Job directory written by run_hfss_frequency.py")
    p.add_argument("--solve", action="store_true",
                   help="solve the copy before probing (uses the HFSS solver licence); needed when the job's "
                        "project holds no saved solution")
    p.add_argument("--key", nargs=2, type=float, metavar=("PHI", "THETA"), default=None,
                   help="incidence key for the probes, in degrees (default: the config's first key)")
    p.add_argument("--out", default=None,
                   help="name of the copy folder inside the job directory (default: diagnose_<utc>); must not exist")
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

    probes = Probes()
    report = probes.report
    report["meta"] = {"job": str(job), "design": design, "setup": setup, "solution": solution,
                      "frequency_label": label, "solve": args.solve}
    print("== 1. results folder of the job (no AEDT)")
    report["meta"]["job_results_folder"] = list_results(results)

    if args.out:
        dest = job / args.out
        if dest.exists():
            print(f"{dest} exists; choose another --out")
            return 2
    else:
        dest = job / f"diagnose_{utc_stamp()}"
        n = 0
        while dest.exists():  # a second diagnostic within the same second
            n += 1
            dest = job / f"diagnose_{utc_stamp()}_{n}"
    dest.mkdir()
    scratch = dest / "scratch"
    scratch.mkdir()
    shutil.copy2(src, dest / src.name)
    if results.is_dir():
        shutil.copytree(results, dest / results.name)
    copy_results = dest / results.name
    report["meta"]["copy"] = str(dest)
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
    from bbsim.hfss_setup import EXIT_CS, EXIT_FACE_LIST, OUTGOING_POWER_EXPRESSION, export_convergence_text, solve, solved_variations
    from bbsim.session import aedt_session, aedt_versions  # lazy PyAEDT import

    exc = cfg.excitation
    ff = cfg.far_field
    phi_values = incident_angle_values(exc.phi_lower_deg, exc.phi_upper_deg, exc.phi_step_deg)
    theta_values = incident_angle_values(exc.theta_lower_deg, exc.theta_upper_deg, exc.theta_step_deg)
    phi0, theta0 = (float(args.key[0]), float(args.key[1])) if args.key else (float(phi_values[0]), float(theta_values[0]))
    swept = len(phi_values) * len(theta_values) > 1
    solved = None
    extraction_ok = False
    working_power: list[str] = []
    try:
        with aedt_session(dest / src.name, design, args.aedt_version) as hfss:
            probes.hfss = hfss
            P = probes.run

            print("== 2. what AEDT holds")
            report["meta"]["versions"] = aedt_versions(hfss)
            print("versions:", report["meta"]["versions"])
            P("designs", "designs", lambda: list(hfss.design_list))
            P("setups", "setups", lambda: list(hfss.odesign.GetModule("AnalysisSetup").GetSetups()))
            P("solutions", "solutions", lambda: list(hfss.existing_analysis_sweeps))
            solved = P("solved_variations", f"solved variations of {solution}", lambda: solved_variations(hfss, solution))
            P("variables", "design variables",
              lambda: {str(v): str(hfss.odesign.GetVariableValue(v)) for v in hfss.odesign.GetVariables()})
            P("excitations", "excitations", lambda: list(hfss.odesign.GetModule("BoundarySetup").GetExcitations()))
            P("plane_wave", "plane wave definition",
              lambda: {k: str(v) for b in hfss.boundaries if "Incident" in str(b.type) for k, v in dict(b.props).items()})

            if args.solve:
                print(f"== 2b. solve the copy ({cfg.solver.cores} cores; HFSS solver licence)")
                t0 = time.perf_counter()
                P("solve", "solve", lambda: solve(hfss, cfg.solver.cores))
                print(f"  solve took {time.perf_counter() - t0:.1f} s")
                solved = P("solved_after_solve", f"solved variations of {solution} after the solve",
                           lambda: solved_variations(hfss, solution))
                print("-- the copy's results folder before saving")
                report["meta"]["copy_results_before_save"] = list_results(copy_results, limit=30, tails=False)
                P("save", "save the copy", lambda: hfss.save_project())
                print("-- the copy's results folder after saving")
                report["meta"]["copy_results_after_save"] = list_results(copy_results, limit=30, tails=False)

            print("== 3. the calculator's inputs")
            fields = hfss.odesign.GetModule("FieldsReporter")
            for name in (OUTGOING_POWER_EXPRESSION, "Vector_RealPoynting"):
                P(f"expr_{name}", f"named expression {name} exists", lambda name=name: fields.DoesNamedExpressionExists(name))
            P("face_lists", "user lists in the modeler", lambda: [str(lst.name) for lst in hfss.modeler.user_lists])
            P("face_list_faces", f"faces of list {EXIT_FACE_LIST}",
              lambda: [lst.props for lst in hfss.modeler.user_lists if str(lst.name) == EXIT_FACE_LIST])

            faces = P("faces", "entrance and exit faces", lambda: select_faces(hfss, cfg.geometry))
            ingoing = incoming_power_w(exc.ei_v_per_m, faces[0].area_mm2) if faces else None
            report["meta"]["incoming_power_w"] = ingoing

            phase = ["Phase:=", "0deg"]
            ang = ["IWavePhi:=", f"{phi0}deg", "IWaveTheta:=", f"{theta0}deg"]
            print(f"== 3b. polarization as integer or text, outgoing_power at key ({phi0}, {theta0})"
                  f"{' with the incident angles' if swept else ''}")
            # The order separates "an integer is ignored" (every integer 1 gives the integer-0 value) from "the last
            # evaluated variation is reused" (an integer 1 after a text 1 gives the Ephi=1 value).
            sequence = [("int 0", 0), ("int 1", 1), ("text 1", "1"), ("int 1 after text 1", 1), ("text 0", "0"),
                        ("int 1 after text 0", 1)]
            for i, (text, value) in enumerate(sequence):
                path = dest / f"pol_{i:02d}.fld"
                variation = ["Ephi:=", value, "Freq:=", label] + (ang if swept else []) + phase

                def evaluate_pol(variation=variation, path=path):
                    if path.exists():
                        path.unlink()
                    fields.CalcStack("clear")
                    fields.CopyNamedExprToStack(OUTGOING_POWER_EXPRESSION)
                    fields.CalculatorWrite(str(path), ["Solution:=", solution], variation)
                    wait_for_file(path, PROBE_WAIT_S)
                    value = read_calculator_scalar(path)
                    return {"value": value, "T": value / ingoing} if ingoing else {"value": value}

                print(f"  arguments: {variation}")
                P(f"pol_{i}", f"outgoing_power, Ephi {text}", evaluate_pol)
            print("== 4. calculator evaluations, each written with CalculatorWrite and read back")
            ang_int = ["IWavePhi:=", f"{int(phi0)}deg", "IWaveTheta:=", f"{int(theta0)}deg"]
            matrix = [
                ("const_freq_phase", "constant 1, Freq and Phase only", "scalar", ["Freq:=", label] + phase),
                ("const_runner", "constant 1, the runner's arguments", "scalar", intrinsics(0, label, phi0, theta0) + phase),
                ("power_runner", "outgoing_power, the runner's arguments (Ephi int, x.0deg angles, Phase)", "power",
                 intrinsics(0, label, phi0, theta0) + phase),
                ("power_ephi_text", "outgoing_power, Ephi as text", "power", ["Ephi:=", "0", "Freq:=", label] + ang + phase),
                ("power_whole_deg", "outgoing_power, Ephi as text, whole-degree angles", "power",
                 ["Ephi:=", "0", "Freq:=", label] + ang_int + phase),
                ("power_no_angles_text", "outgoing_power, Ephi as text, no incident angles", "power",
                 ["Ephi:=", "0", "Freq:=", label] + phase),
                ("power_no_angles_int", "outgoing_power, Ephi int, no incident angles", "power",
                 ["Ephi:=", 0, "Freq:=", label] + phase),
                ("power_freq_phase", "outgoing_power, Freq and Phase only", "power", ["Freq:=", label] + phase),
                ("power_ephi1_text", "outgoing_power, Ephi=1 as text, no incident angles", "power",
                 ["Ephi:=", "1", "Freq:=", label] + phase),
            ]
            for i, (key, text, kind, variation) in enumerate(matrix):
                path = dest / f"calc_{i:02d}_{key}.fld"

                def evaluate(kind=kind, variation=variation, path=path):
                    if path.exists():
                        path.unlink()
                    fields.CalcStack("clear")
                    if kind == "scalar":
                        fields.EnterScalar(1)
                    else:
                        fields.CopyNamedExprToStack(OUTGOING_POWER_EXPRESSION)
                    fields.CalculatorWrite(str(path), ["Solution:=", solution], variation)
                    wait_for_file(path, PROBE_WAIT_S)
                    value = read_calculator_scalar(path)
                    out = {"value": value, "file_head": file_head(path)}
                    if kind == "power" and ingoing:
                        out["T"] = value / ingoing
                    return out

                print(f"  arguments: {variation}")
                result = P(key, text, evaluate)
                if result and kind == "power":
                    working_power.append(text)
            P("power_pyaedt_evaluate", "outgoing_power via PyAEDT fields_calculator.evaluate (its own variation)",
              lambda: hfss.post.fields_calculator.evaluate(OUTGOING_POWER_EXPRESSION, setup=solution))
            try:
                fields.CalcStack("clear")
            except Exception:  # noqa: BLE001
                pass

            print("== 5. the field exports alone")
            if faces:
                entrance, exit_face = faces
                grid = far_field_grid(ff.theta_lower_deg, ff.theta_upper_deg, ff.phi_lower_deg, ff.phi_upper_deg,
                                      ff.a_mm, ff.b_mm, freq, ff.fineness, ff.min_coarseness_deg, ff.max_coarseness_deg)
                ctx = ExtractionContext(
                    frequency_ghz=freq, solution=solution, sphere_name=radiation_sphere_name(freq),
                    phi_values=phi_values, theta_values=theta_values, grid=grid, exit_field=cfg.exit_field,
                    exit_face=exit_face, exit_cs_name=EXIT_CS, exit_cs_x=cfg.geometry.exit_cs_x, exit_cs_y=cfg.geometry.exit_cs_y,
                    scratch_dir=scratch, timeout_s=cfg.output.export_timeout_s,
                    incoming_power_w=ingoing,
                )
                face_geometry = P("exit_face", "exit-face description", lambda: describe_exit_face(hfss, ctx))
                lattice = exit_grid(face_geometry, cfg.exit_field) if face_geometry is not None and cfg.exit_field.manual else None
                if lattice is not None:
                    range_min, range_max = lattice.range_strings(hfss.modeler.model_units)
                    resolution = [f"{v}mm" for v in cfg.exit_field.resolution_mm]
                    grid_variants = [
                        ("grid_int1", "ExportOnGrid, Ephi integer 1", ["Ephi:=", 1, "Freq:=", label] + (ang if swept else []) + phase),
                        ("grid_text1", "ExportOnGrid, Ephi text 1", ["Ephi:=", "1", "Freq:=", label] + (ang if swept else []) + phase),
                        ("grid_runner", "ExportOnGrid, the runner's arguments", intrinsics(0, label, phi0, theta0) + phase),
                        ("grid_ephi_text", "ExportOnGrid, Ephi as text", ["Ephi:=", "0", "Freq:=", label] + ang + phase),
                        ("grid_no_angles", "ExportOnGrid, Ephi as text, no incident angles", ["Ephi:=", "0", "Freq:=", label] + phase),
                    ]
                    for key, text, variation in grid_variants:
                        path = dest / f"{key}.fld"

                        def export(variation=variation, path=path):
                            if path.exists():
                                path.unlink()
                            fields.CalcStack("clear")
                            fields.EnterQty("E")
                            fields.CalcOp("Smooth")
                            fields.ExportOnGrid(
                                str(path), range_min, range_max, resolution, solution, variation,
                                ["NAME:ExportOption", "IncludePtInOutput:=", True, "RefCSName:=", EXIT_CS,
                                 "PtInSI:=", True, "FieldInRefCS:=", False],
                                "Cartesian", ["0mm", "0mm", "0mm"], False,
                            )
                            wait_for_file(path, cfg.output.export_timeout_s)
                            df = read_exit_field_fld(path)
                            mag = (df[["Ex_real", "Ex_imag", "Ey_real", "Ey_imag", "Ez_real", "Ez_imag"]] ** 2).sum(axis=1) ** 0.5
                            return f"{len(df)} points (lattice {lattice.lattice_points}), max |E| {float(mag.max()):.4e} V/m"

                        print(f"  arguments: {variation}")
                        P(key, text, export)
                radiation = hfss.odesign.GetModule("RadField")
                far_variants = [
                    ("far_runner", "ExportFieldsToFile, the runner's keys", intrinsic_variation_key(label, phi0, theta0)),
                    ("far_no_angles", "ExportFieldsToFile, Freq only in the intrinsic key", f"Freq='{label}'"),
                ]
                for key, text, intrinsic_key in far_variants:
                    path = dest / f"{key}.ffd"

                    def export_far(intrinsic_key=intrinsic_key, path=path):
                        if path.exists():
                            path.unlink()
                        radiation.ExportFieldsToFile([
                            "ExportFileName:=", str(path), "SetupName:=", ctx.sphere_name,
                            "IntrinsicVariationKey:=", intrinsic_key, "DesignVariationKey:=", "Ephi='0'",
                            "SolutionName:=", solution, "Quantity:=", "",
                        ])
                        wait_for_file(path, cfg.output.export_timeout_s)
                        df = read_far_field_ffd(path, grid.theta_step_deg, grid.phi_step_deg)
                        return f"{len(df)} points (grid {grid.points_per_angle})"

                    print(f"  intrinsic key: {intrinsic_key!r}")
                    P(key, text, export_far)

                print("== 6. the runner's own extraction on the copy")
                n_angles = len(phi_values) * len(theta_values)
                expect_wg = lattice.lattice_points * n_angles if lattice is not None else None
                expect_ff = grid.points_per_angle * n_angles
                print(f"  expected rows: waveguide {expect_wg} at most, far_field {expect_ff}")
                extraction_ok = face_geometry is not None
                for ephi in cfg.output.polarizations:
                    power = P(f"runner_power_{ephi}", f"runner: outgoing power Ephi={ephi}",
                              lambda e=ephi: evaluate_outgoing_power(hfss, ctx, e))
                    if power is not None:
                        for (phi, theta), value in power.items():
                            print(f"    Ephi={ephi} phi={phi} theta={theta}: P_out={value:.6e} W  T={value / ingoing:.6f}")
                    out_dir = dest / dataset_dir_name(cfg.project.dataset_id, freq, ephi)
                    out_dir.mkdir(exist_ok=True)
                    if face_geometry is not None and power is not None:
                        wg = P(f"runner_waveguide_{ephi}", f"runner: waveguide table Ephi={ephi}",
                               lambda e=ephi: extract_waveguide(hfss, ctx, e, face_geometry))
                        if wg is not None:
                            wg.to_csv(out_dir / "waveguide.csv", index=False)
                            print(f"    waveguide.csv: {len(wg)} rows (expected {expect_wg} at most) -> {out_dir}")
                        else:
                            extraction_ok = False
                    else:
                        extraction_ok = False
                    far = P(f"runner_far_field_{ephi}", f"runner: far-field table Ephi={ephi}",
                            lambda e=ephi: extract_far_field(hfss, ctx, e))
                    if far is not None:
                        far.to_csv(out_dir / "far_field.csv", index=False)
                        print(f"    far_field.csv: {len(far)} rows (expected {expect_ff}) -> {out_dir}")
                    else:
                        extraction_ok = False
            else:
                print("  skipped: the faces could not be selected")

            print("== 7. legacy read-back, for the record")

            def legacy():
                fields.CalcStack("clear")
                fields.CopyNamedExprToStack(OUTGOING_POWER_EXPRESSION)
                fields.ClcEval(solution, intrinsics(0, label, phi0, theta0), "Fields")
                return fields.GetTopEntryValue(solution, intrinsics(0, label, phi0, theta0))

            P("legacy_readback", "ClcEval then GetTopEntryValue", legacy)
            try:
                fields.CalcStack("clear")
            except Exception:  # noqa: BLE001
                pass

            print("== 8. convergence, solver profile and mesh statistics")
            for ephi in (0, 1):
                variation = f"Ephi='{ephi}'"
                text = P(f"convergence_{ephi}", f"convergence table {variation}",
                         lambda v=variation, e=ephi: export_convergence_text(hfss, setup, v, dest / f"convergence_Ephi{e}.txt"))
                if text:
                    for line in text.splitlines()[-8:]:
                        print("     " + line)
                P(f"profile_{ephi}", f"solver profile {variation}",
                  lambda v=variation, e=ephi: file_tail(Path(hfss.export_profile(setup, v, str(dest / f"profile_Ephi{e}.prof")))))
                P(f"mesh_stats_{ephi}", f"mesh statistics {variation}",
                  lambda v=variation, e=ephi: file_tail(Path(hfss.export_mesh_stats(setup, v, str(dest / f"mesh_Ephi{e}.ms")))))
    finally:
        report["meta"]["extraction_ok"] = extraction_ok
        report["meta"]["working_power_variants"] = working_power
        (dest / "diagnose.json").write_text(json.dumps(report, indent=2, default=str) + "\n")
        probes.summary()
        print(f"report written to {dest / 'diagnose.json'}")

    print(f"FIELDS: {'present' if solved else 'absent'} (solved variations: {solved})")
    print(f"WORKING POWER VARIANTS: {working_power or 'none'}")
    print(f"EXTRACTION {'PASS' if extraction_ok else 'FAIL'}: tables in {dest}")
    return 0 if extraction_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
