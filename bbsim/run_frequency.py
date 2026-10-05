"""CLI: solve one frequency headlessly from a config file and export the two CSV tables."""
from __future__ import annotations

import argparse
import json
import logging
import os
import shutil
import sys
import time
from collections.abc import Mapping
from dataclasses import asdict
from pathlib import Path

import numpy as np

from bbsim.config import RunConfig, config_to_dict, load_config
from bbsim.dataset_record import build_dataset_record
from bbsim.manifest import build_manifest, sha256_file, write_json
from bbsim.naming import dataset_record_name, frequency_design_name, frequency_label, last_adaptive_solution
from bbsim.paths import JobDirs, create_job_dirs, default_job_id
from bbsim.sampling import incident_angle_values, incoming_power_w
from bbsim.schema import KEY_COLUMNS

log = logging.getLogger("bbsim.run")


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Headless HFSS gap simulation for one frequency.")
    p.add_argument("--config", type=Path, required=True, help="TOML run configuration (see configs/)")
    p.add_argument("--frequency-ghz", type=float, default=None)
    p.add_argument("--cores", type=int, default=None)
    p.add_argument("--output-root", default=None)
    p.add_argument("--project", default=None, help="Prepared .aedt project path")
    p.add_argument("--polarizations", type=int, nargs="+", default=None, choices=[0, 1])
    p.add_argument("--aedt-version", default=None)
    p.add_argument("--job-id", default=None, help="Defaults to slurm<SLURM_JOB_ID> or a UTC timestamp")
    p.add_argument("--log-level", default="INFO")
    return p.parse_args(argv)


def overrides_from_args(args: argparse.Namespace) -> dict[str, object]:
    return {
        "solver.frequency_ghz": args.frequency_ghz,
        "solver.cores": args.cores,
        "output.root": args.output_root,
        "project.path": args.project,
        "output.polarizations": args.polarizations,
        "project.aedt_version": args.aedt_version,
    }


def allocation_warning(cores: int, env: Mapping[str, str] | None = None) -> str | None:
    """Message when the configured solver cores disagree with the Slurm allocation, else None."""
    env = os.environ if env is None else env
    allocated = env.get("SLURM_CPUS_PER_TASK")
    if allocated and allocated.isdigit() and int(allocated) != cores:
        return f"solver.cores = {cores} but SLURM_CPUS_PER_TASK = {allocated}; pass --cores to match the allocation"
    return None


def _axis_summary(values, step: float) -> dict[str, float | int]:
    """The emitted far-field grid along one angle: first and last value, point count, step (degrees)."""
    return {"min": float(values[0]), "max": float(values[-1]), "count": int(len(values)), "step": float(step)}



def _setup_logging(logfile: Path, level: str) -> None:
    logging.basicConfig(level=level.upper(), format="%(asctime)s %(levelname)s %(name)s: %(message)s",
                        handlers=[logging.FileHandler(logfile), logging.StreamHandler()])


def run_job(cfg: RunConfig, dirs: JobDirs, job_id: str) -> Path:
    from bbsim.extract import ExtractionContext, describe_exit_face, exit_grid, extract_far_field, extract_waveguide
    from bbsim.geometry import create_exit_coordinate_system, create_exit_face_list, select_faces
    from bbsim.hfss_setup import (
        EXIT_CS,
        EXIT_FACE_LIST,
        add_outgoing_power_expression,
        assign_radiation_boundaries,
        clear_boundaries_and_excitations,
        create_plane_wave,
        create_setup_and_polarization_sweep,
        export_convergence_text,
        initialize_variables,
        insert_far_field_sphere,
        reset_initial_mesh_settings,
        set_total_fields,
        solve,
    )
    from bbsim.inventory import object_record
    from bbsim.session import aedt_session, aedt_versions

    t_start = time.perf_counter()
    source = cfg.project.path
    if not source.is_file():
        raise FileNotFoundError(f"prepared project not found: {source}")
    job_project = dirs.project / source.name
    shutil.copy2(source, job_project)
    source_sha, job_sha = sha256_file(source), sha256_file(job_project)
    freq = cfg.solver.frequency_ghz
    label = frequency_label(freq)
    run_design = frequency_design_name(cfg.project.design, freq)
    exc = cfg.excitation
    phi_values = incident_angle_values(exc.phi_lower_deg, exc.phi_upper_deg, exc.phi_step_deg)
    theta_values = incident_angle_values(exc.theta_lower_deg, exc.theta_upper_deg, exc.theta_step_deg)
    outputs: dict[str, str] = {}
    convergence: dict[str, str] = {}
    files: dict[str, dict] = {}

    with aedt_session(job_project, cfg.project.design, cfg.project.aedt_version) as hfss:
        versions = aedt_versions(hfss)
        log.info("versions: %s", versions)
        if cfg.project.design not in hfss.design_list:
            raise RuntimeError(f"design {cfg.project.design!r} not in {job_project}; designs: {list(hfss.design_list)}")
        if run_design in hfss.design_list:
            raise RuntimeError(f"design {run_design!r} already exists in the job copy")
        hfss.set_active_design(cfg.project.design)
        hfss.duplicate_design(run_design)
        if hfss.design_name != run_design:
            hfss.set_active_design(run_design)
        if hfss.modeler.model_units != "mm":
            raise RuntimeError(f"model units are {hfss.modeler.model_units!r}; the configuration assumes mm")
        # duplicate_design copies the base design's global mesh settings (manual 0.1 mm model resolution in the
        # GUI-made base designs); the legacy script solved a fresh design with Auto. Restore the reference settings.
        mesh_inherited, mesh_requested = reset_initial_mesh_settings(hfss)
        log.info("initial mesh settings reset to %s (inherited: %s)", mesh_requested, mesh_inherited)

        # Same order as legacy main(): variables, faces/boundaries, exit CS, calculator, setup, sphere, plane wave.
        clear_boundaries_and_excitations(hfss)
        entrance, exit_face = select_faces(hfss, cfg.geometry)
        solid = object_record(hfss, cfg.geometry.object)
        log.info("entrance face %s, exit face %s, material %s", entrance, exit_face, solid["material"])
        initialize_variables(hfss, exc.ei_v_per_m)
        create_exit_face_list(hfss, exit_face, name=EXIT_FACE_LIST)
        assign_radiation_boundaries(hfss, entrance, exit_face)
        create_exit_coordinate_system(hfss, exit_face, cfg.geometry.exit_cs_x, cfg.geometry.exit_cs_y, name=EXIT_CS)
        add_outgoing_power_expression(hfss)
        setup, parametric = create_setup_and_polarization_sweep(hfss, freq, cfg.solver.max_delta_e, cfg.solver.max_passes)
        sphere, grid = insert_far_field_sphere(hfss, cfg.far_field, freq)
        plane_wave = create_plane_wave(hfss, entrance, exc, freq)
        set_total_fields(hfss)
        hfss.save_project()

        log.info("solving %s with %d cores (%d incident angles x 2 polarizations)", setup, cfg.solver.cores,
                 len(phi_values) * len(theta_values))
        t_solve = time.perf_counter()
        if not solve(hfss, cfg.solver.cores):
            raise RuntimeError("hfss.analyze returned False")
        solve_s = time.perf_counter() - t_solve
        log.info("solve finished in %.1f s", solve_s)

        ctx = ExtractionContext(
            frequency_ghz=freq, solution=last_adaptive_solution(freq), sphere_name=sphere,
            phi_values=phi_values, theta_values=theta_values, grid=grid, exit_field=cfg.exit_field,
            exit_face=exit_face, exit_cs_name=EXIT_CS, exit_cs_x=cfg.geometry.exit_cs_x, exit_cs_y=cfg.geometry.exit_cs_y,
            scratch_dir=dirs.scratch, timeout_s=cfg.output.export_timeout_s,
            incoming_power_w=incoming_power_w(exc.ei_v_per_m, entrance.area_mm2),
        )
        face_geometry = describe_exit_face(hfss, ctx)
        lattice = exit_grid(face_geometry, cfg.exit_field) if cfg.exit_field.manual else None
        exit_point_counts: dict[int, dict] = {}
        for ephi in cfg.output.polarizations:
            waveguide = extract_waveguide(hfss, ctx, ephi, face_geometry)
            far_field = extract_far_field(hfss, ctx, ephi)
            out_dir = dirs.dataset_dir(cfg.project.dataset_id, freq, ephi)
            out_dir.mkdir(parents=True, exist_ok=False)
            for name, table in (("waveguide.csv", waveguide), ("far_field.csv", far_field)):
                path = out_dir / name
                table.to_csv(path, index=False)
                files[f"{out_dir.name}/{name}"] = {"sha256": sha256_file(path), "bytes": path.stat().st_size,
                                                   "rows": len(table)}
            outputs[str(ephi)] = str(out_dir)
            exit_point_counts[ephi] = waveguide.groupby(KEY_COLUMNS).size().to_dict()
            convergence[str(ephi)] = export_convergence_text(hfss, setup, f"Ephi='{ephi}'", dirs.logs / f"convergence_Ephi{ephi}.txt")
            log.info("wrote %s (%d waveguide rows, %d far-field rows)", out_dir, len(waveguide), len(far_field))
        # Geant4 pairs the Ephi=0 and Ephi=1 exit points of each incident angle by row index.
        first_ephi, first_counts = next(iter(exit_point_counts.items()))
        for ephi, counts in exit_point_counts.items():
            if counts != first_counts:
                raise RuntimeError(
                    f"exit-point counts per incident angle differ between polarization Ephi={first_ephi} and "
                    f"Ephi={ephi}; Geant4 pairs them by row index. Ephi={first_ephi}: {first_counts}; Ephi={ephi}: {counts}"
                )
        retained = int(next(iter(first_counts.values())))
        lattice_points = lattice.lattice_points if lattice is not None else None
        outside_points = "omitted" if lattice_points is not None and retained < lattice_points else "none"
        hfss.save_project()

    exit_cs_z = [float(v) + 0.0 for v in np.cross(cfg.geometry.exit_cs_x, cfg.geometry.exit_cs_y)]
    manifest = build_manifest(
        config=config_to_dict(cfg),
        job={"id": job_id, "root": str(dirs.root), "wall_time_s": time.perf_counter() - t_start, "solve_time_s": solve_s},
        project={"source_path": str(source), "source_sha256": source_sha, "job_copy": str(job_project),
                 "job_copy_sha256_before_run": job_sha, "base_design": cfg.project.design, "run_design": run_design},
        geometry={"object": cfg.geometry.object, "units": "mm", "material": solid["material"],
                  "bounding_box_mm": solid["bounding_box_mm"],
                  "entrance_face": asdict(entrance), "exit_face": asdict(exit_face),
                  "exit_cs": {"name": EXIT_CS, "x": list(cfg.geometry.exit_cs_x), "y": list(cfg.geometry.exit_cs_y),
                              "z": exit_cs_z},
                  "pose": cfg.geometry.pose, "symmetry": list(cfg.geometry.symmetry),
                  "cross_section": face_geometry.cross_section, "bounds_method": face_geometry.bounds_method,
                  "exit_face_vertex_count": face_geometry.vertex_count,
                  "edge_sampling_sagitta_mm": face_geometry.sagitta_mm},
        excitation={"plane_wave": plane_wave, "ei_v_per_m": exc.ei_v_per_m, "incident_phi_deg": phi_values.tolist(),
                    "incident_theta_deg": theta_values.tolist(), "polarizations": list(cfg.output.polarizations),
                    "polarization_convention": "Ephi=0: E_theta=1; Ephi=1: E_phi=1", "incoming_power_w": ctx.incoming_power_w},
        solver={"setup": setup, "parametric_sweep": parametric, "frequency_label": label, "formulation": "TotalFields",
                "max_delta_e": cfg.solver.max_delta_e, "max_passes": cfg.solver.max_passes, "cores": cfg.solver.cores,
                "initial_mesh_settings": {"requested": mesh_requested, "inherited_from_base_design": mesh_inherited},
                "boundaries": {"rbin": "Radiation on entrance", "rbout": "Radiation on exit", "other faces": "default PEC"},
                "convergence": convergence},
        far_field={"sphere": sphere, "theta_step_deg": grid.theta_step_deg, "phi_step_deg": grid.phi_step_deg,
                   "points_per_angle": grid.points_per_angle, "coordinate_system": EXIT_CS, "radiation_surface": EXIT_FACE_LIST,
                   "theta_deg": _axis_summary(grid.theta_values, grid.theta_step_deg),
                   "phi_deg": _axis_summary(grid.phi_values, grid.phi_step_deg)},
        exit_field={"manual": cfg.exit_field.manual, "resolution_mm": list(cfg.exit_field.resolution_mm),
                    "points_in_si": True, "field_in_ref_cs": False, "coordinate_system": EXIT_CS,
                    "grid": lattice.axis_record() if lattice is not None else None, "lattice_points": lattice_points,
                    "retained_points_per_key": retained, "outside_points": outside_points},
        outputs=outputs,
        versions=versions,
    )
    manifest_path = dirs.root / "manifest.json"
    write_json(manifest_path, manifest)
    for out_dir in outputs.values():
        shutil.copy2(manifest_path, Path(out_dir) / "manifest.json")
    # Built from the manifest exactly as written, so the sidecar cannot disagree with it.
    record = build_dataset_record(json.loads(manifest_path.read_text()), sha256_file(manifest_path), files)
    write_json(dirs.root / dataset_record_name(cfg.project.dataset_id, freq), record)
    return manifest_path


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    cfg = load_config(args.config, overrides_from_args(args))
    if not cfg.project.path.is_file():
        print(f"prepared project not found: {cfg.project.path} (run prepare_hfss_project.py first, or pass --project)",
              file=sys.stderr)
        return 2
    job_id = args.job_id or default_job_id()
    dirs = create_job_dirs(cfg.output.root, cfg.project.dataset_id, cfg.solver.frequency_ghz, job_id)
    _setup_logging(dirs.logs / "run.log", args.log_level)
    write_json(dirs.root / "effective_config.json", config_to_dict(cfg))
    log.info("job %s -> %s", job_id, dirs.root)
    warning = allocation_warning(cfg.solver.cores)
    if warning:
        log.warning(warning)
    try:
        manifest = run_job(cfg, dirs, job_id)
    except Exception:  # noqa: BLE001 - report and exit non-zero; the session is already released
        log.exception("job failed")
        return 1
    log.info("done; manifest at %s", manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
