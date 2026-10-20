"""CLI: solve one frequency headlessly from a config file and export the two CSV tables.

Two phases share one AEDT session: ``prepare_and_solve`` sets up the run design in the job's own project copy and
solves it; ``extract_and_write`` extracts both polarizations into memory, checks them, and only then writes the
CSVs, the manifest and the dataset sidecar (schema-2 design, sections 2.4 and 4). ``output.schema`` selects the
writer: 2 by default (HF-050), 1 for the 1.0 tables until the ``_v2`` trees pass BBRsim's tree check.

``--from-solved JOB_DIR`` replaces the first phase: a finished job's project and ``.aedtresults`` are copied into a
new job beside it, the solved run design is opened, and ``extract_and_write`` writes schema 2 from the copy with no
solve, under the source manifest's configuration (HF-050).
"""
from __future__ import annotations

import argparse
import copy
import json
import logging
import os
import shutil
import sys
import time
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

from bbsim import dataset_v2
from bbsim.config import RunConfig, build_config, config_to_dict, load_config
from bbsim.dataset_record import REEXPORT_TOOL, build_dataset_record, build_dataset_record_v2, modes_record_v2
from bbsim.dataset_v2 import DatasetV2Error
from bbsim.manifest import REPO_ROOT, build_manifest, git_commit, sha256_file, write_json
from bbsim.naming import dataset_record_name, frequency_design_name, frequency_label, last_adaptive_solution
from bbsim.paths import JobDirs, create_job_dirs, default_job_id
from bbsim.sampling import incident_angle_values, incoming_power_w
from bbsim.schema import KEY_COLUMNS

log = logging.getLogger("bbsim.run")


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Headless HFSS gap simulation for one frequency.")
    p.add_argument("--config", type=Path, default=None,
                   help="TOML run configuration (see configs/); required unless --from-solved")
    p.add_argument("--from-solved", type=Path, default=None, metavar="JOB_DIR",
                   help="Schema 2 re-export with no solve: copy this finished job's project and .aedtresults into a "
                        "new job beside it (--job-id required) and extract from the copy; the configuration is the "
                        "job's manifest's, so --config and the overrides below are refused")
    p.add_argument("--frequency-ghz", type=float, default=None)
    p.add_argument("--cores", type=int, default=None)
    p.add_argument("--output-root", default=None)
    p.add_argument("--project", default=None, help="Prepared .aedt project path")
    p.add_argument("--polarizations", type=int, nargs="+", default=None, choices=[0, 1])
    p.add_argument("--aedt-version", default=None)
    p.add_argument("--job-id", default=None, help="Defaults to slurm<SLURM_JOB_ID> or a UTC timestamp")
    p.add_argument("--column-order-reference", type=Path, default=None, metavar="JSON",
                   help="Schema 2: the exporter column-order reference (HF-041) that the crack1 500 GHz reference "
                        "job keeps as scratch/column_order_reference.json; required by every other schema-2 job")
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


@contextmanager
def aedt_messages_on_error(hfss: Any, seen: set[str]) -> Iterator[None]:
    """On an exception in the block, log AEDT's message-manager entries not in ``seen``, then re-raise.

    PyAEDT turns a failed AEDT command into a bare "Failed to execute gRPC AEDT command"; AEDT's reason is
    only in its message manager (HPC step 4, baseline2).
    """
    from bbsim.hfss_setup import aedt_messages

    try:
        yield
    except Exception:
        new = [m for m in aedt_messages(hfss) if m not in seen]
        for message in new or ["(no new AEDT messages)"]:
            log.warning("aedt: %s", message)
        raise


@dataclass(frozen=True)
class StagedProject:
    """The job's own project copy, made before any AEDT session, and what was checked and hashed on the way."""

    source: Path
    source_sha: str
    job_project: Path
    job_sha: str
    recorded_id: str | None
    record: Path | None


@dataclass(frozen=True)
class SolvedContext:
    """What ``prepare_and_solve`` (or a reuse of a solved job) hands to ``extract_and_write``: the extraction context,
    the exit face and its lattice, the names of everything created in the run design, the solve's evidence, and the
    project provenance. ``solve_s`` is None when no solve ran."""

    ctx: Any
    face_geometry: Any
    lattice: Any
    entrance: Any
    exit_face: Any
    solid: dict
    setup: str
    parametric: str
    sphere: str
    grid: Any
    plane_wave: str
    wall_facets: dict | None
    mesh_inherited: Any
    mesh_requested: dict
    variations: list
    messages: list
    solve_s: float | None
    source: Path
    source_sha: str
    job_project: Path
    job_sha: str
    recorded_id: str | None
    record: Path | None
    run_design: str
    versions: dict


def is_reference_job(cfg: RunConfig) -> bool:
    """True for the one dataset and frequency that measures the exporter column-order reference itself (HF-041)."""
    return (cfg.project.dataset_id == dataset_v2.REFERENCE_DATASET
            and float(cfg.solver.frequency_ghz) == dataset_v2.REFERENCE_GHZ)


def load_reference(path: Path) -> dict:
    """The column-order reference object (field list, section 2) from a JSON file, such as the reference job's
    ``scratch/column_order_reference.json``. Its content is validated against the run's versions later."""
    reference = json.loads(Path(path).read_text())
    if not isinstance(reference, dict):
        raise ValueError(f"{path} must hold the column-order reference object, got {type(reference).__name__}")
    return reference


def check_before_aedt(cfg: RunConfig, reference: Mapping | None) -> None:
    """Refuse, from the configuration alone and before any job copy or AEDT session, a run whose output cannot be
    written: under schema 1 a column-order reference (schema 1 neither needs nor applies it); under schema 2
    anything but both polarizations, a declared theta grid without a 90 deg key and one above it (``cell_90``,
    Ruling 4), and a missing reference on any job but the reference job (HF-041)."""
    if cfg.output.schema == 1:
        if reference is not None:
            raise ValueError("a column-order reference applies to schema 2 only; this run writes schema 1 "
                             "(output.schema = 1), which neither requires nor applies it")
        return
    if sorted(cfg.output.polarizations) != [0, 1]:
        raise DatasetV2Error(f"schema 2 needs both polarizations (output.polarizations = [0, 1]), got "
                             f"{list(cfg.output.polarizations)}; a one-polarization diagnostic writes schema 1")
    exc = cfg.excitation
    theta_values = incident_angle_values(exc.theta_lower_deg, exc.theta_upper_deg, exc.theta_step_deg)
    try:
        dataset_v2.cell_90_band_deg(theta_values)
    except DatasetV2Error as err:
        raise DatasetV2Error(f"schema 2: {err}; a diagnostic without that key writes schema 1") from err
    if reference is None and not is_reference_job(cfg):
        raise DatasetV2Error(
            f"schema 2 needs the column-order reference (--column-order-reference JSON): only "
            f"{dataset_v2.REFERENCE_DATASET} at {dataset_v2.REFERENCE_GHZ:g} GHz measures its own, and "
            f"{cfg.project.dataset_id} at {float(cfg.solver.frequency_ghz):g} GHz must embed that one (HF-041)")
    if reference is None:
        # The reference job measures its two cases at fixed keys; a grid without them would fail after the solve.
        phi_values = incident_angle_values(exc.phi_lower_deg, exc.phi_upper_deg, exc.phi_step_deg)
        missing = [(phi, theta) for _, phi, theta in dataset_v2.REFERENCE_CASES
                   if not (np.any(np.abs(phi_values - phi) <= 1e-9) and np.any(np.abs(theta_values - theta) <= 1e-9))]
        if missing:
            raise DatasetV2Error(f"the column-order reference job needs the keys (phi, theta) = "
                                 f"{[c[1:] for c in dataset_v2.REFERENCE_CASES]} deg in its grid; missing {missing}")
    if reference is not None and not isinstance(reference, Mapping):
        raise DatasetV2Error(f"the column-order reference must be an object, got {type(reference).__name__}")


def stage_project(cfg: RunConfig, dirs: JobDirs) -> StagedProject:
    """Check the config's dataset_id against the record beside the prepared project, then copy the project into the
    job directory; both before any AEDT session, which opens the copy."""
    from bbsim.mapping import recorded_dataset_id

    source = cfg.project.path
    if not source.is_file():
        raise FileNotFoundError(f"prepared project not found: {source}")
    # The run config names the dataset; the record beside the project says which dataset its design is.
    recorded_id, record = recorded_dataset_id(source, cfg.project.design, cfg.geometry.object)
    if recorded_id is not None and recorded_id != cfg.project.dataset_id:
        raise ValueError(f"config dataset_id {cfg.project.dataset_id!r} disagrees with {record}, which gives "
                         f"{recorded_id!r} for design {cfg.project.design!r} object {cfg.geometry.object!r}")
    if recorded_id is None:
        log.warning("dataset_id %r not confirmed: %s", cfg.project.dataset_id,
                    f"{record} has no entry for {cfg.project.design}/{cfg.geometry.object}" if record
                    else f"no .mapping.json or .build.json beside {source}")
    job_project = dirs.project / source.name
    shutil.copy2(source, job_project)
    return StagedProject(source=source, source_sha=sha256_file(source), job_project=job_project,
                         job_sha=sha256_file(job_project), recorded_id=recorded_id, record=record)


def run_job(cfg: RunConfig, dirs: JobDirs, job_id: str, reference: dict | None = None) -> Path:
    """One frequency: stage the project, open one AEDT session, solve, extract both polarizations, write.

    ``reference`` is the column-order reference object every schema-2 job but the reference job embeds."""
    from bbsim.session import aedt_session, aedt_versions

    t_start = time.perf_counter()
    check_before_aedt(cfg, reference)
    staged = stage_project(cfg, dirs)
    with aedt_session(staged.job_project, cfg.project.design, cfg.project.aedt_version) as hfss:
        if cfg.output.schema == 2:
            # A stale or unreadable version would stop the run after its licensed solve; it costs seconds here.
            _check_versions_for_reference(reference, aedt_versions(hfss))
        solved = prepare_and_solve(cfg, dirs, job_id, hfss, staged)
        return extract_and_write(cfg, dirs, job_id, hfss, solved, t_start, reference, None)


def _require_v2_exit_normal(cfg: RunConfig, exit_face: Any) -> None:
    """Schema 2: exit_flux integrates ScalarZ of Re S, which is the flux Re S.n (HF-039(1), CDX-009) only on an exit
    face whose outward normal is +z; any other pose is refused, a pose the config leaves unchecked included."""
    if tuple(exit_face.outward_normal) != (0.0, 0.0, 1.0):
        raise DatasetV2Error(
            f"schema 2 stores the exit flux Integrate(ScalarZ(Re S)), which is Re S.n only when the exit face's "
            f"outward normal is (0, 0, 1) (pose canonical-z); exit face {exit_face.id} has outward normal "
            f"{list(exit_face.outward_normal)} (config pose {cfg.geometry.pose!r})")


def _require_v2_lattice(cross_section: Mapping, lattice: Any) -> None:
    """Schema 2: the quadrature weights and the modes need a rectangle or a disc on the manual exit lattice."""
    shape = cross_section["shape"]
    if shape not in dataset_v2.QUADRATURE_BY_SHAPE or lattice is None:
        raise DatasetV2Error(f"schema 2 needs a rectangle or disc exit cross-section on the manual exit lattice "
                             f"(quadrature weights, modes); got {shape!r}"
                             + ("" if lattice is not None else " without a lattice"))


def _require_v2_exit_face(cfg: RunConfig, solved: SolvedContext) -> None:
    """Both schema-2 exit-face checks on a solved context, whichever path built it (fresh solve or --from-solved)."""
    _require_v2_exit_normal(cfg, solved.exit_face)
    _require_v2_lattice(solved.face_geometry.cross_section, solved.lattice)


def _check_versions_for_reference(reference: Mapping | None, versions: Mapping[str, Any]) -> None:
    """Schema 2: the reference job's versions must be evidence; a supplied reference must match them (HF-041)."""
    aedt, pyaedt = versions.get("aedt"), versions.get("pyaedt")
    if reference is None:
        dataset_v2.check_versions(aedt, pyaedt)
    else:
        dataset_v2.check_reference_compatible(dict(reference), aedt, pyaedt, REPO_ROOT)


def prepare_and_solve(cfg: RunConfig, dirs: JobDirs, job_id: str, hfss: Any, staged: StagedProject) -> SolvedContext:
    """Set up the run design in the open job copy as the legacy script did, check what can be checked before the
    solve, solve, and verify that every requested polarization was solved. Returns what extraction needs.

    ``staged`` is ``stage_project``'s result: the copy that ``hfss`` opened, made and hashed before the session."""
    from bbsim.extract import ExtractionContext, describe_exit_face, exit_grid
    from bbsim.geometry import create_exit_coordinate_system, create_exit_face_list, select_faces
    from bbsim.hfss_setup import (
        EXIT_CS,
        EXIT_FACE_LIST,
        add_exit_flux_expression,
        add_outgoing_power_expression,
        aedt_messages,
        assign_radiation_boundaries,
        assign_wall_facets,
        clear_boundaries_and_excitations,
        create_plane_wave,
        create_setup_and_polarization_sweep,
        initialize_variables,
        insert_far_field_sphere,
        missing_polarizations,
        reset_initial_mesh_settings,
        set_total_fields,
        solve,
        solved_variations,
    )
    from bbsim.inventory import object_record
    from bbsim.session import aedt_versions

    job_project = staged.job_project
    freq = cfg.solver.frequency_ghz
    run_design = frequency_design_name(cfg.project.design, freq)
    exc = cfg.excitation
    phi_values = incident_angle_values(exc.phi_lower_deg, exc.phi_upper_deg, exc.phi_step_deg)
    theta_values = incident_angle_values(exc.theta_lower_deg, exc.theta_upper_deg, exc.theta_step_deg)

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
    wall_facets = None
    if cfg.solver.wall_normal_deviation_deg is not None:
        wall_facets = assign_wall_facets(hfss, cfg.geometry.object, cfg.solver.wall_normal_deviation_deg)
        log.info("wall facets: %s", wall_facets)

    # Same order as legacy main(): variables, faces/boundaries, exit CS, calculator, setup, sphere, plane wave.
    clear_boundaries_and_excitations(hfss)
    entrance, exit_face = select_faces(hfss, cfg.geometry)
    if cfg.output.schema == 2:
        _require_v2_exit_normal(cfg, exit_face)   # before the solve; extract_and_write checks again
    solid = object_record(hfss, cfg.geometry.object)
    log.info("entrance face %s, exit face %s, material %s", entrance, exit_face, solid["material"])
    initialize_variables(hfss, exc.ei_v_per_m)
    create_exit_face_list(hfss, exit_face, name=EXIT_FACE_LIST)
    assign_radiation_boundaries(hfss, entrance, exit_face)
    create_exit_coordinate_system(hfss, exit_face, cfg.geometry.exit_cs_x, cfg.geometry.exit_cs_y, name=EXIT_CS)
    # add_outgoing_power_expression clears every named expression, so the flux is added after it (Ruling 14).
    add_outgoing_power_expression(hfss)
    add_exit_flux_expression(hfss)
    setup, parametric = create_setup_and_polarization_sweep(hfss, freq, cfg.solver.max_delta_e, cfg.solver.max_passes)
    sphere, grid = insert_far_field_sphere(hfss, cfg.far_field, freq)
    plane_wave = create_plane_wave(hfss, entrance, exc, freq)
    set_total_fields(hfss)

    # The export lattice needs only geometry: describe the exit face before the licensed solve, so a
    # geometry problem costs seconds rather than a solve.
    ctx = ExtractionContext(
        frequency_ghz=freq, solution=last_adaptive_solution(freq), sphere_name=sphere,
        phi_values=phi_values, theta_values=theta_values, grid=grid, exit_field=cfg.exit_field,
        exit_face=exit_face, exit_cs_name=EXIT_CS, exit_cs_x=cfg.geometry.exit_cs_x, exit_cs_y=cfg.geometry.exit_cs_y,
        scratch_dir=dirs.scratch, timeout_s=cfg.output.export_timeout_s,
        incoming_power_w=incoming_power_w(exc.ei_v_per_m, entrance.area_mm2),
    )
    face_geometry = describe_exit_face(hfss, ctx)
    lattice = exit_grid(face_geometry, cfg.exit_field) if cfg.exit_field.manual else None
    if lattice is not None:
        log.info("exit lattice: %d points per incident angle (%s)", lattice.lattice_points, face_geometry.bounds_method)
    if cfg.output.schema == 2:
        _require_v2_lattice(face_geometry.cross_section, lattice)
    hfss.save_project()

    log.info("solving %s with %d cores (%d incident angles x 2 polarizations)", setup, cfg.solver.cores,
             len(phi_values) * len(theta_values))
    t_solve = time.perf_counter()
    if not solve(hfss, cfg.solver.cores):
        raise RuntimeError("hfss.analyze returned False")
    solve_s = time.perf_counter() - t_solve
    log.info("solve finished in %.1f s", solve_s)
    # analyze() returning True proves nothing (HPC step 4, baseline1: 66 s, no adaptive pass, no fields).
    # Ask AEDT what it solved, and keep its own messages: the solver's reason lives only there.
    messages = aedt_messages(hfss)
    for message in messages:
        # AEDT tags its entries [info], [warning], [error]; its errors are warnings here because the
        # runner decides below what is fatal.
        (log.warning if "[error]" in message or "[warning]" in message else log.info)("aedt: %s", message)
    variations = solved_variations(hfss, ctx.solution)
    log.info("solved variations of %s: %s", ctx.solution, variations)
    missing = missing_polarizations(variations, cfg.output.polarizations)
    if missing:
        raise RuntimeError(
            f"no solved variations for Ephi={missing} in {ctx.solution!r}: AEDT holds {variations}; "
            f"the solve ran for {solve_s:.1f} s without producing them. AEDT messages: " + " | ".join(messages)
        )
    # AEDT writes solution data into the project's results folder only when the project is saved; a run
    # that fails during extraction closes the project unsaved (HPC steps 4 and 4a lost their solves that way).
    hfss.save_project()
    return SolvedContext(
        ctx=ctx, face_geometry=face_geometry, lattice=lattice, entrance=entrance, exit_face=exit_face, solid=solid,
        setup=setup, parametric=parametric, sphere=sphere, grid=grid, plane_wave=plane_wave, wall_facets=wall_facets,
        mesh_inherited=mesh_inherited, mesh_requested=mesh_requested, variations=variations, messages=messages,
        solve_s=solve_s, source=staged.source, source_sha=staged.source_sha, job_project=job_project,
        job_sha=staged.job_sha, recorded_id=staged.recorded_id, record=staged.record, run_design=run_design,
        versions=versions,
    )


# BBRsim pairs the two polarizations' exit points by row and requires equal X, Y and Z to this (BBR012 at load).
POSITION_MATCH_M = 1e-12


def _require_shared_exit_points(waveguide0: Any, waveguide1: Any) -> None:
    """Ruling 12: row by row, the Ephi=0 and Ephi=1 tables must hold the same key and the same exit point to
    ``POSITION_MATCH_M`` in X, Y and Z; ``dataset_v2.build`` compares lattice nodes only. Names the first key."""
    keys0 = waveguide0[KEY_COLUMNS].to_numpy(dtype=float)
    keys1 = waveguide1[KEY_COLUMNS].to_numpy(dtype=float)
    if keys0.shape != keys1.shape or not np.array_equal(keys0, keys1):
        raise DatasetV2Error("the Ephi=0 and Ephi=1 waveguide tables do not list the same keys row by row: BBRsim "
                             "pairs the polarizations' exit points by row (BBR012)")
    xyz0 = waveguide0[["X", "Y", "Z"]].to_numpy(dtype=float)
    xyz1 = waveguide1[["X", "Y", "Z"]].to_numpy(dtype=float)
    diff = np.abs(xyz0 - xyz1)
    bad = ~(diff <= POSITION_MATCH_M)
    if bad.any():
        i = int(np.flatnonzero(bad.any(axis=1))[0])
        j = int(np.flatnonzero(bad[i])[0])
        phi, theta = keys0[i]
        row = i - int(np.flatnonzero(np.all(keys0 == keys0[i], axis=1))[0])
        raise DatasetV2Error(
            f"at (phi, theta) = ({phi:g}, {theta:g}) deg the Ephi=0 and Ephi=1 exit points differ in {'XYZ'[j]} by "
            f"{diff[i, j]!r} m at row {row} of the key (limit {POSITION_MATCH_M:g} m): BBRsim pairs them by row (BBR012)")


def _column_order_check(cfg: RunConfig, dirs: JobDirs, solved: SolvedContext,
                        extracted: Mapping[int, dataset_v2.Extracted], reference: Mapping | None
                        ) -> tuple[dict, dict, dict | None]:
    """(block, reference, raw tables to keep as evidence or None). The reference job measures the reference on its
    own raw tables (HF-041); every other job validates the supplied one against its versions and exporter first."""
    aedt, pyaedt = solved.versions.get("aedt"), solved.versions.get("pyaedt")
    if reference is None:
        if not is_reference_job(cfg):
            raise DatasetV2Error(f"schema 2 needs the column-order reference for {cfg.project.dataset_id} at "
                                 f"{float(cfg.solver.frequency_ghz):g} GHz (--column-order-reference JSON; HF-041)")
        raw = {e: extracted[e].waveguide for e in (0, 1)}
        measured = dataset_v2.build_reference(raw, cfg.project.dataset_id, cfg.solver.frequency_ghz, str(dirs.root),
                                              git_commit(), aedt, pyaedt, REPO_ROOT)
        return dataset_v2.order_check_block("reference", measured), measured, raw
    reference = dict(reference)
    dataset_v2.check_reference_compatible(reference, aedt, pyaedt, REPO_ROOT)
    return dataset_v2.order_check_block("shared_exporter", reference), reference, None


def _write_reference_evidence(scratch: Path, reference: Mapping, raw: Mapping[int, Any]) -> None:
    """Ruling 6: the reference job keeps the two raw tables it hashed, byte for byte, and the reference object."""
    for e in (0, 1):
        path = scratch / f"column_order_reference_Ephi{e}.csv"
        path.write_bytes(raw[e].to_csv(index=False).encode())
        if sha256_file(path) != reference["raw_table_sha256"][str(e)]:
            raise RuntimeError(f"{path} does not hash to the reference's raw_table_sha256[{e!r}]")
    write_json(scratch / "column_order_reference.json", reference)


def extract_and_write(cfg: RunConfig, dirs: JobDirs, job_id: str, hfss: Any, solved: SolvedContext, t_start: float,
                      reference: dict | None = None, reexport: dict | None = None) -> Path:
    """Extract every polarization into memory, check, then write the CSVs, the manifest and the sidecar.

    Order (schema-2 design, section 2.4): both polarizations are extracted and checked for finite values; under
    schema 2 the column-order reference is measured (the reference job) or validated (every other job), the 2.0
    tables and sidecar blocks are built (``dataset_v2.build``); then the pairing and lattice checks. Only then is
    anything written, so a failure leaves no CSV, manifest or sidecar. Schema 1 writes the 1.0 tables and sidecar
    exactly as before. ``reexport`` (schema 2 only) is the --from-solved provenance for the sidecar; ``t_start`` is
    the job's start on ``time.perf_counter``. Returns the manifest path.
    """
    from bbsim.convergence import qualify
    from bbsim.extract import RIM_BAND_REL, evaluate_exit_powers, extract_far_field, extract_waveguide
    from bbsim.hfss_setup import EXIT_CS, EXIT_FACE_LIST, export_convergence_text
    from bbsim.readers import require_finite_table

    schema = cfg.output.schema
    if schema == 1 and reexport is not None:
        raise ValueError("a re-export record applies to schema 2 only; this run writes schema 1")
    check_before_aedt(cfg, reference)   # pure; repeated here for callers that did not come through run_job
    if schema == 2:
        _require_v2_exit_face(cfg, solved)   # every path in, a --from-solved context included
    ctx = solved.ctx
    freq = cfg.solver.frequency_ghz
    dataset_id = cfg.project.dataset_id
    extracted: dict[int, dataset_v2.Extracted] = {}
    convergence: dict[str, str] = {}
    rim_log: dict[str, list] = {}

    with aedt_messages_on_error(hfss, set(solved.messages)):
        for ephi in cfg.output.polarizations:
            if schema == 2:
                abs_w, flux_w = evaluate_exit_powers(hfss, ctx, ephi)
                waveguide = extract_waveguide(hfss, ctx, ephi, solved.face_geometry, rim_log,
                                              exit_powers=(abs_w, flux_w))
            else:   # schema 1 evaluates |Re S| inside extract_waveguide and needs no flux
                abs_w, flux_w = {}, {}
                waveguide = extract_waveguide(hfss, ctx, ephi, solved.face_geometry, rim_log)
            far_field = extract_far_field(hfss, ctx, ephi)
            for name, table in (("waveguide.csv", waveguide), ("far_field.csv", far_field)):
                require_finite_table(table, f"{dataset_id} Ephi={ephi} {name}")
            extracted[ephi] = dataset_v2.Extracted(waveguide=waveguide, far_field=far_field, flux_w=flux_w, abs_w=abs_w)
            convergence[str(ephi)] = export_convergence_text(hfss, solved.setup, f"Ephi='{ephi}'",
                                                             dirs.logs / f"convergence_Ephi{ephi}.txt")
            log.info("extracted Ephi=%s (%d waveguide rows, %d far-field rows)", ephi, len(waveguide), len(far_field))
        hfss.save_project()

    # CDX-025: whether each polarization met the configured target. Recorded, not yet a gate (OR-027): the
    # job array's assembly will refuse unqualified datasets.
    qualification = {e: qualify(text, cfg.solver.max_delta_e) for e, text in convergence.items()}
    qualified = bool(qualification) and all(q["qualified"] for q in qualification.values())
    if not qualified:
        log.warning("convergence target %g not met or not shown (CDX-025; recorded in the manifest, not yet a "
                    "gate): %s", cfg.solver.max_delta_e,
                    "; ".join(f"Ephi={e}: {', '.join(q['reasons'])}" for e, q in qualification.items()
                              if not q["qualified"]))

    result = None
    evidence = None
    if schema == 2:
        order_check, reference, evidence = _column_order_check(cfg, dirs, solved, extracted, reference)
        modes = modes_record_v2(solved.face_geometry.cross_section)
        if modes is None:
            raise DatasetV2Error(f"schema 2 has no modes for a {solved.face_geometry.cross_section['shape']!r} "
                                 "cross-section")
        result = dataset_v2.build(
            extracted, grid=solved.lattice, cross_section=solved.face_geometry.cross_section,
            incoming_w=ctx.incoming_power_w, phi_values=ctx.phi_values, theta_values=ctx.theta_values,
            ff_steps_deg=(solved.grid.theta_step_deg, solved.grid.phi_step_deg),
            polarization_filter_limit_ghz=modes["polarization_filter_limit_ghz"], order_check=order_check)
        _require_shared_exit_points(extracted[0].waveguide, extracted[1].waveguide)

    # Geant4 pairs the Ephi=0 and Ephi=1 exit points of each incident angle by row index.
    exit_point_counts = {e: ex.waveguide.groupby(KEY_COLUMNS).size().to_dict() for e, ex in extracted.items()}
    first_ephi, first_counts = next(iter(exit_point_counts.items()))
    for ephi, counts in exit_point_counts.items():
        if counts != first_counts:
            raise RuntimeError(
                f"exit-point counts per incident angle differ between polarization Ephi={first_ephi} and "
                f"Ephi={ephi}; Geant4 pairs them by row index. Ephi={first_ephi}: {first_counts}; Ephi={ephi}: {counts}"
            )
    retained = int(next(iter(first_counts.values())))
    lattice = solved.lattice
    lattice_points = lattice.lattice_points if lattice is not None else None
    if lattice_points is not None and retained > lattice_points:
        # HFSS can omit lattice points outside the solid but never add any: the lattice model is wrong,
        # and so would be the grid recorded in the manifest and the sidecar.
        raise RuntimeError(
            f"exit field has {retained} points per incident angle, more than the {lattice_points} points of "
            f"the ExportOnGrid lattice {lattice.counts}; the lattice model disagrees with the HFSS export"
        )
    outside_points = "omitted" if lattice_points is not None and retained < lattice_points else "none"

    # --- write phase: nothing above has written a CSV, manifest or sidecar --------------------------------------
    if evidence is not None:
        _write_reference_evidence(dirs.scratch, reference, evidence)
    outputs: dict[str, str] = {}
    files: dict[str, dict] = {}
    for ephi in cfg.output.polarizations:
        tables = ((result.waveguide[ephi], result.far_field[ephi]) if result is not None
                  else (extracted[ephi].waveguide, extracted[ephi].far_field))
        out_dir = dirs.dataset_dir(dataset_id, freq, ephi)
        out_dir.mkdir(parents=True, exist_ok=False)
        for name, table in zip(("waveguide.csv", "far_field.csv"), tables):
            path = out_dir / name
            table.to_csv(path, index=False)
            files[f"{out_dir.name}/{name}"] = {"sha256": sha256_file(path), "bytes": path.stat().st_size,
                                               "rows": len(table)}
        outputs[str(ephi)] = str(out_dir)
        log.info("wrote %s (%d waveguide rows, %d far-field rows)", out_dir, len(tables[0]), len(tables[1]))

    exc = cfg.excitation
    entrance, exit_face, face_geometry, grid = solved.entrance, solved.exit_face, solved.face_geometry, solved.grid
    exit_cs_z = [float(v) + 0.0 for v in np.cross(cfg.geometry.exit_cs_x, cfg.geometry.exit_cs_y)]
    solver = {"setup": solved.setup, "parametric_sweep": solved.parametric, "frequency_label": frequency_label(freq),
              "formulation": "TotalFields", "max_delta_e": cfg.solver.max_delta_e,
              "max_passes": cfg.solver.max_passes, "cores": cfg.solver.cores,
              "initial_mesh_settings": {"requested": solved.mesh_requested,
                                        "inherited_from_base_design": solved.mesh_inherited},
              "wall_facets": solved.wall_facets,
              "boundaries": {"rbin": "Radiation on entrance", "rbout": "Radiation on exit", "other faces": "default PEC"},
              "solved_variations": solved.variations, "convergence": convergence, "qualified": qualified,
              "convergence_qualification": {e: {k: q[k] for k in ("achieved", "target", "passes_completed",
                                                                  "max_passes", "converged", "qualified", "reasons")}
                                            for e, q in qualification.items()}}
    if reexport is not None:
        solver["reused_solution"] = reexport.get("source_job")
    exit_field = {"manual": cfg.exit_field.manual, "resolution_mm": list(cfg.exit_field.resolution_mm),
                  "step_mm": list(lattice.step_mm) if lattice is not None else None,
                  "points_in_si": True, "field_in_ref_cs": False, "coordinate_system": EXIT_CS,
                  "grid": lattice.axis_record() if lattice is not None else None, "lattice_points": lattice_points,
                  "retained_points_per_key": retained, "outside_points": outside_points,
                  "rim": {"rule": "every lattice point inside the disc kept; inside points without field set to zero; "
                                  "points outside dropped (ledger HF-029)",
                          "band_rel": RIM_BAND_REL, "per_key": rim_log} if rim_log else None}
    extra: dict[str, Any] = {}
    if result is not None:
        # Fresh pre-clamp powers (field list section 3), and the diagnostics that are not sidecar fields.
        exit_field["raw_powers_per_key"] = result.raw_powers_per_key
        extra["transmittance_checks"] = {"flux_over_abs_re_s": result.flux_over_abs_re_s,
                                         "rho_per_key": result.rho_per_key}
    manifest = build_manifest(
        config=config_to_dict(cfg),
        job={"id": job_id, "root": str(dirs.root), "wall_time_s": time.perf_counter() - t_start,
             "solve_time_s": solved.solve_s},
        project={"source_path": str(solved.source), "source_sha256": solved.source_sha,
                 "job_copy": str(solved.job_project),
                 "dataset_id_record": {"path": str(solved.record) if solved.record else None,
                                       "dataset_id": solved.recorded_id},
                 "job_copy_sha256_before_run": solved.job_sha, "base_design": cfg.project.design,
                 "run_design": solved.run_design},
        geometry={"object": cfg.geometry.object, "units": "mm", "material": solved.solid["material"],
                  "bounding_box_mm": solved.solid["bounding_box_mm"],
                  "entrance_face": asdict(entrance), "exit_face": asdict(exit_face),
                  "exit_cs": {"name": EXIT_CS, "x": list(cfg.geometry.exit_cs_x), "y": list(cfg.geometry.exit_cs_y),
                              "z": exit_cs_z},
                  "pose": cfg.geometry.pose, "symmetry": list(cfg.geometry.symmetry),
                  "cross_section": face_geometry.cross_section, "bounds_method": face_geometry.bounds_method,
                  "exit_face_vertex_count": face_geometry.vertex_count,
                  "edge_sampling_sagitta_mm": face_geometry.sagitta_mm},
        excitation={"plane_wave": solved.plane_wave, "ei_v_per_m": exc.ei_v_per_m,
                    "incident_phi_deg": ctx.phi_values.tolist(), "incident_theta_deg": ctx.theta_values.tolist(),
                    "polarizations": list(cfg.output.polarizations),
                    "polarization_convention": "Ephi=0: E_theta=1; Ephi=1: E_phi=1",
                    "incoming_power_w": ctx.incoming_power_w},
        solver=solver,
        far_field={"sphere": solved.sphere, "theta_step_deg": grid.theta_step_deg, "phi_step_deg": grid.phi_step_deg,
                   "points_per_angle": grid.points_per_angle, "coordinate_system": EXIT_CS,
                   "radiation_surface": EXIT_FACE_LIST,
                   "theta_deg": _axis_summary(grid.theta_values, grid.theta_step_deg),
                   "phi_deg": _axis_summary(grid.phi_values, grid.phi_step_deg)},
        exit_field=exit_field,
        outputs=outputs,
        versions=solved.versions,
        **extra,
    )
    manifest_path = dirs.root / "manifest.json"
    write_json(manifest_path, manifest)
    for out_dir in outputs.values():
        shutil.copy2(manifest_path, Path(out_dir) / "manifest.json")
    # Built from the manifest exactly as written, so the sidecar cannot disagree with it. Never edited afterwards:
    # the record shares objects with the parsed manifest.
    written, manifest_sha = json.loads(manifest_path.read_text()), sha256_file(manifest_path)
    if result is None:
        record = build_dataset_record(written, manifest_sha, files)
    else:
        record = build_dataset_record_v2(written, manifest_sha, files, result.blocks, reexport)
    write_json(dirs.root / dataset_record_name(dataset_id, freq), record)
    log.info("wrote the schema-%d sidecar %s", schema, dataset_record_name(dataset_id, freq))
    return manifest_path


# --- --from-solved: a fresh schema-2 extraction on a copy of a solved job, no solve (schema-2 design 2.4, HF-050) ----

@dataclass(frozen=True)
class SourceJob:
    """A finished job to re-export: its directory, its manifest exactly as read (never modified), the run
    configuration built from that manifest's ``config`` block with ``output.schema`` set to 2, its solved project
    and ``.aedtresults`` folder, and the manifest's sha256."""

    job_dir: Path
    manifest: dict
    config: RunConfig
    project: Path
    results: Path
    manifest_sha256: str


def _manifest_field(manifest: Mapping, section: str, key: str) -> Any:
    block = manifest.get(section)
    if not isinstance(block, Mapping) or key not in block:
        raise ValueError(f"the source manifest has no {section}.{key}")
    return block[key]


def _check_source_records(cfg: RunConfig, manifest: Mapping) -> None:
    """The re-export asks AEDT for the source's run design, solution, far-field sphere and incidence keys by the names
    and grids this code derives from the source's config; refuse a source that recorded different ones, before any
    copy, rather than halfway through the extraction."""
    from bbsim.hfss_setup import far_field_sphere
    from bbsim.naming import setup_name

    freq = cfg.solver.frequency_ghz
    exc = cfg.excitation
    sphere, grid = far_field_sphere(cfg.far_field, freq)
    derived = {
        ("project", "run_design"): frequency_design_name(cfg.project.design, freq),
        ("solver", "setup"): setup_name(freq),
        ("far_field", "sphere"): sphere,
        ("far_field", "theta_step_deg"): grid.theta_step_deg,
        ("far_field", "phi_step_deg"): grid.phi_step_deg,
        ("far_field", "points_per_angle"): grid.points_per_angle,
        ("excitation", "incident_phi_deg"): incident_angle_values(exc.phi_lower_deg, exc.phi_upper_deg,
                                                                  exc.phi_step_deg).tolist(),
        ("excitation", "incident_theta_deg"): incident_angle_values(exc.theta_lower_deg, exc.theta_upper_deg,
                                                                    exc.theta_step_deg).tolist(),
    }

    def same(recorded: Any, value: Any) -> bool:
        if isinstance(value, str) or isinstance(recorded, str):
            return recorded == value
        try:
            return bool(np.shape(recorded) == np.shape(value) and np.allclose(recorded, value, rtol=1e-12, atol=1e-9))
        except (TypeError, ValueError):
            return False

    differ = [f"{section}.{key}: recorded {_manifest_field(manifest, section, key)!r}, derived {value!r}"
              for (section, key), value in derived.items() if not same(_manifest_field(manifest, section, key), value)]
    if differ:
        raise ValueError("the source job's manifest records what this code would not reproduce from its config: "
                         + "; ".join(differ))


def load_source_job(job_dir: Path) -> SourceJob:
    """Read and check a finished job for --from-solved, without changing it: its ``manifest.json``, its solved
    project (``project/<name>.aedt``, the manifest's ``job_copy``) and that project's ``.aedtresults`` folder must
    exist (design section 4). The configuration is a deep copy of the manifest's ``config`` block with
    ``output.schema`` = 2, whether the source recorded 1 or nothing (HF-050); every other setting is the source's.
    The job must sit in its dataset folder, where the fresh job goes."""
    job_dir = Path(os.path.abspath(job_dir))
    manifest_path = job_dir / "manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"no manifest.json in {job_dir}: --from-solved needs a finished job")
    manifest = json.loads(manifest_path.read_text())
    if not isinstance(manifest, dict):
        raise ValueError(f"{manifest_path} must hold a manifest object, got {type(manifest).__name__}")
    project = job_dir / "project" / Path(_manifest_field(manifest, "project", "job_copy")).name
    if not project.is_file():
        raise FileNotFoundError(f"the source job's solved project {project} does not exist")
    results = project.with_suffix(".aedtresults")
    if not results.is_dir():
        raise FileNotFoundError(f"the source job's solution folder {results} (.aedtresults) does not exist")
    _manifest_field(manifest, "project", "run_design")
    raw = manifest.get("config")
    if not isinstance(raw, Mapping):
        raise ValueError(f"{manifest_path} has no config block")
    raw = copy.deepcopy(dict(raw))
    output = raw.get("output")
    raw["output"] = {**(output if isinstance(output, Mapping) else {}), "schema": 2}
    cfg = build_config(raw)
    folder = f"{cfg.project.dataset_id}_{frequency_label(cfg.solver.frequency_ghz)}"
    if job_dir.parent.name != folder:
        raise ValueError(f"the source job {job_dir} is not in its dataset folder {folder}/, where the fresh job goes")
    record = manifest["project"].get("dataset_id_record") or {}
    if record.get("dataset_id") is not None and record["dataset_id"] != cfg.project.dataset_id:
        raise ValueError(f"the source manifest's config dataset_id {cfg.project.dataset_id!r} disagrees with its "
                         f"dataset_id record {record!r}")
    _check_source_records(cfg, manifest)
    return SourceJob(job_dir=job_dir, manifest=manifest, config=cfg, project=project, results=results,
                     manifest_sha256=sha256_file(manifest_path))


def reexport_job_dirs(source: SourceJob, job_id: str) -> JobDirs:
    """The fresh job's directories, in the source's dataset folder (``<dataset>_<label>/job_<job_id>``); an existing
    job directory, the source's own included, is refused."""
    return create_job_dirs(source.job_dir.parent.parent, source.config.project.dataset_id,
                           source.config.solver.frequency_ghz, job_id)


def stage_solved_job(source: SourceJob, dirs: JobDirs) -> StagedProject:
    """Copy the source's solved project and its ``.aedtresults`` into the fresh job, before any AEDT session, as
    ``hpc/diagnose_poynting.py`` copies a solved job; AEDT opens only the copy, so the source is never written."""
    job_project = dirs.project / source.project.name
    shutil.copy2(source.project, job_project)
    shutil.copytree(source.results, dirs.project / source.results.name)
    source_sha, job_sha = sha256_file(source.project), sha256_file(job_project)
    if job_sha != source_sha:
        raise RuntimeError(f"the copy {job_project} does not hash to its source {source.project}")
    record = source.manifest["project"].get("dataset_id_record") or {}
    return StagedProject(source=source.project, source_sha=source_sha, job_project=job_project, job_sha=job_sha,
                         recorded_id=record.get("dataset_id"),
                         record=Path(record["path"]) if record.get("path") else None)


def prepare_from_solved(source: SourceJob, dirs: JobDirs, job_id: str, hfss: Any, staged: StagedProject
                        ) -> SolvedContext:
    """Open the solved run design in the job's copy and rebuild the extraction context from its geometry, as
    ``prepare_and_solve`` builds it, without setting up or solving anything. Both polarizations must be solved
    already. ``staged`` is ``stage_solved_job``'s copy, which ``hfss`` opened. ``solve_s`` is None; the solve's mesh
    and facet settings are the source manifest's."""
    from bbsim.extract import ExtractionContext, describe_exit_face, exit_grid
    from bbsim.geometry import select_faces
    from bbsim.hfss_setup import (
        EXIT_CS,
        add_exit_flux_expression,
        add_outgoing_power_expression,
        aedt_messages,
        far_field_sphere,
        missing_polarizations,
        set_total_fields,
        solved_variations,
    )
    from bbsim.inventory import object_record
    from bbsim.naming import parametric_setup_name, plane_wave_name, setup_name
    from bbsim.session import aedt_versions

    cfg = source.config
    freq = cfg.solver.frequency_ghz
    run_design = source.manifest["project"]["run_design"]
    exc = cfg.excitation
    phi_values = incident_angle_values(exc.phi_lower_deg, exc.phi_upper_deg, exc.phi_step_deg)
    theta_values = incident_angle_values(exc.theta_lower_deg, exc.theta_upper_deg, exc.theta_step_deg)

    versions = aedt_versions(hfss)
    log.info("versions: %s", versions)
    if run_design not in hfss.design_list:
        raise RuntimeError(f"design {run_design!r} not in {staged.job_project}; designs: {list(hfss.design_list)}")
    if hfss.design_name != run_design:
        hfss.set_active_design(run_design)
    if hfss.modeler.model_units != "mm":
        raise RuntimeError(f"model units are {hfss.modeler.model_units!r}; the configuration assumes mm")
    solution = last_adaptive_solution(freq)
    variations = solved_variations(hfss, solution)
    log.info("solved variations of %s: %s", solution, variations)
    missing = missing_polarizations(variations, cfg.output.polarizations)
    if missing:
        raise RuntimeError(f"no solved variations for Ephi={missing} in {solution!r} of {staged.job_project} (copied "
                           f"from {source.job_dir}): AEDT holds {variations}, and a re-export never solves")

    entrance, exit_face = select_faces(hfss, cfg.geometry)
    _require_v2_exit_normal(cfg, exit_face)
    solid = object_record(hfss, cfg.geometry.object)
    log.info("entrance face %s, exit face %s, material %s", entrance, exit_face, solid["material"])
    # The field type is a post-processing setting: HPC step 8b set it on copies of these solved projects and then
    # evaluated their fields. The source run set the same before its solve.
    set_total_fields(hfss)
    # HF-052: both named expressions are re-defined exactly as a fresh run defines them, the clear in
    # add_outgoing_power_expression first. This supersedes the schema-2 design's "only exit_flux is added": every
    # solve since the schema-2 runner (b2af196) saves exit_flux in its run design, and what AEDT does when a name
    # is added twice is unverified. Named expressions are calculator definitions; they do not touch the solution.
    add_outgoing_power_expression(hfss)
    add_exit_flux_expression(hfss)
    sphere, grid = far_field_sphere(cfg.far_field, freq)   # the solved design holds this sphere (load_source_job)

    ctx = ExtractionContext(
        frequency_ghz=freq, solution=solution, sphere_name=sphere,
        phi_values=phi_values, theta_values=theta_values, grid=grid, exit_field=cfg.exit_field,
        exit_face=exit_face, exit_cs_name=EXIT_CS, exit_cs_x=cfg.geometry.exit_cs_x, exit_cs_y=cfg.geometry.exit_cs_y,
        scratch_dir=dirs.scratch, timeout_s=cfg.output.export_timeout_s,
        incoming_power_w=incoming_power_w(exc.ei_v_per_m, entrance.area_mm2),
    )
    face_geometry = describe_exit_face(hfss, ctx)
    lattice = exit_grid(face_geometry, cfg.exit_field) if cfg.exit_field.manual else None
    _require_v2_lattice(face_geometry.cross_section, lattice)
    log.info("exit lattice: %d points per incident angle (%s)", lattice.lattice_points, face_geometry.bounds_method)
    solver = source.manifest.get("solver") or {}
    mesh = solver.get("initial_mesh_settings") or {}
    return SolvedContext(
        ctx=ctx, face_geometry=face_geometry, lattice=lattice, entrance=entrance, exit_face=exit_face, solid=solid,
        setup=setup_name(freq), parametric=parametric_setup_name(freq), sphere=sphere, grid=grid,
        plane_wave=plane_wave_name(freq), wall_facets=solver.get("wall_facets"),
        mesh_inherited=mesh.get("inherited_from_base_design"), mesh_requested=mesh.get("requested"),
        variations=variations, messages=aedt_messages(hfss), solve_s=None, source=staged.source,
        source_sha=staged.source_sha, job_project=staged.job_project, job_sha=staged.job_sha,
        recorded_id=staged.recorded_id, record=staged.record, run_design=run_design, versions=versions,
    )


def run_from_solved(source: SourceJob, dirs: JobDirs, job_id: str, reference: dict | None = None) -> Path:
    """--from-solved: copy the source's solved project and results, open the copy, extract and write schema 2 with
    no solve. The convergence tables are exported fresh from the opened design and qualified against the source's
    target (CDX-025); the sidecar records ``provenance.reexport`` and the manifest ``solver.reused_solution``."""
    from bbsim.session import aedt_session, aedt_versions

    t_start = time.perf_counter()
    cfg = source.config
    check_before_aedt(cfg, reference)
    staged = stage_solved_job(source, dirs)
    reexport = {"source_job": str(source.job_dir), "source_manifest_sha256": source.manifest_sha256,
                "tool": REEXPORT_TOOL, "git_commit": git_commit()}
    with aedt_session(staged.job_project, source.manifest["project"]["run_design"], cfg.project.aedt_version) as hfss:
        _check_versions_for_reference(reference, aedt_versions(hfss))
        solved = prepare_from_solved(source, dirs, job_id, hfss, staged)
        return extract_and_write(cfg, dirs, job_id, hfss, solved, t_start, reference, reexport)


_FROM_SOLVED_EXCLUDES = (("--config", "config"), ("--frequency-ghz", "frequency_ghz"), ("--cores", "cores"),
                         ("--output-root", "output_root"), ("--project", "project"), ("--polarizations", "polarizations"),
                         ("--aedt-version", "aedt_version"))


def _main_from_solved(args: argparse.Namespace) -> int:
    given = [flag for flag, attr in _FROM_SOLVED_EXCLUDES if getattr(args, attr) is not None]
    if given:
        print(f"refused: --from-solved takes its configuration from the source job's manifest; drop {', '.join(given)}",
              file=sys.stderr)
        return 2
    if not args.job_id:
        print("refused: --from-solved needs --job-id, the fresh job beside the source (for example tight1_v2)",
              file=sys.stderr)
        return 2
    reference = None
    try:
        source = load_source_job(args.from_solved)
        if args.column_order_reference is not None:
            reference = load_reference(args.column_order_reference)
        check_before_aedt(source.config, reference)
        dirs = reexport_job_dirs(source, args.job_id)
    except (OSError, ValueError) as exc:
        print(f"refused before the job started: {exc}", file=sys.stderr)
        return 2
    _setup_logging(dirs.logs / "run.log", args.log_level)
    write_json(dirs.root / "effective_config.json", config_to_dict(source.config))
    log.info("job %s -> %s: schema-2 re-export of %s, no solve", args.job_id, dirs.root, source.job_dir)
    if reference is not None:
        log.info("column-order reference from %s (%s, %s)", args.column_order_reference, reference.get("source_job"),
                 reference.get("aedt_version"))
    try:
        manifest = run_from_solved(source, dirs, args.job_id, reference)
    except Exception:  # noqa: BLE001 - report and exit non-zero; the session is already released
        log.exception("re-export failed")
        return 1
    log.info("done; manifest at %s", manifest)
    return 0


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.from_solved is not None:
        return _main_from_solved(args)
    if args.config is None:
        print("refused: one of --config or --from-solved is required", file=sys.stderr)
        return 2
    cfg = load_config(args.config, overrides_from_args(args))
    if not cfg.project.path.is_file():
        print(f"prepared project not found: {cfg.project.path} (run prepare_hfss_project.py first, or pass --project)",
              file=sys.stderr)
        return 2
    reference = None
    try:
        if args.column_order_reference is not None:
            if cfg.output.schema != 2:
                raise ValueError(f"--column-order-reference applies to schema 2 only; {args.config} writes schema "
                                 f"{cfg.output.schema}")
            reference = load_reference(args.column_order_reference)
        check_before_aedt(cfg, reference)
    except (OSError, ValueError) as exc:
        print(f"refused before the job started: {exc}", file=sys.stderr)
        return 2
    job_id = args.job_id or default_job_id()
    dirs = create_job_dirs(cfg.output.root, cfg.project.dataset_id, cfg.solver.frequency_ghz, job_id)
    _setup_logging(dirs.logs / "run.log", args.log_level)
    write_json(dirs.root / "effective_config.json", config_to_dict(cfg))
    log.info("job %s -> %s (schema %d)", job_id, dirs.root, cfg.output.schema)
    if reference is not None:
        log.info("column-order reference from %s (%s, %s)", args.column_order_reference, reference.get("source_job"),
                 reference.get("aedt_version"))
    warning = allocation_warning(cfg.solver.cores)
    if warning:
        log.warning(warning)
    try:
        manifest = run_job(cfg, dirs, job_id, reference=reference)
    except Exception:  # noqa: BLE001 - report and exit non-zero; the session is already released
        log.exception("job failed")
        return 1
    log.info("done; manifest at %s", manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
