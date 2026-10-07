"""HFSS design setup, replicating legacy/bbsim1freq.py call for call.

The *_arguments/*_kwargs builders are pure; the functions taking ``hfss`` need a live AEDT session.
"""
from __future__ import annotations

import inspect
from collections.abc import Collection
from pathlib import Path
from typing import Any

from bbsim.config import ExcitationConfig, FarFieldConfig
from bbsim.geometry import FaceInfo
from bbsim.naming import (
    frequency_label,
    parametric_setup_name,
    plane_wave_name,
    radiation_sphere_name,
    setup_name,
)
from bbsim.sampling import FarFieldGrid, far_field_grid, incident_point_count

EXIT_FACE_LIST = "outgoing"
EXIT_CS = "outgoing_cs"
OUTGOING_POWER_EXPRESSION = "outgoing_power"


def setup_arguments(frequency_ghz: float, max_delta_e: float, max_passes: int) -> list:
    return [f"NAME:{setup_name(frequency_ghz)}", "Frequency:=", frequency_label(frequency_ghz),
            "MaxDeltaE:=", max_delta_e, "MaximumPasses:=", max_passes]


def parametric_arguments(frequency_ghz: float) -> list:
    return [
        f"NAME:{parametric_setup_name(frequency_ghz)}",
        "IsEnabled:=", True,
        ["NAME:ProdOptiSetupDataV2", "SaveFields:=", True],
        "Sim. Setups:=", [setup_name(frequency_ghz)],
        ["NAME:Sweeps", ["NAME:SweepDefinition", "Variable:=", "Ephi", "Data:=", "LIN 0 1 1"]],
    ]


def plane_wave_arguments(entrance: FaceInfo, exc: ExcitationConfig, frequency_ghz: float) -> dict[str, Any]:
    """Keyword arguments for Hfss.plane_wave; polarization is symbolic in the design variable Ephi."""
    theta_num = incident_point_count(exc.theta_lower_deg, exc.theta_upper_deg, exc.theta_step_deg)
    phi_num = incident_point_count(exc.phi_lower_deg, exc.phi_upper_deg, exc.phi_step_deg)
    return dict(
        assignment=[entrance.id],
        vector_format="Spherical",
        origin=list(entrance.center_mm),
        polarization=["Ephi", "1-Ephi"],
        propagation_vector=[
            [f"{float(exc.phi_lower_deg)}deg", f"{float(exc.phi_upper_deg)}deg", phi_num],
            [f"{float(exc.theta_lower_deg)}deg", f"{float(exc.theta_upper_deg)}deg", theta_num],
        ],
        name=plane_wave_name(frequency_ghz),
    )


def far_field_sphere_kwargs(accepted: Collection[str], grid: FarFieldGrid, name: str) -> dict[str, Any]:
    """Keyword arguments for Hfss.insert_infinite_sphere, for whichever PyAEDT naming is installed.

    PyAEDT >= 1.0 uses theta_*/phi_*; the release the legacy script used had x_*/y_*, where for the
    Theta-Phi definition x is theta and y is phi (legacy passed x_start=rad_theta_lower).
    """
    common = dict(definition="Theta-Phi", custom_radiation_faces=EXIT_FACE_LIST,
                  custom_coordinate_system=EXIT_CS, name=name)
    theta = (grid.theta_start_deg, grid.theta_stop_deg, grid.theta_step_deg)
    phi = (grid.phi_start_deg, grid.phi_stop_deg, grid.phi_step_deg)
    if {"theta_start", "theta_stop", "theta_step", "phi_start", "phi_stop", "phi_step"} <= set(accepted):
        return {**common, **dict(zip(("theta_start", "theta_stop", "theta_step"), theta)),
                **dict(zip(("phi_start", "phi_stop", "phi_step"), phi))}
    if {"x_start", "x_stop", "x_step", "y_start", "y_stop", "y_step"} <= set(accepted):
        return {**common, **dict(zip(("x_start", "x_stop", "x_step"), theta)),
                **dict(zip(("y_start", "y_stop", "y_step"), phi))}
    raise RuntimeError(f"unrecognised Hfss.insert_infinite_sphere signature: {sorted(accepted)}")


def initial_mesh_settings_kwargs() -> dict[str, Any]:
    """Global mesh settings of the script-generated reference designs (crack1Rohan_500GHz, crack2_500GHz).

    The GUI-made base designs carry ``UseAutoLength=false`` with a 0.1 mm defeature length, which a
    duplicated design would inherit; the legacy script created a fresh design and therefore solved with
    ``UseAutoLength=true``. These kwargs reproduce the reference ``MeshSettings`` block exactly.
    """
    return dict(level=5, method="Auto", dynamic_surface=False, flex_mesh=False, curvilinear=False,
                fallback=True, phi=True, auto_model_resolution=True)


WALL_FACETS = "wall_facets"


def assign_wall_facets(hfss: Any, object_name: str, normal_deviation_deg: float) -> dict[str, Any]:
    """Facet the object's curved surfaces to ``normal_deviation_deg`` with a manual surface approximation.

    HFSS otherwise meshes the round gap's 50 um circle as a coarse polygon: HPC step 7a (2026-10-06) found its
    wall up to 0.76 um inside the circle. 5 deg gives 72 facets around a circle, the wall within 0.05 um of it.
    Rohan, 2026-10-06 (ledger HF-030): measured against the default before step 8.
    """
    op = hfss.mesh.assign_surface_mesh_manual(assignment=[object_name], surface_deviation=None,
                                              normal_dev=f"{normal_deviation_deg}deg", aspect_ratio=None,
                                              name=WALL_FACETS)
    if not op:
        raise RuntimeError(f"assign_surface_mesh_manual returned {op!r} for {object_name!r}")
    # PyAEDT returns the operation object whether or not AEDT created it, so read the design's mesh tree:
    # a silent failure would solve the "fine facets" run on default facets and fake an HF-030 ratio of 1.
    name = str(getattr(op, "name", WALL_FACETS))
    present = [str(n) for n in hfss.mesh.meshoperation_names]
    if name not in present:
        raise RuntimeError(f"mesh operation {name!r} is not in the design after assign_surface_mesh_manual; "
                           f"its mesh operations are {present}")
    return {"object": object_name, "normal_deviation_deg": normal_deviation_deg, "mesh_operation": name}


def reset_initial_mesh_settings(hfss: Any) -> tuple[Any, dict[str, Any]]:
    """Apply the reference global mesh settings; return (inherited settings, requested kwargs) for the manifest."""
    try:
        inherited = dict(hfss.mesh.initial_mesh_settings.props)
    except Exception as exc:  # noqa: BLE001 - informational only
        inherited = f"unavailable: {exc}"
    requested = initial_mesh_settings_kwargs()
    if not hfss.mesh.assign_initial_mesh_from_slider(**requested):
        raise RuntimeError("assign_initial_mesh_from_slider returned False")
    return inherited, requested


def initialize_variables(hfss: Any, ei_v_per_m: float) -> None:
    hfss["Ephi"] = 0
    hfss.variable_manager.set_variable("Ei", str(ei_v_per_m), sweep=False)


def clear_boundaries_and_excitations(hfss: Any) -> None:
    module = hfss.odesign.GetModule("BoundarySetup")
    module.DeleteAllBoundaries()
    module.DeleteAllExcitations()


def assign_radiation_boundaries(hfss: Any, entrance: FaceInfo, exit_face: FaceInfo) -> None:
    hfss.assign_radiation_boundary_to_faces(assignment=[entrance.id], name="rbin")
    hfss.assign_radiation_boundary_to_faces(assignment=[exit_face.id], name="rbout")


def add_outgoing_power_expression(hfss: Any) -> None:
    """|Re(Poynting)| integrated over the exit face list."""
    fields = hfss.odesign.GetModule("FieldsReporter")
    fields.ClearAllNamedExpr()
    fields.CalcStack("Clear")
    fields.CopyNamedExprToStack("Vector_RealPoynting")
    fields.CalcOp("Mag")
    fields.EnterSurf(EXIT_FACE_LIST)
    fields.CalcOp("Integrate")
    fields.AddNamedExpression(OUTGOING_POWER_EXPRESSION, "Fields")


def create_setup_and_polarization_sweep(hfss: Any, frequency_ghz: float, max_delta_e: float,
                                        max_passes: int) -> tuple[str, str]:
    hfss.odesign.GetModule("AnalysisSetup").InsertSetup("HfssDriven", setup_arguments(frequency_ghz, max_delta_e, max_passes))
    hfss.odesign.GetModule("Optimetrics").InsertSetup("OptiParametric", parametric_arguments(frequency_ghz))
    return setup_name(frequency_ghz), parametric_setup_name(frequency_ghz)


def insert_far_field_sphere(hfss: Any, ff: FarFieldConfig, frequency_ghz: float) -> tuple[str, FarFieldGrid]:
    grid = far_field_grid(ff.theta_lower_deg, ff.theta_upper_deg, ff.phi_lower_deg, ff.phi_upper_deg,
                          ff.a_mm, ff.b_mm, frequency_ghz, ff.fineness, ff.min_coarseness_deg, ff.max_coarseness_deg)
    name = radiation_sphere_name(frequency_ghz)
    accepted = inspect.signature(hfss.insert_infinite_sphere).parameters
    hfss.insert_infinite_sphere(**far_field_sphere_kwargs(accepted, grid, name))
    return name, grid


def create_plane_wave(hfss: Any, entrance: FaceInfo, exc: ExcitationConfig, frequency_ghz: float) -> str:
    kwargs = plane_wave_arguments(entrance, exc, frequency_ghz)
    hfss.plane_wave(**kwargs)
    return kwargs["name"]


def set_total_fields(hfss: Any) -> None:
    hfss.odesign.GetModule("Solutions").EditSources(["FieldType:=", "TotalFields"])


def solve(hfss: Any, cores: int) -> bool:
    return bool(hfss.analyze(cores=cores))


def aedt_messages(hfss: Any) -> list[str]:
    """Every entry of AEDT's message manager for the open design (all severities), or one line saying why not.

    The solver's own reasons (licence, validation, mesh) live only here; PyAEDT's analyze() does not
    surface them and reports "solved correctly" from the elapsed time alone.
    """
    try:
        return [str(m) for m in hfss.odesktop.GetMessages(hfss.project_name, hfss.design_name, 0)]
    except Exception as exc:  # noqa: BLE001 - diagnostics must not fail the run
        return [f"(AEDT message manager unavailable: {type(exc).__name__}: {exc})"]


def solved_variations(hfss: Any, solution: str) -> list[str]:
    """AEDT's solved-variation strings for a solution, e.g. ``["Ephi='0'", "Ephi='1'"]``; empty when nothing was solved."""
    return [str(v) for v in hfss.odesign.GetModule("Solutions").GetAvailableVariations(solution)]


def missing_polarizations(variations: Collection[str], polarizations: Collection[int]) -> list[int]:
    """Requested Ephi values with no solved variation (``Ephi='n'`` appears in none of the strings)."""
    return [int(p) for p in polarizations if not any(f"Ephi='{int(p)}'" in v for v in variations)]


def export_convergence_text(hfss: Any, setup: str, variation: str, path: Path) -> str:
    """Convergence table for one variation (e.g. ``"Ephi='0'"``); never fails the run."""
    try:
        hfss.export_convergence(setup, variation, str(path))
        return Path(path).read_text()
    except Exception as exc:  # noqa: BLE001
        return f"unavailable: {exc}"
