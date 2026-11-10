"""AEDT calls of the ports-only mode study of a hollow round bore (ledger OR-053).

Driven Modal solution type, the lossy wall, the port-face mesh, two N-mode wave ports, one ports-only setup per
frequency, the solution lookup, the ports' Gamma (AEDT 2025 R2 scripting guide 28-18; PyAEDT 1.7.0 names) and a
Touchstone export. PyAEDT is never imported here: every call goes through the ``hfss`` handle.
"""
from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any

from bbsim.hfss_setup import assign_wall_facets, initial_mesh_settings_kwargs
from bbsim.modes import GammaRow
from bbsim.naming import frequency_label, setup_name

WALL_NAME = "wall"
PORT_FACES = "port_faces"
FREQUENCY_RTOL = 1e-9


def set_driven_modal(hfss: Any) -> str:
    """Select SolutionsHfss.DrivenModal and return what the design reads back; anything but "Modal" is refused."""
    hfss.solution_type = "Modal"
    got = hfss.solution_type
    if str(got) != "Modal":
        raise RuntimeError(f"the design's solution type reads {got!r} after setting 'Modal'")
    return str(got)


def assign_wall(hfss: Any, side_face_id: int, wall: Any) -> dict[str, Any]:
    """Assign the bore's side face its wall. PEC is the default boundary, so nothing is assigned."""
    if wall.kind == "pec":
        return {"kind": "pec", "assigned": None}
    kwargs: dict[str, Any] = dict(assignment=[side_face_id], conductivity=wall.conductivity_s_per_m,
                                  use_thickness=wall.thickness_um is not None, name=WALL_NAME)
    if wall.thickness_um is not None:
        kwargs["thickness"] = f"{wall.thickness_um}um"
    boundary = hfss.assign_finite_conductivity(**kwargs)
    if not boundary:
        raise RuntimeError(f"assign_finite_conductivity returned {boundary!r}")
    out: dict[str, Any] = {"kind": "finite_conductivity", "name": str(getattr(boundary, "name", WALL_NAME)),
                           "conductivity_s_per_m": wall.conductivity_s_per_m, "thickness_um": wall.thickness_um}
    props = getattr(boundary, "props", None)
    if isinstance(props, dict):
        out["properties"] = {str(k): _plain(v) for k, v in props.items()}
    return out


def _plain(value: Any) -> Any:
    return value if isinstance(value, (bool, int, float, str, type(None))) else str(value)


def assign_mesh(hfss: Any, cfg: Any, end_face_ids: list[int]) -> dict[str, Any]:
    """Reference slider settings (curvilinear from the config), faceted wall, and a length mesh on the port faces."""
    requested = dict(initial_mesh_settings_kwargs(), curvilinear=cfg.solver.curvilinear)
    if not hfss.mesh.assign_initial_mesh_from_slider(**requested):
        raise RuntimeError("assign_initial_mesh_from_slider returned False")
    facets = assign_wall_facets(hfss, cfg.object, cfg.solver.wall_normal_deviation_deg)
    out: dict[str, Any] = {"initial_mesh_settings": requested, "wall_facets": facets, "port_faces": None}
    length = cfg.solver.port_face_max_length_mm
    if end_face_ids and length is not None:
        # maximum_elements=None: PyAEDT's default cap of 1000 elements would stop a port face short of the length.
        op = hfss.mesh.assign_length_mesh(list(end_face_ids), inside_selection=False, maximum_length=f"{length}mm",
                                          maximum_elements=None, name=PORT_FACES)
        if not op:
            raise RuntimeError(f"assign_length_mesh returned {op!r}")
        out["port_faces"] = {"faces": list(end_face_ids), "maximum_length_mm": length,
                             "name": str(getattr(op, "name", PORT_FACES))}
    return out


def assign_ports(hfss: Any, entrance_face_id: int, exit_face_id: int, modes: int) -> list[str]:
    """Two wave ports with ``modes`` modes each, not renormalized, no de-embedding: port1 at z = 0, port2 at the far end."""
    names = ["port1", "port2"]
    for face, name in zip((entrance_face_id, exit_face_id), names):
        port = hfss.wave_port(assignment=face, modes=modes, name=name, renormalize=False, deembed=0)
        if not port:
            raise RuntimeError(f"wave_port returned {port!r} for {name} on face {face}")
    return names


def insert_ports_only_setup(hfss: Any, frequency_ghz: float, port_accuracy_percent: int) -> str:
    """A driven setup that solves the ports only (no 3-D field solve) at one frequency."""
    name = setup_name(frequency_ghz)
    hfss.odesign.GetModule("AnalysisSetup").InsertSetup(
        "HfssDriven", ["NAME:" + name, "Frequency:=", frequency_label(frequency_ghz), "PortsOnly:=", True,
                       "PortAccuracy:=", port_accuracy_percent])
    return name


def solution_of(hfss: Any, setup: str) -> str:
    """The one valid solution of ``setup`` ("<setup> : <sweep>")."""
    listed = [str(s) for s in hfss.odesign.GetModule("Solutions").GetValidISolutionList()]
    mine = [s for s in listed if s.startswith(setup + " :")]
    if len(mine) != 1:
        raise RuntimeError(f"expected exactly one solution of setup {setup!r}, found {len(mine)}; "
                           f"the design's solutions are {listed}")
    return mine[0]


def _port_number(name: str) -> int:
    m = re.fullmatch(r"port(\d+)", name)
    if m is None:
        raise RuntimeError(f"port name {name!r} is not of the form 'port<N>'")
    return int(m.group(1))


def read_gamma(hfss: Any, setup: str, ports: list[str], modes: int, frequency_ghz: float) -> list[GammaRow]:
    """Gamma(port:m) of every port mode: rows ordered by port, then mode 1..modes (alpha = Re, beta = Im)."""
    expressions = [(_port_number(p), m, f"Gamma({p}:{m})") for p in ports for m in range(1, modes + 1)]
    data = hfss.post.get_solution_data(expressions=[e for _, _, e in expressions],
                                       setup_sweep_name=solution_of(hfss, setup))
    swept = [float(x) for x in data.primary_sweep_values]
    if len(swept) != 1 or not math.isclose(swept[0], frequency_ghz, rel_tol=FREQUENCY_RTOL, abs_tol=0.0):
        raise RuntimeError(f"setup {setup!r} solved at {swept} GHz, not the requested {frequency_ghz} GHz")
    rows: list[GammaRow] = []
    for port, mode, expr in expressions:
        values = []
        for formula in ("real", "imag"):
            x, y = data.get_expression_data(expr, formula=formula)
            if len(y) != 1 or len(x) != 1:
                raise RuntimeError(f"AEDT returned no {formula} data for {expr} in setup {setup!r} "
                                   f"({len(y)} values)")
            if not math.isclose(float(x[0]), frequency_ghz, rel_tol=FREQUENCY_RTOL, abs_tol=0.0):
                raise RuntimeError(f"{expr} was read at {float(x[0])} GHz, not {frequency_ghz} GHz")
            values.append(float(y[0]))
        rows.append(GammaRow(frequency_ghz, port, mode, values[0], values[1]))
    return rows


def export_touchstone_with_gamma(hfss: Any, setup: str, path: Path) -> str | None:
    """Evidence for the first HPC run: a Touchstone export that also carries Gamma. None on success, else the error text."""
    try:
        hfss.odesign.GetModule("Solutions").ExportNetworkData(
            "", [solution_of(hfss, setup)], 3, str(path), ["all"], False, 50, "S", -1, 1, 15, False, True)
    except Exception as exc:  # noqa: BLE001 - reported as evidence, not raised
        return f"{type(exc).__name__}: {exc}"
    return None
