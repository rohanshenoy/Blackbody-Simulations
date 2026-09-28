"""Inventory designs, objects, faces, boundaries and setups; diff retained designs after cleanup."""
from __future__ import annotations

from typing import Any


def _native_pairs(flat: Any) -> list[dict[str, str]]:
    """BoundarySetup.GetBoundaries()/GetExcitations() return [name, type, name, type, ...]."""
    items = [str(x) for x in (flat or [])]
    return [{"name": items[i].split(":")[0], "type": items[i + 1]} for i in range(0, len(items) - 1, 2)]


def object_record(hfss: Any, object_name: str) -> dict[str, Any]:
    obj = hfss.modeler[object_name]
    faces = []
    for face in obj.faces:
        try:
            normal = [float(v) for v in face.normal] if face.normal else None
        except Exception:  # noqa: BLE001 - normal is informational only
            normal = None
        faces.append({
            "id": int(face.id),
            "center_mm": [float(c) for c in face.center],
            "area_mm2": float(hfss.modeler.get_face_area(face.id)),
            "normal": normal,
        })
    return {
        "name": object_name,
        "material": str(obj.material_name),
        "bounding_box_mm": [float(v) for v in obj.bounding_box],
        "faces": sorted(faces, key=lambda f: f["id"]),
    }


def _mesh_settings(hfss: Any) -> Any:
    try:
        return dict(hfss.mesh.initial_mesh_settings.props)
    except Exception as exc:  # noqa: BLE001 - informational only
        return f"unavailable: {exc}"


def inventory_design(hfss: Any) -> dict[str, Any]:
    odesign = hfss.odesign
    boundary_module = odesign.GetModule("BoundarySetup")
    analysis_module = odesign.GetModule("AnalysisSetup")
    variables = {str(name): str(odesign.GetVariableValue(name)) for name in odesign.GetVariables()}
    return {
        "name": str(hfss.design_name),
        "solution_type": str(hfss.solution_type),
        "units": str(hfss.modeler.model_units),
        "variables": variables,
        "objects": [object_record(hfss, name) for name in hfss.modeler.object_names],
        "boundaries": _native_pairs(boundary_module.GetBoundaries()),
        "excitations": _native_pairs(boundary_module.GetExcitations()),
        "setups": [str(s) for s in analysis_module.GetSetups()],
        "coordinate_systems": [str(cs.name) for cs in hfss.modeler.coordinate_systems],
        "initial_mesh_settings": _mesh_settings(hfss),
    }


def inventory_project(hfss: Any) -> dict[str, Any]:
    designs: dict[str, Any] = {}
    for name in list(hfss.design_list):
        hfss.set_active_design(name)
        designs[str(name)] = inventory_design(hfss)
    return {"project_file": str(hfss.project_file), "project_name": str(hfss.project_name), "designs": designs}


def _face_signature(faces: list[dict[str, Any]]) -> list[tuple]:
    return sorted((tuple(round(c, 9) for c in f["center_mm"]), round(f["area_mm2"], 9)) for f in faces)


def compare_retained_designs(source: dict[str, Any], prepared: dict[str, Any], design_map: dict[str, str],
                             object_map: dict[str, dict[str, str]], tol_mm: float = 1e-9) -> list[str]:
    """Differences between each source design and its renamed copy. Face IDs are ignored; geometry is not."""
    diffs: list[str] = []
    for old, new in design_map.items():
        s, p = source["designs"].get(old), prepared["designs"].get(new)
        if s is None:
            diffs.append(f"source design {old!r} missing")
            continue
        if p is None:
            diffs.append(f"prepared design {new!r} missing")
            continue
        label = f"{old} -> {new}"
        for key in ("solution_type", "units", "variables", "boundaries", "excitations", "setups", "initial_mesh_settings"):
            if s.get(key) != p.get(key):
                diffs.append(f"{label}: {key} differ: {s.get(key)!r} vs {p.get(key)!r}")
        renames = object_map.get(old, {})
        s_objs = {renames.get(o["name"], o["name"]): o for o in s["objects"]}
        p_objs = {o["name"]: o for o in p["objects"]}
        if set(s_objs) != set(p_objs):
            diffs.append(f"{label}: objects differ: {sorted(s_objs)} vs {sorted(p_objs)}")
        for name in sorted(set(s_objs) & set(p_objs)):
            so, po = s_objs[name], p_objs[name]
            if so["material"] != po["material"]:
                diffs.append(f"{label}: object {name} material {so['material']!r} vs {po['material']!r}")
            if max(abs(a - b) for a, b in zip(so["bounding_box_mm"], po["bounding_box_mm"])) > tol_mm:
                diffs.append(f"{label}: object {name} bounding box {so['bounding_box_mm']} vs {po['bounding_box_mm']}")
            if _face_signature(so["faces"]) != _face_signature(po["faces"]):
                diffs.append(f"{label}: object {name} face centres/areas differ")
    return diffs
