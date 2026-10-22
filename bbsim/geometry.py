"""Identify entrance/exit faces geometrically (never by stored face ID) and build the exit frame."""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from bbsim.config import FaceSelector, GeometryConfig

_AXIS_INDEX = {"x": 0, "y": 1, "z": 2}


class GeometryError(RuntimeError):
    """The model does not match the configured geometry."""


@dataclass(frozen=True)
class FaceInfo:
    id: int
    center_mm: tuple[float, float, float]
    area_mm2: float
    outward_normal: tuple[float, float, float]


def check_box_dimensions(bounding_box_mm: Sequence[float], expected_size_mm: Sequence[float] | None,
                         tol_mm: float = 1e-6) -> None:
    if expected_size_mm is None:
        return
    size = [bounding_box_mm[3 + i] - bounding_box_mm[i] for i in range(3)]
    if any(abs(s - e) > tol_mm for s, e in zip(size, expected_size_mm)):
        raise GeometryError(f"object size {size} mm does not match expected size {list(expected_size_mm)} mm")


def check_face_area(face: FaceInfo, expected_area_mm2: float | None, rel_tol: float = 1e-6) -> None:
    if expected_area_mm2 is None:
        return
    if abs(face.area_mm2 - expected_area_mm2) > rel_tol * expected_area_mm2:
        raise GeometryError(f"face {face.id} area {face.area_mm2} mm^2 does not match expected {expected_area_mm2} mm^2")


def select_face_on_bounding_plane(faces: Sequence[Mapping[str, Any]], bounding_box_mm: Sequence[float],
                                  selector: FaceSelector, tol_mm: float = 1e-6) -> FaceInfo:
    """Pick the unique planar face whose centre lies on the bounding-box plane ``axis = min|max``.

    Non-planar faces are never candidates. AEDT's GetFaceCenter fails on a curved face, and PyAEDT 1.7.0 then
    reports the centroid of the face's vertices if it has more than one (a face without AEDT vertices gets one
    per edge, at the edge's start), else the centroid of four samples on its first edge. The side of a cylinder
    from z = 0 to z = L, bounded by its end circles, thus reports a mid-length point such as (R, 0, L/2), on
    the plane x = R but on neither end plane. A curved face with at most one vertex whose first edge is an end
    circle would report that circle's centre, on the end plane.
    The outward normal is defined by construction (-axis for ``min``, +axis for ``max``); a normal
    reported by PyAEDT must agree with it.
    """
    i = _AXIS_INDEX[selector.axis]
    plane = bounding_box_mm[i] if selector.side == "min" else bounding_box_mm[3 + i]
    on_plane = [f for f in faces if abs(float(f["center_mm"][i]) - plane) <= tol_mm]
    curved = [f for f in on_plane if f.get("is_planar") is False]
    hits = [f for f in on_plane if f.get("is_planar") is not False]
    if len(hits) != 1:
        ignored = f"; ignored non-planar faces on that plane: ids {[f['id'] for f in curved]}" if curved else ""
        raise GeometryError(
            f"expected exactly one planar face with centre on {selector.axis} = {plane} mm ({selector.side}), "
            f"found {len(hits)}: ids {[h['id'] for h in hits]}{ignored}"
        )
    normal = [0.0, 0.0, 0.0]
    normal[i] = -1.0 if selector.side == "min" else 1.0
    face = hits[0]
    reported = face.get("normal")
    if reported is not None and sum(float(a) * b for a, b in zip(reported, normal)) < 1.0 - 1e-6:
        raise GeometryError(f"face {face['id']} normal {list(reported)} does not point along the outward "
                            f"{selector.side} {selector.axis} direction {normal}")
    return FaceInfo(int(face["id"]), tuple(float(c) for c in face["center_mm"]), float(face["area_mm2"]), tuple(normal))


def check_exit_frame(exit_face: FaceInfo, exit_cs_x: Sequence[float], tol: float = 1e-9) -> None:
    """The exit CS x axis must be the exit face's outward normal: the far-field hemisphere relies on it."""
    if any(abs(float(a) - b) > tol for a, b in zip(exit_cs_x, exit_face.outward_normal)):
        raise GeometryError(f"exit_cs_x {list(exit_cs_x)} must equal the exit face's outward normal "
                            f"{list(exit_face.outward_normal)}")


def check_pose(entrance: FaceInfo, exit_face: FaceInfo, pose: str) -> None:
    """Pose rule canonical-z: propagation along +z, the HFSS global frame in which BBRsim reads incident angles."""
    if pose == "unchecked":
        return
    if pose != "canonical-z":
        raise GeometryError(f"unknown pose rule {pose!r}")
    if tuple(entrance.outward_normal) != (0.0, 0.0, -1.0) or tuple(exit_face.outward_normal) != (0.0, 0.0, 1.0):
        raise GeometryError(
            "pose rule canonical-z needs the entrance outward normal (0, 0, -1) and the exit outward normal "
            f"(0, 0, 1); got entrance {list(entrance.outward_normal)}, exit {list(exit_face.outward_normal)}"
        )


def select_faces(hfss: Any, geom: GeometryConfig) -> tuple[FaceInfo, FaceInfo]:
    """Entrance and exit faces of ``geom.object`` in the active design, validated against the config."""
    from bbsim.inventory import object_record

    if geom.object not in hfss.modeler.object_names:
        raise GeometryError(f"object {geom.object!r} not in design; objects are {list(hfss.modeler.object_names)}")
    record = object_record(hfss, geom.object)
    check_box_dimensions(record["bounding_box_mm"], geom.expected_box_size_mm)
    entrance = select_face_on_bounding_plane(record["faces"], record["bounding_box_mm"], geom.entrance_face)
    exit_face = select_face_on_bounding_plane(record["faces"], record["bounding_box_mm"], geom.exit_face)
    if entrance.id == exit_face.id:
        raise GeometryError("entrance and exit selectors resolve to the same face")
    check_face_area(entrance, geom.expected_face_area_mm2)
    check_face_area(exit_face, geom.expected_face_area_mm2)
    check_exit_frame(exit_face, geom.exit_cs_x)
    check_pose(entrance, exit_face, geom.pose)
    return entrance, exit_face


def create_exit_coordinate_system(hfss: Any, exit_face: FaceInfo, x_pointing: Sequence[float],
                                  y_pointing: Sequence[float], name: str = "outgoing_cs") -> None:
    """Legacy create_local_coordinate_system: axis mode at the exit-face centre."""
    hfss.modeler.create_coordinate_system(
        origin=list(exit_face.center_mm),
        reference_cs="Global",
        name=name,
        mode="axis",
        x_pointing=list(x_pointing),
        y_pointing=list(y_pointing),
    )


def create_exit_face_list(hfss: Any, exit_face: FaceInfo, name: str = "outgoing") -> None:
    """Legacy: face list named ``outgoing`` used by the calculator and the infinite sphere."""
    hfss.modeler.create_face_list([exit_face.id], name=name)
