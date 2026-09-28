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


def select_face_on_bounding_plane(faces: Sequence[Mapping[str, Any]], bounding_box_mm: Sequence[float],
                                  selector: FaceSelector, tol_mm: float = 1e-6) -> FaceInfo:
    """Pick the unique face whose centre lies on the bounding-box plane ``axis = min|max``.

    The outward normal is defined by construction (-axis for ``min``, +axis for ``max``).
    """
    i = _AXIS_INDEX[selector.axis]
    plane = bounding_box_mm[i] if selector.side == "min" else bounding_box_mm[3 + i]
    hits = [f for f in faces if abs(float(f["center_mm"][i]) - plane) <= tol_mm]
    if len(hits) != 1:
        raise GeometryError(
            f"expected exactly one face with centre on {selector.axis} = {plane} mm ({selector.side}), "
            f"found {len(hits)}: ids {[h['id'] for h in hits]}"
        )
    normal = [0.0, 0.0, 0.0]
    normal[i] = -1.0 if selector.side == "min" else 1.0
    face = hits[0]
    return FaceInfo(int(face["id"]), tuple(float(c) for c in face["center_mm"]), float(face["area_mm2"]), tuple(normal))


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
