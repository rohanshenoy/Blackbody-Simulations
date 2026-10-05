"""Geometry spec: TOML -> validated dataclasses for building a bare base design from scratch. Pure: no AEDT.

A cylinder spec fixes the canonical pose by construction: axis along +z from z = 0 to z = length, centred
on x = y = 0. A "cad" kind will add a file and its two end planes later, through the same loader.
"""
from __future__ import annotations

import math
import re
import tomllib
from collections.abc import Mapping
from dataclasses import MISSING, dataclass, fields
from pathlib import Path
from typing import Any


class GeometrySpecError(ValueError):
    """The geometry spec is missing, malformed, or violates a constraint."""


DATASET_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]*$")


def _non_empty(where: str, name: str, value: Any) -> None:
    if not isinstance(value, str) or not value.strip():
        raise GeometrySpecError(f"[{where}] {name} must be a non-empty string, got {value!r}")


@dataclass(frozen=True)
class SpecProject:
    output: Path
    design: str
    object: str
    dataset_id: str
    aedt_version: str = "2025.2"

    def __post_init__(self) -> None:
        for name in ("design", "object", "dataset_id", "aedt_version"):
            _non_empty("project", name, getattr(self, name))
        if (not DATASET_ID_PATTERN.match(self.dataset_id) or "GHz" in self.dataset_id
                or "Ephi=" in self.dataset_id):
            raise GeometrySpecError(
                f"[project] dataset_id {self.dataset_id!r} must match {DATASET_ID_PATTERN.pattern} and contain "
                "neither 'GHz' nor 'Ephi=' (BBRsim parses <id>_<freq>GHz_Ephi=<n>)")
        if self.output.suffix != ".aedt":
            raise GeometrySpecError(f"[project] output must end in .aedt, got {str(self.output)!r}")


@dataclass(frozen=True)
class CylinderSolid:
    kind: str
    material: str
    radius_mm: float
    length_mm: float

    def __post_init__(self) -> None:
        _non_empty("solid", "material", self.material)
        for name in ("radius_mm", "length_mm"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not value > 0:
                raise GeometrySpecError(f"[solid] {name} must be a positive number, got {value!r}")

    @property
    def expected_bounding_box_mm(self) -> list[float]:
        r, length = float(self.radius_mm), float(self.length_mm)
        return [-r, -r, 0.0, r, r, length]

    @property
    def end_face_area_mm2(self) -> float:
        return math.pi * float(self.radius_mm) ** 2

    @property
    def face_count(self) -> int:
        return 3


SOLIDS: dict[str, type] = {"cylinder": CylinderSolid}


@dataclass(frozen=True)
class GeometrySpec:
    project: SpecProject
    solid: CylinderSolid


def _build(cls: type, data: Mapping[str, Any], where: str) -> Any:
    names = {f.name for f in fields(cls)}
    unknown = sorted(set(data) - names)
    if unknown:
        raise GeometrySpecError(f"[{where}] unknown keys: {unknown}")
    kwargs: dict[str, Any] = {}
    for f in fields(cls):
        if f.name in data:
            value = data[f.name]
            if f.name == "output":
                if not isinstance(value, str):
                    raise GeometrySpecError(f"[{where}] output must be a string path")
                value = Path(value).expanduser()
            kwargs[f.name] = value
        elif f.default is MISSING:
            raise GeometrySpecError(f"[{where}] missing required key {f.name!r}")
    try:
        return cls(**kwargs)
    except TypeError as exc:
        raise GeometrySpecError(f"[{where}] {exc}") from None


def build_geometry_spec(raw: Mapping[str, Any]) -> GeometrySpec:
    unknown = sorted(set(raw) - {"project", "solid"})
    if unknown:
        raise GeometrySpecError(f"unknown sections: {unknown}")
    for section in ("project", "solid"):
        if section not in raw:
            raise GeometrySpecError(f"missing section [{section}]")
    solid_raw = dict(raw["solid"])
    kind = solid_raw.get("kind")
    if kind not in SOLIDS:
        raise GeometrySpecError(f"[solid] kind must be one of {sorted(SOLIDS)}, got {kind!r}")
    return GeometrySpec(project=_build(SpecProject, dict(raw["project"]), "project"),
                        solid=_build(SOLIDS[kind], solid_raw, "solid"))


def load_geometry_spec(path: Path | str) -> GeometrySpec:
    with open(path, "rb") as fh:
        return build_geometry_spec(tomllib.load(fh))


def plan_text(spec: GeometrySpec) -> str:
    """What the build will create and check, for --dry-run and the start of every build log."""
    p, s = spec.project, spec.solid
    bbox = ", ".join(f"{v:g}" for v in s.expected_bounding_box_mm)
    return "\n".join([
        f"Geometry spec -> new project {p.output}",
        f"  design {p.design}, object {p.object}, dataset_id {p.dataset_id}, AEDT {p.aedt_version}",
        f"  solid: {s.kind}, material {s.material}, radius {s.radius_mm:g} mm, length {s.length_mm:g} mm",
        f"  pose: axis +z from z = 0 to z = {s.length_mm:g} mm, centred on x = y = 0 (canonical-z)",
        f"  checks: bounding box [{bbox}] mm; {s.face_count} faces; one planar end face on z = 0 and one on "
        f"z = {s.length_mm:g} mm, area {s.end_face_area_mm2:.9g} mm^2 each; no boundaries, excitations or setups",
    ])
