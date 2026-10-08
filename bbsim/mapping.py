"""Old-to-new design/object names for project cleanup, validated against a source inventory."""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class DesignRename:
    source_design: str
    source_object: str
    new_design: str
    new_object: str
    dataset_id: str
    gap_um: float


@dataclass(frozen=True)
class DesignMapping:
    renames: tuple[DesignRename, ...]

    @property
    def design_map(self) -> dict[str, str]:
        return {r.source_design: r.new_design for r in self.renames}

    @property
    def object_map(self) -> dict[str, dict[str, str]]:
        return {r.source_design: {r.source_object: r.new_object} for r in self.renames}


def load_mapping(path: Path) -> DesignMapping:
    data = json.loads(Path(path).read_text())
    renames = tuple(DesignRename(**entry) for entry in data["renames"])
    if not renames:
        raise ValueError(f"{path}: no renames")
    ids = [r.dataset_id for r in renames]
    repeated = sorted({i for i in ids if ids.count(i) > 1})
    if repeated:
        # Two designs sharing a dataset_id would write into one dataset directory and tree entry.
        raise ValueError(f"{path}: dataset_id repeated across renames: {repeated}")
    return DesignMapping(renames)


def validate_mapping(mapping: DesignMapping, source_inventory: dict[str, Any], tol_um: float = 0.5) -> list[str]:
    """Every source design/object must exist; the micrometre label must equal the smallest box dimension."""
    errors: list[str] = []
    designs = source_inventory["designs"]
    new_names = [r.new_design for r in mapping.renames]
    if len(set(new_names)) != len(new_names):
        errors.append(f"new design names are not unique: {new_names}")
    for r in mapping.renames:
        design = designs.get(r.source_design)
        if design is None:
            errors.append(f"design {r.source_design!r} not found in source project")
            continue
        objects = {o["name"]: o for o in design["objects"]}
        obj = objects.get(r.source_object)
        if obj is None:
            errors.append(f"design {r.source_design!r}: object {r.source_object!r} not found; objects are {sorted(objects)}")
            continue
        bb = obj["bounding_box_mm"]
        smallest_um = 1000.0 * min(bb[3 + i] - bb[i] for i in range(3))
        if abs(smallest_um - r.gap_um) > tol_um:
            errors.append(f"design {r.source_design!r}: gap_um={r.gap_um} but smallest box dimension is {smallest_um:.3f} um")
    return errors


def plan_deletions(mapping: DesignMapping, source_inventory: dict[str, Any]) -> list[str]:
    retained = set(mapping.design_map)
    return sorted(name for name in source_inventory["designs"] if name not in retained)
