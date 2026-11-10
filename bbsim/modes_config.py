"""Config for the ports-only mode study of a hollow round bore (ledger OR-053). Pure: no AEDT.

Sections: [project] geometry_spec, [wall], [ports], [solver], [output]. The bore's material, dataset ID, design
and object come from the geometry spec, whose cylinder already refuses a non-vacuum bore (CDX-029).
"""
from __future__ import annotations

import math
import tomllib
from collections.abc import Mapping
from dataclasses import MISSING, dataclass, fields
from pathlib import Path
from typing import Any

from bbsim.geometry_spec import GeometrySpec, load_geometry_spec
from bbsim.modes import cutoff_distance, propagating_count


class ModesConfigError(ValueError):
    """The modes configuration is missing, malformed, or violates a constraint."""


WALL_KINDS = ("pec", "finite_conductivity")
MIN_CUTOFF_DISTANCE = 0.01
EXTRA_MODES = 2


def _finite(name: str, value: Any) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ModesConfigError(f"{name} must be a finite number, got {value!r}")
    return float(value)


def _positive(name: str, value: Any) -> float:
    v = _finite(name, value)
    if not v > 0:
        raise ModesConfigError(f"{name} must be > 0, got {value!r}")
    return v


def _integer(name: str, value: Any) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ModesConfigError(f"{name} must be an integer, got {value!r}")
    return value


@dataclass(frozen=True)
class WallConfig:
    kind: str
    conductivity_s_per_m: float | None = None
    thickness_um: float | None = None

    def __post_init__(self) -> None:
        if self.kind not in WALL_KINDS:
            raise ModesConfigError(f"kind must be one of {list(WALL_KINDS)}, got {self.kind!r}")
        if self.kind == "pec":
            for name in ("conductivity_s_per_m", "thickness_um"):
                if getattr(self, name) is not None:
                    raise ModesConfigError(f"{name} is not allowed with kind = 'pec'")
            return
        if self.conductivity_s_per_m is None:
            raise ModesConfigError("conductivity_s_per_m is required with kind = 'finite_conductivity'")
        _positive("conductivity_s_per_m", self.conductivity_s_per_m)
        if self.thickness_um is not None:
            _positive("thickness_um", self.thickness_um)


@dataclass(frozen=True)
class PortsConfig:
    modes: int
    frequencies_ghz: tuple[float, ...]

    def __post_init__(self) -> None:
        if _integer("modes", self.modes) < 1:
            raise ModesConfigError(f"modes must be >= 1, got {self.modes!r}")
        if not isinstance(self.frequencies_ghz, tuple) or not self.frequencies_ghz:
            raise ModesConfigError("frequencies_ghz must be a non-empty list of numbers")
        for i, f in enumerate(self.frequencies_ghz):
            _positive(f"frequencies_ghz[{i}]", f)
        if len(set(self.frequencies_ghz)) != len(self.frequencies_ghz):
            raise ModesConfigError(f"frequencies_ghz repeats a frequency: {list(self.frequencies_ghz)}")


@dataclass(frozen=True)
class ModesSolverConfig:
    cores: int = 4
    curvilinear: bool = True
    wall_normal_deviation_deg: float = 5.0
    port_accuracy_percent: int = 1
    port_face_max_length_mm: float | None = None

    def __post_init__(self) -> None:
        if _integer("cores", self.cores) < 1:
            raise ModesConfigError(f"cores must be >= 1, got {self.cores!r}")
        if not isinstance(self.curvilinear, bool):
            raise ModesConfigError(f"curvilinear must be true or false, got {self.curvilinear!r}")
        dev = _finite("wall_normal_deviation_deg", self.wall_normal_deviation_deg)
        if not 0 < dev < 90:
            raise ModesConfigError(f"wall_normal_deviation_deg must be in (0, 90), got {dev!r}")
        if _integer("port_accuracy_percent", self.port_accuracy_percent) < 1:
            raise ModesConfigError(f"port_accuracy_percent must be >= 1, got {self.port_accuracy_percent!r}")
        if self.port_face_max_length_mm is not None:
            _positive("port_face_max_length_mm", self.port_face_max_length_mm)


@dataclass(frozen=True)
class ModesOutputConfig:
    root: Path


@dataclass(frozen=True)
class ModesConfig:
    geometry_spec: Path
    spec: GeometrySpec
    wall: WallConfig
    ports: PortsConfig
    solver: ModesSolverConfig
    output: ModesOutputConfig

    @property
    def dataset_id(self) -> str:
        return self.spec.project.dataset_id

    @property
    def design(self) -> str:
        return self.spec.project.design

    @property
    def object(self) -> str:
        return self.spec.project.object

    @property
    def radius_m(self) -> float:
        return float(self.spec.solid.radius_mm) * 1e-3


_SECTIONS: dict[str, type] = {"wall": WallConfig, "ports": PortsConfig, "solver": ModesSolverConfig,
                              "output": ModesOutputConfig}


def _build(cls: type, data: Mapping[str, Any], where: str) -> Any:
    names = {f.name for f in fields(cls)}
    unknown = sorted(set(data) - names)
    if unknown:
        raise ModesConfigError(f"[{where}] unknown keys: {unknown}")
    kwargs: dict[str, Any] = {}
    for f in fields(cls):
        if f.name in data:
            value = data[f.name]
            if f.name == "root":
                if not isinstance(value, str):
                    raise ModesConfigError(f"[{where}] root must be a string path")
                value = Path(value).expanduser()
            elif isinstance(value, list):
                value = tuple(value)
            kwargs[f.name] = value
        elif f.default is MISSING:
            raise ModesConfigError(f"[{where}] missing required key {f.name!r}")
    try:
        return cls(**kwargs)
    except ModesConfigError as e:
        raise ModesConfigError(f"[{where}] {e}") from None
    except TypeError as e:
        raise ModesConfigError(f"[{where}] {e}") from None


def load_modes_config(path: Path | str) -> ModesConfig:
    path = Path(path)
    with open(path, "rb") as fh:
        raw = tomllib.load(fh)
    unknown = sorted(set(raw) - {"project", *_SECTIONS})
    if unknown:
        raise ModesConfigError(f"unknown sections: {unknown}")
    for section in ("project", *_SECTIONS):
        if section not in raw:
            raise ModesConfigError(f"missing section [{section}]")
    project = dict(raw["project"])
    extra = sorted(set(project) - {"geometry_spec"})
    if extra:
        raise ModesConfigError(f"[project] unknown keys: {extra}")
    if not isinstance(project.get("geometry_spec"), str) or not project["geometry_spec"].strip():
        raise ModesConfigError("[project] geometry_spec must be a path string")
    spec_path = Path(project["geometry_spec"]).expanduser()
    if not spec_path.is_absolute():
        spec_path = path.parent / spec_path
    spec = load_geometry_spec(spec_path)
    parts = {s: _build(cls, dict(raw[s]), s) for s, cls in _SECTIONS.items()}
    ports: PortsConfig = parts["ports"]
    radius_m = float(spec.solid.radius_mm) * 1e-3
    top = max(ports.frequencies_ghz)
    need = propagating_count(radius_m, top) + EXTRA_MODES
    if ports.modes < need:
        raise ModesConfigError(f"[ports] modes = {ports.modes} is below {need}: {propagating_count(radius_m, top)} "
                               f"modes propagate at {top:g} GHz and the port needs {EXTRA_MODES} more")
    for f in ports.frequencies_ghz:
        d = cutoff_distance(radius_m, f)
        if d < MIN_CUTOFF_DISTANCE:
            raise ModesConfigError(f"[ports] frequencies_ghz: {f:g} GHz is {d:.4f} of an ideal cutoff, closer "
                                   f"than {MIN_CUTOFF_DISTANCE}")
    return ModesConfig(geometry_spec=spec_path, spec=spec, **parts)
