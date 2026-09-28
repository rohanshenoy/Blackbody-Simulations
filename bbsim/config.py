"""Run configuration: TOML file -> validated frozen dataclasses, with dotted CLI overrides."""
from __future__ import annotations

import tomllib
from collections.abc import Mapping
from dataclasses import MISSING, asdict, dataclass, fields
from pathlib import Path
from typing import Any


class ConfigError(ValueError):
    """The configuration is missing, malformed, or violates a constraint."""


@dataclass(frozen=True)
class FaceSelector:
    axis: str
    side: str

    def __post_init__(self) -> None:
        if self.axis not in ("x", "y", "z"):
            raise ConfigError(f"face axis must be 'x', 'y' or 'z', got {self.axis!r}")
        if self.side not in ("min", "max"):
            raise ConfigError(f"face side must be 'min' or 'max', got {self.side!r}")


@dataclass(frozen=True)
class ProjectConfig:
    path: Path
    design: str
    dataset_id: str
    aedt_version: str = "2025.2"


@dataclass(frozen=True)
class GeometryConfig:
    object: str
    entrance_face: FaceSelector
    exit_face: FaceSelector
    exit_cs_x: tuple[float, float, float]
    exit_cs_y: tuple[float, float, float]
    expected_box_size_mm: tuple[float, float, float] | None = None
    expected_face_area_mm2: float | None = None

    def __post_init__(self) -> None:
        for name in ("exit_cs_x", "exit_cs_y"):
            if len(getattr(self, name)) != 3:
                raise ConfigError(f"{name} must have 3 components")
        if self.expected_box_size_mm is not None and len(self.expected_box_size_mm) != 3:
            raise ConfigError("expected_box_size_mm must have 3 components")


@dataclass(frozen=True)
class ExcitationConfig:
    ei_v_per_m: float
    theta_lower_deg: float
    theta_upper_deg: float
    theta_step_deg: float
    phi_lower_deg: float
    phi_upper_deg: float
    phi_step_deg: float

    def __post_init__(self) -> None:
        if self.ei_v_per_m != 1.0:
            raise ConfigError(
                "ei_v_per_m must be 1.0: the HFSS plane wave amplitude is fixed at 1 V/m by the polarization "
                "expressions ['Ephi', '1-Ephi'] (as in the reference run); Ei only scales IngoingPower, so any "
                "other value would silently rescale every transmission ratio"
            )
        for name in ("theta", "phi"):
            lo, hi, step = (getattr(self, f"{name}_{k}_deg") for k in ("lower", "upper", "step"))
            if hi < lo:
                raise ConfigError(f"{name}_upper_deg must not be below {name}_lower_deg")
            if step <= 0:
                raise ConfigError(f"{name}_step_deg must be positive")


@dataclass(frozen=True)
class SolverConfig:
    frequency_ghz: float
    max_delta_e: float
    max_passes: int
    cores: int
    sweep: str = "discrete"

    def __post_init__(self) -> None:
        if self.frequency_ghz <= 0:
            raise ConfigError("frequency_ghz must be positive")
        if self.max_passes < 1:
            raise ConfigError("max_passes must be at least 1")
        if self.cores < 1:
            raise ConfigError("cores must be at least 1")
        if self.sweep != "discrete":
            raise ConfigError(
                f"sweep={self.sweep!r} is not supported; only 'discrete' is ported. "
                "Adaptive refinement remains in legacy/bbsim1freq.py."
            )


@dataclass(frozen=True)
class FarFieldConfig:
    theta_lower_deg: float
    theta_upper_deg: float
    phi_lower_deg: float
    phi_upper_deg: float
    a_mm: float
    b_mm: float
    fineness: float
    min_coarseness_deg: float
    max_coarseness_deg: float

    def __post_init__(self) -> None:
        if self.a_mm <= 0 or self.b_mm <= 0 or self.fineness <= 0:
            raise ConfigError("a_mm, b_mm and fineness must be positive")
        if not 0 < self.min_coarseness_deg <= self.max_coarseness_deg:
            raise ConfigError("need 0 < min_coarseness_deg <= max_coarseness_deg")
        for name in ("theta", "phi"):
            if getattr(self, f"{name}_upper_deg") < getattr(self, f"{name}_lower_deg"):
                raise ConfigError(f"{name}_upper_deg must not be below {name}_lower_deg")


@dataclass(frozen=True)
class ExitFieldConfig:
    manual: bool = True
    resolution_mm: tuple[float, float, float] = (0.0, 0.1, 0.001)
    boundary_mm: tuple[tuple[float, float, float], tuple[float, float, float]] | None = None

    def __post_init__(self) -> None:
        res = self.resolution_mm
        if len(res) != 3 or any(v < 0 for v in res) or not any(v > 0 for v in res):
            raise ConfigError(f"resolution_mm must be 3 non-negative steps with at least one positive, got {list(res)}")
        if self.boundary_mm is not None:
            ok = (len(self.boundary_mm) == 2 and all(len(corner) == 3 for corner in self.boundary_mm)
                  and all(hi >= lo for lo, hi in zip(*self.boundary_mm)))
            if not ok:
                raise ConfigError("boundary_mm must be [[xmin, ymin, zmin], [xmax, ymax, zmax]] with max >= min")


@dataclass(frozen=True)
class OutputConfig:
    root: Path
    polarizations: tuple[int, ...] = (0, 1)
    export_timeout_s: float = 300.0

    def __post_init__(self) -> None:
        if not self.polarizations or any(p not in (0, 1) for p in self.polarizations):
            raise ConfigError(f"polarizations must be a non-empty subset of [0, 1], got {list(self.polarizations)}")
        if len(set(self.polarizations)) != len(self.polarizations):
            raise ConfigError(f"polarizations must not repeat, got {list(self.polarizations)}")
        if self.export_timeout_s <= 0:
            raise ConfigError("export_timeout_s must be positive")


@dataclass(frozen=True)
class RunConfig:
    project: ProjectConfig
    geometry: GeometryConfig
    excitation: ExcitationConfig
    solver: SolverConfig
    far_field: FarFieldConfig
    exit_field: ExitFieldConfig
    output: OutputConfig


_SECTIONS: dict[str, type] = {
    "project": ProjectConfig,
    "geometry": GeometryConfig,
    "excitation": ExcitationConfig,
    "solver": SolverConfig,
    "far_field": FarFieldConfig,
    "exit_field": ExitFieldConfig,
    "output": OutputConfig,
}
_PATH_FIELDS = {"path", "root"}
_FACE_FIELDS = {"entrance_face", "exit_face"}


def _coerce(name: str, value: Any, where: str) -> Any:
    if name in _PATH_FIELDS:
        if not isinstance(value, str):
            raise ConfigError(f"[{where}] {name} must be a string path")
        return Path(value).expanduser()
    if name in _FACE_FIELDS:
        if not isinstance(value, Mapping):
            raise ConfigError(f"[{where}] {name} must be a table like {{ axis = 'z', side = 'min' }}")
        return _build(FaceSelector, value, f"{where}.{name}")
    if isinstance(value, list):
        return tuple(tuple(v) if isinstance(v, list) else v for v in value)
    return value


def _build(cls: type, data: Mapping[str, Any], where: str) -> Any:
    names = {f.name for f in fields(cls)}
    unknown = sorted(set(data) - names)
    if unknown:
        raise ConfigError(f"[{where}] unknown keys: {unknown}")
    kwargs: dict[str, Any] = {}
    for f in fields(cls):
        if f.name in data:
            kwargs[f.name] = _coerce(f.name, data[f.name], where)
        elif f.default is MISSING and f.default_factory is MISSING:
            raise ConfigError(f"[{where}] missing required key {f.name!r}")
    try:
        return cls(**kwargs)
    except ConfigError as e:
        raise ConfigError(f"[{where}] {e}") from None
    except TypeError as e:
        raise ConfigError(f"[{where}] {e}") from None


def build_config(raw: Mapping[str, Any], overrides: Mapping[str, Any] | None = None) -> RunConfig:
    data: dict[str, dict[str, Any]] = {k: dict(v) for k, v in raw.items()}
    unknown_sections = sorted(set(data) - set(_SECTIONS))
    if unknown_sections:
        raise ConfigError(f"unknown sections: {unknown_sections}")
    for key, value in (overrides or {}).items():
        if value is None:
            continue
        section, _, name = key.partition(".")
        if section not in _SECTIONS or not name:
            raise ConfigError(f"unknown override {key!r}")
        if name not in {f.name for f in fields(_SECTIONS[section])}:
            raise ConfigError(f"unknown override {key!r}")
        data.setdefault(section, {})[name] = value
    return RunConfig(**{s: _build(cls, data.get(s, {}), s) for s, cls in _SECTIONS.items()})


def load_config(path: Path | str, overrides: Mapping[str, Any] | None = None) -> RunConfig:
    with open(path, "rb") as fh:
        raw = tomllib.load(fh)
    return build_config(raw, overrides)


def config_to_dict(cfg: RunConfig) -> dict[str, Any]:
    """JSON-serializable copy: Paths become strings, tuples become lists."""
    def convert(obj: Any) -> Any:
        if isinstance(obj, Path):
            return str(obj)
        if isinstance(obj, dict):
            return {k: convert(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [convert(v) for v in obj]
        return obj
    return convert(asdict(cfg))
