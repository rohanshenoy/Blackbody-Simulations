"""Run configuration: TOML file -> validated frozen dataclasses, with dotted CLI overrides."""
from __future__ import annotations

import math
import os
import tomllib
from collections.abc import Mapping
from dataclasses import MISSING, asdict, dataclass, fields
from pathlib import Path
from typing import Any

from bbsim import acceptance_gate


class ConfigError(ValueError):
    """The configuration is missing, malformed, or violates a constraint."""


POSE_RULES = ("canonical-z", "unchecked")
SYMMETRIES = ("mirror_l", "mirror_g", "end_to_end", "rotational")


def _finite(name: str, value: Any) -> float:
    """A real, finite number (CDX-024): nan passes every inequality guard and inf passes "> 0", so both are refused
    here, before any job directory or AEDT session exists. Booleans and strings are refused too. The value itself
    is left as given (an integer frequency keeps its label)."""
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigError(f"{name} must be a finite number, got {value!r}")
    return float(value)


def _finite_vector(name: str, values: Any) -> tuple[float, ...]:
    if isinstance(values, (str, bytes)) or not hasattr(values, "__iter__"):
        raise ConfigError(f"{name} must be a list of numbers, got {values!r}")
    return tuple(_finite(f"{name}[{i}]", v) for i, v in enumerate(values))


def _integer(name: str, value: Any) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ConfigError(f"{name} must be an integer, got {value!r}")
    return value


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
    pose: str = "canonical-z"
    symmetry: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        for name in ("exit_cs_x", "exit_cs_y"):
            if len(_finite_vector(name, getattr(self, name))) != 3:
                raise ConfigError(f"{name} must have 3 components")
        x = [float(c) for c in self.exit_cs_x]
        y = [float(c) for c in self.exit_cs_y]
        for name, v in (("exit_cs_x", x), ("exit_cs_y", y)):
            if abs(math.sqrt(sum(c * c for c in v)) - 1.0) > 1e-9:
                raise ConfigError(f"{name} must be a unit vector, got {v}")
        if abs(sum(a * b for a, b in zip(x, y))) > 1e-9:
            raise ConfigError(f"exit_cs_x and exit_cs_y must be orthogonal, got {x} and {y}")
        if self.expected_box_size_mm is not None and len(_finite_vector("expected_box_size_mm",
                                                                         self.expected_box_size_mm)) != 3:
            raise ConfigError("expected_box_size_mm must have 3 components")
        if self.expected_face_area_mm2 is not None:
            _finite("expected_face_area_mm2", self.expected_face_area_mm2)
        if self.pose not in POSE_RULES:
            raise ConfigError(f"pose must be one of {list(POSE_RULES)}, got {self.pose!r}")
        if isinstance(self.symmetry, str):
            raise ConfigError(f"symmetry must be a list such as [\"mirror_l\"], got the string {self.symmetry!r}")
        unknown = [s for s in self.symmetry if s not in SYMMETRIES]
        if unknown:
            raise ConfigError(f"symmetry values must be from {list(SYMMETRIES)}, got {unknown}")
        if len(set(self.symmetry)) != len(self.symmetry):
            raise ConfigError(f"symmetry must not repeat, got {list(self.symmetry)}")


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
        for f in fields(self):
            _finite(f.name, getattr(self, f.name))
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
    # A curved wall's facets: the normal deviation of a manual surface approximation on the geometry object
    # (Rohan, 2026-10-06, ledger HF-030). None keeps HFSS's own faceting, as every config did before.
    wall_normal_deviation_deg: float | None = None

    def __post_init__(self) -> None:
        _finite("frequency_ghz", self.frequency_ghz)
        _finite("max_delta_e", self.max_delta_e)
        _integer("max_passes", self.max_passes)
        _integer("cores", self.cores)
        if self.wall_normal_deviation_deg is not None:
            _finite("wall_normal_deviation_deg", self.wall_normal_deviation_deg)
        if self.frequency_ghz <= 0:
            raise ConfigError("frequency_ghz must be positive")
        if self.max_delta_e <= 0:
            raise ConfigError("max_delta_e must be positive")
        if self.wall_normal_deviation_deg is not None and not 0 < self.wall_normal_deviation_deg < 90:
            raise ConfigError("wall_normal_deviation_deg must lie strictly between 0 and 90 degrees")
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
        for f in fields(self):
            _finite(f.name, getattr(self, f.name))
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
    edge_samples: int = 64

    def __post_init__(self) -> None:
        if self.manual is not True:
            raise ConfigError(
                f"manual={self.manual!r} is not supported; only manual = true (ExportOnGrid on the exit-frame "
                "lattice) is ported. The mesh export (CalculatorWrite) writes points in HFSS global coordinates "
                "with no lattice, so its CSVs and sidecar would fail BBRsim's exit-plane and grid checks; it "
                "remains in legacy/bbsim1freq.py."
            )
        res = _finite_vector("resolution_mm", self.resolution_mm)
        if len(res) != 3 or any(v < 0 for v in res) or not any(v > 0 for v in res):
            raise ConfigError(f"resolution_mm must be 3 non-negative steps with at least one positive, got {list(res)}")
        if self.boundary_mm is not None:
            for i, corner in enumerate(self.boundary_mm):
                _finite_vector(f"boundary_mm[{i}]", corner)
            ok = (len(self.boundary_mm) == 2 and all(len(corner) == 3 for corner in self.boundary_mm)
                  and all(hi >= lo for lo, hi in zip(*self.boundary_mm)))
            if not ok:
                raise ConfigError("boundary_mm must be [[xmin, ymin, zmin], [xmax, ymax, zmax]] with max >= min")
        if isinstance(self.edge_samples, bool) or not isinstance(self.edge_samples, int) or self.edge_samples < 8:
            raise ConfigError(f"edge_samples must be an integer >= 8, got {self.edge_samples!r}")


SCHEMAS = (1, 2)


@dataclass(frozen=True)
class OutputConfig:
    root: Path
    polarizations: tuple[int, ...] = (0, 1)
    export_timeout_s: float = 300.0
    # The dataset schema the runner writes (HF-050): 2 for every new dataset; 1 keeps the 1.0 writer until the _v2
    # trees pass BBRsim's tree check, and for diagnostics whose grid schema 2 cannot describe.
    schema: int = 2

    def __post_init__(self) -> None:
        _finite("export_timeout_s", self.export_timeout_s)
        if isinstance(self.schema, bool) or not isinstance(self.schema, int) or self.schema not in SCHEMAS:
            raise ConfigError(f"schema must be the integer 1 or 2, got {self.schema!r}")
        if not self.polarizations or any(p not in (0, 1) for p in self.polarizations):
            raise ConfigError(f"polarizations must be a non-empty subset of [0, 1], got {list(self.polarizations)}")
        if len(set(self.polarizations)) != len(self.polarizations):
            raise ConfigError(f"polarizations must not repeat, got {list(self.polarizations)}")
        if self.export_timeout_s <= 0:
            raise ConfigError("export_timeout_s must be positive")


@dataclass(frozen=True)
class AcceptanceConfig:
    """The acceptance the job's dataset needs before assembly (CDX-039, OR-042): ``kind`` is "none",
    "grid_acceptance" (a 35-key solve judged against the 15-key ``reference_job`` of the same geometry, at
    ``t_rel_tol``, default the checker's 0.005) or "reexport_check" (set by --from-solved from its source job, never
    in a config file). The manifest records it, and assembly requires exactly it."""
    kind: str
    reference_job: Path | None = None
    t_rel_tol: float | None = None

    def __post_init__(self) -> None:
        try:
            acceptance_gate.validate_requirement_fields(self.kind, self.reference_job, self.t_rel_tol)
        except acceptance_gate.GateError as exc:
            raise ConfigError(str(exc)) from None


@dataclass(frozen=True)
class RunConfig:
    project: ProjectConfig
    geometry: GeometryConfig
    excitation: ExcitationConfig
    solver: SolverConfig
    far_field: FarFieldConfig
    exit_field: ExitFieldConfig
    output: OutputConfig
    acceptance: AcceptanceConfig


_SECTIONS: dict[str, type] = {
    "project": ProjectConfig,
    "geometry": GeometryConfig,
    "excitation": ExcitationConfig,
    "solver": SolverConfig,
    "far_field": FarFieldConfig,
    "exit_field": ExitFieldConfig,
    "output": OutputConfig,
    "acceptance": AcceptanceConfig,
}
_PATH_FIELDS = {"path", "root"}
_FACE_FIELDS = {"entrance_face", "exit_face"}


def _coerce(name: str, value: Any, where: str) -> Any:
    if name == "reference_job":
        # Optional: None stays None. A path is made absolute from the loading directory (the checkout, for the
        # bridge scripts), not the TOML's, and kept as written otherwise; comparisons resolve both sides (CDX-039).
        if value is None:
            return None
        if not isinstance(value, str) or not value.strip():
            raise ConfigError(f"[{where}] reference_job must be a nonblank string path, got {value!r}")
        return Path(os.path.abspath(os.path.expanduser(value)))
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
    if "acceptance" not in data:
        raise ConfigError('[acceptance] is required: kind = "none", or kind = "grid_acceptance" with reference_job '
                          "(CDX-039, OR-042)")
    return RunConfig(**{s: _build(cls, data.get(s, {}), s) for s, cls in _SECTIONS.items()})


def load_config(path: Path | str, overrides: Mapping[str, Any] | None = None) -> RunConfig:
    with open(path, "rb") as fh:
        raw = tomllib.load(fh)
    cfg = build_config(raw, overrides)
    if cfg.acceptance.kind == "reexport_check":
        raise ConfigError("[acceptance] kind reexport_check: a re-export's requirement is set by --from-solved from its "
                          "source job; a config file declares none or grid_acceptance")
    return cfg


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
