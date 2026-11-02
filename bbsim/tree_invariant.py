"""BBRsim's cross-frequency invariant of a dataset sidecar, with the standard library only (G4-012, S65, HF-057).

Every frequency of one dataset ID in a tree must carry the same frequency-independent physics, or BBRsim's reader
refuses the tree. ``invariant(sidecar)`` builds the representation its sidecar parser builds for that comparison, and
``invariant_diff`` compares two of them block by block as its reader does, in ``INVARIANT_BLOCKS`` order:

- ``schema``: the integer major of ``schema_version`` ("1" or "1.x" is 1, "2" or "2.x" is 2);
- ``excitation``: only ``incident_phi_deg`` and ``incident_theta_deg``;
- ``frames``: only ``hfss_global_axes_in_canonical`` and ``exit_cs_axes_in_canonical``;
- ``symmetry`` and ``boundaries``: their blocks, whose flags must be booleans and whose boundary names strings, as
  BBRsim's loader requires (BBR024), so a value it refuses is refused here even when every frequency agrees on it;
- ``geometry``: ``extent_mm``, and ``shape`` when the sidecar has one;
- ``modes``: the block without ``basis`` and ``propagating_count``;
- ``cross_section``: ``exit_field.cross_section``, which is not a top-level key of a sidecar.

Recorded-only fields (origins, pose text, face selectors, bounding boxes, provenance) are outside it, and so is
``exit_field.grid``, so a frequency whose exit step was fitted finer (1010 GHz on crack1, 136 x 51) agrees with the
others. Before anything is compared, every input the projection reads must exist with a usable type (``InvariantError``
names the first that does not), so a sidecar lacking a block never skips that block's comparison.

``close(a, b)`` is the reader's comparison: objects over the keys both carry, arrays element by element with equal
lengths, numbers to 1e-9 relative, strings and null exactly; and one rule the reader's Python mirror lacks: booleans are
equal only to booleans of the same value, so ``True`` and ``1`` differ. Shared by the assembly and ``--verify``; it
imports nothing outside the standard library.
"""
from __future__ import annotations

import copy
import math
from collections.abc import Mapping
from typing import Any

INVARIANT_BLOCKS = ("schema", "excitation", "frames", "symmetry", "boundaries", "geometry", "modes", "cross_section")
INVARIANT_REL = 1e-9
SYMMETRY_FLAGS = ("mirror_l", "mirror_g", "end_to_end", "rotational")
BOUNDARY_NAMES = ("entrance", "exit", "walls")
_DROPPED_MODES = ("basis", "propagating_count")


class InvariantError(ValueError):
    """A sidecar lacks an input of the projection, or holds it with an unusable type; the message names it."""


def _is_number(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _block(parent: Mapping, key: str, where: str) -> Mapping:
    value = parent.get(key) if isinstance(parent, Mapping) else None
    if not isinstance(value, Mapping):
        raise InvariantError(f"{where}{key} is missing or not an object (got {type(value).__name__})")
    return value


def _present(block: Mapping, keys: tuple[str, ...], where: str, kind: type, what: str) -> None:
    missing = [k for k in keys if k not in block]
    if missing:
        raise InvariantError(f"{where} lacks {', '.join(missing)}")
    for key in keys:
        if not isinstance(block[key], kind):
            raise InvariantError(f"{where}.{key} must be {what}, got {block[key]!r} (BBRsim's loader refuses it, "
                                 "BBR024)")


def _numbers(block: Mapping, key: str, where: str) -> list:
    value = block.get(key)
    if not (isinstance(value, list) and value and all(_is_number(v) for v in value)):
        raise InvariantError(f"{where}{key} must be a nonempty list of numbers, got {value!r}")
    return value


def _axes(frames: Mapping, key: str) -> Mapping:
    axes = _block(frames, key, "frames.")
    for name in ("x", "y", "z"):
        vector = axes.get(name)
        if not (isinstance(vector, list) and len(vector) == 3 and all(_is_number(v) for v in vector)):
            raise InvariantError(f"frames.{key}.{name} must be a 3-vector of numbers, got {vector!r}")
    return axes


def _finite(block: Mapping, keys: tuple[str, ...], where: str) -> None:
    for key in keys:
        value = block.get(key)
        if not (_is_number(value) and math.isfinite(value)):
            raise InvariantError(f"{where}{key} must be a finite number, got {value!r}")


def schema_major(sidecar: Mapping) -> int:
    """1 for ``schema_version`` "1" or "1.x", 2 for "2" or "2.x"; anything else raises ``InvariantError``."""
    version = sidecar.get("schema_version") if isinstance(sidecar, Mapping) else None
    if not isinstance(version, str):
        raise InvariantError(f"schema_version must be a string, got {version!r}")
    for major in (1, 2):
        if version == str(major) or version.startswith(f"{major}."):
            return major
    raise InvariantError(f"schema_version {version!r} is neither 1.x nor 2.x")


def invariant(sidecar: Mapping) -> dict:
    """The sidecar's cross-frequency invariant, as BBRsim's sidecar parser builds it (the module docstring); a deep
    copy, so editing it leaves the sidecar alone. Raises ``InvariantError`` for a missing or mistyped input."""
    if not isinstance(sidecar, Mapping):
        raise InvariantError(f"a sidecar is a JSON object, got {type(sidecar).__name__}")
    major = schema_major(sidecar)
    excitation = _block(sidecar, "excitation", "")
    phis = _numbers(excitation, "incident_phi_deg", "excitation.")
    thetas = _numbers(excitation, "incident_theta_deg", "excitation.")
    frames = _block(sidecar, "frames", "")
    hfss_axes, exit_axes = _axes(frames, "hfss_global_axes_in_canonical"), _axes(frames, "exit_cs_axes_in_canonical")
    symmetry = _block(sidecar, "symmetry", "")
    _present(symmetry, SYMMETRY_FLAGS, "symmetry", bool, "true or false")
    boundaries = _block(sidecar, "boundaries", "")
    _present(boundaries, BOUNDARY_NAMES, "boundaries", str, "a string")
    geometry = _block(sidecar, "geometry", "")
    extent = _block(geometry, "extent_mm", "geometry.")
    _finite(extent, ("p", "l", "g"), "geometry.extent_mm.")
    modes = _block(sidecar, "modes", "")
    if not isinstance(modes.get("mode"), str):
        raise InvariantError(f"modes.mode must be a string, got {modes.get('mode')!r}")
    _finite(modes, ("cutoff_ghz", "polarization_filter_limit_ghz"), "modes.")
    if major == 2 and not isinstance(modes.get("regime_cutoffs_ghz"), list):
        raise InvariantError("modes.regime_cutoffs_ghz must be a list under schema 2.x (OR-031)")
    section = _block(_block(sidecar, "exit_field", ""), "cross_section", "exit_field.")
    if not isinstance(section.get("shape"), str):
        raise InvariantError(f"exit_field.cross_section.shape must be a string, got {section.get('shape')!r}")
    shape = {"shape": geometry["shape"]} if "shape" in geometry else {}
    projected: dict[str, Any] = {
        "schema": major,
        "excitation": {"incident_phi_deg": phis, "incident_theta_deg": thetas},
        "frames": {"hfss_global_axes_in_canonical": hfss_axes, "exit_cs_axes_in_canonical": exit_axes},
        "symmetry": symmetry,
        "boundaries": boundaries,
        "geometry": {"extent_mm": extent, **shape},
        "modes": {k: v for k, v in modes.items() if k not in _DROPPED_MODES},
        "cross_section": section,
    }
    return copy.deepcopy(projected)


def close(a: Any, b: Any) -> bool:
    """BBRsim's comparison of invariant values (the module docstring), booleans told from numbers."""
    if isinstance(a, bool) or isinstance(b, bool):
        return isinstance(a, bool) and isinstance(b, bool) and a == b
    if _is_number(a) and _is_number(b):
        return abs(a - b) <= INVARIANT_REL * max(abs(a), abs(b))
    if isinstance(a, Mapping) and isinstance(b, Mapping):
        return all(close(a[k], b[k]) for k in a.keys() & b.keys())
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(close(x, y) for x, y in zip(a, b))
    if _is_number(a) or _is_number(b):
        return False
    return a == b


def invariant_diff(a: Mapping, b: Mapping) -> list[str]:
    """The ``INVARIANT_BLOCKS`` on which two projections (``invariant``) differ, in that order; empty when they
    agree."""
    return [k for k in INVARIANT_BLOCKS if k in a and k in b and not close(a[k], b[k])]
