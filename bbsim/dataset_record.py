"""The HFSS dataset sidecar, schema 1.0 (agreed with the BBRSimulation session on 2026-10-04) and schema 2.0.

Built only from the run manifest as written to disk, so the sidecar cannot disagree with it. Pure:
no AEDT; scipy is imported only for the schema-1 circular-guide Bessel zeros (the 2.0 path needs none).
Canonical labels p, l, g (propagation, long, gap) are the Geant4 crack-local x, y, z and equal the exit-CS
axes X, Y, Z. Schema 2.0 (the sidecar field list of ledger HF-047 and HF-048) is built from
the same shared core as 1.0, so the blocks the two schemas have in common cannot drift apart.
"""
from __future__ import annotations

import copy
import math
from collections.abc import Mapping
from typing import Any

import numpy as np

from bbsim.sampling import SPEED_OF_LIGHT_M_PER_S as C
from bbsim.schema import FAR_FIELD_COLUMNS, WAVEGUIDE_COLUMNS, WAVEGUIDE_COLUMNS_V2

SCHEMA_VERSION = "1.0"
LIST_LIMIT_GHZ = 20000.0   # top of BBRsim's Planck band; fixed, so the modes block agrees across frequencies
TIE_REL = 1e-12
CANONICAL = "p,l,g; p x l = g; Geant4 crack-local (x,y,z) = (p,l,g)"
SYMMETRY_FLAGS = ("mirror_l", "mirror_g", "end_to_end", "rotational")
TRANSMITTANCE = {
    "definition": "outgoing_power_w / incoming_power_w, first row per key",
    "outgoing_power": "integral |Re S| over exit face",
    "incoming_includes_cos_theta": False,
}
_SHAPES = {"rectangle": "box", "disc": "cylinder"}

SCHEMA_VERSION_V2 = "2.0"
MODEL_KIND_V2 = "screen_free_abc"
CONVENTIONS_V2 = {"time_dependence": "exp(+j omega t)", "phasor": "peak"}
PLANE_TOLERANCE_M_V2 = 1e-12   # X = 0 on every row; BBRsim's tolerance (CDX-023)
REEXPORT_TOOL = "bbsim reexport_v2"
# Circular-guide cutoff zeros, to the digits the schema-2.0 field list fixes (HF-047, OR-031; no scipy on the 2.0 path):
X11P = 1.8411837813   # first zero of J1'  -> TE11
X01 = 2.4048255577    # first zero of J0   -> TM01
X21P = 3.0542369282   # first zero of J2'  -> TE21
_CONVERGENCE_FIELDS = ("achieved", "target", "passes_completed", "max_passes", "converged", "reasons")
_REEXPORT_FIELDS = ("source_job", "source_manifest_sha256", "git_commit")
_EXIT_FIELD_BLOCKS = ("column_order_check", "tangential_components", "quadrature", "scale", "sampler")
_CROSS_TERM_REL_TOL = 1e-6   # BBRsim's tolerance for valid_below_ghz against its re-derived limit


def mode_label(family: str, i: int, j: int) -> str:
    """``TE10``; a comma once either index reaches 10 (``TE0,10``, ``TE18,1``)."""
    return f"{family}{i}{j}" if i < 10 and j < 10 else f"{family}{i},{j}"


def _ordered(entries: list[tuple]) -> list[tuple]:
    """Ascending cutoff; ties within TIE_REL relative ordered TE before TM, then first, then second index.

    Entries are (cutoff_ghz, family, i, j, ...). Ties are mathematical (x'_0p = x_1p; TE_mn = TM_mn), so
    their order must not depend on the last bit a library returns.
    """
    def tie_key(entry):
        return (entry[1] != "TE", entry[2], entry[3])

    out: list[tuple] = []
    group: list[tuple] = []
    for entry in sorted(entries, key=lambda e: e[0]):
        if group and abs(entry[0] - group[0][0]) > TIE_REL * group[0][0]:
            out += sorted(group, key=tie_key)
            group = []
        group.append(entry)
    return out + sorted(group, key=tie_key)


def _rect_cutoff_ghz(a_m: float, b_m: float, m: int, n: int) -> float:
    return C / 2 * math.hypot(m / a_m, n / b_m) / 1e9


def rectangle_mode_count(a_m: float, b_m: float, limit_ghz: float) -> int:
    """TE (m, n) != (0, 0) plus TM m, n >= 1 index pairs with cutoff <= limit."""
    count = 0
    for m in range(int(2 * a_m * limit_ghz * 1e9 / C) + 2):
        for n in range(int(2 * b_m * limit_ghz * 1e9 / C) + 2):
            if (m, n) != (0, 0) and _rect_cutoff_ghz(a_m, b_m, m, n) <= limit_ghz:
                count += 2 if m >= 1 and n >= 1 else 1
    return count


def rectangle_modes(a_m: float, b_m: float, frequency_ghz: float, limit_ghz: float = LIST_LIMIT_GHZ) -> dict:
    """Closed ideal-PEC rectangular guide; a along l (index m) and b along g (index n), taken literally."""
    lowest = (1, 0) if a_m >= b_m else (0, 1)
    onsets = [{"n": 0, "mode": "TE10", "cutoff_ghz": _rect_cutoff_ghz(a_m, b_m, 1, 0)}]
    k = 1
    while _rect_cutoff_ghz(a_m, b_m, 0, k) <= limit_ghz:
        onsets.append({"n": k, "mode": mode_label("TE", 0, k), "cutoff_ghz": _rect_cutoff_ghz(a_m, b_m, 0, k)})
        k += 1
    return {
        "mode": mode_label("TE", *lowest),
        "cutoff_ghz": _rect_cutoff_ghz(a_m, b_m, *lowest),
        "basis": (f"closed ideal-PEC rectangular guide, a = {a_m!r} m along l, b = {b_m!r} m along g; "
                  "f_mn = (c/2) sqrt((m/a)^2 + (n/b)^2); no TEM mode"),
        "list_limit_ghz": limit_ghz,
        "polarization_filter_limit_ghz": _rect_cutoff_ghz(a_m, b_m, 0, 1),
        "gap_family_onsets": onsets,
        "mode_count_below_limit": rectangle_mode_count(a_m, b_m, limit_ghz),
        "propagating_count": rectangle_mode_count(a_m, b_m, frequency_ghz),
    }


def disc_entries(radius_m: float, limit_ghz: float) -> list[tuple[float, str, int, int, int]]:
    """Ordered (cutoff_ghz, family, n, p, degeneracy) with cutoff <= limit: TE_np from J_n' zeros, TM_np from J_n."""
    from scipy.special import jn_zeros, jnp_zeros

    scale = C / (2 * math.pi * radius_m) / 1e9
    x_max = limit_ghz / scale
    entries: list[tuple[float, str, int, int, int]] = []
    n = 0
    while True:
        k = int(x_max / math.pi) + 3
        found = False
        for family, zeros in (("TE", jnp_zeros(n, k)), ("TM", jn_zeros(n, k))):
            for p, x in enumerate(zeros, 1):
                cutoff = float(x) * scale
                if cutoff <= limit_ghz:
                    found = True
                    entries.append((cutoff, family, n, p, 1 if n == 0 else 2))
        # Order 0 can be empty while TE11 lies below the limit (x'_11 < j_01 < x'_01). From n = 1 on, the
        # smallest zero x'_n1 grows with n, so the first empty order ends the search.
        if not found and n >= 1:
            break
        n += 1
    return _ordered(entries)


def disc_modes(radius_m: float, frequency_ghz: float, limit_ghz: float = LIST_LIMIT_GHZ) -> dict:
    """Closed ideal-PEC circular guide: TE_np = x'_np c/(2 pi R), TM_np = x_np c/(2 pi R)."""
    from scipy.special import jnp_zeros

    cutoffs = [{"mode": mode_label(family, n, p), "cutoff_ghz": cutoff, "degeneracy": degeneracy}
               for cutoff, family, n, p, degeneracy in disc_entries(radius_m, limit_ghz)]
    # TE11 is always the lowest circular mode; take it from the list so the file is self-consistent.
    lowest = (cutoffs[0]["cutoff_ghz"] if cutoffs
              else float(jnp_zeros(1, 1)[0]) * C / (2 * math.pi * radius_m) / 1e9)
    return {
        "mode": "TE11",
        "cutoff_ghz": lowest,
        "basis": f"closed ideal-PEC circular guide, R = {radius_m!r} m; f = x c/(2 pi R), x the J_n / J'_n zeros",
        "list_limit_ghz": limit_ghz,
        "polarization_filter_limit_ghz": lowest,
        "cutoffs": cutoffs,
        "propagating_count": len(disc_entries(radius_m, frequency_ghz)),
    }


def modes_record(cross_section: Mapping[str, Any], frequency_ghz: float) -> dict | None:
    """Closed-PEC-guide modes from the declared cross-section; None (field omitted) for a polygon."""
    shape = cross_section["shape"]
    if shape == "rectangle":
        return rectangle_modes(2 * cross_section["y_e_half_m"], 2 * cross_section["z_e_half_m"], frequency_ghz)
    if shape == "disc":
        return disc_modes(cross_section["radius_m"], frequency_ghz)
    return None


def _rectangle_modes_v2(a_m: float, b_m: float) -> dict:
    """Schema 2.0 rectangle block: the lowest mode as in 1.0, the polarization-filter limit TE01 = c/(2b), and
    the regimes TE10 plus the gap-family onsets TE0n up to the 20 THz limit, ascending (ties as ``_ordered``)."""
    lowest = (1, 0) if a_m >= b_m else (0, 1)
    entries: list[tuple] = [(_rect_cutoff_ghz(a_m, b_m, 1, 0), "TE", 1, 0)]
    n = 1
    while _rect_cutoff_ghz(a_m, b_m, 0, n) <= LIST_LIMIT_GHZ:
        entries.append((_rect_cutoff_ghz(a_m, b_m, 0, n), "TE", 0, n))
        n += 1
    return {
        "mode": mode_label("TE", *lowest),
        "cutoff_ghz": _rect_cutoff_ghz(a_m, b_m, *lowest),
        "polarization_filter_limit_ghz": _rect_cutoff_ghz(a_m, b_m, 0, 1),
        "regime_cutoffs_ghz": [{"mode": mode_label(family, i, j), "cutoff_ghz": cutoff}
                               for cutoff, family, i, j in _ordered(entries)],
    }


def _disc_modes_v2(radius_m: float) -> dict:
    """Schema 2.0 disc block: TE11 is the lowest mode and the polarization-filter limit (a disc has no filter);
    the regimes are TE11, TM01 and TE21 from the closed-form zeros."""
    scale = C / (2 * math.pi * radius_m) / 1e9
    te11 = X11P * scale
    return {
        "mode": "TE11",
        "cutoff_ghz": te11,
        "polarization_filter_limit_ghz": te11,
        "regime_cutoffs_ghz": [{"mode": mode_label(family, n, p), "cutoff_ghz": x * scale}
                               for family, n, p, x in (("TE", 1, 1, X11P), ("TM", 0, 1, X01), ("TE", 2, 1, X21P))],
    }


def modes_record_v2(cross_section: Mapping[str, Any]) -> dict | None:
    """The schema-2.0 ``modes`` block from the declared cross-section, closed forms only; None (the field is
    omitted) for a polygon. Identical across all frequencies of one geometry (OR-031)."""
    shape = cross_section["shape"]
    if shape == "rectangle":
        return _rectangle_modes_v2(2 * cross_section["y_e_half_m"], 2 * cross_section["z_e_half_m"])
    if shape == "disc":
        return _disc_modes_v2(cross_section["radius_m"])
    return None


def exit_axes(x, y) -> np.ndarray:
    """Rows X, Y, Z of the exit frame in global coordinates (Z = X x Y)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    return np.vstack((x, y, np.cross(x, y)))


def _rows(matrix) -> list[list[float]]:
    return [[round(float(v), 12) + 0.0 for v in row] for row in matrix]


def _core_record(manifest: Mapping[str, Any], manifest_sha256: str, files: Mapping[str, Mapping[str, Any]],
                 schema_version: str, columns: list[str], transmittance: Mapping[str, Any]) -> dict[str, Any]:
    """Everything the two schemas share, from the run manifest and the CSV checksums; ``modes`` is the caller's.

    The caller supplies what differs: the schema version, the ``waveguide.csv`` column list and the transmittance
    block. ``exit_field.grid`` is the manifest's own object (callers that change it copy it first).
    """
    cfg = manifest["config"]
    geo = manifest["geometry"]
    exc = manifest["excitation"]
    ff = manifest["far_field"]
    ef = manifest["exit_field"]
    prj = manifest["project"]
    axes = exit_axes(geo["exit_cs"]["x"], geo["exit_cs"]["y"])
    exit_rows = _rows(axes)
    bbox = [float(v) for v in geo["bounding_box_mm"]]
    size = [hi - lo for lo, hi in zip(bbox[:3], bbox[3:])]
    frequency_ghz = float(cfg["solver"]["frequency_ghz"])
    # Arrival-direction convention: (theta, phi) are the angles of -k. Normal entry travels into the gap,
    # k = -n_entrance, so -k is the entrance outward normal: theta = acos(n_z), 180 under canonical-z.
    entrance_normal_z = float(geo["entrance_face"]["outward_normal"][2])
    normal_entry_theta_deg = round(math.degrees(math.acos(max(-1.0, min(1.0, entrance_normal_z)))), 12) + 0.0
    record: dict[str, Any] = {
        "schema_version": schema_version,
        "dataset_id": cfg["project"]["dataset_id"],
        "frequency_ghz": frequency_ghz,
        "frequency_label": manifest["solver"]["frequency_label"],
        "provenance": {
            "producer": "bbsim run_frequency", "hand_written": False,
            "git_commit": manifest.get("git_commit"), "created_utc": manifest["created_utc"],
            "job_id": manifest["job"]["id"], "manifest_sha256": manifest_sha256,
            "aedt_version": manifest["versions"].get("aedt"), "pyaedt_version": manifest["versions"].get("pyaedt"),
            "project": prj["source_path"], "project_sha256": prj["source_sha256"],
            "design": prj["base_design"], "run_design": prj["run_design"],
            "object": geo["object"], "material": geo["material"], "model_units": geo["units"],
        },
        "frames": {
            "canonical": CANONICAL,
            "pose_rule": geo["pose"],
            "hfss_global_axes_in_canonical": dict(zip("xyz", _rows(axes.T))),
            "entrance_face": cfg["geometry"]["entrance_face"],
            "exit_face": cfg["geometry"]["exit_face"],
            "entrance_outward_normal_global": [float(v) for v in geo["entrance_face"]["outward_normal"]],
            "exit_outward_normal_global": [float(v) for v in geo["exit_face"]["outward_normal"]],
            "exit_cs": {"name": geo["exit_cs"]["name"], "origin": "exit_face_center",
                        "origin_mm_global": [float(v) for v in geo["exit_face"]["center_mm"]],
                        "x": exit_rows[0], "y": exit_rows[1], "z": exit_rows[2]},
            "exit_cs_axes_in_canonical": dict(zip("xyz", _rows(axes @ axes.T))),
        },
        "excitation": {
            "coordinate_system": "global", "incidence_convention": "arrival_direction",
            "normal_entry_theta_deg": normal_entry_theta_deg,
            "plane_wave_origin": "entrance_face_center",
            "origin_mm_global": [float(v) for v in geo["entrance_face"]["center_mm"]],
            "incident_phi_deg": exc["incident_phi_deg"], "incident_theta_deg": exc["incident_theta_deg"],
            "ei_v_per_m": exc["ei_v_per_m"], "incoming_power_w": exc["incoming_power_w"],
            "polarization_convention": exc["polarization_convention"],
        },
        "far_field": {
            "coordinate_system": "exit_cs", "definition": "Theta-Phi", "component_basis": "spherical_in_exit_cs",
            "radiation_surface": "exit_face", "theta_deg": ff["theta_deg"], "phi_deg": ff["phi_deg"],
            "points_per_key": ff["points_per_angle"], "columns": list(FAR_FIELD_COLUMNS),
        },
        "exit_field": {
            "coordinate_system": "exit_cs", "points_in_si": True, "field_in_ref_cs": False,
            "field_components_frame": "hfss_global", "plane": "x_e=0", "resolution_mm": ef["resolution_mm"],
            "grid": ef["grid"], "bounds_method": geo["bounds_method"], "cross_section": geo["cross_section"],
            "outside_points": ef["outside_points"], "rim_points": "included",
            "points_per_key_retained": ef["retained_points_per_key"], "columns": list(columns),
        },
        "transmittance": dict(transmittance),
        "symmetry": {flag: flag in cfg["geometry"]["symmetry"] for flag in SYMMETRY_FLAGS},
        "boundaries": {"entrance": "radiation", "exit": "radiation", "walls": "PEC"},
        "geometry": {
            "shape": _SHAPES.get(geo["cross_section"]["shape"], "prism"),
            "extent_mm": {label: round(float(sum(abs(a) * s for a, s in zip(row, size))), 12) + 0.0
                          for label, row in zip("plg", axes)},
            "bounding_box_mm": bbox,
        },
        "files": {name: dict(entry) for name, entry in sorted(files.items())},
    }
    return record


def build_dataset_record(manifest: Mapping[str, Any], manifest_sha256: str,
                         files: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    """The schema-1.0 sidecar for one solved frequency, from the run manifest and the CSV checksums."""
    record = _core_record(manifest, manifest_sha256, files, SCHEMA_VERSION, WAVEGUIDE_COLUMNS, TRANSMITTANCE)
    modes = modes_record(manifest["geometry"]["cross_section"], float(manifest["config"]["solver"]["frequency_ghz"]))
    if modes is not None:
        record["modes"] = modes
    return record


def _require(mapping: Any, name: str, where: str) -> Any:
    if not isinstance(mapping, Mapping) or name not in mapping:
        raise ValueError(f"schema-2 sidecar: {where} lacks {name!r}")
    return mapping[name]


def _convergence_record(solver: Mapping[str, Any]) -> dict[str, Any]:
    """provenance.convergence from the manifest's ``solver.qualified`` and ``convergence_qualification``."""
    qualified = _require(solver, "qualified", "the manifest's solver block")
    if not isinstance(qualified, bool):
        raise ValueError(f"schema-2 sidecar: solver.qualified must be a bool, got {qualified!r}")
    per_pol = _require(solver, "convergence_qualification", "the manifest's solver block")
    # JSON turns the integer polarization keys of the in-memory manifest into strings; accept either.
    per_pol = {str(key): value for key, value in per_pol.items()}
    if set(per_pol) != {"0", "1"}:
        raise ValueError("schema-2 sidecar: solver.convergence_qualification must hold exactly polarizations 0 "
                         f"and 1, got {sorted(per_pol)}")
    per_polarization = {}
    for e in ("0", "1"):
        where = f"solver.convergence_qualification[{e}]"
        per_polarization[e] = {field: copy.deepcopy(_require(per_pol[e], field, where))
                               for field in _CONVERGENCE_FIELDS}
    return {"qualified": qualified, "per_polarization": per_polarization}


def _reexport_record(reexport: Mapping[str, Any]) -> dict[str, Any]:
    """provenance.reexport: exactly source_job, source_manifest_sha256, tool and git_commit; the tool name is the
    builder's, so a caller may leave it out but may not name another."""
    if not isinstance(reexport, Mapping):
        raise ValueError(f"schema-2 sidecar: reexport must be a mapping or None, got {type(reexport).__name__}")
    unexpected = sorted(set(reexport) - set(_REEXPORT_FIELDS) - {"tool"})
    missing = [field for field in _REEXPORT_FIELDS if field not in reexport]
    if unexpected or missing:
        raise ValueError(f"schema-2 sidecar: reexport must hold exactly {list(_REEXPORT_FIELDS)} (and tool "
                         f"{REEXPORT_TOOL!r}); missing {missing}, unexpected {unexpected}")
    if reexport.get("tool", REEXPORT_TOOL) != REEXPORT_TOOL:
        raise ValueError(f"schema-2 sidecar: reexport tool must be {REEXPORT_TOOL!r}, got {reexport['tool']!r}")
    return {"source_job": str(reexport["source_job"]),
            "source_manifest_sha256": reexport["source_manifest_sha256"],
            "tool": REEXPORT_TOOL, "git_commit": reexport["git_commit"]}


def build_dataset_record_v2(manifest: Mapping[str, Any], manifest_sha256: str,
                            files: Mapping[str, Mapping[str, Any]], blocks: Mapping[str, Any],
                            reexport: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """The schema-2.0 sidecar for one solved frequency (HF-047, HF-048).

    The 1.0 core (provenance, frames, excitation, far_field, exit_field, symmetry, boundaries, geometry, files)
    plus the 2.0 additions. ``blocks`` is ``dataset_v2.DatasetV2Result.blocks``: its ``transmittance``,
    ``exit_field`` (the five blocks), ``cross_term`` and ``far_field`` (``normalization``, ``ff_over_flux``)
    are copied to their field-list locations, never recomputed. ``reexport`` is None for a fresh solve, else the
    source job, its manifest hash and the git commit of the re-export. The 1.0 ``incoming_includes_cos_theta``
    and the 1.0 ``modes`` lists do not appear. The manifest and ``blocks`` are not modified.
    """
    exit_blocks = _require(blocks, "exit_field", "blocks")
    far_blocks = _require(blocks, "far_field", "blocks")
    transmittance = _require(blocks, "transmittance", "blocks")
    cross_term = _require(blocks, "cross_term", "blocks")
    chosen = {name: _require(exit_blocks, name, "blocks['exit_field']") for name in _EXIT_FIELD_BLOCKS}
    normalization = _require(far_blocks, "normalization", "blocks['far_field']")
    ff_over_flux = _require(far_blocks, "ff_over_flux", "blocks['far_field']")
    ef = manifest["exit_field"]
    step_mm = _require(ef, "step_mm", "the manifest's exit_field block")
    if step_mm is None or len(step_mm) != 3:
        raise ValueError("schema-2 sidecar: the manifest's exit_field.step_mm must be [x, y, z] in mm, "
                         f"got {step_mm!r}")
    convergence = _convergence_record(manifest["solver"])
    reexport_record = None if reexport is None else _reexport_record(reexport)

    record = _core_record(manifest, manifest_sha256, files, SCHEMA_VERSION_V2, WAVEGUIDE_COLUMNS_V2,
                          copy.deepcopy(transmittance))
    record["model"] = {"kind": MODEL_KIND_V2}
    record["conventions"] = dict(CONVENTIONS_V2)
    exit_field = record["exit_field"]
    grid = copy.deepcopy(exit_field["grid"])      # the manifest's object: copy before adding the step
    for axis, index in (("y_e", 1), ("z_e", 2)):
        grid[axis]["step"] = step_mm[index] * 1e-3   # the step the runner used, copied, never recomputed
    exit_field["grid"] = grid
    exit_field.update({name: copy.deepcopy(value) for name, value in chosen.items()})
    exit_field["plane_tolerance_m"] = PLANE_TOLERANCE_M_V2
    record["cross_term"] = copy.deepcopy(cross_term)
    record["far_field"]["normalization"] = normalization
    record["far_field"]["ff_over_flux"] = copy.deepcopy(ff_over_flux)
    record["provenance"]["convergence"] = convergence
    if reexport_record is not None:
        record["provenance"]["reexport"] = reexport_record
    modes = modes_record_v2(manifest["geometry"]["cross_section"])
    if modes is not None:
        limit = modes["polarization_filter_limit_ghz"]
        valid_below = record["cross_term"].get("valid_below_ghz")
        if (not isinstance(valid_below, (int, float)) or isinstance(valid_below, bool)
                or not math.isclose(valid_below, limit, rel_tol=_CROSS_TERM_REL_TOL)):
            raise ValueError(f"schema-2 sidecar: cross_term.valid_below_ghz {valid_below!r} must equal "
                             f"modes.polarization_filter_limit_ghz {limit!r} (BBRsim compares them to 1e-6)")
        record["modes"] = modes
    return record
