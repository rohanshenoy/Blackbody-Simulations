"""The HFSS dataset sidecar, schema 1.0, agreed with the BBRSimulation session on 2026-10-04.

Built only from the run manifest as written to disk, so the sidecar cannot disagree with it. Pure:
no AEDT; scipy is imported only for the circular-guide Bessel zeros. Canonical labels p, l, g
(propagation, long, gap) are the Geant4 crack-local x, y, z and equal the exit-CS axes X, Y, Z.
"""
from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numpy as np

from bbsim.sampling import SPEED_OF_LIGHT_M_PER_S as C
from bbsim.schema import FAR_FIELD_COLUMNS, WAVEGUIDE_COLUMNS

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


def exit_axes(x, y) -> np.ndarray:
    """Rows X, Y, Z of the exit frame in global coordinates (Z = X x Y)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    return np.vstack((x, y, np.cross(x, y)))


def _rows(matrix) -> list[list[float]]:
    return [[round(float(v), 12) + 0.0 for v in row] for row in matrix]


def build_dataset_record(manifest: Mapping[str, Any], manifest_sha256: str,
                         files: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    """The schema-1.0 sidecar for one solved frequency, from the run manifest and the CSV checksums."""
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
        "schema_version": SCHEMA_VERSION,
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
            "points_per_key_retained": ef["retained_points_per_key"], "columns": list(WAVEGUIDE_COLUMNS),
        },
        "transmittance": dict(TRANSMITTANCE),
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
    modes = modes_record(geo["cross_section"], frequency_ghz)
    if modes is not None:
        record["modes"] = modes
    return record
