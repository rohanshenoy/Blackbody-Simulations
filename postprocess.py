"""Per-key photon sampling from an HFSS dataset pair, following BBRsim's sampler (CDX-004, HF-044(2)).

Schema 1.0 (default): coherent polarization cross term, load-time normalization, sin(Theta)-weighted outgoing
directions. Schema 2.x (read with its dataset sidecar): T = OutgoingPower / (IngoingPower x InterceptFactor), a cap
on each photon's own mixture instead of the rescale, exit positions and the cross term from the quadrature-weighted
tangential field (Ex, Ey, read by position), and a dark basis (a numerical-zero flux) contributing nothing to any
mixture. Key blocks are paired by key and by position within the key, never by DataFrame label or block order."""
import json
import math
import os
import re

import pandas as pd
import numpy as np

# BBRsim's sampler (BBRHFSSData.cc, mirrored by BBRSimulation/tools/python/bbrsim/hfss.py) is the one supported
# physics; this script follows it without importing it, and without importing bbsim (CDX-004, HF-044(2)): it is a
# third implementation, so a drift in either of the others shows up as a failing conformance test.
K_MIN_NORMAL_COMPONENT = 1e-6   # kMinNormalComponent in BBRHFSSData.cc
EXIT_VALUES = ["Ex_real", "Ey_real", "Ez_real", "Ex_imag", "Ey_imag", "Ez_imag"]   # 1.0: read by label
# 2.0 header: the stored slots are interleaved (Re Ex, Im Ex, Re Ey, Im Ey, Re Ez, Im Ez, columns 9-14) and are
# read by position, Ex = (9, 10) and Ey = (11, 12), never by label; InterceptFactor is the last column.
V2_COLUMNS = ["Freq", "Ephi", "IWavePhi", "IWaveTheta", "OutgoingPower", "IngoingPower", "X", "Y", "Z",
              "Ex_real", "Ex_imag", "Ey_real", "Ey_imag", "Ez_real", "Ez_imag", "InterceptFactor"]
V2_SLOTS = slice(9, 15)
KEY_TOL_DEG = 1e-9      # an incidence key matches a row or a sidecar record within this, in degrees
COORD_TOL_M = 1e-12     # paired exit points agree within this (the producer's Ephi equality, BBRsim's BBR012)
ANGLE_TOL_DEG = 1e-9    # paired far-field angles agree within this
NODE_TOL = 1e-6         # an exit point lies within this fraction of the lattice step of its node
FACTOR_TOL = 1e-6       # InterceptFactor equals the sidecar's c_k within this, relative (BBRsim's tolerance)
_EPHI_DIR = re.compile(r"^(?P<stem>.+)_Ephi=(?P<ephi>[01])$")


def _rows(df, phi, theta):
    """The rows of one incidence key (within KEY_TOL_DEG), in frame order, with a fresh 0..n-1 index."""
    near = ((np.abs(df["IWavePhi"].to_numpy(dtype=float) - phi) <= KEY_TOL_DEG)
            & (np.abs(df["IWaveTheta"].to_numpy(dtype=float) - theta) <= KEY_TOL_DEG))
    return df[near].reset_index(drop=True)


def _key(df, phi, theta):
    rows = _rows(df, phi, theta)
    if rows.empty:
        raise ValueError(f"No data for IWavePhi={phi}, IWaveTheta={theta}")
    return rows


def _pair(df_0, df_1, phi, theta, coords, tol, what):
    """The key's rows of both polarizations, paired by position within the key as BBRsim pairs them (CDX-032): the
    same count and every paired coordinate within ``tol``, else ValueError."""
    a, b = _key(df_0, phi, theta), _key(df_1, phi, theta)
    if len(a) != len(b) or not np.all(np.abs(a[coords].to_numpy(dtype=float) - b[coords].to_numpy(dtype=float)) <= tol):
        raise ValueError(f"Mismatch in {what} between df_0 and df_1.")
    return a, b


def _major(version, what):
    """The major (1 or 2) of a schema version text such as "1.0", "2.0" or "2": the text before the first ".".
    Anything else (not a string, no numeric major, a major other than 1 or 2) is a ValueError."""
    major = version.split(".")[0] if isinstance(version, str) else None
    if major not in ("1", "2"):
        raise ValueError(f"{what} {version!r} is not a 1.x or 2.x schema version")
    return int(major)


def _is_v2(schema, sidecar):
    """True for a 2.x table, which needs a 2.x sidecar; False for 1.x, whose sidecar (if any) must be 1.x too."""
    major = _major(schema, "schema")
    if sidecar is None:
        if major == 2:
            raise ValueError(f"schema {schema!r} is read with its 2.x dataset sidecar; none was given")
        return False
    declared = _major(sidecar.get("schema_version") if isinstance(sidecar, dict) else None,
                      "the sidecar's schema_version")
    if declared != major:
        raise ValueError(f"schema {schema!r} contradicts the sidecar's schema_version {sidecar['schema_version']!r}")
    return major == 2


def _v1_frame(df):
    if "InterceptFactor" in df.columns:
        raise ValueError("a 2.0 table (it has InterceptFactor) needs schema='2.0' and its dataset sidecar")


def _v2_frame(df):
    if list(df.columns) != V2_COLUMNS:
        raise ValueError(f"a 2.0 waveguide table has the columns {V2_COLUMNS}, got {list(df.columns)}")


def _num(value, where):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{where} must be a finite number, got {value!r}")
    return float(value)


def _sidecar_keys(sidecar):
    """The 2.x sidecar's per-key records, checked: [{key, c_k, t_flux (Ephi 0, 1), dark (Ephi 0, 1)}] from
    transmittance.per_key and exit_field.scale.per_key, which must hold the same keys once each."""
    try:
        per_key = sidecar["transmittance"]["per_key"]
        scale = sidecar["exit_field"]["scale"]["per_key"]
    except (KeyError, TypeError) as exc:
        raise ValueError(f"the 2.x sidecar lacks transmittance.per_key or exit_field.scale.per_key: {exc!r}") from exc
    if not isinstance(per_key, list) or not isinstance(scale, list) or not per_key:
        raise ValueError("transmittance.per_key and exit_field.scale.per_key must be non-empty lists")
    records = []
    for i, rec in enumerate(per_key):
        at = f"transmittance.per_key[{i}]"
        try:
            records.append({"key": (_num(rec["phi_deg"], f"{at}.phi_deg"), _num(rec["theta_deg"], f"{at}.theta_deg")),
                            "c_k": _num(rec["c_k"], f"{at}.c_k"),
                            "t_flux": [_num(rec["polarizations"][str(e)]["t_flux"], f"{at}.polarizations.{e}.t_flux")
                                       for e in (0, 1)],
                            "dark": None})
        except (KeyError, TypeError) as exc:
            raise ValueError(f"{at} is malformed: {exc!r}") from exc
    for i, rec in enumerate(scale):
        at = f"exit_field.scale.per_key[{i}]"
        try:
            key = (_num(rec["phi_deg"], f"{at}.phi_deg"), _num(rec["theta_deg"], f"{at}.theta_deg"))
            dark = [rec["dark0"], rec["dark1"]]
        except (KeyError, TypeError) as exc:
            raise ValueError(f"{at} is malformed: {exc!r}") from exc
        if not all(isinstance(d, bool) for d in dark):
            raise ValueError(f"{at}.dark0 and .dark1 must be true or false, got {dark!r}")
        match = [r for r in records
                 if abs(r["key"][0] - key[0]) <= KEY_TOL_DEG and abs(r["key"][1] - key[1]) <= KEY_TOL_DEG]
        if len(match) != 1 or match[0]["dark"] is not None:
            raise ValueError(f"{at}: (phi, theta) = ({key[0]:g}, {key[1]:g}) deg has no unique transmittance.per_key "
                             "record, or a second scale record")
        match[0]["dark"] = dark
    if any(r["dark"] is None for r in records):
        raise ValueError("exit_field.scale.per_key lacks a record for a key of transmittance.per_key")
    return records


def _check_dark(rec, rows):
    """OR-038: a dark record needs t_flux exactly 0 in the sidecar and, when the key's two exit tables are given,
    OutgoingPower exactly 0 on every row of the key in its polarization's table. Nothing is repaired."""
    for e in (0, 1):
        if not rec["dark"][e]:
            continue
        at = f"(phi, theta) = ({rec['key'][0]:g}, {rec['key'][1]:g}) deg, Ephi={e}"
        if rec["t_flux"][e] != 0.0:
            raise ValueError(f"{at}: dark{e} is true but the sidecar's t_flux is {rec['t_flux'][e]!r}, not exactly 0")
        if rows is not None:
            power = rows[e]["OutgoingPower"].to_numpy(dtype=float)
            bad = np.flatnonzero(power != 0.0)
            if bad.size:
                raise ValueError(f"{at}: dark{e} is true but OutgoingPower is {float(power[bad[0]])!r} on row "
                                 f"{int(bad[0]) + 1} of the key's {power.size} ({bad.size} nonzero), not exactly 0")


def _v2_key(sidecar, phi, theta, rows=None):
    """The sidecar's record of one key, its dark flags validated against the key's exit rows (when given) first."""
    for rec in _sidecar_keys(sidecar):
        if abs(rec["key"][0] - phi) <= KEY_TOL_DEG and abs(rec["key"][1] - theta) <= KEY_TOL_DEG:
            _check_dark(rec, rows)
            return rec
    raise ValueError(f"The sidecar has no record of IWavePhi={phi}, IWaveTheta={theta}")


def _quadrature_weights(y_m, z_m, sidecar):
    """Dimensionless quadrature weight of each exit point: a disc 1 everywhere; a rectangle the trapezoid rule on the
    declared lattice, node = round((coordinate - min) / step) on Y and on Z, 1/2 at node 0 or count - 1 per axis,
    multiplied (so 1/4 at a corner). The points must be the whole lattice, each on its own node."""
    shape = sidecar.get("exit_field", {}).get("cross_section", {}).get("shape")
    y, z = np.asarray(y_m, dtype=float), np.asarray(z_m, dtype=float)
    if shape == "disc":
        return np.ones(y.size)
    if shape != "rectangle":
        raise ValueError(f"exit_field.cross_section.shape must be 'rectangle' or 'disc', got {shape!r}")
    if not (np.all(np.isfinite(y)) and np.all(np.isfinite(z))):
        raise ValueError("exit-point Y and Z must be finite")
    weights, nodes, counts = np.ones(y.size), [], []
    for axis, coord in (("y_e", y), ("z_e", z)):
        try:
            grid = sidecar["exit_field"]["grid"][axis]
            lo, step, count = _num(grid["min"], f"{axis}.min"), _num(grid["step"], f"{axis}.step"), grid["count"]
        except (KeyError, TypeError) as exc:
            raise ValueError(f"exit_field.grid.{axis} is malformed: {exc!r}") from exc
        if isinstance(count, bool) or not isinstance(count, int) or count < 2:
            raise ValueError(f"exit_field.grid.{axis}.count must be an integer of at least 2, got {count!r}")
        if not step > 0.0:
            raise ValueError(f"exit_field.grid.{axis}.step must be positive, got {step!r}")
        offset = (coord - lo) / step
        node = np.rint(offset)
        if np.any((node < 0) | (node > count - 1)):
            i = int(np.flatnonzero((node < 0) | (node > count - 1))[0])
            raise ValueError(f"exit point {i} lies outside the {count}-node {axis} lattice (node {node[i]:g})")
        if np.any(np.abs(offset - node) > NODE_TOL):
            i = int(np.argmax(np.abs(offset - node)))
            raise ValueError(f"exit point {i} is off its {axis} node by {abs(offset[i] - node[i]):g} of the step")
        weights = weights * np.where((node == 0) | (node == count - 1), 0.5, 1.0)
        nodes.append(node.astype(int))
        counts.append(count)
    if y.size != counts[0] * counts[1]:
        raise ValueError(f"a {counts[0]} x {counts[1]} lattice has {counts[0] * counts[1]} points, got {y.size}")
    if np.bincount(nodes[0] * counts[1] + nodes[1], minlength=y.size).max() > 1:
        raise ValueError("two exit points sit at the same node: a lattice holds each node once")
    return weights


def _e_t(values):
    """Complex (Ex, Ey) of the six stored slots, by position."""
    return values[:, 0] + 1j * values[:, 1], values[:, 2] + 1j * values[:, 3]


def _mix(c0, x0, c1, x1, dark):
    """c0 x0 + c1 x1, a dark basis contributing nothing (its stored values are not read, NaN included)."""
    zero = np.zeros_like(x0)
    return (zero if dark[0] else c0 * x0) + (zero if dark[1] else c1 * x1)


def key_transmittances(phi_in, theta_in, df_0, df_1, schema="1.0", sidecar=None):
    """(T0, T1, Re rho) of one incidence key.

    1.0, as BBRsim's loader computes them: T_e = OutgoingPower / IngoingPower; rho = sum E0.E1* /
    sqrt(sum|E0|^2 sum|E1|^2) over the exit points, rows paired by order (its real part is the sum of products of
    the six stored values, whatever their column order); when the key's largest transmittance over linear
    polarizations, the top eigenvalue of [[T0, c], [c, T1]] with c = sqrt(T0 T1) Re rho, exceeds 1, both T are
    divided by it.

    2.x (``schema="2.0"`` and the sidecar): T_e = OutgoingPower / (IngoingPower x InterceptFactor) from the key's
    first row, uncapped and not rescaled (the cap belongs to each photon's own mixture, ``get_S21``); Re rho on the
    quadrature-weighted tangential field (Ex, Ey by position), 0 if either basis is dark. A dark basis is validated
    first: t_flux and every OutgoingPower of the key exactly 0.
    """
    if _is_v2(schema, sidecar):
        return _transmittances_v2(phi_in, theta_in, df_0, df_1, sidecar)
    _v1_frame(df_0)
    _v1_frame(df_1)
    a, b = _pair(df_0, df_1, phi_in, theta_in, ["X", "Y", "Z"], COORD_TOL_M, "exit points")
    t = []
    for rows in (a, b):
        p_in = float(rows["IngoingPower"].iloc[0])
        if p_in == 0:
            raise ZeroDivisionError("IngoingPower is zero.")
        t.append(float(rows["OutgoingPower"].iloc[0]) / p_in)
    v0, v1 = a[EXIT_VALUES].to_numpy(float), b[EXIT_VALUES].to_numpy(float)
    p0, p1 = float(np.sum(v0 * v0)), float(np.sum(v1 * v1))
    rho = float(np.sum(v0 * v1) / np.sqrt(p0 * p1)) if p0 > 0 and p1 > 0 else 0.0
    c = np.sqrt(t[0] * t[1]) * rho
    lam = 0.5 * (t[0] + t[1]) + float(np.hypot(0.5 * (t[0] - t[1]), c))
    if lam > 1.0:
        t = [x / lam for x in t]
    return t[0], t[1], rho


def _transmittances_v2(phi_in, theta_in, df_0, df_1, sidecar):
    _v2_frame(df_0)
    _v2_frame(df_1)
    a, b = _pair(df_0, df_1, phi_in, theta_in, ["X", "Y", "Z"], COORD_TOL_M, "exit points")
    rec = _v2_key(sidecar, phi_in, theta_in, (a, b))
    t = []
    for e, rows in enumerate((a, b)):
        p_in, out = float(rows["IngoingPower"].iloc[0]), float(rows["OutgoingPower"].iloc[0])
        factor = rows["InterceptFactor"].to_numpy(dtype=float)
        if not (math.isfinite(p_in) and p_in > 0.0 and math.isfinite(out) and np.all(np.isfinite(factor))
                and np.all(factor > 0.0)):
            raise ValueError(f"Ephi={e}: IngoingPower, OutgoingPower and InterceptFactor must be finite and positive")
        if not np.all(np.abs(factor - rec["c_k"]) <= FACTOR_TOL * rec["c_k"]):
            raise ValueError(f"Ephi={e}: InterceptFactor differs from the sidecar's c_k {rec['c_k']!r} on some row of "
                             f"the key (limit {FACTOR_TOL:g} relative)")
        t.append(out / (p_in * float(factor[0])))
        if t[-1] < 0.0:
            raise ValueError(f"Ephi={e}: negative transmittance {t[-1]!r}")
    if rec["dark"][0] or rec["dark"][1]:
        return t[0], t[1], 0.0
    weights = _quadrature_weights(a["Y"].to_numpy(float), a["Z"].to_numpy(float), sidecar)
    e0, e1 = _e_t(a.iloc[:, V2_SLOTS].to_numpy(float)), _e_t(b.iloc[:, V2_SLOTS].to_numpy(float))

    def power(e):
        return float(sum(np.sum(weights * np.abs(c) ** 2) for c in e))

    p0, p1 = power(e0), power(e1)
    cross = float(np.real(sum(np.sum(weights * x * np.conj(y)) for x, y in zip(e0, e1))))
    if not (math.isfinite(p0) and math.isfinite(p1) and p0 > 0.0 and p1 > 0.0 and math.isfinite(cross)):
        raise ValueError(f"rho is undefined for a non-dark pair: weighted tangential powers {p0!r} and {p1!r}")
    return t[0], t[1], cross / (math.sqrt(p0) * math.sqrt(p1))


def _stem(directory, ephi):
    """The ``<ID>_<label>`` stem of a directory named ``<stem>_Ephi=<ephi>``; any other name is a ValueError."""
    match = _EPHI_DIR.match(os.path.basename(os.path.normpath(directory)))
    if match is None or int(match["ephi"]) != ephi:
        raise ValueError(f"{directory} must be a directory named <ID>_<label>_Ephi={ephi}")
    return match["stem"]


def _sidecar_of(directory, stem):
    """The sidecar ``<stem>.dataset.json`` beside an Ephi directory, or None when there is none."""
    path = os.path.join(os.path.dirname(os.path.normpath(directory)), f"{stem}.dataset.json")
    return path if os.path.isfile(path) else None


def load_datasets(ephi0_dir, ephi1_dir):
    """Load waveguide and far-field tables for the Ephi=0 and Ephi=1 datasets of one frequency, and the dataset
    sidecar beside them (``None`` for a 1.0 dataset without one).

    The directories must be named ``<stem>_Ephi=0`` and ``<stem>_Ephi=1`` with the same ``<ID>_<label>`` stem, and
    name the same sidecar (ValueError otherwise). A 2.x sidecar needs two 16-column 2.0 tables, a 16-column table
    needs a 2.x sidecar, and an unknown ``schema_version`` is refused. Every dark record of a 2.x sidecar is
    validated against both tables here, before any helper masks anything: t_flux and every OutgoingPower of a dark
    key exactly 0.
    """
    stems = (_stem(ephi0_dir, 0), _stem(ephi1_dir, 1))
    if stems[0] != stems[1]:
        raise ValueError(f"the directories belong to different datasets: stem {stems[0]!r} against {stems[1]!r}")
    paths = [_sidecar_of(d, stems[0]) for d in (ephi0_dir, ephi1_dir)]
    if paths[0] != paths[1]:
        raise ValueError(f"the two directories name different sidecars: {paths[0]} and {paths[1]}")
    wg0, wg1 = (pd.read_csv(os.path.join(d, "waveguide.csv")) for d in (ephi0_dir, ephi1_dir))
    ff0, ff1 = (pd.read_csv(os.path.join(d, "far_field.csv")) for d in (ephi0_dir, ephi1_dir))
    sidecar = None
    if paths[0] is not None:
        try:
            with open(paths[0]) as handle:
                sidecar = json.load(handle)
        except (OSError, ValueError) as exc:
            raise ValueError(f"the sidecar {paths[0]} is unreadable: {exc}") from exc
        if not isinstance(sidecar, dict):
            raise ValueError(f"the sidecar {paths[0]} is not a JSON object")
    v2 = sidecar is not None and _major(sidecar.get("schema_version"), "the sidecar's schema_version") == 2
    for df in (wg0, wg1):
        if v2:
            _v2_frame(df)
        elif "InterceptFactor" in df.columns:
            raise ValueError("a 2.0 waveguide table (it has InterceptFactor) needs its 2.x dataset sidecar beside "
                             "the Ephi directories")
    if v2:
        for rec in _sidecar_keys(sidecar):
            rows = tuple(_rows(df, *rec["key"]) for df in (wg0, wg1))
            _check_dark(rec, rows)
    return wg0, wg1, ff0, ff1, sidecar


#%% Stage 1: Finding the transmission power fraction of the waveguide

def get_S21(phi_in, theta_in, polarization_0, polarization_1, df_0, df_1, schema="1.0", sidecar=None):
    """
    Retrieve the power transmission fraction |S21|^2 from HFSS data,
    weighted by the incident wave's polarization.

    Parameters:
        phi_in (float): IWavePhi angle (in degrees)
        theta_in (float): IWaveTheta angle (in degrees)
        polarization_0 (float): Amplitude of polarization mode 0 (Ephi=0)
        polarization_1 (float): Amplitude of polarization mode 1 (Ephi=1)
        df_0 (pd.DataFrame): DataFrame for polarization mode 0 (Ephi=0)
        df_1 (pd.DataFrame): DataFrame for polarization mode 1 (Ephi=1)
        schema (str): "1.0" (default) or "2.x", the latter with its dataset ``sidecar``

    Returns:
        float: BBRsim's transmission probability for this linear polarization, in [0, 1]. The mixture
        T_mix = a^2 T0 + b^2 T1 + 2 a b sqrt(T0 T1) Re rho for the normalized coefficients (a, b) is capped on
        this photon's own value: min(1, max(0, T_mix)); under 1.0 (T0, T1) are already rescaled by the key's top
        eigenvalue, under 2.x they are not (``key_transmittances``).
    """

    norm = float(np.hypot(polarization_0, polarization_1))
    if norm == 0:
        raise ValueError("polarization amplitudes are both zero")
    e_t, e_p = polarization_0 / norm, polarization_1 / norm
    t0, t1, rho = key_transmittances(phi_in, theta_in, df_0, df_1, schema, sidecar)
    t = e_t ** 2 * t0 + e_p ** 2 * t1 + 2.0 * e_t * e_p * np.sqrt(t0 * t1) * rho
    return min(1.0, max(0.0, float(t)))
#%% Stage 2: Determining a probability distribution over the location of emission from the outgoing face.

def get_face_emission_distribution(phi_in, theta_in, polarization_0, polarization_1,
                              df_0, df_1, schema="1.0", sidecar=None):
    """
    Compute normalized emission probability distribution over the output face
    for a given incoming angle and polarization, using linear superposition of E-fields.
    Matches BBRsim's exit-position weights: 1.0 samples |E|^2 over the six stored values with uniform weights; 2.x
    samples q_i = w_i |a E_t0,i + b E_t1,i|^2 with E_t = (Ex, Ey) read by position, the dimensionless quadrature
    weights w_i of the sidecar's lattice, and a dark basis (numerical-zero flux) contributing nothing.

    Parameters:
        phi_in (float): IWavePhi in degrees
        theta_in (float): IWaveTheta in degrees
        polarization_0 (float): Amplitude of polarization mode 0 (Ephi=0)
        polarization_1 (float): Amplitude of polarization mode 1 (Ephi=1)
        df_0 (pd.DataFrame): E-field data for (Etheta=1, Ephi=0)
        df_1 (pd.DataFrame): E-field data for (Etheta=0, Ephi=1)
        schema (str): "1.0" (default) or "2.x", the latter with its dataset ``sidecar``

    Returns:
        pd.DataFrame: With columns [X, Y, Z, E2, Probability]; E2 is the mixed field's |E|^2 (|E_t|^2 under 2.x),
        Probability the normalized weight. A zero total is a ValueError: the caller does not sample at zero
        transmission, and a positive transmission with nothing to sample is an error.
    """
    if _is_v2(schema, sidecar):
        return _face_distribution_v2(phi_in, theta_in, polarization_0, polarization_1, df_0, df_1, sidecar)
    _v1_frame(df_0)
    _v1_frame(df_1)
    f0, f1 = _pair(df_0, df_1, phi_in, theta_in, ["X", "Y", "Z"], COORD_TOL_M, "spatial coordinates")

    # Reconstruct complex E-fields
    E0 = {
        "Ex": f0["Ex_real"] + 1j * f0["Ex_imag"],
        "Ey": f0["Ey_real"] + 1j * f0["Ey_imag"],
        "Ez": f0["Ez_real"] + 1j * f0["Ez_imag"],
    }

    E1 = {
        "Ex": f1["Ex_real"] + 1j * f1["Ex_imag"],
        "Ey": f1["Ey_real"] + 1j * f1["Ey_imag"],
        "Ez": f1["Ez_real"] + 1j * f1["Ez_imag"],
    }

    # Linearly combine fields according to polarization
    E_total = {
        comp: polarization_0 * E0[comp] + polarization_1 * E1[comp]
        for comp in ["Ex", "Ey", "Ez"]
    }

    # Compute |E|^2
    E2 = np.abs(E_total["Ex"])**2 + np.abs(E_total["Ey"])**2 + np.abs(E_total["Ez"])**2

    # Normalize to get probability distribution
    total_E2 = E2.sum()
    if total_E2 == 0:
        raise ValueError("Total field energy is zero — can't normalize.")

    probabilities = E2 / total_E2

    # Return results
    result = f0[["X", "Y", "Z"]].copy()
    result["E2"] = E2
    result["Probability"] = probabilities

    return result


def _face_distribution_v2(phi_in, theta_in, polarization_0, polarization_1, df_0, df_1, sidecar):
    """2.x exit positions: q_i = w_i |a E_t0,i + b E_t1,i|^2 normalized, E_t = (Ex, Ey) by position, w_i the
    dimensionless quadrature weights of the sidecar's lattice, a dark basis contributing nothing (E2 is |E_t|^2)."""
    _v2_frame(df_0)
    _v2_frame(df_1)
    f0, f1 = _pair(df_0, df_1, phi_in, theta_in, ["X", "Y", "Z"], COORD_TOL_M, "spatial coordinates")
    rec = _v2_key(sidecar, phi_in, theta_in, (f0, f1))
    ex0, ey0 = _e_t(f0.iloc[:, V2_SLOTS].to_numpy(float))
    ex1, ey1 = _e_t(f1.iloc[:, V2_SLOTS].to_numpy(float))
    E2 = (np.abs(_mix(polarization_0, ex0, polarization_1, ex1, rec["dark"])) ** 2
          + np.abs(_mix(polarization_0, ey0, polarization_1, ey1, rec["dark"])) ** 2)
    weights = _quadrature_weights(f0["Y"].to_numpy(float), f0["Z"].to_numpy(float), sidecar) * E2
    if not weights.sum() > 0:
        raise ValueError("Total field energy is zero (or not finite) — can't normalize.")
    result = f0[["X", "Y", "Z"]].copy()
    result["E2"] = E2
    result["Probability"] = weights / weights.sum()
    return result

#%% Stage 3: Determining a probability distribution over the angle of emission parametrized by (theta, phi)

def get_angular_emission_distribution(phi_in, theta_in, polarization_0, polarization_1,
                                      df_0, df_1, schema="1.0", sidecar=None):
    """
    Compute normalized emission probability distribution over outgoing angles (theta, phi)
    for given incoming angle and polarization, using linear superposition of E-fields
    in spherical coordinates (rEtheta, rEphi).

    Each row stands for sin(Theta) dTheta dPhi of solid angle (the constant steps cancel on normalization). 1.0
    keeps the rows with sin(Theta) cos(Phi) >= 1e-6 and 2.x keeps the rows above 1e-6 (as the field list says); the
    others (the polar axis, the exit plane, directions behind it) carry none.
    Under 2.x the far field is the unscaled CSV rE and a dark basis contributes nothing; exit-area weights and the
    exit E_t are not direction weights. Called directly, it checks only the sidecar half of the dark rule (t_flux
    exactly 0 for a dark basis), because far_field.csv carries no OutgoingPower; the waveguide helpers
    (``key_transmittances``, ``get_face_emission_distribution``) and ``load_datasets`` check the table half.

    Parameters:
        phi_in (float): Incoming IWavePhi in degrees
        theta_in (float): Incoming IWaveTheta in degrees
        polarization_0 (float): Amplitude of polarization mode 0 (Ephi=0)
        polarization_1 (float): Amplitude of polarization mode 1 (Ephi=1)
        df_0 (pd.DataFrame): Far-field E data for (Etheta=1, Ephi=0)
        df_1 (pd.DataFrame): Far-field E data for (Etheta=0, Ephi=1)
        schema (str): "1.0" (default) or "2.x", the latter with its dataset ``sidecar``

    Returns:
        pd.DataFrame: With columns [Theta, Phi, E2, Probability, Polarization]
                     Polarization column contains normalized complex vectors [Etheta, Ephi].
        A zero total is a ValueError (the caller does not sample at zero transmission).
    """
    v2 = _is_v2(schema, sidecar)
    f0, f1 = _pair(df_0, df_1, phi_in, theta_in, ["Theta", "Phi"], ANGLE_TOL_DEG, "angular coordinates")

    # Reconstruct spherical components
    Etheta_0 = (f0["rEtheta_real"] + 1j * f0["rEtheta_imag"]).to_numpy()
    Ephi_0   = (f0["rEphi_real"]   + 1j * f0["rEphi_imag"]).to_numpy()

    Etheta_1 = (f1["rEtheta_real"] + 1j * f1["rEtheta_imag"]).to_numpy()
    Ephi_1   = (f1["rEphi_real"]   + 1j * f1["rEphi_imag"]).to_numpy()

    # Linear superposition in spherical basis
    dark = _v2_key(sidecar, phi_in, theta_in)["dark"] if v2 else (False, False)
    Etheta_total = _mix(polarization_0, Etheta_0, polarization_1, Etheta_1, dark)
    Ephi_total = _mix(polarization_0, Ephi_0, polarization_1, Ephi_1, dark)

    # |E|^2 = |Etheta|^2 + |Ephi|^2 is power per unit solid angle; on the (Theta, Phi) grid each row stands for
    # sin(Theta) dTheta dPhi of solid angle, and rows in or behind the exit plane carry none (BBRsim's sampler).
    E2 = np.abs(Etheta_total)**2 + np.abs(Ephi_total)**2
    sin_t = np.sin(np.radians(f0["Theta"].to_numpy(dtype=float)))
    outward = sin_t * np.cos(np.radians(f0["Phi"].to_numpy(dtype=float)))
    if v2:
        if not (np.all(np.isfinite(outward)) and np.all(np.isfinite(E2))):
            raise ValueError("the far-field angles and the non-dark rE values must be finite")
        weights = np.where(outward > K_MIN_NORMAL_COMPONENT, sin_t * E2, 0.0)
    else:
        weights = np.where(outward < K_MIN_NORMAL_COMPONENT, 0.0, sin_t * E2)

    total = weights.sum()
    if (not total > 0) if v2 else total == 0:
        raise ValueError("No outgoing far-field power at this key — can't normalize.")

    probabilities = weights / total

    # Normalize outgoing polarization vectors at each point
    polarization_vectors = np.vstack([Etheta_total, Ephi_total]).T  # shape (N, 2)
    norms = np.linalg.norm(polarization_vectors, axis=1)
    # Avoid division by zero (set zero norm to 1 to avoid NaNs, those rows will have zero probability anyway)
    norms_safe = np.where(norms == 0, 1, norms)
    normalized_polarizations = polarization_vectors / norms_safe[:, None]

    # Construct result DataFrame
    result = f0[["Theta", "Phi"]].copy()
    result["E2"] = E2
    result["Probability"] = probabilities
    result["Polarization"] = [vec for vec in normalized_polarizations]

    return result
#%%
# Sample run
def sample_emission_point(face_distribution):
    """Randomly sample a point from the outgoing face distribution."""
    probabilities = face_distribution["Probability"].values
    indices = np.arange(len(probabilities))
    sampled_index = np.random.choice(indices, p=probabilities)
    sampled_point = face_distribution.iloc[sampled_index][["X", "Y", "Z"]]
    return sampled_point, sampled_index

def sample_emission_angle(angular_distribution):
    """Randomly sample an emission angle (Theta, Phi) from the angular distribution."""
    probabilities = angular_distribution["Probability"].values
    indices = np.arange(len(probabilities))
    sampled_index = np.random.choice(indices, p=probabilities)
    sampled_angle = angular_distribution.iloc[sampled_index][["Theta", "Phi"]]
    return sampled_angle, sampled_index

def simulate_photon_emission(
    phi_in,
    theta_in,
    polarization_0,
    polarization_1,
    waveguide_0_df,
    waveguide_1_df,
    far_field_0_df,
    far_field_1_df,
    schema="1.0",
    sidecar=None
):
    """
    Simulate a single photon emission event through a waveguide and into free space.

    The simulation proceeds in three stages:
    1. Compute the transmission probability |S21|² for the given incident angles and polarizations,
       and decide probabilistically if the photon reaches the waveguide's output face.
    2. If the photon reaches the end face, sample a spatial emission point from the normalized
       outgoing face distribution.
    3. Sample an angular emission direction and polarization from the normalized far-field distribution.

    Parameters
    ----------
    phi_in : float
        Incident azimuthal angle (degrees).
    theta_in : float
        Incident polar angle (degrees).
    polarization_0 : str
        Polarization state label for waveguide 0 input.
    polarization_1 : str
        Polarization state label for waveguide 1 input.
    waveguide_0_df : pandas.DataFrame
        Data for waveguide 0 containing S-parameters and spatial emission info.
    waveguide_1_df : pandas.DataFrame
        Data for waveguide 1 containing S-parameters and spatial emission info.
    far_field_0_df : pandas.DataFrame
        Far-field angular distribution data for waveguide 0.
    far_field_1_df : pandas.DataFrame
        Far-field angular distribution data for waveguide 1.
    schema : str
        "1.0" (default) or "2.x"; a 2.x table is read with its dataset ``sidecar``.

    Returns
    -------
    dict
        A dictionary with the following keys:
        - "S21_squared" : float
            Transmission probability |S21|².
        - "reached_end" : bool
            Whether the photon reached the output face.
        - "sampled_point" : pandas.Series, optional
            Coordinates (X, Y, Z) of the sampled emission point (only if reached_end is True).
        - "point_probability" : float, optional
            Probability density at the sampled point.
        - "face_distribution" : pandas.DataFrame, optional
            Full normalized spatial emission distribution.
        - "sampled_angle" : tuple, optional
            Emission angles (Theta, Phi) of the photon.
        - "angle_probability" : float, optional
            Probability density at the sampled angle.
        - "angular_distribution" : pandas.DataFrame, optional
            Full normalized angular emission distribution.
        - "polarization_at_sampled_angle" : str, optional
            Polarization state corresponding to the sampled angle.
    """
    # Stage 1: Compute |S21|^2 (power transmission)
    S21_squared = get_S21(
        phi_in, theta_in, polarization_0, polarization_1, waveguide_0_df, waveguide_1_df, schema, sidecar
    )

    # Decide if the photon reaches the end face; at zero transmission nothing is sampled
    reached_end = np.random.rand() < S21_squared
    if not reached_end or S21_squared <= 0.0:
        # Photon is lost: nothing else to simulate
        return {
            "S21_squared": S21_squared,
            "reached_end": False
        }

    # Stage 2: Get normalized spatial emission distribution on outgoing face
    face_distribution = get_face_emission_distribution(
        phi_in, theta_in, polarization_0, polarization_1, waveguide_0_df, waveguide_1_df, schema, sidecar
    )

    # Stage 3: Get normalized angular emission distribution with polarization
    angular_distribution = get_angular_emission_distribution(
        phi_in, theta_in, polarization_0, polarization_1, far_field_0_df, far_field_1_df, schema, sidecar
    )

    # Sample an emission point on the outgoing face
    sampled_point, sampled_point_index = sample_emission_point(face_distribution)

    # Sample an emission angle (Theta, Phi)
    sampled_angle, sampled_angle_index = sample_emission_angle(angular_distribution)

    # Extract emission point probability
    prob_point = face_distribution.loc[
        (face_distribution["X"] == sampled_point["X"]) &
        (face_distribution["Y"] == sampled_point["Y"]) &
        (face_distribution["Z"] == sampled_point["Z"]),
        "Probability"
    ].values[0]

    # Extract angular emission probability and polarization
    angle_row = angular_distribution.iloc[sampled_angle_index]
    prob_angle = angle_row["Probability"]
    polarization = angle_row["Polarization"]

    return {
        "S21_squared": S21_squared,
        "reached_end": True,
        "sampled_point": sampled_point,
        "point_probability": prob_point,
        "face_distribution": face_distribution,
        "sampled_angle": sampled_angle,
        "angle_probability": prob_angle,
        "angular_distribution": angular_distribution,
        "polarization_at_sampled_angle": polarization,
    }

def simulate_multiple_emissions(
    N,
    phi_in,
    theta_in,
    polarization_0,
    polarization_1,
    waveguide_0_df,
    waveguide_1_df,
    far_field_0_df,
    far_field_1_df,
    schema="1.0",
    sidecar=None
):
    """
    Run multiple photon emission simulations for a fixed incoming wave configuration.

    This function calls `simulate_photon_emission` N times and stores the results
    for statistical analysis or Monte Carlo estimation.

    Parameters
    ----------
    N : int
        Number of photons to simulate.
    phi_in : float
        Incident azimuthal angle (degrees).
    theta_in : float
        Incident polar angle (degrees).
    polarization_0 : str
        Polarization state label for waveguide 0 input.
    polarization_1 : str
        Polarization state label for waveguide 1 input.
    waveguide_0_df : pandas.DataFrame
        Data for waveguide 0 containing S-parameters and spatial emission info.
    waveguide_1_df : pandas.DataFrame
        Data for waveguide 1 containing S-parameters and spatial emission info.
    far_field_0_df : pandas.DataFrame
        Far-field angular distribution data for waveguide 0.
    far_field_1_df : pandas.DataFrame
        Far-field angular distribution data for waveguide 1.
    schema : str
        "1.0" (default) or "2.x", the latter with its dataset ``sidecar``.

    Returns
    -------
    list of dict
        A list containing the results of each simulation, as returned by
        `simulate_photon_emission`.
    """
    results = []
    for _ in range(N):
        result = simulate_photon_emission(
            phi_in,
            theta_in,
            polarization_0,
            polarization_1,
            waveguide_0_df,
            waveguide_1_df,
            far_field_0_df,
            far_field_1_df,
            schema,
            sidecar
        )
        results.append(result)
    return results


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(description="Sample photon emission through a gap from HFSS exports.")
    parser.add_argument("ephi0_dir", help="<id>_<freq>GHz_Ephi=0 directory with waveguide.csv and far_field.csv")
    parser.add_argument("ephi1_dir", help="<id>_<freq>GHz_Ephi=1 directory")
    parser.add_argument("--phi-in", type=float, default=0.0)
    parser.add_argument("--theta-in", type=float, default=180.0)
    parser.add_argument("--polarization", type=float, nargs=2, default=(1.0, 0.0), metavar=("P0", "P1"),
                        help="Amplitudes of (E_theta=1, E_phi=0) and (E_theta=0, E_phi=1); squares must sum to 1")
    parser.add_argument("--n", type=int, default=1, help="Number of photons to sample")
    args = parser.parse_args(argv)

    waveguide_0_df, waveguide_1_df, far_field_0_df, far_field_1_df, sidecar = load_datasets(args.ephi0_dir,
                                                                                            args.ephi1_dir)
    schema = sidecar["schema_version"] if sidecar is not None else "1.0"
    polarization_0, polarization_1 = args.polarization
    magnitude_squared = polarization_0**2 + polarization_1**2
    assert np.isclose(magnitude_squared, 1.0), f"Magnitude squared is {magnitude_squared}, expected 1."
    phi_in, theta_in = args.phi_in, args.theta_in

    # Stage 1
    power_transmission_fraction = get_S21(
        phi_in, theta_in, polarization_0, polarization_1, waveguide_0_df, waveguide_1_df, schema, sidecar
    )
    print(power_transmission_fraction)

    if power_transmission_fraction > 0:
        # Stage 2
        face_emission_distribution = get_face_emission_distribution(
            phi_in, theta_in, polarization_0, polarization_1, waveguide_0_df, waveguide_1_df, schema, sidecar
        )
        print(face_emission_distribution.sort_values("Probability", ascending=False).head())

        # Stage 3
        angular_emission_distribution = get_angular_emission_distribution(
            phi_in, theta_in, polarization_0, polarization_1, far_field_0_df, far_field_1_df, schema, sidecar
        )
        print(
            angular_emission_distribution
            .assign(
                Polarization=lambda df: df["Polarization"].apply(
                    lambda arr: np.array2string(arr, precision=4, separator=',', suppress_small=True))
            )[["Theta", "Phi", "Probability", "Polarization"]]
            .sort_values("Probability", ascending=False)
            .head()
        )
    else:
        print("No transmission at this key and polarization: nothing is sampled.")

    N = args.n
    all_results = simulate_multiple_emissions(
        N, phi_in, theta_in, polarization_0, polarization_1,
        waveguide_0_df, waveguide_1_df, far_field_0_df, far_field_1_df, schema, sidecar
    )
    num_reached = sum(1 for r in all_results if r["reached_end"])
    print(f"Out of {N} photons, {num_reached} reached the end face.")
    print(f"Transmission rate: {num_reached / N:.4f}")
    reached = [r for r in all_results if r["reached_end"]]
    points = pd.DataFrame([{
        "X": r["sampled_point"]["X"], "Y": r["sampled_point"]["Y"], "Z": r["sampled_point"]["Z"],
        "Probability": r["point_probability"]
    } for r in reached])
    angles = pd.DataFrame([{
        "Theta": r["sampled_angle"]["Theta"], "Phi": r["sampled_angle"]["Phi"],
        "Probability": r["angle_probability"]
    } for r in reached])
    print(points.head())
    print(angles.head())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
