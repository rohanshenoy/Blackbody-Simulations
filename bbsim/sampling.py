"""Angular sampling rules and incident power, ported unchanged from legacy/bbsim1freq.py."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

SPEED_OF_LIGHT_M_PER_S = 299792458.0
# CODATA 2022 (scipy >= 1.15). The Windows reference IngoingPower was computed with this value;
# CODATA 2018 (1.25663706212e-6) differs by 7e-10 relative.
VACUUM_PERMEABILITY_H_PER_M = 1.25663706127e-6


def incident_point_count(lower_deg: float, upper_deg: float, step_deg: float) -> int:
    """Number of plane-wave sweep points for one angle (legacy get_incoming_phi_theta_num)."""
    span = float(upper_deg) - float(lower_deg)
    if span < 0:
        raise ValueError(f"upper ({upper_deg}) must not be below lower ({lower_deg})")
    if span == 0:
        return 1
    if step_deg <= 0:
        raise ValueError(f"step must be positive, got {step_deg}")
    return int(np.ceil(span / step_deg)) + 1


def incident_angle_values(lower_deg: float, upper_deg: float, step_deg: float) -> np.ndarray:
    return np.linspace(float(lower_deg), float(upper_deg), incident_point_count(lower_deg, upper_deg, step_deg))


def far_field_step_deg(range_deg: float, length_mm: float, frequency_ghz: float, fineness: float,
                       min_coarseness_deg: float, max_coarseness_deg: float) -> float:
    """Angular step = range/pi * lambda / (fineness * length), clipped to [min, max] coarseness."""
    wavelength_m = SPEED_OF_LIGHT_M_PER_S / (float(frequency_ghz) * 1e9)
    step = min(max_coarseness_deg, range_deg / np.pi * wavelength_m / (fineness * (length_mm * 1e-3)))
    return float(max(min_coarseness_deg, step))


@dataclass(frozen=True)
class FarFieldGrid:
    theta_start_deg: float
    theta_stop_deg: float
    theta_step_deg: float
    phi_start_deg: float
    phi_stop_deg: float
    phi_step_deg: float

    @property
    def theta_values(self) -> np.ndarray:
        # Same arange rule as the .ffd reader; HFSS emits exactly these points.
        return np.arange(self.theta_start_deg, self.theta_stop_deg + 1e-8, self.theta_step_deg)

    @property
    def phi_values(self) -> np.ndarray:
        return np.arange(self.phi_start_deg, self.phi_stop_deg + 1e-8, self.phi_step_deg)

    @property
    def points_per_angle(self) -> int:
        return len(self.theta_values) * len(self.phi_values)


def far_field_grid(theta_lower_deg: float, theta_upper_deg: float, phi_lower_deg: float, phi_upper_deg: float,
                   a_mm: float, b_mm: float, frequency_ghz: float, fineness: float,
                   min_coarseness_deg: float, max_coarseness_deg: float) -> FarFieldGrid:
    """Theta resolution follows the short dimension b, phi resolution the long dimension a."""
    theta_step = far_field_step_deg(theta_upper_deg - theta_lower_deg, b_mm, frequency_ghz, fineness,
                                    min_coarseness_deg, max_coarseness_deg)
    phi_step = far_field_step_deg(phi_upper_deg - phi_lower_deg, a_mm, frequency_ghz, fineness,
                                  min_coarseness_deg, max_coarseness_deg)
    return FarFieldGrid(theta_lower_deg, theta_upper_deg, theta_step, phi_lower_deg, phi_upper_deg, phi_step)


def incoming_power_w(ei_v_per_m: float, face_area_mm2: float) -> float:
    """Plane-wave power through the entrance face: area * Ei^2 / (2 c mu0), area converted from mm^2."""
    return 1e-6 * face_area_mm2 * ei_v_per_m ** 2 / (2 * SPEED_OF_LIGHT_M_PER_S * VACUUM_PERMEABILITY_H_PER_M)
