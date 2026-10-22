"""Names shared between the runner, the exports and the Geant4 directory contract."""
from __future__ import annotations


def frequency_label(frequency_ghz: float) -> str:
    """``500`` and ``500.0`` -> ``"500GHz"``; ``512.5`` -> ``"512.5GHz"``; ``1234.5678`` -> ``"1234.5678GHz"``.

    The label is also the HFSS solve frequency, so it keeps every digit (shortest round-trip repr).
    """
    value = float(frequency_ghz)
    if value.is_integer():
        return f"{int(value)}GHz"
    return f"{value!r}GHz"


def frequency_design_name(base_design: str, frequency_ghz: float) -> str:
    return f"{base_design}_{frequency_label(frequency_ghz)}"


def setup_name(frequency_ghz: float) -> str:
    return frequency_label(frequency_ghz)


def parametric_setup_name(frequency_ghz: float) -> str:
    return f"E_phi_sweep_{frequency_label(frequency_ghz)}"


def plane_wave_name(frequency_ghz: float) -> str:
    return f"plane_wave_{frequency_label(frequency_ghz)}"


def radiation_sphere_name(frequency_ghz: float) -> str:
    return f"radiation_sphere_{frequency_label(frequency_ghz)}"


def last_adaptive_solution(frequency_ghz: float) -> str:
    return f"{setup_name(frequency_ghz)} : LastAdaptive"


def dataset_dir_name(dataset_id: str, frequency_ghz: float, ephi: int) -> str:
    """Directory stem read by BBRCrackLibrary: ``<id>_<freq>GHz_Ephi=<n>``."""
    return f"{dataset_id}_{frequency_label(frequency_ghz)}_Ephi={int(ephi)}"


def dataset_record_name(dataset_id: str, frequency_ghz: float) -> str:
    """Sidecar beside the two dataset directories: ``<id>_<freq>GHz.dataset.json``."""
    return f"{dataset_id}_{frequency_label(frequency_ghz)}.dataset.json"
