"""CSV column contracts.

Both tables are read positionally by the Geant4 loader
(BBRSimulation/library/src/BBRHFSSData.cc). Never reorder these lists.
"""
from __future__ import annotations

import pandas as pd

WAVEGUIDE_COLUMNS = [
    "Freq", "Ephi", "IWavePhi", "IWaveTheta", "OutgoingPower", "IngoingPower",
    "X", "Y", "Z", "Ex_real", "Ey_real", "Ez_real", "Ex_imag", "Ey_imag", "Ez_imag",
]
FAR_FIELD_COLUMNS = [
    "Freq", "Ephi", "IWavePhi", "IWaveTheta", "Phi", "Theta",
    "rEphi_real", "rEphi_imag", "rEtheta_real", "rEtheta_imag",
]
KEY_COLUMNS = ["IWavePhi", "IWaveTheta"]
EXIT_FIELD_COLUMNS = ["X", "Y", "Z", "Ex_real", "Ey_real", "Ez_real", "Ex_imag", "Ey_imag", "Ez_imag"]
FAR_FIELD_POINT_COLUMNS = ["Phi", "Theta", "rEphi_real", "rEphi_imag", "rEtheta_real", "rEtheta_imag"]


class SchemaError(ValueError):
    """A table does not carry the required columns."""


def order_columns(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    """Return ``df`` restricted to ``columns`` in that order; raise if any is missing."""
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise SchemaError(f"missing columns: {missing}")
    return df.loc[:, list(columns)]
