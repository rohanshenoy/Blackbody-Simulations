import pandas as pd
import numpy as np

def load_datasets(ephi0_dir, ephi1_dir):
    """Load waveguide and far-field tables for the Ephi=0 and Ephi=1 datasets of one frequency."""
    import os
    return (
        pd.read_csv(os.path.join(ephi0_dir, "waveguide.csv")),
        pd.read_csv(os.path.join(ephi1_dir, "waveguide.csv")),
        pd.read_csv(os.path.join(ephi0_dir, "far_field.csv")),
        pd.read_csv(os.path.join(ephi1_dir, "far_field.csv")),
    )

#%% Stage 1: Finding the transmission power fraction of the waveguide

def get_S21(phi_in, theta_in, polarization_0, polarization_1, df_0, df_1):
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

    Returns:
        float: Weighted and clipped |S21|^2 (P_out / P_in, capped at 1)
    """

    def lookup_power(df, phi, theta):
        row = df[(df["IWavePhi"] == phi) & (df["IWaveTheta"] == theta)]
        if row.empty:
            raise ValueError(f"No data for IWavePhi={phi}, IWaveTheta={theta}")
        P_out = row["OutgoingPower"].values[0]
        P_in = row["IngoingPower"].values[0]
        if P_in == 0:
            raise ZeroDivisionError("IngoingPower is zero.")
        return min(P_out / P_in, 1)

    transmission_0 = lookup_power(df_0, phi_in, theta_in)
    transmission_1 = lookup_power(df_1, phi_in, theta_in)

    weighted_S21_sq = (
        polarization_0**2 * transmission_0 +
        polarization_1**2 * transmission_1
    )

    return weighted_S21_sq
#%% Stage 2: Determining a probability distribution over the location of emission from the outgoing face.

def get_face_emission_distribution(phi_in, theta_in, polarization_0, polarization_1,
                              df_0, df_1):
    """
    Compute normalized emission probability distribution over the output face
    for a given incoming angle and polarization, using linear superposition of E-fields.

    Parameters:
        phi_in (float): IWavePhi in degrees
        theta_in (float): IWaveTheta in degrees
        polarization_0 (float): Amplitude of polarization mode 0 (Ephi=0)
        polarization_1 (float): Amplitude of polarization mode 1 (Ephi=1)
        df_0 (pd.DataFrame): E-field data for (Etheta=1, Ephi=0)
        df_1 (pd.DataFrame): E-field data for (Etheta=0, Ephi=1)

    Returns:
        pd.DataFrame: With columns [X, Y, Z, E2, Probability]
    """

    # Filter to the relevant angle in both datasets
    f0 = df_0[(df_0["IWavePhi"] == phi_in) & (df_0["IWaveTheta"] == theta_in)].copy()
    f1 = df_1[(df_1["IWavePhi"] == phi_in) & (df_1["IWaveTheta"] == theta_in)].copy()

    if f0.empty or f1.empty:
        raise ValueError(f"No data for IWavePhi={phi_in}, IWaveTheta={theta_in}")

    # Ensure same spatial points
    if not f0[["X", "Y", "Z"]].equals(f1[["X", "Y", "Z"]]):
        raise ValueError("Mismatch in spatial coordinates between df_0 and df_1.")

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

#%% Stage 3: Determining a probability distribution over the angle of emission parametrized by (theta, phi)

def get_angular_emission_distribution(phi_in, theta_in, polarization_0, polarization_1,
                                      df_0, df_1):
    """
    Compute normalized emission probability distribution over outgoing angles (theta, phi)
    for given incoming angle and polarization, using linear superposition of E-fields
    in spherical coordinates (rEtheta, rEphi).

    Parameters:
        phi_in (float): Incoming IWavePhi in degrees
        theta_in (float): Incoming IWaveTheta in degrees
        polarization_0 (float): Amplitude of polarization mode 0 (Ephi=0)
        polarization_1 (float): Amplitude of polarization mode 1 (Ephi=1)
        df_0 (pd.DataFrame): Far-field E data for (Etheta=1, Ephi=0)
        df_1 (pd.DataFrame): Far-field E data for (Etheta=0, Ephi=1)

    Returns:
        pd.DataFrame: With columns [Theta, Phi, E2, Probability, Polarization]
                     Polarization column contains normalized complex vectors [Etheta, Ephi]
    """

    # Filter to the relevant incoming angles
    f0 = df_0[(df_0["IWavePhi"] == phi_in) & (df_0["IWaveTheta"] == theta_in)].copy()
    f1 = df_1[(df_1["IWavePhi"] == phi_in) & (df_1["IWaveTheta"] == theta_in)].copy()

    if f0.empty or f1.empty:
        raise ValueError(f"No data for IWavePhi={phi_in}, IWaveTheta={theta_in}")

    # Ensure angle match
    if not f0[["Theta", "Phi"]].equals(f1[["Theta", "Phi"]]):
        raise ValueError("Mismatch in angular coordinates between df_0 and df_1.")

    # Reconstruct spherical components
    Etheta_0 = f0["rEtheta_real"] + 1j * f0["rEtheta_imag"]
    Ephi_0   = f0["rEphi_real"]   + 1j * f0["rEphi_imag"]

    Etheta_1 = f1["rEtheta_real"] + 1j * f1["rEtheta_imag"]
    Ephi_1   = f1["rEphi_real"]   + 1j * f1["rEphi_imag"]

    # Linear superposition in spherical basis
    Etheta_total = polarization_0 * Etheta_0 + polarization_1 * Etheta_1
    Ephi_total   = polarization_0 * Ephi_0   + polarization_1 * Ephi_1

    # Compute |E|^2 = |Etheta|^2 + |Ephi|^2
    E2 = np.abs(Etheta_total)**2 + np.abs(Ephi_total)**2

    # Normalize to get probability distribution
    total_E2 = E2.sum()
    if total_E2 == 0:
        raise ValueError("Total E-field magnitude squared is zero — can't normalize.")

    probabilities = E2 / total_E2

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
    far_field_1_df
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
        phi_in, theta_in, polarization_0, polarization_1, waveguide_0_df, waveguide_1_df
    )

    # Decide if the photon reaches the end face
    reached_end = np.random.rand() < S21_squared
    if not reached_end:
        # Photon is lost: nothing else to simulate
        return {
            "S21_squared": S21_squared,
            "reached_end": False
        }

    # Stage 2: Get normalized spatial emission distribution on outgoing face
    face_distribution = get_face_emission_distribution(
        phi_in, theta_in, polarization_0, polarization_1, waveguide_0_df, waveguide_1_df
    )

    # Stage 3: Get normalized angular emission distribution with polarization
    angular_distribution = get_angular_emission_distribution(
        phi_in, theta_in, polarization_0, polarization_1, far_field_0_df, far_field_1_df
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
    far_field_1_df
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
            far_field_1_df
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

    waveguide_0_df, waveguide_1_df, far_field_0_df, far_field_1_df = load_datasets(args.ephi0_dir, args.ephi1_dir)
    polarization_0, polarization_1 = args.polarization
    magnitude_squared = polarization_0**2 + polarization_1**2
    assert np.isclose(magnitude_squared, 1.0), f"Magnitude squared is {magnitude_squared}, expected 1."
    phi_in, theta_in = args.phi_in, args.theta_in

    # Stage 1
    power_transmission_fraction = get_S21(
        phi_in, theta_in, polarization_0, polarization_1, waveguide_0_df, waveguide_1_df
    )
    print(power_transmission_fraction)

    # Stage 2
    face_emission_distribution = get_face_emission_distribution(
        phi_in, theta_in, polarization_0, polarization_1, waveguide_0_df, waveguide_1_df
    )
    print(face_emission_distribution.sort_values("Probability", ascending=False).head())

    # Stage 3
    angular_emission_distribution = get_angular_emission_distribution(
        phi_in, theta_in, polarization_0, polarization_1, far_field_0_df, far_field_1_df
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

    N = args.n
    all_results = simulate_multiple_emissions(
        N, phi_in, theta_in, polarization_0, polarization_1,
        waveguide_0_df, waveguide_1_df, far_field_0_df, far_field_1_df
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
