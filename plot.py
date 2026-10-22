import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def load_and_pivot(csv_path):
    """
    Load a CSV file containing waveguide simulation data and create a pivot table
    of the outgoing power indexed by incident angles.

    Parameters:
        csv_path (str): Path to the CSV file to load.

    Returns:
        pandas.DataFrame: A pivot table with 'IWaveTheta' as the index, 'IWavePhi' as the columns,
                          and 'OutgoingPower' as the values. Duplicate angle pairs are removed
                          before pivoting.
    """
    df = pd.read_csv(csv_path)
    df_unique = df.drop_duplicates(subset=["IWavePhi", "IWaveTheta"])
    pivot = df_unique.pivot(index="IWaveTheta", columns="IWavePhi", values="OutgoingPower")
    return pivot

def plot_outgoing_power_heatmap(pivot_table, title, cmap="viridis"):
    """
    Plot heatmap with centered ticks every 15 degrees and larger fonts.
    """
    plt.figure(figsize=(8, 6))

    ax = sns.heatmap(
        pivot_table,
        annot=False,
        fmt=".2e",
        cmap=cmap,
        cbar_kws={"label": r"$P_{out}$ (W)"}
    )

    theta_vals = pivot_table.index.astype(float).to_numpy()
    phi_vals   = pivot_table.columns.astype(float).to_numpy()

    # uniform 15° ticks
    theta_ticks = np.arange(np.floor(theta_vals.min()/15)*15,
                            np.ceil(theta_vals.max()/15)*15 + 1, 15)
    phi_ticks   = np.arange(np.floor(phi_vals.min()/15)*15,
                            np.ceil(phi_vals.max()/15)*15 + 1, 15)

    theta_tick_positions = np.interp(theta_ticks, theta_vals, np.arange(len(theta_vals)))
    phi_tick_positions   = np.interp(phi_ticks, phi_vals, np.arange(len(phi_vals)))

    # center ticks
    ax.set_yticks(theta_tick_positions + 0.5)
    ax.set_xticks(phi_tick_positions + 0.5)

    ax.set_yticklabels(theta_ticks)
    ax.set_xticklabels(phi_ticks)

    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)

    ax.set_title(title, fontsize=18)
    ax.set_xlabel(r"$\phi_{in}$ (degrees)", fontsize=16)
    ax.set_ylabel(r"$\theta_{in}$ (degrees)", fontsize=16)

    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=15)
    cbar.set_label(r"$P_{out}$ (W)", fontsize=16)
    cbar.ax.yaxis.get_offset_text().set_fontsize(16)

    plt.tight_layout()
    
def plot_outgoing_power_by_incoming_angle(csv_path, fixed_param, fixed_values):
    """
    Plot outgoing power as a function of one incoming angle (theta_in or phi_in)
    for fixed values of the other angle.

    Parameters:
        csv_path (str): Path to the CSV file containing the waveguide data.
        fixed_param (str): The angle to fix, either r"$\\theta_{in}$" or r"$\\phi_{in}$".
        fixed_values (float or list of floats): One or more fixed values for the fixed_param angle (in degrees).

    Raises:
        ValueError: If fixed_param is not one of the allowed strings.

    Returns:
        None. Displays a line plot.
    """
    param_map = {
        r"$\theta_{in}$": "IWaveTheta",
        r"$\phi_{in}$": "IWavePhi"
    }
    if fixed_param not in param_map:
        raise ValueError("fixed_param must be r'$\\theta_{in}$' or r'$\\phi_{in}$'")
    
    col_fix = param_map[fixed_param]
    col_var = "IWavePhi" if col_fix == "IWaveTheta" else "IWaveTheta"
    xlabel_units = r"$\phi_{in}$ (degrees)" if col_fix == "IWaveTheta" else r"$\theta_{in}$ (degrees)"
    xlabel_no_units = r"$\phi_{in}$" if col_fix == "IWaveTheta" else r"$\theta_{in}$"

    df = pd.read_csv(csv_path)
    if not isinstance(fixed_values, (list, tuple)):
        fixed_values = [fixed_values]

    plt.figure(figsize=(8, 6))
    for val in fixed_values:
        df_filt = df[df[col_fix] == val]
        plt.plot(df_filt[col_var], df_filt["OutgoingPower"], marker='o', markersize=4, label=f"{fixed_param} = {val}°")

    plt.xlabel(xlabel_units)
    plt.ylabel(r"$P_{out}$ (W)")
    plt.title(rf"$P_{{out}}$ vs {xlabel_no_units}")
    plt.grid(True)
    if len(fixed_values) > 1:
        plt.legend()
    plt.tight_layout()

def plot_exit_field_by_incoming_angle(csv_path, theta_in, phi_in,
                                     fixed_coord="Z", fixed_value=0.001,
                                     x_axis="X", y_axis="Y",
                                     s=10, cmap="viridis"):
    """
    Plot the electric field magnitude on a fixed exit plane slice as a function of position,
    for specified incoming angles θ_in and φ_in.

    Parameters:
        csv_path (str): Path to the CSV file containing the exit field data.
        theta_in (float): Incoming polar angle θ_in in degrees.
        phi_in (float): Incoming azimuthal angle φ_in in degrees.
        fixed_coord (str, optional): The coordinate to fix (e.g., 'Z'). Default is 'Z'.
        fixed_value (float, optional): The fixed value of the coordinate (e.g., plane location). Default is 0.001 (meters).
        x_axis (str, optional): Coordinate to plot on x-axis. Default is 'X'.
        y_axis (str, optional): Coordinate to plot on y-axis. Default is 'Y'.
        s (int, optional): Marker size for scatter plot. Default is 10.
        cmap (str, optional): Matplotlib colormap name. Default is 'viridis'.

    Returns:
        None. Displays a scatter plot of |E| magnitude on the exit plane.
    """
    df = pd.read_csv(csv_path)
    # Filter data for the given incoming angles and fixed coordinate slice (e.g., z=0.001)
    df_filt = df[(df["IWaveTheta"] == theta_in) & 
                 (df["IWavePhi"] == phi_in) & 
                 (np.isclose(df[fixed_coord], fixed_value))]

    if df_filt.empty:
        print(f"No data for θ={theta_in}°, φ={phi_in}° at {fixed_coord.lower()} = {fixed_value}")
        return

    E_mag = np.sqrt(df_filt["Ex_real"]**2 + df_filt["Ey_real"]**2 + df_filt["Ez_real"]**2 + df_filt["Ex_imag"]**2 + df_filt["Ey_imag"]**2 + df_filt["Ez_imag"]**2)

    # Scatter plot of the E-field over the output plane
    plt.figure(figsize=(7, 6))
    sc = plt.scatter(df_filt[x_axis], df_filt[y_axis], c=E_mag, cmap=cmap, s=s)
    plt.xlabel(f"{x_axis.lower()} (m)")
    plt.ylabel(f"{y_axis.lower()} (m)")
    plt.title(rf"$|E|$ at {fixed_coord.lower()} = {fixed_value}m for $\theta_{{in}}$ = {theta_in}°, $\phi_{{in}}$ = {phi_in}°")
    plt.colorbar(sc, label=r"$|E|$ (V/m)")
    plt.tight_layout()

def plot_far_field_by_incoming_angle_fixed(
    csv_path,
    theta_in_list,
    phi_in_list,
    fixed_param,
    fixed_value,
    title_prefix="Far-field E-Field Magnitude"
):
    """
    Plot 1D slices of far-field electric field magnitude at a fixed outgoing angle,
    with improved font sizing and tick spacing every 15 degrees.
    """

    param_map = {
        r"$\theta_{out}$": "Theta",
        r"$\phi_{out}$": "Phi"
    }
    if fixed_param not in param_map:
        raise ValueError("fixed_param must be r'$\\theta_{out}$' or r'$\\phi_{out}$'")

    fixed_col = param_map[fixed_param]
    varying_col = "Phi" if fixed_col == "Theta" else "Theta"
    xlabel_units = r"$\phi_{out}$ (degrees)" if fixed_col == "Theta" else r"$\theta_{out}$ (degrees)"

    df = pd.read_csv(csv_path)

    if not isinstance(theta_in_list, (list, tuple)):
        theta_in_list = [theta_in_list]
    if not isinstance(phi_in_list, (list, tuple)):
        phi_in_list = [phi_in_list]

    rounded_fixed_value = round(fixed_value)
    plt.figure(figsize=(10, 7))

    for theta_in in theta_in_list:
        for phi_in in phi_in_list:

            df_filt = df[
                (df["IWaveTheta"] == theta_in) &
                (df["IWavePhi"] == phi_in) &
                (np.isclose(df[fixed_col], fixed_value))
            ]
            if df_filt.empty:
                continue

            # E field magnitude
            Ephi_mag = np.sqrt(df_filt["rEphi_real"]**2 + df_filt["rEphi_imag"]**2)
            Etheta_mag = np.sqrt(df_filt["rEtheta_real"]**2 + df_filt["rEtheta_imag"]**2)
            Etot = np.sqrt(Ephi_mag**2 + Etheta_mag**2)

            label = (
                fr"{fixed_param} = {rounded_fixed_value}°, "
                fr"$\theta_{{in}}$ = {theta_in}°, "
                fr"$\phi_{{in}}$ = {phi_in}°"
            )

            plt.plot(df_filt[varying_col], Etot, label=label)

    plt.xlabel(xlabel_units, fontsize=16)
    plt.ylabel(r"$|E|$ (V/m)", fontsize=16)
    plt.title(f"{title_prefix} at {fixed_param} = {rounded_fixed_value}°", fontsize=18)

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)


    var_vals = df[varying_col].unique().astype(float)
    tick_min = np.floor(var_vals.min() / 15) * 15
    tick_max = np.ceil(var_vals.max() / 15) * 15
    ticks = np.arange(tick_min, tick_max + 1, 15)
    plt.xticks(ticks)

    plt.legend(fontsize=12, loc='best')

    plt.grid(True)
    plt.tight_layout()

def plot_far_field_by_incoming_angle(csv_path,
                                     theta_in,
                                     phi_in,
                                     cmap="viridis",
                                     title_prefix="Far-field E-Field Magnitude",
                                     tick_interval=10):
    """
   Plot a 2D heatmap of far-field E-field magnitude as a function of outgoing angles
   for a specific incoming angle (theta_in, phi_in).

   Parameters:
       csv_path (str): Path to CSV file with far-field data.
       theta_in (float): Incoming polar angle in degrees.
       phi_in (float): Incoming azimuthal angle in degrees.
       cmap (str): Matplotlib colormap name.
       title_prefix (str): Title prefix for the plot.
       tick_interval (int): Interval between axis ticks in degrees.

   Returns:
       None. Displays a heatmap plot.
   """
    df = pd.read_csv(csv_path)

    # Filter for specified incoming angles
    df_filt = df[(df["IWaveTheta"] == theta_in) & (df["IWavePhi"] == phi_in)]

    if df_filt.empty:
        print(f"No data for θ_in={theta_in}°, φ_in={phi_in}°")
        return

    # Compute total E-field magnitude from real/imag components
    Ephi_mag = np.sqrt(df_filt["rEphi_real"]**2 + df_filt["rEphi_imag"]**2)
    Etheta_mag = np.sqrt(df_filt["rEtheta_real"]**2 + df_filt["rEtheta_imag"]**2)
    Etot = np.sqrt(Ephi_mag**2 + Etheta_mag**2)
    df_filt = df_filt.copy()
    df_filt["|E|"] = Etot

    # Pivot to get 2D matrix: Theta_out vs Phi_out
    pivot = df_filt.pivot(index="Theta", columns="Phi", values="|E|")

    plt.figure(figsize=(8, 6))
    ax = sns.heatmap(pivot, cmap=cmap, cbar_kws={"label": r"$|E|$ (V/m)"})

    plt.xlabel(r"$\phi_{out}$ (degrees)")
    plt.ylabel(r"$\theta_{out}$ (degrees)")
    plt.title(rf"{title_prefix} for $\theta_{{in}}$ = {theta_in}°, $\phi_{{in}}$ = {phi_in}°")

    # Set neat tick intervals on axes, even if raw data has lots of decimal places
    x_ticks = np.arange(pivot.columns.min(), pivot.columns.max() + 1, tick_interval)
    y_ticks = np.arange(pivot.index.min(), pivot.index.max() + 1, tick_interval)

    x_tick_positions = [np.abs(pivot.columns - val).argmin() for val in x_ticks]
    y_tick_positions = [np.abs(pivot.index - val).argmin() for val in y_ticks]

    ax.set_xticks(x_tick_positions)
    ax.set_xticklabels([f"{int(round(val))}" for val in x_ticks])
    ax.set_yticks(y_tick_positions)
    ax.set_yticklabels([f"{int(round(val))}" for val in y_ticks])

    plt.tight_layout()


def main(argv=None):
    import argparse
    import os
    parser = argparse.ArgumentParser(description="Plot HFSS gap exports (waveguide.csv, far_field.csv).")
    parser.add_argument("dataset_dir", help="Directory containing waveguide.csv and far_field.csv")
    parser.add_argument("--out-dir", default=None, help="Save PNGs here instead of opening windows")
    parser.add_argument("--theta-in", type=float, default=180.0)
    parser.add_argument("--phi-in", type=float, default=0.0)
    args = parser.parse_args(argv)
    if args.out_dir:
        plt.switch_backend("Agg")

    waveguide_data = os.path.join(args.dataset_dir, "waveguide.csv")
    far_field_data = os.path.join(args.dataset_dir, "far_field.csv")

    pivot = load_and_pivot(waveguide_data)
    plot_outgoing_power_heatmap(pivot, r"$P_{out}$ vs $\theta_{in}$ and $\phi_{in}$")
    plot_outgoing_power_by_incoming_angle(waveguide_data, fixed_param=r"$\theta_{in}$", fixed_values=[args.theta_in])
    plot_exit_field_by_incoming_angle(csv_path=waveguide_data, theta_in=args.theta_in, phi_in=args.phi_in,
                                      fixed_coord="X", x_axis="Z", y_axis="Y", fixed_value=0)
    plot_far_field_by_incoming_angle_fixed(csv_path=far_field_data, theta_in_list=[args.theta_in],
                                           phi_in_list=[args.phi_in], fixed_param=r"$\theta_{out}$", fixed_value=90)
    plot_far_field_by_incoming_angle(csv_path=far_field_data, theta_in=args.theta_in, phi_in=args.phi_in)

    if args.out_dir:
        os.makedirs(args.out_dir, exist_ok=True)
        for i, num in enumerate(plt.get_fignums(), start=1):
            plt.figure(num).savefig(os.path.join(args.out_dir, f"figure_{i:02d}.png"), dpi=150)
    else:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
