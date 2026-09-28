# -*- coding: utf-8 -*-
"""
Created on Fri Aug 15 13:08:44 2025

@author: Jason Wang
"""

import os
import pandas as pd
import numpy as np
from scipy.interpolate import RegularGridInterpolator, griddata
import matplotlib.pyplot as plt

def filter_by_variables(predictor_variables, response_variables, df):
    # Keep only predictor + response columns
    df_filtered = df[predictor_variables + response_variables].drop_duplicates()
    
    return df_filtered    

def interpolate_and_save(df, predictor_variables, response_variables, output_folder, resolutions, method='linear', save_to_csv=False):
    """
    Interpolate response(s) over all predictor variables and optionally save as CSV.
    
    Parameters:
        df: pandas DataFrame with predictors + responses
        predictor_variables: list of predictor column names
        response_variables: list of response column names to interpolate
        output_folder: folder to save CSV
        resolutions: dict mapping predictor names to grid resolution (same units as df)
        method: interpolation method ('linear', 'nearest', 'cubic')
        save_to_csv: if True, save the interpolated data to CSV
    """

    # Create 1D arrays for each predictor
    grid_axes = []
    for p in predictor_variables:
        min_val = df[p].min()
        max_val = df[p].max()
        res = resolutions.get(p, 1)
        axis = np.arange(min_val, max_val + res, res)
        grid_axes.append(axis)
        print(f"{p}: {len(axis)} points from {min_val} to {max_val} with resolution {res}")


    # Create N-D grid
    mesh = np.meshgrid(*grid_axes, indexing='ij')
    mesh_points = np.stack([m.ravel() for m in mesh], axis=-1)

    # Interpolate each response
    interpolated_data = {}
    points = df[predictor_variables].values
    for resp in response_variables:
        values = df[resp].values
        grid_values = griddata(points, values, mesh_points, method=method)
        interpolated_data[resp] = grid_values

    # Build DataFrame
    export_dict = {p: mesh_points[:, i] for i, p in enumerate(predictor_variables)}
    for resp in response_variables:
        export_dict[resp] = interpolated_data[resp]

    export_df = pd.DataFrame(export_dict)

    # Save CSV if requested
    if save_to_csv:
        output_csv_path = os.path.join(output_folder, "interpolated.csv")
        export_df.to_csv(output_csv_path, index=False)
        print(f"Saved interpolated data to {output_csv_path}")

    return export_df  # return the interpolated DataFrame for further use


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(description="Interpolate refined S21 grids onto a fine angular grid.")
    parser.add_argument("--sim-data-root",
                        default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "HFSSSimData"))
    parser.add_argument("--project", default="InfParallelPlate")
    parser.add_argument("--design", default="interpolate")
    parser.add_argument("--frequencies", type=int, nargs="+", default=[500, 550, 600])
    args = parser.parse_args(argv)

    for freq in args.frequencies:
        for Ephi in [0, 1]:
            freq_sim_data_location = f"{args.project}_{args.design}_{freq}GHz_Ephi={Ephi}"
            freq_sim_data_folder = os.path.join(args.sim_data_root, freq_sim_data_location)
            waveguide_path = os.path.join(freq_sim_data_folder, "refined_waveguide.csv")
            far_field_path = os.path.join(freq_sim_data_folder, "refined_far_field.csv")

            if not os.path.exists(waveguide_path):
                print(f"Missing file: {waveguide_path}")
                continue  # the legacy loop crashed on an undefined DataFrame here
            waveguide_df = pd.read_csv(waveguide_path)
            print(f"Loaded refined_waveguide.csv for freq={freq} GHz, Ephi={Ephi}")
            # Add a outgoing power fraction column
            waveguide_df["OutgoingPowerFraction"] = (
                waveguide_df["OutgoingPower"] / waveguide_df["IngoingPower"]
            ).clip(upper=1.0)

            if os.path.exists(far_field_path):
                far_field_df = pd.read_csv(far_field_path)
                print(f"Loaded refined_far_field.csv for freq={freq} GHz, Ephi={Ephi}")
            else:
                print(f"Missing file: {far_field_path}")
                far_field_df = None

            S21_predictor_variables = ["IWaveTheta", "IWavePhi"]
            S21_response_variables = ["OutgoingPowerFraction"]
            S21_df = filter_by_variables(S21_predictor_variables, S21_response_variables, waveguide_df)

            S21_resolutions = {
                "IWaveTheta": 0.1,  # degrees
                "IWavePhi": 0.1,  # degrees
            }
            interpolate_and_save(
                df=S21_df,
                predictor_variables=S21_predictor_variables,
                response_variables=S21_response_variables,
                output_folder=freq_sim_data_folder,
                resolutions=S21_resolutions,
                method='linear',
                save_to_csv=True,
            )
            # The exit-field and far-field interpolations in the original script were commented out
            # as too slow ("POD + modal coefficients perhaps"); far_field_df is loaded for that future work.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
