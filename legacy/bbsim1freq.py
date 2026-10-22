# -*- coding: utf-8 -*-
"""
Created on 03/20/2025

@author: Jason Wang
"""

# Notes:
# (1) The code takes in a command-line argument for the frequency.
# (2) The code creates a copy of the base design for each frequency. If a mistake is made and the code needs to be
# rerun for a frequency for a particular design, it is advisable to rename the original design due to a technicality in HFSS
# (3) If is desired, it is possible to Keyboard Interrupt in the middle of a simulation and save to csv. Then,
# it is possible to begin a new simulation with this existing data. 
# (4) The output .csv files are large, so to use the data, import directly using pandas.
# (5) Make sure there are no additional variables defined in the HFSS GUI! This can lead to errors in exporting far field data due
# to a technicality

# Packages: pyaedt, scipy, pandas, numpy, seaborn (for plotting)

import sys
import threading

# Make sure to install the package pyaedt!
from ansys.aedt.core.hfss import Hfss
from scipy import constants
import pandas as pd
import numpy as np
import os
import re
import time
from collections import defaultdict
from typing import Optional

#%%
frequency = int(sys.argv[1]) # Frequency [GHz] (command-line argument)
print("Frequency: " + str(frequency) + "GHz") # first argument after the script name; the frequency
c = constants.c # Speed of light
mu_0 = constants.mu_0 # permeability of free space

project_name = "InfParallelPlate" # Name of the HFSS project
design_name = "adaptive_plot_6" # Name of the HFSS design
repo_root = os.path.dirname(os.path.abspath(__file__))
output_file_location = os.path.join(repo_root, "HFSSSimData") # The folder to output all data files
os.makedirs(output_file_location, exist_ok=True)

import_from_existing_csv = False # Whether to import waveguide and far_field_data from existing CSV (default is false)
waveguide_data_csv = os.path.join(output_file_location, "InfParallelPlate_adaptive_plot_500GHz_Ephi=1/waveguide.csv")
far_field_data_csv = os.path.join(output_file_location, "InfParallelPlate_adaptive_plot_500GHz_Ephi=1/far_field.csv")

Ei = 1 # strength of incident electric field [V/m]
max_delta_E = 0.02 # HFSS analysis sweep parameter
max_passes = 10 # HFSS analysis sweep parameter
num_cores = 4

ingoing_face_id = 8 # Check face id's of the plane wave ingoing face and outgoing face by clicking select object/by name
outgoing_face_id = 7

# Define x and y directions of outgoing coordinate systems (vectors relative to global coordinate system)
# x direction points outward from face. The z direction is automatic from the right-hand rule.
# It is helpful to redefine a coordinate system so that the theta and phi sweep correspond to sweeps corresponding
# to the two length scales. The 4 angular variables refer to sweeps over the far field radiation, with respect to the user-defined CS
outgoing_face_cs_x = [0, 0, 1]
outgoing_face_cs_y = [0, 1, 0]
rad_theta_lower, rad_theta_upper, rad_phi_lower, rad_phi_upper = 0, 180, -90, 90

# Whether to specify manually the boundary and resolution of the outgoing face. If manual_field is set to False,
# HFSS will infer the appropriate sampling resolution, but it is coarse. The coordinates are relative to the global CS.
# If manual_field is set to True, the coordinates are relative to outgoing CS. If outgoing_face_boundary is set to None,
# the bounding region for the outgoing face is inferred, but it only works for rectilinear faces.
# outgoing_face_boundary can also be manually specified for non-rectilinear geometries.
manual_field = True
outgoing_face_boundary = None # Or specify ([x1, y1, z1], [x2, y2, z2]), e.g.([0, -5, -0.025], [0, 5, 0.025]) Boundary of the outgoing face (mm)
outgoing_face_field_resolution = ["0mm", "0.1mm", "0.001mm"] # Resolution in outputting the electric field at the outgoing face

# The following 4 variables refer to sweeps over incident plane wave. These angles are with respect to the global
# coordinate system shfited to the face. Symmetry can be used to make these sweeps less wide
i_theta_lower, i_theta_upper, i_phi_lower, i_phi_upper = 90, 180, 0, 90
i_theta_step = 15 # Initial step size over theta and phi (adaptive), or step size over theta and phi (discrete)
i_phi_step = 15

a = 10 # length scale of the dimension coinciding with the sweep over phi [mm]
b = 0.05 # length scale of the dimension coinciding with the sweep over theta [mm]
fineness = 10 #1/fineness is the fraction of lambda/a swept over each radiation step in theta and phi

# Maximum coarseness is the maximum coarseness of the angular sweeps, in degrees.
# For small widths (a or b small), phi or theta will have only 1 main lobe with angular width lambda/a very high.
# The general guideline of sampling is 10 points across the narrowest feature
minimum_coarseness = 0.1 # degrees
maximum_coarseness = 5 # degrees

# Collect radiation parameters into one variable for cleanness
rad_params = rad_theta_lower, rad_theta_upper, rad_phi_lower, rad_phi_upper, a, b, fineness, maximum_coarseness, minimum_coarseness

# Adaptive or discrete sweep. For adaptive sweep, max difference is maximum fractional difference allowed between
# any two points in the sweep, relative to the total maximum value of the outgoing power.
sweep = "adaptive"
max_difference = 0.25

#%%
# Instantiate HFSS session and set the project and design appropriately
hfss = Hfss(project=project_name, design=design_name, non_graphical=False)
oDesktop = hfss.odesktop
oProject = oDesktop.SetActiveProject(project_name)
oDesign = oProject.SetActiveDesign(design_name)
oEditor = oDesign.SetActiveEditor("3D Modeler")

# For each frequency, a copy of the existing geometry is made, but HFSS sometimes renames
# face ids when copying a geometry from one design to another. Thus we ensure that the
# correct face id is acquired by checking the position of the face and finding the
# face with the same position in the copied design
ingoing_face = hfss.modeler.get_face_by_id(ingoing_face_id)
outgoing_face = hfss.modeler.get_face_by_id(outgoing_face_id)
old_ingoing_center = ingoing_face.center
old_outgoing_center = outgoing_face.center

new_name = f"{design_name}_{frequency}GHz"
refine_name = f"refined_{design_name}_{frequency}GHz"

objects = hfss.modeler.object_names
selection_string = ",".join(objects)
print(f"Objects found: {objects}")

# Create a copy of the geometry in a new design, labeling it appropriately
if not import_from_existing_csv:
    selection_string = ",".join(objects)
    oEditor.Copy([
        "NAME:Selections",
        "Selections:=", selection_string
    ])
    oProject.InsertDesign("HFSS", new_name, "HFSS Terminal Network", "")
    oDesign = oProject.SetActiveDesign(new_name)
    oEditor = oDesign.SetActiveEditor("3D Modeler")
    oEditor.Paste()
    hfss = Hfss(project=project_name, design=new_name, non_graphical=False)
else:    
    selection_string = ",".join(objects)
    oEditor.Copy([
        "NAME:Selections",
        "Selections:=", selection_string
    ])
    oProject.InsertDesign("HFSS", refine_name, "HFSS Terminal Network", "")
    oDesign = oProject.SetActiveDesign(refine_name)
    oEditor = oDesign.SetActiveEditor("3D Modeler")
    oEditor.Paste()
    hfss = Hfss(project=project_name, design=refine_name, non_graphical=False)

tol = 1e-6
ingoing_face_id = None
outgoing_face_id = None

for obj in hfss.modeler.object_names:
    for fid in hfss.modeler.get_object_faces(obj):
        face = hfss.modeler.get_face_by_id(fid)
        cx, cy, cz = face.center
        # Compare ingoing
        if all(abs(a - b) < tol for a, b in zip([cx, cy, cz], old_ingoing_center)):
            ingoing_face_id = fid
        # Compare outgoing
        if all(abs(a - b) < tol for a, b in zip([cx, cy, cz], old_outgoing_center)):
            outgoing_face_id = fid

if ingoing_face_id is None or outgoing_face_id is None:
    raise RuntimeError("Could not find ingoing or outgoing face in new copied design")

oEditor.SetModelUnits(
    [
        "NAME:Units Parameter",
        "Units:=", "mm",
        "Rescale:="	, True,
        "Max Model Extent:="	, 10000
    ])
# Get HFSS modules that we will use throughout
oModuleSolution = oDesign.GetModule("Solutions")
oModuleAnalysis = oDesign.GetModule("AnalysisSetup")
oModuleParametric = oDesign.GetModule("Optimetrics")
oModuleBoundary = oDesign.GetModule("BoundarySetup")
oModuleFields = oDesign.GetModule("FieldsReporter")
oModuleRad = oDesign.GetModule("RadField")

# Clear the simulation of existing radiation boundaries and plane wave excitations
oModuleBoundary.DeleteAllBoundaries()
oModuleBoundary.DeleteAllExcitations()
#%%

def initialize_variables(Ei):
    """
    Initialize project variables for HFSS simulation.

    This function creates the polarization Ephi and sets the incident electric 
    field amplitude Ei to the specified value without sweeping it.

    Parameters:
        Ei (float): The amplitude of the incident electric field to set in the simulation [V/m]
    """
    hfss["Ephi"] = 0
    hfss.variable_manager.set_variable("Ei", str(Ei), sweep=False)

def get_faces_from_face_id(ingoing_face_id, outgoing_face_id):
    """
    Retrieve ingoing and outgoing faces by their IDs and assign radiation boundaries.

    This function:
    - Retrieves the face corresponding to the outgoing_face_id and stores it in a named face list "outgoing".
    - Retrieves the ingoing face by ingoing_face_id.
    - Assigns radiation boundaries to both the ingoing and outgoing faces, allowing them to behave as open-space/vacuum boundaries.

    Parameters:
        ingoing_face_id (int): The ID of the face where the incident plane wave enters the structure.
        outgoing_face_id (int): The ID of the face from which waves exit the structure.

    Returns:
        tuple: A tuple (plane_wave_face, outgoing_face), where each element is a Face object from the HFSS model.
    """
    try:
        outgoing_face = hfss.modeler.get_face_by_id(outgoing_face_id)
        hfss.modeler.create_face_list([outgoing_face], name="outgoing")
    except Exception as e:
        raise RuntimeError(f"Error retrieving 'outgoing' facelist: {e}")

    try:
        plane_wave_face = hfss.modeler.get_face_by_id(ingoing_face_id)
    except Exception as e:
        raise RuntimeError(f"Error retrieving 'ingoing' facelist: {e}")
    hfss.assign_radiation_boundary_to_faces(assignment=[plane_wave_face.id], name="rbin")
    hfss.assign_radiation_boundary_to_faces(assignment=[outgoing_face.id], name="rbout")
    
    return plane_wave_face, outgoing_face

def create_local_coordinate_system(outgoing_face, outgoing_face_cs_x, outgoing_face_cs_y):
    """
    Create a local coordinate system aligned with the outgoing face.

    This function creates a new coordinate system named "outgoing_cs" at the center 
    of the specified outgoing face. The new coordinate system is defined in "axis" 
    mode, with custom x and y directions provided.

    Parameters:
        outgoing_face (Face): The face object where the local coordinate system will be centered.
        outgoing_face_cs_x (list or tuple): A 3-element vector specifying the x-direction.
        outgoing_face_cs_y (list or tuple): A 3-element vector specifying the y-direction.
    """
    hfss.modeler.create_coordinate_system(
        origin=outgoing_face.center,
        reference_cs="Global",
        name="outgoing_cs",
        mode="axis",
        x_pointing=outgoing_face_cs_x,
        y_pointing=outgoing_face_cs_y
    )

# Add the expression for outgoing power to the fields calculator in HFSS
def add_outgoing_power_to_calculator():
    """
    Add an expression to the HFSS fields calculator to compute outgoing power.

    This function constructs and stores a named calculator expression "outgoing_power" 
    that computes the magnitude of the real part of the Poynting vector integrated 
    over the surface named "outgoing". This corresponds to the total power exiting 
    through the output face of the waveguide.

    Notes:
        - In some cases (e.g., during early iterations or with certain boundary conditions), 
          HFSS may report outgoing power greater than incoming power. 
        - For stable simulations, the outgoing power should eventually converge to the correct value.
    """
    oModuleFields.ClearAllNamedExpr()
    oModuleFields.CalcStack("Clear")
    oModuleFields.CopyNamedExprToStack("Vector_RealPoynting")
    oModuleFields.CalcOp("Mag")
    oModuleFields.EnterSurf("outgoing")
    oModuleFields.CalcOp("Integrate")
    oModuleFields.AddNamedExpression("outgoing_power", "Fields")

def get_incoming_phi_theta_num(theta_lower, theta_upper, phi_lower, phi_upper, theta_step, phi_step):
    """
    Calculate the number of points in the sweep for incoming theta and phi angles,
    based on their ranges and step sizes.

    Parameters:
        theta_lower (float): Lower bound of theta angle (degrees).
        theta_upper (float): Upper bound of theta angle (degrees).
        phi_lower (float): Lower bound of phi angle (degrees).
        phi_upper (float): Upper bound of phi angle (degrees).
        theta_step (float): Step size for theta angle (degrees).
        phi_step (float): Step size for phi angle (degrees).

    Returns:
        tuple: (theta_num, phi_num), number of points in theta and phi directions.
    """
    theta_range = theta_upper - theta_lower
    phi_range = phi_upper - phi_lower
    if theta_range == 0:
        theta_num = 1
    else:
        theta_num = int(np.ceil(theta_range / theta_step)) + 1
    if phi_range == 0:
        phi_num = 1
    else:
        phi_num = int(np.ceil(phi_range / phi_step)) + 1
    return theta_num, phi_num

def get_radiation_phi_theta_num(theta_lower, theta_upper, phi_lower, phi_upper, a, b, freq, fineness, maximum_coarseness,
                                minimum_coarseness):
    """
    Determine angular resolution and number of points to sweep in theta and phi for radiation pattern.

    This function calculates the angular step sizes and total number of angular points 
    in both theta and phi directions based on the specified waveguide dimensions, frequency, 
    and desired sampling fineness. The step sizes are bounded between specified minimum 
    and maximum angular coarseness.

    Parameters:
        theta_lower (float): Lower bound of theta sweep (in degrees).
        theta_upper (float): Upper bound of theta sweep (in degrees).
        phi_lower (float): Lower bound of phi sweep (in degrees).
        phi_upper (float): Upper bound of phi sweep (in degrees).
        a (float): Width of the waveguide (in mm), used to scale phi resolution.
        b (float): Height of the waveguide (in mm), used to scale theta resolution.
        freq (float): Frequency in GHz.
        fineness (float): A scaling factor controlling the desired angular resolution.
        maximum_coarseness (float): Upper bound on angular step size (in degrees).
        minimum_coarseness (float): Lower bound on angular step size (in degrees).

    Returns:
        tuple:
            - theta_step (float): Computed angular step size in theta (degrees).
            - phi_step (float): Computed angular step size in phi (degrees).
            - theta_num (int): Number of theta points in the sweep.
            - phi_num (int): Number of phi points in the sweep.
    """
    theta_range = theta_upper - theta_lower
    phi_range = phi_upper - phi_lower

    # Find wavelength and steps (using appropriate unit conversions)
    wavelength = c / (freq * 1e9)

    # Compute angular resolution for theta
    theta_step = np.minimum(maximum_coarseness, theta_range / np.pi * wavelength / (fineness * (b * 1e-3)))
    theta_step = np.maximum(minimum_coarseness, theta_step)

    # Compute angular resolution for phi
    phi_step = np.minimum(maximum_coarseness, phi_range / np.pi * wavelength / (fineness * (a * 1e-3)))
    phi_step = np.maximum(minimum_coarseness, phi_step)

    # Compute number of steps
    theta_num = int(np.ceil(theta_range / theta_step)) + 1
    phi_num = int(np.ceil(phi_range / phi_step)) + 1

    return theta_step, phi_step, theta_num, phi_num


# Helper methods for naming conventions

def get_plane_wave_name(frequency):
    """
    Generate a name for a plane wave excitation at a given frequency.

    Parameters:
        frequency (float): Frequency in GHz.

    Returns:
        str: Name of the plane wave excitation (e.g., "plane_wave_150GHz").
    """
    return f"plane_wave_{frequency}GHz"

def get_radiation_sphere_name(frequency):
    """
    Generate a name for the radiation sphere at a given frequency.

    Parameters:
        frequency (float): Frequency in GHz.

    Returns:
        str: Name of the radiation sphere (e.g., "radiation_sphere_150GHz").
    """
    return f"radiation_sphere_{frequency}GHz"

def get_setup_name(frequency):
    """
    Generate a standard simulation setup name for a given frequency.

    Parameters:
        frequency (float): Frequency in GHz.

    Returns:
        str: Setup name (e.g., "150GHz").
    """
    return f"{frequency}GHz"

def get_parametric_setup_name(frequency):
    """
    Generate a name for a parametric setup (e.g., Ephi sweep) at a given frequency.

    Parameters:
        frequency (float): Frequency in GHz.

    Returns:
        str: Name of the parametric setup (e.g., "E_phi_sweep_150GHz").
    """
    return f"E_phi_sweep_{frequency}GHz"

def setup_radiation(rad_params, frequency):
    """
    Set up radiation boundaries and insert a far-field infinite sphere for HFSS simulation.

    This function:
    - Computes angular step sizes and number of points for theta and phi based on the given
      waveguide geometry and frequency.
    - Constructs a far-field infinite radiation sphere using the computed resolution.
    - Associates the radiation sphere with the "outgoing" face and the local coordinate system "outgoing_cs".

    Parameters:
        rad_params (tuple): A tuple containing the following values:
            - rad_theta_lower (float): Lower bound of theta (degrees)
            - rad_theta_upper (float): Upper bound of theta (degrees)
            - rad_phi_lower (float): Lower bound of phi (degrees)
            - rad_phi_upper (float): Upper bound of phi (degrees)
            - a (float): Waveguide width (mm)
            - b (float): Waveguide height (mm)
            - fineness (float): Desired resolution control factor
            - maximum_coarseness (float): Maximum allowable angular step (degrees)
            - minimum_coarseness (float): Minimum allowable angular step (degrees)
        frequency (float): Frequency in GHz at which the radiation sphere is defined.

    Returns:
        None
    """
    rad_theta_lower, rad_theta_upper, rad_phi_lower, rad_phi_upper, a, b, fineness, maximum_coarseness, minimum_coarseness = rad_params

    # Compute angular resolution and number of points
    rad_theta_step, rad_phi_step, rad_theta_num, rad_phi_num = get_radiation_phi_theta_num(
        rad_theta_lower, rad_theta_upper,
        rad_phi_lower, rad_phi_upper,
        a, b, frequency,
        fineness, maximum_coarseness, minimum_coarseness
    )

    sphere_name = get_radiation_sphere_name(frequency)

    # Insert infinite radiation sphere for far-field observation
    hfss.insert_infinite_sphere(
        definition="Theta-Phi",
        x_start=rad_theta_lower, x_stop=rad_theta_upper, x_step=rad_theta_step,
        y_start=rad_phi_lower, y_stop=rad_phi_upper, y_step=rad_phi_step,
        custom_radiation_faces="outgoing",
        custom_coordinate_system="outgoing_cs",
        name=sphere_name
    )
    
# Setup plane wave
def setup_plane_wave(plane_wave_face, i_theta_lower, i_theta_upper, i_phi_lower, i_phi_upper,
                     frequency, i_theta_step, i_phi_step, refine=False):

    """
    Set up an incident plane wave excitation in HFSS using spherical angle sweeps.

    This function:
    - Computes the number of theta and phi points based on the given step sizes.
    - Constructs a plane wave source on the specified input face using spherical coordinates.
    - Uses symbolic polarization components: `Ephi` for phi-polarized and `1 - Ephi` for theta-polarized.
    - Optionally adds a "_refine" suffix to the excitation name.

    Note:
        HFSS does not allow the plane wave excitation to be defined in a relative coordinate system; it must use the Global CS.
        Polarization vectors:
            - θ̂ = (cosθcosϕ, cosθsinϕ, -sinθ)
            - ϕ̂ = (-sinϕ, cosϕ, 0)

    Parameters:
        plane_wave_face (Face): The face where the plane wave excitation is applied.
        i_theta_lower (float): Lower bound of theta sweep (in degrees).
        i_theta_upper (float): Upper bound of theta sweep (in degrees).
        i_phi_lower (float): Lower bound of phi sweep (in degrees).
        i_phi_upper (float): Upper bound of phi sweep (in degrees).
        frequency (float): Frequency in GHz for naming and setup.
        i_theta_step (float): Angular step size for theta sweep (in degrees).
        i_phi_step (float): Angular step size for phi sweep (in degrees).
        refine (bool, optional): If True, appends "_refine" to the plane wave excitation name.

    Returns:
        None
    """
    
    # Compute number of angular steps
    i_theta_num, i_phi_num = get_incoming_phi_theta_num(i_theta_lower, i_theta_upper, i_phi_lower,
                                                                                  i_phi_upper, i_theta_step, i_phi_step)
    plane_wave_name = get_plane_wave_name(frequency)
    if refine:
        plane_wave_name += "_refine"
    
    # Define plane wave excitation using symbolic polarization
    hfss.plane_wave(
        assignment=[plane_wave_face],
        vector_format="Spherical",
        origin=plane_wave_face.center,
        polarization=["Ephi", "1-Ephi"],
        propagation_vector=[
            [f"{i_phi_lower}deg", f"{i_phi_upper}deg", i_phi_num],
            [f"{i_theta_lower}deg", f"{i_theta_upper}deg", i_theta_num]
        ],
        name=plane_wave_name
    )

def create_setup(frequency, max_delta_E, max_passes):
    """
    Create an HFSS simulation setup for a given frequency

    This function:
    - Creates a new solution setup for the specified frequency.
      from that frequency's setup to speed up convergence.
    - Adds a parametric sweep over the `Ephi` variable from 0 to 1 in steps of 1 
      (to sweep between phi and theta polarization).

    Parameters:
        frequency (float): Frequency in GHz for this setup.
        max_delta_E (float): Maximum allowed change in electric field between adaptive passes (convergence criteria).
        max_passes (int): Maximum number of adaptive passes allowed for convergence.

    Returns:
        None

    Notes:
        - The parametric sweep is stored under a name like "E_phi_sweep_150GHz".
    """
    freq_str = f"{frequency}GHz"
    setup_name = get_setup_name(frequency)

    oModuleAnalysis.InsertSetup("HfssDriven", [
        f"NAME:{setup_name}",
        "Frequency:=", freq_str,
        "MaxDeltaE:=", max_delta_E,
        "MaximumPasses:=", max_passes
    ])

    # Insert a parametric sweep over the two electric field polarizations
    parametric_setup_name = get_parametric_setup_name(frequency)
    oModuleParametric.InsertSetup("OptiParametric",
                                  [
                                      f"NAME:{parametric_setup_name}",
                                      "IsEnabled:=", True,
                                      [
                                          "NAME:ProdOptiSetupDataV2",
                                          "SaveFields:=", True,
                                      ],
                                      "Sim. Setups:=", [setup_name],
                                      [
                                          "NAME:Sweeps",
                                          [
                                              "NAME:SweepDefinition",
                                              "Variable:=", "Ephi",
                                              "Data:=", "LIN 0 1 1",
                                          ]
                                      ],
                                  ])

def run_analysis(num_cores, max_delta_E, max_passes, plane_wave_face, Ei, output_file_location, i_theta_lower, i_theta_upper, i_phi_lower, i_phi_upper,
                 frequency, i_theta_step, i_phi_step, rad_params, sweep = "adaptive", max_difference = 0.02, import_from_existing_csv = False, E_phi = None,
                 waveguide_data_csv = None, far_field_data_csv = None):
    """
    Run a full electromagnetic simulation in HFSS, including plane wave setup, field analysis, 
    adaptive or discrete refinement sweeps, and optional data import/export.

    This function supports multiple sweep modes:
      - "discrete": Runs a single setup and extracts results without refinement.
      - "adaptive": Iteratively refines regions based on difference threshold until convergence.
      - "zoom": Placeholder mode, intended for future implementation.
    
    It also allows importing previously computed waveguide and far-field data to bypass the initial simulation.

    Parameters:
        num_cores (int): Number of processor cores to allocate for the analysis.
        max_delta_E (float): Maximum acceptable change in electric field between adaptive passes.
        max_passes (int): Maximum number of adaptive passes allowed.
        plane_wave_face (object): HFSS face object where the plane wave is applied.
        Ei (float): Incident electric field magnitude.
        output_file_location (str): Path to the directory where result CSV files should be saved.
        i_theta_lower (float): Lower bound of the incident theta angle (degrees).
        i_theta_upper (float): Upper bound of the incident theta angle (degrees).
        i_phi_lower (float): Lower bound of the incident phi angle (degrees).
        i_phi_upper (float): Upper bound of the incident phi angle (degrees).
        frequency (float): Frequency of simulation in GHz.
        i_theta_step (float): Step size in theta for incident plane wave sweep.
        i_phi_step (float): Step size in phi for incident plane wave sweep.
        rad_params (tuple): Tuple containing radiation configuration parameters:
            (rad_theta_lower, rad_theta_upper, rad_phi_lower, rad_phi_upper, a, b, fineness, maximum_coarseness, minimum_coarseness)
        sweep (str, optional): Type of sweep ("adaptive", "discrete", "zoom"). Default is "adaptive".
        max_difference (float, optional): Tolerance for convergence when refining fields. Default is 0.02.
        import_from_existing_csv (bool, optional): Whether to skip simulation and import existing waveguide/far-field data. Default is False.
        E_phi (float, optional): Optional specific polarization value (0 or 1) to sweep only one polarization. Default is both.
        waveguide_data_csv (str, optional): Path to CSV file containing precomputed waveguide data (used only if importing).
        far_field_data_csv (str, optional): Path to CSV file containing precomputed far-field data (used only if importing).

    Behavior:
        - Creates and applies a plane wave excitation.
        - Runs HFSS analysis on total field configuration.
        - Extracts waveguide and far-field data for E_phi = 0 and 1, or a specified value.
        - If `adaptive` sweep is chosen, performs iterative refinement until convergence criteria is met.
        - Stores final refined data in CSV files for each polarization and frequency.
        - If importing CSVs, skips simulation and loads precomputed results.

    Returns:
        None

    Notes:
        - The “zoom” sweep option is currently a placeholder for future development.
    """
    
    if sweep == "discrete":
        print(f"Discrete sweep has begun for frequency {frequency}GHz")
    elif sweep == "adaptive":
        print(f"Adaptive sweep has begun for frequency {frequency}GHz")
    elif sweep == "zoom":
        print(f"Zoom sweep has begun for frequency {frequency}GHz")

    # Declare waveguide_data and far_field_data variables
    waveguide_data_0 = waveguide_data_1 = None
    far_field_data_0 = far_field_data_1 = None
    waveguide_data_from_csv = far_field_data_from_csv = None

    if not import_from_existing_csv:
        print("Proceeding with default setup")
        setup_plane_wave(
            plane_wave_face, i_theta_lower, i_theta_upper, i_phi_lower, i_phi_upper,
            frequency, i_theta_step, i_phi_step
        )

        # Change from default scattered fields to total fields formulation
        oModuleSolution.EditSources(["FieldType:=", "TotalFields"])

        hfss.analyze(cores=num_cores)

        rad_theta_lower, rad_theta_upper, rad_phi_lower, rad_phi_upper, a, b, fineness, maximum_coarseness, minimum_coarseness = rad_params
        rad_theta_step, rad_phi_step, rad_theta_num, rad_phi_num = get_radiation_phi_theta_num(rad_theta_lower, rad_theta_upper,
            rad_phi_lower, rad_phi_upper, a, b, frequency, fineness, maximum_coarseness, minimum_coarseness)

        # Extract initial waveguide and far field data (final for discrete sweep, to be adjusted for adaptive sweep)
        waveguide_data_0 = extract_waveguide_data(frequency, Ei, plane_wave_face, i_theta_lower, i_theta_upper,
            i_phi_lower, i_phi_upper, i_theta_step, i_phi_step, output_file_location, 0, csv=True, refine=False)
        far_field_data_0 = extract_far_field_data(frequency, i_theta_lower, i_theta_upper, i_phi_lower, i_phi_upper,
            i_theta_step, i_phi_step, rad_theta_step, rad_phi_step, output_file_location, 0, csv=True, refine=False)
        waveguide_data_1 = extract_waveguide_data(frequency, Ei, plane_wave_face, i_theta_lower, i_theta_upper,
            i_phi_lower, i_phi_upper, i_theta_step, i_phi_step, output_file_location, 1, csv=True, refine=False)
        far_field_data_1 = extract_far_field_data(frequency, i_theta_lower, i_theta_upper, i_phi_lower, i_phi_upper,
            i_theta_step, i_phi_step, rad_theta_step, rad_phi_step, output_file_location, 1, csv=True, refine=False)
    else:
        print("Reading in waveguide data and far-field data from existing csv files")
        waveguide_data_from_csv = pd.read_csv(waveguide_data_csv)
        far_field_data_from_csv = pd.read_csv(far_field_data_csv)

    if sweep == "adaptive" or sweep == "zoom":
        # Disable analysis and parametric setup before adaptive sweep begins, and insert new refinement setup.
        parametric_setup_name = get_parametric_setup_name(frequency)
        try:
            oModuleParametric.EnableSetup(parametric_setup_name, False)
        except Exception:
            pass
        setup_name = get_setup_name(frequency)
        try:
            oModuleAnalysis.EditSetup(setup_name,
                                  [
                                      f"NAME:{setup_name}",
                                      "IsEnabled:=", False,
                                  ])
        except Exception:
            pass

        Ephi_values = [E_phi] if E_phi is not None else [0, 1]
        # Sweep over the polarizations
        for Ephi in Ephi_values:
            # Set the Ephi variable to the desired polarization
            hfss["Ephi"] = Ephi
            print(f"Beginning refinement process for frequency {frequency}GHz, Ephi={Ephi}")

            # Select correct dataset for the polarization
            if not import_from_existing_csv:
                waveguide_data = waveguide_data_0 if Ephi == 0 else waveguide_data_1
                far_field_data = far_field_data_0 if Ephi == 0 else far_field_data_1
            else:
                waveguide_data = waveguide_data_from_csv
                far_field_data = far_field_data_from_csv

            if sweep == "adaptive":
                # Repeat until converge: (1) find a new region to refine (2) create a plane wave to refine this region
                # and analyze
                incoming_power = get_incoming_power(Ei, plane_wave_face)
                converged = False
                while not converged:
                    result = find_regions_to_refine(incoming_power, waveguide_data, max_difference)
                    if result is not False:
                        refine_regions = result
                        # Update the waveguide and far-field data frames every cycle of refinement (refining all regions). You can change the 
                        # csv=False to csv=True if you want the .csv file to be updated every cycle
                        waveguide_data, far_field_data = run_refined_plane_wave(plane_wave_face, frequency, num_cores, max_delta_E, max_passes, Ei, output_file_location,
                        Ephi, waveguide_data, far_field_data, rad_params, refine_regions, csv=False)
                    else:
                        print("Adaptive sweep has converged")
                        converged = True

                # Output two csv files at the end of the sweep for each polarization for each frequency
                print(f"Exporting waveguide and far field data to csv for frequency {frequency}GHz, Ephi={Ephi}")
                folder_name = f"{project_name}_{design_name}_{frequency}GHz_Ephi={Ephi}"
                folder_path = os.path.join(output_file_location, folder_name)
                os.makedirs(folder_path, exist_ok=True)

                waveguide_output_filename = "refined_waveguide.csv"
                full_path = os.path.join(folder_path, waveguide_output_filename)
                waveguide_data.to_csv(full_path, index=False)
                
                far_field_output_filename = "refined_far_field.csv"
                full_path = os.path.join(folder_path, far_field_output_filename)
                far_field_data.to_csv(full_path, index=False)
                
                print(f"Waveguide and far field data to csv for frequency {frequency}GHz, Ephi={Ephi} exported")
            if sweep == "zoom":
                print("Temporary filler")

def run_refined_plane_wave(plane_wave_face, frequency, num_cores, max_delta_E, max_passes, Ei, output_file_location, Ephi, waveguide_data, far_field_data, rad_params,
    refine_regions, csv=False):
    """
    Perform refined plane wave sweeps over specified angular regions to improve simulation resolution,
    updating waveguide and far-field data accordingly.

    This function:
    - Iterates over a list of refinement regions, each defining angular bounds and step sizes.
    - For each region, deletes existing excitations and inserts a new adaptive setup with refined plane wave sweeps.
    - Extracts waveguide and far-field data for each refined region.
    - Merges new data with existing dataframes, optionally saving merged data to CSV files.
    - Handles user interruption gracefully, allowing partial data export.

    Parameters:
        plane_wave_face (Face): HFSS face object where plane wave excitation is applied.
        frequency (float): Simulation frequency in GHz.
        num_cores (int): Number of CPU cores to use for HFSS analysis.
        max_delta_E (float): Convergence criterion for maximum allowed electric field change.
        max_passes (int): Maximum number of adaptive passes per setup.
        Ei (float): Incident electric field magnitude.
        output_file_location (str): Directory path to save output CSV files.
        Ephi (float): Polarization parameter value for the sweep (0 or 1).
        waveguide_data (pd.DataFrame): Existing waveguide data to update.
        far_field_data (pd.DataFrame): Existing far-field data to update.
        rad_params (tuple): Radiation parameters tuple:
            (rad_theta_lower, rad_theta_upper, rad_phi_lower, rad_phi_upper, a, b, fineness, maximum_coarseness, minimum_coarseness)
        refine_regions (list of dict): List of regions to refine, each dict contains:
            - "phi_range": (phi_lower, phi_upper) tuple
            - "theta_range": (theta_lower, theta_upper) tuple
            - "phi_step": angular step size in phi (degrees)
            - "theta_step": angular step size in theta (degrees)
        csv (bool, optional): If True, merged dataframes are saved to CSV files after refinement. Default is False.

    Returns:
        tuple:
            - waveguide_data (pd.DataFrame): Updated waveguide data including refined regions.
            - far_field_data (pd.DataFrame): Updated far-field data including refined regions.

    Notes:
        - Uses a nested helper function `merge_fast` to efficiently combine existing and new data.
        - Deletes adaptive setups after each region refinement to keep the project clean.
        - Handles KeyboardInterrupt allowing the user to save partial results before exiting.
    """
    rad_theta_lower, rad_theta_upper, rad_phi_lower, rad_phi_upper, a, b, fineness, maximum_coarseness, minimum_coarseness = rad_params
    rad_theta_step, rad_phi_step, rad_theta_num, rad_phi_num = get_radiation_phi_theta_num(rad_theta_lower, rad_theta_upper,
        rad_phi_lower, rad_phi_upper, a, b, frequency, fineness, maximum_coarseness, minimum_coarseness)
    key_cols = ["IWavePhi", "IWaveTheta"]

    def merge_fast(existing_df, new_dfs, name, merge_csv=False):
        """
    Efficiently merge new chunks of data with an existing DataFrame, avoiding duplicates based on key columns.

    This function:
    - Concatenates a list of new DataFrame chunks into a single DataFrame.
    - Removes rows from the existing DataFrame that have keys overlapping with the new data.
    - Combines the filtered existing data with the new data to produce an updated DataFrame.
    - Optionally saves the combined DataFrame to a CSV file in a structured folder path.

    Parameters:
        existing_df (pd.DataFrame): The original DataFrame to be updated.
        new_dfs (list of pd.DataFrame): List of new DataFrame chunks to merge into the existing DataFrame.
        name (str): Identifier name for the dataset being merged (e.g., "waveguide" or "far_field").
        merge_csv (bool, optional): If True, saves the merged DataFrame to a CSV file after merging. Default is False.

    Returns:
        pd.DataFrame: The combined DataFrame containing existing and new data without duplicate keys.

    Notes:
        - Assumes that the key columns used for identifying duplicates are ["IWavePhi", "IWaveTheta"].
        - When saving to CSV, the file is named "refined_{name}.csv" and saved under a folder structure:
          {output_file_location}/{project_name}_{design_name}_{frequency}GHz_Ephi={Ephi}/
        - Prints timing and progress information during the merge process.
    """
        print(f"[{name}] Starting merge of {len(new_dfs)} chunks...")
        t0 = time.perf_counter()

        new_df = pd.concat(new_dfs, ignore_index=True)
        existing_idx = existing_df.set_index(key_cols, drop=False)
        new_keys = new_df.set_index(key_cols, drop=False).index.unique()
        filtered_existing = existing_idx.drop(new_keys, errors="ignore").reset_index(drop=True)

        combined = pd.concat([filtered_existing, new_df], ignore_index=True)

        if merge_csv:
            print("Updating CSV file following this set of refining")
            folder_name = f"{project_name}_{design_name}_{frequency}GHz_Ephi={Ephi}"
            folder_path = os.path.join(output_file_location, folder_name)
            os.makedirs(folder_path, exist_ok=True)

            output_filename = f"refined_{name}.csv"
            full_path = os.path.join(folder_path, output_filename)
            combined.to_csv(full_path, index=False)
            
        elapsed = time.perf_counter() - t0
        print(f"[{name}] Merge took {elapsed:.3f}s")
        return combined

    print(f"Running {len(refine_regions)} refined plane wave regions at {frequency} GHz")

    waveguide_chunks = []
    far_field_chunks = []

    try:
        for idx, region in enumerate(refine_regions):
            freq_str = f"{frequency}GHz"

            print(f"\n--- Refining region {idx + 1} of {len(refine_regions)} ---")
            plane_phi_lower, plane_phi_upper = region["phi_range"]
            plane_theta_lower, plane_theta_upper = region["theta_range"]
            plane_phi_step = region["phi_step"]
            plane_theta_step = region["theta_step"]

            oModuleBoundary.DeleteAllExcitations()
            # Insert a new setup for refinement
            setup_name = get_setup_name(frequency) + "_refine"
            oModuleAnalysis.InsertSetup("HfssDriven", [
                f"NAME:{setup_name}",
                "Frequency:=", freq_str,
                "MaxDeltaE:=", max_delta_E,
                "MaximumPasses:=", max_passes,
            ])

            setup_plane_wave(
                plane_wave_face,
                plane_theta_lower, plane_theta_upper,
                plane_phi_lower, plane_phi_upper,
                frequency, plane_theta_step, plane_phi_step,
                refine=True
            )

            # Change from default scattered fields formulation to total fields formulation
            oModuleSolution.EditSources(["FieldType:=", "TotalFields"])
            hfss.analyze(cores=num_cores)

            waveguide_chunks.append(
                extract_waveguide_data(
                    frequency, Ei, plane_wave_face,
                    plane_theta_lower, plane_theta_upper,
                    plane_phi_lower, plane_phi_upper,
                    plane_theta_step, plane_phi_step,
                    output_file_location, Ephi, csv=False, refine=True
                )
            )

            far_field_chunks.append(
                extract_far_field_data(
                    frequency, plane_theta_lower, plane_theta_upper,
                    plane_phi_lower, plane_phi_upper,
                    plane_theta_step, plane_phi_step, rad_theta_step, rad_phi_step,
                    output_file_location, Ephi, csv=False, refine=True
                )
            )
            # Delete the analysis setup after refining each region
            oModuleAnalysis.DeleteSetups([setup_name])


    # Add safe exits for keyboard interruptions (just in case the user desires to export the CSV early)
    except KeyboardInterrupt:
        print(f"\n[!] User interrupted refining process")

        user_input = input(
            "Type 'csv' to save partial results to file, or press Enter to skip saving: ").strip().lower()
        if user_input == "csv":
            print("Saving last updated results to a CSV")
            csv = True
        else:
            csv = False

        waveguide_data = merge_fast(waveguide_data, waveguide_chunks, "waveguide", merge_csv=csv)
        far_field_data = merge_fast(far_field_data, far_field_chunks, "far_field", merge_csv=csv)
        if csv:
            print("Latest csv files have been successfully saved")
        sys.exit(1)

    # If no interruption, continue here
    waveguide_data = merge_fast(waveguide_data, waveguide_chunks, "waveguide", merge_csv=csv)
    far_field_data = merge_fast(far_field_data, far_field_chunks, "far_field", merge_csv=csv)
    return waveguide_data, far_field_data

def find_regions_to_refine(incoming_power, waveguide_data, max_difference):
    """
    Analyze outgoing power data on a theta-phi grid to identify angular regions requiring finer sampling.

    This function compares adjacent points in the 2D outgoing power distribution (theta vs phi)
    and detects intervals where the fractional difference in power exceeds a given threshold.
    It generates refinement intervals along the phi and theta axes, merges contiguous intervals,
    and returns a list of ordered regions prioritized by power contrast.

    Parameters:
        incoming_power (float): Reference power (usually the incident power) used for normalization.
        waveguide_data (pd.DataFrame): DataFrame containing waveguide measurement points with columns:
            - "IWaveTheta": Incident theta angles (degrees)
            - "IWavePhi": Incident phi angles (degrees)
            - "OutgoingPower": Measured outgoing power at each (theta, phi) point
        max_difference (float): Fractional difference threshold above which a region is marked for refinement.

    Returns:
        list of dict or bool: 
            - If regions needing refinement are found, returns a list of dictionaries, each describing:
                - "phi_range": Tuple (phi_lower, phi_upper) in degrees defining phi bounds
                - "theta_range": Tuple (theta_lower, theta_upper) in degrees defining theta bounds
                - "phi_step": Suggested phi step size for refinement (degrees)
                - "theta_step": Suggested theta step size for refinement (degrees)
                - "frac_diff": The fractional power difference motivating refinement
                - "pair": Tuple of index pairs ((i,j), (i2,j2)) indicating data points defining the interval
            - Returns False if no regions exceed the threshold (no refinement needed).

    Behavior:
        - Removes duplicate angle entries.
        - Reshapes data into a 2D grid indexed by theta and phi.
        - Compares each point only with its immediate next valid neighbor in phi and theta directions.
        - Generates finer subdivision steps to split intervals into 10 parts for refinement.
        - Merges contiguous intervals that share boundaries and step sizes to reduce redundancy.
        - Sorts the output list by descending fractional power difference (priority for refinement).

    Notes:
        - Prints informative messages about found regions and merging progress.
        - Uses inner helper function `merge_contiguous_regions_by_group` to group and merge intervals.

    Example output item:
        {
            "phi_range": (10.0, 15.0),
            "theta_range": (20.0, 20.0),
            "phi_step": 0.5,
            "theta_step": 0.0,
            "frac_diff": 0.12,
            "pair": ((i1, j1), (i2, j2))
        }
    """
    max_frac_diff_seen = 0.0

    print("Finding regions to refine!")

    # Focus only on unique angular combinations
    key_cols = ["IWaveTheta", "IWavePhi"]
    unique_points = waveguide_data.drop_duplicates(subset=key_cols)

    # Pivot data to form a 2D power grid: rows = theta, columns = phi
    pivot = unique_points.pivot(index="IWaveTheta", columns="IWavePhi", values="OutgoingPower")

    theta_vals = pivot.index.values
    phi_vals = pivot.columns.values
    power_grid = pivot.values
    max_power = incoming_power

    print(f"Theta values ({len(theta_vals)}): {theta_vals}")
    print(f"Phi values ({len(phi_vals)}): {phi_vals}")
    print(f"OutgoingPower max: {max_power}")

    refine_regions = []

    # Loop over each point in the grid to find high-contrast neighbors
    for i in range(power_grid.shape[0]):
        for j in range(power_grid.shape[1]):
            current_power = power_grid[i, j]
            if np.isnan(current_power):
                continue

            # Check neighbor in +φ direction
            for jj in range(j + 1, power_grid.shape[1]):
                neighbor_power = power_grid[i, jj]
                if not np.isnan(neighbor_power):
                    frac_diff = abs(current_power - neighbor_power) / max_power
                    max_frac_diff_seen = max(max_frac_diff_seen, frac_diff)

                    if frac_diff > max_difference:
                        phi1, phi2 = phi_vals[j], phi_vals[jj]
                        theta1 = theta_vals[i]

                        # Define φ-refinement region at fixed θ
                        plane_phi_lower = phi1
                        plane_phi_upper = phi2
                        plane_theta_lower = theta1
                        plane_theta_upper = theta1

                        # Break interval up until 10 distinct pieces (10 times finer)
                        plane_phi_num = 11
                        plane_phi_step = (plane_phi_upper - plane_phi_lower) / (plane_phi_num - 1) if plane_phi_num > 1 else 0.0
                        plane_theta_step = 0.0

                        refine_regions.append({
                            "phi_range": (plane_phi_lower, plane_phi_upper),
                            "theta_range": (plane_theta_lower, plane_theta_upper),
                            "phi_step": plane_phi_step,
                            "theta_step": plane_theta_step,
                            "frac_diff": frac_diff,
                            "pair": ((i, j), (i, jj))
                        })
                    break  # Only check immediate next valid neighbor

            # Check neighbor in +θ direction
            for ii in range(i + 1, power_grid.shape[0]):
                neighbor_power = power_grid[ii, j]
                if not np.isnan(neighbor_power):
                    frac_diff = abs(current_power - neighbor_power) / max_power
                    max_frac_diff_seen = max(max_frac_diff_seen, frac_diff)

                    if frac_diff > max_difference:
                        theta1, theta2 = theta_vals[i], theta_vals[ii]
                        phi1 = phi_vals[j]

                        # Define θ-refinement region at fixed φ
                        plane_phi_lower = phi1
                        plane_phi_upper = phi1
                        plane_theta_lower = theta1
                        plane_theta_upper = theta2

                        # Break interval up until 10 distinct pieces (10 times finer)
                        plane_theta_num = 11
                        plane_theta_step = (plane_theta_upper - plane_theta_lower) / (plane_theta_num - 1) if plane_theta_num > 1 else 0.0
                        plane_phi_step = 0.0

                        refine_regions.append({
                            "phi_range": (plane_phi_lower, plane_phi_upper),
                            "theta_range": (plane_theta_lower, plane_theta_upper),
                            "phi_step": plane_phi_step,
                            "theta_step": plane_theta_step,
                            "frac_diff": frac_diff,
                            "pair": ((i, j), (ii, j))
                        })
                    break  # Only check immediate next valid neighbor

    if not refine_regions:
        print("No regions found that exceed the fractional difference threshold.")
        print(f"Maximum fractional difference seen: {max_frac_diff_seen:.5f}")
        return False

    # Sort regions by strength of power contrast (descending)
    refine_regions.sort(key=lambda x: x["frac_diff"], reverse=True)

    print(f"Found {len(refine_regions)} regions to refine.")
    for region in refine_regions:
        (phi_low, phi_high) = region["phi_range"]
        (theta_low, theta_high) = region["theta_range"]
        print(f"φ∈[{phi_low:.3f}, {phi_high:.3f}] θ∈[{theta_low:.3f}, {theta_high:.3f}] "
              f"φ_step={region['phi_step']:.3f}, θ_step={region['theta_step']:.3f}, "
              f"Δ={region['frac_diff']:.5f}")

    # Merge regions along each axis
    def merge_contiguous_regions_by_group(regions, orthogonal_range_key, step_key, merge_range_key):
        if not regions:
            print("No regions to merge.")
            return []

        print(f"\nStarting merge on {len(regions)} regions grouped by '{orthogonal_range_key}'")

        # Group by orthogonal axis (e.g. θ constant for φ merge)
        groups = defaultdict(list)
        for r in regions:
            groups[r[orthogonal_range_key]].append(r)

        merged_all = []

        for ortho_val, group in sorted(groups.items()):
            print(f"\nMerging group with {orthogonal_range_key} = {ortho_val} ({len(group)} regions)")

            group_sorted = sorted(group, key=lambda reg: reg[merge_range_key][0])
            current = group_sorted[0].copy()

            for k, next_region in enumerate(group_sorted[1:], start=1):
                curr_start, curr_end = current[merge_range_key]
                next_start, next_end = next_region[merge_range_key]

                adjacent = (curr_end == next_start)
                same_step = (current[step_key] == next_region[step_key])

                if adjacent and same_step:
                    current[merge_range_key] = (curr_start, next_end)
                    current["frac_diff"] = max(current["frac_diff"], next_region["frac_diff"])
                else:
                    merged_all.append(current)
                    current = next_region.copy()

            merged_all.append(current)

        print(f"\nFinished merging. Total merged regions: {len(merged_all)}")
        return merged_all

    # Split regions by sweep direction: φ-sweeps (fixed θ) vs θ-sweeps (fixed φ)
    phi_sweep = [r for r in refine_regions if r["theta_range"][0] == r["theta_range"][1]]
    theta_sweep = [r for r in refine_regions if r["phi_range"][0] == r["phi_range"][1]]

    # Merge contiguous φ-sweep regions at constant θ
    merged_phi = merge_contiguous_regions_by_group(
        phi_sweep, orthogonal_range_key="theta_range",
        step_key="phi_step", merge_range_key="phi_range"
    )

    # Merge contiguous θ-sweep regions at constant φ
    merged_theta = merge_contiguous_regions_by_group(
        theta_sweep, orthogonal_range_key="phi_range",
        step_key="theta_step", merge_range_key="theta_range"
    )

    # Combine all merged regions and sort by importance
    merged = merged_phi + merged_theta
    merged.sort(key=lambda r: r["frac_diff"], reverse=True)

    print(f"\nMerged into {len(merged)} refined regions.")
    for region in merged:
        (phi_low, phi_high) = region["phi_range"]
        (theta_low, theta_high) = region["theta_range"]
        print(f"φ∈[{phi_low:.3f}, {phi_high:.3f}] θ∈[{theta_low:.3f}, {theta_high:.3f}] "
              f"φ_step={region['phi_step']:.3f}, θ_step={region['theta_step']:.3f}, "
              f"Δ={region['frac_diff']:.5f}")

    return merged

def get_incoming_power(Ei, plane_wave_face):
    """
    Calculate the incoming electromagnetic power on a given face due to an incident plane wave.

    The incoming power is computed using the incident electric field magnitude, the area of the specified face,
    and fundamental constants (speed of light and vacuum permeability). The face area is converted from mm² to m².

    Parameters:
        Ei (float): Magnitude of the incident electric field (in V/m).
        plane_wave_face (Face): HFSS face object representing the incoming wave boundary.

    Returns:
        float: Incoming power in watts (W) incident on the specified face.
    """
    return 1e-6 *  hfss.modeler.get_face_area(plane_wave_face.id) * Ei ** 2 / (2 * c * mu_0)

def extract_waveguide_data(frequency, Ei, plane_wave_face, i_theta_lower, i_theta_upper, i_phi_lower, i_phi_upper,
                           i_theta_step, i_phi_step, output_file_location, E_phi, csv=True, refine=False):
    """
    Extracts waveguide output data by sweeping over incident wave angles and electric field polarization,
    generating a Pandas DataFrame containing outgoing power and electric field magnitude at the waveguide exit.

    This function performs two passes:
    1. Computes the outgoing power at the waveguide exit for each angle and polarization.
    2. Exports and reads the electric field data on the output face grid for each angle and polarization.

    The resulting data is concatenated into a DataFrame with columns including frequency, polarization, 
    incident angles, outgoing power, ingoing power, and field magnitudes. Optionally exports the data as CSV.

    Parameters:
        frequency (float): Frequency in GHz at which data is extracted.
        Ei (float): Magnitude of the incident electric field (V/m).
        plane_wave_face (Face): HFSS face object representing the input wave boundary.
        i_theta_lower (float): Lower bound of incident theta angle sweep (degrees).
        i_theta_upper (float): Upper bound of incident theta angle sweep (degrees).
        i_phi_lower (float): Lower bound of incident phi angle sweep (degrees).
        i_phi_upper (float): Upper bound of incident phi angle sweep (degrees).
        i_theta_step (float): Step size for theta angle sweep (degrees).
        i_phi_step (float): Step size for phi angle sweep (degrees).
        output_file_location (str): Directory path to save output files.
        E_phi (float or None): Electric field polarization parameter; if None, default value is used.
        csv (bool, optional): If True, exports the resulting DataFrame to a CSV file. Defaults to True.
        refine (bool, optional): If True, uses refined solution data. Defaults to False.

    Returns:
        pandas.DataFrame: DataFrame containing waveguide output power, field magnitudes, and metadata for all swept angles and polarizations.

    Notes:
        - The function prints timing info for debugging and progress.
        - The waveguide exit field export method can be toggled by the 'manual_field' global flag.
    """
    freq_str = f"{frequency}GHz"
    setup_name = get_setup_name(frequency) + ("_refine : LastAdaptive" if refine else " : LastAdaptive")
    exit_field_path = os.path.join(output_file_location, f"{project_name}_exitfield_{freq_str}.fld")
        
    i_theta_num, i_phi_num = get_incoming_phi_theta_num(
        i_theta_lower, i_theta_upper, i_phi_lower, i_phi_upper, i_theta_step, i_phi_step)
    i_phi_values = np.linspace(i_phi_lower, i_phi_upper, i_phi_num)
    i_theta_values = np.linspace(i_theta_lower, i_theta_upper, i_theta_num)
    total = len(i_phi_values) * len(i_theta_values)

    power_data = {}
    all_data = []
    count = 0
    start_time = time.perf_counter()

    print(f"Extracting waveguide data for frequency {freq_str}" + (f" and Ephi = {E_phi}" if E_phi is not None else ""))

    # First Pass: Outgoing Power
    print("Starting first pass: outgoing power")
    t1 = time.perf_counter()
    oModuleFields.CalcStack("clear")
    oModuleFields.CopyNamedExprToStack("outgoing_power")
    for i_phi in i_phi_values:
        for i_theta in i_theta_values:
            i_phi_str = f"{i_phi}deg"
            i_theta_str = f"{i_theta}deg"
            t_power = time.perf_counter()
            oModuleFields.ClcEval(setup_name,
                ["Ephi:=", E_phi, "Freq:=", freq_str, "IWavePhi:=", i_phi_str, "IWaveTheta:=", i_theta_str], "Fields")
            result = oModuleFields.GetTopEntryValue(setup_name,
                ["Ephi:=", E_phi, "Freq:=", freq_str, "IWavePhi:=", i_phi_str, "IWaveTheta:=", i_theta_str])
            oModuleFields.CalcStack("pop")
            power_data[(E_phi, i_phi, i_theta)] = float(result[0])
            count += 1
            print(f"[w1] [{count}/{total}] φ={i_phi}deg, θ={i_theta}deg, Ephi={E_phi}, {time.perf_counter() - t_power:.3f}s")

    print(f"First pass time: {(time.perf_counter() - t1):.3f}s")

    count = 0
    # Second Pass: Field Export + Read
    print("Starting second pass: waveguide exit field export")
    t2 = time.perf_counter()

    if manual_field:
        oModuleFields.CalcStack("clear")
        oModuleFields.EnterQty("E")
        oModuleFields.CalcOp("Smooth")
        unit = oEditor.GetModelUnits()
        if outgoing_face_boundary is not None:
            range_min = ['{}{}'.format(i, unit) for i in outgoing_face_boundary[0]]
            range_max = ['{}{}'.format(i, unit) for i in outgoing_face_boundary[1]]
        else:
            # Convert the vertices from the GLOBAL CS to the OUTGOING CS
            outgoing_face = hfss.modeler.get_face_by_id(outgoing_face_id)
            outgoing_cs_origin = np.array(outgoing_face.center)
            
            x_axis_new = np.array(outgoing_face_cs_x)  
            y_axis_new = np.array(outgoing_face_cs_y)
            z_axis_new = np.cross(x_axis_new, y_axis_new)
            
            # Transformation matrix: UDCS -> global
            M = np.column_stack((x_axis_new, y_axis_new, z_axis_new))
            Minv = np.linalg.inv(M)  # global -> UDCS
            
            # Get vertex positions in GLOBAL CS
            vertex_ids = oEditor.GetVertexIDsFromFace(outgoing_face_id)
            vertex_positions = [oEditor.GetVertexPosition(i) for i in vertex_ids]
            
            # Convert to floats and transform to UDCS using outgoing_cs_origin as the center
            vertex_positions_udcs = [
                Minv @ (np.array([float(vx), float(vy), float(vz)]) - outgoing_cs_origin)
                for vx, vy, vz in vertex_positions
            ]
            
            # Find bounding box in UDCS
            x_udcs, y_udcs, z_udcs = zip(*vertex_positions_udcs)
            range_min = ['{}{}'.format(i, unit) for i in [min(x_udcs), min(y_udcs), min(z_udcs)]]
            range_max = ['{}{}'.format(i, unit) for i in [max(x_udcs), max(y_udcs), max(z_udcs)]]
    else:
        oModuleFields.CalcStack("clear")
        oModuleFields.EnterQty("E")
        oModuleFields.CalcOp("Smooth")
        oModuleFields.EnterSurf("outgoing")
        oModuleFields.CalcOp("Value") 

    for i_phi in i_phi_values:
        for i_theta in i_theta_values:
            i_phi_str = f"{i_phi}deg"
            i_theta_str = f"{i_theta}deg"

            t_export = time.perf_counter()
        
            if(manual_field):  
                # If the grid is manually specified, it should be specified with respect to the outgoing CS
                oModuleFields.ExportOnGrid(exit_field_path, range_min, range_max, outgoing_face_field_resolution, setup_name, ["Ephi:=", E_phi, "Freq:=", freq_str, "IWavePhi:=", i_phi_str,
    		    "IWaveTheta:=", i_theta_str, "Phase:=", "0deg"], ["NAME:ExportOption", "IncludePtInOutput:=", True, "RefCSName:=", "outgoing_cs", "PtInSI:=", True, "FieldInRefCS:=", False], "Cartesian", 
            ["0mm","0mm","0mm"], False)
            else:
                # If no grid is specified, HFSS automatically writes with respect to the global coordinate system
                oModuleFields.CalculatorWrite(exit_field_path,
                    ["Solution:=", setup_name],
                    ["Ephi:=", E_phi, "Freq:=", freq_str, "IWavePhi:=", i_phi_str, "IWaveTheta:=", i_theta_str])
                print(f"Exported field in {time.perf_counter() - t_export:.3f}s")

            t_read = time.perf_counter()
            df = read_hfss_field(exit_field_path)
            print(f"Read field in {time.perf_counter() - t_read:.3f}s")

            df["Freq"] = freq_str
            df["Ephi"] = E_phi
            df["IWavePhi"] = i_phi
            df["IWaveTheta"] = i_theta
            df["OutgoingPower"] = power_data[(E_phi, i_phi, i_theta)]
            all_data.append(df)

            count += 1
            print(f"[w2] [{count}/{total}] φ={i_phi}deg, θ={i_theta}deg, Ephi={E_phi}, {time.perf_counter() - t_export:.3f}s")

    print(f"Second pass time: {(time.perf_counter() - t2):.3f}s")

    # Concatenate everything
    waveguide_df = pd.concat(all_data, ignore_index=True)
    waveguide_df["IngoingPower"] = get_incoming_power(Ei, plane_wave_face)

    priority_cols = ["Freq", "Ephi", "IWavePhi", "IWaveTheta", "OutgoingPower", "IngoingPower"]
    waveguide_df = waveguide_df[[*priority_cols, *[col for col in waveguide_df.columns if col not in priority_cols]]]

    if csv:
        folder_name = f"{project_name}_{design_name}_{frequency}GHz_Ephi={E_phi}"
        folder_path = os.path.join(output_file_location, folder_name)
        os.makedirs(folder_path, exist_ok=True)

        waveguide_filename = "waveguide.csv"
        waveguide_df.to_csv(os.path.join(folder_path, waveguide_filename), index=False)

    print(f"[{freq_str}] Total time: {(time.perf_counter() - start_time) / 60:.3f} min")
    return waveguide_df

def extract_far_field_data(frequency, i_theta_lower, i_theta_upper, i_phi_lower, i_phi_upper,
                           i_theta_step, i_phi_step, rad_theta_step, rad_phi_step, output_file_location, E_phi, csv=True, refine=False):
    """
    Extracts the far field total electric field data from HFSS simulations over specified incident angles,
    returning a Pandas DataFrame with complex electric field components for each combination of input parameters.

    Parameters:
        frequency (float): Frequency in GHz for which the data is extracted.
        i_theta_lower (float): Lower bound of the incident theta angle sweep (degrees).
        i_theta_upper (float): Upper bound of the incident theta angle sweep (degrees).
        i_phi_lower (float): Lower bound of the incident phi angle sweep (degrees).
        i_phi_upper (float): Upper bound of the incident phi angle sweep (degrees).
        i_theta_step (float): Step size for incident theta angle sweep (degrees).
        i_phi_step (float): Step size for incident phi angle sweep (degrees).
        rad_theta_step (float): Step size for far field theta sampling (degrees).
        rad_phi_step (float): Step size for far field phi sampling (degrees).
        output_file_location (str): Directory path to save exported data and CSV files.
        E_phi (float or int): Electric field polarization parameter for the incident wave.
        csv (bool, optional): If True, exports the resulting DataFrame to a CSV file. Defaults to True.
        refine (bool, optional): If True, uses refined adaptive solution data. Defaults to False.

    Returns:
        pandas.DataFrame: DataFrame containing the far field data with columns:
            - Freq: Frequency string (e.g., "100GHz")
            - Ephi: Electric field polarization parameter
            - IWavePhi: Incident wave phi angle (degrees)
            - IWaveTheta: Incident wave theta angle (degrees)
            - Phi: Far field observation phi angle (degrees)
            - Theta: Far field observation theta angle (degrees)
            - rEphi_real, rEphi_imag: Real and imaginary parts of far field E_phi component
            - rEtheta_real, rEtheta_imag: Real and imaginary parts of far field E_theta component

    Notes:
        - The function internally launches a monitoring thread to track the export file creation and completion,
          ensuring that the file is fully written before reading.
        - Timing printouts are included to help with performance debugging.
    """
    freq_str = f"{frequency}GHz"
    setup_name = get_setup_name(frequency) + ("_refine : LastAdaptive" if refine else " : LastAdaptive")
    sphere_name = get_radiation_sphere_name(frequency)
    far_field_path = os.path.join(output_file_location, f"{project_name}_farfield_{freq_str}.ffd")

    i_theta_num, i_phi_num = get_incoming_phi_theta_num(
        i_theta_lower, i_theta_upper, i_phi_lower, i_phi_upper, i_theta_step, i_phi_step
    )
    i_phi_values = np.linspace(i_phi_lower, i_phi_upper, i_phi_num)
    i_theta_values = np.linspace(i_theta_lower, i_theta_upper, i_theta_num)
    total = len(i_phi_values) * len(i_theta_values)
    all_data = []
    count = 0
    start_time = time.perf_counter()

    print(f"Extracting far field data for frequency {freq_str}" + (f" and Ephi = {E_phi}" if E_phi is not None else ""))

    # Internal helper function to extract far field data
    def evaluate_fields_and_export(Ephi, i_phi, i_theta):
        i_phi_str = f"{i_phi}deg"
        i_theta_str = f"{i_theta}deg"
        t_export_start = time.perf_counter()

        # We add timing here for debugging purposes, monitored via threading
        file_state: dict[str, Optional[float]] = {"created": None, "written": None}

        if os.path.exists(far_field_path):
            os.remove(far_field_path)

        def monitor_file():
            # Wait for file creation
            while not os.path.exists(far_field_path):
                time.sleep(0.001)
            file_state["created"] = time.perf_counter()
            print(f"File created at {file_state['created'] - t_export_start:.3f} s")

            stable_count = 0
            last_size = -1
            started_writing = False
            sleep_interval = 0.005  # 5 ms

            while stable_count < 2:
                time.sleep(sleep_interval)
                size = os.path.getsize(far_field_path)

                if size > 0 and not started_writing:
                    file_state["started_writing"] = time.perf_counter()
                    print(
                        f"File started writing at {file_state['started_writing'] - t_export_start:.3f} s, size: {size} bytes")
                    started_writing = True

                # Only consider stable if size > 0 and unchanged
                if size == last_size and size > 0:
                    stable_count += 1
                else:
                    stable_count = 0
                    last_size = size

            file_state["written"] = time.perf_counter()
            print(f"File write stabilized at {file_state['written'] - t_export_start:.3f} s, size: {last_size} bytes")

        monitor_thread = threading.Thread(target=monitor_file)
        monitor_thread.start()

        oModuleRad.ExportFieldsToFile([
            "ExportFileName:=", far_field_path,
            "SetupName:=", sphere_name,
            "IntrinsicVariationKey:=", f"Freq='{freq_str}' IWavePhi='{i_phi_str}' IWaveTheta='{i_theta_str}'",
            "DesignVariationKey:=", f"Ephi='{Ephi}'",
            "SolutionName:=", setup_name,
            "Quantity:=", ""
        ])

        monitor_thread.join()
        t_export_end = time.perf_counter()

        print(f"Total ExportFieldsToFile() time: {(t_export_end - t_export_start):.3f} s")

        t_read = time.perf_counter()
        df = read_hfss_far_field(far_field_path, rad_theta_step, rad_phi_step)
        print(f"Read far field in {time.perf_counter() - t_read:.3f}s")

        df["Freq"] = freq_str
        df["Ephi"] = int(Ephi)
        df["IWavePhi"] = i_phi
        df["IWaveTheta"] = i_theta
        return df

    # Loop through all input combinations and evaluate the far field
    for i_phi in i_phi_values:
        for i_theta in i_theta_values:
            t1 = time.perf_counter()
            df = evaluate_fields_and_export(E_phi, i_phi, i_theta)
            all_data.append(df)
            count += 1
            print(f"[f1] [{count}/{total}] φ={i_phi}deg, θ={i_theta}deg, Ephi={E_phi}, {time.perf_counter() - t1:.3f}s")

    far_field_df = pd.concat(all_data, ignore_index=True)

    far_field_df = far_field_df[[
        "Freq", "Ephi", "IWavePhi", "IWaveTheta",
        "Phi", "Theta", "rEphi_real", "rEphi_imag", "rEtheta_real", "rEtheta_imag"
    ]]

    if csv:
        folder_name = f"{project_name}_{design_name}_{frequency}GHz_Ephi={E_phi}"
        folder_path = os.path.join(output_file_location, folder_name)
        os.makedirs(folder_path, exist_ok=True)

        far_field_filename = "far_field.csv"
        far_field_df.to_csv(os.path.join(folder_path, far_field_filename), index=False)

    print(f"[{freq_str}] Done in {(time.perf_counter() - start_time ) /60:.3f} min")
    return far_field_df

def read_hfss_field(filepath):
    """
   Reads a HFSS .fld file containing the electric field data at the exit of the waveguide
   and converts it into a pandas DataFrame.

   Assumes the file format has two header lines, followed by data lines each containing
   nine floating-point numbers representing:
   X, Y, Z coordinates and the real and imaginary parts of the electric field components:
   Ex_real, Ey_real, Ez_real, Ex_imag, Ey_imag, Ez_imag.

   Lines containing 'nan' values in the last column are skipped.

   After reading, the function attempts to delete the file from disk to save memory.

   Parameters:
       filepath (str): Path to the .fld file to be read.

   Returns:
       pandas.DataFrame: DataFrame with columns ['X', 'Y', 'Z',
       'Ex_real', 'Ey_real', 'Ez_real', 'Ex_imag', 'Ey_imag', 'Ez_imag'] containing
       the field data.
   """
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Skip the first 2 header lines
    data_lines = lines[2:]

    # Read numeric values into a (N, 6) NumPy array: X, Y, Z, Ex_real, Ey_real, Ez_real, Ex_imag, Ey_imag, Ez_imag
    data = [
        list(map(float, line.strip().split()))
        for line in data_lines
        if line.strip() and line.strip().split()[-1].lower() != 'nan'
    ]

    df = pd.DataFrame(data, columns=['X', 'Y', 'Z', 'Ex_real', 'Ey_real', 'Ez_real', 'Ex_imag', 'Ey_imag', 'Ez_imag'])

    # Delete the file from the computer for memory purposes
    try:
        os.remove(filepath)
    except FileNotFoundError:
        print(f"File not found, cannot delete: {filepath}")

    return df

def read_hfss_far_field(filepath, rad_theta_step, rad_phi_step):
    """
    Reads a HFSS .ffd far field file and returns a DataFrame with far field electric field
    components, adding columns for real and imaginary parts of rEphi and rEtheta.

    The file format is expected as follows:
    - Line 0: theta_start, theta_end, number of theta points
    - Line 1: phi_start, phi_end, number of phi points
    - Lines 4 onward: field data for each (theta, phi) point, with four floats per line:
      rEtheta_real, rEtheta_imag, rEphi_real, rEphi_imag

    The function constructs a meshgrid of theta and phi values based on start, end, and step sizes,
    and flattens the grid alongside the field data into a DataFrame.

    Parameters:
        filepath (str): Path to the .ffd far field file.
        rad_theta_step (float): Angular step size for theta (degrees).
        rad_phi_step (float): Angular step size for phi (degrees).

    Returns:
        pandas.DataFrame: DataFrame with columns ["Phi", "Theta", "rEphi_real", "rEphi_imag",
        "rEtheta_real", "rEtheta_imag"], where Phi and Theta are angles in degrees.
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Parse angle metadata
    theta_start, theta_end, n_theta = map(float, lines[0].split())
    phi_start, phi_end, n_phi = map(float, lines[1].split())
    n_theta = int(n_theta)
    n_phi = int(n_phi)

    theta_vals = np.arange(theta_start, theta_end + 1e-8, rad_theta_step)
    phi_vals = np.arange(phi_start, phi_end + 1e-8, rad_phi_step)

    theta_grid, phi_grid = np.meshgrid(theta_vals, phi_vals, indexing="ij")

    # Parse field data lines (line 5 onward, due to the format of the file)
    data_lines = lines[4:]
    assert len(data_lines) == n_theta * n_phi, "Mismatch in expected vs actual data lines"

    rEtheta_real = []
    rEtheta_imag = []
    rEphi_real = []
    rEphi_imag = []

    for line in data_lines:
        parts = list(map(float, line.strip().split()))
        rEtheta_real.append(parts[0])
        rEtheta_imag.append(parts[1])
        rEphi_real.append(parts[2])
        rEphi_imag.append(parts[3])

    # Build DataFrame
    df = pd.DataFrame({
        "Phi": phi_grid.ravel(),
        "Theta": theta_grid.ravel(),
        "rEphi_real": rEphi_real,
        "rEphi_imag": rEphi_imag,
        "rEtheta_real": rEtheta_real,
        "rEtheta_imag": rEtheta_imag,
    })

    return df
#%%
def main():
    """
    Main entry point for the waveguide simulation workflow.

    This function orchestrates the entire simulation process, including:
    - Initializing key variables and geometry references.
    - Creating local coordinate systems and adding necessary calculator expressions.
    - For the given frequency:
        - Creating HFSS simulation setup
        - Setting up far field radiation boundaries.
        - Running the analysis including field extractions and data exports.

    This function serves as the high-level control flow to execute the full simulation, data extraction, 
    and cleanup process.

    Returns:
        None
    """

    # Simulation begins here
    initialize_variables(Ei)
    plane_wave_face, outgoing_face = get_faces_from_face_id(ingoing_face_id, outgoing_face_id)
    create_local_coordinate_system(outgoing_face, outgoing_face_cs_x, outgoing_face_cs_y)
    add_outgoing_power_to_calculator()

    if not import_from_existing_csv:
        print(f"Frequency: {frequency}GHz")
        create_setup(frequency, max_delta_E, max_passes)
        setup_radiation(rad_params, frequency)
        run_analysis(num_cores, max_delta_E, max_passes, plane_wave_face, Ei, output_file_location, i_theta_lower, i_theta_upper, i_phi_lower,
                     i_phi_upper, frequency, i_theta_step, i_phi_step, rad_params, sweep, max_difference)
    else:
        # We use regex to find the Ephi and frequency from the file path
        parent_folder = os.path.basename(os.path.dirname(waveguide_data_csv))
        m = re.search(r"_(\d+)GHz_Ephi=(\d+)$", parent_folder)
        if not m:
            raise ValueError(f"Cannot parse frequency/Ephi from '{parent_folder}'")

        freq = int(m.group(1))
        E_phi = int(m.group(2))

        hfss["Ephi"] = E_phi
        setup_radiation(rad_params, frequency)
        run_analysis(num_cores, max_delta_E, max_passes, plane_wave_face, Ei, output_file_location, i_theta_lower,
                     i_theta_upper, i_phi_lower, i_phi_upper, freq, i_theta_step, i_phi_step, rad_params,
                     sweep, max_difference, True, E_phi, waveguide_data_csv,
                     far_field_data_csv)

    hfss.release_desktop(close_projects=False, close_desktop=False)

if __name__ == "__main__":
    main()


