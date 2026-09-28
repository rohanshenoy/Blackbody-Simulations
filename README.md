# HFSS Simulations – Blackbody Radiation in Cryogenic Detector Gaps
This repository contains HFSS models and scripts for simulating how blackbody photons from warmer components propagate through small mechanical gaps ($$O(10^2–10^3 μm)$$) in the dark matter detector cryostats, especially in relation to the Cryogenic Dark Matter Search (CDMS).

These gaps behave like semi-open parallel-plate waveguides that can support low-frequency TEM modes with negligible cutoff, allowing photons to reach ultra-cold detectors and cause photoionization or leakage currents.

**We model:**
- Gap geometries and dimensions from mechanical designs
- Plane-wave excitation over blackbody-relevant frequencies
- |S21| parameters, field distributions, and far-field emission

**Workflow (Linux/HPC, headless):**

1. Prepare a cleaned project once from the reference `InfParallelPlate.aedt`:
   `python prepare_hfss_project.py --input InfParallelPlate.aedt --output <dir>/ParallelPlateGaps.aedt --mapping configs/design_mapping.json`
2. Run one frequency from a config file, inside a Slurm allocation:
   `python run_hfss_frequency.py --config configs/crack1_500GHz_reference.toml`
3. Compare against a reference dataset:
   `python compare_hfss_exports.py <job>/InfParallelPlate_crack1Rohan_500GHz_Ephi=0 <reference dir>`
4. Post-process: `python postprocess.py <Ephi=0 dir> <Ephi=1 dir>`, `python plot.py <dataset dir> --out-dir figs`.

See `hpc/README.md` for the cluster environment.
For the original geometry and simulation guide, see [here](https://github.com/ModerJason/Blackbody-Simulations/blob/main/Blackbody%20Simulations%20Usage%20Guide.pdf).
The original Windows script is preserved in `legacy/`. The usage guide PDF refers to `bbsim.py`;
the actual legacy script is `legacy/bbsim1freq.py`.

These results will later integrate with Geant4 open-space photon simulations for an end-to-end model.

For more details on analytical calculations, validation for HFSS simulations, and plots, see [here](https://github.com/ModerJason/Blackbody-Simulations/blob/main/Blackbody%20Simulation%20Report.pdf).

**Requirements:**
- Ansys Electronics Desktop 2025 R2 (HFSS) with a license; PyAEDT (`ansys.aedt.core`)
- Python 3.11 with `numpy`, `pandas`, `scipy`, `matplotlib`, `seaborn`
