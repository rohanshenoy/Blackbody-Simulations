# HPC verification log

Fill each section with the exact command, date, Slurm job ID, and pasted output.
A section without pasted output is not done. Steps are defined in `hpc/README.md`.

## Status (2026-10-04)

- Steps 0 to 3: PASS (2026-09-28). The step 1 smoke test also passes under the
  merged `bb_env.sh` that loads Geant4 (Slurm job 3693291, 2026-09-30).
- Next: step 4. First update the HPC checkout to the fork tip and capture
  `requirements-hpc.txt` with `python -m pip freeze --exclude-editable`.
- After that: step 5.
- BBRsim (Geant4) on HPC is tracked in `ANSYS_HPC_INSTALL_README.md` section 8.

## 0. Get the branch onto HPC

2026-09-28. Branch fetched from the git bundle; checkout at 07b387c. Python 3.11.6 (GCC 13.2.0),
PyAEDT 1.7.0. `requirements-hpc.txt` still to be captured with `python -m pip freeze` on the login node.

## 1. AEDT smoke test

2026-09-28 18:19, Slurm job 3614903 (`srun -A golwala -p expansion -q debug -N 1 -c 4 --mem=16G -t 00:30:00`),
node hpc-91-17. Command: `python hpc/aedt_smoke_test.py --workdir /home/rshenoy/BBRSim/scratch/smoke`.
Result: PASS.

    pyaedt version: 1.7.0
    Hfss.__init__ signature: (self, project, design, solution_type, setup, version, non_graphical,
                              new_desktop, close_on_exit, student_version, machine, port, aedt_process_id, remove_lock)
    PyAEDT INFO: AEDT version 2025.2.
    PyAEDT INFO: Launching AEDT server with gRPC transport mode: TransportMode.UDS
    PyAEDT INFO: Electronics Desktop started on gRPC port 49125 after 42.9 seconds.
    PyAEDT INFO: AEDT installation Path /resnick/groups/golwala/rshenoy/software/ansys_em/v252/AnsysEM
    PyAEDT WARNING: Service Pack is not detected. PyAEDT is currently connecting in Insecure Mode.
    PyAEDT INFO: Non-graphical mode detected. Disabling Desktop logs.
    PyAEDT INFO: Added design 'smoke' of type HFSS.
    desktop started in 77.2s
    versions: {'pyaedt': '1.7.0', 'aedt': '2025.2.0'}
    created object: probe model units: mm
    saved: /home/rshenoy/BBRSim/scratch/smoke/smoke.aedt True
    PyAEDT INFO: Desktop has been released and closed.
    closed cleanly

Verified: headless (non-graphical, gRPC over UDS) startup on a compute node, base AEDT license checkout,
project creation and save, clean release. No arguments were dropped from `Hfss()`. Not yet verified:
HFSS solver license and a solve (step 4). Cosmetic: PyAEDT's logger and the script's logger both print,
so lines appear twice.

## 2. Dry-run inventory

2026-09-28 18:22, same allocation (job 3614903, hpc-91-17). Command:
`python prepare_hfss_project.py --input InfParallelPlate.aedt --output /home/rshenoy/BBRSim/projects/ParallelPlateGaps.aedt --mapping configs/design_mapping.json --dry-run`.
Result: inventory PASS; the 2023 R2 project opened headlessly in 2025 R2 with no conversion prompt.
AEDT started in 5.2 s (warm); the 18-design inventory took about 1.5 s.

    Designs found: 18
      crack2          units=mm  objects: crack2[vacuum]  setups=[] boundaries=0
      crack1Rohan     units=mm  objects: crack1[vacuum]  setups=[] boundaries=0
      (16 others, all HFSS, e.g. crack1Rohan_500GHz setups=['500GHz'] boundaries=2,
       backagain2_500GHz setups=['500GHz', '500GHz_refine'], infrarmitigation boundaries=3)
    Renames:
      crack1Rohan/crack1 -> parallel_plate_gap_50um/gap  (dataset_id InfParallelPlate_crack1Rohan, 50 um)
      crack2/crack2 -> parallel_plate_gap_100um/gap  (dataset_id InfParallelPlate_crack2, 100 um)
    Deletions (16): [backagain, backagain2, backagain2_500GHz, backagain_500GHz, crack1Rohan_500GHz,
      crack1again_500GHz, crack2_500GHz, cylindrical2, height_perturb, hollow, hollow_500GHz,
      infrarmitigation, rectangle, right_angle, stub, stub2]
    Dry run: no project written.

Both base designs confirmed bare (no setups, boundaries or excitations), as the runner assumes.
Two harmless PyAEDT warnings "Not enough vertices or non-planar face" while reading cylindrical2.

Defect found: after the desktop was released, Python's TemporaryDirectory cleanup raised
`OSError: [Errno 39] Directory not empty: .../projects/.prepare_tmp_vxeyq9hv` (AEDT/NFS left entries
behind), so the successful dry run exited 1. Fixed in f6ececa: the working copy now lives in the
node-local $TMPDIR, removal is retried for 60 s, and a final failure is logged rather than raised.
The leftover directory on HPC must be removed by hand (`rm -rf /home/rshenoy/BBRSim/projects/.prepare_tmp_*`).

## 3. Prepared project

2026-09-28, same allocation (job 3614903, hpc-91-17), code at f6ececa. Command:
`python prepare_hfss_project.py --input InfParallelPlate.aedt --output /home/rshenoy/BBRSim/projects/ParallelPlateGaps.aedt --mapping configs/design_mapping.json`.
Result: PASS.

    Prepared project verified: /resnick/home/rshenoy/BBRSim/projects/ParallelPlateGaps.aedt
    a264705100ab80ff2a3b2315b0529a8505c2903b996c98cd0d65d7144097482f  InfParallelPlate.aedt

The cleaned project holds parallel_plate_gap_50um and parallel_plate_gap_100um (object `gap` in each),
saved by AEDT 2025 R2; the reopen inventory matched the source designs (geometry, materials, units,
mesh settings) and the reference project's hash is unchanged. Records next to the output:
ParallelPlateGaps.inventory.source.json, ParallelPlateGaps.inventory.json, ParallelPlateGaps.mapping.json.

## 4. Single-angle baseline

## 5. Full reference run and comparison

## Discrepancies and decisions
