# HPC verification log

Fill each section with the exact command, date, Slurm job ID, and pasted output.
A section without pasted output is not done. Steps are defined in `hpc/README.md`.

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

## 3. Prepared project

## 4. Single-angle baseline

## 5. Full reference run and comparison

## Discrepancies and decisions
