# Running on the Caltech HPC cluster

Every HPC command below runs in `/home/rshenoy/BBRSim/Blackbody-Simulations` after
`source /home/rshenoy/BBRSim/bb_env.sh`. Nothing here needs a display, X forwarding,
or Open OnDemand. Paste each step's output into `hpc/RESULTS.md` before moving on.

## 0. Get the branch onto HPC

`origin` is Jason's repository, which we can only read, so the branch travels as a git bundle
(no publishing) or via a fork under your own GitHub account.

Bundle route. On the Mac (`/Users/rohanshenoy/supercdms/Blackbody-Simulations`):

    git bundle create ~/Desktop/bbsim-linux-hpc-migration.bundle main..linux-hpc-migration
    scp ~/Desktop/bbsim-linux-hpc-migration.bundle rshenoy@login.hpc.caltech.edu:/home/rshenoy/BBRSim/

On HPC (login node is fine; no solver runs here):

    cd /home/rshenoy/BBRSim/Blackbody-Simulations
    git status --short                      # expect nothing
    git fetch origin                        # main must be present (the bundle builds on it)
    git bundle verify ../bbsim-linux-hpc-migration.bundle
    git fetch ../bbsim-linux-hpc-migration.bundle linux-hpc-migration:linux-hpc-migration
    git checkout linux-hpc-migration

To pick up later commits, rebuild and copy the bundle, then on HPC
`git fetch ../bbsim-linux-hpc-migration.bundle linux-hpc-migration:linux-hpc-migration && git checkout linux-hpc-migration && git reset --hard linux-hpc-migration`.

Fork route. On the Mac, `gh repo fork ModerJason/Blackbody-Simulations --remote --remote-name fork`
then `git push -u fork linux-hpc-migration`; on HPC, `git remote add fork <fork url>`,
`git fetch fork`, `git checkout -b linux-hpc-migration fork/linux-hpc-migration`.

Then, on HPC:

    source /home/rshenoy/BBRSim/bb_env.sh
    python -m pip freeze > requirements-hpc.txt      # commit this: it pins the tested combination
    python -c "import ansys.aedt.core as c; print(c.__version__)"
    python -c "import bbsim.run_frequency, bbsim.prepare_project, bbsim.compare_exports; print('bbsim imports ok')"

Expected: the PyAEDT version is printed and the imports succeed.

## 1. Headless AEDT smoke test (debug QOS, on a compute node)

    srun -A golwala -p expansion -q debug -N 1 -c 4 --mem=16G -t 00:30:00 --pty bash
    cd /home/rshenoy/BBRSim/Blackbody-Simulations && source /home/rshenoy/BBRSim/bb_env.sh
    python hpc/aedt_smoke_test.py --workdir /home/rshenoy/BBRSim/scratch/smoke

Expected: the PyAEDT version, the `Hfss.__init__` signature, `desktop started in ...s`,
`versions: {...'aedt': '2025.2...'}`, `saved: ... True`, `closed cleanly`.
If a warning says arguments were dropped from `Hfss()`, record them in RESULTS.md.
A licensing failure shows up here, before any solve is attempted.

## 2. Inventory the reference project (dry run, same allocation)

    mkdir -p /home/rshenoy/BBRSim/projects
    python prepare_hfss_project.py --input InfParallelPlate.aedt \
        --output /home/rshenoy/BBRSim/projects/ParallelPlateGaps.aedt \
        --mapping configs/design_mapping.json --dry-run

Expected: 18 designs listed, `crack1Rohan/crack1 -> parallel_plate_gap_50um/gap`,
`crack2/crack2 -> parallel_plate_gap_100um/gap`, 16 deletions, no ERRORS, exit 0.
`ParallelPlateGaps.inventory.source.json` is written next to the output path.

## 3. Prepare the cleaned project (same allocation)

    python prepare_hfss_project.py --input InfParallelPlate.aedt \
        --output /home/rshenoy/BBRSim/projects/ParallelPlateGaps.aedt \
        --mapping configs/design_mapping.json

Expected: `Prepared project verified`, exit 0, and `ParallelPlateGaps.mapping.json` with an
empty `verification_differences` list. Confirm the reference is untouched:

    sha256sum InfParallelPlate.aedt   # a264705100ab80ff2a3b2315b0529a8505c2903b996c98cd0d65d7144097482f

## 4. Single-angle baseline (same allocation, or a new debug one)

    python run_hfss_frequency.py --config configs/crack1_500GHz_single_angle.toml --job-id baseline1

Expected: exit 0, and
`/home/rshenoy/BBRSim/outputs/InfParallelPlate_crack1Rohan_500GHz/job_baseline1/` containing
`InfParallelPlate_crack1Rohan_500GHz_Ephi=0/{waveguide.csv,far_field.csv,manifest.json}`,
`logs/run.log` and `logs/convergence_Ephi0.txt`. Check the grid sizes and transmission:

    python - <<'PY'
    import pandas as pd
    d = "/home/rshenoy/BBRSim/outputs/InfParallelPlate_crack1Rohan_500GHz/job_baseline1/InfParallelPlate_crack1Rohan_500GHz_Ephi=0"
    w = pd.read_csv(f"{d}/waveguide.csv"); f = pd.read_csv(f"{d}/far_field.csv")
    print(len(w), len(f), w.OutgoingPower.iloc[0] / w.IngoingPower.iloc[0])
    PY

Expected: `5151 19388 <T>`, with T close to the Windows reference 1.0545 (AEDT 2023 R2).
Record T, the number of passes and final delta E from the convergence file, and the solve time.

## 5. Full reference run (batch) and comparison

On the Mac, copy the reference datasets to HPC:

    scp -r /Users/rohanshenoy/geant4/BBRSimulation/data/waveguides/InfParallelPlate_crack1Rohan_500GHz_Ephi=* \
        rshenoy@login.hpc.caltech.edu:/home/rshenoy/BBRSim/reference/

(Create `/home/rshenoy/BBRSim/reference` on HPC first; use your usual login host if it differs.)

On HPC:

    mkdir -p /home/rshenoy/BBRSim/outputs      # Slurm cannot create its --output directory
    sbatch hpc/run_frequency.sbatch configs/crack1_500GHz_reference.toml
    squeue -u rshenoy

When it finishes, with `<job>` the job directory named in `/home/rshenoy/BBRSim/outputs/slurm-<id>.out`:

    for e in 0 1; do
      python compare_hfss_exports.py "<job>/InfParallelPlate_crack1Rohan_500GHz_Ephi=$e" \
          "/home/rshenoy/BBRSim/reference/InfParallelPlate_crack1Rohan_500GHz_Ephi=$e" \
          --json "<job>/comparison_Ephi$e.json"
    done

Expected: `RESULT: PASS` for both polarizations. Paste both tables into RESULTS.md. If
transmission fails only at a few angles, report the values and decide the tolerance
explicitly rather than loosening it silently. Field-distribution lines are warnings only.

Memory (32G) and time (8h) in the batch script are first guesses; set them from the
baseline's `sacct -j <id> --format=MaxRSS,Elapsed` before larger runs.
