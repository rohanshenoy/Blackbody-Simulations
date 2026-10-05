# Running on the Caltech HPC cluster

Every HPC command below runs in `/home/rshenoy/BBRSim/Blackbody-Simulations` after
`source /home/rshenoy/BBRSim/bb_env.sh`. Nothing here needs a display, X forwarding,
or Open OnDemand. Paste each step's output into `hpc/RESULTS.md` before moving on.

## 0. Get the branch onto HPC

`origin` is Jason's repository, which we can only read. The branch lives on the fork
https://github.com/rohanshenoy/Blackbody-Simulations (remote name `fork` on the Mac).

First time on HPC:

    cd /home/rshenoy/BBRSim/Blackbody-Simulations
    git status --short                      # expect nothing
    git remote add fork https://github.com/rohanshenoy/Blackbody-Simulations.git
    git fetch fork
    git checkout -B linux-hpc-migration fork/linux-hpc-migration

Later updates on HPC (branch checked out):

    git pull --ff-only fork linux-hpc-migration

Publishing from the Mac: `git push fork linux-hpc-migration`. Never push to `origin`.

Offline alternative (no GitHub): `git bundle create ~/Desktop/bbsim-linux-hpc-migration.bundle main..linux-hpc-migration`
on the Mac, `scp` it to `/home/rshenoy/BBRSim/`, then on HPC
`git fetch ../bbsim-linux-hpc-migration.bundle linux-hpc-migration && git merge --ff-only FETCH_HEAD`.

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

    scp -r /Users/rohanshenoy/BBRsim/BBRSimulation/data/waveguides/InfParallelPlate_crack1Rohan_500GHz_Ephi=* \
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

## 5b. Incident-direction check (batch, no new solve)

After step 5, with `<job>` the step-5 job directory named in its `slurm-<id>.out`:

    sbatch -A golwala -p expansion -q debug -N 1 -c 4 --mem=16G -t 00:30:00 \
        -o /home/rshenoy/BBRSim/outputs/slurm-%j.out \
        --wrap 'source /home/rshenoy/BBRSim/bb_env.sh && cd /home/rshenoy/BBRSim/Blackbody-Simulations && python hpc/check_incident_direction.py --project <job>/project/ParallelPlateGaps.aedt --design parallel_plate_gap_50um_500GHz --out <job>/incident_direction.json'

Expected: exit 0, `"verdict": "PASS: arrival direction, k = -r_hat(theta, phi)"`, `ratio_ky_over_kz`
near -1 and `magnitude_over_k0` near 1. `ALTERNATIVE` (ratio near +1) means HFSS uses the other
reading: stop and tell the BBRsim side, whose azimuth formula depends on it. If `EditSources` rejects
`ScatteredFields`, paste the error; Rohan then reads the incident wave direction once in the GUI
(Open OnDemand) and records it here.

## 6. Build the round-gap project (batch, debug QOS)

    sbatch -A golwala -p expansion -q debug -N 1 -c 4 --mem=16G -t 00:30:00 \
        -o /home/rshenoy/BBRSim/outputs/slurm-%j.out \
        --wrap 'source /home/rshenoy/BBRSim/bb_env.sh && cd /home/rshenoy/BBRSim/Blackbody-Simulations && python build_hfss_project.py --spec configs/geometries/round_gap_r50um.toml'

Expected: `Build verified: /resnick/home/rshenoy/BBRSim/projects/RoundGap.aedt` (or the same path under
`/home`), exit 0, `RoundGap.build.json` with `"verification_differences": []`, and
`RoundGap.inventory.json` listing one design `round_gap_r50um` with one object `gap` [vacuum],
bounding box [-0.05, -0.05, 0, 0.05, 0.05, 0.4] and 3 faces. Record the PyAEDT warnings, if any, about
non-planar faces: they come from the curved side and are expected.

## 7. Round gap, single angle at 2000 GHz

    sbatch -t 02:00:00 hpc/run_frequency.sbatch configs/round_gap_r50um_2000GHz_single_angle.toml --job-id roundgap1

When it finishes:

    python - <<'PY'
    import json, pandas as pd
    job = "/home/rshenoy/BBRSim/outputs/RoundGap_r50um_2000GHz/job_roundgap1"
    d = f"{job}/RoundGap_r50um_2000GHz_Ephi=0"
    w = pd.read_csv(f"{d}/waveguide.csv"); f = pd.read_csv(f"{d}/far_field.csv")
    s = json.load(open(f"{job}/RoundGap_r50um_2000GHz.dataset.json"))
    print(len(w), len(f), w.OutgoingPower.iloc[0] / w.IngoingPower.iloc[0])
    print(s["exit_field"]["points_per_key_retained"], s["exit_field"]["cross_section"], s["modes"]["mode"], s["modes"]["propagating_count"])
    PY

Expected: `7845 1369 <T>` or `7825 1369 <T>` (record which: it tells whether HFSS evaluates the four
rim points on the axes), T strictly between 0 and 1, then `7845` or `7825`,
`{'shape': 'disc', 'radius_m': 5e-05}` (to rounding), `TE11 1`. Record T as the first round-gap
reference value (no Windows reference exists), the passes and final delta E from
`logs/convergence_Ephi0.txt`, and `sacct -j <id> --format=MaxRSS,Elapsed`.

## 8. Round gap, full sweep

    sbatch hpc/run_frequency.sbatch configs/round_gap_r50um_2000GHz_reference.toml

When it finishes, with `<job>` named in its `slurm-<id>.out`:

    python - <<'PY'
    import pandas as pd
    job = "<job>"
    for e in (0, 1):
        w = pd.read_csv(f"{job}/RoundGap_r50um_2000GHz_Ephi={e}/waveguide.csv")
        r = w[(w.IWavePhi == 0) & (w.IWaveTheta == 180)].iloc[0]
        print(e, r.OutgoingPower / r.IngoingPower, len(w))
    PY

Expected: two lines with 15 x (7845 or 7825) rows each, and the two T values at normal incidence equal
to within the convergence tolerance (a few percent at MaxDeltaE 0.02): a round gap cannot prefer a
polarization at normal incidence. A large difference means a frame or polarization error; stop and
report. Record memory and wall time; they size the job array (next-steps item 5).
