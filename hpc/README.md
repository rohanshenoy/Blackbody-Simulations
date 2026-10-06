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

## 4. Single-angle baseline (batch)

The first use of the HFSS solver licence. One paste on the login node:

    mkdir -p /home/rshenoy/BBRSim/outputs      # Slurm cannot create its --output directory
    cd /home/rshenoy/BBRSim/Blackbody-Simulations && git pull --ff-only fork linux-hpc-migration
    sbatch -t 02:00:00 hpc/run_frequency.sbatch configs/crack1_500GHz_single_angle.toml --job-id baseline3
    squeue -u rshenoy

Expected: the job ends `COMPLETED`; `/home/rshenoy/BBRSim/outputs/slurm-<id>.out` ends with
`done; manifest at .../job_baseline3/manifest.json`; and
`/home/rshenoy/BBRSim/outputs/InfParallelPlate_crack1Rohan_500GHz/job_baseline3/` holds
`InfParallelPlate_crack1Rohan_500GHz_Ephi=0/{waveguide.csv,far_field.csv,manifest.json}`,
`InfParallelPlate_crack1Rohan_500GHz.dataset.json`, `logs/run.log` and `logs/convergence_Ephi0.txt`
(this config extracts Ephi=0 only; the solve covers both polarizations).
A job id is used once: resubmitting with an existing `--job-id` stops with `FileExistsError` before AEDT
starts, because the job directory exists. Pick a new id (`baseline4`, ...).

Then check the grid sizes and the transmission:

    python - <<'PY'
    import pandas as pd
    d = "/home/rshenoy/BBRSim/outputs/InfParallelPlate_crack1Rohan_500GHz/job_baseline3/InfParallelPlate_crack1Rohan_500GHz_Ephi=0"
    w = pd.read_csv(f"{d}/waveguide.csv"); f = pd.read_csv(f"{d}/far_field.csv")
    print(len(w), len(f), w.OutgoingPower.iloc[0] / w.IngoingPower.iloc[0])
    PY

Expected: `5151 19388 <T>`, with T near 1.0: the step 4a solve of the same setup gave 1.0017. The
Windows reference (AEDT 2023 R2, 15 directions) gives 1.0545 at this key and 0.9961 at the physically
identical key (90 deg, 180 deg) for Ephi=1; its two polarizations come from separate adaptive meshes, so
a spread of about 6 % is the method's own noise at MaxDeltaE 0.02. Record T, the number of passes and
final delta E from the convergence file, the solve time from `logs/run.log`,
`sacct -j <id> -o JobID,State,ExitCode,Elapsed,MaxRSS` and `du -sh <job> <job>/project/*.aedtresults`.

History: `baseline1` (2026-10-05) and `baseline2` (2026-10-06) solved and failed at the first field
evaluation: a single-direction plane wave defines no `IWavePhi`/`IWaveTheta`, and the runner named them
(RESULTS.md, step 4). Fixed; use a new job id for every attempt.

### 4a. Solve-and-probe diagnostic

When a run fails after its solve, `hpc/diagnose_calculator.py` copies the job's project and results,
lists the results folder, asks AEDT what it holds, then evaluates `outgoing_power` with the runner's
arguments and with variants, exports the exit and far fields alone, and runs the runner's own
extraction on the copy, writing the CSV tables beside it. With `--solve` it first solves the copy
(HFSS solver licence) and saves it; use that when the job's project holds no saved solution (runs
made before the runner saved right after its solve, baseline1 and baseline2 among them). Debug QOS:

    sbatch -A golwala -p expansion -q debug -N 1 -c 4 --mem=16G -t 00:30:00 \
      -o /home/rshenoy/BBRSim/outputs/diag-%j.out \
      --wrap 'bash -lc "source /home/rshenoy/BBRSim/bb_env.sh && cd /home/rshenoy/BBRSim/Blackbody-Simulations && python hpc/diagnose_calculator.py --solve --job /home/rshenoy/BBRSim/outputs/InfParallelPlate_crack1Rohan_500GHz/job_baseline2"'

Expected at the end of `/home/rshenoy/BBRSim/outputs/diag-<id>.out`: a summary with one OK or FAILED
line per probe, `FIELDS: present`, `WORKING POWER VARIANTS: [...]` naming the argument forms that
evaluate (each with T near 1.054 at normal incidence), and `EXTRACTION PASS` or `FAIL`. The AEDT
message printed under a failed probe is the reason.

## 5. Full reference run (batch) and comparison

The 15-direction crack1 run at 500 GHz, then the comparison with the Windows reference as a dependent
job. The reference CSVs are the ones tracked in the HPC BBRSimulation checkout (`main` holds the same
bytes as the Mac copies); `hpc/compare_reference.sbatch` checks their sha256 against
`configs/legacy_reference_500GHz.sha256` before comparing. One paste on the login node:

    cd /home/rshenoy/BBRSim/Blackbody-Simulations && git pull --ff-only fork linux-hpc-migration
    A=$(sbatch --parsable -t 02:00:00 hpc/run_frequency.sbatch configs/crack1_500GHz_reference.toml --job-id reference2); A=${A%%;*}; echo "run job $A"
    sbatch --dependency=afterok:$A --kill-on-invalid-dep=yes hpc/compare_reference.sbatch \
        /home/rshenoy/BBRSim/outputs/InfParallelPlate_crack1Rohan_500GHz/job_reference2
    squeue -u rshenoy

When both jobs have left the queue:

    tail -8 /home/rshenoy/BBRSim/outputs/slurm-<run id>.out
    cat /home/rshenoy/BBRSim/outputs/slurm-<compare id>.out

Expected: the run ends with `done; manifest at .../job_reference2/manifest.json`; the comparison prints
four reference files `OK`, a table and `RESULT: PASS` for each polarization, and
`COMPARISON PASS`. The JSON records are `<job>/comparison_Ephi{0,1}.json`. If the run fails, Slurm
cancels the comparison (`--kill-on-invalid-dep`). For per-key detail on the login node:
`python hpc/compare_detail.py <job>`. History: `reference1` (job 4069263) carried Ephi=0's power in its
Ephi=1 table (RESULTS.md step 5); fixed, so the run is `reference2`.

Transmission tolerance: |dT| <= max(1e-3, 0.05 T). The reference's own two physically identical keys,
(0 deg, 180 deg) for Ephi=0 and (90 deg, 180 deg) for Ephi=1, differ by 5.9 % (separate adaptive meshes per
polarization; RESULTS.md step 4), so a miss of a few percent at a few keys is within the method's noise.
Report the values and decide the tolerance explicitly rather than loosening it silently. Field
distribution lines are warnings only.

Record what this run used, the multi-frequency sizing rule depends on it:

    sacct -j <run id> -o JobID,State,ExitCode,Elapsed,MaxRSS
    du -sh <job> <job>/project/*.aedtresults

## 5b. Incident-direction check (batch, no new solve)

After step 5, on a copy of its solved project (the script switches the field type while exporting and
never saves, but the job's own project stays untouched this way). One paste on the login node:

    J=/home/rshenoy/BBRSim/outputs/InfParallelPlate_crack1Rohan_500GHz/job_reference2
    mkdir $J/incident_check && cp -r $J/project/ParallelPlateGaps.aedt $J/project/ParallelPlateGaps.aedtresults $J/incident_check/
    sbatch -A golwala -p expansion -q debug -N 1 -c 4 --mem=16G -t 00:30:00 \
        -o /home/rshenoy/BBRSim/outputs/slurm-%j.out \
        --wrap "bash -lc 'source /home/rshenoy/BBRSim/bb_env.sh && cd /home/rshenoy/BBRSim/Blackbody-Simulations && python hpc/check_incident_direction.py --project $J/incident_check/ParallelPlateGaps.aedt --design parallel_plate_gap_50um_500GHz --out $J/incident_direction.json'"

Then `cat $J/incident_direction.json`.

Expected: exit 0, `"verdict": "PASS: arrival direction, k = -r_hat(theta, phi)"`, `ratio_ky_over_kz`
near -1 and `magnitude_over_k0` near 1. `ALTERNATIVE` (ratio near +1) means HFSS uses the other
reading: stop and tell the BBRsim side, whose azimuth formula depends on it. If `EditSources` rejects
`ScatteredFields`, paste the error; Rohan then reads the incident wave direction once in the GUI
(Open OnDemand) and records it here.

Also expected: `"slope_z_negative_as_expected": true`. It is a self-check outside the verdict and the
exit code: the ratio cannot tell k from -k, and under HFSS's e^{+jwt} a wave entering the gap needs a
negative z slope. `false` means the wave runs out of the gap or the phasors use e^{-jwt}: stop and
tell the BBRsim side, as for `ALTERNATIVE`.

## 6. Build the round-gap project (batch, debug QOS)

    sbatch -A golwala -p expansion -q debug -N 1 -c 4 --mem=16G -t 00:30:00 \
        -o /home/rshenoy/BBRSim/outputs/slurm-%j.out \
        --wrap 'bash -lc "source /home/rshenoy/BBRSim/bb_env.sh && cd /home/rshenoy/BBRSim/Blackbody-Simulations && python build_hfss_project.py --spec configs/geometries/round_gap_r50um.toml"'

Expected: `Build verified: /resnick/home/rshenoy/BBRSim/projects/RoundGap.aedt` (or the same path under
`/home`), exit 0, `RoundGap.build.json` with `"verification_differences": []`, and
`RoundGap.inventory.json` listing one design `round_gap_r50um` with one object `gap` [vacuum],
bounding box [-0.05, -0.05, 0, 0.05, 0.05, 0.4] and 3 faces. PyAEDT warnings "Not enough vertices or
non-planar face" are expected: the two vertex-free end caps and the curved side have no computable normal.

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
    print(sum("all-nan rows" in line for line in open(f"{job}/logs/run.log")))
    PY

Expected: `N 1369 <T>` with N from 7825 to 7845 (twenty lattice points lie exactly on the rim, four on
the axes and sixteen off them; 7845 means HFSS kept them all, 7825 none; record the exact N), T strictly
between 0 and 1, then the same N,
`{'shape': 'disc', 'radius_m': 5e-05}` (to rounding), `TE11 1`, then `0`. A non-zero last line means
`logs/run.log` holds the reader's warning that HFSS wrote outside lattice points as rows other than
nine tokens; those points are still skipped and the run stands, but copy the warning into RESULTS.md.
Record T as the first round-gap reference value (no Windows reference exists), the passes and final
delta E from `logs/convergence_Ephi0.txt`, and `sacct -j <id> --format=MaxRSS,Elapsed`.

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

Expected: two lines with 15 x N rows each (N as in step 7), and the two T values at normal incidence equal
to within the convergence tolerance (a few percent at MaxDeltaE 0.02): a round gap cannot prefer a
polarization at normal incidence. A large difference means a frame or polarization error; stop and
report. Record memory and wall time; they size the job array (next-steps item 5).
