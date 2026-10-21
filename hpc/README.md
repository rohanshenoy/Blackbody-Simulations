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

## 5c. crack2 reference run and comparison (batch)

The 15-direction crack2 run at 500 GHz (`configs/crack2_500GHz_reference.toml`, the Windows baseline's
settings), compared with the Windows crack2 reference exactly as in step 5. Rohan decided on 2026-10-06
that BBRsim moves both cracks to 2025 R2 data in one commit if crack2 compares cleanly. One paste:

    cd /home/rshenoy/BBRSim/Blackbody-Simulations && git pull --ff-only fork linux-hpc-migration
    B=$(sbatch --parsable -t 02:00:00 hpc/run_frequency.sbatch configs/crack2_500GHz_reference.toml --job-id reference1); B=${B%%;*}; echo "crack2 run job $B"
    sbatch --dependency=afterok:$B --kill-on-invalid-dep=yes hpc/compare_reference.sbatch \
        /home/rshenoy/BBRSim/outputs/InfParallelPlate_crack2_500GHz/job_reference1 InfParallelPlate_crack2_500GHz

When both jobs have left the queue:

    grep -v "^ \|^Loading" /home/rshenoy/BBRSim/outputs/slurm-<compare id>.out
    python hpc/compare_detail.py /home/rshenoy/BBRSim/outputs/InfParallelPlate_crack2_500GHz/job_reference1 --stem InfParallelPlate_crack2_500GHz

Expected: the four crack2 reference files `OK`, structure checks PASS for both polarizations, and
transmission within 5 % wherever light passes; a constant offset of one Windows polarization, as for
crack1, is the case to look for in the detail. Record `sacct` and `du` as in step 5.

## 5d. Convergence check: both cracks at MaxDeltaE 0.005 (batch)

Steps 5 and 5c showed per-solve offsets of 3.5 % to 5.5 % in T at the reference target (Delta Mag
Energy 0.02), in both directions. Rohan decided on 2026-10-06 to rerun both cracks once at 0.005
(`configs/crack{1,2}_500GHz_reference_tight.toml`, up to 20 passes; nothing else differs) before the
two-crack data tree is built. One paste:

    cd /home/rshenoy/BBRSim/Blackbody-Simulations && git pull --ff-only fork linux-hpc-migration
    A=$(sbatch --parsable -t 03:00:00 hpc/run_frequency.sbatch configs/crack1_500GHz_reference_tight.toml --job-id tight1); A=${A%%;*}; echo "crack1 run job $A"
    B=$(sbatch --parsable -t 03:00:00 hpc/run_frequency.sbatch configs/crack2_500GHz_reference_tight.toml --job-id tight1); B=${B%%;*}; echo "crack2 run job $B"
    squeue -u rshenoy

When both have left the queue, on the login node (each tight run against its 0.02 run, then against
the Windows reference; sections 1 to 3 of the detail report):

    D=/home/rshenoy/BBRSim/outputs; C1=InfParallelPlate_crack1Rohan_500GHz; C2=InfParallelPlate_crack2_500GHz
    python hpc/compare_detail.py $D/$C1/job_tight1 --reference-root $D/$C1/job_reference2 | sed -n '/== 1/,/== 4/p'
    python hpc/compare_detail.py $D/$C1/job_tight1 | sed -n '/== 2/,/== 4/p'
    python hpc/compare_detail.py $D/$C2/job_tight1 --stem $C2 --reference-root $D/$C2/job_reference1 | sed -n '/== 1/,/== 4/p'
    python hpc/compare_detail.py $D/$C2/job_tight1 --stem $C2 | sed -n '/== 2/,/== 4/p'
    sacct -j <crack1 run job>,<crack2 run job> -o JobID,State,ExitCode,Elapsed,MaxRSS

Expected: both runs end `done`; passes completed below 20 with a last delta under 0.005; the T change
from the 0.02 run to the 0.005 run at the transmitting keys is the 0.02 run's convergence error, and
the two polarizations at normal incidence stay equal. The tree then takes the better-converged runs.

## 6. Build the round-gap project (batch, debug QOS)

Steps 6 to 8 can go in as one chain, each starting only when the previous one succeeded
(`--kill-on-invalid-dep` cancels the rest after a failure). Both round-gap configs use the production
convergence target, MaxDeltaE 0.005 with up to 20 passes (Rohan, 2026-10-06, ledger HF-027):

    cd /home/rshenoy/BBRSim/Blackbody-Simulations && git pull --ff-only fork linux-hpc-migration
    S6=$(sbatch --parsable -A golwala -p expansion -q debug -N 1 -c 4 --mem=16G -t 00:30:00 -o /home/rshenoy/BBRSim/outputs/slurm-%j.out \
      --wrap 'bash -lc "source /home/rshenoy/BBRSim/bb_env.sh && cd /home/rshenoy/BBRSim/Blackbody-Simulations && python build_hfss_project.py --spec configs/geometries/round_gap_r50um.toml"'); S6=${S6%%;*}; echo "step 6 job $S6"
    S7=$(sbatch --parsable --dependency=afterok:$S6 --kill-on-invalid-dep=yes -t 02:00:00 hpc/run_frequency.sbatch configs/round_gap_r50um_2000GHz_single_angle.toml --job-id roundgap1); S7=${S7%%;*}; echo "step 7 job $S7"
    S8=$(sbatch --parsable --dependency=afterok:$S7 --kill-on-invalid-dep=yes -t 03:00:00 hpc/run_frequency.sbatch configs/round_gap_r50um_2000GHz_reference.toml --job-id reference1); S8=${S8%%;*}; echo "step 8 job $S8"
    squeue -u rshenoy

Each step's own command and its checks follow.

    sbatch -A golwala -p expansion -q debug -N 1 -c 4 --mem=16G -t 00:30:00 \
        -o /home/rshenoy/BBRSim/outputs/slurm-%j.out \
        --wrap 'bash -lc "source /home/rshenoy/BBRSim/bb_env.sh && cd /home/rshenoy/BBRSim/Blackbody-Simulations && python build_hfss_project.py --spec configs/geometries/round_gap_r50um.toml"'

Expected: `Build verified: /resnick/home/rshenoy/BBRSim/projects/RoundGap.aedt` (or the same path under
`/home`), exit 0, `RoundGap.build.json` with `"verification_differences": []`, and
`RoundGap.inventory.json` listing one design `round_gap_r50um` with one object `gap` [vacuum],
bounding box [-0.05, -0.05, 0, 0.05, 0.05, 0.4] and 3 faces. PyAEDT warnings "Not enough vertices or
non-planar face" are expected: the two vertex-free end caps and the curved side have no computable normal.

## 7. Round gap, single angle at 2000 GHz

The first run (job 4077742, 2026-10-06) solved and then stopped at the rim of the disc; step 7a found why,
and the rim rule (ledger HF-029) fixed it. Its job directory `job_roundgap1` exists and holds no data, so the
runner refuses that job id again. Step 7b runs this config as `roundgap2`, beside a fine-facet variant, and its
check prints the facts below; this section keeps the single-run check, pointed at `roundgap2`:

    python - <<'PY'
    import json, pandas as pd
    job = "/home/rshenoy/BBRSim/outputs/RoundGap_r50um_2000GHz/job_roundgap2"
    d = f"{job}/RoundGap_r50um_2000GHz_Ephi=0"
    w = pd.read_csv(f"{d}/waveguide.csv"); f = pd.read_csv(f"{d}/far_field.csv")
    s = json.load(open(f"{job}/RoundGap_r50um_2000GHz.dataset.json"))
    print(len(w), len(f), w.OutgoingPower.iloc[0] / w.IngoingPower.iloc[0])
    print(s["exit_field"]["points_per_key_retained"], s["exit_field"]["cross_section"], s["modes"]["mode"], s["modes"]["propagating_count"])
    print(sum("all-nan rows" in line for line in open(f"{job}/logs/run.log")))
    PY

Expected: `7845 1369 <T>`, T strictly between 0 and 1: every lattice point inside the disc, twenty of them
exactly on the rim (rim rule, Rohan 2026-10-06, ledger HF-029: HFSS meshes the circle as a polygon, so some
inside points get no field and are written with zero field, and some just outside get field and are dropped;
`manifest.json` `exit_field.rim` counts both per key and polarization). Then the same 7845,
`{'shape': 'disc', 'radius_m': 5e-05}` (to rounding), `TE11 1`, then `0`. A non-zero last line means
`logs/run.log` holds the reader's warning that HFSS wrote unsolved lattice points as rows other than nine
tokens; under the rim rule those inside the disc are kept with zero field (counted in `exit_field.rim`) and
those outside dropped, so the run stands, but copy the warning into RESULTS.md.
Record T as the first round-gap reference value (no Windows reference exists), the passes and final
delta E from `logs/convergence_Ephi0.txt`, and `sacct -j <id> --format=MaxRSS,Elapsed`.

## 7a. Exit-lattice diagnostic (batch, debug QOS)

Step 7 (job 4077742, 2026-10-06) solved and then stopped at its first exit-field export: HFSS had written field
values at 50 lattice points outside the declared 50 um disc, which the runner refuses, and the runner deletes
each raw export after reading it. `hpc/diagnose_exit_grid.py` reruns the runner unchanged on the single-angle
config with both polarizations (two separately adapted meshes), keeps every raw exit-field export, reports the
containment check instead of raising, and maps the valued and nan lattice points against their distance from
the rim:

    cd /home/rshenoy/BBRSim/Blackbody-Simulations && git pull --ff-only fork linux-hpc-migration
    D=$(sbatch --parsable -A golwala -p expansion -q debug -N 1 -c 4 --mem=16G -t 00:30:00 -o /home/rshenoy/BBRSim/outputs/slurm-%j.out \
      --wrap 'bash -lc "source /home/rshenoy/BBRSim/bb_env.sh && cd /home/rshenoy/BBRSim/Blackbody-Simulations && python hpc/diagnose_exit_grid.py --config configs/round_gap_r50um_2000GHz_single_angle.toml --polarizations 0 1 --cores 4 --job-id diag_exitgrid1"'); D=${D%%;*}; echo "diagnostic job $D"

When it has left the queue (same shell, or put the job number in place of `$D`):

    sacct -j $D -o JobID,State,ExitCode,Elapsed,MaxRSS
    grep -E "^(containment|rim rule)" /home/rshenoy/BBRSim/outputs/slurm-$D.out
    sed -n '/^== runner exit status/,$p' /home/rshenoy/BBRSim/outputs/slurm-$D.out

Expected: one `containment` line per polarization; then per export the valued and nan lattice points inside
and outside the disc (7845 lattice points lie inside it, 20 of them on the rim), how far beyond the rim the
valued outside points lie, |E| there against the ring just inside the rim, and whether the two polarizations'
valued points inside the disc are the same set; a closing `EXIT GRID:` line. The job directory
`job_diag_exitgrid1` keeps the outside points in its tables and is not a dataset.

## 7b. Wall facets: the round gap as configured and with fine facets (batch)

Ran 2026-10-07 (jobs 4155776 and 4155777, RESULTS 7b): with 5 deg facets T came out 1.0698 times the default
for both polarizations, so the round gap's production configs now carry `wall_normal_deviation_deg = 5.0`
(Rohan, ledger HF-032) and `round_gap_r50um_2000GHz_single_angle_facets5.toml` was removed. The commands below
reproduce the step at commit `fa6669d`.

Step 7a showed HFSS's default faceting puts the wall up to 0.76 um inside the 50 um circle, and 2000 GHz is
only 14 % above the TE11 cutoff, where T may depend on the effective radius. Rohan, 2026-10-06 (ledger
HF-030): measure it before step 8. Two single-angle runs at normal incidence with both polarizations,
identical except `solver.wall_normal_deviation_deg = 5.0` (72 facets around the circle) in the second:

    cd /home/rshenoy/BBRSim/Blackbody-Simulations && git pull --ff-only fork linux-hpc-migration
    A=$(sbatch --parsable -t 02:00:00 hpc/run_frequency.sbatch configs/round_gap_r50um_2000GHz_single_angle.toml --polarizations 0 1 --job-id roundgap2); A=${A%%;*}; echo "default facets job $A"
    B=$(sbatch --parsable -t 02:00:00 hpc/run_frequency.sbatch configs/round_gap_r50um_2000GHz_single_angle_facets5.toml --polarizations 0 1 --job-id roundgap2_facets5); B=${B%%;*}; echo "fine facets job $B"
    squeue -u rshenoy

When both have left the queue (same shell, or put the job numbers in place of `$A` and `$B`):

    sacct -j $A,$B -o JobID,State,ExitCode,Elapsed,MaxRSS
    python - <<'PY'
    import json, re
    from pathlib import Path
    import pandas as pd

    root = Path("/home/rshenoy/BBRSim/outputs/RoundGap_r50um_2000GHz")
    FIELDS = ["Ex_real", "Ey_real", "Ez_real", "Ex_imag", "Ey_imag", "Ez_imag"]
    T = {}
    for name in ("roundgap2", "roundgap2_facets5"):
        job = root / f"job_{name}"
        print("==", name)
        if not (job / "manifest.json").is_file():
            log = job / "logs" / "run.log"
            lines = log.read_text(errors="replace").splitlines() if log.is_file() else ["no logs/run.log"]
            print("   no manifest; last lines of run.log:", *lines[-4:], sep="\n   ")
            continue
        m = json.loads((job / "manifest.json").read_text())
        print(f"   wall facets {m['solver']['wall_facets']}; solve {m['job']['solve_time_s']:.0f} s")
        s = json.loads((job / "RoundGap_r50um_2000GHz.dataset.json").read_text())
        short = sum("all-nan rows" in line for line in (job / "logs" / "run.log").read_text(errors="replace").splitlines())
        print(f"   sidecar: {s['exit_field']['points_per_key_retained']} per key, {s['exit_field']['cross_section']}, "
              f"{s['modes']['mode']} {s['modes']['propagating_count']}; short-row warnings {short}")
        for e in (0, 1):
            wg = pd.read_csv(job / f"RoundGap_r50um_2000GHz_Ephi={e}" / "waveguide.csv")
            T[name, e] = float(wg.OutgoingPower.iloc[0] / wg.IngoingPower.iloc[0])
            zero = int(wg[FIELDS].eq(0.0).all(axis=1).sum())
            rim = m["exit_field"]["rim"]["per_key"][str(e)][0]
            conv = (job / "logs" / f"convergence_Ephi{e}.txt").read_text(errors="replace")
            done = re.search(r"Completed\s*:\s*(\S+)", conv)
            rows = [line.strip() for line in conv.splitlines() if re.match(r"^\s*\d+\s*\|", line)]
            print(f"   Ephi={e}: T {T[name, e]:.5f}; {len(wg)} rows, {zero} zero; rim {rim['unsolved_inside_zeroed']} zeroed "
                  f"(deepest {rim['deepest_unsolved_inside_m'] * 1e9:.0f} nm), {rim['valued_outside_dropped']} dropped "
                  f"(farthest {rim['farthest_valued_outside_m'] * 1e9:.0f} nm); passes {done.group(1) if done else '?'}, "
                  f"last {rows[-1] if rows else '-'}")
        print(f"   T Ephi=1 / Ephi=0: {T[name, 1] / T[name, 0]:.5f}")
    for e in (0, 1):
        if ("roundgap2", e) in T and ("roundgap2_facets5", e) in T:
            print(f"FACETS Ephi={e}: T fine / default = {T['roundgap2_facets5', e] / T['roundgap2', e]:.5f}")
    PY

Expected: both COMPLETED; per run a `wall facets` line (None, then the `wall_facets` operation at 5 deg, which
the runner has read back from the design's mesh tree before the solve), a sidecar line `7845 per key,
{'shape': 'disc', 'radius_m': 5e-05}, TE11 1; short-row warnings 0` (a non-zero count: see step 7); per
polarization 7845 rows, as many zero rows as the rim line zeroed, at most 2.5 um deep or beyond (5 % of R),
passes under 20 with a last delta under 0.005; within each run `T Ephi=1 / Ephi=0` within about 0.5 % of 1
(normal incidence: a round gap cannot prefer a polarization). With fine facets the rim counts should fall
sharply. The `FACETS` lines are the result: a ratio within about 0.5 % of 1 means HFSS's default faceting
is good enough; a larger one means step 8 and the round-gap production configs need the fine facets.
Rohan chooses; record both runs, their memory and wall time.

## 8. Round gap, full sweep (5 deg facets) and the faceting check

Both round-gap configs carry the fine wall facets (Rohan, 2026-10-07, ledger HF-032). Alongside the sweep, a
single-angle run at 2.5 deg checks that 5 deg has converged, against step 7b's 5 deg run `roundgap2_facets5`:

    cd /home/rshenoy/BBRSim/Blackbody-Simulations && git pull --ff-only fork linux-hpc-migration
    C=$(sbatch --parsable -t 03:00:00 hpc/run_frequency.sbatch configs/round_gap_r50um_2000GHz_reference.toml --job-id reference1); C=${C%%;*}; echo "step 8 job $C"
    F=$(sbatch --parsable -t 02:00:00 hpc/run_frequency.sbatch configs/round_gap_r50um_2000GHz_single_angle_facets2p5.toml --polarizations 0 1 --job-id roundgap3_facets2p5); F=${F%%;*}; echo "2.5 deg check job $F"
    squeue -u rshenoy

When both have left the queue (same shell, or put the job numbers in place of `$C` and `$F`):

    sacct -j $C,$F -o JobID,State,ExitCode,Elapsed,MaxRSS
    python - <<'PY'
    import json, re
    from pathlib import Path
    import pandas as pd

    root = Path("/home/rshenoy/BBRSim/outputs/RoundGap_r50um_2000GHz")
    KEYS = ["IWavePhi", "IWaveTheta"]


    def converged(job, e):
        text = (job / "logs" / f"convergence_Ephi{e}.txt").read_text(errors="replace")
        done = re.search(r"Completed\s*:\s*(\S+)", text)
        rows = [line.strip() for line in text.splitlines() if re.match(r"^\s*\d+\s*\|", line)]
        return f"passes {done.group(1) if done else '?'}, last {rows[-1] if rows else '-'}"


    def manifest(job):
        if (job / "manifest.json").is_file():
            return json.loads((job / "manifest.json").read_text())
        log = job / "logs" / "run.log"
        lines = log.read_text(errors="replace").splitlines() if log.is_file() else ["no logs/run.log"]
        print("   no manifest; last lines of run.log:", *lines[-4:], sep="\n   ")
        return None


    job = root / "job_reference1"
    print("== step 8:", job.name)
    m = manifest(job)
    if m:
        s = json.loads((job / "RoundGap_r50um_2000GHz.dataset.json").read_text())
        print(f"   wall facets {m['solver']['wall_facets']}; solve {m['job']['solve_time_s']:.0f} s")
        print(f"   sidecar: {s['exit_field']['points_per_key_retained']} per key, {s['exit_field']['cross_section']}, "
              f"{s['modes']['mode']} {s['modes']['propagating_count']}")
        six = []
        for e in (0, 1):
            wg = pd.read_csv(job / f"RoundGap_r50um_2000GHz_Ephi={e}" / "waveguide.csv")
            k = wg.drop_duplicates(KEYS).set_index(KEYS)
            t = k.OutgoingPower / k.IngoingPower
            vals = [float(t[(p, 180.0)]) for p in (0.0, 45.0, 90.0)]
            six += vals
            rim = m["exit_field"]["rim"]["per_key"][str(e)]
            print(f"   Ephi={e}: {len(wg)} rows, {len(k)} keys x {sorted(set(wg.groupby(KEYS).size()))}; "
                  f"T at theta 180, phi 0/45/90: {[round(v, 5) for v in vals]}; rim zeroed <= "
                  f"{max(r['unsolved_inside_zeroed'] for r in rim)}, dropped <= {max(r['valued_outside_dropped'] for r in rim)}; "
                  f"{converged(job, e)}")
        print(f"   spread of the six normal-incidence T: {max(six) / min(six) - 1:.3%}")

    print("== faceting check: 2.5 deg against 5 deg (step 7b's roundgap2_facets5)")
    T = {}
    for name in ("roundgap2_facets5", "roundgap3_facets2p5"):
        job = root / f"job_{name}"
        m = manifest(job)
        if not m:
            continue
        for e in (0, 1):
            wg = pd.read_csv(job / f"RoundGap_r50um_2000GHz_Ephi={e}" / "waveguide.csv")
            T[name, e] = float(wg.OutgoingPower.iloc[0] / wg.IngoingPower.iloc[0])
        print(f"   {name}: {m['solver']['wall_facets']['normal_deviation_deg']} deg, T {T[name, 0]:.5f} / {T[name, 1]:.5f}, "
              f"solve {m['job']['solve_time_s']:.0f} s; {converged(job, 0)}")
    for e in (0, 1):
        if ("roundgap2_facets5", e) in T and ("roundgap3_facets2p5", e) in T:
            print(f"FACETS Ephi={e}: T 2.5 deg / 5 deg = {T['roundgap3_facets2p5', e] / T['roundgap2_facets5', e]:.5f}")
    PY

Expected: both COMPLETED. Step 8: the 5 deg `wall_facets` operation; the sidecar line `7845 per key, disc, TE11
1`; per polarization 15 x 7845 = 117675 rows, 15 keys of 7845; the six T values at normal incidence (three
azimuths, two polarizations) near step 7b's 0.5096 and equal to within about 0.5 %: at normal incidence a round
gap cannot prefer a polarization or an azimuth, and step 7b's polarizations agreed to 0.02 %. A large spread
means a frame or polarization error; stop and report. The `FACETS` lines: 2.5 deg against 5 deg within about
0.5 % of 1 means 5 deg has converged and step 8 is the round gap's first dataset; a larger change means step 8
reruns at 2.5 deg. Record memory and wall time; they size the job array.

## 8b. Power diagnostic on the solved projects (batch, no solve)

The power normalization is under review (Codex queue CDX-009, 013, 014, 015; ledger G4-020, OR-019). This step
measures what the tables' T is made of, on copies of the three solved projects: the exit integrals of |Re S|,
Re S_z, |Re S_z| and Re S . n; the entrance flux with total and with scattered fields; the incident amplitude
|E_total - E_scattered| just inside the entrance; and the far-field power over both hemispheres, with total and
with scattered fields. No solve and no solver licence, one AEDT session per job. Each job gets
`diagnose_poynting/poynting.json`; nothing else in the job changes, and the copy is removed after a successful
run. One paste:

    cd /home/rshenoy/BBRSim/Blackbody-Simulations && git pull --ff-only fork linux-hpc-migration
    O=/home/rshenoy/BBRSim/outputs
    for J in $O/InfParallelPlate_crack1Rohan_500GHz/job_tight1 $O/InfParallelPlate_crack2_500GHz/job_tight1 $O/RoundGap_r50um_2000GHz/job_reference1; do
      du -sh $J/project/*.aedtresults
      sbatch -A golwala -p expansion -q normal -N 1 -c 4 --mem=16G -t 02:00:00 -o $O/slurm-%j.out \
        --wrap "bash -lc 'source /home/rshenoy/BBRSim/bb_env.sh && cd /home/rshenoy/BBRSim/Blackbody-Simulations && python hpc/diagnose_poynting.py --job $J'"
    done
    squeue -u rshenoy

When the three jobs have left the queue, for each job id:

    sed -n '/== diagnose_poynting/,/POYNTING:/p' /home/rshenoy/BBRSim/outputs/slurm-<id>.out
    sacct -j <id> -o JobID,State,ExitCode,Elapsed,MaxRSS

Expected: exit 0 and `POYNTING: ok` for each, which means the exit |Re S| integral reproduced the recorded
OutgoingPower at every key (the copy, the variation and the expression are the run's). `out/in`, the exit flux
over the entrance flux with total fields, is 1.00 at every transmitting key to the mesh's accuracy (about 1 %),
since the walls are PEC. `|E_inc|` is 1.000 V/m at every key, as IngoingPower assumes. `|S|/S_z` is 1.000 at
normal incidence. Measured rather than predicted: `|S|/S_z` at the oblique keys (CDX-009); `in/geom`, the power
entering per unit of geometric incident flux I A |cos alpha| (above 1 means capture beyond the geometric
opening; blank at 90 deg); `scat/geom`; and `ff/flux`, `ff_fwd/flux`, `ff_scat/flux`, the far field against the
exit flux (Codex found 0.18 and 0.34 on the stored outward hemisphere, CDX-015). Exit 2 with "nothing started"
names its reason; exit 1 lists what failed in the table's header and keeps the copy under
`diagnose_poynting/project` for inspection (remove it after reporting). Paste the three blocks and the `sacct`
lines.

## 8c. Power diagnostic: the cracks again, and the round gap's exit code (batch, no solve)

Step 8b's crack jobs failed only on the far field, which the shared `.ffd` reader refused (RESULTS 8b; fixed in
`9fccac4`). This reruns the two cracks and reads why the round-gap job exited 1 after printing `ok`. Each run 1
record stays as `diagnose_poynting_run1/poynting.json`; its project copy and exports (the size of the job's
`.aedtresults`) are removed first, and the jobs' own `project/` directories are not touched. One paste:

    cd /home/rshenoy/BBRSim/Blackbody-Simulations && git pull --ff-only fork linux-hpc-migration
    O=/home/rshenoy/BBRSim/outputs
    tail -15 $O/slurm-4264788.out
    R=$O/RoundGap_r50um_2000GHz/job_reference1/diagnose_poynting; ls $R; rm -rf "$R/project" "$R/scratch"
    for J in $O/InfParallelPlate_crack1Rohan_500GHz/job_tight1 $O/InfParallelPlate_crack2_500GHz/job_tight1; do
      D=$J/diagnose_poynting
      test -f "$D/poynting.json" && rm -rf "$D/project" "$D/scratch" && mv "$D" "${D}_run1" && ls "${D}_run1"
      sbatch -A golwala -p expansion -q normal -N 1 -c 4 --mem=16G -t 02:00:00 -o $O/slurm-%j.out \
        --wrap "bash -lc 'source /home/rshenoy/BBRSim/bb_env.sh && cd /home/rshenoy/BBRSim/Blackbody-Simulations && python hpc/diagnose_poynting.py --job $J'"
    done
    squeue -u rshenoy

When both jobs have left the queue, for each job id:

    sed -n '/== diagnose_poynting/,/POYNTING:/p' /home/rshenoy/BBRSim/outputs/slurm-<id>.out
    tail -3 /home/rshenoy/BBRSim/outputs/slurm-<id>.out
    sacct -j <id> -o JobID,State,ExitCode,Elapsed,MaxRSS

Expected: the round gap's log ends in a traceback from `shutil.rmtree` (an `OSError` such as "Directory not
empty" or "Device or resource busy", from AEDT's `.nfs` files) after its `POYNTING: ok` line, and `ls` shows
`project` beside `poynting.json`; anything else, stop and report. Each `ls ..._run1` shows only
`poynting.json`. Both crack jobs: exit 0 and `POYNTING: ok`; the T, `|S|/S_z`, `out/in`, `in/geom` and
`|E_inc|` columns equal run 1's (same project, same expressions); `ff/flux`, `ff_fwd/flux` and `ff_scat/flux`
filled at every transmitting key and blank at the keys that carry no power. Predicted, not assumed: `ff/flux`
near 0.2 for crack1 and 0.4 for crack2 at normal incidence, the fraction of the exit flux that equivalent
currents on a slot about a tenth of a wavelength wide radiate (Codex found 0.18 and 0.34 on the outward
hemisphere, CDX-015). If a copy cannot be removed, the run now says so and still exits 0; remove it by hand.
Paste the round gap's `tail` and `ls`, the two blocks with their `tail` lines, and the `sacct` lines.

## 8d. OR-035 check on the solved projects (batch, no solve)

OR-035 changed two things the fake AEDT models but real AEDT has not yet shown. HF-045: a rectangle whose extent is
not a whole number of exit-lattice steps gets the fitted step, with the ExportOnGrid stop padded by half a step, so
both walls are nodes. CDX-023: exit points must lie on the exit plane x_e = 0 to 1e-12 m, and lattice bounds off it
stop a run before its solve. `hpc/check_or035.sbatch` checks both on copies of the three solved projects, with no
solve and one AEDT seat per project. (a) On crack1's copy, the runner's own lattice and export over the exit face's
y lattice and z in [-0.01875, 0.01875] mm at the configured 1 um: 37.5 steps, so `fit_axis` gives 38 intervals of
0.0375/38 mm and 39 nodes. (b) On all three copies, the exit-face outline in the exit frame, from the runner's own
face description. It runs inside Rohan's one HPC batch on its own (B4), independent of step 9. crack1's copy carries
its `.aedtresults` (about 2.5 GB) under `<job>/check_or035/` while the job runs; with step 9b's re-exports (about
7.5 GB) running at the same time, the peak on `/home` is about 10 GB, so check `quota -s` or `df -h /home/rshenoy`
first. One paste:

    cd /home/rshenoy/BBRSim/Blackbody-Simulations && git pull --ff-only fork linux-hpc-migration
    mkdir -p /home/rshenoy/BBRSim/outputs
    B4=$(sbatch --parsable hpc/check_or035.sbatch); B4=${B4%%;*}; echo "B4 job $B4"

When it has left the queue (same shell, or put the job number in place of `$B4`):

    sacct -j $B4 -o JobID,State,ExitCode,Elapsed,MaxRSS
    grep -E '^(== |   \((a|b)\)|       |   HF-045|   copies|   note)|: stopped: |Error:|Exception:' /home/rshenoy/BBRSim/outputs/check_or035-$B4.out

Expected: (a) the fitted step 0.0375/38 mm (0.000986842105263158) with that spacing string in the ExportOnGrid
call; `expect 101 x 39 = 3939 rows per key, z end nodes at -1.875e-05 and 1.875e-05 m to 1e-12 m`; 15 keys exported
with `rows per key [3939]`; `HF-045: PASS`. (b) `max |x|` at most 1e-09 mm and PASS for crack1 and crack2 (4 vertices
each) and the round gap (edge samples). The last line `== RESULT: PASS` with exit 0, and the copies removed. A FAIL
in (a) means AEDT did not put a node on both walls, or gave another count, despite the padded stop; the runner's
refusal is printed and the raw export of that key is kept in `<crack1 job>/check_or035/scratch/`. A FAIL in (b)
names the project whose exit face lies off the plane. After a failure the copies are kept; remove each
`<job>/check_or035/` by hand before a rerun (exit 2, "nothing started", if one exists). Paste the `grep` output and
the `sacct` lines.

## 9. Schema-2 re-export of the three solved jobs, no solve (batch)

The schema-2 tables (HF-050) for the three solved 2025 R2 jobs, crack1 and crack2 `job_tight1` (step 5d) and the round
gap's `job_reference1` (step 8), come from a fresh extraction on copies of their solved projects: no solve, one AEDT
seat per project, no solver licence. Each fresh job is accepted against its source, then the `_v2` trees are
assembled for BBRsim's schema-2 reader. These steps run inside Rohan's one HPC batch, in chain B: B1 (9b-9c) -> B2
(9d) -> B3 (9e); B4 (step 8d) alone; chain B never waits on chain C, and a failure in one chain never blocks another.
The pastes below are what the batch submits; each also works on its own, in this order.

With schema 2 the default (HF-050): a run that is not crack1 at 500 GHz needs `--column-order-reference`, the file
crack1's job keeps as `scratch/column_order_reference.json` (HF-041); a `--polarizations 1` run is refused; the
single-angle pastes (steps 4, 7, 7b, and step 8's 2.5 deg check) still run schema 1, as their configs say
(`output.schema = 1`).

### 9a. Pull the fork at the tip

    cd /home/rshenoy/BBRSim/Blackbody-Simulations && git pull --ff-only fork linux-hpc-migration
    git rev-parse HEAD
    quota -s; df -h /home/rshenoy

Expected: `git rev-parse HEAD` prints the tip pushed for the batch. Each re-export copies its source's
`.aedtresults` (about 2.5 GB) onto `/home` NFS, about 7.5 GB for the three (B1). Step 8d (B4) is an independent
chain that may run at the same time and holds its own copy of crack1's `.aedtresults` (about 2.5 GB) under
`check_or035/` while it runs, so the peak on `/home` is about 10 GB: `quota -s` or `df -h /home/rshenoy` must leave
that much free.

### 9b. Re-export, with its acceptance (9c) in the same job (batch B1)

`hpc/reexport_v2.sbatch` checks before any AEDT session that each source job and its step-8b record exist, that the
record was made on that source (its `manifest_sha256`) and covers every key and both polarizations
(`bbsim.reexport_check --preflight`), and that no fresh job exists yet. It then re-exports crack1 into
`job_tight1_v2`, with no reference (it measures its own), and stops if that fails; then crack2 into `job_tight1_v2`
and the round gap into `job_reference1_v2`, both with crack1's reference, a failure of one not stopping the other.
The records it compares against are step 8c's `<job_tight1>/diagnose_poynting/poynting.json` for the cracks and step
8b's `<job_reference1>/diagnose_poynting/poynting.json` for the round gap. One paste:

    mkdir -p /home/rshenoy/BBRSim/outputs
    B1=$(sbatch --parsable hpc/reexport_v2.sbatch); B1=${B1%%;*}; echo "B1 job $B1"

When it has left the queue (same shell, or put the job number in place of `$B1`):

    sacct -j $B1 -o JobID,State,ExitCode,Elapsed,MaxRSS
    grep -E '^(== |PREFLIGHT|preflight: |REEXPORT|report: |(PASS|FAIL) |   (manifest|Ephi=|sidecar|ff_over_flux|reference|column-order|\(no summary|\(the reference))|wrote the schema-2 sidecar|done; manifest at|failed|refused|earlier attempt|reexport_check: |bbsim.run: aedt: |Error:|Exception:' /home/rshenoy/BBRSim/outputs/reexport_v2-$B1.out
    sed -n '/^largest deviations:/,/^REEXPORT/p' /home/rshenoy/BBRSim/outputs/reexport_v2-$B1.out

Expected: three `PREFLIGHT: PASS` lines. Each re-export's log ends with the schema-2 sidecar written and
`done; manifest at ...`, followed by its summary: `qualified true` for both polarizations, the achieved delta
against the source's target 0.005 (HF-046: crack1 0.0036596 / 0.0037117, crack2 0.004126 / 0.00434, the round gap
about 0.0025 / 0.0029), `solve_time_s null` (no solve, no `analyze`) and the source as `reused_solution`. crack1's two
column-order cases print `passed true` with at least one `discriminating true`, and their P2/P (the legacy tables
gave 0.713 / 0.725; a value far from them is worth a look, not a failure). `ff_over_flux` at normal incidence is near
0.179 for crack1 and 0.348 for crack2 (step 8c's outward hemisphere over the exit flux), and near 0.66 for the round
gap (step 8b's forward fraction).

9c, inside the same job: `REEXPORT: PASS` three times, each table with its largest deviations printed, its
`powers.poynting` row passing against the step-8b record, and the last line `== RESULT: PASS` with exit 0. What
this run shows for the first time, beyond PASS: re-defining both named expressions (`ClearAllNamedExpr`, then
`outgoing_power` and `exit_flux`), `export_convergence` and `EditSources(TotalFields)` work on a reopened solved
design (an AEDT failure there names itself in the log before anything is written); the `convergence` row passes,
so the reopened design's convergence tables equal the source's; and the dark keys that carry only noise (theta 0,
T about 1e-44) are judged by the `exit_field.dark` row (W = 1 and the source values kept within eps), not by
`exit_field.values`, which covers the records with a fitted W. Each table is printed before its
`reexport_check.json` is written into the fresh job, on a `report:` line; a report that cannot be written (for
example a full `/home`) fails that dataset on that line.

Anything else ends `== RESULT: FAIL <what failed>`: exit 2 with "nothing started" names the preflight's reason;
exit 1 names the failed re-export or acceptance (a failed crack1 re-export means the other two did not start).
Which job may be removed: a failed `job_*_v2` holds only a copy of its source's solved project, re-creatable by
rerunning the re-export, so remove it by hand before a rerun (the runner refuses an existing one, and the script
never deletes anything). Never remove the sources `job_tight1` and `job_reference1`: their solves are the
originals. Paste the `grep` and `sed` output and the `sacct` lines.

### 9d. Assemble the two `_v2` trees (batch B2)

`hpc/assemble_tree.sbatch <tree_dir> <expected_files> <job_dir>...` gates every job before copying anything:
`solver.qualified` true in its manifest (OR-027: assembly takes qualified runs only), printed with one line per
polarization; one sidecar named by its dataset; both `_Ephi` directories with `waveguide.csv`, `far_field.csv` and
`manifest.json`; the acceptance verdicts (HF-052, CDX-026), each printed as an `acceptance:` line: a re-export
(its sidecar records `provenance.reexport`) must hold `reexport_check.json` with `passed` true (step 9c), and a job
holding `grid_acceptance.json` (step 10b) must record `passed` true there, so the tree never rests on the Slurm
dependency alone; no dataset twice; one schema major per tree; the expected file count. It refuses an existing
non-empty tree directory. It lays the tree out as steps 5e and 8a did, and writes and checks `SHA256SUMS`. The
schema-1 trees `cracks_500GHz_2025r2` and `testworld_2025r2` stay as they are. One paste (after 9b; the
dependency waits for it; same shell as 9b, or put B1's job number in place of `$B1`):

    T=/resnick/groups/golwala/rshenoy/bbsim/trees; O=/home/rshenoy/BBRSim/outputs
    C1=$O/InfParallelPlate_crack1Rohan_500GHz/job_tight1_v2; C2=$O/InfParallelPlate_crack2_500GHz/job_tight1_v2
    R=$O/RoundGap_r50um_2000GHz/job_reference1_v2
    A=$(sbatch --parsable --dependency=afterok:$B1 --kill-on-invalid-dep=yes hpc/assemble_tree.sbatch $T/cracks_500GHz_2025r2_v2 14 $C1 $C2); A=${A%%;*}
    W=$(sbatch --parsable --dependency=afterok:$B1 --kill-on-invalid-dep=yes hpc/assemble_tree.sbatch $T/testworld_2025r2_v2 21 $C1 $C2 $R); W=${W%%;*}
    echo "B2 jobs $A $W"

When both have left the queue:

    sacct -j $A,$W -o JobID,State,ExitCode,Elapsed
    for j in $A $W; do grep -v ': OK$' /home/rshenoy/BBRSim/outputs/assemble_tree-$j.out; done

Expected, per tree: every job's two polarization lines read `qualified true` with the achieved delta against 0.005;
every dataset `schema_version 2.0`; every job's `acceptance: reexport_check.json passed true (REEXPORT: PASS)`;
`copied 14` (cracks) and `copied 21` (test world) files;
`SHA256SUMS: 14 files, all OK` and `21 files`; the last line `== RESULT: PASS` with exit 0. A refusal names the job
and its reason, and nothing is copied; a tree directory that exists and is not empty is never touched. If B1 failed,
Slurm cancels both (`--kill-on-invalid-dep`).

### 9e. BBRsim's tree check on the `_v2` trees (BBRsim side, batch B3)

BBRsim's `g4_tree_check.sbatch`, after its full check, with the BBRSimulation checkout at the schema-2 reader
(`feat/schema-2-reader`). That branch carries OR-024's BBR028 whitelist, which `testworld_2025r2_v2` needs: the
round gap's 2000 GHz set lies above its 1757 GHz polarization-filter limit, so BBR028 fires once. Its header
states the expected results; record them with the BBRsim side's. The paste, with `<A1>` the full check's job:

    sbatch --dependency=afterok:<A1>:$A:$W /home/rshenoy/BBRSim/g4_tree_check.sbatch $T/cracks_500GHz_2025r2_v2 $T/testworld_2025r2_v2

## 10. The first 35-key solves (OR-036; licensed solves, Rohan's separate go)

Licensed solves: run only on Rohan's separate go (OR-036). Each solve holds an AEDT seat and an HFSS solver
increment for its whole run; the three run one after another.

OR-036's incidence grid (Rohan, 2026-10-19) is a set of incoming directions in the crack-local frame: theta 90 to
180 deg in 15 deg steps and phi 0 to 90 deg in 22.5 deg steps, 35 keys, all reachable. Its first solves are crack1
and crack2 at 500 GHz and the round gap at 2000 GHz, from `configs/crack1_500GHz_grid35.toml`,
`configs/crack2_500GHz_grid35.toml` and `configs/round_gap_r50um_2000GHz_grid35.toml`: each is its 15-key source
(`*_reference_tight.toml` for the cracks, `round_gap_r50um_2000GHz_reference.toml`) with only the six incidence lines
changed and `schema = 2` in `[output]`, so MaxDeltaE 0.005, 20 passes, and for the round gap 5 deg wall facets and
2000 GHz. They run in Rohan's one HPC batch as chain C: C1 (10a) -> C2 (10b) -> C3 (10c) -> C4 (10d); chain C never
waits on chain B, and a failure in one chain never blocks another.

With schema 2 the default (HF-050), every job but the crack1 500 GHz reference needs `--column-order-reference`
(HF-041). crack1's 35-key solve is itself that reference job (crack1 at 500 GHz, a grid holding (0, 180) and
(90, 180) deg), so it runs with none and keeps its own as `job_grid35_1/scratch/column_order_reference.json`;
crack2 and the round gap take that file, not step 9's, so chain C does not wait for chain B.

Space on `/home`: each job keeps its solved project and `.aedtresults`, at least the 15-key runs' sizes, about
2.8 / 4.8 / 2.5 GB for crack1 / crack2 / the round gap. Field solutions are stored per excitation, so a 35-key
`.aedtresults` may approach 35/15 of the 15-key one, roughly 24 GB for the three. Not measured yet: record it with
`du -sh` (it sizes the job array). With chain B's copies (about 10 GB at its peak, step 9a) on `/home` at the same
time, check `quota -s` or `df -h /home/rshenoy` first.

### 10a. Solve the three (batch C1)

`hpc/solve_grid35.sbatch <config> [extra runner args]` runs the runner with `--job-id grid35_1`, so each job is
`<dataset>_<label>/job_grid35_1`. Before any AEDT session it checks that the config declares the 35 keys, that the
runner's own pre-AEDT check passes, that a `--column-order-reference` file exists, and that the job does not exist
yet; a refusal ends `== RESULT: FAIL nothing started: <why>` with exit 2. The three are chained with `afterany`, so
crack2 and the round gap start even if crack1 failed; without crack1's reference they stop there, before AEDT (the
runner refuses a missing reference before the job starts too). One paste:

    cd /home/rshenoy/BBRSim/Blackbody-Simulations && git pull --ff-only fork linux-hpc-migration
    mkdir -p /home/rshenoy/BBRSim/outputs; quota -s; df -h /home/rshenoy
    O=/home/rshenoy/BBRSim/outputs
    REF=$O/InfParallelPlate_crack1Rohan_500GHz/job_grid35_1/scratch/column_order_reference.json
    G1=$(sbatch --parsable hpc/solve_grid35.sbatch configs/crack1_500GHz_grid35.toml); G1=${G1%%;*}
    G2=$(sbatch --parsable --dependency=afterany:$G1 hpc/solve_grid35.sbatch configs/crack2_500GHz_grid35.toml --column-order-reference $REF); G2=${G2%%;*}
    G3=$(sbatch --parsable --dependency=afterany:$G2 hpc/solve_grid35.sbatch configs/round_gap_r50um_2000GHz_grid35.toml --column-order-reference $REF); G3=${G3%%;*}
    echo "C1 jobs $G1 $G2 $G3"

When all three have left the queue (same shell, or put the job numbers in place of `$G1`, `$G2`, `$G3`):

    sacct -j $G1,$G2,$G3 -o JobID,State,ExitCode,Elapsed,MaxRSS
    for j in $G1 $G2 $G3; do grep -E '^(== |host=|preflight: |   (manifest|Ephi=|sidecar|cell_90|column-order|reference passed))|wrote the schema-2 sidecar|solve finished|convergence target|done; manifest at|job failed|refused|bbsim.run: aedt: |Error:|Exception:' /home/rshenoy/BBRSim/outputs/solve_grid35-$j.out; done
    du -sh /home/rshenoy/BBRSim/outputs/*/job_grid35_1 /home/rshenoy/BBRSim/outputs/*/job_grid35_1/project/*.aedtresults

Expected, for each: the `host=` line with the commit that ran; `preflight: ok; 5 x 7 = 35 keys, schema 2`; the run
logs `wrote the schema-2 sidecar` and ends `done; manifest at .../job_grid35_1/manifest.json` and
`run_hfss_frequency.py exited 0`, with no `convergence target ... not met` warning; the summary shows
`qualified true` for both polarizations with the achieved delta at most the target 0.005, the sidecar with
`schema_version 2.0; 35 keys` and `cell_90 0.036633054690581604` (0.0366330547: the band 82.5 to 90 deg, HF-049's
rule on this grid); crack1's two column-order cases `passed true`, at least one `discriminating true`, P2/P printed
(the legacy tables gave 0.713 / 0.725; printed, not pinned); the last line `== RESULT: PASS` with exit 0. Cost, from
the 15-key runs times 35/15 (time grows with the keys, memory does not): crack1 about 20 min and 2.8 GB (15 keys:
8.7 min, 2.8 GB), crack2 about 34 min and 4.8 GB (14.4 min, 4.8 GB), the round gap about 10 min and 2.5 GB (4 min
16 s, 2.5 GB). Record `Elapsed`, `MaxRSS` and the `du -sh` sizes: they size the job array. Anything else ends
`== RESULT: FAIL <reasons>` with exit 1 (the runner's exit code, `solver.qualified`, the key count, cell_90,
crack1's reference), or exit 2 with "nothing started". The script never deletes anything. A failed job whose
`logs/run.log` shows the solve finished (`solve finished in ... s`; its FAIL line then ends by saying the solved
project is kept) holds a licensed solve that a later re-extraction can reuse: do not delete it. Move it aside
(`mv job_grid35_1 job_grid35_1_failed`), keep it and report the failure; rerun only on purpose, because a rerun
solves again. Paste the `grep` and `du` output and the `sacct` lines.

### 10b. Key-by-key acceptance against the 15-key solves (batch C2)

`hpc/grid_acceptance.sbatch` runs `python -m bbsim.grid_acceptance --new <job_grid35_1> --reference <15-key job>` on
the three pairs: crack1 and crack2 against `job_tight1` (step 5d), the round gap against `job_reference1` (step 8).
Each compares, at the nine keys both grids declare (phi 0, 45, 90 x theta 90, 135, 180 deg), the new job's
magnitude-power ratio, its sidecar's `outgoing_power_abs_re_s_w` over IngoingPower (A I), with the 15-key job's
OutgoingPower / IngoingPower. A key passes within 0.5 % of T_ref; one whose T_ref is below 1e-3 of its
polarization's largest is reported as noise (with a note if its T_new is not also below that band). Neither the
stored schema-2 flux nor c_k enters: the flux differs from |Re S| by 0.2 % at normal incidence and 10 to 18 % on
oblique keys (step 8b), and c_k at theta 90 deg follows the grid, so the comparison isolates mesh and solver
differences from the decided numerator and divisor changes. 0.5 % sits above the mesh noise at MaxDeltaE 0.005
(0.1 % between the separately meshed polarizations, step 5d; 0.3 % for the round gap's facet check, step 8). It also
requires the new job to declare schema 2.x and exactly the 35-key grid, the reference schema 1.x and all nine keys,
the same dataset and frequency label, and IngoingPower constant within each table and equal between the jobs to
1e-9; a non-finite number fails the row that reads it. Each pair prints its table, then writes
`<job_grid35_1>/grid_acceptance.json`; a report that cannot be written fails that pair. A failing pair does not stop
the others, and each is judged by its exit code. One paste, after 10a (same shell):

    GACC=$(sbatch --parsable --dependency=afterok:$G1:$G2:$G3 --kill-on-invalid-dep=yes hpc/grid_acceptance.sbatch); GACC=${GACC%%;*}
    echo "C2 job $GACC"

When it has left the queue:

    sacct -j $GACC -o JobID,State,ExitCode,Elapsed
    cat /home/rshenoy/BBRSim/outputs/grid_acceptance-$GACC.out

Expected: `GRID35: PASS` three times at 0.5 % on the nine shared keys, each table with its 18 rows (nine keys, two
polarizations), and the last line `== RESULT: PASS` with exit 0. Noise rows, as guidance and not a pin: on the
legacy 500 GHz tables 6 per crack, the polarization the slot cuts off (Ephi=0 at phi 90 deg and Ephi=1 at phi 0 deg,
three theta each, T about 1e-10); the round gap's count is whatever its reference gives. A pair that fails names its
key, polarization and relative difference, and the last line `== RESULT: FAIL <datasets>` names it with its exit
code (exit 1): the new solve does not reproduce the old one there. Paste the tables; do not loosen the tolerance
silently. If a solve failed, Slurm cancels this job (`--kill-on-invalid-dep`).

### 10c. Assemble the two 35-key trees (batch C3)

The same `hpc/assemble_tree.sbatch` as step 9d, with the same gate (qualified runs only, the acceptance verdicts,
one schema major per tree, the expected file count; here each job's `grid_acceptance.json` from step 10b must
record `passed` true), into two new trees beside the schema-1 and `_v2` ones: one grid per dataset per tree
(OR-036), which BBRsim's schema-2 reader enforces once G4-032 lands (one grid and one schema major per dataset ID).
One paste (same shell as 10b):

    T=/resnick/groups/golwala/rshenoy/bbsim/trees; O=/home/rshenoy/BBRSim/outputs
    N1=$O/InfParallelPlate_crack1Rohan_500GHz/job_grid35_1; N2=$O/InfParallelPlate_crack2_500GHz/job_grid35_1
    NR=$O/RoundGap_r50um_2000GHz/job_grid35_1
    GT1=$(sbatch --parsable --dependency=afterok:$GACC --kill-on-invalid-dep=yes hpc/assemble_tree.sbatch $T/cracks_500GHz_grid35 14 $N1 $N2); GT1=${GT1%%;*}
    GT2=$(sbatch --parsable --dependency=afterok:$GACC --kill-on-invalid-dep=yes hpc/assemble_tree.sbatch $T/testworld_grid35 21 $N1 $N2 $NR); GT2=${GT2%%;*}
    echo "C3 jobs $GT1 $GT2"

When both have left the queue:

    sacct -j $GT1,$GT2 -o JobID,State,ExitCode,Elapsed
    for j in $GT1 $GT2; do grep -v ': OK$' /home/rshenoy/BBRSim/outputs/assemble_tree-$j.out; done

Expected, as step 9d: every job's two polarization lines `qualified true` against 0.005; every dataset
`schema_version 2.0`; every job's `acceptance: grid_acceptance.json passed true (GRID35: PASS)`; `copied 14`
(cracks) and `copied 21` (test world); `SHA256SUMS: 14 files, all OK` and
`21 files`; the last line `== RESULT: PASS` with exit 0. A refusal names the job and its reason, and nothing is
copied. If the acceptance failed, Slurm cancels both (`--kill-on-invalid-dep`).

### 10d. BBRsim's tree check on the 35-key trees (BBRsim side, batch C4)

BBRsim's `g4_tree_check.sbatch`, after its full check, with the BBRSimulation checkout at the schema-2 reader
(`feat/schema-2-reader`, as step 9e; BBRsim `main` refuses a schema-2 dataset with BBR024 until that reader is
merged). Its header states the expected results; record them with the BBRsim side's. The paste, with `<A1>` the full
check's job:

    sbatch --dependency=afterok:<A1>:$GT1:$GT2 /home/rshenoy/BBRSim/g4_tree_check.sbatch $T/cracks_500GHz_grid35 $T/testworld_grid35

Expected: as step 9e (batch B3), with 35-key manifests: checksums 14 / 21 OK, BBRsim's sidecar check PASS on schema
2, every row PASS.
