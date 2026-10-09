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
