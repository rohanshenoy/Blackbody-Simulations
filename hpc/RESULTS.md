# HPC verification log

Fill each section with the exact command, date, Slurm job ID, and pasted output.
A section without pasted output is not done. Steps are defined in `hpc/README.md`.

## Status (2026-10-05)

- Steps 0 to 3: PASS (2026-09-28). The step 1 smoke test also passes under the
  merged `bb_env.sh` that loads Geant4 (Slurm job 3693291, 2026-09-30).
- Next: step 4. First update the HPC checkout to the fork tip and capture
  `requirements-hpc.txt` with `python -m pip freeze --exclude-editable`.
- After that: step 5, then the round-gap steps 5b, 6, 7 and 8 (code for them is on
  the branch and tested on the Mac against a fake AEDT only).
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

2026-10-05, batch (`hpc/run_frequency.sbatch`: 4 cores, 32G, `-t 02:00:00`), code at 271c280, job
`baseline1`. Result: FAIL after the solve.

    13:09:37 solving 500GHz with 4 cores (1 incident angles x 2 polarizations)
    13:10:44 Design setup None solved correctly in 0.0h 1.0m 6.0s      (PyAEDT: AnalyzeAll, blocking)
    13:10:44 solve finished in 66.6 s
    GrpcApiError: Failed to execute gRPC AEDT command: ClcEval        (extract.py, evaluate_outgoing_power)

Diagnostic (`hpc/diagnose_calculator.py` at 38a17f4, 13:17, on a copy of the job's project and results):
the three designs, the setup `500GHz`, the solutions `500GHz : LastAdaptive` and the Ephi table, the
variable `Ephi='0'` and the excitation `plane_wave_500GHz` (PhiPoints 1, ThetaPoints 1) are all present.
Every `ClcEval` failed with AEDT's message "Script macro error: Error in performing operation", a bare
constant evaluated with `Freq` alone included, so neither the expression nor the variation arguments are
at fault. (`odesign.ListVariations` is not reachable over gRPC; the call belongs to the Solutions module.)

Root cause: `ClcEval` followed by `GetTopEntryValue` reads the GUI calculator stack back. That stateful
round trip exists over COM (Windows, where the legacy script ran) but not over the gRPC transport PyAEDT
uses on Linux. PyAEDT 1.7.0 evaluates calculator expressions by `CalculatorWrite` to a file and reading
it back (`visualization/post/fields_calculator.py`: "ClcEval does not return any value"), and pyEPR's
notes on its PyAEDT gRPC backend name this read-back as the one operation that does not survive gRPC.

Fix (this branch, with this entry): `evaluate_outgoing_power` writes each value with `CalculatorWrite`
(the intrinsics plus `Phase='0deg'`, as PyAEDT passes them) and reads it with `read_calculator_scalar`.
The fake AEDT used by the tests now refuses `ClcEval` and `GetTopEntryValue` the way gRPC does.

Step 4a (2026-10-06, job 4064670, hpc-22-18, code at `893f4f6`), the extraction diagnostic on a copy
of `job_baseline1`: **the solve produced no solution.** The results folder holds 20 files, 0.0 MB
(bookkeeping `.asol` files, geometry caches, one `opti906_0.profile` of 1 KB written at 13:10:42, no
mesh, no fields); `Solutions.GetAvailableVariations("500GHz : LastAdaptive")` is empty; the Ephi=0
convergence table reads `Completed : N/A` with no pass rows; the profile is empty ("There is no profile
data to be exported"). `AnalyzeAll(blocking)` returned after 66 s without having run an adaptive pass,
and PyAEDT logged "solved correctly" from the elapsed time alone. The setup (`500GHz`, MaxDeltaE 0.02,
10 passes) and the parametric sweep (`E_phi_sweep_500GHz`, SaveFields true, enabled) were saved as
intended. So the ClcEval failure of baseline1 followed from having no data; whether the stack
read-back also fails over gRPC with data present (what PyAEDT and pyEPR document) is untested, and the
`CalculatorWrite` path stays as the sanctioned one. `sacct`: job 4005538 `FAILED 1:0`, 13:08:22 to
13:10:49, MaxRSS 1.72 GB; it ended two seconds after Python did (the "still running" impression was
those seconds).

Two defects in the diagnostic itself: a PyAEDT-wrapped `export_mesh_stats` call raised, and PyAEDT's
default `settings.release_on_exception` then closed the desktop, so the extraction test never ran; and
`.profile` files were not printed. Fixed with the runner change below.

Runner change (this entry's commit): after the solve the runner reads AEDT's message manager into
the log, asks the Solutions module for the solved variations, and stops with a RuntimeError naming the
missing polarizations and quoting AEDT's messages when any requested `Ephi='n'` is absent; the manifest
records `solver.solved_variations`; `aedt_session` sets `settings.release_on_exception = False` so
PyAEDT never closes the desktop behind the runner's back.

`opti906_0.profile` (2026-10-06): the Optimetrics run's record. Host hpc-25-06, HFSS 2025.2.0; the
sweep started at 13:10:19 PDT, ran two phases (13:10:19 to 13:10:24, progress 0; 13:10:24 to 13:10:42,
progress 1.0) and was marked `Finished` at 13:10:42 without solving a variation. Between 13:09:37 and
13:10:19 `Analyze All` spent 42 s on the nominal setup and left no mesh and no pass. `slurm-4005538.out`
holds no AEDT message at all: PyAEDT prints "Non-graphical mode detected. Disabling Desktop logs." at
startup and stops forwarding AEDT's message window, so the solver's reason was never written down.

Licensing client log `~/.ansys/ansyscl.hpc-25-06.1731008.39413.log` (2026-10-06): the licence is not
the cause. At 13:09:42 the solver engine (`HFSSCOMENGINE.EXE`, pid 1732508, "Using HPC Parametric
context") checked out `elec_solve_hfss` (1/1/1/25: granted, 1 of 25 in use) and the HPC pack features
`elec_solve_level2` and `elec_solve_level1` (1 of 30 each); it checked them in at 13:10:19 after 37 s,
having produced no mesh and no pass. A second engine (pid 1735543, the sweep) checked the same features
out and in twice between 13:10:24 and 13:10:42 and shut down at 13:10:43; the final check-ins show
0 in use. No `DENIED` anywhere. So the solver started twice with its licences and aborted each time
before writing anything: a meshing or validation error, or a failure to launch the solver's worker
processes on the node; the reason was in AEDT's message window, which nothing recorded.

`baseline2` (2026-10-06, job 4065548, hpc-24-18, code at `29cd12a`): **the solve completed.** AEDT's
messages, now logged by the runner:

    [warning] Excitation 'plane_wave_500GHz' is an analytical incident wave with geometry assignment and
              it requires total field formulation. (09:05:03, at creation; the runner then selects total fields)
    [info]    Normal completion of simulation on server: hpc-24-18. (09:05:45)
    [info]    Parametric Analysis on E_phi_sweep_500GHz has been started. (09:05:45)
    [info]    A variation (Ephi='0') has been requested using the following machines: hpc-24-18. (09:05:50)
    [info]    A variation (Ephi='1') has been requested using the following machines: hpc-24-18. (09:05:51)
    [info]    Parametric Analysis is done. (09:06:09)
    solved variations of 500GHz : LastAdaptive: ["Ephi='0'", "Ephi='1'"]

Then the first `CalculatorWrite` of `outgoing_power` failed with a bare
`GrpcApiError: Failed to execute gRPC AEDT command: CalculatorWrite`, and the run closed the project
unsaved again (the runner saved it only before the solve).

**Correction to the entries above** (2026-10-06). baseline1 solved as well: the same 42 s nominal
solve and 24 s sweep, the same licence pattern, and baseline2's messages for an identical setup. The
step 4a copy of baseline1 held no solution because the runner saved the project only before the solve
and the failed run closed it unsaved; AEDT evidently puts solution data into the results folder when
the project is saved (the `--solve` diagnostic lists the folder before and after a save to confirm).
Withdrawn: "the solve produced no solution", and "ClcEval does not work over gRPC" as the established
cause (the gRPC limitation is documented by pyEPR but our evidence does not show it: the step 4a probes
ran on a copy without a solution). What is established: evaluating `outgoing_power` with the runner's
arguments fails on a solved design by both `ClcEval` (baseline1) and `CalculatorWrite` (baseline2).
Candidates: the variation arguments (`Ephi` passed as an integer, the incident-angle intrinsics
`IWavePhi`/`IWaveTheta` of a single-direction plane wave) or the expression's inputs (the face list
`outgoing`, `Vector_RealPoynting` after `ClearAllNamedExpr`).

Runner changes (this entry's commit): the project is saved right after the verified solve, so a failed
extraction keeps the solve; any failure after the solve logs AEDT's new messages before the session
closes. Next: step 4a with `--solve` on `job_baseline2`: solves a copy of its saved pre-solve project,
saves it, and evaluates `outgoing_power` with the runner's arguments and with variants in one session.

Step 4a with `--solve` (2026-10-06, job 4066881, code at `e9ea19a`) on a copy of `job_baseline2`:
**root cause found.** The copy solved (Ephi=0 converged in 4 passes, Ephi=1 in 2) and saved. Then:

| Probe | Result |
|---|---|
| constant 1, Freq and Phase only | OK, 1.0 |
| constant 1, the runner's arguments (with IWavePhi, IWaveTheta) | FAILED |
| outgoing_power with the incident angles, any form (Ephi int or text, 0.0deg or 0deg) | FAILED |
| outgoing_power without the incident angles (Ephi int or text, or Freq and Phase only) | OK, 6.6471e-10 W, T = 1.00167 |
| outgoing_power, Ephi=1, without the incident angles | OK, 1.08e-19 W, T = 1.6e-10 |
| ExportOnGrid with the incident angles | FAILED |
| ExportOnGrid without the incident angles | OK, 5151 points |
| ExportFieldsToFile with the runner's keys, and with Freq only | OK, 19388 points each |
| named expressions `outgoing_power`, `Vector_RealPoynting`; face list `outgoing` = [7] | present |
| PyAEDT's own `fields_calculator.evaluate` | FAILED (its variation lists every design variable) |

A plane wave with a single incident direction (PhiPoints = ThetaPoints = 1) defines no incident-wave
variables, so naming `IWavePhi`/`IWaveTheta` fails every calculator evaluation and grid export, even a
constant; the one direction is then implied. The legacy script always swept 15 directions, which is why
this never appeared on Windows. `Ephi` as an integer is fine, so neither gRPC marshalling nor the
read-back method was the cause. Fix (this entry's commit): `field_variation` names the angles only when
the plane wave has more than one direction; the far-field key is unchanged (AEDT accepts it either way).
A sweep along one axis only is unverified.

Transmission at (0 deg, 180 deg), Ephi=0: 1.00167 here (single-direction solve, 2025 R2) against 1.05451
in the Windows reference (15-direction solve, 2023 R2): -5.0 %. The reference itself gives 0.99607 at
(90 deg, 180 deg) for Ephi=1, the physically identical incident field (E along the gap, normal incidence),
because each polarization is a separate parametric variation with its own adaptive mesh: its internal
spread between the two is 5.9 %, so this baseline lies within the method's mesh noise at MaxDeltaE 0.02.
`compare_hfss_exports.py` allows |dT| <= max(1e-3, 0.05 T), which this key misses by 0.0001; step 5 is a
15-direction solve with its own mesh, so the tolerance is left unchanged and taken to Rohan as a decision
if step 5 misses it.

**Step 4: PASS** (2026-10-06, job 4068495, hpc-22-11, code at `9b9f1db`, `--job-id baseline3`). The solve
took 62.5 s (AnalyzeAll; "Normal completion of simulation"; parametric sweep done; solved variations
Ephi='0' and Ephi='1'); the project was saved right after it; then:

    [power 1/1] phi=0.0 theta=180.0 Ephi=0 -> 6.647119e-10 W        T = 1.001668
    exit-field grid in outgoing_cs: 5151 lattice points; [exit 1/1] 5151 points
    [far 1/1] 19388 points
    wrote .../job_baseline3/InfParallelPlate_crack1Rohan_500GHz_Ephi=0 (5151 waveguide rows, 19388 far-field rows)
    done; manifest at /home/rshenoy/BBRSim/outputs/InfParallelPlate_crack1Rohan_500GHz/job_baseline3/manifest.json

T equals the step 4a value to every printed digit (same setup, same mesh). Wall time on the node about
1 min 51 s, of which about 45 s AEDT start-up. MaxRSS and the results-folder size: recorded with step 5.

## 5. Full reference run and comparison

2026-10-06, run job 4069263 (hpc-89-25, code at `a9df6e4`, `--job-id reference1`, `-t 02:00:00`) and
comparison job 4069264 (`hpc/compare_reference.sbatch`, dependent `afterok`). **Run: PASS; comparison:
FAIL on transmission only.**

    run:   COMPLETED 0:0, Elapsed 00:06:00, MaxRSS 1995148K (batch step); job directory 345M
           Ephi=0 and Ephi=1: 77265 waveguide rows, 290820 far-field rows each; done; manifest written
    reference CSVs (HPC BBRSimulation checkout): all four sha256 OK
    Ephi=0: columns, Freq, Ephi, keys (15), ingoing power (rel 0), waveguide and far-field points: PASS
            transmission FAIL: 3 pass, 7 below noise floor, 5 fail at (0,45) (0,90) (0,135) (0,180) (45,180);
            max rel diff 0.055 at (0,135)
            fields WARN: waveguide min |E| correlation -0.2472, max nRMS 0.3720; far field -0.2101, 0.2723
    Ephi=1: same structural PASSes
            transmission FAIL: 4 pass, 3 below noise floor, 8 fail at (0,45) (0,90) (0,135) (0,180)
            (90,45) (90,90) (90,135) (90,180); max rel diff 7.3e9 at (0,180), where the reference is 1.4e-10
            fields WARN: waveguide min corr -0.0571, max nRMS 0.3437; far field 0.0358, 0.2207

The schema, keys, grids and ingoing power are identical to the Windows reference: the pipeline
reproduces the dataset's structure exactly on AEDT 2025 R2. The transmission values differ by up to
5.5 % at the transmitting keys, against the reference's own 5.9 % spread between physically identical
keys (step 4). At Ephi=1, phi=0 the candidate exceeds the 1e-6 noise floor where the reference holds
1e-12 to 1e-10 (E along the long side, below the TE01 cutoff). The field minima are aggregates over all
15 keys, cut-off ones included. Per-key detail: `hpc/compare_detail.py` (next).

Per-key detail (`hpc/compare_detail.py`, 2026-10-06, login node):

    convergence  Ephi=0: 7 passes, last delta 0.032014, 0.030064, 0.015993 (26008 elements)
                 Ephi=1: 7 passes, last delta 0.036096, 0.022741, 0.0056319 (23001 elements); solve 189.5 s
    T (phi, theta)   ref E0      cand E0     ref E1      cand E1
      (0, 45)       2.236e-02   2.112e-02   2.963e-12   2.112e-02
      (0, 90)       2.577e-01   2.434e-01   3.467e-11   2.434e-01
      (0, 135)      7.595e-01   7.174e-01   1.006e-10   7.174e-01
      (0, 180)      1.0545      0.99605     1.367e-10   0.99605
      (45, 45)      1.199e-02   1.160e-02   1.160e-02   1.160e-02
      (45, 90)      1.337e-01   1.312e-01   1.311e-01   1.312e-01
      (45, 135)     4.073e-01   3.941e-01   3.939e-01   3.941e-01
      (45, 180)     0.52725     0.49802     0.49803     0.49802
      (90, 45..180) ~1e-11      ~1e-11      2.276e-02 .. 0.99607   ~1e-11
    fields at the transmitting keys: exit |E| correlation 0.9990 (E0) and 0.9999 (E1), far field 0.9999 / 1.0000

**Bug found: the candidate's Ephi=1 OutgoingPower is a copy of its Ephi=0 OutgoingPower** at all 15 keys,
to every printed digit, while its Ephi=1 exit and far fields match the reference's Ephi=1 fields
(correlation 0.9999 at the transmitting keys). The runner passed the polarization as an integer
(`"Ephi:=", 1`); over gRPC, CalculatorWrite evidently ignored it and evaluated the nominal Ephi=0 without
an error, while ExportOnGrid honoured the same integer (step 4a: text "1" gave the Ephi=1 value).
Fix (this entry's commit): `field_variation` passes the polarization as text, as PyAEDT passes design
variables; the fake AEDT now ignores an integer Ephi in CalculatorWrite, and the existing test that
each polarization carries its own power reproduced the bug before the fix. The reference1 dataset must
not be used. The diagnostic gained an integer-versus-text probe sequence at a chosen key to confirm the
mechanism on the saved reference1 solution.

Against the reference, the candidate's Ephi=0 agrees with the reference's Ephi=1 at the physically
identical keys (normal incidence: 0.99605 against 0.99607; (45, theta): within 0.05 %), and both lie a
constant 5.54 % below the reference's Ephi=0 at phi = 0: the reference's Ephi=0 variation is the
outlier. Once the Ephi=1 column is right, the Ephi=0 keys at phi = 0 and (45, 180) will still miss the
5 % tolerance by 0.5 %; that is Rohan's decision, with this evidence.

Polarization probe (2026-10-06, job 4071971, diagnostic on the saved reference1 solution, key (90, 180),
no new solve): integer 0 and integer 1 both gave T = 3.534e-11 (the Ephi=0 value, cut off), also right
after a text 1; text "1" gave T = 0.99937; text "0" gave 3.534e-11. ExportOnGrid gave the same field for
integer 1 and text 1 (max |E| 1.2935 V/m). **Root cause confirmed:** over gRPC, CalculatorWrite ignores
an integer design-variable value and evaluates the nominal variation; no caching is involved.

**reference2** (2026-10-06, run job 4071972, code at `470bd48`; comparison job 4071973): run done
(solve 198.6 s; Ephi=0 and Ephi=1 converged in 7 passes, last delta 0.0160 and 0.0056). Comparison:

    reference CSVs: all four sha256 OK
    Ephi=1: RESULT: PASS. Transmission 8 pass, 7 below noise floor, 0 fail; max rel diff 0.003 at (45,180).
    Ephi=0: RESULT: FAIL. Transmission 3 pass, 7 below noise floor, 5 fail at (0,45) (0,90) (0,135) (0,180)
            (45,180), each exactly -5.544 %; (45,45) -3.24 %, (45,90) -1.87 %, (45,135) -3.24 % pass.
    Structure (columns, Freq, Ephi, 15 keys, ingoing power, waveguide and far-field grids): PASS, both.
    Fields: |E| correlation >= 0.999 (exit) and >= 0.9999 (far) at every transmitting key; the WARN
            minima come from cut-off keys (T <= 1e-10), whose fields are numerical noise on both sides.

Self-consistency at normal incidence (E along the gap; physics makes (0,180) Ephi=0 and (90,180) Ephi=1
the same incident field): ours 0.99605 and 0.99937, ratio 0.9967; the Windows reference 1.05451 and
0.99607, ratio 1.0587. Our Ephi=0 equals the reference's Ephi=1 at these keys, and our Ephi=1 matches
the reference's Ephi=1 within 0.33 %. The five failing keys are the Windows Ephi=0 variation's constant
5.54 % offset (its own adaptive mesh), not a difference in the 2025 R2 dataset. Acceptance of step 5 with
this documented exception, a tolerance change, or a convergence study first: Rohan's decision.
Comparison records: `/home/rshenoy/BBRSim/outputs/InfParallelPlate_crack1Rohan_500GHz/job_reference2/comparison_Ephi{0,1}.json`.

## 5b. Incident-direction check

## 6. Round-gap project

## 7. Round gap, single angle at 2000 GHz

## 8. Round gap, full sweep

## Discrepancies and decisions
