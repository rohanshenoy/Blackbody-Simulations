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

**Step 5: PASS with a recorded exception** (decided by Rohan, 2026-10-06). The five Ephi=0 keys at
phi = 0 and (45, 180) differ from the Windows reference by its own constant 5.54 % Ephi=0 offset; the
transmission tolerance stays at 5 % for later comparisons. Note from the BBRsim side: its loader
normalizes a key only when the coherent T exceeds 1, so after loading the two datasets differ at
(0, 180) by 0.4 %, not zero (legacy 1.05451 becomes 1; reference2's 0.99605 stays).

## 5b. Incident-direction check

2026-10-06, job 4073526 (`hpc/check_incident_direction.py` on a copy of `job_reference2`'s solved
project, code at `ce9fd30`). **PASS.**

    (IWavePhi, IWaveTheta) = (90, 135) deg, Ephi=0, 500 GHz
    slope_y = 7409.05 rad/m, slope_z = -7408.35 rad/m, k0 = 10479.23 rad/m
    ratio_ky_over_kz = -1.000095   magnitude_over_k0 = 0.999834
    verdict: PASS: arrival direction, k = -r_hat(theta, phi)
    slope_z_negative_as_expected: true

HFSS applies one direction to both angles, and the incident wave travels along -r_hat(theta, phi): the
incidence definition in the dataset sidecar and BBRsim's azimuth formula are confirmed by data. This
was the frame contract's last unproven assumption (ledger G4-004, HF-012).

## 5c. crack2 reference run and comparison

2026-10-06, run job 4073527 (`configs/crack2_500GHz_reference.toml`, `--job-id reference1`, code at
`ce9fd30`) and comparison job 4073528. Run done (solve 95.1 s; Ephi=0 and Ephi=1 converged in 4 passes,
last delta 0.0087 and 0.0131, about 11.6k elements). Comparison against the Windows crack2 reference:

    reference CSVs: all four sha256 OK; structure (columns, Freq, Ephi, 15 keys, ingoing power
    1.327209365e-09 W, waveguide and far-field grids): PASS, both polarizations
    Ephi=0: RESULT: PASS. Transmission 8 pass, 7 below noise floor; max rel diff 0.035 at (0, 90).
    Ephi=1: RESULT: FAIL. Transmission 7 pass, 7 below noise floor, 1 fail at (90, 90): 0.02906 against
            0.03169, -8.3 % (grazing incidence); elsewhere within +4.1 %.
    Fields at the transmitting keys: exit |E| correlation >= 0.998 (0.994 at (90, 90)), far >= 0.998.

    T (phi, theta)  ref E0    cand E0   ref E1    cand E1
      (0, 45..180)  cand/ref = 1.03501 at every theta for Ephi=0
      (0, 180)      1.00873   1.04403   7.2e-10   5.5e-09
      (45, 180)     0.50436   0.52201   0.50140   0.52199
      (90, 90)      6.1e-10   1.3e-09   0.03169   0.02906
      (90, 180)     1.5e-09   5.5e-09   1.00280   1.04399

Self-consistency at normal incidence: ours 1.04403 (Ephi=0, (0,180)) and 1.04399 (Ephi=1, (90,180)),
ratio 1.0000; the Windows reference 1.00873 and 1.00280, ratio 1.0059. Our crack2 transmission is
3.5 % to 4.1 % above the reference's, uniformly within each polarization, the opposite sign of crack1,
where the Windows Ephi=0 was 5.54 % high. Per-variation offsets of a few percent in both directions
point at the convergence target: Delta Mag Energy 0.02 leaves T uncertain by several percent per solve.

## 5d. Convergence check (MaxDeltaE 0.005)

2026-10-06, jobs 4075000 (crack1, `configs/crack1_500GHz_reference_tight.toml`) and 4075001 (crack2,
`configs/crack2_500GHz_reference_tight.toml`), `--job-id tight1`, code at `c465069`. Both COMPLETED.

    crack1: Elapsed 00:08:42, MaxRSS 2.84 GB, solve 340.8 s
            Ephi=0 9 passes, last delta 0.0036596 (38898 elements); Ephi=1 11 passes, 0.0037117 (64396)
    crack2: Elapsed 00:14:25, MaxRSS 4.84 GB, solve 613.2 s
            Ephi=0 13 passes, last delta 0.004126 (107406 elements); Ephi=1 12 passes, 0.00434 (81044)

    T, E along the gap at normal incidence      Ephi=0 (0,180)   Ephi=1 (90,180)   ratio
      crack1 tight                              0.99234          0.99120           1.0011
      crack1 at 0.02 (reference2)               0.99605          0.99937           0.9967
      crack1 Windows (2023 R2)                  1.05451          0.99607           1.0587
      crack2 tight                              0.99226          0.99305           0.9992
      crack2 at 0.02 (reference1)               1.04403          1.04399           1.0000
      crack2 Windows (2023 R2)                  1.00873          1.00280           1.0059

    tight against the 0.02 run, transmitting keys: crack1 -0.82 % to +0.13 %; crack2 -4.96 % to -1.47 %,
      except the grazing key (90, 90) Ephi=1, +10.3 %
    tight against Windows, transmitting keys: crack1 Ephi=1 -0.49 % to +0.18 %, crack1 Ephi=0 -5.90 % at
      phi = 0 (the Windows Ephi=0 offset) and -3.5 % to -2.0 % at phi = 45; crack2 -1.63 % to +1.19 %

At the reference target 0.02 a solve can carry about 5 % error in T (crack2 at 0.02, and the Windows
crack1 Ephi=0 set); at 0.005 both cracks converge, the two separately meshed polarizations agree to
0.1 %, and the 2025 R2 data agree with the Windows data within 1.6 % everywhere except the Windows crack1
Ephi=0 set. Both cracks give T = 0.992 at normal incidence. Cost at 500 GHz against 0.02: crack1 1.5x
the wall time and 1.4x the memory, crack2 about 3.5x and 2.4x. Per HF-025 the two-crack tree takes the
tight runs.

## 5e. Two-crack 500 GHz data tree

2026-10-06, login node, from the step 5d runs (`job_tight1` of both cracks), per Rohan's data decisions
(ledger G4-016, HF-023, HF-025; record HF-028): `/resnick/groups/golwala/rshenoy/bbsim/trees/cracks_500GHz_2025r2`,
the data root BBRsim reads through `/bbr/dataDir` or `BBRSIMDATA`.

    sha256sum: all OK
    waveguides/: InfParallelPlate_crack1Rohan_500GHz{.dataset.json, _Ephi=0, _Ephi=1}
                 InfParallelPlate_crack2_500GHz{.dataset.json, _Ephi=0, _Ephi=1}
    SHA256SUMS: 14 files; du -sh 252M (108.6 + 143.0 MiB, the per-crack sizes measured before)

The tight runs against the Windows reference (jobs 4077456 crack1 and 4077457 crack2, each writing
`comparison_Ephi{0,1}.json` into its `job_tight1`):

    crack1 Ephi=0: FAIL, transmission 3 pass, 7 below noise floor, 5 fail at (0,45) (0,90) (0,135) (0,180)
                   (45,180), max rel diff 0.059 at (0,135); fields WARN only
    crack1 Ephi=1: PASS (fields WARN only)
    crack2 Ephi=0 and Ephi=1: PASS (fields WARN only)

The five crack1 Ephi=0 failures are the keys of step 5's accepted exception, the Windows Ephi=0 set's constant
offset (HF-021); the field warnings come from the keys below the noise floor, as in step 5. The tree went to the
BBRsim side (HF-028), whose HPC validation is `g4_full_check.sbatch` followed by `g4_tree_check.sbatch` (G4-018).

## 6. Round-gap project

2026-10-06, job 4077741 (debug QOS, chained with steps 7 and 8, code `c8519be`): COMPLETED in 00:01:18,
MaxRSS 1.42 GB. PASS.

    Build verified: /resnick/home/rshenoy/BBRSim/projects/RoundGap.aedt
    RoundGap.build.json: verification_differences [], error None
    RoundGap.inventory.json: design round_gap_r50um, units mm, object gap [vacuum],
                             bounding box [-0.05, -0.05, 0.0, 0.05, 0.05, 0.4], 3 faces

## 7. Round gap, single angle at 2000 GHz

2026-10-06, job 4077742 (`roundgap1`, Ephi=0): FAILED after 00:02:12 (MaxRSS 1.84 GB) at the first exit-field
export after the solve:

    FieldFileError: exit field for phi=0.0 theta=180.0 Ephi=0: 50 exit points lie outside the declared disc
    cross-section; HFSS must write nan, not values, outside the solid

Step 8 (job 4077743) was cancelled by its dependency, as designed. The runner had already deleted the raw
export; step 7a found the cause.

## 7a. Exit-lattice diagnostic

2026-10-06, job 4089551 (`diag_exitgrid1`, code `e173f1e`): the runner unchanged on the single-angle config
with both polarizations, every raw export kept, the containment check reported instead of raised. COMPLETED
in 00:02:12, MaxRSS 1.81 GB.

                                              Ephi=0           Ephi=1
    lattice rows (101 x 101)                  10201            10201
    lattice points inside the disc            7845 (20 on the rim)
    inside with field                         7783             7786
    inside without field (holes)              62, to 756 nm    59, to 523 nm inside the rim
    outside with field (spill)                50, to 912 nm    56, to 1078 nm beyond it
    rim points with field                     7 of 20          7 of 20
    |E| at spill / median within 2 um inside  0.999 (0.03-1.38) 0.984 (0.03-1.46)

The two exports share the lattice but not the pattern (3 holes and 8 spill points differ). HFSS meshes the
circle as a polygon (the runner keeps curvilinear elements off, as in the crack reference designs), so its
wall cuts up to about 0.76 um inside the circle, and the export gives points up to about 1 um beyond it the
field of the edge elements: the spill is the field, not junk. Even without the containment check the run
could not have been written: the polarizations' valued counts differ (7833 and 7842), and BBRsim pairs their
rows by position.

Decisions (Rohan, 2026-10-06):
- Rim rule (ledger HF-029): for a disc, every lattice point inside it is kept for every key and both
  polarizations (7845 at R = 50 um and a 1 um step), in export order; valued points outside are dropped;
  inside points without field get zero field, which BBRsim never samples (its loader accepts them, confirmed
  by the Geant4 session); anything more than 5 % of R beyond or inside the rim stops the run; `manifest.json`
  `exit_field.rim` records the counts per key and polarization. Rectangles are unchanged.
- Wall facets (ledger HF-030): measure first (step 7b), the single-angle run as configured and with
  `solver.wall_normal_deviation_deg = 5.0`, both polarizations; then choose for step 8.

## 7b. Wall facets

2026-10-07, code `fa6669d`, single angle at 2000 GHz with both polarizations, MaxDeltaE 0.005: job 4155776
(`roundgap2`, HFSS's default faceting) and job 4155777 (`roundgap2_facets5`, `wall_normal_deviation_deg = 5.0`,
read back from the design before the solve). Both COMPLETED.

                                        default facets          5 deg facets
    T, Ephi=0 / Ephi=1                  0.47632 / 0.47639       0.50958 / 0.50969
    T Ephi=1 / Ephi=0                   1.00015                 1.00021
    rows per polarization               7845                    7845
    holes zeroed (deepest)              62 / 59 (756 / 523 nm)  1 / 1 (0 nm, on the rim)
    outside values dropped (farthest)   50 / 56 (912 / 1078 nm) 57 / 61 (606 nm)
    passes, last delta                  3, 0.0049 (6911 tets)   4, 0.0025 / 0.0031 (15258 / 15058 tets)
    solve, wall time, MaxRSS            57 s, 00:02:19, 1.78 GB 97 s, 00:02:59, 2.09 GB

    FACETS: T 5 deg / default = 1.06983 (Ephi=0), 1.06989 (Ephi=1)

HFSS's default faceting makes the guide effectively narrower: 2000 GHz is only 14 % above the TE11 cutoff, so
the raised cutoff, and the smaller faceted exit face over which the outgoing power is integrated, put T 7 %
low. The default run reproduced step 7a's rim counts exactly (the meshing is deterministic). The holes came from
the polygonal wall and vanish at 5 deg; the values beyond the rim persist at 5 deg (up to 0.6 um, with the wall
within 0.05 um of the circle), so they come from the export's point location, and the rim rule (HF-029) stays.
Decision (Rohan, 2026-10-07, ledger HF-032): both round-gap configs carry 5 deg; step 8 runs at 5 deg alongside
a single-angle 2.5 deg check that 5 deg has converged.

## 8. Round gap, full sweep

2026-10-07, code `18ff009`, both at MaxDeltaE 0.005 and both polarizations: job 4176305 (`reference1`, the
15-direction sweep at 5 deg facets) and job 4176308 (`roundgap3_facets2p5`, single angle at 2.5 deg). Both
COMPLETED. PASS: the round gap's first dataset, `/home/rshenoy/BBRSim/outputs/RoundGap_r50um_2000GHz/job_reference1`.

    step 8: 00:04:16, MaxRSS 2.47 GB, solve 171 s; wall facets 5 deg (read back)
    sidecar: 7845 per key, disc radius 5e-05 m, TE11, 1 propagating mode
    Ephi=0: 117675 rows, 15 keys x 7845; T at theta 180, phi 0/45/90: 0.51013 0.51031 0.51001;
            rim zeroed <= 1, dropped <= 58 per key; 5 passes, last delta 0.0025 (21999 tets)
    Ephi=1: 117675 rows, 15 keys x 7845; T at theta 180, phi 0/45/90: 0.50977 0.50976 0.50986;
            rim zeroed <= 1, dropped <= 61 per key; 5 passes, last delta 0.0029 (21956 tets)
    spread of the six normal-incidence T: 0.108 %

    faceting check, single angle: 5 deg (step 7b) T 0.50958 / 0.50969, 97 s;
                                  2.5 deg 00:05:03, MaxRSS 3.14 GB, T 0.51139 / 0.51111, 216 s, 4 passes,
                                  last delta 0.0010 (33559 tets)
    FACETS: T 2.5 deg / 5 deg = 1.00354 (Ephi=0), 1.00278 (Ephi=1)

At normal incidence the round gap prefers no polarization and no azimuth: the six values agree to 0.11 %, and
the sweep's normal-incidence T matches the single-angle 5 deg run to 0.1 %. 2.5 deg moves T by 0.3 %, inside the
0.5 % criterion, so 5 deg stands (HF-032); the change from the default (+7.0 %) to 5 deg and on to 2.5 deg
(+0.3 %) is consistent with an error falling as the square of the facet angle, which puts 5 deg about 0.3 to
0.5 % below the converged T, within the 1 % that MaxDeltaE 0.005 aims for (HF-027).

## Discrepancies and decisions
