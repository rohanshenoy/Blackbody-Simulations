#!/bin/bash -l
# Step 10b (hpc/README.md): the key-by-key acceptance of the three 35-key solves against the 15-key solves of the same
# geometry (OR-036), batch C2. No AEDT and no licence: bbsim.grid_acceptance reads the jobs' sidecars and tables.
#
# Usage (HPC):  sbatch --dependency=afterok:<crack1>:<crack2>:<round gap> --kill-on-invalid-dep=yes \
#                 hpc/grid_acceptance.sbatch
# Pairs (O = /home/rshenoy/BBRSim/outputs; BBSIM_OUTPUTS overrides it), each new job made by step 10a with
# --job-id grid35_1:
#   crack1     $O/InfParallelPlate_crack1Rohan_500GHz/job_grid35_1  against  .../job_tight1      (step 5d)
#   crack2     $O/InfParallelPlate_crack2_500GHz/job_grid35_1       against  .../job_tight1      (step 5d)
#   round gap  $O/RoundGap_r50um_2000GHz/job_grid35_1               against  .../job_reference1  (step 8)
# python -m bbsim.grid_acceptance runs on each pair, printing its table and writing <new job>/grid_acceptance.json; a
# failing pair does not stop the others. Each is judged by its exit code, never by its last line.
# It runs from this checkout (BBSIM_CHECKOUT overrides it), with bb_env.sh's python (BBSIM_ENV, BBSIM_PYTHON).
#
# Expect: three tables, each with its 18 rows (nine shared keys, two polarizations) and "GRID35: PASS"; the last line
# "== RESULT: PASS" and exit 0. Otherwise the last line is "== RESULT: FAIL <datasets>" naming each failing pair with
# its exit code, exit 1; arguments (it takes none) exit 2.
#SBATCH --job-name=bbsim_grid35_acceptance
#SBATCH --account=golwala
#SBATCH --partition=expansion
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/home/rshenoy/BBRSim/outputs/grid_acceptance-%j.out

BBSIM_ENV=${BBSIM_ENV:-/home/rshenoy/BBRSim/bb_env.sh}
# bb_env.sh is sourced before strict mode, because module and the env scripts are not proven under set -u.
if ! source "$BBSIM_ENV"; then
  echo "ERROR: sourcing $BBSIM_ENV failed" >&2
  echo "== RESULT: FAIL environment: sourcing $BBSIM_ENV failed; nothing started"
  exit 1
fi
# Slurm exports the submitting shell: no PyAEDT is started here, but a stale PYAEDT_* never reaches the children.
unset "${!PYAEDT_@}"
# No set -e: every pair's exit status is recorded, and the last line decides the job's.
set -uo pipefail
if [ $# -ne 0 ]; then
  echo "usage: sbatch hpc/grid_acceptance.sbatch   (no arguments; BBSIM_OUTPUTS sets the outputs root)"
  echo "== RESULT: FAIL usage: it takes no arguments, got $#; nothing started"
  exit 2
fi
CHECKOUT=${BBSIM_CHECKOUT:-/home/rshenoy/BBRSim/Blackbody-Simulations}
if ! cd "$CHECKOUT"; then
  echo "== RESULT: FAIL checkout $CHECKOUT not found; nothing started"
  exit 1
fi
PYTHON=${BBSIM_PYTHON:-python}
O=${BBSIM_OUTPUTS:-/home/rshenoy/BBRSim/outputs}
PAIRS="InfParallelPlate_crack1Rohan_500GHz/job_tight1 InfParallelPlate_crack2_500GHz/job_tight1 RoundGap_r50um_2000GHz/job_reference1"
echo "host=$(hostname) job=${SLURM_JOB_ID:-none} checkout=$CHECKOUT commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"

FAILED=""
for P in $PAIRS; do
  DATASET=${P%/*}
  echo "== grid_acceptance $DATASET: $O/$DATASET/job_grid35_1 against $O/$P"
  "$PYTHON" -m bbsim.grid_acceptance --new "$O/$DATASET/job_grid35_1" --reference "$O/$P"
  rc=$?
  echo "== grid_acceptance $DATASET exited $rc"
  [ "$rc" -eq 0 ] || FAILED="$FAILED${FAILED:+; }$DATASET (exit $rc)"
done

if [ -n "$FAILED" ]; then
  echo "== RESULT: FAIL $FAILED"
  exit 1
fi
echo "== RESULT: PASS"
exit 0
