#!/bin/bash -l
# Step 8d (hpc/README.md): the OR-035 checks on copies of the three solved projects, with no solve (HF-045, CDX-023).
# One AEDT session per project, no solver licence. hpc/check_or035.py does the work; its docstring states both checks.
#
# Usage (HPC):  mkdir -p /home/rshenoy/BBRSim/outputs   # Slurm needs the --output directory to exist
#               sbatch hpc/check_or035.sbatch [--crack1 JOB] [--crack2 JOB] [--round-gap JOB]
# Defaults: crack1 and crack2 job_tight1 (step 5d), the round gap job_reference1 (step 8), under
# /home/rshenoy/BBRSim/outputs. Each job's project is copied into <job>/check_or035/ (crack1's with its .aedtresults,
# about 2.5 GB, for the export); AEDT opens only the copies, and the source jobs are never written. The copies are
# removed after a passing run and kept after a failure, which names where.
# Resources are first guesses: a few minutes per project, the copy included.
#
# Expect: (a) "expect 101 x 39 = 3939 rows per key, z end nodes at -1.875e-05 and 1.875e-05 m to 1e-12 m", the
# fitted step 0.0375/38 mm, 15 keys exported with rows per key [3939], "HF-045: PASS"; (b) max |x| <= 1e-09 mm and
# PASS for crack1, crack2 and the round gap; the last line "== RESULT: PASS" and exit 0. Otherwise the last line is
# "== RESULT: FAIL <what failed>" with exit 1, or exit 2 when nothing was started (a job that cannot be read, or an
# existing <job>/check_or035/ from an earlier attempt, which is removed by hand before a rerun).
#SBATCH --job-name=bbsim_check_or035
#SBATCH --account=golwala
#SBATCH --partition=expansion
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/home/rshenoy/BBRSim/outputs/check_or035-%j.out

BBSIM_ENV=${BBSIM_ENV:-/home/rshenoy/BBRSim/bb_env.sh}
# bb_env.sh is sourced before strict mode, because module and the env scripts are not proven under set -u.
if ! source "$BBSIM_ENV"; then
  echo "ERROR: sourcing $BBSIM_ENV failed" >&2
  echo "== RESULT: FAIL environment: sourcing $BBSIM_ENV failed; nothing started"
  exit 1
fi
# Slurm exports the submitting shell: a stale PYAEDT_DESKTOP_PORT or PYAEDT_PROCESS_ID would make PyAEDT attach to
# another desktop instead of starting its own; bbsim.session.aedt_session also removes them.
unset "${!PYAEDT_@}"
# No set -e: the script's exit status is recorded, and its verdict line decides the job's.
set -uo pipefail
CHECKOUT=${BBSIM_CHECKOUT:-/home/rshenoy/BBRSim/Blackbody-Simulations}
if ! cd "$CHECKOUT"; then
  echo "== RESULT: FAIL checkout $CHECKOUT not found; nothing started"
  exit 1
fi
PYTHON=${BBSIM_PYTHON:-python}
echo "host=$(hostname) job=${SLURM_JOB_ID:-none} checkout=$CHECKOUT commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"

LOG=$(mktemp "${TMPDIR:-/tmp}/check_or035.XXXXXX")
"$PYTHON" -u hpc/check_or035.py "$@" 2>&1 | tee "$LOG"
rc=$?
last=$(tail -n 1 "$LOG")
rm -f "$LOG"
# The script ends with its own verdict line; a crash before it (an import error, a killed process) is a FAIL here.
case $last in
  "== RESULT: PASS"*)
    [ "$rc" -eq 0 ] && exit 0
    echo "== RESULT: FAIL hpc/check_or035.py printed PASS but exited $rc"
    exit 1 ;;
  "== RESULT: FAIL"*)
    [ "$rc" -eq 2 ] && exit 2
    exit 1 ;;
esac
echo "== RESULT: FAIL hpc/check_or035.py exited $rc without a verdict"
[ "$rc" -eq 2 ] && exit 2
exit 1
