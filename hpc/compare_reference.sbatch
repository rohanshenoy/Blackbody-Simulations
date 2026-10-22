#!/bin/bash -l
# Usage (HPC): sbatch hpc/compare_reference.sbatch <job directory> [<stem>] [<reference root>]
#   Compares a job's two polarizations against the Windows reference datasets (AEDT 2023 R2) with
#   compare_hfss_exports.py, writing <job>/comparison_Ephi{0,1}.json. <stem> defaults to
#   InfParallelPlate_crack1Rohan_500GHz; the reference root defaults to the CSVs tracked in the HPC
#   BBRSimulation checkout, whose sha256 must match configs/legacy_reference_500GHz.sha256 first.
#   Step 5 submits it with --dependency=afterok on the run (hpc/README.md).
#SBATCH --job-name=bbsim_compare
#SBATCH --account=golwala
#SBATCH --partition=expansion
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/home/rshenoy/BBRSim/outputs/slurm-%j.out
JOB=${1:?job directory required}
STEM=${2:-InfParallelPlate_crack1Rohan_500GHz}
REF=${3:-/home/rshenoy/BBRSim/BBRSimulation/data/waveguides}

BBSIM_ENV=${BBSIM_ENV:-/home/rshenoy/BBRSim/bb_env.sh}
source "$BBSIM_ENV" || { echo "ERROR: sourcing $BBSIM_ENV failed" >&2; exit 1; }
unset "${!PYAEDT_@}"
set -eo pipefail
CODE=/home/rshenoy/BBRSim/Blackbody-Simulations
cd "$CODE"

echo "host=$(hostname) job=${SLURM_JOB_ID} candidate=${JOB} stem=${STEM} reference=${REF}"
# Refuse to compare against anything but the reference bytes (fails when the stem has no recorded hashes).
grep -F "${STEM}_Ephi=" "$CODE/configs/legacy_reference_500GHz.sha256" | (cd "$REF" && sha256sum -c -)

status=0
for e in 0 1; do
  echo "== Ephi=$e"
  python compare_hfss_exports.py "$JOB/${STEM}_Ephi=$e" "$REF/${STEM}_Ephi=$e" --json "$JOB/comparison_Ephi$e.json" || status=1
done
if [ "$status" -eq 0 ]; then verdict=PASS; else verdict=FAIL; fi
echo "COMPARISON ${verdict}: ${JOB}/comparison_Ephi0.json and comparison_Ephi1.json"
exit "$status"
