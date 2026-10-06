#!/bin/bash -l
# Usage (HPC):  mkdir -p /home/rshenoy/BBRSim/outputs   # Slurm needs the --output directory to exist
#               sbatch hpc/run_frequency.sbatch configs/crack1_500GHz_reference.toml [extra run_hfss_frequency.py args]
# Memory and time are first guesses for the 15-angle reference run; adjust after the baseline solve.
# Shell pattern proven on this cluster (g4_full_check.sbatch, job 3996706): a login shell, and bb_env.sh
# sourced before strict mode, because module and the env scripts are not proven under set -u.
#SBATCH --job-name=bbsim_hfss
#SBATCH --account=golwala
#SBATCH --partition=expansion
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=/home/rshenoy/BBRSim/outputs/slurm-%j.out
CONFIG=${1:?config path required}
shift || true

BBSIM_ENV=${BBSIM_ENV:-/home/rshenoy/BBRSim/bb_env.sh}
source "$BBSIM_ENV" || { echo "ERROR: sourcing $BBSIM_ENV failed" >&2; exit 1; }
# Slurm exports the submitting shell: a stale PYAEDT_DESKTOP_PORT or PYAEDT_PROCESS_ID would make PyAEDT
# attach to another desktop instead of starting its own, and the runner refuses to start with them set.
unset "${!PYAEDT_@}"
set -eo pipefail
cd /home/rshenoy/BBRSim/Blackbody-Simulations

echo "host=$(hostname) job=${SLURM_JOB_ID} cpus=${SLURM_CPUS_PER_TASK} config=${CONFIG}"
python run_hfss_frequency.py --config "$CONFIG" --cores "$SLURM_CPUS_PER_TASK" "$@"
