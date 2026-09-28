#!/bin/bash
# Usage (HPC):  mkdir -p /home/rshenoy/BBRSim/outputs   # Slurm needs the --output directory to exist
#               sbatch hpc/run_frequency.sbatch configs/crack1_500GHz_reference.toml [extra run_hfss_frequency.py args]
# Memory and time are first guesses for the 15-angle reference run; adjust after the baseline solve.
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
set -euo pipefail

CONFIG=${1:?config path required}
shift || true

source /home/rshenoy/BBRSim/bb_env.sh
cd /home/rshenoy/BBRSim/Blackbody-Simulations

echo "host=$(hostname) job=${SLURM_JOB_ID} cpus=${SLURM_CPUS_PER_TASK} config=${CONFIG}"
python run_hfss_frequency.py --config "$CONFIG" --cores "$SLURM_CPUS_PER_TASK" "$@"
