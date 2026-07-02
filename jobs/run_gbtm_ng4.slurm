#!/bin/bash
#SBATCH --job-name=q1_gbtm_biweekly_ng4_v3plus
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=/userhome/cs3/u3011656/hypertension/hypertension/analysis/outputs/slurm_logs/q1_gbtm_biweekly_ng4_v3plus_%j.log

set -euo pipefail

PROJECT_ROOT="/userhome/cs3/u3011656/hypertension/hypertension"
OUT_DIR="${PROJECT_ROOT}/analysis/outputs/q1_gbtm_relative_visit_only_biweekly_full_ng4_gpu_minvis3_parallel_2026-04-17"
RSCRIPT="/userhome/cs3/u3011656/miniconda3/envs/r_env/bin/Rscript"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

mkdir -p "${PROJECT_ROOT}/analysis/outputs/slurm_logs"
cd "${PROJECT_ROOT}"

echo "Starting full biweekly relative-time visit-only GBTM run at $(date)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-}"
echo "SLURM_JOB_PARTITION: ${SLURM_JOB_PARTITION:-}"
echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK:-1}"
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES:-}"
echo "OMP_NUM_THREADS: ${OMP_NUM_THREADS:-}"
echo "OPENBLAS_NUM_THREADS: ${OPENBLAS_NUM_THREADS:-}"
echo "MKL_NUM_THREADS: ${MKL_NUM_THREADS:-}"
nvidia-smi || true

"${RSCRIPT}" analysis/scripts/07i_q1_gbtm_relative_visit_only_shorttest.R \
  --base-dir "${PROJECT_ROOT}" \
  --out-dir "${OUT_DIR}" \
  --ng 4 \
  --bin-days 14 \
  --min-visits 3 \
  --sample-n 0 \
  --nstart 4 \
  --ncores "${SLURM_CPUS_PER_TASK:-1}" \
  --maxit 150

echo "Finished full biweekly relative-time visit-only GBTM run at $(date)"
