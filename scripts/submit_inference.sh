#!/bin/bash
#SBATCH --account=rrg-mkoz_gpu
#SBATCH --job-name=rnp_full
#SBATCH --time=30:00:00
#SBATCH --gres=gpu:h100:1
#SBATCH --mem=40G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/rnp_full_%j.out

module load python/3.11 rdkit/2024.09.6 cuda/12.6
source venv/bin/activate

protenix pred \
    --input /scratch/wvirany/boltz-UQ/inputs/rnp_full \
    --out_dir /scratch/wvirany/boltz-UQ/output/rnp_full \
    --use_default_params true \
    --use_msa true \
    --sample 1 \
    --trimul_kernel cuequivariance \
    --triatt_kernel cuequivariance