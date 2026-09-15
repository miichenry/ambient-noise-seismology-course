#!/bin/bash
#SBATCH --job-name=2Dbayes
#SBATCH --partition=public-cpu,public-bigmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=10G
#SBATCH --time=12:00:00
#SBATCH --output="out/%x-%j.out"

# Initialize Conda (ensure this path is correct for your system)
source ~/.bashrc  
conda activate bayesbay

# Fix for Illegal Instruction (AVX/CPU mismatch)
# Forces libraries to use compatible instruction sets
export OPENBLAS_CORETYPE=Nehalem
export MKL_CBWR=COMPATIBLE

python -u create_observed_depth_inv.py
