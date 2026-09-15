#!/bin/bash
#SBATCH --job-name=2Dbayes
#SBATCH --partition=public-cpu,public-bigmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=5G
#SBATCH --time=24:00:00
#SBATCH --output="out_vul/%x-%a.out"
#SBATCH --array=0-29 ##19
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=michail.henry@invert-geoscience.com

# Initialize Conda (ensure this path is correct for your system)
source ~/.bashrc  
conda activate bayesbay

# Fix for Illegal Instruction (AVX/CPU mismatch)
# Forces libraries to use compatible instruction sets
export OPENBLAS_CORETYPE=Nehalem
export MKL_CBWR=COMPATIBLE

# Ignore user-local packages (~/.local) to ensure we use the Conda environment's libraries
# NOTE: You must run 'sbatch fix_seislib.sh' once before running this array!
##export PYTHONNOUSERSITE=1

# Define the list of periods
plist=(0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0)
#plist=(0.2 1.0)
# Get the specific period for this task ID
# SLURM_ARRAY_TASK_ID will be 0, 1, 2... up to 19
current_period=${plist[$SLURM_ARRAY_TASK_ID]}

echo "Running array task $SLURM_ARRAY_TASK_ID for period $current_period"

# Run the python script passing the period as an argument
python -u 2D_v2.py $current_period
