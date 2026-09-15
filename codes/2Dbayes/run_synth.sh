#!/bin/bash
#SBATCH --job-name=synth_checker
#SBATCH --partition=public-cpu,public-bigmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=10G
#SBATCH --time=24:00:00
#SBATCH --output="out_synth/%x-%a.out"
#SBATCH --array=0-18
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=michail.henry@invert-geoscience.com

# Initialize Conda (ensure this path is correct for your system)
source ~/.bashrc
conda activate bayesbay

# Fix for Illegal Instruction (AVX/CPU mismatch)
# Forces libraries to use compatible instruction sets
export OPENBLAS_CORETYPE=Nehalem
export MKL_CBWR=COMPATIBLE

# Define the list of periods
plist=(0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0)

# Get the specific period for this task ID
current_period=${plist[$SLURM_ARRAY_TASK_ID]}

echo "Running array task $SLURM_ARRAY_TASK_ID for period $current_period"

mkdir -p out_synth

# Run the python script passing the period as an argument
python -u SYNTH.py $current_period
