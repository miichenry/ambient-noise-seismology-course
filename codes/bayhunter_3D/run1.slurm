#!/bin/bash

#SBATCH --job-name=Bayhunter_Campiglia
#SBATCH --partition=shared-cpu,shared-bigmem,public-cpu,public-bigmem
#SBATCH --time=01:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30  # This is equal to the number of chains
#SBATCH --mem-per-cpu=30GB   # memory per cpu-core
#SBATCH --array=1-200 # from 1 to the number of lines in input_files.list
#SBATCH --output="outslurm_failedmodels/slurm-%A_%a_bayhunter.out"  # output path

#unset DISPLAY  # command to allow plots 

fname="input_files.list"  # File with for each line on different input file

dispfile=`cat $fname | sed -n "${SLURM_ARRAY_TASK_ID}p"`
echo "Processing this input file: $dispfile"

python bayhunter_1input.py $dispfile #adapt_bayhunter.py 
