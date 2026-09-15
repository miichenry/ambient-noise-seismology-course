#!/bin/bash
# NB: When defining Slurm SBATCH parameters, there is only one "#" before "SBATCH", if there are 2, this line is commented.
#SBATCH --job-name=S2_ML_165_2  # create a name for your job
#SBATCH --partition=public-cpu,shared-cpu,shared-bigmem,public-bigmem #,public-cpu,public-bigmem # list of partitions to use (remember each partition has their own constraints!)
#SBATCH --ntasks=20              # total number of parallel tasks in MPI (increase for more speed, but see below)
#SBATCH --cpus-per-task=1        # cpu-cores per task (always choose 1 here)
###SBATCH --mem-per-cpu=20G               # Total memory requested, shared among all CPUs (divide mem 
#SBATCH --mem=250G                # by ntasks to get memory allocated per task) also option to use mem = 800G (total mem needed for all ntasks in 1 node)
#SBATCH --time=12:00:00          # total run time limit (HH:MM:SS) (keep in mind limits of partition requested)
#SBATCH --output="outslurm/slurm-%j-%x.out"   # Path where to write output files. Change as you want. %j is job ID, %x is job-name defined above.
##SBATCH --mail-type=END         # if you want to get an email when your job is done
##SBATCH --mail-user=elliot.jiwani-brown@unige.ch  # Your email

# Copy this slurm files and the param files in your project directory and modify them there. Adjust the path to the script below

# *** Choice of NTASKS (NPROCS): ***
# Step 2: Do not specify more than the number of station pairs.
# Step 1: Do not specify more than the number of .h5 files in RAW_DATA
# Step 0: Do not specify more than the is number of time chunks (i.e. total dataset duration in hours divided by inc_hours)

# *** Define the name of the input parameter file. ***
#PARAMS="S0B_params.yaml" #ntasks = 80, mem = 800G
#PARAMS="S1_params.yaml" #ntasks = 40, mem = 1000G
PARAMS="S2_params.yaml"

# *** Define the name of the script to use ***
#SCRIPT=S0B_to_ASDF.py   # Convert mseed dataset to H5 files in time chunks
#SCRIPT=S1_fft_cc_MPI.py # Get cross-correlations  
#SCRIPT=S1_fft_cc_MPI_normZ.py # Get cross-correlations, with 3C data normalized by Z component (preserves relative amplitudes)
SCRIPT=S2_stacking.py  # Stacking of the cross-correlations
#SCRIPT=S2_stacking_acorr.py # Stacking for autocorrelations only

*** Now call run the script with MPI! ***
mpiexec -n $SLURM_NPROCS python $SCRIPT $PARAMS
