""" 
Script to merge all the group dispersion picks for all the .csv files corresponding to each station pair into one big table.
Written by Genevieve Savard @UniGe (2023)
"""

import pandas as pd
import os
import glob
import numpy as np
import time
import sys

t0 = time.time()

# Input/Output
dispsearch = "/srv/beegfs/scratch/users/h/henrymi/campiglia_data/campiglia_STACK/STACK_phase_only_rma/dispersion_V2/*.csv"

# Define some filtering threshold
snr_nbG_thresh = 10. # Force SNR threshold for the filtered waveform at each period (narrowband Gaussian)
d_lambda_thresh = 1.0 # Retain only picks for a minimum distance/wavelength multiple

output_file = f"/srv/beegfs/scratch/users/h/henrymi/campiglia_data/postprocessing/picking_csv/phaseonly_rma_snr{snr_nbG_thresh}_dlambda{d_lambda_thresh}_3Jul.csv"


# Get file list
dispfiles = glob.glob(dispsearch)
print(f"Number of pick files (.csv): {len(dispfiles)} files")

# Merge all picks and do some quality control
dflist = []
nfiles = len(dispfiles)
for k, dfile in enumerate(dispfiles):
    if k % 1000 == 0: 
        print(f"{k+1}/{nfiles}")
    
    try:
        df = pd.read_csv(dfile) # Read dispersion picks
    except:
        print(f"WARNING: Error when reading file {dfile}. Is the file empty? Skipping.")
        continue
    df["stasrc"] = os.path.split(dfile)[1].split("_")[0] # source station
    df["starcv"] = os.path.split(dfile)[1].split("_")[1] # receiver station
    
    # Apply QC criteria
    df = df.loc[(df.snr_nbG > snr_nbG_thresh) & (df.ratio_d_lambda > d_lambda_thresh), :]
    
    if len(df.inst_period.values) > 0: # If there are picks left after filtering
        # print(dfile)
        dflist.append(df)

# Merge the pandas dataframes        
picks = pd.concat(dflist)
del dflist

# Output to new .csv file
picks.to_csv(output_file)
print(f"Output file written: {output_file}")
print(f"Time elapsed: {time.time()-t0:.1f} s")
