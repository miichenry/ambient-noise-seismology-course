import re
import glob
import os
 
path = os.path.expanduser("~/scratch/campiglia_data/postprocessing/disp/disp_ZZ_*_td.dat")
 
files = glob.glob(path)
 
stations = set()
for f in files:
    m = re.search(r'disp_ZZ_(\d+)_(\d+)_td\.dat', f)
    if m:
        stations.update(m.groups())
 
print(f"Unique stations: {len(stations)}")
print(sorted(stations))
 
