import numpy as np
from glob import glob
import os
from multiprocessing import Pool, cpu_count
from scipy.ndimage import gaussian_filter1d

#=================================================
#                   PARAMETERS
#=================================================
work_path = '/home/users/h/henrymi/jectpro/campiglia/100mcell_90stations'
data_file_path = os.path.join(work_path, 'results')
output_dir = os.path.join(work_path, "observed")
os.makedirs(output_dir, exist_ok=True)
bbox_raw_dir = os.path.join(work_path, "bbox_raw")
os.makedirs(bbox_raw_dir, exist_ok=True)

# Lat/lon bounding box — set to None to disable
LAT_MIN = None
LAT_MAX = None
LON_MIN = None
LON_MAX = None
print(f"Processing with bounding box:  {LAT_MIN} <= lat <= {LAT_MAX}, {LON_MIN} <= lon <= {LON_MAX}")  

files = sorted(glob(os.path.join(data_file_path, 'cell100_model_*')))
if len(files) == 0:
    raise RuntimeError(f"No raw_model_* files found in {data_file_path}")

#=================================================
# 1. ESTABLISH MASTER TEMPLATE FROM FIRST FILE
#=================================================
print(f"Establishing master grid template from: {os.path.basename(files[0])}")

# Using comments='#' to handle headers safely
ref_data = np.loadtxt(files[0], comments='#')
if ref_data.ndim == 1: ref_data = ref_data.reshape(1, -1)

lons_ref = ref_data[:, 0]
lats_ref = ref_data[:, 1]
n_index = len(lats_ref)

# Create lookup for the master 550 points using 6-decimal precision
coord_to_idx = {(round(lo, 6), round(la, 6)): i for i, (lo, la) in enumerate(zip(lons_ref, lats_ref))}

print(f"Master grid established with {n_index} points.")

#=================================================
# 2. LOAD AND MAP DATA (WITH NaN FILLING)
#=================================================
print("Loading files and mapping to master template...")

ALL_VELS_list, TIMES_list = [], []

for f in files:
    try:
        # Read header comment lines to preserve original format
        with open(f) as fh:
            header_lines = [ln for ln in fh if ln.startswith('#')]

        data = np.loadtxt(f, comments='#')
        if data.size == 0: continue
        if data.ndim == 1: data = data.reshape(1, -1)

        # Extract time from filename (assumes raw_model_TIME_.dat)
        t_val = float(os.path.basename(f).split('_')[2])

        # --- Write bbox-filtered copy in original format ---
        lons_f, lats_f = data[:, 0], data[:, 1]
        inbox_mask = np.ones(len(data), dtype=bool)
        if LAT_MIN is not None: inbox_mask &= lats_f >= LAT_MIN
        if LAT_MAX is not None: inbox_mask &= lats_f <= LAT_MAX
        if LON_MIN is not None: inbox_mask &= lons_f >= LON_MIN
        if LON_MAX is not None: inbox_mask &= lons_f <= LON_MAX
        bbox_data = data[inbox_mask]
        if bbox_data.size > 0:
            raw_out = os.path.join(bbox_raw_dir, os.path.basename(f))
            with open(raw_out, 'w') as fout:
                fout.writelines(header_lines)
            with open(raw_out, 'ab') as fout:
                np.savetxt(fout, bbox_data, fmt='%s')

        # Initialize a grid of NaNs based ONLY on the 550 master points
        full_vel_grid = np.full(n_index, np.nan)

        points_found = 0
        for row in data:
            lon, lat, vel = row[0], row[1], row[2]
            key = (round(lon, 6), round(lat, 6))

            # Only map if the point exists in the master template
            if key in coord_to_idx:
                full_vel_grid[coord_to_idx[key]] = vel
                points_found += 1

        ALL_VELS_list.append(full_vel_grid)
        TIMES_list.append(t_val)

        if points_found < n_index:
            print(f"File {os.path.basename(f)}: Only {points_found}/{n_index} points matched (rest filled with NaN).")
        else:
            print(f"File {os.path.basename(f)}: Full match ({points_found} points).")

    except Exception as e:
        print(f"ERROR: Could not process {f}: {e}")

ALL_VELS = np.array(ALL_VELS_list)
TIMES = np.array(TIMES_list)

# Indices where at least one file provided a value
has_data = ~np.all(np.isnan(ALL_VELS), axis=0)

# Mask points outside the lat/lon bounding box
in_box = np.ones(n_index, dtype=bool)
if LAT_MIN is not None: in_box &= lats_ref >= LAT_MIN
if LAT_MAX is not None: in_box &= lats_ref <= LAT_MAX
if LON_MIN is not None: in_box &= lons_ref >= LON_MIN
if LON_MAX is not None: in_box &= lons_ref <= LON_MAX

indexs_to_process = np.where(has_data & in_box)[0]
print(f"Bounding box filter: {in_box.sum()}/{n_index} points inside box, {len(indexs_to_process)} with data.")

#=================================================
# 3. WORKER FUNCTION & PROCESSING
#=================================================
def process_index(nx):
    vels = ALL_VELS[:, nx]
    mask = ~np.isnan(vels)

    # Need at least 2 points for smoothing
    if np.sum(mask) < 2: 
        return nx

    ts = TIMES[mask]
    vels_valid = vels[mask]

    # Sort chronologically
    sort_idx = np.argsort(ts)
    ts_sorted = ts[sort_idx]
    vels_sorted = vels_valid[sort_idx]

    # Apply Gaussian smoothing
    vg_gau = gaussian_filter1d(vels_sorted, 1)

    # Filename with 6 decimal places
    outname = os.path.join(output_dir, f"disp_{lats_ref[nx]:.6f}_{lons_ref[nx]:.6f}.dat")
    
    # Save content with 6 decimal places
    output = np.column_stack((ts_sorted, vg_gau))
    np.savetxt(outname, output, fmt='%.6f')

    return nx

if __name__ == "__main__":
    print(f"Processing {len(indexs_to_process)} dispersion curves using {cpu_count()} cores...")
    with Pool(cpu_count()) as pool:
        for _ in pool.imap_unordered(process_index, indexs_to_process):
            pass 
    print("Done!")
