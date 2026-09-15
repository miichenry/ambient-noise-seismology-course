#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import logging
from glob import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
from pyproj import Geod
import geopy.distance
import platform

# Specialized seismic/inversion libraries
import bayesbay as bb
import seislib
from seislib.tomography import SeismicTomography
from bayesbay.discretization import Voronoi2D
from bayesbay.prior import UniformPrior

# =============================================================================
# LOGGING CONFIGURATION
# =============================================================================
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

# =============================================================================
# CONSTANTS & PARAMETERS
# =============================================================================
WORK_PATH = '/home/users/h/henrymi/jectpro/campiglia'
STATION_FILE = os.path.join(WORK_PATH, 'campiglia_station_v3.csv')
DISP_FOLDER = os.path.expanduser('~/scratch/campiglia_data/postprocessing/disp/')
OUTPUT_DIR = '/home/users/h/henrymi/jectpro/campiglia/100mcell_90stations'

# Inversion Hyperparameters
PERTURBATION_STD_INV = 0.2
NUMBER_CHAINS = 20
ITER_BURNING_PHASE = 50000 #100000
ITER_MAIN_PHASE = 150000 #100000
CELL_SIZE_METERS = 100
VORONOI_DIM_MIN = 50
VORONOI_DIM_MAX = 1500

apply_gaussian_smoothing = False
GEOD = Geod(ellps="WGS84")

def setup_directories():
    """Ensure output directories exist."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for folder in ['fig', 'results']:
        os.makedirs(os.path.join(OUTPUT_DIR, folder), exist_ok=True)

def create_source_receivers(disp_folder, stations_file, period, all_pairs=False):
    """Processes dispersion data to create source-receiver observation files."""
    logger.info(f"Creating source-receivers for T={period}s (all_pairs={all_pairs})")
    
    df = pd.read_csv(stations_file, dtype={'station': str})
    sta_names = df['station'].values
    sta_lats = df['latitude'].values
    sta_lons = df['longitude'].values
    sta_elevs = df['elevation'].values

    src_lat, src_lon, rcv_lat, rcv_lon, travel_time, vels = [], [], [], [], [], []

    if not all_pairs:
        files = glob(os.path.join(disp_folder, '*'))
        for f in files:
            data = np.loadtxt(f, skiprows=1) #sometimes row 1 is a single number = to the interstation distance
            if data.size == 0:
                print(f"DEBUG {os.path.basename(f)}: empty file"); continue
            if data.ndim == 1: data = data.reshape(1, -1)
            t_axis, vg_axis = data[:, 0], data[:, 1]


            if len(t_axis) != len(np.unique(t_axis)):
                unique_t, indices = np.unique(t_axis, return_inverse=True)
                vg_axis = np.array([np.mean(vg_axis[indices == i]) for i in range(len(unique_t))])
                t_axis = unique_t

            if np.min(t_axis) - 1e-6 <= period <= np.max(t_axis) + 1e-6:
                parts = os.path.basename(f).split('_')
                s1_name, s2_name = parts[2], parts[3]
                idx1 = np.where(sta_names == s1_name)[0]
                idx2 = np.where(sta_names == s2_name)[0]

                if len(idx1) > 0 and len(idx2) > 0:
                    i1, i2 = idx1[0], idx2[0]
                    vg_interp = np.interp(period, t_axis, vg_axis)
                    if vg_interp <= 1e-9: continue
                    _, _, dist2d = GEOD.inv(sta_lons[i1], sta_lats[i1], sta_lons[i2], sta_lats[i2])
                    dz = sta_elevs[i2] - sta_elevs[i1]
                    dist3d_km = np.sqrt(dist2d**2 + dz**2) / 1000.0

                    src_lat.append(sta_lats[i1]); src_lon.append(sta_lons[i1])
                    rcv_lat.append(sta_lats[i2]); rcv_lon.append(sta_lons[i2])
                    travel_time.append(1.0 / vg_interp)
                    vels.append(vg_interp)  
    else:
        print(f"DEBUG period {period} outside range!")
        # Implementation for all_pairs=True omitted for brevity
        pass

    output_data = np.column_stack((src_lat, src_lon, rcv_lat, rcv_lon, travel_time, vels))
    if len(vels) == 0:
        raise ValueError(f"No velocity data found for period {period}s")
    return output_data, np.min(vels), np.max(vels)

def plot_convergence(results):
    """
    Plots Log-Likelihood and Model Complexity (Number of Voronoi Cells) 
    across all chains to check for convergence.
    """
    fig, ax = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    
    # 1. Plot Log-Likelihood
    # 'd_obs.log_likelihood' is the standard naming convention in bayesbay results
    log_likeli = results['d_obs.log_likelihood']
    for i in range(log_likeli.shape[0]):
        ax[0].plot(log_likeli[i, :], alpha=0.5, label=f'Chain {i}' if i < 5 else None)
    
    ax[0].set_ylabel('Log-Likelihood')
    ax[0].set_title('Convergence: Log-Likelihood per Chain')
    ax[0].grid(True, alpha=0.3)

    # 2. Plot Number of Voronoi Cells (Complexity)
    n_cells = np.array(results['voronoi.n_dimensions'])                                                                                                      
    if n_cells.ndim == 1:                                                                                                                                    
        n_cells = n_cells[np.newaxis, :]  # treat as single chain
    for i in range(n_cells.shape[0]):
        ax[1].plot(n_cells[i, :], alpha=0.5)
    
    ax[1].set_ylabel('Number of Voronoi Cells')
    ax[1].set_xlabel('Saved Iterations (x SAVE_EVERY)')
    ax[1].set_title('Model Complexity: Number of Voronoi Cells')
    ax[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'convergence_diagnostics.png'))
    plt.show()

def run_inversion(obs_data, period, vmin, vmax):
    """Performs the Bayesian Tomographic Inversion."""
    logger.info(f"Starting Inversion for Period: {period}s")
    obs = obs_data
    s_lat, s_lon, r_lat, r_lon, slow_init = obs[:, 0], obs[:, 1], obs[:, 2], obs[:, 3], obs[:, 4]
    
    lon_min, lon_max = min(s_lon.min(), r_lon.min()) - 0.01, max(s_lon.max(), r_lon.max()) + 0.01
    lat_min, lat_max = min(s_lat.min(), r_lat.min()) - 0.01, max(s_lat.max(), r_lat.max()) + 0.01

    lat_mean = 0.5 * (lat_min + lat_max)
    deg_per_m_lat = 1.0 / 111320.0
    deg_per_m_lon = 1.0 / (111320.0 * np.cos(np.deg2rad(lat_mean)))
    cell_size_deg = 0.5 * (CELL_SIZE_METERS * (deg_per_m_lat + deg_per_m_lon))

    tomo = SeismicTomography(cell_size=cell_size_deg, lonmin=lon_min, lonmax=lon_max, latmin=lat_min, latmax=lat_max, regular_grid=True)
    tomo.data_coords = np.column_stack((s_lat, s_lon, r_lat, r_lon))
    tomo.compile_coefficients()

    jacobian = scipy.sparse.csr_matrix(tomo.A)
    grid_points = np.column_stack(tomo.grid.midpoints_lon_lat())

    vel_prior = UniformPrior('vel', vmin=vmin, vmax=vmax, perturb_std=PERTURBATION_STD_INV)
    voronoi = Voronoi2D(name='voronoi', vmin=[lon_min, lat_min], vmax=[lon_max, lat_max], perturb_std=0.001, n_dimensions_min=VORONOI_DIM_MIN, n_dimensions_max=VORONOI_DIM_MAX, parameters=[vel_prior], compute_kdtree=True)

    def forward(state):
        vor = state['voronoi']
        kdtree = vor.load_from_cache('kdtree')
        nearest = kdtree.query(grid_points)[1]
        interp_vel = vor.get_param_values('vel')[nearest]
        state.save_to_extra_storage('interp_vel', interp_vel)
        return jacobian @ (1.0 / interp_vel)

    target = bb.likelihood.Target('d_obs', slow_init, std_min=0, std_max=1.0, std_perturb_std=0.005, noise_is_correlated=False)
    log_likelihood = bb.likelihood.LogLikelihood(targets=target, fwd_functions=forward)
    
    inversion = bb.BayesianInversion(parameterization=bb.parameterization.Parameterization(voronoi), log_likelihood=log_likelihood, n_chains=NUMBER_CHAINS)
    inversion.run(n_iterations=ITER_MAIN_PHASE + ITER_BURNING_PHASE, burnin_iterations=ITER_BURNING_PHASE, save_every=500, verbose=False)

    results = inversion.get_results()
    print(f"Keys in results: {list(results.keys())}")

    # Plot Complexity
    plot_complexity_only(results, period)
    
    # Plot Data Fit (using slow_init as the observed_data)
    plot_data_fit(results, slow_init, period)
    save_results(results, grid_points, vmin, vmax, period, tomo.A)

    # Terminal Stats
    n_cells = np.array(results['voronoi.n_dimensions'])
    if n_cells.ndim == 1:
        n_cells = n_cells[np.newaxis, :]
    print(f"Mean cells at start: {np.mean(n_cells[:, 0]):.1f}")
    print(f"Mean cells at end:   {np.mean(n_cells[:, -1]):.1f}")
    print(f"Std dev of cells at end: {np.std(n_cells[:, -1]):.1f}")

def plot_complexity_only(results, period):
    # Convert list to numpy array to use .shape and indexing
    n_cells = np.array(results['voronoi.n_dimensions'])
    # Ensure 2D shape (n_chains, n_iterations)
    if n_cells.ndim == 1:
        n_cells = n_cells[np.newaxis, :]
    plt.figure(figsize=(10, 5))

    # n_cells has shape (n_chains, n_iterations)
    for i in range(n_cells.shape[0]):
        plt.plot(n_cells[i, :], alpha=0.5, label=f'Chain {i}' if i < 3 else None)
    
    plt.title(f'Model Complexity (T={period}s)')
    plt.xlabel('Saved Iterations')
    plt.ylabel('Number of Voronoi Cells')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.savefig(os.path.join(OUTPUT_DIR, f'fig/complexity_T{period:.2f}.png'), dpi=300)
    plt.close()

def plot_data_fit(results, observed_data, period):
    # Convert list to numpy array
    d_pred = np.array(results['d_obs.dpred'])
    
    # d_pred is (n_chains, n_obs) — mean across chains
    final_pred = np.mean(d_pred, axis=0)
    residuals = observed_data - final_pred
    
    plt.figure(figsize=(10, 5))
    plt.hist(residuals, bins=30, color='skyblue', edgecolor='black')
    plt.axvline(0, color='red', linestyle='--')
    plt.title(f'Final Data Residuals (T={period}s)')
    plt.xlabel('Residual (Seconds)')
    plt.ylabel('Frequency')
    plt.savefig(os.path.join(OUTPUT_DIR, f'fig/data_fit_T{period:.2f}.png'), dpi=300)
    plt.close()

def save_results(results, grid_points, vmin, vmax, period, A_mat):
    """Generates plots and saves masked data files."""
    inferred_vel = np.mean(results['interp_vel'], axis=0)
    
    # 1. CALCULATE MASK BASED ON RAY COVERAGE
    # A_mat.sum(axis=0) gives the total sensitivity/path length in each cell
    col_sums = np.array(A_mat.sum(axis=0)).ravel()
    raw_mask = col_sums > 0  # True where we have ray paths
    
    # 2. APPLY MASK TO RAW VELOCITY (1D)
    # This ensures areas without data are NaN in the .dat output
    masked_inferred_vel = np.where(raw_mask, inferred_vel, np.nan)

    # 3. INTERPOLATION FOR 2D PLOTTING
    xi = np.linspace(grid_points[:, 0].min(), grid_points[:, 0].max(), 1000)
    yi = np.linspace(grid_points[:, 1].min(), grid_points[:, 1].max(), 1000)
    xx, yy = np.meshgrid(xi, yi)
    
    # Standard deviation propagation for plotting
    d_pred_std = np.std(1.0 / np.array(results['d_obs.dpred']), axis=0)
    # Use np.divide to handle 0-coverage cells safely
    cell_std = np.divide((A_mat.T @ d_pred_std), col_sums, out=np.zeros_like(col_sums), where=col_sums!=0)
    
    # Interpolate both velocity and uncertainty
    zi = griddata(grid_points, masked_inferred_vel, (xx, yy), method='linear')
    std_grid = griddata(grid_points, cell_std, (xx, yy), method='linear')

    if apply_gaussian_smoothing:
        zi = gaussian_filter(zi, sigma=0.2)


    # 4. PLOTTING voronoi grid + interpolation
    plt.figure(figsize=(10, 8))
    plt.scatter(grid_points[:, 0], grid_points[:, 1], c=masked_inferred_vel, cmap='RdYlBu', s=10, vmin=vmin, vmax=vmax)
    plt.colorbar(label='Velocity [km/s]')
    plt.title(f'Mean Velocity (T={period}s)')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.axis('equal')
    plt.savefig(os.path.join(OUTPUT_DIR, f'fig/grid_T{period:.2f}_.png'), dpi=300)
    plt.close()

    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    im1 = ax[0].imshow(zi, origin='lower', extent=[xi.min(), xi.max(), yi.min(), yi.max()], vmin=vmin, vmax=vmax, cmap='RdYlBu', aspect='auto')
    ax[0].set_title(f'Inferred velocity [km/s] (T={period}s)')
    
    im2 = ax[1].imshow(std_grid, origin='lower', extent=[xi.min(), xi.max(), yi.min(), yi.max()], vmin=0, cmap='viridis', aspect='auto')
    ax[1].set_title('Measurement std (Uncertainty)')

    plt.colorbar(im1, ax=ax[0], label='Velocity [km/s]')
    plt.colorbar(im2, ax=ax[1], label='d_obs [std]')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f'fig/inter_T{period:.2f}_.png'), dpi=500)
    plt.close()
    
    ############### Hit Count of the Jacobian ##############
    # 1. CALCULATE HIT COUNT (Total ray length per cell)
    # A_mat is (n_observations x n_grid_cells)
    # Summing over axis 0 gives the total path length in each cell
    hit_count = np.array(A_mat.sum(axis=0)).ravel()

    # 2. INTERPOLATE HIT COUNT FOR PLOTTING
    hit_grid = griddata(grid_points, hit_count, (xx, yy), method='linear')

    # 3. PLOT THE HIT COUNT MAP
    plt.figure(figsize=(8, 6))
    im = plt.imshow(hit_grid, origin='lower', 
                    extent=[xi.min(), xi.max(), yi.min(), yi.max()], 
                    cmap='hot_r', aspect='auto')
    plt.colorbar(im, label='Total Ray Path Length [m]')
    
    # Plot the stations if you have their coordinates available
    # plt.scatter(sta_lons, sta_lats, c='blue', marker='^', label='Stations')
    
    plt.title(f'Ray Density / Hit Count (T={period}s)')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f'fig/hitcount_T{period:.2f}.png'))
    plt.close()

    # 5. RANDOM VORONOI TESSELLATIONS (6-panel)
    fig, axes = plt.subplots(2, 3, figsize=(10, 6))
    random_indexes = np.random.choice(range(len(results['voronoi.vel'])), size=6, replace=False)
    for ipanel, (ax, irandom) in enumerate(zip(axes.ravel(), random_indexes)):
        voronoi_sites = results['voronoi.discretization'][irandom]
        vel_sample = results['voronoi.vel'][irandom]
        ax, cbar = Voronoi2D.plot_tessellation(
            voronoi_sites,
            vel_sample,
            ax=ax,
            voronoi_sites_kwargs=dict(markersize=0)
        )
        ax.tick_params(labelleft=False, labelbottom=False)
        ax.set_xlabel('')
        ax.set_ylabel('')
        cbar.set_label('Velocity [km/s]')
        if ipanel in [0, 3]:
            ax.tick_params(labelleft=True)
            ax.set_ylabel('Latitude')
        if ipanel not in [2, 5]:
            cbar.set_ticklabels('')
            cbar.set_label('')
        if ipanel > 2:
            ax.set_xlabel('Longitude')
            ax.tick_params(labelbottom=True)
    plt.suptitle(f'Random Voronoi Realizations (T={period}s)')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f'fig/voronoi_T{period:.2f}_.png'), dpi=300)
    plt.close()

    # 6. NOISE STD HISTOGRAM
    plt.figure()
    plt.hist(results['d_obs.std'], density=True, bins=20, ec='w', zorder=100, label='Posterior')
    plt.xlabel('Noise standard deviation')
    plt.ylabel('Density')
    plt.title(f'Inferred Data Noise (T={period}s)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f'fig/noise_T{period:.2f}_.png'), dpi=300)
    plt.close()

    # 7. SAVE MASKED RAW DATA
    # We use 'masked_inferred_vel' which contains NaNs for empty cells
    raw_data = np.column_stack((grid_points[:, 0], grid_points[:, 1], masked_inferred_vel))
    raw_results_path = os.path.join(OUTPUT_DIR, f"results/cell{CELL_SIZE_METERS}_model_{period:.6f}_.dat")
    np.savetxt(raw_results_path, raw_data, fmt="%10.6f", header="Lon Lat Vel")

    logger.info(f"Masked raw data saved to {raw_results_path}")

# =============================================================================
# MAIN EXECUTION
# =============================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Seismic Tomography Inversion Script")
    parser.add_argument("period", type=float, help="The period (T) to process")
    args = parser.parse_args()

    setup_directories()
    try:
        obs_data, vmin, vmax = create_source_receivers(disp_folder=DISP_FOLDER, stations_file=STATION_FILE, period=args.period)
        run_inversion(obs_data, args.period, vmin, vmax)
    except Exception as e:
        logger.exception(f"An error occurred: {e}")
