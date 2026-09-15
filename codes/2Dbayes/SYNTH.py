#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import logging
from glob import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.sparse import csr_matrix
from pyproj import Geod
from seislib.tomography import SeismicTomography
from bayesbay.discretization import Voronoi2D
from bayesbay.prior import UniformPrior
import bayesbay as bb

GEOD = Geod(ellps="WGS84")

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
WORK_PATH = '/home/users/h/henrymi/jectpro/vulcano/'
STATION_FILE = os.path.join(WORK_PATH, 'stations_vulcano.csv')
DISP_FOLDER = os.path.join(WORK_PATH, 'disp')
OUTPUT_DIR = os.path.join(WORK_PATH, '100mcell_synth_checkerboard_400')

# Inversion hyperparameters
VMIN_INV = 1.9
VMAX_INV = 2.2
PERTURBATION_STD_INV = 0.2
NUMBER_CHAINS = 10
ITER_BURNING_PHASE = 30_000
ITER_MAIN_PHASE = 20_000
CELL_SIZE_METERS = 100
VORONOI_DIM_MIN = 50
VORONOI_DIM_MAX = 1500

# Checkerboard parameters
CHECKER_REF_VEL = 2.0   # reference velocity [km/s]
CHECKER_ANOM_AMP = 0.1  # anomaly amplitude [km/s]
CHECKER_KX = 16         # number of alternations in x direction
CHECKER_KY = 26         # number of alternations in y direction
NOISE_STD_OBS = 0.0     # noise std added to synthetic data (0 = no noise)

BBOX_LAT_MIN = 38.395161892417796
BBOX_LAT_MAX = 38.41207025207962
BBOX_LON_MIN = 14.947286879484098
BBOX_LON_MAX = 14.958928854610402

def create_source_receivers(disp_folder, stations_file, period, output_filename, all_pairs=False):
    """Process dispersion data and write source-receiver geometry to a CSV file.

    With all_pairs=True: uses all N*(N-1) ordered station pairs.  The travel
    times written are placeholders because SYNTH.py replaces them with
    synthetic forward-modelled data anyway.

    With all_pairs=False: only includes pairs that have a measured dispersion
    file covering the requested period (mirrors 2D_v2.py behaviour).
    """
    df = pd.read_csv(stations_file, dtype={'station': str})
    sta_names = df['station'].values
    sta_lats = df['latitude'].values
    sta_lons = df['longitude'].values
    sta_elevs = df['elevation'].values if 'elevation' in df.columns else np.zeros(len(sta_lats))

    src_lat, src_lon, rcv_lat, rcv_lon, travel_time, vels = [], [], [], [], [], []

    if all_pairs:
        n = len(sta_names)
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                src_lat.append(sta_lats[i])
                src_lon.append(sta_lons[i])
                rcv_lat.append(sta_lats[j])
                rcv_lon.append(sta_lons[j])
                travel_time.append(1.0)  # placeholder; replaced by synthetic d_obs
                vels.append(1.0)
    else:
        for f in glob(os.path.join(disp_folder, '*')):
            data = np.loadtxt(f, skiprows=0)
            if data.size == 0:
                continue
            if data.ndim == 1:
                data = data.reshape(1, -1)
            t_axis, vg_axis = data[:, 0], data[:, 1]

            # Average duplicate period entries
            if len(t_axis) != len(np.unique(t_axis)):
                unique_t, indices = np.unique(t_axis, return_inverse=True)
                vg_axis = np.array([np.mean(vg_axis[indices == i]) for i in range(len(unique_t))])
                t_axis = unique_t

            if not (np.min(t_axis) - 1e-6 <= period <= np.max(t_axis) + 1e-6):
                continue

            parts = os.path.basename(f).split('_')
            s1_name, s2_name = parts[2], parts[3]
            idx1 = np.where(sta_names == s1_name)[0]
            idx2 = np.where(sta_names == s2_name)[0]
            if len(idx1) == 0 or len(idx2) == 0:
                continue

            i1, i2 = idx1[0], idx2[0]
            vg_interp = np.interp(period, t_axis, vg_axis)
            if vg_interp <= 1e-9:
                continue

            _, _, dist2d = GEOD.inv(sta_lons[i1], sta_lats[i1], sta_lons[i2], sta_lats[i2])
            dz = sta_elevs[i2] - sta_elevs[i1]
            dist3d_km = np.sqrt(dist2d**2 + dz**2) / 1000.0

            src_lat.append(sta_lats[i1]); src_lon.append(sta_lons[i1])
            rcv_lat.append(sta_lats[i2]); rcv_lon.append(sta_lons[i2])
            travel_time.append(dist3d_km / vg_interp)
            vels.append(vg_interp)

    if len(vels) == 0:
        raise ValueError(f"No station pairs found for period {period}s")

    output = np.column_stack((src_lat, src_lon, rcv_lat, rcv_lon, travel_time, vels))
    np.savetxt(output_filename, output, fmt="%12.6f",
               header="src_lat src_lon rcv_lat rcv_lon travel_time vel")
    return np.min(vels), np.max(vels)


def setup_directories():
    """Ensure output directories exist."""
    for folder in ['fig', 'results']:
        os.makedirs(os.path.join(OUTPUT_DIR, folder), exist_ok=True)


def load_geometry(period):
    """Load source-receiver geometry for a given period."""
    logger.info(f"Loading geometry for T={period}s")

    # Use period-specific filename to avoid collisions in parallel array jobs
    obs_file = os.path.join(OUTPUT_DIR, f"observations_{period:.2f}.csv")
    create_source_receivers(
        disp_folder=DISP_FOLDER,
        stations_file=STATION_FILE,
        period=period,
        output_filename=obs_file,
        all_pairs=True
    )

    sta_lat = np.loadtxt(STATION_FILE, usecols=1, skiprows=1, delimiter=',')
    sta_lon = np.loadtxt(STATION_FILE, usecols=2, skiprows=1, delimiter=',')

    src_lat = np.loadtxt(obs_file, usecols=0, comments='#')
    src_lon = np.loadtxt(obs_file, usecols=1, comments='#')
    rcv_lat = np.loadtxt(obs_file, usecols=2, comments='#')
    rcv_lon = np.loadtxt(obs_file, usecols=3, comments='#')

    logger.info(f"Number of observations: {len(src_lat)}")
    return src_lat, src_lon, rcv_lat, rcv_lon, sta_lat, sta_lon


def build_tomo_grid(src_lat, src_lon, rcv_lat, rcv_lon):
    """Build SeismicTomography grid from source-receiver geometry."""
    lon_min = min(np.min(src_lon), np.min(rcv_lon)) - 0.01
    lon_max = max(np.max(src_lon), np.max(rcv_lon)) + 0.01
    lat_min = min(np.min(src_lat), np.min(rcv_lat)) - 0.01
    lat_max = max(np.max(src_lat), np.max(rcv_lat)) + 0.01

    # Convert cell size from meters to approximate degrees
    cell_size_deg_lat = (CELL_SIZE_METERS / 1000.0) / 111.32
    lat_mean = 0.5 * (lat_min + lat_max)
    meters_per_deg_lon = 111320.0 * np.cos(np.deg2rad(lat_mean))
    cell_size_deg_lon = CELL_SIZE_METERS / meters_per_deg_lon
    cell_size_deg = 0.5 * (cell_size_deg_lat + cell_size_deg_lon)

    logger.info(f"Cell size: {CELL_SIZE_METERS} m -> {cell_size_deg:.6f} deg")

    tomo = SeismicTomography(
        cell_size=cell_size_deg,
        lonmin=lon_min,
        lonmax=lon_max,
        latmin=lat_min,
        latmax=lat_max,
        regular_grid=True
    )

    nobs = len(src_lat)
    data_coords = np.zeros((nobs, 4))
    data_coords[:, 0] = src_lat
    data_coords[:, 1] = src_lon
    data_coords[:, 2] = rcv_lat
    data_coords[:, 3] = rcv_lon
    tomo.data_coords = data_coords

    logger.info(
        f"Grid: lon [{tomo.grid.lonmin:.4f}, {tomo.grid.lonmax:.4f}], "
        f"lat [{tomo.grid.latmin:.4f}, {tomo.grid.latmax:.4f}], "
        f"cells: {np.column_stack(tomo.grid.midpoints_lon_lat()).shape[0]}"
    )
    return tomo


def generate_checkerboard_data(tomo):
    """Generate synthetic checkerboard model and compute forward travel times."""
    grid_points = np.column_stack(tomo.grid.midpoints_lon_lat())

    logger.info("Compiling Jacobian matrix (A)...")
    tomo.compile_coefficients()
    A = csr_matrix(tomo.A)
    logger.info(f"Jacobian shape: {A.shape}, nnz = {A.nnz}")
    if A.nnz == 0:
        raise RuntimeError("Empty Jacobian: check cell_size or geometry.")

    vel_true_fun = tomo.checkerboard(
        ref_value=CHECKER_REF_VEL,
        kx=CHECKER_KX,
        ky=CHECKER_KY,
        lonmin=tomo.grid.lonmin,
        lonmax=tomo.grid.lonmax,
        latmin=tomo.grid.latmin,
        latmax=tomo.grid.latmax,
        anom_amp=CHECKER_ANOM_AMP
    )
    vel_true = vel_true_fun(grid_points[:, 0], grid_points[:, 1])

    # Synthetic data: d_obs = A @ slowness_true + noise
    d_obs_clean = A @ (1.0 / vel_true)
    noise = np.random.normal(0, NOISE_STD_OBS, size=d_obs_clean.shape)
    d_obs = d_obs_clean + noise
    logger.info(f"Synthetic checkerboard data generated (noise std = {NOISE_STD_OBS})")

    return A, grid_points, vel_true, d_obs


def run_inversion(A, grid_points, d_obs, period):
    """Run the Bayesian tomographic inversion."""
    logger.info(f"Starting Bayesian inversion for T={period}s")

    lon_min, lon_max = grid_points[:, 0].min(), grid_points[:, 0].max()
    lat_min, lat_max = grid_points[:, 1].min(), grid_points[:, 1].max()

    vel_prior = UniformPrior('vel', vmin=VMIN_INV, vmax=VMAX_INV, perturb_std=PERTURBATION_STD_INV)
    voronoi = Voronoi2D(
        name='voronoi',
        vmin=[lon_min, lat_min],
        vmax=[lon_max, lat_max],
        perturb_std=0.03,
        n_dimensions_min=VORONOI_DIM_MIN,
        n_dimensions_max=VORONOI_DIM_MAX,
        parameters=[vel_prior],
        compute_kdtree=True
    )

    def forward(state):
        vor = state['voronoi']
        kdtree = vor.load_from_cache('kdtree')
        nearest = kdtree.query(grid_points)[1]
        interp_vel = vor.get_param_values('vel')[nearest]
        state.save_to_extra_storage('interp_vel', interp_vel)
        return A @ (1.0 / interp_vel)

    target = bb.likelihood.Target(
        'd_obs', d_obs,
        std_min=0, std_max=0.1,
        std_perturb_std=0.001,
        noise_is_correlated=False
    )
    log_likelihood = bb.likelihood.LogLikelihood(targets=target, fwd_functions=forward)

    inversion = bb.BayesianInversion(
        parameterization=bb.parameterization.Parameterization(voronoi),
        log_likelihood=log_likelihood,
        n_chains=NUMBER_CHAINS
    )

    inversion.run(
        sampler=None,
        n_iterations=ITER_MAIN_PHASE + ITER_BURNING_PHASE,
        burnin_iterations=ITER_BURNING_PHASE,
        save_every=250,
        verbose=False,
        print_every=2000
    )

    return inversion.get_results()


def save_results(results, grid_points, vel_true, sta_lat, sta_lon, period):
    """Generate true vs. recovered comparison plots and save data files."""
    inferred_vel = np.mean(results['interp_vel'], axis=0)

    n_grid = 200
    xi = np.linspace(grid_points[:, 0].min(), grid_points[:, 0].max(), n_grid)
    yi = np.linspace(grid_points[:, 1].min(), grid_points[:, 1].max(), n_grid)
    XX, YY = np.meshgrid(xi, yi)

    true_grid = griddata(grid_points, vel_true, (XX, YY), method='linear')
    recovered_grid = griddata(grid_points, inferred_vel, (XX, YY), method='linear')

    vmin_plot = CHECKER_REF_VEL - CHECKER_ANOM_AMP
    vmax_plot = CHECKER_REF_VEL + CHECKER_ANOM_AMP

    def _make_checkerboard_fig():
        fig, axs = plt.subplots(1, 2, figsize=(12, 5))
        im0 = axs[0].imshow(
            true_grid, origin='lower',
            extent=(xi.min(), xi.max(), yi.min(), yi.max()),
            cmap='RdBu', vmin=vmin_plot, vmax=vmax_plot
        )
        axs[0].set_title('Checkerboard True Model')
        axs[0].plot(sta_lon, sta_lat, 'k.', markersize=2)
        plt.colorbar(im0, ax=axs[0], label='Velocity [km/s]')
        im1 = axs[1].imshow(
            recovered_grid, origin='lower',
            extent=(xi.min(), xi.max(), yi.min(), yi.max()),
            cmap='RdBu', vmin=vmin_plot, vmax=vmax_plot
        )
        axs[1].set_title('Recovered Model')
        axs[1].plot(sta_lon, sta_lat, 'k.', markersize=2)
        plt.colorbar(im1, ax=axs[1], label='Velocity [km/s]')
        plt.tight_layout()
        return fig, axs

    # Without bounding box
    fig, axs = _make_checkerboard_fig()
    fig_path = os.path.join(OUTPUT_DIR, f'fig/checkerboard_T{period:.2f}.pdf')
    plt.savefig(fig_path, dpi=300)
    plt.close()
    logger.info(f"Figure saved to {fig_path}")

    # With bounding box
    fig, axs = _make_checkerboard_fig()
    bbox_x = [BBOX_LON_MIN, BBOX_LON_MAX, BBOX_LON_MAX, BBOX_LON_MIN, BBOX_LON_MIN]
    bbox_y = [BBOX_LAT_MIN, BBOX_LAT_MIN, BBOX_LAT_MAX, BBOX_LAT_MAX, BBOX_LAT_MIN]
    for ax in axs:
        ax.plot(bbox_x, bbox_y, 'k-', linewidth=1.5)
    fig_bbox_path = os.path.join(OUTPUT_DIR, f'fig/checkerboard_T{period:.2f}_bbox.pdf')
    plt.savefig(fig_bbox_path, dpi=300)
    plt.close()
    logger.info(f"Figure with bbox saved to {fig_bbox_path}")

    # # Random Voronoi tessellations (6-panel)
    # fig, axes = plt.subplots(2, 3, figsize=(10, 6))
    # random_indexes = np.random.choice(range(len(results['voronoi.vel'])), size=6, replace=False)
    # for ipanel, (ax, irandom) in enumerate(zip(axes.ravel(), random_indexes)):
    #     voronoi_sites = results['voronoi.discretization'][irandom]
    #     vel_sample = results['voronoi.vel'][irandom]
    #     ax, cbar = Voronoi2D.plot_tessellation(
    #         voronoi_sites, vel_sample, ax=ax,
    #         voronoi_sites_kwargs=dict(markersize=0)
    #     )
    #     ax.tick_params(labelleft=False, labelbottom=False)
    #     ax.set_xlabel(''); ax.set_ylabel('')
    #     cbar.set_label('Velocity [km/s]')
    #     if ipanel in [0, 3]:
    #         ax.tick_params(labelleft=True); ax.set_ylabel('Latitude')
    #     if ipanel not in [2, 5]:
    #         cbar.set_ticklabels(''); cbar.set_label('')
    #     if ipanel > 2:
    #         ax.set_xlabel('Longitude'); ax.tick_params(labelbottom=True)
    # plt.suptitle(f'Random Voronoi Realizations (T={period}s)')
    # plt.tight_layout()
    # vor_path = os.path.join(OUTPUT_DIR, f'fig/voronoi_T{period:.2f}.png')
    # plt.savefig(vor_path, dpi=300)
    # plt.close()
    # logger.info(f"Voronoi figure saved to {vor_path}")

    data = np.column_stack((XX.ravel(), YY.ravel(), recovered_grid.ravel(), true_grid.ravel()))
    dat_path = os.path.join(OUTPUT_DIR, f'results/synth_model_{period:.2f}.dat')
    np.savetxt(dat_path, data, fmt="%10.5f", header="X Y recovered synthetic")
    logger.info(f"Results saved to {dat_path}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Checkerboard Synthetic Test - Bayesian Seismic Tomography")
    parser.add_argument("period", type=float, help="The period (T) to process")
    args = parser.parse_args()

    setup_directories()
    try:
        src_lat, src_lon, rcv_lat, rcv_lon, sta_lat, sta_lon = load_geometry(args.period)
        tomo = build_tomo_grid(src_lat, src_lon, rcv_lat, rcv_lon)
        A, grid_points, vel_true, d_obs = generate_checkerboard_data(tomo)
        results = run_inversion(A, grid_points, d_obs, args.period)
        save_results(results, grid_points, vel_true, sta_lat, sta_lon, args.period)
        logger.info("Checkerboard test complete.")
    except Exception as e:
        logger.exception(f"An error occurred: {e}")
