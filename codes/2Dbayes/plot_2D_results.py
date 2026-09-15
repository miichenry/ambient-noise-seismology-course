#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
from pathlib import Path

# =============================================================================
# PARAMETERS
# =============================================================================
results_dir   = Path('/home/users/h/henrymi/jectpro/campiglia/100mcell_90stations/results')
stations_file = Path('/home/users/h/henrymi/jectpro/campiglia/campiglia_station_v3.csv')
output_fig_dir = Path('/home/users/h/henrymi/jectpro/campiglia/100mcell_90stations/fig')
show = False

grid_size = 5000
interpolation_method = 'linear'

apply_gaussian_smoothing = False
gaussian_sigma_deg = 0.001  # smoothing radius in degrees (~500 m at Vulcano)

tile_zoom_level = 13
colormap = 'RdBu'
cbar_pad = 0.05
marker_shape = '^'
marker_color = 'darkgrey'
marker_edge_color = 'black'
alpha_value = 0.7
x_stretch = 2.0  # >1 makes X wider relative to geographic scale; 1 = true scale
_fig_h = 8.0
fig_size = (12,8)  # ~12 × 8 in
base_font_size = fig_size[0] * 2
marker_size = 10

output_fig_dir.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'font.size': base_font_size})

_TILE_URLS = {
    'esri_topo':  'https://server.arcgisonline.com/ArcGIS/rest/services/World_Topo_Map/MapServer/tile/{z}/{y}/{x}.jpg',
    'google_sat': 'https://mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z}',
    'google_hybrid': 'https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}',
}
tile_source = 'google_sat'  # 'esri_topo' | 'google_sat' | 'google_hybrid'

use_utm_gridlines = False
utm_zone = 33  # UTM zone for Vulcano
utm_x_ticks = list(range(494000, 501000, 1000))   # Easting  (m)
utm_y_ticks = list(range(4248000, 4257000, 1000))  # Northing (m)


def _parse_period(filename: str) -> float:
    for part in Path(filename).stem.split('_'):
        try:
            return float(part)
        except ValueError:
            continue
    raise ValueError(f"Cannot parse period from filename: {filename}")


def _load_dat(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r'\s+', comment='#', names=['lon', 'lat', 'vel'])
    return df.apply(pd.to_numeric, errors='coerce').dropna()


def _map_extent(station_df: pd.DataFrame, pad: float = 0.005):
    return (
        station_df['longitude'].min() - pad,
        station_df['longitude'].max() + pad,
        station_df['latitude'].min() - pad,
        station_df['latitude'].max() + pad,
    )


def plot_one_period(data: pd.DataFrame, period: float, station_df: pd.DataFrame,
                    ax, tile, vmin, vmax):
    min_lon, max_lon, min_lat, max_lat = _map_extent(station_df)

    ax.set_title(f"T = {period:.2f} s", fontsize=base_font_size, loc='right')
    ax.add_image(tile, tile_zoom_level)
    ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.Geodetic())
    ax.set_aspect('auto')

    lon_grid, lat_grid = np.mgrid[
        data['lon'].min():data['lon'].max():complex(0, grid_size),
        data['lat'].min():data['lat'].max():complex(0, grid_size),
    ]
    vel_grid = griddata(
        (data['lon'], data['lat']), data['vel'],
        (lon_grid, lat_grid), method=interpolation_method,
    )
    if apply_gaussian_smoothing:
        px_lon = (data['lon'].max() - data['lon'].min()) / grid_size
        px_lat = (data['lat'].max() - data['lat'].min()) / grid_size
        sigma_px = [gaussian_sigma_deg / px_lon, gaussian_sigma_deg / px_lat]
        nan_mask = np.isnan(vel_grid)
        vel_filled = np.where(nan_mask, np.nanmean(vel_grid), vel_grid)
        vel_grid = np.where(nan_mask, np.nan, gaussian_filter(vel_filled, sigma=sigma_px))

    im = ax.pcolormesh(
        lon_grid, lat_grid, vel_grid,
        cmap=colormap, shading='auto',
        transform=ccrs.PlateCarree(),
        vmin=vmin, vmax=vmax,
    )
    ax.scatter(
        station_df['longitude'], station_df['latitude'],
        marker=marker_shape, color=marker_color, edgecolor=marker_edge_color,
        s=marker_size, transform=ccrs.PlateCarree(), alpha=alpha_value, zorder=5,
    )
    ax.text(
        0.05, 0.05, f"Vmean = {data['vel'].mean():.3f} km/s",
        transform=ax.transAxes, fontsize=10,
        verticalalignment='bottom', horizontalalignment='left',
        bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.4'),
    )
    return im


def process_results(results_dir: Path):
    dat_files = sorted(results_dir.glob('cell*_model_*.dat'))
    if not dat_files:
        raise FileNotFoundError(f"No .dat files found in {results_dir}")

    station_df = pd.read_csv(stations_file)
    tile = cimgt.GoogleTiles(url=_TILE_URLS[tile_source])
    suffix = f'_{gaussian_sigma_deg}deg_smooth' if apply_gaussian_smoothing else ''

    for dat_path in dat_files:
        period = _parse_period(dat_path.name)
        data = _load_dat(dat_path)
        vmin, vmax = data['vel'].min(), data['vel'].max()
        print(f"{dat_path.name}  T={period:.4f}s  Vmin={vmin:.4f}  Vmax={vmax:.4f}")

        fig, ax = plt.subplots(1, 1, subplot_kw={'projection': tile.crs}, figsize=fig_size)
        im = plot_one_period(data, period, station_df, ax, tile, vmin, vmax)

        cbar = fig.colorbar(im, ax=ax, orientation='vertical', pad=cbar_pad, fraction=0.046)
        cbar.ax.ticklabel_format(useOffset=False, style='plain')
        cbar.set_label('Vel [km/s]', fontsize=8)

        if use_utm_gridlines:
            utm_crs = ccrs.UTM(zone=utm_zone)
            ax.gridlines(crs=utm_crs, draw_labels=False, alpha=0.2,
                         xlocs=utm_x_ticks, ylocs=utm_y_ticks)
            proj = ax.projection
            mid_n = np.mean(utm_y_ticks)
            mid_e = np.mean(utm_x_ticks)
            x_native = [proj.transform_point(x, mid_n, utm_crs)[0] for x in utm_x_ticks]
            y_native = [proj.transform_point(mid_e, y, utm_crs)[1] for y in utm_y_ticks]
            ax.set_xticks(x_native)
            ax.set_xticklabels([f'{x/1e5:.2f}' for x in utm_x_ticks], fontsize=7)
            ax.set_yticks(y_native)
            ax.set_yticklabels([f'{y/1e6:.3f}' for y in utm_y_ticks], fontsize=7)
            ax.set_xlabel('Easting (×10⁵ m)', fontsize=8)
            ax.set_ylabel('Northing (×10⁶ m)', fontsize=8)
        else:
            gl = ax.gridlines(draw_labels=True, alpha=0.2)
            gl.top_labels = False
            gl.right_labels = False

        out_path = output_fig_dir / f'2D_T{period:.4f}_{interpolation_method}{suffix}.png'
        plt.savefig(out_path, dpi=300, bbox_inches='tight')
        print(f"  → {out_path}")

        if show:
            plt.show()
        else:
            plt.close()


if __name__ == '__main__':
    process_results(results_dir)
