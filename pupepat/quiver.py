import os
from pupepat.surface import SurfaceFitter
from pupepat.utils import get_inliers, offsets, estimate_scatter
from pupepat.plot import plot_quiver
import numpy as np
from astropy.io import ascii


def make_quiver_plot(output_dir, output_table, output_plot):
    """
    Make a quiver plot for a full M2 focus grid
    :param output_dir: Output Directory for the plot
    :param output_table: Name of the ascii file with fit results
    :param output_plot: Filename for the plot
    """
    # Read in the data table
    data = ascii.read(os.path.join(output_dir, output_table))
    data.sort('FOCDMD')

    for demanded_focus in np.unique(data['FOCDMD']):

        focus_set = data[data['FOCDMD'] == demanded_focus]
        # Calculate outliers using plain centroid and center of the circles
        centroid_inliers = get_inliers(offsets(focus_set['x0_centroid'], focus_set['y0_centroid'], focus_set['x0_outer'], focus_set['y0_outer']), 6.0)
        focus_set = focus_set[centroid_inliers]

        # Calculate outliers using nearby points
        neighbor_inliers = get_neighbor_inliers(focus_set)
        focus_set = focus_set[neighbor_inliers]

        # Fit a smooth surface to inner and outer offsets
        focus_surface = fit_focus_surface(focus_set)

        # Make the actual plot
        plot_quiver(focus_set, focus_surface, '{basename}_{focdmd: +d}.pdf'.format(basename=output_plot, focdmd=demanded_focus))


def fit_focus_surface(data):
    fitter = SurfaceFitter()
    dx_surface = fitter.fit(data['M2ROLL'], data['M2PITCH'], data['x0_outer'] - data['x0_inner'], 1)
    dy_surface = fitter.fit(data['M2ROLL'], data['M2PITCH'], data['y0_outer'] - data['y0_inner'], 1)

    return dx_surface, dy_surface


def get_neighbor_inliers(data, min_m2_offset=2.0, max_m2_offset=150.0, min_neighbors=10, threshold=6.0):

    is_outlier = np.zeros(len(data), dtype=bool)
    for i, row in enumerate(data):
        m2_offset = offsets(data['M2ROLL'], data['M2PITCH'], row['M2ROLL'], row['M2PITCH'])
        is_neighbor = np.logical_and(m2_offset >= min_m2_offset, m2_offset <= max_m2_offset)
        # Short circuit if there are not enough neighbors
        if is_neighbor.sum() - 1 < min_neighbors:
            continue
        center_offsets = offsets(data['x0_inner'], data['y0_inner'], data['x0_outer'], data['y0_outer'])
        neighbor_scatter = estimate_scatter(center_offsets[is_neighbor])

        this_offset = offsets(row['x0_inner'], row['y0_inner'], row['x0_outer'], row['y0_outer'])

        if this_offset > threshold * neighbor_scatter:
            is_outlier[i] = True

    return np.logical_not(is_outlier)
