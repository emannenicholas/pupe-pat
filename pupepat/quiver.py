"""
pupepat.quiver - Tools to make a quiver plot to visualize the best M2 tip and tilt

Author
    Curtis McCully (cmccully@lco.global)

License
    GPL v3.0

February 2018
"""
import os
from pupepat.surface import SurfaceFitter
from pupepat.utils import get_inliers, offsets, estimate_scatter, prettify_focdmd
from pupepat.plot import plot_quiver
import numpy as np
from astropy.io import ascii


def make_quiver_plot(output_dir, output_table, output_plot='pupe-pat'):
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

        # If there are less than 5 images taken, don't bother making a quiver plot. The quiver plot requires
        # 2x2 free parameters so the fit will fail for less than 5.
        if len(focus_set) < 5:
            continue

        # Fit a smooth surface to inner and outer offsets
        focus_surface = fit_focus_surface(focus_set)

        # Make the actual plot
        output_filename = os.path.join(output_dir, '{basename}-{focdmd}-quiver.pdf'.format(basename=output_plot,
                                                                                      focdmd=prettify_focdmd(demanded_focus)))
        plot_quiver(focus_set, focus_surface, output_filename)


def fit_focus_surface(data):
    # Fit a first degree polynomial in x and y
    dx_fitter = SurfaceFitter(1)
    dx_fitter.fit(data['M2ROLL'], data['M2PITCH'], data['x0_outer'] - data['x0_inner'])
    dy_fitter = SurfaceFitter(1)
    dy_fitter.fit(data['M2ROLL'], data['M2PITCH'], data['y0_outer'] - data['y0_inner'])

    return dx_fitter, dy_fitter


def get_neighbor_inliers(data, min_m2_offset=2.0, max_m2_offset=150.0, min_neighbors=10, threshold=6.0):

    is_outlier = np.zeros(len(data), dtype=bool)
    for i, row in enumerate(data):
        m2_offset = offsets(data['M2ROLL'], data['M2PITCH'], row['M2ROLL'], row['M2PITCH'])
        is_neighbor = np.logical_and(m2_offset >= min_m2_offset, m2_offset <= max_m2_offset)
        # Short circuit if there are not enough neighbors
        # Note, this was originally is_neighbor.sum() - 1 (to remove the current data point) < min_neighbors
        if is_neighbor.sum() <= min_neighbors:
            continue
        center_offsets = offsets(data['x0_inner'], data['y0_inner'], data['x0_outer'], data['y0_outer'])
        neighbor_scatter = estimate_scatter(center_offsets[is_neighbor])

        this_offset = offsets(row['x0_inner'], row['y0_inner'], row['x0_outer'], row['y0_outer'])

        if this_offset > threshold * neighbor_scatter:
            is_outlier[i] = True

    return np.logical_not(is_outlier)
