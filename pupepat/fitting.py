"""
pupepat.fitting - Tools to fit models to pupil plate images

Author
    Curtis McCully (cmccully@lco.global)

License
    GPL v3.0

February 2018
"""
from astropy.modeling import custom_model
import numpy as np
from astropy.io import fits

from pupepat.ellipse import inside_ellipse
from pupepat.utils import estimate_bias_level, make_cutout, cutout_coordinates, run_sep
from pupepat.plot import plot_best_fit_ellipse
from astropy.modeling import fitting
import os
import emcee
import corner

import logging

logger = logging.getLogger('pupepat')


@custom_model
def elliptical_annulus(x, y, x0_inner=0.0, y0_inner=0.0, a_inner=1.0, b_inner=1.0, theta_inner=0.0, amplitude_inner=1.0,
                       x0_outer=0.0, y0_outer=0.0, a_outer=1.0, b_outer=1.0, theta_outer=0.0, amplitude_outer=1.0,
                       x_slope=0.0, y_slope=0.0, background=0.0):
    """
    2D Elliptical Annulus

    :param x: X positions to evaluate the model
    :param y: Y positions to evaluate the model
    :param x0_inner: X center of the inner boundary of the annulus
    :param y0_inner: Y center of the inner boundary of the annulus
    :param a_inner: Semi-major axis length of the inner boundary of the annulus
    :param b_inner: Semi-minor axis length of the inner boundary of the annulus
    :param theta_inner: The rotation angle of the semi-major axis in radians of the inner boundary of the annulus
    :param amplitude_inner: Amplitude of the region inside the inner boundary of the annulus
    :param x0_outer: X center of the outer boundary of the annulus
    :param y0_outer: Y center of the outer boundary of the annulus
    :param a_outer:  Semi-major axis length of the outer boundary of the annulus
    :param b_outer:  Semi-minor axis length of the outer boundary of the annulus
    :param theta_outer: The rotation angle of the semi-major axis in radians of the outer boundary of the annulus
    :param amplitude_outer: Amplitude of the region inside the annulus
    :param x_slope: Slope in the X-direction for the flux in the annulus
    :param y_slope: Slope in the Y-direction for the flux in the annulus
    :param background: Background level

    :return: numpy array with amplitude_inner inside the inner boundary of the annulus and amplitude_outer in the
             annulus
    """

    result = np.zeros(x.shape)
    # Include a gradient in the flux
    inside_outer = inside_ellipse(x, y, x0_outer, y0_outer, a_outer, b_outer, theta_outer)
    result[inside_outer] = amplitude_outer + x_slope * x[inside_outer] + y_slope * y[inside_outer]
    result[inside_ellipse(x, y, x0_inner, y0_inner, a_inner, b_inner, theta_inner)] = amplitude_inner
    result[~inside_outer] = background
    return result


def fit_defocused_image(filename, plot_basename):
    logger.info('Extracting sources', extra={'tags': {'filename': os.path.basename(filename)}})
    hdu = fits.open(filename)
    data = float(hdu[0].header['GAIN']) * (hdu[0].data - estimate_bias_level(hdu))
    sources = run_sep(data)
    best_fit_models = [fit_cutout(data, source, plot_basename, os.path.basename(filename),
                                  hdu[0].header, i)
                       for i, source in enumerate(sources)]
    return best_fit_models


def fit_cutout(data, source, plot_filename, image_filename, header, id, fit_circle=True, fit_gradient=False):
    logger.info('Fitting source', extra={'tags': {'filename': image_filename, 'x': source['x'], 'y': source['y']}})

    background = np.median(data)

    # Short circuit if either the source is in focus or if we just picked up a couple of hot columns
    got_bad_columns = (source['xmax'] - source['xmin']) < 200 and (source['ymax'] - source['ymin']) > 400
    got_bad_columns |= (source['xmax'] - source['xmin']) > 400 and (source['ymax'] - source['ymin']) < 200
    in_focus = np.abs(data[int(source['y']), int(source['x'])] - background) > 200.0 * np.sqrt(background)
    too_close_to_edge = source['x'] - 150 < 0 or source['y'] - 150 < 0
    too_close_to_edge |= source['x'] + 150 > data.shape[1] or source['y'] + 150 > data.shape[0]
    if got_bad_columns or in_focus or too_close_to_edge:
        if got_bad_columns:
            error_message = 'Star did not have enough columns. Likely an artifact'
        elif in_focus:
            error_message = 'Star was not a donut.'
        elif too_close_to_edge:
            error_message = 'Star was too close to the edge of the image.'
        logger.error(error_message, extra={'tags': {'filename': image_filename, 'x': source['x'], 'y': source['y']}})
        return

    cutout = make_cutout(data, source['x'], source['y'], 150)
    r = cutout_coordinates(cutout, 151, 151)

    inner_brightness_guess = np.median(cutout[r < 5])
    inside_donut_guess = cutout > (20.0 * np.sqrt(np.median(inner_brightness_guess)) + inner_brightness_guess)
    inner_radius_guess = np.percentile(r[inside_donut_guess].flatten(), 5)

    # Remake the cutout with +- 4 * inner_radius_guess
    cutout_radius = 4.0 * inner_radius_guess
    cutout = make_cutout(data, source['x'], source['y'], cutout_radius)
    x0 = y0 = cutout_radius + 1.0
    r = cutout_coordinates(cutout, x0, y0)

    # Run sep again and make sure there is only one source in the cutout
    cutout_sources = run_sep(cutout, mask_threshold=25.0 * np.sqrt(background) + background)

    if len(cutout_sources) > 1:
        logger.error('Too many sources detected in cutout. Likely source crowding.',
                     extra={'tags': {'filename': image_filename, 'x': source['x'], 'y': source['y']}})
        return

    inside_donut_guess = cutout > (20.0 * np.sqrt(np.median(inner_brightness_guess)) + background)
    outer_radius_guess = np.percentile(r[inside_donut_guess].flatten(), 90)

    brightness_guess = np.median(cutout[np.logical_and(r < outer_radius_guess, r > inner_radius_guess)])
    inner_brightness_guess = np.median(cutout[r < inner_radius_guess])

    initial_model = elliptical_annulus(x0_inner=x0, y0_inner=y0,
                                       x0_outer=x0, y0_outer=y0,
                                       amplitude_inner=inner_brightness_guess, amplitude_outer=brightness_guess,
                                       a_inner=inner_radius_guess, b_inner=inner_radius_guess,
                                       a_outer=outer_radius_guess,
                                       b_outer=outer_radius_guess,
                                       background=np.median(data))

    if not fit_gradient:
        initial_model.x_slope.fixed = True
        initial_model.y_slope.fixed = True

    if fit_circle:
        initial_model.theta_inner.fixed = True
        initial_model.theta_outer.fixed = True
        initial_model.b_inner.tied = lambda x: x.a_inner
        initial_model.b_outer.tied = lambda x: x.a_outer

    x, y = np.meshgrid(np.arange(cutout.shape[1]), np.arange(cutout.shape[0]))
    fitter = fitting.SimplexLSQFitter()
    best_fit_model = fitter(initial_model, x, y, cutout, weights=1.0 / np.abs(cutout), maxiter=20000, acc=1e-6)
    plot_best_fit_ellipse(plot_filename + '_{id}.pdf'.format(id=id), cutout, best_fit_model, header)

    logging_tags = {parameter: getattr(best_fit_model, parameter).value for parameter in best_fit_model.param_names}
    logging_tags['filename'] = image_filename
    for keyword in ['xmin', 'xmax', 'ymin', 'ymax']:
        logging_tags[keyword] = float(source[keyword])

    logger.info('Best fit parameters for PUPE-PAT model',
                extra={'tags': logging_tags})

    fit_results = {param: getattr(best_fit_model, param).value for param in best_fit_model.param_names}
    fit_results['x0_centroid'], fit_results['y0_centroid'] = x0, y0
    fit_results['sourceid'] = id
    return fit_results


def ln_likelihood_ellipse(theta, x, y, yerr):
    x0_inner, y0_inner, a_inner, b_inner, theta_inner, amplitude_inner, \
    x0_outer, y0_outer, a_outer, b_outer, theta_outer, amplitude_outer, \
    x_slope, y_slope, background = theta
    model = elliptical_annulus(x0_inner=x0_inner, y0_inner=y0_inner, a_inner=a_inner, b_inner=b_inner,
                               theta_inner=theta_inner,
                               amplitude_inner=amplitude_inner, x0_outer=x0_outer, y0_outer=y0_outer, a_outer=a_outer,
                               b_outer=b_outer, theta_outer=theta_outer, amplitude_outer=amplitude_outer,
                               x_slope=x_slope,
                               y_slope=y_slope, background=background)
    inv_sigma2 = yerr ** -2.0
    return -0.5 * (np.sum((y - model(x[0], x[1])) ** 2 * inv_sigma2 - np.log(inv_sigma2)))


def run_mcmc(x, y, yerr, theta0, likelihood_function, labels, output_plot='pupe-pat-mcmc.pdf', nwalkers=50):
    n_dim = len(theta0)
    starting_positions = [np.array(theta0) * (1.0 + np.random.uniform(-0.1, 0.1, len(theta0))) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, n_dim, likelihood_function, args=(x, y, yerr), threads=4)
    sampler.run_mcmc(starting_positions, 1000)
    samples = sampler.chain[:, 200:, :].reshape((-1, n_dim))
    fig = corner.corner(samples, labels=labels)
    fig.savefig(output_plot)


def ln_likelihood_circle(theta, x, y, yerr):
    x0_inner, y0_inner, r_inner, amplitude_inner, \
    x0_outer, y0_outer, r_outer, amplitude_outer, background = theta
    model = elliptical_annulus(x0_inner=x0_inner, y0_inner=y0_inner, a_inner=r_inner, b_inner=r_inner, theta_inner=0.0,
                               amplitude_inner=amplitude_inner, x0_outer=x0_outer, y0_outer=y0_outer, a_outer=r_outer,
                               b_outer=r_outer, theta_outer=0.0, amplitude_outer=amplitude_outer, x_slope=0.0,
                               y_slope=0.0, background=background)
    inv_sigma2 = yerr ** -2.0
    return -0.5 * (np.sum((y - model(x[0], x[1])) ** 2 * inv_sigma2 - np.log(inv_sigma2)))


ellipse_labels = ['x0_inner', 'y0_inner', 'a_inner', 'b_inner', 'theta_inner', 'amplitude_inner',
                  'x0_outer', 'y0_outer', 'a_outer', 'b_outer', 'theta_outer', 'amplitude_outer', 'x_slope', 'y_slope',
                  'background']

circle_labels = ['x0_inner', 'y0_inner', 'r_inner', 'amplitude_inner', 'x0_outer', 'y0_outer', 'r_outer',
                 'amplitude_outer', 'background']
