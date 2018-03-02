from astropy.modeling import custom_model
import numpy as np
from astropy.io import fits
import sep
from pupepat.ellipse import inside_ellipse
from pupepat.utils import estimate_bias_level
from pupepat.plot import plot_best_fit_ellipse
from astropy.modeling import fitting
import os

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
    :param a_inner: Semimajor axis length of the inner boundary of the annulus
    :param b_inner: Semiminor axis length of the inner boundary of the annulus
    :param theta_inner: The rotation angle of the semimajor axis in radians of the inner boundary of the annulus
    :param amplitude_inner: Amplitude of the region inside the inner boundary of the annulus
    :param x0_outer: X center of the outer boundary of the annulus
    :param y0_outer: Y center of the outer boundary of the annulus
    :param a_outer:  Semimajor axis length of the outer boundary of the annulus
    :param b_outer:  Semiminor axis length of the outer boundary of the annulus
    :param theta_outer: The rotation angle of the semimajor axis in radians of the outer boundary of the annulus
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
    result += background
    return result


def fit_defocused_image(filename, plot_basename):
    logger.info('Extracting sources', extra={'tags':{'filename': os.path.basename(filename)}})
    hdu = fits.open(filename)
    data = float(hdu[0].header['GAIN']) * (hdu[0].data - estimate_bias_level(hdu))
    background = sep.Background(data)
    sources = sep.extract(data - background, 5.0, err=np.sqrt(data), minarea=5000, deblend_cont=1.0, filter_kernel=None)
    best_fit_models = [fit_cutout(data, source, plot_basename+'_{id}.pdf'.format(id=i), os.path.basename(filename))
                       for i, source in enumerate(sources)]
    return best_fit_models


def fit_cutout(data, source, plot_filename, image_filename):
    logger.info('Fitting source', extra={'tags': {'filename': image_filename, 'x': source['x'], 'y':source['y']}})
    cutout = data[source['y'] - 150:source['y'] + 151, source['x']- 150:source['x'] + 151]
    x, y = np.meshgrid(np.arange(cutout.shape[1]), np.arange(cutout.shape[0]))

    # Short circuit if either the source is in focus or if we just picked up a couple of hot columns
    got_bad_columns = (source['xmax'] - source['xmin']) < 100 or (source['ymax'] - source['ymin']) < 100
    background = np.median(data)
    in_focus = np.abs(data[int(source['y']), int(source['x'])] - background) > 200.0 * np.sqrt(background)
    if got_bad_columns or in_focus:
        if got_bad_columns:
            error_message = 'Star did not have enough columnns. Likely an artifact'
        elif in_focus:
            error_message = 'Star was not a donut.'
        logger.error(error_message, extra={'tags': {'filename': image_filename, 'x': source['x'], 'y':source['y']}})
        return

    x0 = source['x'] - source['xmin']
    y0 = source['y'] - source['ymin']
    r = np.sqrt((x - x0) ** 2.0 + (y - y0) ** 2.0)
    brightness_guess = np.median(cutout[r < 100])
    initial_model = elliptical_annulus(x0_inner=x0, y0_inner=y0,
                                       x0_outer=x0, y0_outer=y0,
                                       amplitude_inner=np.median(data), amplitude_outer=brightness_guess,
                                       a_inner=30, b_inner=30, a_outer=100, b_outer=100, background=np.median(data))

    fitter = fitting.SimplexLSQFitter()
    best_fit_model = fitter(initial_model, x, y, cutout, weights=1.0 / np.abs(cutout), maxiter=20000, acc=1e-6)
    plot_best_fit_ellipse(plot_filename, cutout, best_fit_model)

    logging_tags = {parameter: getattr(best_fit_model, parameter).value for parameter in best_fit_model.param_names}
    logging_tags['filename'] = image_filename
    for keyword in ['xmin', 'xmax', 'ymin', 'ymax']:
        logging_tags[keyword] = float(source[keyword])

    logger.info('Best fit parameters for PUPE-PAT model',
                extra={'tags': logging_tags})
    return best_fit_model
