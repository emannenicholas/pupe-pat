from astropy.stats import median_absolute_deviation
import numpy as np
from astropy.table import Table
from astropy.io import fits
import os
import sep
import logging
logger = logging.getLogger('pupepat')


def estimate_bias_level(hdu):
    """
    Estimate the bias level of image
    :param hdu: fits hdu with the image in question

    :return: float bias_level
    """
    # 1.48 here goes from median absolute deviation to standard deviation
    noise = 1.48 * median_absolute_deviation(hdu[0].data - sep.Background(hdu[0].data.astype(np.float)))
    gain = float(hdu[0].header['GAIN'])
    read_noise = float(hdu[0].header['RDNOISE'])
    bias_level = gain * np.median(hdu[0].data) - gain * gain * noise * noise + read_noise * read_noise
    bias_level /= gain
    return bias_level


def save_results(input_filename: str, best_fit_models: list,
                 output_table: Table, output_filename: str):
    """
    Save the best fit model to a text file and into the table

    :param input_filename: Filename of the image that was fit
    :param best_fit_models: list of astropy.modeling.Model for the best fit elliptical annulus
    :param output_table: astropy.tables.Table with the previous results
    :param output_filename: filename to save the results
    :return: output_table
    """
    header_keywords = ['M2ROLL', 'MOLTYPE', 'EXPTIME', 'FOCDMD', 'M2PITCH']
    hdu = fits.open(input_filename)

    if output_table is None:
        param_names = [param for param in best_fit_models[0].param_names]
        output_table = {param: [getattr(best_fit_model, param).value for best_fit_model in best_fit_models]
                        for param in param_names}
        output_table['filename'] = [os.path.basename(input_filename)] * len(best_fit_models)
        output_table['sourceid'] = range(len(best_fit_models))
        for keyword in header_keywords:
            output_table[keyword] = [hdu[0].header[keyword]] * len(best_fit_models)
        output_table = Table(output_table)
    else:
        for i, best_fit_model in enumerate(best_fit_models):
            best_fit_parameters = {param: getattr(best_fit_model, param).value for param in best_fit_model.param_names}
            best_fit_parameters['filename'] = input_filename
            best_fit_parameters['sourceid'] = i
            for keyword in header_keywords:
                best_fit_parameters[keyword] = hdu[0].header[keyword]
            logger.info(best_fit_parameters)
            output_table.add_row(best_fit_parameters)

    output_table.write(output_filename, format='ascii', overwrite=True)
    return output_table
