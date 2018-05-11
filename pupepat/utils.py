from astropy.stats import median_absolute_deviation
import numpy as np
from astropy.table import Table
from astropy.io import fits, ascii
import os
import sep
import logging
from PyPDF2 import PdfFileReader, PdfFileWriter

logger = logging.getLogger('pupepat')


def offsets(x1, y1, x2, y2):
    return ((x2 - x1) ** 2.0 + (y2 - y1) ** 2.0) ** 0.5

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
        param_names = [param for param in best_fit_models[0].keys()]
        output_table = {param: [best_fit_model[param] for best_fit_model in best_fit_models]
                        for param in param_names}
        output_table['filename'] = [os.path.basename(input_filename)] * len(best_fit_models)
        output_table['sourceid'] = range(len(best_fit_models))
        for keyword in header_keywords:
            output_table[keyword] = [hdu[0].header[keyword]] * len(best_fit_models)
        output_table = Table(output_table)
    else:
        for i, best_fit_model in enumerate(best_fit_models):
            best_fit_model['filename'] = os.path.basename(input_filename)
            best_fit_model['sourceid'] = i
            for keyword in header_keywords:
                best_fit_model[keyword] = hdu[0].header[keyword]
            logger.info(best_fit_model)
            output_table.add_row(best_fit_model)

    output_table.write(output_filename, format='ascii', overwrite=True)
    return output_table


def make_cutout(data, x0, y0, r):
    return data[int(y0 - r):int(y0 + r + 1), int(x0 - r):int(x0 + r + 1)]


def cutout_coordinates(cutout, x0, y0):
    x, y = np.meshgrid(np.arange(cutout.shape[1]), np.arange(cutout.shape[0]))
    return np.sqrt((x - x0) ** 2.0 + (y - y0) ** 2.0)


def run_sep(data, mask_threshold=None):
    if mask_threshold is None:
        mask_threshold = 20.0 * np.sqrt(np.median(data)) + np.median(data)
    background = sep.Background(np.ascontiguousarray(data), mask=np.ascontiguousarray(data > mask_threshold))
    return sep.extract(data - background, 20.0, err=np.sqrt(data), minarea=1000, deblend_cont=1.0, filter_kernel=None)


def prettify_focdmd(demanded_focus):
    return '{focdmd:+.0f}'.format(focdmd=demanded_focus).replace('-', 'm').replace('+', 'p')


def merge_pdfs(output_directory, output_table, output_pdf='pupe-pat'):
    data = ascii.read(os.path.join(output_directory, output_table))
    data.sort(['FOCDMD', 'filename'])
    pdf_files = np.array([os.path.join(output_directory,
                                       '{basename}_{sourceid}.pdf'.format(basename=os.path.splitext(row['filename'])[0],
                                                                          sourceid=row['sourceid']))
                          for row in data])
    for demanded_focus in np.unique(data['FOCDMD']):
        focus_set_indexes = data['FOCDMD'] == demanded_focus

        pdf_writer = PdfFileWriter()
        for pdf_file in pdf_files[focus_set_indexes]:
            pdf_writer.appendPagesFromReader(PdfFileReader(pdf_file))

        output_filename = os.path.join(output_directory,
                                       '{basename}-{focdmd}-donut.pdf'.format(basename=output_pdf,
                                                                              focdmd=prettify_focdmd(demanded_focus)))
        with open(output_filename, 'wb') as output_stream:
            pdf_writer.write(output_stream)


def estimate_scatter(a, axis=None):
    """
    Caclulate the inter-quartile range of *a* (scaled to normal).
    :param a: samples for the distribution
    :param axis: axis to calculate the scatter
    :return: scatter
    """
    # 2 Φ^−1(0.75)σ ≈ 1.349σ ≈ (27/20)σ for the interquartile range
    # 0.741301109 = 1.0 / (2**0.5 * special.erfinv(1.5 - 1.0) * 2.0)
    lower_quartile, upper_quartile = np.percentile(a, [25, 75], axis=axis)
    return 0.741301109 * (upper_quartile - lower_quartile)


def get_inliers(data, threshold):
    """
    Get the indexes of data that are within the threshold of the median
    :param data: Data to check for inliers
    :param threshold: n sigma threshold to consider a sample an inlier
    :return: Boolean array with True for inliers
    """
    scatter = estimate_scatter(data)
    offsets = np.abs(data - np.median(data)) / scatter
    return  offsets <= threshold


def image_is_defocused(filename):
    header = fits.getheader(filename)
    return header['FOCDMD'] != 0.0
