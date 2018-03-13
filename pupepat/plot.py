from matplotlib import pyplot
from pupepat.ellipse import generate_ellipse
import os

def plot_best_fit_ellipse(output_filename, data_cutout, best_fit_model, header):
    """
    Plot the best fit elliptical annulus on top of the data.

    :param output_filename: Filename to save the plot
    :param data_cutout: numpy array with the cutout of the data used for the fit
    :param best_fit_model: astropy.modeling.Model corresponding to the best fit ellipse.
    """
    pyplot.clf()
    pyplot.imshow(data_cutout, origin='lower')
    x_inner, y_inner = generate_ellipse(best_fit_model.x0_inner, best_fit_model.y0_inner, best_fit_model.a_inner,
                                        best_fit_model.b_inner, best_fit_model.theta_inner)
    pyplot.plot(x_inner, y_inner, color='cyan')
    x_outer, y_outer = generate_ellipse(best_fit_model.x0_outer, best_fit_model.y0_outer, best_fit_model.a_outer,
                                        best_fit_model.b_outer, best_fit_model.theta_outer)
    pyplot.plot(x_outer, y_outer, color='cyan')
    pyplot.title(os.path.basename(os.path.splitext(output_filename)[0]))
    pyplot.text(data_cutout.shape[1] * 0.01, data_cutout.shape[0] * 0.95,
                'Inner: {x: 0.3f}, {y: 0.3f}'.format(x=best_fit_model.x0_inner.value, y=best_fit_model.y0_inner.value), color='white')
    pyplot.text(data_cutout.shape[1] * 0.01, data_cutout.shape[0] * 0.915,
                'Outer: {x: 0.3f}, {y: 0.3f}'.format(x=best_fit_model.x0_outer.value, y=best_fit_model.y0_outer.value), color='white')
    pyplot.text(data_cutout.shape[1] * 0.6, data_cutout.shape[0] * 0.95,
                'M2 Pitch: {pitch: 0.3f}"'.format(pitch=header['M2PITCH']), color='white')
    pyplot.text(data_cutout.shape[1] * 0.6, data_cutout.shape[0] * 0.915,
                'M2 Roll: {roll: 0.3f}"'.format(roll=header['M2ROLL']), color='white')
    pyplot.text(data_cutout.shape[1] * 0.6, data_cutout.shape[0] * 0.88,
                'FOCDMD: {focdmd: 0.3f} mm'.format(focdmd=header['FOCDMD']), color='white')
    pyplot.tight_layout()
    pyplot.savefig(output_filename)
