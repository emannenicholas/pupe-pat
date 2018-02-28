from matplotlib import pyplot
from pupepat.ellipse import generate_ellipse

def plot_best_fit_ellipse(output_filename, data_cutout, best_fit_model):
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
    pyplot.savefig(output_filename)
