from matplotlib import pyplot
from pupepat.ellipse import generate_ellipse
import numpy as np
import os
from pupepat.surface import SurfaceFitter


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


def plot_quiver(data, focus_surface, output_plot, dpi=200):
    pyplot.clf()
    pyplot.grid(True)
    pyplot.gcf().set_size_inches(11, 10, forward=True)
    fit_evaluator = SurfaceFitter()
    X, Y = np.meshgrid(1.1 * np.linspace(min(data['M2ROLL']), max(data['M2ROLL']), 10 * dpi),
                       1.1 * np.linspace(min(data['M2PITCH']), max(data['M2PITCH']), 10 * dpi))
    center_offset_surface = fit_evaluator.eval(X, Y, focus_surface[0], 1) ** 2.0
    center_offset_surface += fit_evaluator.eval(X, Y, focus_surface[1], 1) ** 2.0
    center_offset_surface **= 0.5
    pyplot.contourf(X, Y, center_offset_surface, 15)
    pyplot.quiver(data['M2ROLL'], data['M2PITCH'], data['x0_inner'] - data['x0_outer'],
                  data['y0_inner'] - data['y0_outer'], headlength=3, headaxislength=3.0)
    best_roll = X[np.argmin(center_offset_surface)]
    best_pitch = Y[np.argmin(center_offset_surface)]
    pyplot.scatter(best_roll, best_pitch, marker='X', s=100, c='r', lw=0,
                   label='ROLL={roll:+d}, PITCH={pitch:+d}'.format(roll=best_pitch, pitch=best_pitch))
    pyplot.title("Donut Decenter (Inner - Outer)")
    pyplot.xlabel("M2 Roll tilt (arcsec)")
    pyplot.ylabel("M2 Pitch tilt (arcsec)")
    pyplot.legend(loc='upper right', framealpha=0.9)
    pyplot.tight_layout()
    pyplot.savefig(output_plot)
