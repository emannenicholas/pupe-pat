import numpy as np


class SurfaceFitter(object):

    def _get_powers(self, degree):
        """
        Make a list of exponents for the polynomial terms
        :param degree: Maximum of the sum of the x and y powers
        :return: list of x,y powers to use for the polynomial
        """
        x_powers, y_powers = np.meshgrid(range(degree + 1), range(degree + 1))
        powers_to_use = (x_powers + y_powers <= degree)
        return zip(x_powers[powers_to_use], y_powers[powers_to_use])

    def fit(self, x, y, z, degree):
        """
        Fit a surface
        :param x: x coordinate grid for the surface
        :param y: y coordinate grid for the surface
        :param z: 2D surface data
        :param degree: maximum of the sum of x and y powers for the polynomial
        :return: Least squares coefficients (best fit parameters)
        """
        exponents = self._get_powers(degree)
        design_matrix = np.zeros((x.size, len(exponents)), dtype=np.float)
        for k, (i, j) in enumerate(exponents):
            design_matrix[:, k] = x.flatten() ** i * y.flatten() ** j
        coefficients, _, _, _ = np.linalg.lstsq(design_matrix, z.flatten())
        return coefficients

    def eval(self, x, y, coefficients, degree):
        """
        Evaluate a surface given coefficients
        :param x: x coordinate grid for the surface
        :param y: y coordinate grid for the surface
        :param coefficients: coefficients of the polynomial
        :param degree: maximum of the sum of x and y powers for the polynomial
        :return:
        """
        exponents = self._get_powers(degree)
        z = np.zeros(x.shape)
        for a, (i, j) in zip(coefficients, exponents):
            z += a * x ** i * y ** j
        return z
