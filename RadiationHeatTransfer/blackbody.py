# Jack C. Cook
# Sunday, January 31, 2021

"""
Blackbody Radiation:
This file will include all functions associated with black body radiation
"""

from math import pi, exp
from scipy.integrate import quad
import numpy as np
np.warnings.filterwarnings('ignore')
from scipy.interpolate import interp1d
from itertools import count
from itertools import takewhile
from scipy.integrate import cumtrapz


def Eblambda(lmbda: float, T: float, n: float = 1):
    """
    The relatiion for the spectral blackbody emissive power Eb,lambda was developed by Max Planck in
    1901. This formulation is also referred to as Planck's law.

    Equation 12-4 on page 720.

    E_{b\lambda}(\lambda,T) = \dfrac{C_1}{\lambda^5[exp(C_2/\lambda T)-1]}\;\;\;(W/m^2\cdot\mu m)

    :param T: the temperature of the surface (K)
    :param lmbda: the wavelength of the radiation emitted (micrometers)
    :param n: refractive index, assumed to be 1 (air)
    :return: the spectral black body emissive power at some temperature and wavelength (W/m^2/micrometer)
    """
    # Constants
    # c_0: the speed of light in a vacuum
    c_0: float = 2.9979 * 10**8  # (m/s)
    # k: Boltzmann's constant
    k: float = 1.38065 * 10**(-23)  # (J/K)
    # h: Planck's constant
    h: float = 6.626 * 10**(-34)

    # Convert micrometers to meters, for every one million micrometers there is one meter
    mm_to_m: float = 10 ** 6
    lmbda /= mm_to_m

    # Compute additional constants
    C_1: float = 2 * pi * h * c_0**2
    C_2: float = h * c_0 / k

    numerator: float = C_1

    # Theres some floating point issues with this function at some values
    denominator: float = lmbda**5 / n**2 * (np.exp(C_2/(n * lmbda * T))-1)

    return numerator / denominator / mm_to_m


def Eb(T: float):
    """
    Equation 12-6 on page 722.

    The total black body emissive power integrated over the whole range
    :param T: the temperature of the surface (K)
    :return: the total black body emissive power (W/m^2)
    """
    y, erf = quad(Eblambda, 0, np.inf, args=(T,))
    return y, erf


def Eb_vectorized(T: np.ndarray):
    """
    The Eb(T) function specifically takes in one temperature float value.
    To make use of the function by passing a numpy array, the function must be vectorized.


    Parameters
    ----------
    T: np.ndarray
        A numpy matrix of size (m, n)

    Returns
    -------
    The black body emissive power matrix that is the same size as the input T matrix
    """

    def call_Eb(Temp):
        y, erf = Eb(Temp)
        return y

    _Eb = np.vectorize(call_Eb)
    return _Eb(T)


def Eb_0_lambda(lmbda: float, T: float):
    """
    Integrate over some wavelength band.

    Equation 12-7

    :param lmbda: the wavelength of the radiation emitted (micrometers)
    :param T: the temperature of the surface (K)
    :return: the black body emissive power over some band (W/m^2)
    """
    return quad(Eblambda, 0, lmbda, args=(T,))


def f_lambda(lmbda: float, T: float):
    """
    The fraction of radiation emitted from a blackbody at temperature T
    in the wavelength band from lambda = 0 to lambda.

    Equation 12-8 on page 724

    :param lmbda: the wavelength of the radiation emitted (micrometers)
    :param T: the temperature of the surface (K)
    :return: The blackbody radiation function
    """
    y_num, _ = Eb_0_lambda(lmbda, T)
    y_den, _ = Eb(T)

    return y_num / y_den


def effective_spectral(lmbda_1: float, lmbda_2: float, T: float, f: interp1d, step: float=.01):
    """
    Find the effective or average
    Parameters
    ----------
    lmbda_1: float
        The $\lambda_1$ in micrometers
    lmbda_2: float
        The $\lambda_2$ in micrometers
    T: float
        Temperature of the surface of the blackbody
    f: interp1d
        A 1d scipy interpolation of the $\alpha$, $\rho$, $\tau$

    Returns
    -------
    The effective or average spectral hemispherical __
    """
    def any_range(start, stop, step):
        start = type(start + step)(start)
        return takewhile(lambda n: n < stop, count(start, step))

    y_den, _ = Eb(T)
    # lambda_1 may need to be 0.1
    wavelengths: list = list(any_range(lmbda_1, lmbda_2, step))
    eff_Eblambdas: list = []
    Eblambdas: list = []
    for _, lmbda in enumerate(wavelengths):
        Eblmda = Eblambda(lmbda, T)
        eff_Eblambdas.append(Eblmda * f(lmbda))
        Eblambdas.append(Eblmda)

    # integrate over y(x) using the composite trapezoidal rule
    y = cumtrapz(eff_Eblambdas, wavelengths)
    y_num = y.tolist()[-1]  # take the last value as the total integral

    eff = y_num / y_den

    return eff, wavelengths, eff_Eblambdas, Eblambdas
