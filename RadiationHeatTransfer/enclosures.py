# Jack C. Cook
# Sunday, February 28, 2021

import numpy as np


def F_matrix(Fij: np.ndarray, eps: np.ndarray):
    """
    Using the view factors and the epsilon values, compute the F_matrix

    F_i = 1 / eps_i + (eps_i - 1) / eps_1 * summation_{j=1}{N}(Fij)

    Parameters
    ----------
    Fij: np.ndarray
        View factor matrix
    eps: np.ndarray
        Emissive vector

    Returns
    -------

    """
    # each position in the view factor matrix is multiplied by (eps_i - 1) / eps_i
    eps_frac = (eps - 1) / eps
    # F = np.zeros_like(Fij)
    # for i in range(len(Fij)):
    #     F[:, i] = Fij[:, i] * eps_frac
    F = Fij * eps_frac
    # the diagonal (i==j) has an additional 1 / eps_i term
    eye = np.identity(4)
    eps_inv = 1 / eps
    eps_diag = eye * eps_inv
    F += eps_diag
    return F


def heat_flux(J: np.ndarray, Eb: np.ndarray, eps: np.ndarray):
    """
    After having found the radiosity leaving each surface, we can now find the
    heat flux.

    Parameters
    ----------
    J: np.ndarray
        The radiosity of each surface (matrix)
    Eb: np.ndarray
        The blackbody emissive power of each surface (column vector)
    eps: np.ndarray
        The emissivity of each surface (column vector)

    Returns
    -------

    """
    denomenator = (1 - eps) / eps
    Eb_i = Eb / denomenator
    J_i = J / denomenator
    flux = Eb_i - J_i
    return flux
