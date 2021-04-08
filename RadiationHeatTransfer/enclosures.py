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
    F = Fij * eps_frac
    # the diagonal (i==j) has an additional 1 / eps_i term
    m, n = eps.shape
    eye = np.identity(m)
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


def T_matrix(Fij, eps, A):
    # Calculate the reflectivity's for convenience
    rho = 1 - eps

    T = -rho.T * A * Fij / eps.T / A.T

    m, n = eps.shape
    eye = np.identity(m)
    eps_inv = 1 / eps
    eps_diag = eye * eps_inv
    T += eps_diag

    return T


def S_matrix(Fij, eps, A):
    return Fij * A * eps.T


def heat_flow(SS, eb):
    Q = np.zeros_like(eb)
    for i in range(4):
        for j in range(4):
            Q[i] += SS[i, j] * (eb[i] - eb[j])
    return Q
