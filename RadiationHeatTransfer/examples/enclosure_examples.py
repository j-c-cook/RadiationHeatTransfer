# Jack C. Cook
# Sunday, February 28, 2021

# A method for calling and testing the functions in the enclosures.py module

import RadiationHeatTransfer as RHT
import numpy as np


def main():
    # First we need the view factors, the view factor matrix is 2D
    # create a 2D list of view factors
    Fij = [[0., 0.25, 0.5, 0.25],
           [0.3333333, 0., 0.3333333, 0.3333333],
           [0.5, 0.25, 0., 0.25],
           [0.3333333, 0.3333333, 0.3333333, 0.]]
    # it will be easier to solve these using the linear algebra in numpy
    Fij = np.array(Fij)

    # the following are column vectors
    # the areas of each surface are necessary
    A = np.array([0.4, 0.3, 0.4, 0.3])
    # make the shape of the area matrix be (n, 1)
    A = A[:, np.newaxis]  # correspond to the ith surface
    s = A.shape
    print('The shape of the area vector is: {}'.format(s))

    # The temperature of each surface, of shape (n, 1)
    T = np.array([1000., 600., 1000., 600.])
    T = T[:, np.newaxis]

    # the epsilon values for each surface, of shape (1, n)
    eps = np.array([0.3, 0.8, 0.3, 0.8])  # these correspond to the ith surface
    eps = eps[:, np.newaxis]

    # Compute the black body emissive power of each surface
    # Note: to make use of the Eb(T) function in RHT.blackbody.py, we need to vectorize the function
    eb = RHT.blackbody.Eb_vectorized(T)
    s = eb.shape
    print('The shape of the eb vector should be the same as the temp: {}'.format(s))
    print('Eb: ')
    print(eb)

    F = RHT.enclosures.F_matrix(Fij, eps)
    print('F: ')
    print(F)

    # F J = Eb  --  Ax = b
    J = np.linalg.solve(F, eb)
    print('J: ')
    print(J)

    heat_flux = RHT.enclosures.heat_flux(J, eb, eps)
    print('heat_flux: ')
    print(heat_flux)

    Q = heat_flux * A
    print('Q: ')
    print(Q)

    a = 1


if __name__ == '__main__':
    main()
