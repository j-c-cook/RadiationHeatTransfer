# Jack C. Cook
# Monday, March 29, 2021

# An example of the zonal method

import numpy as np
import RadiationHeatTransfer as rht


def main():
    third = 1. / 3.
    Fij = np.array([[0., 0.25, .5, .25],
                    [third, 0., third, third],
                    [0.5, 0.25, 0., 0.25],
                    [third, third, third, 0.]], dtype=np.float64)

    # the following are column vectors
    # the areas of each surface are necessary
    A = np.array([0.4, 0.3, 0.4, 0.3])
    # make the shape of the area matrix be (n, 1)
    A = A[:, np.newaxis]  # correspond to the ith surface
    s = A.shape
    print('The shape of the area vector is: {}'.format(s))

    # The temperature of each surface, of shape (n, 1)
    Temp = np.array([1000., 600., 1000., 600.])
    Temp = Temp[:, np.newaxis]

    # the epsilon values for each surface, of shape (1, n)
    eps = np.array([0.3, 0.8, 0.3, 0.8])  # these correspond to the ith surface
    eps = eps[:, np.newaxis]

    T = rht.enclosures.T_matrix(Fij, eps, A)
    print('T: ')
    print(T)
    S = rht.enclosures.S_matrix(Fij, eps, A)
    print('S: ')
    print(S)

    # T SS = S, Ax=b,  SS = inv(T) S
    SS = np.linalg.solve(T, S)
    print('SS: ')
    print(SS)

    eb = rht.blackbody.Eb_vectorized(Temp)

    Q = rht.enclosures.heat_flow(SS, eb)

    qflux = Q / A
    print('Heat flux: ')
    print(qflux)



if __name__ == '__main__':
    main()
