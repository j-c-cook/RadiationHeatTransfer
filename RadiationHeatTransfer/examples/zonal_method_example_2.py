# Jack C. Cook
# Monday, March 29, 2021

# An example of the zonal method using the oven

import numpy as np
import RadiationHeatTransfer as rht
from copy import deepcopy
from numpy.linalg import inv
from numpy.linalg import solve


def main():
    path_to_output = 'oven_output.txt'

    A, Fij, eps = rht.view3d.read_output(path_to_output)

    s, r = A.shape
    print('The shape of the area vector is: {}'.format(s))

    U = np.array([.15] * s)  # U-value, W/m^2.K to back side
    h = np.array([10.] * s)  # internal convection coef., W/m2.K
    # Heat input by surface W/m2, e.g. by embedded electric resistance
    # heating cable
    qf_in = np.array([200., 0., 0., 0., 0., 0.])

    U = U[:, np.newaxis]
    h = h[:, np.newaxis]
    qf_in = qf_in[:, np.newaxis]

    TOS = 300.  # back side temperature

    Ts = np.array([300.] * s)  # Initial guess of surface temperatures
    Ts = Ts[:, np.newaxis]

    T = rht.enclosures.T_matrix(Fij, eps, A)
    S = rht.enclosures.S_matrix(Fij, eps, A)

    # T SS = S, Ax=b,  SS = inv(T) S
    SS = np.linalg.solve(T, S)

    sigma = 5.67e-08

    diff = 99999999
    count = 0

    # Begin while loop here
    while diff > 1.0E-06:

        Tsold = deepcopy(Ts)  # store Ts for convergence check

        UA = np.zeros_like(SS)
        for i in range(s):
            for j in range(s):
                UA[i, j] = SS[i, j] * sigma * (Ts[i]**2 + Ts[j]**2) \
                           * (Ts[i] + Ts[j])

        # initialize terms used in filling matrices
        UAjsum = np.zeros_like(A)
        RHS = np.zeros((s+1, r))
        for j in range(s):
            for i in range(s):
                UAjsum[j] = UAjsum[j] + UA[i, j]
        # fill conductance array, rows 1-N, columns 1-N
        UCA = np.zeros((s+1, s+1))
        for i in range(s):
            for j in range(s):
                UCA[i, j] = UA[i, j]
        # Correct diagonal through N
        for j in range(s):
            UCA[j, j] = UCA[j, j] - h[j] * A[j] - U[j] * A[j] - UAjsum[j]
        # fill in the N+1th row and N+1th column
        SumhA = 0
        for j in range(s):
            UCA[s, j] = h[j] * A[j]  # row
            UCA[j, s] = h[j] * A[j]  # column
            SumhA += h[j] * A[j]
        UCA[s, s] = -SumhA
        # Fill in RHS array
        for j in range(s):
            RHS[j] = -U[j] * A[j] * TOS - qf_in[j] * A[j]
        RHS[s] = 0.

        Tsnext = solve(UCA, RHS)
        Ts_n = Tsnext[0:s, :]

        diff = np.sqrt(np.mean((Ts_n - Tsold)**2))  # rmse

        Ts = deepcopy(Ts_n)

        count += 1

    print('Number of iterations: {}'.format(count))

    print(Ts)


if __name__ == '__main__':
    main()
