# Jack C. Cook
# Monday, March 29, 2021

# An example showing how to read the output of view3d into numpy arrays

import RadiationHeatTransfer as rht


def main():
    path_to_output = 'output.txt'

    A, Fij, eps = rht.view3d.read_output(path_to_output)

    # Compute the summation of each row to see if the values are equal to 1.
    summation_rule = Fij.sum(axis=1)
    print('The result of the summation rule: ')
    print(summation_rule)

    # Check the symmetry of A*Fij to ensure reciprocity exists
    A_Fij = Fij * A
    reciprocity_check = rht.view3d.check_symmetric(A_Fij)
    print('The result of checking symmetry in the A * F_ij matrix: ')
    print(reciprocity_check)


if __name__ == '__main__':
    main()
