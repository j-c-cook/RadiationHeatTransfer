# Jack C. Cook
# Monday, March 29, 2021

import numpy as np


def read_output(path_to_file, delimiter=' '):

    file = open(path_to_file, 'r+')
    data = file.read().split('\n')
    file.close()
    # remove blank lines
    data = [data[i] for i in range(len(data)) if data[i] != '']

    areas = data[1].split(delimiter)  # make a list of areas by splittingby the delimiter
    _view_factors = data[2:len(data)-1]
    view_factors = []
    for i in range(len(_view_factors)):  # loop through the list and build a 2D list by splitting by delimiter
        vf = _view_factors[i].split(delimiter)
        view_factors.append(vf)
    emissivities = data[-1].split(delimiter)

    areas = np.array(areas, dtype=np.float64)
    view_factors = np.array(view_factors, dtype=np.float64)
    emissivities = np.array(emissivities, dtype=np.float64)

    areas = areas[:, np.newaxis]  # correspond to the ith surface
    emissivities = emissivities[:, np.newaxis]

    return areas, view_factors, emissivities


def check_symmetric(a, rtol=1e-05, atol=1e-08):
    # Check the symmetry of a matrix
    # https://stackoverflow.com/a/42913743/11637415
    return np.allclose(a, a.T, rtol=rtol, atol=atol)
