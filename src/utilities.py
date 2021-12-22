import numpy as np
from numpy import cos, arctan
from numpy.linalg import norm


def get_feature_func(ff_spec):
    ff_funcs = {'relative_orientation': lambda c1, c2: abs(cos((c1[-1] - c2[-1]) % np.pi)),
                'euclidean_distance': lambda c1, c2: norm(c1[:2] - c2[:2])}
    assert ff_spec in ff_funcs, f'Feature function specifier must be one of {list(ff_funcs.keys())}'
    return ff_funcs[ff_spec]


def get_nearest_neighbor_func(nn_spec):
    nn_funcs = {'euclidean_distance': lambda c1, c2: (c1[:2] - c2[:2]) @ (c1[:2] - c2[:2]),
                'annulus_distance': lambda c1, c2: abs(c1[:2] @ c1[:2] - c2[:2] @ c2[:2])
                                                   + abs(arctan(c1[0] / c1[1]) - arctan(c2[0] / c1[1]))}
    assert nn_spec in nn_funcs, f'NN function specifier must be one of {list(nn_funcs.keys())}'
    return nn_funcs[nn_spec]


if __name__ == "__main__":
    pass
