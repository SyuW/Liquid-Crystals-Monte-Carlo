from argparse import Namespace
from datetime import date
from joblib import Parallel, delayed
import math
from multiprocessing import Pool
import numpy as np
import os

# custom imports
from utilities import get_feature_func, get_nearest_neighbor_func


def nematic_order_param(coord_list):
    """
    Calculate the nematic order parameter at a given area fraction
    :param coord_list:
    :return:
    """
    angles = coord_list[:, -1]
    avg_sines = np.mean(np.sin(angles))
    avg_coses = np.mean(np.cos(angles))

    return np.sqrt(avg_sines**2 + avg_coses**2)


