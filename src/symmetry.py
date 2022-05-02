import argparse
import os
import time

import numpy as np
import matplotlib.pyplot as plt

from datetime import date
import pickle
from sklearn.decomposition import PCA
from time import gmtime, strftime
import time

# custom imports
from utilities import get_feature_func, get_nearest_neighbor_func
from order_params import nematic_order_param


def polar_system(coordinates, R):
    num_polar_steps = 30
    num_radial_steps = 30
    angle_divs, angle_step = np.linspace(0, 2 * np.pi, num_polar_steps, retstep=True)
    radial_divs, radial_step = np.linspace(0, R, num_polar_steps, retstep=True)

    r_mesh, theta_mesh = np.meshgrid(radial_divs, angle_divs)

    x_coords = coordinates[:, 0]
    y_coords = coordinates[:, 1]

    # convert to polar
    rlist = np.sqrt(x_coords ** 2 + y_coords ** 2)
    thetalist = np.arctan2(x_coords / y_coords)

    polar_map = dict()

    fig, ax = plt.subplots(dpi=120, subplot_kw=dict(projection='polar'))
    ax.contourf()


def in_quadrant(angle, x, y):
    """
    Check that (x,y) and ray at angle are within same quadrant
    :param angle:
    :param x:
    :param y:
    :return:
    """
    if 0 <= angle <= np.pi / 2:
        quadrant = 1
    elif np.pi / 2 <= angle <= np.pi:
        quadrant = 2
    elif np.pi <= angle <= 3 * np.pi / 2:
        quadrant = 3
    elif 3 * np.pi / 2 <= angle <= 2 * np.pi:
        quadrant = 4
    else:
        raise ValueError("Angle out of bounds")
    if x >= 0 and y >= 0:
        coord_quadrant = 1
    elif x < 0 < y:
        coord_quadrant = 2
    elif x < 0 and y < 0:
        coord_quadrant = 3
    else:
        coord_quadrant = 4
    if coord_quadrant == quadrant:
        return True
    else:
        return False


def symmetry_group_detection(coordinates, outer_radius):
    """
    sweep a rotating radial band around 'snapshot' and
    bin points by angle

    :param coordinates:
    :param outer_radius:
    :return:
    """
    # parameters for number of angular divisions
    num_polar_steps = 100
    num_radial_steps = 100
    angle_divs, angle_step = np.linspace(0, 2 * np.pi, num_polar_steps, retstep=True)
    radial_divs, radial_step = np.linspace(0, outer_radius, num_polar_steps, retstep=True)
    width = 4.5

    # add points to each radial band
    radial_bands = dict()
    for theta in angle_divs:
        radial_bands[theta] = []
        # bounding lines of radial band
        line1 = lambda x: (np.sin(theta) * x + (width / 2)) / np.cos(theta)
        line2 = lambda x: (np.sin(theta) * x - (width / 2)) / np.cos(theta)
        for coord in coordinates:
            xc, yc = coord
            # check that xc, yc are in the same quadrant as the radial band
            if in_quadrant(theta, xc, yc):
                a = line1(xc)
                b = line2(xc)
                if b < a:
                    a = line2(xc)
                    b = line1(xc)
                if a <= yc <= b:
                    radial_bands[theta] += [np.sqrt(xc ** 2 + yc ** 2)]
        # check if empty
        if not radial_bands[theta]:
            radial_bands[theta] += [outer_radius]

    # histogram of angle vs points within band
    # counts = {t: len(radial_bands[t]) for t in radial_bands}
    # t_counts = list(counts.values())

    return radial_bands


if __name__ == "__main__":
    from features import load_dataset, create_data_matrix

    # start the timer
    start = time.perf_counter()

    sys_config = {"R": 30.1, "r": 0, "b": 5, "a": 0.25, "nf": 480, "ns": 15, "start": 200000, "end": 800000}

    project_path = "C:\\Users\\Sam Yu\\Desktop\\School\\4B\\Phys_437B_RP"
    data_path = os.path.join(project_path, f"datasets\\r={sys_config['r']}\\2_17_2022")
    current_time = strftime("%Y-%m-%d_%H-%M-%S", gmtime())
    results_path = f"results\\r={sys_config['r']}\\{current_time}_nf_{sys_config['nf']}_ns_{sys_config['ns']}" \
                   f"_start_{sys_config['start']}_end_{sys_config['end']}"
    results_path = os.path.join(project_path, results_path)
    cache_path = os.path.join(project_path, f"cache\\r={sys_config['r']}\\lcs.pickle")

    """
    # cache the file so that we don't to load the dataset again each time
    if not os.path.isfile(os.path.join(cache_path, "lcs.pickle")):
        print("Cached dataset file not found, making a new one...")
        lcs = load_dataset(dataset_path=data_path)
        with open(cache_path, "wb") as out_f:
            pickle.dump(lcs, out_f)
    else:
        print("Cached dataset found, loading...")
        with open(cache_path, "rb") as in_f:
            lcs = pickle.load(in_f)
    """

    test_datapoints = np.loadtxt(os.path.join(project_path, "datasets\\symmetry_data\\11fold.txt"), delimiter=" ")

    symmetry_group_detection(test_datapoints, outer_radius=5)

    print(f"Symmetry group analysis time: {time.perf_counter() - start} seconds")

    pass
