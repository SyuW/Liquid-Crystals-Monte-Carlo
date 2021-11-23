import argparse
import math
import os

import numpy as np

# custom imports
from LiquidCrystalSystem import LCSystem
from utilities import get_feature_func, get_nearest_neighbor_func


def load_dataset(dataset_path):

    # hard-coded values for dimensions, move to config file later
    R = 25
    r = 0
    b = 5
    a = 0.25

    # load a circular or annular confinement dataset
    if r == 0:
        particle_no_index = 3
        confinement = "Circle"
    else:
        particle_no_index = 4
        confinement = "Annulus"

    systems = dict()

    # iterate over datasets on path
    for _path_ in os.listdir(dataset_path):
        full_path = os.path.join(dataset_path, _path_, 'instanceRun')

        # Simulation summary notes exists
        if os.path.exists(os.path.join(full_path, f"MonteCarlo_{confinement}_SimNotes.txt")):
            lc = LCSystem(lc_data_path=full_path, confinement=confinement)
            num_of_particles = lc.sim_params['# of Ellipse']
            systems[num_of_particles] = lc

        # Simulation summary notes DNE
        else:
            lc = LCSystem(lc_data_path=full_path, confinement=confinement)
            lc.sim_params["R"] = R
            lc.sim_params["r"] = r
            lc.sim_params["Semi Major Axis"] = b
            lc.sim_params["Semi Minor Axis"] = a

            num_of_particles = lc.sim_params["# of Ellipse"]

            systems[num_of_particles] = lc

    return systems


def create_data_matrix(systems, num_of_features, num_of_samples, start=1000000, end=1500000):
    """

    :param systems:
    :param num_of_features:
    :param num_of_samples:
    :param start:
    :param end:
    :return:
    """

    # function for computing features
    feature_func = get_feature_func('orientation')
    data_matrix = []
    samples = dict()
    # iterate over densities
    for particle_number in systems.keys():
        samples[particle_number] = []
        system_state_at_mc_step = systems[particle_number].snapshots
        # iterate over Monte Carlo steps for chosen range
        for mc_step in system_state_at_mc_step.keys():
            if start <= mc_step <= end:
                # Get snapshot of system at Monte Carlo step
                snapshot = system_state_at_mc_step[mc_step]
                # Create the feature vectors
                fvs, _ = create_feature_vectors_from_snapshot(snapshot,
                                                              num_features=num_of_features,
                                                              num_samples=num_of_samples,
                                                              feature_func=feature_func)
                # Add to data matrix
                data_matrix = data_matrix + fvs
                for fv in fvs:
                    # subtract out the mean
                    # fv = fv - np.mean(fv)
                    samples[particle_number].append(fv)
        data_matrix = np.stack(data_matrix, axis=0)

    return data_matrix, samples


def create_feature_vectors_from_snapshot(coordinates, num_features, num_samples,
                                         feature_func=get_feature_func('relative_orientation'),
                                         nn_func=get_nearest_neighbor_func('euclidean_distance')):
    """
    Create the feature vectors from provided coordinates for all particles

    :param coordinates: list of (x, y, theta) coordinates
    :param num_features: number of features in each vector
    :param num_samples: total number of vectors to create
    :param feature_func: function for calculating features from coordinates
    :param nn_func: function for determining nearest neighbors

    :return: feature_vectors, feature_particle_coordinates
    """
    assert (num_features < len(coordinates)), \
        f"Number of features {num_features} cannot be greater than number of particles {len(coordinates)}"

    N = len(coordinates)

    # set the sampling rate for nearest neighbors
    if N % num_features == 0:
        nn_sampling_rate = N / num_features - 1
    else:
        nn_sampling_rate = math.floor(N / num_features)

    # set a random seed for reproducibility
    rng = np.random.default_rng(666)
    probe_indices = rng.choice(N, size=num_samples, replace=False)

    # initialize stuff
    feature_vectors = []
    neighbor_coords = dict()

    # select random probe particle for nearest neighbor distance
    for probe_index in probe_indices:

        # list of distances relative to probe particle
        probe_coord = coordinates[probe_index]
        neighbor_coords[tuple(probe_coord)] = []

        # sort based on nearest distance to probe
        nn_sorted = sorted(coordinates, key=nn_func)

        fv = []
        feature_particle_coords = [probe_coord]

        # add feature based on nearest neighbor distance
        for i, c in enumerate(nn_sorted):

            # don't add the probe particle
            if (i > 0) and (i % nn_sampling_rate) == 0:
                feature = feature_func(c, probe_coord)
                feature_particle_coords.append(c)
                fv.append(feature)

            # stop adding features if total number is met
            if len(fv) == num_features:
                break

        feature_vectors.append(fv)

    return feature_vectors, feature_particle_coords
