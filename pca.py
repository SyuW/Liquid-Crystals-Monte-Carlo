import argparse
import math
import os

import numpy as np
from numpy import array as arr

from sklearn.decomposition import PCA

# custom imports
from LiquidCrystalSystem import LCSystem


def create_feature_vectors_from_snapshot(coordinates, num_features, num_samples,
                                         feature_func=(lambda x, y: abs(np.cos(x - y)))):

    assert (num_features < len(coordinates)), \
        f"Number of features {num_features} cannot be greater than number of particles {len(coordinates)}"

    N = len(coordinates)
    feature_vectors = []

    # set the sampling rate for nearest neighbors
    if num_samples % num_features == 0:
        nn_sampling_rate = math.floor(N / num_features) - 1
    else:
        nn_sampling_rate = math.floor(N / num_features)

    # set a random seed for reproducibility
    rng = np.random.default_rng(666)
    probe_indices = rng.choice(N, size=num_samples, replace=False)

    # select random probe particle for nearest neighbor distance
    for probe_index in probe_indices:
        # list of distances relative to probe particle
        chosen_coord = coordinates[probe_index]

        # special functions for nearest neighbor sort and feature calculation
        norm_squared = lambda x, y: (arr(x) - arr(y)) @ (arr(x) - arr(y))
        distance_to_probe = lambda x: norm_squared(x[:2], chosen_coord[:2])

        # sort based on nearest distance to probe
        nn_sorted = sorted(coordinates, key=distance_to_probe)

        fv = []
        # add feature based on nearest neighbor distance
        for i, c in enumerate(nn_sorted):
            if (i > 0) and (i % nn_sampling_rate) == 0:
                feature = feature_func(c[-1], chosen_coord[-1])
                fv.append(feature)

            if len(fv) == num_features:
                break

        feature_vectors.append(fv)

    # concatenate into a matrix
    feature_vectors = np.stack(feature_vectors, axis=0)

    return feature_vectors


def create_data_matrix(lc_systems_to_density):

    number_of_features = 10
    datapoints_per_step = 5

    # number of rows AKA number of data-points
    # = number of densities * number of captures per density >= 1e6 * number of data-points per capture
    number_of_densities = len(lc_systems_to_density)
    number_of_captures = len([step for step in
                              lc_systems_to_density[0.2913752913752914].system_state_at_step
                              if step >= 1e6])

    print(f"Number of columns: {number_of_densities * number_of_captures * datapoints_per_step}")
    print(f"Number of features: {number_of_features}")

    data_matrix = []

    # iterate over densities
    for density in lc_systems_to_density.keys():

        system_state_over_steps = lc_systems_to_density[density].system_state_at_step

        # iterate over Monte Carlo steps
        for mc_step in system_state_over_steps:

            if mc_step >= 1e6:
                snapshot = system_state_over_steps[mc_step]

                feature_vecs = create_feature_vectors_from_snapshot(snapshot, number_of_features, datapoints_per_step)

                data_matrix.append(feature_vecs)

    data_matrix = np.stack(data_matrix, axis=0)

    return data_matrix


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="apply principal component analysis")
    parser.add_argument("--data-path", help="path to dataset")

    args = parser.parse_args()

    lc_systems_to_density = dict()
    for _path_ in os.listdir(args.data_path):
        full_path = os.path.join(args.data_path, _path_, "instanceRun")

        if os.path.exists(os.path.join(full_path, "MonteCarlo_Annulus_SimNotes.txt")):
            lc = LCSystem(lc_data_path=full_path)
            global_packing_fraction = lc.sim_params['reduced density']
        else:
            continue

        lc_systems_to_density[global_packing_fraction] = lc.system_state_at_step

    data = create_data_matrix(lc_systems_to_density)
    print(data.shape)
