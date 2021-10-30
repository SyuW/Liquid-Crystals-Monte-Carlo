import argparse
import math
import os

import numpy as np
from numpy import array as arr

from sklearn.decomposition import PCA

# custom imports
from LiquidCrystalSystem import LCSystem
from utilities import get_feature_func, get_nearest_neighbor_func


def create_feature_vectors_from_snapshot(coordinates, num_features, num_samples,
                                         feature_func=get_feature_func('relative_orientation'),
                                         nn_func=get_nearest_neighbor_func('euclidean_distance')):

    assert (num_features < len(coordinates)), \
        f"Number of features {num_features} cannot be greater than number of particles {len(coordinates)}"

    N = len(coordinates)
    feature_vectors = []

    # set the sampling rate for nearest neighbors
    if N % num_features == 0:
        nn_sampling_rate = N / num_features - 1
    else:
        nn_sampling_rate = math.floor(N / num_features)

    # set a random seed for reproducibility
    rng = np.random.default_rng(666)
    probe_indices = rng.choice(N, size=num_samples, replace=False)

    # select random probe particle for nearest neighbor distance
    for probe_index in probe_indices:
        # list of distances relative to probe particle
        chosen_coord = coordinates[probe_index]

        # sort based on nearest distance to probe
        nn_sorted = sorted(coordinates, key=nn_func)

        fv = []
        feature_particle_coordinates = [chosen_coord]
        # add feature based on nearest neighbor distance
        for i, c in enumerate(nn_sorted):

            if (i > 0) and (i % nn_sampling_rate) == 0:
                feature = feature_func(c, chosen_coord)
                feature_particle_coordinates.append(c)
                fv.append(feature)

            # stop adding features if total number is met
            if len(fv) == num_features:
                break

        feature_vectors.append(fv)

    return feature_vectors, feature_particle_coordinates


def create_data_matrix(systems):
    n_features = 10
    pts_per_snap = 5

    # number of rows AKA number of data-points
    # = number of densities * number of captures per density >= 1e6 * number of data-points per capture
    n_densities = len(systems)
    n_snaps = len([step for step in
                   systems[0.2913752913752914].snapshots
                   if step >= 1e6])

    print(f"Number of datapoints: {n_densities * n_snaps * pts_per_snap}")
    print(f"Number of features: {n_features}")

    # function for computing features
    feature_func = get_feature_func('orientation')

    data_matrix = []

    # iterate over densities
    for density in systems.keys():

        snapshot_at_step = systems[density].snapshots

        # iterate over Monte Carlo steps
        for mc_step in snapshot_at_step:

            # choose step where system has equilibrated
            if mc_step >= 1e6:
                snapshot = snapshot_at_step[mc_step]
                feature_vecs, _ = create_feature_vectors_from_snapshot(snapshot, n_features,
                                                                       pts_per_snap, feature_func)
                data_matrix = data_matrix + feature_vecs

    data_matrix = np.stack(data_matrix, axis=0)

    return data_matrix


def run_pca(pca_data):
    """
    run principal component analysis
    :return:
    """
    pca = PCA()
    pca.fit(pca_data)

    return pca.components_, pca.explained_variance_ratio_


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

        lc_systems_to_density[global_packing_fraction] = lc

    data = create_data_matrix(lc_systems_to_density)
    print(data.shape)
