import argparse
import math
import os

import numpy as np
from numpy import array as arr

from sklearn.decomposition import PCA

# custom imports
from LiquidCrystalSystem import LCSystem


def create_feature_vectors_from_snapshot(coordinates, num_features):
    assert (num_features < len(coordinates)), \
        f"Number of features {num_features} cannot be greater than number of particles {len(coordinates)}"

    N = len(coordinates)
    nn_sampling_number = math.floor(N / num_features)
    print(f"Nearest neighbor sampling number: {nn_sampling_number}")

    # set a random seed for reproducibility
    rng = np.random.default_rng(666)
    probe_indices = rng.choice(N, size=num_features, replace=False)

    # x, y positions for particle center of masses
    spatial = [c[:2] for c in coordinates]
    # angle of long axis with respect to x for particles
    angular = [(c[-1] % np.pi) for c in coordinates]

    feature_vectors = []

    for probe_index in probe_indices:
        print(f"Probe index: {probe_index}")

        # list of distances relative to probe particle
        chosen_coord = coordinates[probe_index]

        # special functions for nearest neighbor sort and feature calculation
        norm_squared = lambda x, y: (arr(x) - arr(y)) @ (arr(x) - arr(y))
        distance_to_probe = lambda x: norm_squared(x[:2], chosen_coord[:2])

        # may abstract as function param later
        feature_func = lambda x, y: abs(np.cos(x - y))

        nn_sorted = sorted(coordinates, key=distance_to_probe)

        fv = []
        # add feature based on nearest neighbor distance
        for i, c in enumerate(nn_sorted):
            if (i > 0) and (i % nn_sampling_number) == 0:
                feature = feature_func(c[-1], chosen_coord[-1])
                fv.append(feature)

            if len(fv) == num_features:
                break

        feature_vectors.append(fv)

    return feature_vectors


def pca_phase_transition(lc_systems_to_density):
    for density in lc_systems_to_density:
        print(density)

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a liquid crystal system dataset")
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

    print(lc_systems_to_density.keys())
