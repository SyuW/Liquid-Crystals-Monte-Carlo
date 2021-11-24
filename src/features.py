import math
import os

import numpy as np

from datetime import date
from multiprocessing import Pool

# custom imports
from LiquidCrystalSystem import LCSystem
from utilities import get_feature_func, get_nearest_neighbor_func


def load_dataset(dataset_path, verbose=False):
    """

    :param verbose:
    :param dataset_path:
    :return:
    """

    # hard-coded values for dimensions, move to config file later
    R = 25
    r = 0
    b = 5
    a = 0.25
    # load a circular or annular confinement dataset
    if r == 0:
        confinement = "Circle"
    else:
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
    if verbose:
        print(f"Found particle numbers: {sorted([int(x) for x in systems.keys()])}")
    # return particle number: monte carlo step: positions array
    return systems


def create_data_matrix(systems, num_of_features, num_of_samples,
                       start=1000000, end=1500000, save_path=None, verbose=False):
    """
    Construct the data matrix for PCA input
    :param verbose:
    :param base_save_path:
    :param systems:
    :param num_of_features:
    :param num_of_samples:
    :param start:
    :param end:
    :return: data_matrix, samples
    """

    # function for computing features
    feature_func = get_feature_func('relative_orientation')
    data_matrix = []
    samples = dict()
    # for particle number
    for particle_number, system_at_n in systems.items():
        samples[particle_number] = []
        system_state_at_mc_step = systems[particle_number].snapshots
        # for monte carlo step, position array
        for mc_step, pos_array in system_at_n.snapshots.items():
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
    # stack feature vecs to form matrix
    data_matrix = np.stack(data_matrix, axis=0)
    # remove empty elements from samples
    samples = {n: sample for n, sample in samples.items() if sample != []}
    # print out useful info
    if verbose:
        print(f"Number of features used: {num_of_features}")
        print(f"Number of samples used: {num_of_samples}")
        print(f"Shape of data matrix: {data_matrix.shape}")
        for n in sorted([int(x) for x in samples.keys()]):
            print(f"N: {n}, total # of samples: {len(samples[n])}")
    # save all the arrays so that they can be loaded in later stages
    if save_path:
        save_path = os.path.join(save_path, "data")
        os.makedirs(save_path, exist_ok=True)
        # save the data matrix
        np.save(os.path.join(save_path, "data.npy"), data_matrix)
        # save the samples individually
        for particle_number, sample in samples.items():
            np.save(os.path.join(save_path, f"{particle_number}.npy"), np.array(sample))

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
    feature_particle_coords = []
    # select random probe particle for nearest neighbor distance
    for probe_index in probe_indices:
        # list of distances relative to probe particle
        probe_coord = coordinates[probe_index]
        feature_particle_coords = [probe_coord]
        # sort based on nearest distance to probe
        dist_to_P = lambda x: nn_func(x[:2], probe_coord[:2])
        nn_sorted = sorted(coordinates, key=dist_to_P)
        fv = []
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


if __name__ == "__main__":
    project_path = "C:\\Users\\Sam Yu\\Desktop\\School\\4A\\Phys_437A_Research_Project"
    data_path = os.path.join(project_path, "datasets\\r=0")
    s_path = os.path.join(project_path, f"results\\{date.today()}")
    lcs = load_dataset(dataset_path=data_path, verbose=True)

    # m, s = create_data_matrix(lcs, 10, 5, base_save_path=s_path, verbose=True)
