import pickle
from argparse import Namespace
from datetime import date
from joblib import Parallel, delayed
import math
from multiprocessing import Pool
import numpy as np
import os

# custom imports
from LiquidCrystalSystem import LCSystem
from utilities import get_feature_func, get_nearest_neighbor_func
from visualizations import create_feature_vector_visualization


def retrieve_system_data(run_path):
    """
    Get system data from files at path
    :param run_path:
    :return:
    """

    with open(os.path.join(run_path, "checkpoint.pickle"), "rb") as in_f:
        params = pickle.load(in_f)

    _system_ = Namespace()
    # get the main parameters for the simulation run
    _system_.outer_radius = params["R"]
    _system_.inner_radius = params["r"]
    _system_.num_of_particles = len(params["pos_array"])
    _system_.total_steps = params["end_step"]
    _system_.minor_axis = params["a"]
    _system_.major_axis = params["b"]
    _system_.area_fraction = len(params["pos_array"]) * params["a"] * params["b"] / \
                             (params["R"] ** 2 - params["r"] ** 2)
    # get paths of all data files and retrieve system state information for each Monte Carlo step
    pos_array_paths = [os.path.join(run_path, p) for p in os.listdir(run_path) if p.endswith(".csv")]
    pos_array_paths = [os.path.normpath(p) for p in pos_array_paths]
    states = dict()
    for p in pos_array_paths:
        basename = os.path.basename(p)
        mc_step = int(basename.strip("step_").strip(".csv"))
        pos_array = np.loadtxt(os.path.join(run_path, p), delimiter=",", dtype=np.float32)
        states[mc_step] = pos_array
    _system_.states = states

    return _system_


def load_dataset(dataset_path):
    """
    Load dataset from path
    :param dataset_path: path containing all the run folders produced from simulations
    :return:
    """

    systems = Parallel(n_jobs=4)(delayed(retrieve_system_data)(os.path.join(dataset_path, _p_))
                                 for _p_ in os.listdir(dataset_path))
    systems = {s.area_fraction: s for s in systems}
    return systems


def create_data_matrix(systems, num_of_features, num_of_samples, feature_func, neighbor_func,
                       start=1000000, end=1500000, save_path=None, verbose=False):
    """
    Construct the data matrix for PCA input
    :param neighbor_func:
    :param feature_func:
    :param save_path:
    :param verbose:
    :param systems:
    :param num_of_features:
    :param num_of_samples:
    :param start:
    :param end:
    :return: data_matrix, samples
    """

    data_matrix = []
    samples = {af: [] for af in systems.keys()}
    # for particle number
    for area_fraction in sorted(samples.keys()):
        system_at_area_fraction = systems[area_fraction]
        # for monte carlo step, position array
        for mc_step, pos_array in system_at_area_fraction.states.items():
            # if monte carlo step of snapshot is between start and end
            if start <= mc_step <= end:
                # get snapshot of system at Monte Carlo step
                snapshot = system_at_area_fraction.states[mc_step]
                # create the feature vectors
                fvs, _ = create_feature_vectors_from_snapshot(snapshot,
                                                              num_features=num_of_features,
                                                              num_samples=num_of_samples,
                                                              feature_func=feature_func,
                                                              nn_func=neighbor_func)

                # add to data matrix
                data_matrix = data_matrix + fvs

                # add to samples
                samples[area_fraction] += fvs

    # stack feature vecs to form matrix
    data_matrix = np.stack(data_matrix, axis=0)
    # remove empty elements from samples
    samples = {n: sample for n, sample in samples.items() if sample != []}
    # print out useful info
    if verbose:
        print(f"Number of features used per feature vector: {num_of_features}")
        print(f"Number of samples used per area fraction: {num_of_samples}")
        print(f"Shape of data matrix: {data_matrix.shape}")
        # for n in sorted([int(x) for x in samples.keys()]):
        #    print(f"N: {n}, # of samples: {len(samples[n])}")
        # print(f"Total # of samples: {sum([len(samples[int(n)]) for n in samples.keys()])}")

    return data_matrix, samples


def create_feature_vectors_from_snapshot(coordinates, num_features, num_samples,
                                         feature_func=get_feature_func('relative_orientation'),
                                         nn_func=get_nearest_neighbor_func('euclidean_distance')):
    """
    Create the feature vectors from provided particle coordinates

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

    # choosing probe particle using normal distribution
    gaussian_points = np.random.multivariate_normal((0, 0), np.identity(2)*4, num_samples)
    probe_indices = []
    for p in gaussian_points:
        dists = [nn_func(c[:2], p) for c in coordinates]
        probe_indices.append(np.argmin(dists))

    # initialize stuff
    feature_vectors = []
    feature_particle_coords = dict()
    # select random probe particle for nearest neighbor distance
    for probe_index in probe_indices:
        # list of distances relative to probe particle
        probe_coord = coordinates[probe_index]
        feature_particle_coords[tuple(probe_coord)] = []
        # sort based on nearest distance to probe
        dist_to_P = lambda x: nn_func(x[:2], probe_coord[:2])
        nn_sorted = sorted(coordinates, key=dist_to_P)
        fv = []
        # add feature based on nearest neighbor distance
        for i, c in enumerate(nn_sorted):
            # don't add the probe particle
            if (i > 0) and (i % nn_sampling_rate) == 0:
                feature = feature_func(c, probe_coord)
                feature_particle_coords[tuple(probe_coord)].append(c)
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
