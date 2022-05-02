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
from symmetry import symmetry_group_detection
from visualizations import plot_single_state, plot_single_state_centers, plot_ellipse_ends_only


def load_feature_data(path):
    """
    Load feature data
    :param path:
    :return:
    """
    samples = dict()
    for _path_ in os.listdir(os.path.join(path, "data")):
        basename = os.path.basename(_path_)
        if _path_ == "data.npy":
            matrix = np.load(os.path.join(path, "data", _path_))
        else:
            particle_number = int(basename.strip(".npy"))
            samples[particle_number] = np.load(os.path.join(path, "data", _path_))
    return matrix, samples


def pca_procedure():
    return


def run_pca(pca_data, save_path, create_plots=True, verbose=False):
    """
    Run the PCA procedure on feature data collected from liquid crystal sims
    :param save_path:
    :param create_plots:
    :param pca_data:
    :param verbose:
    :return:
    """

    # pca_whiten = PCA()
    # pca_whiten.fit(pca_data)
    # whitened_data = pca_whiten.transform(pca_data)
    pca = PCA()
    pca.fit(pca_data)

    if verbose:
        print(f"First principal component: {pca.components_[0]}")
        print(f"Explained variance ratios: {pca.explained_variance_ratio_}")

    if not os.path.exists(save_path):
        os.makedirs(save_path, exist_ok=True)

    num_of_features = pca_data.shape[1]
    if create_plots:
        with plt.ioff():
            fig, ax = plt.subplots()
        # explained variances
        ax.plot(range(1, num_of_features + 1), pca.explained_variance_ratio_)
        ax.scatter(range(1, num_of_features + 1), pca.explained_variance_ratio_, marker="v")
        ax.grid()
        ax.set_xlabel("k-th principal component")
        ax.set_ylabel("Explained variance")
        ax.set_title("Explained variance ratios from PCA")
        fig.savefig(os.path.join(save_path, "explained_variances.png"))
        ax.cla()
        # Cumulative explained variance ratios
        ax.plot(range(1, num_of_features + 1), np.cumsum(pca.explained_variance_ratio_))
        ax.scatter(range(1, num_of_features + 1), np.cumsum(pca.explained_variance_ratio_), marker="v")
        ax.grid()
        ax.set_xlabel("number of components")
        ax.set_ylabel("cumulative explained variance")
        ax.set_title("Cumulative explained variance ratios from PCA")
        fig.savefig(os.path.join(save_path, "cumulative_explained_variances.png"))
        ax.cla()
        # weights of first principal component
        w1 = pca.components_[0]
        w1 = np.abs(w1)
        ax.plot(range(1, num_of_features + 1), w1)
        ax.scatter(range(1, num_of_features + 1), w1, marker="s")
        ax.grid()
        ax.set_xlabel(r"component weight, $[\vec{w_1}]_k$")
        ax.set_ylabel(r"$[\mathbf{w}_1]^k$")
        ax.set_title(r"Weights of $w_1$")
        fig.savefig(os.path.join(save_path, "1st_pc_weights.png"))
        ax.cla()
        plt.close(fig)

    return pca.components_, pca.explained_variance_ratio_


def create_phase_diagram(phase_boundaries):
    """

    :param phase_boundaries:
    :return:
    """

    return


def find_phase_transition(w1, data_matrix, config, samples, save_path, verbose=False):
    """

    :param data_matrix:
    :param config:
    :param save_path:
    :param verbose:
    :param w1:
    :param samples:
    :return:
    """

    num_samples = config["ns"]
    num_features = config["nf"]
    means = []
    stds = []
    # scores_pca = pca.transform(data_matrix)[:, 0].reshape(len(samples.keys()), len(fvs))
    # compute scores for each area fraction value
    for area_fraction in sorted(samples):
        sample = samples[area_fraction]
        if verbose:
            print(f"Area Fraction: {area_fraction}, No. feature vectors: {len(sample)}")
        # scores = scores_pca[i]
        scores = np.array([w1 @ vec for vec in sample])
        # print(scores)
        avg = np.mean(scores)
        std = np.std(scores)
        means.append(avg)
        stds.append(std)

    # convert to numpy
    means = np.array(means)
    stds = np.array(stds)
    # normalize the values
    norm_means = (means - min(means)) / max(means - min(means))
    norm_stds = stds / max(stds)
    densities = sorted(samples.keys())
    # find the critical density for phase transition
    critical_density = densities[np.argmax(norm_stds)]
    # bar for marking critical density
    particle_diff = 5

    # create plots
    with plt.ioff():
        fig, ax = plt.subplots()

    # plot the row means
    row_means = [np.mean(row) for row in data_matrix][::50]
    ax.plot(range(len(row_means)), row_means)
    ax.set_xlabel(r"Row of $X$")
    ax.set_ylabel("Mean of row")
    # ax.set_ylim(0, 1)
    ax.grid()
    ax.set_title("Plot of row means of data matrix")
    fig.savefig(os.path.join(save_path, "row_means.png"))
    ax.cla()
    # plot the means of samples
    fv_means = [np.mean(np.mean(samples[N])) for N in sorted(samples.keys())]
    ax.plot(densities, fv_means)
    ax.set_xlabel(r"Area fraction, $\phi$")
    ax.set_ylabel("Mean of samples for N")
    # ax.set_ylim(0, 1)
    ax.grid()
    ax.set_title("Plot of mean of samples for particle numbers")
    fig.savefig(os.path.join(save_path, "samples_means.png"))
    ax.cla()
    # plotting normalized mean and normalized standard deviation together
    ax.set_xlabel(r"Area fraction $\phi$")
    ax.set_ylabel(r"Order parameter ${P_1}$")
    ax.grid()
    # plot the mean response
    ax.plot(densities, means, color="orange")
    ax.scatter(densities, means, color="orange", marker="<", label=r"$p_1$")
    # plot the std response
    ax.plot(densities, stds, color="blue")
    ax.scatter(densities, stds, color="blue", marker=">", label=r"$\sigma_1$")
    ax.set_title(f"Means and Stds plot (s={num_samples}, f={num_features})")
    ax.legend()
    fig.savefig(os.path.join(save_path, "means_stds.png"))
    ax.cla()
    # plotting normalized mean and normalized standard deviation together
    ax.set_xlabel(r"Area fraction $\eta$")
    ax.set_ylabel(r"Normalized order parameter $\bar{P}_1$")
    ax.grid()
    # plot the mean response
    ax.plot(densities, norm_means, color="orange")
    ax.scatter(densities, norm_means, color="orange", marker=">", label=r"$\tilde{p}_1$")
    # plot the std response
    ax.plot(densities, norm_stds, color="blue")
    ax.scatter(densities, norm_stds, color="blue", marker="<", label=r"$\tilde{\sigma}_1$")
    # ax.set_title(f"Normalized Means and Stds plot (s={ns}, f={nf})")
    # mark the critical density
    ax.axvline(x=critical_density, linestyle="--", color="black", label="Critical density")
    # ax.axvspan(xmin=critical_density - spacing, xmax=critical_density + spacing, color="red", alpha=0.2)
    ax.legend()
    fig.savefig(os.path.join(save_path, "norm_means_stds.png"))
    ax.cla()
    plt.close(fig)

    print(fr"Critical density detected at {critical_density:.4f} by PCA")

    return critical_density


if __name__ == "__main__":
    from features import load_dataset, create_data_matrix, create_feature_vectors_from_snapshot
    from visualizations import create_feature_vector_visualization

    # start the timer
    start = time.perf_counter()

    sys_config = {"R": 30.1, "r": 0, "b": 5, "a": 0.25, "nf": 150, "ns": 15, "start": 200000, "end": 800000}

    project_path = "C:\\Users\\Sam Yu\\Desktop\\School\\4B\\Phys_437B_RP"
    data_path = os.path.join(project_path, f"datasets\\r={sys_config['r']}\\3_6_2022")
    current_time = strftime("%Y-%m-%d_%H-%M-%S", gmtime())
    results_path = f"results\\r={sys_config['r']}\\{current_time}_nf_{sys_config['nf']}_ns_{sys_config['ns']}" \
                   f"_start_{sys_config['start']}_end_{sys_config['end']}"
    results_path = os.path.join(project_path, results_path)
    cache_path = os.path.join(project_path, f"cache\\r={sys_config['r']}")

    # make the save path
    os.makedirs(results_path, exist_ok=True)

    # cache the file so that we don't to load the dataset again each time
    load_dataset_from_cache = True
    dataset_cache = os.path.join(cache_path, "lcs.pickle")
    # load dataset from cache option
    if load_dataset_from_cache:
        # check the cached dataset file exists
        if os.path.isfile(dataset_cache):
            print("Cached dataset found, loading...")
            with open(dataset_cache, "rb") as in_f:
                lcs = pickle.load(in_f)
        # else, make a new one
        else:
            print("Cached dataset file not found, making a new one...")
            lcs = load_dataset(dataset_path=data_path)
            with open(dataset_cache, "wb") as out_f:
                pickle.dump(lcs, out_f)
    # make a new one if caching option was not selected
    else:
        print("Cached dataset file not found, making a new one...")
        lcs = load_dataset(dataset_path=data_path)
        with open(dataset_cache, "wb") as out_f:
            pickle.dump(lcs, out_f)

    print(f"Dataset Loading Execution Time: {time.perf_counter() - start} seconds")

    # analyze using nematic order parameter
    lambda_means = []
    lambda_vars = []
    area_fractions = sorted(lcs.keys())
    for a_f in area_fractions:
        lambdas = []
        for mc_step in lcs[a_f].states.keys():
            lambdas_at_step = nematic_order_param(lcs[a_f].states[mc_step])
            lambdas.append(lambdas_at_step)
        lambda_means.append(np.mean(lambdas))
        lambda_vars.append(np.var(lambdas))
    with plt.ioff():
        plt.scatter(area_fractions, lambda_means, color="blue", marker="+", label=r"mean($\Lambda$)")
        plt.plot(area_fractions, lambda_means, color="blue", linestyle="--")
        plt.scatter(area_fractions, lambda_vars, color="red", marker="*", label=r"var($\Lambda$)")
        plt.plot(area_fractions, lambda_vars, color="red", linestyle="--")
        plt.scatter(area_fractions, np.array(lambda_vars) / np.max(lambda_vars), color="orange", marker="*",
                    label=r"var($\Lambda$)/max(var($\Lambda$))")
        plt.plot(area_fractions, np.array(lambda_vars) / np.max(lambda_vars), color="orange", linestyle="--")
        # plt.plot(area_fractions, np.array(lambda_means) * np.array(lambda_vars), color="orange", linestyle="--",
        #          label=r"mean($\Lambda$)*var($\Lambda$)")
        plt.xlabel(r"Area fractions, $\phi$")
        plt.ylabel(r"Order parameter, $\Lambda$")
        plt.title(fr"Order Parameter $\Lambda$ analysis (steps {sys_config['start']}-{sys_config['end']})")
        plt.legend()
        plt.savefig(os.path.join(results_path, "lambda_nematic_order_param.png"))
        plt.close()

    print(f"Order parameter Lambda calculation time: {time.perf_counter() - start} seconds")

    # fold of symmetry experiments
    symmetry_analysis_path = os.path.join(results_path, "symm")
    os.makedirs(symmetry_analysis_path, exist_ok=True)
    chosen_densities = [0.105, 0.27, 0.5816669999999999]
    chosen_mc_step = 800000
    analyze_symm = True
    if analyze_symm:
        for chosen_density in chosen_densities:
            lc_system = lcs[chosen_density]
            snapshot = lc_system.states[chosen_mc_step]
            snapshot_fig = plot_single_state(pos_array=snapshot, mc_step=chosen_mc_step, N=lc_system.num_of_particles,
                                             R=lc_system.outer_radius, r=lc_system.inner_radius,
                                             a=lc_system.minor_axis, b=lc_system.major_axis)
            snapshot_fig.savefig(os.path.join(symmetry_analysis_path, f"chosen_snapshot"
                                                                      f"_phi_{chosen_density}_step_{chosen_mc_step}.png"))
            plt.close(snapshot_fig)
            center_of_masses_fig = plot_single_state_centers(pos_array=snapshot, mc_step=chosen_mc_step,
                                                             N=lc_system.num_of_particles, R=lc_system.outer_radius,
                                                             r=lc_system.inner_radius, a=lc_system.minor_axis,
                                                             b=lc_system.major_axis)
            center_of_masses_fig.savefig(os.path.join(symmetry_analysis_path,
                                                      f"center_of_masses_phi_{chosen_density}_step_{chosen_mc_step}.png"))
            plt.close(center_of_masses_fig)
            shifted_coords_1 = np.zeros_like(snapshot[:, :2])
            shifted_coords_1[:, 1] = snapshot[:, 1] + lc_system.major_axis * np.sin(snapshot[:, -1])
            shifted_coords_1[:, 0] = snapshot[:, 0] + lc_system.major_axis * np.cos(snapshot[:, -1])
            shifted_coords_2 = np.zeros_like(snapshot[:, :2])
            shifted_coords_2[:, 1] = snapshot[:, 1] - lc_system.major_axis * np.sin(snapshot[:, -1])
            shifted_coords_2[:, 0] = snapshot[:, 0] - lc_system.major_axis * np.cos(snapshot[:, -1])
            ellipse_ends_fig = plot_ellipse_ends_only(shifted_coords_1, shifted_coords_2, com_pos_array=snapshot,
                                                      mc_step=chosen_mc_step, N=lc_system.num_of_particles,
                                                      R=lc_system.outer_radius, r=lc_system.inner_radius,
                                                      a=lc_system.minor_axis, b=lc_system.major_axis)
            ellipse_ends_fig.savefig(os.path.join(symmetry_analysis_path, f"ellipse_ends_phi_{chosen_density}"
                                                                          f"_step_{chosen_mc_step}.png"))
            plt.close(ellipse_ends_fig)
            radial_bands_1 = symmetry_group_detection(shifted_coords_1, lc_system.outer_radius)
            radial_bands_2 = symmetry_group_detection(shifted_coords_2, lc_system.outer_radius)
            thetas1 = sorted(radial_bands_1.keys())
            thetas2 = sorted(radial_bands_2.keys())
            r_dist_1 = [lc_system.outer_radius - max(radial_bands_1[t]) for t in thetas1]
            r_dist_2 = [lc_system.outer_radius - max(radial_bands_2[t]) for t in thetas1]
            with plt.ioff():
                plt.plot(thetas1, r_dist_1, color="red", label="red ellipse ends")
                plt.plot(thetas2, r_dist_2, color="blue", label="blue ellipse ends")
                plt.xlabel(r"Theta, $\theta$")
                plt.ylabel(r"$R - $max($\{\sqrt{x^2+y^2}|(x,y)\in{S_{\theta}}\}$)")
                plt.title(fr"for $\phi={chosen_density},step={chosen_mc_step}$")
                plt.legend()
                plt.savefig(os.path.join(symmetry_analysis_path, f"radial_distributions_phi_{chosen_density}"
                                                                 f"_step_{chosen_mc_step}_d_{5}.png"))
                plt.close()
                plt.plot(thetas1, np.abs(np.array(r_dist_1) - np.array(r_dist_2)), color="red")
                plt.xlabel(r"Theta, $\theta$")
                plt.ylabel(r"Absolute Difference")
                plt.title(fr"Absolute difference between red vs blue radial distributions")
                plt.savefig(os.path.join(symmetry_analysis_path, f"distribution_difference_phi_{chosen_density}"
                                                                 f"_step_{chosen_mc_step}.png"))
                plt.close()

    print(f"Time for rotational symmetry method: {time.perf_counter() - start} seconds.")

    # functions for computing features
    ff = get_feature_func('relative_orientation')
    nf = get_nearest_neighbor_func('euclidean_distance')

    # cache the constructed data matrix so we don't have to create it from scratch
    feature_samples_path = f"features_R_{sys_config['R']}_r_{sys_config['r']}_nf_{sys_config['nf']}" \
                           f"_ns_{sys_config['ns']}_start_{sys_config['start']}_end_{sys_config['end']}.pickle"
    features_cache = os.path.join(cache_path, feature_samples_path)
    load_features_from_cache = True
    # load features from cache option
    if load_features_from_cache:
        if os.path.isfile(features_cache):
            print("Cached features found, loading...")
            with open(features_cache, "rb") as in_f:
                sam = pickle.load(in_f)
                mat = []
                for af in sorted(sam.keys()):
                    mat += sam[af]
                mat = np.stack(mat, axis=0)
        else:
            print("Cached features not found, creating from scratch...")
            mat, sam = create_data_matrix(lcs,
                                          num_of_features=sys_config["nf"], num_of_samples=sys_config["ns"],
                                          feature_func=ff, neighbor_func=nf,
                                          start=200000, end=800000, verbose=True)
            with open(features_cache, "wb") as out_f:
                pickle.dump(sam, out_f)
    # load features from cache option not selected -> create from scratch
    else:
        print("Feature caching option not selected, creating from scratch...")
        mat, sam = create_data_matrix(lcs,
                                      num_of_features=sys_config["nf"], num_of_samples=sys_config["ns"],
                                      feature_func=ff, neighbor_func=nf,
                                      start=200000, end=800000, verbose=True)
        with open(features_cache, "wb") as out_f:
            pickle.dump(sam, out_f)

    print(f"Data Matrix Creation Time: {time.perf_counter() - start} seconds")

    # visualize examples of the feature vector selection
    chosen_densities = [0.15083300000000002]
    chosen_mc_steps = [800000]
    fv_vis_path = os.path.join(results_path, "pca", "vis")
    os.makedirs(fv_vis_path, exist_ok=True)
    # option for plotting feature vector visualizations
    plot_fv_visualizations = False
    if plot_fv_visualizations:
        for chosen_density in chosen_densities:
            lc_system = lcs[chosen_density]
            chosen_mc_steps = list(lc_system.states.keys())
            chosen_mc_steps = [800000]
            for mc_step in chosen_mc_steps:
                snapshot = lc_system.states[mc_step]
                fvs, fvcs = create_feature_vectors_from_snapshot(snapshot, num_features=sys_config["nf"], num_samples=1)
                probe_coord = list(fvcs.keys())[0]
                feature_particles = fvcs[probe_coord]
                fv_vis_fig = create_feature_vector_visualization(snapshot, probe=probe_coord,
                                                                 feature_particles=feature_particles,
                                                                 mc_step=mc_step,
                                                                 inner_radius=lc_system.inner_radius,
                                                                 outer_radius=lc_system.outer_radius,
                                                                 a=lc_system.minor_axis,
                                                                 b=lc_system.major_axis)
                fv_vis_fig.savefig(os.path.join(fv_vis_path, f"phi_{chosen_density}_nf_{sys_config['nf']}"
                                                             f"_ns_{sys_config['ns']}_step_{mc_step}.png"))
                plt.close()

    # principal components and explained variance ratios
    ws, evrs = run_pca(pca_data=mat, save_path=os.path.join(results_path, "pca"))

    print(f"PCA Analysis Time: {time.perf_counter() - start} seconds")
    # find phase transition
    # there's an ambiguity in the direction that the principal components point in, up to a global sign
    eta_c = find_phase_transition(w1=-ws[0], data_matrix=mat, samples=sam, save_path=os.path.join(results_path, "pca"),
                                  config=sys_config, verbose=True)
