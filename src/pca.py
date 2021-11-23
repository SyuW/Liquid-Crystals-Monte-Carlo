import argparse
import os

import numpy as np
import matplotlib.pyplot as plt

from datetime import date
from sklearn.decomposition import PCA


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


def run_pca(pca_data, save_path, create_plots=True, verbose=False):
    """
    Run the PCA procedure on feature data collected from liquid crystal sims
    :param save_path:
    :param create_plots:
    :param pca_data:
    :param verbose:
    :return:
    """

    pca = PCA()
    pca.fit(pca_data)

    if verbose:
        print(f"First principal component: {pca.components_[0]}")
        print(f"Explained variance ratios: {pca.explained_variance_ratio_}")

    num_of_features = pca_data.shape[1]
    if create_plots:
        with plt.ioff():
            fig, ax = plt.subplots()
        # explained variances
        ax.plot(range(1, num_of_features + 1), pca.explained_variance_ratio_)
        ax.scatter(range(1, num_of_features + 1), pca.explained_variance_ratio_, marker="v")
        ax.grid()
        ax.set_xlabel("number of components")
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
        ax.set_xlabel("k-th component")
        ax.set_ylabel(r"$[\mathbf{w}_1]^k$")
        ax.set_title(r"Components of $w_1$")
        fig.savefig(os.path.join(save_path, "1st_pc_weights.png"))
        ax.cla()
        plt.close(fig)

    return pca.components_, pca.explained_variance_ratio_


def create_phase_diagram(phase_boundaries):
    """

    :param phase_boundaries:
    :return:
    """
    with plt.ioff():
        fig, ax = plt.subplots()
    return


def find_phase_transition(w1, samples, config, save_path, verbose=False):
    """

    :param config:
    :param save_path:
    :param verbose:
    :param w1:
    :param samples:
    :return:
    """
    R = config["R"]
    r = config["r"]
    b = config["b"]
    a = config["a"]
    nf = config["nf"]
    ns = config["ns"]

    means = []
    stds = []
    # scores_pca = pca.transform(data_matrix)[:, 0].reshape(len(samples.keys()), len(fvs))
    # compute scores for each density value
    for i, particle_number in enumerate(sorted(samples.keys())):
        fvs = samples[particle_number]
        if verbose:
            print(f"No. particles: {particle_number}, No. feature vectors: {len(fvs)}")
        # scores = scores_pca[i]
        scores = np.array([w1 @ vec for vec in fvs])
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
    densities = [N * a * b / (R ** 2 - r ** 2) for N in sorted(samples.keys())]
    critical_density = densities[np.argmax(norm_stds)]
    # bar for marking critical density
    particle_diff = 5
    spacing = particle_diff * a * b / (R ** 2 - r ** 2)

    # create plots
    with plt.ioff():
        fig, ax = plt.subplots()
    # plotting mean and standard deviation together
    ax.set_xlabel(r"Density $\eta$")
    ax.set_ylabel("Response")
    ax.grid()
    # plot the mean response
    ax.plot(densities, norm_means, color="orange", label="Mean")
    ax.scatter(densities, norm_means, color="orange")
    # plot the std response
    ax.plot(densities, norm_stds, color="blue", label="Std")
    ax.scatter(densities, norm_stds, color="blue")
    ax.set_title(f"PCA Summary for orientational features (s={ns}, f={nf})")
    # mark the critical density
    ax.axvline(x=critical_density, linestyle="--", color="black", label="Critical density")
    ax.axvspan(xmin=critical_density - spacing, xmax=critical_density + spacing, color="red", alpha=0.2)
    ax.legend()
    ax.cla()
    fig.savefig(os.path.join(save_path, "all.png"))
    plt.close(fig)

    print(fr"Critical density detected at {critical_density:.4f} by PCA")

    return critical_density


if __name__ == "__main__":
    from features import load_dataset, create_data_matrix

    sys_config = {"R": 25, "r": 0, "b": 5, "a": 0.25, "nf": 10, "ns": 5}

    project_path = "C:\\Users\\Sam Yu\\Desktop\\School\\4A\\Phys_437A_Research_Project"
    data_path = os.path.join(project_path, "datasets\\r=0")
    results_path = os.path.join(project_path, f"results\\{date.today()}",
                                f"features_{sys_config['nf']}_samples_{sys_config['ns']}")

    if os.path.exists(results_path):
        m, s = load_feature_data(results_path)
    else:
        lcs = load_dataset(dataset_path=data_path, verbose=False)
        m, s = create_data_matrix(lcs, 10, 5,
                                  base_save_path=os.path.join(project_path, f"results\\{date.today()}"),
                                  verbose=False)

    # principal components and explained variance ratios
    ws, evrs = run_pca(pca_data=m, save_path=results_path)
    # find phase transition
    eta_c = find_phase_transition(w1=ws[0], samples=s, )
