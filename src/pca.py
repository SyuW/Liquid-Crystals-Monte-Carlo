import argparse
import os

import numpy as np
import matplotlib.pyplot as plt

from datetime import date
from sklearn.decomposition import PCA
from time import gmtime, strftime


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

    return


def find_phase_transition(w1, data_matrix, samples, config, save_path, verbose=False):
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
    for particle_number, sample in samples.items():
        fvs = samples[particle_number]
        if verbose:
            print(f"No. particles: {particle_number}, No. feature vectors: {len(sample)}")
        # scores = scores_pca[i]
        scores = np.array([w1 @ vec for vec in fvs])
        # print(scores)
        avg = np.mean(scores)
        std = np.std(scores)
        means.append(avg)
        stds.append(std)

    pca = PCA(n_components=1)
    pca.fit(data_matrix)
    pca_fit = pca.transform(data_matrix)
    inverse = pca.inverse_transform(pca_fit)

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

    # plot the row means
    row_means = [np.mean(row) for row in m][::50]
    ax.plot(range(len(row_means)), row_means)
    ax.set_xlabel(r"Row of $X$")
    ax.set_ylabel("Mean of row")
    ax.set_ylim(0, 1)
    ax.set_title("Plot of row means of data matrix")
    fig.savefig(os.path.join(save_path, "row_means.png"))
    ax.cla()
    # plot the means of samples
    fv_means = [np.mean(np.mean(samples[N])) for N in samples]
    ax.plot(sorted(samples.keys()), fv_means)
    ax.set_xlabel("Particle Number, N")
    ax.set_ylabel("Mean of samples for N")
    ax.set_ylim(0, 1)
    ax.grid()
    ax.set_title("Plot of mean of samples for particle numbers")
    fig.savefig(os.path.join(save_path, "samples_means.png"))
    ax.cla()
    # plotting normalized mean and normalized standard deviation together
    ax.set_xlabel(r"Density $\eta$")
    ax.set_ylabel(r"Order parameter ${P_1}$")
    ax.grid()
    # plot the mean response
    ax.plot(densities, means, color="orange", label="Mean")
    ax.scatter(densities, means, color="orange")
    # plot the std response
    ax.plot(densities, stds, color="blue", label="Std")
    ax.scatter(densities, stds, color="blue")
    ax.set_title(f"Means and Stds plot (s={ns}, f={nf})")
    ax.legend()
    fig.savefig(os.path.join(save_path, "means_stds.png"))
    ax.cla()
    # plotting normalized mean and normalized standard deviation together
    ax.set_xlabel(r"Density $\eta$")
    ax.set_ylabel(r"Normalized order parameter $\bar{P}_1$")
    ax.grid()
    # plot the mean response
    ax.plot(densities, norm_means, color="orange", label="Mean")
    ax.scatter(densities, norm_means, color="orange")
    # plot the std response
    ax.plot(densities, norm_stds, color="blue", label="Std")
    ax.scatter(densities, norm_stds, color="blue")
    ax.set_title(f"Normalized Means and Stds plot (s={ns}, f={nf})")
    # mark the critical density
    ax.axvline(x=critical_density, linestyle="--", color="black", label="Critical density")
    ax.axvspan(xmin=critical_density - spacing, xmax=critical_density + spacing, color="red", alpha=0.2)
    ax.legend()
    fig.savefig(os.path.join(save_path, "norm_means_stds.png"))
    ax.cla()
    plt.close(fig)

    print(fr"Critical density detected at {critical_density:.4f} by PCA")

    return critical_density


if __name__ == "__main__":
    from features import load_dataset, create_data_matrix

    sys_config = {"R": 25, "r": 0, "b": 5, "a": 0.25, "nf": 15, "ns": 8}

    project_path = "C:\\Users\\Sam Yu\\Desktop\\School\\4A\\Phys_437A_Research_Project"
    data_path = os.path.join(project_path, "datasets\\r=0")
    current_time = strftime("%Y-%m-%d_%H-%M-%S", gmtime())
    results_path = os.path.join(project_path, f"results\\r={sys_config['r']}\\{current_time}")

    lcs = load_dataset(dataset_path=data_path, verbose=False)
    m, s = create_data_matrix(lcs, sys_config["nf"], sys_config["ns"],
                              verbose=True)

    # principal components and explained variance ratios
    ws, evrs = run_pca(pca_data=m, save_path=results_path)
    # find phase transition
    eta_c = find_phase_transition(w1=np.abs(ws[0]), data_matrix=m, samples=s,
                                  config=sys_config, save_path=results_path)
