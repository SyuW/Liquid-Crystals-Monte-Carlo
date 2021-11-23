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


def create_phase_diagram(phase_boundaries, ):
    return


def find_phase_transition(w1, samples, ):
    return


if __name__ == "__main__":
    from features import load_dataset, create_data_matrix

    nf = 10
    ns = 5

    project_path = "C:\\Users\\Sam Yu\\Desktop\\School\\4A\\Phys_437A_Research_Project"
    data_path = os.path.join(project_path, "datasets\\r=0")
    results_path = os.path.join(project_path, f"results\\{date.today()}",
                                f"features_{nf}_samples_{ns}")

    if os.path.exists(results_path):
        m, s = load_feature_data(results_path)
    else:
        lcs = load_dataset(dataset_path=data_path, verbose=False)
        m, s = create_data_matrix(lcs, 10, 5,
                                  base_save_path=os.path.join(project_path, f"results\\{date.today()}"),
                                  verbose=False)

    # principal components and explained variance ratios
    ws, evrs = run_pca(pca_data=m, save_path=results_path)


