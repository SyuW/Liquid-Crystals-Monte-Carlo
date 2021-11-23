import argparse
import os

import numpy as np
import matplotlib.pyplot as plt

from datetime import date
from sklearn.decomposition import PCA


def run_pca(pca_data, base_save_path, create_plots=True,
            load_path=None, save_path=None, verbose=False):
    """

    :param save_path:
    :param load_path:
    :param base_save_path:
    :param create_plots:
    :param pca_data:
    :param verbose:
    :return:
    """
    # create the save path for plots
    save_path = os.path.join(base_save_path, date.today())

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
        ax.set_xlabel("k")
        ax.set_ylabel(r"$[\mathbf{w}_1]^k$")
        ax.set_title(r"Components of $w_1$")
        fig.savefig(os.path.join(save_path, "1st_pc_weights.png"))
        ax.cla()
        fig.close()

    return pca.components_, pca.explained_variance_ratio_


def visualize_feature_selection(systems):
    return


def plot_all_snapshots(systems):
    return


if __name__ == "__main__":
    from features import load_dataset, create_data_matrix

    project_path = "C:\\Users\\Sam Yu\\Desktop\\School\\4A\\Phys_437A_Research_Project"
    data_path = os.path.join(project_path, "datasets\\r=0")
    s_path = os.path.join(project_path, f"results\\{date.today()}")
    lcs = load_dataset(dataset_path=data_path, verbose=True)
    m, s = create_data_matrix(lcs, 10, 5, base_save_path=s_path, verbose=True)
