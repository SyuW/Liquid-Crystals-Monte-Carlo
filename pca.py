import argparse
import os

import numpy as np
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA

# custom imports
from LiquidCrystalSystem import LCSystem
from features import create_data_matrix


def run_pca(pca_data, save_path, create_plots=True, save_arrays=True, verbose=False):
    """

    :param save_arrays:
    :param create_plots:
    :param pca_data:
    :param save_path:
    :param verbose:
    :return:
    """

    pca = PCA()
    pca.fit(pca_data)

    if verbose:
        print(f"First principal component: {pca.components_[0]}")
        print(f"Explained variance ratios: {pca.explained_variance_ratio_}")

    if save_arrays:
        np.save(pca.components_[0])
        np.save(pca.explained_variance_ratio_[0])

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
