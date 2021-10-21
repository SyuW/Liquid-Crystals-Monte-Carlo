import argparse
import numpy as np
import os

from sklearn.decomposition import PCA

# custom imports
from interactions import compute_ellipse_line_intersection
from interactions import determine_ellipse_overlap
from LiquidCrystalSystem import LCSystem


def compute_local_packing_fraction(self, mc_step):
    # get system state at specified Monte Carlo step
    system_state = self.system_state_at_step[mc_step]
    a = self.sim_params['Semi Major Axis']
    b = self.sim_params['Semi Minor Axis']

    # generate rays
    theta_vals = np.arange(0, 2 * np.pi, 0.1)
    rays = [np.tan(theta) for theta in theta_vals]
    print(rays)


def pca_phase_transition():

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a liquid crystal system dataset")
    parser.add_argument("--data-path", help="path to dataset")

    args = parser.parse_args()

    # plot all snapshots for all system sizes in dataset
    for _path_ in os.listdir(args.data_path):
        full_path = os.path.join(args.data_path, _path_, "instanceRun")
        lc = LCSystem(lc_data_path=full_path)
