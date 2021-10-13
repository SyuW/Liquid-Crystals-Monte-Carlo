import argparse
import numpy as np
import os

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

plt.ioff()


class LCDataset:

    def plot_snapshot(self, mc_step):
        assert (mc_step in self.system_state_at_step.keys()), f"Monte Carlo step {mc_step} is not valid"
        # create figure and axes
        fig, ax = plt.subplots()
        fig.set_size_inches(10, 10)

        # add annular boundaries to plot
        outer_circle = plt.Circle((0, 0), self.sim_params["R"], color='black', fill=False, linewidth=2.2)
        inner_circle = plt.Circle((0, 0), self.sim_params["r"], color='black', fill=False, linewidth=2.2)
        ax.add_patch(inner_circle)
        ax.add_patch(outer_circle)

        # plot liquid crystal ellipses
        a = self.sim_params['Semi Major Axis']
        b = self.sim_params['Semi Minor Axis']
        for crystal_pos in self.system_state_at_step[mc_step]:
            crystal = Ellipse(xy=crystal_pos[:-1], angle=(180 / np.pi) * crystal_pos[-1], width=2 * a, height=2 * b,
                              linewidth=1.7, color='black', fill=False)
            ax.add_patch(crystal)

        num_ellipses = self.sim_params["# of Ellipse"]
        outer_radius = self.sim_params["R"]
        inner_radius = self.sim_params["r"]

        circle_pad = 5
        ax.set_xlim(-outer_radius - circle_pad, outer_radius + circle_pad)
        ax.set_ylim(-outer_radius - circle_pad, outer_radius + circle_pad)

        ax.set_title(f'N={num_ellipses} b={b} k={a / b} R={outer_radius} r={inner_radius}', size=20)
        ax.tick_params(axis='both', labelsize=20)

        return fig

    def compute_local_packing_fraction(self, mc_step):
        # get system state at specified Monte Carlo step
        system_state = self.system_state_at_step[mc_step]
        a = self.sim_params['Semi Major Axis']
        b = self.sim_params['Semi Minor Axis']

    def __init__(self, lc_data_path):
        # get parameters from simulation
        self.sim_params = dict()
        with open(os.path.join(lc_data_path, "MonteCarlo_Annulus_SimNotes.txt")) as file:
            for line in file.readlines()[2:]:
                info = line.strip().split(": ")
                # special case of acceptance rate
                if info[0] == "Acceptance Rate":
                    self.sim_params[info[0]] = float(info[1].strip(" %")) / 100
                else:
                    self.sim_params[info[0]] = float(info[1])

        # get paths of all data files and retrieve system state information for each Monte Carlo step
        state_file_paths = [os.path.join(lc_data_path, p) for p in os.listdir(lc_data_path) if p.endswith(".csv")]
        state_file_paths = [os.path.normpath(p) for p in state_file_paths]
        # get system state for each MC step
        self.system_state_at_step = dict()
        for p in state_file_paths:
            basename = os.path.basename(p)
            # final state of MC sim
            if basename == "FinalPosArray.csv":
                MC_step_no = self.sim_params["Monte Carlo steps"]
            # initial state of MC sim
            elif basename == "PosArray.csv":
                MC_step_no = 0
            # intermediate step during MC sim
            else:
                MC_step_no = int(basename.strip("PosArray").strip(".csv"))
            particle_positions = []
            with open(os.path.join(lc_data_path, p)) as data_file:
                for line in data_file:
                    # get position of particle and store in list
                    particle_pos = tuple(float(x) for x in line.split(","))
                    particle_positions.append(particle_pos)
            self.system_state_at_step[MC_step_no] = particle_positions


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a liquid crystal system dataset")
    parser.add_argument("--data-path", help="path to dataset")

    args = parser.parse_args()

    lc_systems = []
    for _path_ in os.listdir(args.data_path):
        full_path = os.path.join(args.data_path, _path_, "instanceRun")
        lc_systems.append(LCDataset(lc_data_path=full_path))
        break
