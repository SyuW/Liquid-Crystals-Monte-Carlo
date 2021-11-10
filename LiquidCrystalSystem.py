import argparse
import numpy as np
import os

from matplotlib.cm import get_cmap
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse


def get_configuration_information(path, confinement="Annulus", verbose=False):
    """
    Get information from the simulation
    :param verbose:
    :param confinement:
    :param path:
    :return:
    """
    params = dict()
    params["Confinement Type"] = confinement
    sim_notes_file = f"MonteCarlo_{confinement}_SimNotes.txt"
    # load information from the simulation notes file if it exists
    if os.path.exists(os.path.join(path, sim_notes_file)):
        with open(os.path.join(path, sim_notes_file)) as file:
            for line in file.readlines()[2:]:
                info = line.strip().split(": ")
                # special case of acceptance rate
                if info[0] == "Acceptance Rate":
                    params[info[0]] = float(info[1].strip(" %")) / 100
                else:
                    params[info[0]] = float(info[1])
    # if it does not exist then move on
    else:
        if verbose:
            print(f"Warning: {sim_notes_file} does not exist at {path}\n")
            pass

    return params


class LCSystem:

    def plot_snapshot(self, mc_step, extra_particles=[], color_angles=False):
        """
        Plot a snapshot of the system at a specific Monte Carlo step
        :param mc_step: Monte Carlo step
        :param extra_particles: extra particles to render
        :return: fig: figure of plot
        """
        assert (mc_step in self.snapshots.keys()), f"Monte Carlo step {mc_step} is not valid"

        # create figure and axes
        with plt.ioff():
            fig, ax = plt.subplots()

        if color_angles:
            fig.set_size_inches(10, 12)
        else:
            fig.set_size_inches(10, 10)

        # add annular boundaries to plot
        outer_circle = plt.Circle((0, 0), self.sim_params["R"], color='black', fill=False, linewidth=2.2)
        inner_circle = plt.Circle((0, 0), self.sim_params["r"], color='black', fill=False, linewidth=2.2)
        ax.add_patch(inner_circle)
        ax.add_patch(outer_circle)

        # plot liquid crystal ellipses
        b = self.sim_params['Semi Major Axis']
        a = self.sim_params['Semi Minor Axis']

        # color map
        color_map = get_cmap("twilight")

        # add ellipses
        p = []
        colors = []
        for particle_pos in self.snapshots[mc_step]:
            ellipse_angle = (180 / np.pi) * (particle_pos[-1] % np.pi)
            if color_angles:
                fill_color = color_map(ellipse_angle / 180)
                particle = Ellipse(xy=particle_pos[:-1], angle=ellipse_angle,
                                   width=2 * b, height=2 * a, ec="black",
                                   linewidth=1.7, fc=fill_color, fill=True)
                p.append(particle)
                colors.append(fill_color)
            else:
                particle = Ellipse(xy=particle_pos[:-1], angle=ellipse_angle,
                                   width=2 * b, height=2 * a,
                                   linewidth=1.7, color='black', fill=False)
            ax.add_patch(particle)

        if color_angles:
            e = PatchCollection(p, cmap=color_map)
            e.set_array(colors)
            cbar = fig.colorbar(e, ticks=[0, 0.5, 1], label="Angle wrt x-axis",
                                orientation="horizontal", pad=0.05)
            cbar.ax.set_xticklabels(["0", r"$\pi/2$", r"$\pi$"])

        # add the extra ellipses
        for i, particle_pos in enumerate(extra_particles):
            if i == 0:
                color = 'red'
            else:
                color = 'cyan'
            particle = Ellipse(xy=particle_pos[:-1], angle=(180 / np.pi) * (particle_pos[-1] % np.pi),
                               width=2 * b, height=2 * a,
                               linewidth=1.7, facecolor=color)
            ax.add_patch(particle)
            ax.annotate(f"{i}", particle_pos[:2], color='g', weight='bold', fontsize=12, ha='center', va='center')

        num_ellipses = self.sim_params["# of Ellipse"]
        outer_radius = self.sim_params["R"]
        inner_radius = self.sim_params["r"]

        # x,y-axis limits
        circle_pad = 5
        ax.set_xlim(-outer_radius - circle_pad, outer_radius + circle_pad)
        ax.set_ylim(-outer_radius - circle_pad, outer_radius + circle_pad)

        # set title, tweaks to font sizes
        ax.set_title(f'N={num_ellipses} b={b} k={b / a} R={outer_radius} r={inner_radius}', size=20)
        ax.tick_params(axis='both', labelsize=20)

        return fig

    def plot_batch_snapshots(self, save_dir, limit=50):
        """
        Plot a batch of snapshots at once
        :param save_dir: directory to save at
        :param limit: snapshot limit
        :return:
        """
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        count = 0
        for step in sorted(self.snapshots.keys()):
            save_path = os.path.join(save_dir, f"system_snapshot_at_step_{step}.png")
            snapshot_plot = self.plot_snapshot(step)
            snapshot_plot.savefig(save_path)

            count += 1
            if count > limit:
                break

    def __init__(self, lc_data_path, confinement="Annulus"):
        # get parameters from simulation
        self.sim_params = get_configuration_information(lc_data_path, confinement=confinement)

        if confinement == "Circle":
            self.sim_params["r"] = 0

        # get paths of all data files and retrieve system state information for each Monte Carlo step
        pos_array_paths = [os.path.join(lc_data_path, p) for p in os.listdir(lc_data_path) if p.endswith(".csv")]
        pos_array_paths = [os.path.normpath(p) for p in pos_array_paths]

        # get system state for each MC step
        self.snapshots = dict()
        for p in pos_array_paths:
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

            self.snapshots[MC_step_no] = np.loadtxt(p, delimiter=",", dtype=np.float32)

            # add number of particles to dictionary if doesn't already exist
            if "# of Ellipse" not in self.sim_params.keys():
                self.sim_params["# of Ellipse"] = len(self.snapshots[MC_step_no])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a liquid crystal system dataset")
    parser.add_argument("--data-path", help="path to dataset")

    data_path = "C:\\Users\\Sam Yu\\Desktop\\School\\4A\\Phys_437A_Research_Project\\datasets\\r=14"

    args = parser.parse_args()

    # plot all snapshots for all system sizes in dataset
    _path_ = os.listdir(args.data_path)[-1]
    full_path = os.path.join(args.data_path, _path_, "instanceRun")
    lc = LCSystem(lc_data_path=full_path)
    '''
    for _path_ in os.listdir(args.data_path):
        full_path = os.path.join(args.data_path, _path_, "instanceRun")
        lc = LCSystem(lc_data_path=full_path)

        # only do something for the first dataset
        break
    '''
    pass
