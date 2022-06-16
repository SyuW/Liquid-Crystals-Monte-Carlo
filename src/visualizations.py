import os
import numpy as np
import pickle
import imageio

from matplotlib.cm import get_cmap
from matplotlib.collections import PatchCollection
from matplotlib.patches import Ellipse, Rectangle, Circle, Wedge
import matplotlib.pyplot as plt

from joblib import Parallel, delayed


def iterate_over_runs(root, gifs, parallel=False):
    sims = [os.path.join(root, d) for d in os.listdir(root)]

    if not gifs:
        if parallel:
            Parallel(n_jobs=4)(delayed(plot_system_states)(run, color_angles=True) for run in sims)
    else:
        if parallel:
            Parallel(n_jobs=4)(delayed(gif_test)(run, 0, 2000000) for run in sims)
        else:
            for run in sims:
                with open(os.path.join(run, "checkpoint.pickle"), "rb") as in_f:
                    params = pickle.load(in_f)
                R = params["R"]
                gif_test(run, 0, 2000000)


def plot_single_state(pos_array, mc_step, N, R, r, a, b, color_angles=False):

    fig, ax = plt.subplots()

    if color_angles:
        fig.set_size_inches(10, 12)
    else:
        fig.set_size_inches(10, 10)

    outer_circle = plt.Circle((0, 0), R, color="black", fill=False, linewidth=2.2)
    inner_circle = plt.Circle((0, 0), r, color="black", fill=False, linewidth=2.2)
    ax.add_patch(inner_circle)
    ax.add_patch(outer_circle)

    # color map
    color_map = get_cmap("twilight")

    p = []
    colors = []

    # add ellipses to figure
    for pos in pos_array:
        ellipse_angle = (180 / np.pi) * (pos[-1] % np.pi)
        if color_angles:
            fill_color = color_map(ellipse_angle / 180)
            particle = Ellipse(xy=pos[:-1], angle=ellipse_angle,
                               width=2 * b, height=2 * a, ec="black",
                               linewidth=1.7, fc=fill_color, fill=True)
            p.append(particle)
            colors.append(fill_color)
        else:
            particle = Ellipse(xy=pos[:-1], angle=ellipse_angle,
                               width=2 * b, height=2 * a,
                               linewidth=0.5, color='black', fill=False)
        ax.add_patch(particle)

    if color_angles:
        e = PatchCollection(p, cmap=color_map)
        e.set_array(colors)
        e.set_clim([0, np.pi])
        cbar = fig.colorbar(e, label="Angle wrt x-axis",
                            orientation="horizontal", ax=plt.gca(), pad=0.05)
        cbar.set_ticks([0, np.pi/2, np.pi])
        cbar.set_ticklabels(["0", r"$\pi/2$", r"$\pi$"])

    # x,y-axis limits
    circle_pad = 5
    ax.set_xlim(-R - circle_pad, R + circle_pad)
    ax.set_ylim(-R - circle_pad, R + circle_pad)
    ax.set_title(f"N={N:} b={b} k={b / a:.2f} R={R:.2f} r={r:.2f} "
                 f"step={mc_step}", size=20)

    return fig


def plot_ellipse_ends_only(pos_array1, pos_array2, com_pos_array, mc_step, N, R, r, a, b):

    fig, ax = plt.subplots()
    fig.set_size_inches(10, 10)

    outer_circle = plt.Circle((0, 0), R, color="black", fill=False, linewidth=2.2)
    inner_circle = plt.Circle((0, 0), r, color="black", fill=False, linewidth=2.2)
    ax.add_patch(inner_circle)
    ax.add_patch(outer_circle)

    ax.scatter(pos_array1[:, 0], pos_array1[:, 1], color="red", marker="*")
    ax.scatter(pos_array2[:, 0], pos_array2[:, 1], color="blue", marker="+")
    # ax.scatter(com_pos_array[:, 0], com_pos_array[:, 1], color="blue", marker="+")

    # x,y-axis limits
    circle_pad = 5
    ax.set_xlim(-R - circle_pad, R + circle_pad)
    ax.set_ylim(-R - circle_pad, R + circle_pad)
    ax.set_title(f"N={N:} b={b} k={b / a:.2f} R={R:.2f} r={r:.2f} "
                 f"step={mc_step}", size=20)

    return fig


def plot_single_state_centers(pos_array, mc_step, N, R, r, a, b):

    fig, ax = plt.subplots()
    fig.set_size_inches(10, 10)

    outer_circle = plt.Circle((0, 0), R, color="black", fill=False, linewidth=2.2)
    inner_circle = plt.Circle((0, 0), r, color="black", fill=False, linewidth=2.2)
    ax.add_patch(inner_circle)
    ax.add_patch(outer_circle)

    ax.scatter(pos_array[:, 0], pos_array[:, 1])

    # x,y-axis limits
    circle_pad = 5
    ax.set_xlim(-R - circle_pad, R + circle_pad)
    ax.set_ylim(-R - circle_pad, R + circle_pad)
    ax.set_title(f"N={N:} b={b} k={b / a:.2f} R={R:.2f} r={r:.2f} "
                 f"step={mc_step}", size=20)

    return fig


def plot_system_states(run_dir, color_angles=False):
    with open(os.path.join(run_dir, "checkpoint.pickle"), "rb") as f:
        params = pickle.load(f)

    a = params["a"]
    b = params["b"]
    outer_radius = params["R"]
    inner_radius = params["r"]
    num_ellipses = len(params["pos_array"])

    # make the directory for plots
    save_dir = os.path.join(run_dir, "plots")
    os.makedirs(save_dir, exist_ok=True)

    # create figure and axes
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 10)

    # plot each system snapshot
    for sf in os.listdir(run_dir):
        if sf.endswith(".csv"):

            # get MC step and coordinates array
            pos_array = np.loadtxt(os.path.join(run_dir, sf), delimiter=",", dtype=np.float32)
            mc_step = int(sf.split(".csv")[0].split('_')[1])

            outer_circle = plt.Circle((0, 0), outer_radius, color="black", fill=False, linewidth=2.2)
            inner_circle = plt.Circle((0, 0), inner_radius, color="black", fill=False, linewidth=2.2)
            ax.add_patch(inner_circle)
            ax.add_patch(outer_circle)

            # add ellipses to figure
            for pos in pos_array:
                ellipse_angle = (180 / np.pi) * (pos[-1] % np.pi)
                particle = Ellipse(xy=pos[:-1], angle=ellipse_angle,
                                   width=2 * b, height=2 * a,
                                   linewidth=0.5, color='black', fill=False)
                ax.add_patch(particle)

            # x,y-axis limits
            circle_pad = 5
            ax.set_xlim(-outer_radius - circle_pad, outer_radius + circle_pad)
            ax.set_ylim(-outer_radius - circle_pad, outer_radius + circle_pad)
            ax.set_title(f"N={num_ellipses:} b={b} k={b / a:.2f} R={outer_radius:.2f} r={inner_radius:.2f} "
                         f"step={mc_step}", size=20)

            fig.savefig(os.path.join(save_dir, f"step_{mc_step}.png"))
            ax.cla()

    plt.close('all')


def create_feature_vector_visualization(pos_array, probe, feature_particles, mc_step,
                                        inner_radius, outer_radius, a, b):
    """
    Create visualizations of feature vector construction from a given snapshot
    """
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 10)

    outer_circle = plt.Circle((0, 0), outer_radius, color="black", fill=False, linewidth=2.2)
    inner_circle = plt.Circle((0, 0), inner_radius, color="black", fill=False, linewidth=2.2)
    ax.add_patch(inner_circle)
    ax.add_patch(outer_circle)

    for pos in pos_array:
        ellipse_angle = (180 / np.pi) * (pos[-1] % np.pi)
        particle = Ellipse(xy=pos[:-1], angle=ellipse_angle,
                           width=2 * b, height=2 * a,
                           linewidth=0.5, color='black', fill=False)
        ax.add_patch(particle)

    # add the extra ellipses
    probe_particle = Ellipse(xy=probe[:-1], angle=(180 / np.pi) * (probe[-1] % np.pi),
                             width=2 * b, height=2 * a, linewidth=1.7, facecolor="red")
    ax.add_patch(probe_particle)
    for i, particle_pos in enumerate(feature_particles):
        color = 'cyan'
        particle = Ellipse(xy=particle_pos[:-1], angle=(180 / np.pi) * (particle_pos[-1] % np.pi),
                           width=2 * b, height=2 * a,
                           linewidth=1.7, facecolor=color)
        ax.add_patch(particle)
        ax.annotate(f"{i}", particle_pos[:2], color='g', weight='bold', fontsize=12, ha='center', va='center')

    # x,y-axis limits
    circle_pad = 5
    ax.set_xlim(-outer_radius - circle_pad, outer_radius + circle_pad)
    ax.set_ylim(-outer_radius - circle_pad, outer_radius + circle_pad)
    ax.set_title(f"N={len(pos_array)} step={mc_step } b={b} k={b / a:.2f} R={outer_radius:.2f} r={inner_radius:.2f}",
                 size=20)

    return fig


def create_gif(start, end, in_dir, save_name):
    step_finder = lambda path: float(path.strip(".png").split(sep="_")[-1])

    # build the gif
    frames = []
    for fname in sorted([f for f in os.listdir(in_dir) if f.endswith(".png")], key=step_finder):
        if start - 1 <= step_finder(fname) <= end - 1:
            frames.append(imageio.imread(os.path.join(in_dir, fname)))
    # save the frames as gif
    imageio.mimsave(save_name, frames, 'GIF', fps=6)


def gif_test(in_dir, start_mc_step, end_mc_step):

    with open(os.path.join(in_dir, "checkpoint.pickle"), "rb") as in_f:
        params = pickle.load(in_f)

    N = len(params["pos_array"])
    a = params["a"]
    b = params["b"]
    R = params["R"]
    r = params["r"]
    phi = N*a*b / (R**2 - r**2)

    # make the output directories for GIFs
    _out_ = f"C:\\Users\\Sam Yu\\Desktop\\School\\4B\\Phys_437B_RP\\visualizations\\r={r/R:.2f}R\\gifs"
    os.makedirs(_out_, exist_ok=True)

    # save path
    _save_ = os.path.join(_out_, f"phi_{phi:.2f}_start_{start_mc_step}_end_{end_mc_step}.gif")

    # create and save gif
    create_gif(start=start_mc_step, end=end_mc_step, in_dir=os.path.join(in_dir, "plots"), save_name=_save_)


if __name__ == "__main__":
    dataset = "C:\\Users\\Sam Yu\\Desktop\\School\\4B\\Phys_437B_RP\\datasets\\r=0.4R\\3_30_2022"
    # gif_test(in_dir="C:\\Users\\Sam Yu\\Desktop\\School\\4B\\Phys_437B_RP\\datasets\\r=0\\3_6_2022\\"
    #                "R_30.1_r_0_N_500_STEPS_2000000_A_0.25_B_5_TARGET_0.0500000",
    #         R=0, start_mc_step=0, end_mc_step=1000000)
    iterate_over_runs(root=dataset, gifs=True, parallel=True)
    # plot_system_states(run_dir=os.path.join(dataset, "R_30.1_r_0_N_500_STEPS_2000000_A_0.25_B_5_TARGET_0.0500000"))
