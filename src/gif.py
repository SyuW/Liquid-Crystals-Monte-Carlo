import os
import numpy as np
import matplotlib.pyplot as plt
import imageio


def create_gif(N, start, end, in_dir, save_name):

    step_finder = lambda path: float(path.strip(".png").split(sep="_")[-1])

    # build the gif
    frames = []
    for fname in sorted(list(os.listdir(in_dir)), key=step_finder):
        if start-1 <= step_finder(fname) <= end-1:
            frames.append(imageio.imread(os.path.join(in_dir, fname)))

    # save the frames as gif
    imageio.mimsave(save_name, frames, 'GIF', fps=3)


if __name__ == "__main__":
    # parameters
    inner_radius = 0
    num_particles = 50
    start_mc_step = 1000000
    end_mc_step = 1500000

    # base
    _base_ = "C:\\Users\\Sam Yu\\Desktop\\School\\4A\\Phys_437A_Research_Project\\visualizations\\"
    # input
    _in_ = os.path.join(_base_, f"r={inner_radius}\\plots\\N={num_particles}")
    # output
    _out_ = os.path.join(_base_, f"r={inner_radius}\\gifs\\N={num_particles}")
    os.makedirs(_out_, exist_ok=True)
    _save_ = os.path.join(_out_, f"r_{inner_radius}_N_{num_particles}_start_{start_mc_step}_end_{end_mc_step}.gif")

    # create and save gif
    create_gif(N=num_particles, start=start_mc_step, end=end_mc_step, in_dir=_in_, save_name=_save_)
