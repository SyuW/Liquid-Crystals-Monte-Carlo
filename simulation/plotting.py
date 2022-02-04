import argparse
import os
import pickle

import numpy as np

import matplotlib
matplotlib.use('Agg')

from matplotlib.patches import Ellipse, Rectangle, Circle, Wedge
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.cm import get_cmap

from joblib import Parallel, delayed


def iterate_over_runs(root, parallel=False):
    
    runs = [os.path.join(root, d) for d in os.listdir(root)]
    
    if parallel:
        Parallel(n_jobs=2)(delayed(plot_system_states)(r, color_angles=True) for r in runs)
        

def plot_system_states(run_dir, color_angles=False):
    
    with open(os.path.join(run_dir, "checkpoint.pickle"), "rb") as f:
        params = pickle.load(f)
    
    a = params["a"]
    b = params["b"]
    outer_radius = params["R"]
    inner_radius = params["r"]
    
    # make the directory for plots
    save_dir = os.path.join(run_dir, "plots")
    os.makedirs(save_dir)
    
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
            
            fig.savefig(os.path.join(save_dir, f"step_{mc_step}.png"))
            ax.cla()
      
    plt.close('all')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plotting functions")
    parser.add_argument("--root-path", help="Path containing simulation results")
    args = parser.parse_args()
    
    import time
   
    start = time.perf_counter()
    plot_system_states(args.root_path, color_angles=True)
    print(f"Plotting execution time: {time.perf_counter() - start} seconds")