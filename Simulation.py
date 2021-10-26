import numpy as np
import os
import sys
import math

R = 25
r = .5
a = 0.25
b = 5
monte_carlo_steps = 2000000
step_xy = 0.5 * R
step_th = np.pi / 2
delta_y = 0  # % of a
delta_x = 0  # % of a


class InitialState:

    def __init__(self):
        return


class MonteCarloSim:

    def run_sim(self):
        moves = 0
        accepted_moves = 0

    def monte_carlo_step(self, positions):

        # minimum shifts in position and angle
        min_position_shift = .05 * a
        min_angular_shift = .025

        # counts
        moves = 0
        accepted_moves = 0
        fix_count = 0
        plot_count = 0

        # single Monte Carlo step
        for i, particle in enumerate(positions):
            pass

        return

    def __init__(self, initial_positions, total_steps):
        out_path = os.path.join(os.getcwd(), "data")
        out_name = ""

        for step in range(total_steps):
            self.monte_carlo_step()




if __name__ == "__main__":
    mc = MonteCarloSim()
