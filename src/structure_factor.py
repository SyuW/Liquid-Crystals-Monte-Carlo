import argparse
import numpy as np
from numpy import sin, cos, heaviside, tan
import os

# custom imports
from interactions import compute_ellipse_line_intersection as el_int
from LiquidCrystalSystem import LCSystem


def sample_on_annulus(num_points, r1, r2, t1=0, t2=2 * np.pi):
    unif = np.random.uniform(low=0, high=1, size=num_points)

    # use inverse transform sampling for polar coords
    r = np.sqrt((r2 ** 2 - r1 ** 2) * unif + r1 ** 2)
    t = np.random.uniform(low=t1, high=t2, size=num_points)

    return r * cos(t), r * sin(t)


def calculate_radial_dist_func(r, snapshots, r_max, step1, step2, volume):
    dr = 0.1
    rs = np.arange(0.001, r_max, dr)
    possible_steps = [s for s in snapshots.keys() if step1 <= step2]

    pair_dist_func = []
    # loop over radii
    for r in rs:
        # loop over monte carlo steps for time average
        for step in possible_steps:
            time_avg_arr = []
            pos_array = snapshots[step]
            N = len(pos_array)
            spatial = np.vstack([pos[:2] for pos in pos_array])
            ref_particles_sum = 0
            for ref_pos in spatial:
                pair_distances = np.linalg.norm(spatial-ref_pos, axis=1)
                # count the number of particles that lie in a circular shell (r, r+dr) centered on ref particle
                contained_in_r = np.sum([1 if r < r_jk < r+dr else 0 for r_jk in pair_distances])
                ref_particles_sum += contained_in_r
            time_avg_arr.append(ref_particles_sum)
        norm_factor = (volume / N) * (1 / (2 * np.pi * r * dr)) * (1 / N)
        pair_dist_func.append(norm_factor * np.mean(time_avg_arr))

    return rs, pair_dist_func


def calculate_local_density(pos_array, inner_radius, outer_radius, a, b, f):
    # number of points per stratified sample
    thetas, tstep = np.linspace(0, 2 * np.pi, num=50, endpoint=False, retstep=True)
    radii, rstep = np.linspace(inner_radius, outer_radius, num=50, endpoint=False, retstep=True)

    # for sanity check
    total_points = 0
    total_count = 0

    local_packing_fractions = dict()

    for r1 in radii:
        r2 = r1 + rstep

        for t1 in thetas:
            t2 = t1 + tstep

            # generate sample
            n = f * int(np.ceil(2 * np.pi * (r1 + r2) * rstep / 2))
            xs, ys = sample_on_annulus(n, r1, r2, t1, t2)

            # slope of rays from theta values
            k1 = 'inf' if abs(t1 - np.pi / 2) < 1e-15 else tan(t1)
            k2 = 'inf' if abs(t2 - np.pi / 2) < 1e-15 else tan(t2)

            # initialize counter
            count = 0

            # check which points fall within an ellipse
            for ellipse_pos in pos_array:

                xc, yc, t = ellipse_pos

                # check if ellipse intersects with line
                #if not (el_int(t, xc, yc, a, b, k1, 0) or el_int(t, xc, yc, a, b, k2, 0)):
                #    continue

                point_check = lambda xp, yp: 1 - heaviside((cos(t) * (xp - xc) - sin(t) * (yp - yc)) ** 2 / a ** 2
                                                           + (sin(t) * (xp - xc) + cos(t) * (yp - yc)) ** 2 / b ** 2
                                                           - 1, 1)
                count += sum(point_check(xs, ys))

            local_packing_fractions[(r1, r2, t1, t2)] = count / n

            total_count += count
            total_points += n

    print(len(pos_array))

    print(f"Estimated: {total_count / total_points}")
    print(f"True: {len(pos_array) * a * b / (outer_radius ** 2 - inner_radius ** 2)}")

    return local_packing_fractions


if __name__ == "__main__":
    data_path = "C:\\Users\\Sam Yu\\Desktop\\School\\4A\\Phys_437A_Research_Project\\datasets\\r=0"
    # plot all snapshots for all system sizes in dataset
    _path_ = os.listdir(data_path)[-1]
    full_path = os.path.join(data_path, _path_, "instanceRun")
    lc = LCSystem(lc_data_path=full_path, confinement="Circle", verbose=True)

    _R = lc.sim_params["R"]
    _r = lc.sim_params["r"]
    _a = lc.sim_params["Semi Major Axis"]
    _b = lc.sim_params["Semi Minor Axis"]
    _pos_array = lc.snapshots[199999]

    calculate_local_density(_pos_array, _r, _R, _a, _b, f=3)

