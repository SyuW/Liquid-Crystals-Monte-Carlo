import argparse
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos, heaviside, tan
from scipy.spatial import distance_matrix
import os

# custom imports
from interactions import compute_ellipse_line_intersection as el_int
from LiquidCrystalSystem import LCSystem


def sample_on_annulus(num_points, r1, r2, t1=0, t2=2 * np.pi, ret_polar=False):
    unif = np.random.uniform(low=0, high=1, size=num_points)

    # use inverse transform sampling for polar coords
    r = np.sqrt((r2 ** 2 - r1 ** 2) * unif + r1 ** 2)
    t = np.random.uniform(low=t1, high=t2, size=num_points)

    if ret_polar:
        return t, r
    else:
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
            pairwise_distances = distance_matrix(spatial, spatial)
            ref_particles_sum = 0
            for ref_pos in spatial:
                pair_distances = np.linalg.norm(spatial - ref_pos, axis=1)
                # count the number of particles that lie in a circular shell (r, r+dr) centered on ref particle
                contained_in_r = np.sum([1 if r < r_jk < r + dr else 0 for r_jk in pair_distances])
                ref_particles_sum += contained_in_r
            time_avg_arr.append(ref_particles_sum)
        norm_factor = (volume / N) * (1 / (2 * np.pi * r * dr)) * (1 / N)
        pair_dist_func.append(norm_factor * np.mean(time_avg_arr))

    return rs, pair_dist_func


def compute_structure_factor(pos_array):
    """

    :param pos_array:
    :return:
    """

    ks = np.arange(0, 12, 0.01)
    directions = np.arange(0, 2 * np.pi, 0.1)
    kx = np.outer(ks, np.cos(directions))
    ky = np.outer(ks, np.sin(directions))
    n = len(pos_array)

    exp_sum = 0
    for x, y in pos_array[:, :2]:
        exp_sum += np.exp(1j * (kx * x + ky * y))

    structure_factor = np.mean(np.real(exp_sum * np.conjugate(exp_sum) / n), axis=1)

    return structure_factor


def get_ensemble_average(lc_system, start=999999, end=1500000):
    avg = 0
    ensemble_size = 0
    for mc_step, pos_array in lc_system.snapshots.items():
        if start <= mc_step <= end:
            avg += compute_structure_factor(pos_array)
            ensemble_size += 1

    return avg / ensemble_size


def monte_carlo_integration_test(pos_array, inner_radius, outer_radius, a, b, f=5):

    # number of divisions for each strata
    t_divs, tstep = np.linspace(0, 2*np.pi, num=50, endpoint=True, retstep=True)
    r_divs, rstep = np.linspace(inner_radius, outer_radius, num=50, endpoint=False, retstep=True)

    # sample a bunch of points
    # ts, rs = sample_on_annulus(num_points, r1=inner_radius, r2=outer_radius, ret_polar=True)
    # xs = rs * cos(ts)
    # ys = rs * sin(ts)

    # for sanity check
    total_ellipse_count = 0
    total_point_count = 0

    # polar coordinate grid
    local_packing_fractions = np.zeros((len(t_divs), len(r_divs)))

    for i, r1 in enumerate(r_divs):
        r2 = r1 + rstep
        print(r1)
        for j, t1 in enumerate(t_divs):
            t1 = t1 % (2 * np.pi)
            t2 = (t1 + tstep) % (2 * np.pi)

            # count number of points that fall within the sector
            # r_sector_inds = np.where((rs > r1) & (rs <= r2))
            # t_sector_inds = np.where((ts > t1) & (ts <= t2))
            # common_inds = np.intersect1d(r_sector_inds, t_sector_inds)
            # sector_count = len(common_inds)
            # if sector_count == 0:
            #     continue

            # r_sector = rs[common_inds]
            # t_sector = ts[common_inds]
            # x_sector = r_sector * cos(t_sector)
            # y_sector = r_sector * sin(t_sector)
            n = f * int(np.ceil(2 * np.pi * (r1 + r2) * rstep / 2))
            xs, ys = sample_on_annulus(n, r1, r2, t1, t2)

            # initialize counter
            ellipse_count = 0

            # check which points fall within an ellipse
            for ellipse_pos in pos_array:
                xc, yc, t = ellipse_pos

                point_check = lambda xp, yp: 1 - heaviside((cos(t) * (xp - xc) - sin(t) * (yp - yc)) ** 2 / b ** 2
                                                           + (sin(t) * (xp - xc) + cos(t) * (yp - yc)) ** 2 / a ** 2
                                                           - 1, 1)
                ellipse_count += sum(point_check(xs, ys))

            local_packing_fraction = ellipse_count / n
            local_packing_fractions[j, i] = local_packing_fraction

            total_ellipse_count += ellipse_count
            total_point_count += n

    print(f"Estimated: {total_ellipse_count / total_point_count}")
    print(f"True: {len(pos_array) * a * b / (outer_radius ** 2 - inner_radius ** 2)}")

    fig, ax = plt.subplots(subplot_kw=dict(projection="polar"))
    r, theta = np.meshgrid(r_divs, t_divs)
    c = ax.contourf(theta, r, local_packing_fractions)
    ax.grid(False)
    fig.colorbar(c, orientation='vertical')

    return fig


if __name__ == "__main__":
    import time
    start = time.perf_counter()
    data_path = "C:\\Users\\Sam Yu\\Desktop\\School\\4A\\Phys_437A_Research_Project\\datasets\\r=0"
    #data_path = "C:\\Users\\Sam Yu\\Desktop\\School\\4A\\Phys_437A_Research_Project\\datasets\\k=1"
    data_path = "C:\\Users\\Sam Yu\\Desktop\\School\\4B\\Phys_437B_RP\\datasets\\DataToSam\\instanceRun2"
    N = 200
    #_path_ = [p for p in os.listdir(data_path) if f"n_{N}_" in p][0]
    #full_path = os.path.join(data_path, _path_)
    full_path = "C:\\Users\\Sam Yu\\Desktop\\School\\4B\\Phys_437B_RP\\" \
                "datasets\\r=0\\3_6_2022\\R_30.1_r_0_N_500_STEPS_2000000_A_0.25_B_5_TARGET_0.1325000"
    lc = LCSystem(lc_data_path=full_path, confinement="Circle", verbose=True)
    lc.sim_params["R"] = 68.68028197434451
    lc.sim_params["r"] = 0
    lc.sim_params['Semi Major Axis'] = 5
    lc.sim_params['Semi Minor Axis'] = 0.25
    lpf = monte_carlo_integration_test(pos_array=lc.snapshots[800000], inner_radius=0, outer_radius=25,
                                       f=2, a=0.25, b=5)
    print(f"{time.perf_counter() - start} seconds to calculate the local packing fraction")

    print(lpf)
    # _R = lc.sim_params["R"]
    #_r = lc.sim_params["r"]
    # _a = lc.sim_params["Semi Major Axis"]
    # _b = lc.sim_params["Semi Minor Axis"]
    # _pos_array = lc.snapshots[199999]

    # monte_carlo_integration_test(_pos_array, _r, _R, _a, _b, f=3)
    # averaged_structure_factor = get_ensemble_average(lc, start=0, end=1000000)
    # subtract out the delta function
    # averaged_structure_factor

    # k_mags = np.arange(0, 12, 0.01)
    # from scipy.signal import argrelextrema

    # local_maxima = argrelextrema(averaged_structure_factor, np.greater)
    # k_peaks = [k_mags[x] for x in local_maxima]
    # print(k_peaks)

    # import matplotlib.pyplot as plt

    # plt.plot(k_mags, averaged_structure_factor)
    # plt.gca().axhline(1, linestyle=':', color='black', label="High k-limit")
    # plt.gca().axvline(0.25, linestyle=':', color="orange", label="b")
    # plt.gca().axvline(0.25*2, linestyle=':', color="green", label="2b")
    # plt.gca().axvline(5, linestyle=':', color="red", label="a")
    # plt.gca().axvline(5*2, linestyle=":", color="blue", label="2a")
    # plt.legend()
    # plt.show()

    print("hello")
    pass
