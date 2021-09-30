import numpy as np


# TODO: eclipse-line intersection, eclipse-eclipse intersection
# TODO: implement in Python an O(N^2) FORTRAN routine for determining closest approach between two ellipses
# TODO: calculate the local packing density
# TODO: calculate Fourier transform, structure factor etc.. for an LC system snapshot,
#       compare to diffraction pattern that Dr. Chen showed: should be a continuous distribution
#       due to lack of crystalline symmetry in the LC distribution


# ellipse params: a -- long axis, b -- short axis, theta -- angle with x-axis
def compute_closest_approach(a1, b1, a2, b2, k1k2, k1d, k2d):
    # eccentricities
    e1 = np.sqrt(1 - (b1 ** 2) / (a1 ** 2))
    e2 = np.sqrt(1 - (b2 ** 2) / (a2 ** 2))
    # anisotropic scaling parameter
    eta = a1 / b1 - 1
    # components of A'
    a11 = ((b1 ** 2) / (b2 ** 2)) * (1 + 0.5 * (1 + k1k2) * (eta * (2 + eta) - (e2 ** 2) * (1 + eta * k1k2) ** 2))
    a22 = ((b1 ** 2) / (b2 ** 2)) * (1 + 0.5 * (1 - k1k2) * (eta * (2 + eta) - (e2 ** 2) * (1 - eta * k1k2) ** 2))
    a12 = ((b1 ** 2) / (b2 ** 2)) * 0.5 * np.sqrt(1 - (k1k2 ** 2)) * (
            eta * (2 + eta) + (e2 ** 2) * (1 - (eta * k1k2) ** 2))
    # eigenvalues of A'
    lambda1 = 0.5 * (a11 + a22) + 0.5 * np.sqrt((a11 - a22) ** 2 + 4 * (a12 ** 2))
    lambda2 = 0.5 * (a11 + a22) + 0.5 * np.sqrt((a11 - a22) ** 2 + 4 * (a12 ** 2))
    # major and minor axes of transformed ellipse
    b2p = 1 / np.sqrt(lambda1)
    a2p = 1 / np.sqrt(lambda2)

    deltap = (a2p / b2p) ** 2 - 1
    # if the angle between the two transformed major axes is small
    if abs(k1k2) == 1:
        if a11 > a22:
            kpmp = (1 / np.sqrt(1 - e1 * k1d ** 2)) * (b1 / a1) * k1d
        elif a11 < a22:
            kpmp = np.sqrt(1 - k1d ** 2) / np.sqrt(1 - e1 * k1d ** 2)
        else:
            pass
    else:
        kpmp = ((a12 / np.sqrt(1 + k1k2)) * (b1 / a1 * k1d + k2d + (b1 / a1 - 1) * k1d * k1k2) \
                + (lambda1 - a11) / np.sqrt(1 - k1k2) * (b1 / a1 * k1d - k2d - (b1 / a1 - 1) * k1d * k1k2)) \
               / (np.sqrt(2 * (a12 ** 2 + (lambda1 - a11) ** 2) * (1 - (e1 * k1d) ** 2)))

    if kpmp == 0 or deltap == 0:
        dp = a2p + 1
    else:
        # coefficients of quartic for q
        t = 1 / (kpmp ** 2) - 1
        A = -1 / (b2p ** 2 * (1 + t))
        B = -2 / (b2p * (1 + t + deltap))
        C = -t - (1 + deltap) ** 2 + (1 / b2p ** 2) * (1 + t + deltap * t)
        D = 2 / b2p * (1 + t) * (1 + deltap)
        E = (1 + t + deltap) * (1 + deltap)

        # solution for quartic
        alpha = -3 / 8 * (B / A) ** 2 + C / A
        beta = (B / A) ** 3 / 8 - (B / A) * (C / A) / 2 + D / A
        gamma = -3 / 256 * (B / A) ** 4 + C / A * (B / A) ** 2 / 16 - (B / A) * (D / A) / 4 + E / A

        if beta == 0:
            qq = -B / (4 * A) + np.sqrt((-alpha + np.sqrt(alpha ** 2 - 4 * gamma)) / 2)
        else:
            P = -alpha ** 2 / 12 - gamma
            Q = -alpha ** 3 / 108 + gamma * alpha / 3 - beta ** 2 / 8
            U = (- 0.5 * Q + np.sqrt(Q ** 2 / 4 + P ** 3 / 27)) ** (1 / 3)

            if abs(U) != 0:
                y = (-5 / 6) * alpha + U - P / (3 * U)
            else:
                y = (-5 / 6) * alpha - Q ** (1 / 3)

            qq = -B / (4 * A) + 0.5 * (np.sqrt(alpha + 2 * y) \
                                       + np.sqrt(-(3 * alpha + 2 * y + 2 * beta / np.sqrt(alpha + 2 * y))))

            # substitute for R'
            dp = np.sqrt((qq ** 2 - 1) / deltap * (1 + b2p * (1 + deltap) / qq) ** 2 \
                         + (1 - (qq ** 2 - 1) / deltap) * (1 + b2p / qq) ** 2)

    dist = dp * b1 / np.sqrt(1 - e1 * k1d ** 2)

    return


def compute_ellipse_wall():
    return


def run_test_cases():
    # a_1, b_1, a_2, b_2, theta_1, theta_2
    test_pair_1 = {"inputs": (3, 1, 5, 1, 0, 0), "output": 8}
    test_pair_2 = {"inputs": (5, 1, 2, 0.25, np.pi / 2, 0), "output": 1.25}
    test_pair_3 = {"inputs": (9, 2, 5, 0.25, 5 * np.pi / 6, np.pi / 3), "output": 7}

    return


def main():
    return


if __name__ == "__main__":
    main()
