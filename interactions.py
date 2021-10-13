import numpy as np
import unittest

import math


# TODO: eclipse-line intersection, eclipse-eclipse intersection
# TODO: implement in Python an O(N^2) FORTRAN routine for determining closest approach between two ellipses
# TODO: calculate the local packing density
# TODO: calculate Fourier transform, structure factor etc.. for an LC system snapshot,
#       compare to diffraction pattern that Dr. Chen showed: should be a continuous distribution
#       due to lack of crystalline symmetry in the LC distribution


def determine_ellipse_overlap(x1, y1, x2, y2, theta1, theta2):
    # vector joining centers of ellipses
    d_vec = np.array([x1 - x2, y1 - y2], dtype=np.float16)
    # normalized major axis vectors
    k1 = np.array([np.cos(theta1), np.sin(theta1)], dtype=np.float16)
    k2 = np.array([np.cos(theta2), np.sin(theta2)], dtype=np.float16)
    # length between centers of ellipses
    d = np.linalg.norm(d_vec)
    # calculate cosines
    d_vec_n = d_vec / d
    k1d = k1 @ d_vec_n
    k2d = k2 @ d_vec_n
    k1k2 = k1 @ k2

    # should get output of 8
    a1, b1, a2, b2, k1d, k2d, k1k2 = (3, 1, 5, 1, 1, -1, 1)
    d_closest = compute_ellipse_ellipse_closest_approach(a1, b1, a2, b2, k1k2, k1d, k2d)

    return d_closest


# ellipse params: a -- long axis, b -- short axis, theta -- angle with x-axis
def compute_ellipse_ellipse_closest_approach(a1, b1, a2, b2, k1k2, k1d, k2d):
    """
    calculate distance of closest approach between two ellipses along direction of line joining their centers
    :param a1: length of major axis of untransformed ellipse 1
    :param b1: length of minor axis of untransformed ellipse 1
    :param a2: length of major axis of untransformed ellipse 2
    :param b2: length of minor axis of untransformed ellipse 2
    :param k1k2: cosine of angle between major axis of ellipse 1 and major axis of ellipse 2
    :param k1d: cosine of angle between major axis of ellipse 1 and direction of line joining centers
    :param k2d: cosine of angle between major axis of ellipse 2 and direction of line joining centers
    :return: distance of closest approach
    """
    # squared eccentricities
    e1 = 1 - b1 ** 2 / a1 ** 2
    e2 = 1 - b2 ** 2 / a2 ** 2
    # anisotropic scaling parameter
    eta = a1 / b1 - 1
    # components of A'
    a11 = (b1 / b2) ** 2 * (1 + 0.5 * (1 + k1k2) * (eta * (2 + eta) - e2 * (1 + eta * k1k2) ** 2))
    a22 = (b1 / b2) ** 2 * (1 + 0.5 * (1 - k1k2) * (eta * (2 + eta) - e2 * (1 - eta * k1k2) ** 2))
    a12 = (b1 / b2) ** 2 * 0.5 * np.sqrt(1 - k1k2 ** 2) * (eta * (2 + eta) + e2 * (1 - (eta * k1k2) ** 2))
    # eigenvalues of A'
    lambda1 = 0.5 * (a11 + a22) + 0.5 * np.sqrt((a11 - a22) ** 2 + 4 * (a12 ** 2))
    lambda2 = 0.5 * (a11 + a22) - 0.5 * np.sqrt((a11 - a22) ** 2 + 4 * (a12 ** 2))
    # major and minor axes lengths of transformed ellipse
    b2p = 1 / np.sqrt(lambda1)
    a2p = 1 / np.sqrt(lambda2)

    delta = (a2p / b2p) ** 2 - 1
    # if the angle between the transformed major axes is small
    if abs(k1k2) == 1:
        if a11 > a22:
            kpmpd = np.sqrt(1 - e1 * k1d ** 2) * b1 / a1 * k1d
        elif a11 < a22:
            kpmpd = np.sqrt(1 - k1d ** 2) / np.sqrt(1 - e1 * k1d ** 2)
        else:
            kpmpd = 0

    # determine k1d in transformed coordinate system
    else:
        kpmpd = (a12 / (1 + k1k2) ** .5 * (b1 / a1 * k1d + k2d + (b1 / a1 - 1) * k1d * k1k2)
                 + (lambda1 - a11) / np.sqrt(1 - k1k2) * (b1 / a1 * k1d - k2d - (b1 / a1 - 1) * k1d * k1k2)) \
                / ((2 * (a12 ** 2 + (lambda1 - a11) ** 2) * (1 - (e1 * k1d) ** 2)) ** .5)

    if abs(kpmpd) <= 1e-7 or abs(delta) <= 1e-7:
        dp = a2p + 1
    else:
        # coefficients of quartic for q
        t = 1 / kpmpd ** 2 - 1
        A = -1 / b2p ** 2 * (1 + t)
        B = -2 / b2p * (1 + t + delta)
        C = -t - (1 + delta) ** 2 + 1 / b2p ** 2 * (1 + t + delta * t)
        D = 2 / b2p * (1 + t) * (1 + delta)
        E = (1 + t + delta) * (1 + delta)

        # solution for quartic
        alpha = -3 / 8 * (B / A) ** 2 + C / A
        beta = (B / A) ** 3 / 8 - (B / A) * (C / A) / 2 + D / A
        gamma = -3 / 256 * (B / A) ** 4 + (C / A) * (B / A) ** 2 / 16 - (B / A) * (D / A) / 4 + E / A

        # beta = 0 degenerate case
        if abs(beta) <= 1e-7:
            qq = -B / (4 * A) + ((-alpha + (alpha ** 2 - 4 * gamma) ** .5) / 2) ** .5
        # beta != 0 general case
        else:
            P = -alpha ** 2 / 12 - gamma
            Q = -alpha ** 3 / 108 + gamma * alpha / 3 - beta ** 2 / 8
            if abs(Q ** 2 / 4 + P ** 3 / 27) < 1e-15:
                U = 0
            else:
                U = np.cbrt(-.5 * Q + (Q ** 2 / 4 + P ** 3 / 27) ** .5)

            if abs(U) != 0:
                y = (-5 / 6) * alpha + U - P / (3 * U)
            else:
                y = (-5 / 6 * alpha - Q) ** (1 / 3)

            qq = -B / (4 * A) + .5 * (np.sqrt(alpha + 2 * y)
                                      + np.sqrt(-(3 * alpha + 2 * y + 2 * beta / np.sqrt(alpha + 2 * y))))

        # substitute for R'
        dp = np.sqrt((qq ** 2 - 1) / delta * (1 + b2p * (1 + delta) / qq) ** 2
                     + (1 - (qq ** 2 - 1) / delta) * (1 + b2p / qq) ** 2)

    dist = dp * b1 / np.sqrt(1 - e1 * k1d ** 2)

    return round(dist, 3)


def compute_ellipse_line_intersection(theta, x_c, y_c, a, b, k, d):
    """
    Determine if an ellipse intersects a line
    :param theta: angle wrt x-axis
    :param x_c: x-coordinate shift in center of mass
    :param y_c: y-coordinate shift in center of mass
    :param a: major axis of ellipse
    :param b: minor axis of ellipse
    :param k: slope of line or 'inf' if line is vertical
    :param d: y-intercept of line or x-intercept if 'inf' flag is specified
    :return: true/false for intersection
    """
    # discriminant for determining intersection with arbitrary line: y = kx + d
    if k.isnumeric():
        disc = (d * k * math.cos(theta) ** 2 / b ** 2 + d * k * math.sin(theta) ** 2 / a ** 2 - k * y_c * math.cos(
            theta) / b ** 2 + k * x_c * math.sin(theta) / a ** 2 - d * math.cos(theta) * math.sin(
            theta) / a ** 2 + d * math.cos(theta) * math.sin(theta) / b ** 2 - x_c * math.cos(
            theta) / a ** 2 - y_c * math.sin(
            theta) / b ** 2) ** 2 - (
                       d ** 2 * math.cos(theta) ** 2 / b ** 2 + d ** 2 * math.sin(
                   theta) ** 2 / a ** 2 - 2 * d * y_c * math.cos(
                   theta) / b ** 2 + 2 * d * x_c * math.sin(
                   theta) / a ** 2 + x_c ** 2 / a ** 2 + y_c ** 2 / b ** 2 - 1) * (
                       k ** 2 * math.cos(theta) ** 2 / b ** 2 + k ** 2 * math.sin(
                   theta) ** 2 / a ** 2 - 2 * k * math.cos(
                   theta) * math.sin(theta) / a ** 2 + 2 * k * math.cos(theta) * math.sin(theta) / b ** 2 + math.cos(
                   theta) ** 2 / a ** 2 + math.sin(theta) ** 2 / b ** 2)
    # degenerate case where line is vertical x = d
    elif k == "inf":
        disc = (d * math.cos(theta) * math.sin(theta) / a ** 2 - d * math.cos(theta) * math.sin(theta) / b ** 2
                + y_c * math.cos(theta) / b ** 2 - x_c * math.sin(theta) / a ** 2) ** 2 \
               - (d ** 2 * math.cos(theta) ** 2 / a ** 2 + d ** 2 * math.sin(theta) ** 2 / b ** 2
                  - 2 * d * x_c * math.cos(theta) / a ** 2 - 2 * d * y_c * math.sin(theta) / b ** 2 + x_c ** 2 / a ** 2
                  + y_c ** 2 / b ** 2 - 1) * (math.cos(theta) ** 2 / b ** 2 + math.sin(theta) ** 2 / a ** 2)
    else:
        raise ValueError("k must be numeric or 'inf'")

    print(disc)

    if disc > 0:
        return True

    return False


def main():
    compute_ellipse_line_intersection(0.5, 0.1, 1, 2.2, 1, 'inf', 2.6)
    return


def cube_root(z):
    z = complex(z)
    x = z.real
    y = z.imag
    r = (x ** 2 + y ** 2) ** (1 / 6)
    theta = math.atan2(y, x) / 3
    cbrt = complex(r * math.cos(theta), r * math.sin(theta))

    return cbrt


class TestClosestApproachDistance(unittest.TestCase):

    def test_cube_root(self):
        test_inputs = [1, 10000, -300000, -7, 5.128193898, 69e5]
        for test_input in test_inputs:
            actual = cube_root(test_input)
            cube_real = round((actual ** 3).real, 8)
            cube_imag = round((actual ** 3).imag, 8)
            self.assertEqual(complex(test_input).real, cube_real)
            self.assertEqual(complex(test_input).imag, cube_imag)

    # inputs: a_1, b_1, a_2, b_2, theta_1, theta_2
    # expected output: distance of closest approach
    def test_closest_approach_distance(self):
        test_cases = [{"inputs": (3, 1, 5, 1, 0, 0), "expected": 8},
                      {"inputs": (5, 1, 2, 0.25, np.pi / 2, np.pi / 2), "expected": 1.25},
                      {"inputs": (9, 2, 5, 0.25, 5 * np.pi / 6, np.pi / 3), "expected": 7}]
        for case in test_cases:
            expected = case["expected"]
            a1, b1, a2, b2, theta1, theta2 = case["inputs"]
            actual = compute_ellipse_ellipse_closest_approach(a1=a1, b1=b1, a2=a2, b2=b2,
                                                              k1d=np.cos(theta1), k2d=np.cos(theta2),
                                                              k1k2=np.cos(theta1 - theta2))
            print(f"Expected: {expected}, Actual: {actual}")
            self.assertEqual(expected, actual)


if __name__ == "__main__":
    main()
