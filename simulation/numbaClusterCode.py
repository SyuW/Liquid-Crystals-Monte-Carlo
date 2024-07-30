import numpy as np
from numpy.linalg.linalg import norm
import numpy.random as rd
import os
from math import trunc
from operator import itemgetter
import sys
import math
import time

#from numba import jit


def dist(x1, y1, x2, y2):
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)


def atan2(x, y):
    if x > 0:
        val = np.arctan(y / x)
    elif x < 0 and y >= 0:
        val = np.arctan(y / x) + np.pi
    elif x < 0 and y < 0:
        val = np.arctan(y / x) - np.pi
    elif x == 0 and y > 0:
        val = np.pi / 2
    elif x == 0 and y < 0:
        val = np.pi / 2

    return val


#@jit(cache=True)
def overlap_Ellipse(x1, y1, theta1, x2, y2, theta2, a, b):
    k = b / a
    xi = (k ** 2 - 1) / (k ** 2 + 1)
    w = 2 * a
    r_val = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
    r_vec = np.array([(x2 - x1) / r_val, (y2 - y1) / r_val])

    u1 = np.array([np.cos(theta1), np.sin(theta1)])
    u2 = np.array([np.cos(theta2), np.sin(theta2)])

    r_dot_u1 = np.dot(r_vec, u1)

    r_dot_u2 = np.dot(r_vec, u2)

    # u1_dot_u2 = np.cos(theta_2_prime)

    u1_dot_u2 = np.dot(u1, u2)

    denom = (1 - 0.5 * xi * (((r_dot_u1 + r_dot_u2) ** 2 / (1 + xi * u1_dot_u2)) + (
                (r_dot_u1 - r_dot_u2) ** 2 / (1 - xi * u1_dot_u2)))) ** (0.5)

    sigma_2D = w / denom

    if r_val >= sigma_2D:
        val = False
    else:
        val = True

        # print('sigma: ',sigma_2D,' r dist: ',r_dist,' xi: ',xi,' u1: ',u1, ' u2: ', u2, ' r vec: ',r)

    return val


#@jit(cache=True)
def overlap_Ellipse2(x1, y1, theta1, x2, y2, theta2, a, b):
    k = b / a
    xi = (k ** 2 - 1) / (k ** 2 + 1)
    w = 2 * a
    r_val = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
    r_vec = np.array([(x2 - x1) / r_val, (y2 - y1) / r_val])

    u1 = np.array([np.cos(theta1), np.sin(theta1)])
    u2 = np.array([np.cos(theta2), np.sin(theta2)])

    r_dot_u1 = np.dot(r_vec, u1)

    r_dot_u2 = np.dot(r_vec, u2)

    # u1_dot_u2 = np.cos(theta_2_prime)

    u1_dot_u2 = np.dot(u1, u2)

    denom = (1 - 0.5 * xi * (((r_dot_u1 + r_dot_u2) ** 2 / (1 + xi * u1_dot_u2)) + (
                (r_dot_u1 - r_dot_u2) ** 2 / (1 - xi * u1_dot_u2)))) ** (0.5)

    sigma_2D = w / denom

    # sigma_2D = contact(theta_r,phi,xi,w)

    return sigma_2D


#@jit(cache=True)
def GeometricPotential(x1, y1, theta1, x2, y2, theta2, LongAxis1, ShortAxis1, LongAxis2, ShortAxis2):
    small = 10 ** (-14)

    k1 = [np.cos(theta1), np.sin(theta1)]
    k2 = [np.cos(theta2), np.sin(theta2)]

    k1k1_dyad = [[np.cos(theta1) ** 2, np.sin(theta1) * np.cos(theta1)],
                 [np.sin(theta1) * np.cos(theta1), np.sin(theta1) ** 2]]
    k2k2_dyad = [[np.cos(theta2) ** 2, np.sin(theta2) * np.cos(theta2)],
                 [np.sin(theta2) * np.cos(theta2), np.sin(theta2) ** 2]]
    identity = [[1, 0], [0, 1]]

    e1 = np.sqrt(1 - (ShortAxis1 / LongAxis1) ** 2)
    e2 = np.sqrt(1 - (ShortAxis2 / LongAxis2) ** 2)

    eta = (LongAxis1 / ShortAxis1) - 1

    a_prime1 = (ShortAxis1 / ShortAxis2) ** 2 * np.add(np.array(identity), (eta * np.array(k1k1_dyad)))
    a_prime2 = np.add(np.array(identity), (-1) * (e2 ** 2) * np.array(k2k2_dyad))
    a_prime3 = np.add(np.array(identity), eta * np.array(k1k1_dyad))

    A_prime = np.matmul(np.matmul(np.array(a_prime1), np.array(a_prime2)), np.array(a_prime3))

    eigenVal, eigenVec = np.linalg.eig(A_prime)

    if eigenVal[0] > eigenVal[1]:
        LongAxis2_prime = 1 / (np.sqrt(eigenVal[1]))
        k_minus = eigenVec[:, 1]
        ShortAxis2_prime = 1 / (np.sqrt(eigenVal[0]))
        k_plus = eigenVec[:, 0]
    else:
        LongAxis2_prime = 1 / (np.sqrt(eigenVal[0]))
        k_minus = eigenVec[:, 0]
        ShortAxis2_prime = 1 / (np.sqrt(eigenVal[1]))
        k_plus = eigenVec[:, 1]

    d_hat = (1 / np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)) * np.array([x1 - x2, y1 - y2])

    d_hat_prime = (1 / (np.sqrt(1 - (e1 ** 2) * (np.dot(np.array(k1), d_hat)) ** 2))) * np.add(d_hat, (
                ShortAxis1 / LongAxis1 - 1) * (np.dot(np.array(k1), d_hat)) * np.array(k1))

    sinPhi = np.dot(k_minus, d_hat_prime)
    cosPhi = np.dot(k_plus, d_hat_prime)

    delta = (LongAxis2_prime / ShortAxis2_prime) ** 2 - 1

    if np.abs(cosPhi) < small:

        d_prime = 1 + LongAxis2_prime
    else:
        tanPhiSq = (sinPhi / cosPhi) ** 2

        A = (-1 / (ShortAxis2_prime ** 2)) * (1 + tanPhiSq)
        B = (-2 / ShortAxis2_prime) * (1 + tanPhiSq + delta)
        C = -tanPhiSq - (1 + delta) ** 2 + (1 / (ShortAxis2_prime ** 2)) * (1 + tanPhiSq * (1 + delta))
        D = (2 / ShortAxis2_prime) * (1 + tanPhiSq) * (1 + delta)
        E = (1 + tanPhiSq + delta) * (1 + delta)

        alpha = -(3 * B ** 2) / (8 * A ** 2) + C / A
        beta = B ** 3 / (8 * A ** 3) - (B * C) / (2 * A ** 2) + D / A
        gamma = -3 * B ** 4 / (256 * A ** 4) + (C * B ** 2) / (16 * A ** 3) - (B * D) / (4 * A ** 2) + (E / A)

        settings = np.seterr(divide='raise')

        P = - (alpha ** 2 / 12) - gamma
        Q = - (alpha ** 3 / 108) + (alpha * gamma / 3) - (beta ** 2 / 8)
        U = (-Q / 2 + np.sqrt((Q ** 2 / 4) + (P ** 3 / 27))) ** (1 / 3)

        if np.abs(U) < small or math.isnan(U) == True:
            y = -(5 * alpha / 6) - (Q ** (1 / 3))
        else:
            y = -(5 * alpha / 6) + U - P / (3 * U)

        q = -(B / (4 * A)) + (1 / 2) * (
                    np.sqrt(alpha + 2 * y) + np.sqrt(-(3 * alpha + 2 * y + (2 * beta / (np.sqrt(alpha + 2 * y))))))

        d_prime = np.sqrt(
            ((q ** 2 - 1) / delta) * (1 + (ShortAxis2_prime * (1 + delta) / q)) ** 2 + (1 - ((q ** 2 - 1) / delta)) * (
                        1 + ShortAxis2_prime / q) ** 2)

        d = d_prime * ShortAxis1 / (np.sqrt(1 - (e1 ** 2) * (np.dot(k1, d_hat) ** 2)))

    return d


#@jit(cache=True)
def HardBoundaryCircle_Disc(R, shortAxis, longAxis, xc, yc, theta):
    overlap = False

    cos = np.cos(theta)
    sin = np.sin(theta)
    sin2 = np.sin(2 * theta)

    A = longAxis ** 2 + xc ** 2 + yc ** 2 - 2 * longAxis * (xc * cos + yc * sin) - R ** 2
    B = 4 * shortAxis * (yc * cos - xc * sin)
    C = 4 * (shortAxis ** 2) - 2 * (longAxis ** 2) + 2 * (xc ** 2) + 2 * (yc ** 2) - 2 * (R ** 2)
    D = 4 * shortAxis * (yc * cos - xc * sin)
    E = longAxis ** 2 + 2 * longAxis * (xc * cos + yc * sin) + xc ** 2 + yc ** 2 - R ** 2

    delta = (256 * (A * E) ** 3
             - 192 * (B * D) * (A * E) ** 2
             - 128 * (A * C * E) ** 2
             + 144 * (C * E) * (A * D) ** 2
             - 27 * (A ** 2) * (D ** 4)
             + 144 * (A * C) * (B * E) ** 2
             - 6 * (A * E) * (B * D) ** 2
             - 80 * (A * B * D * E) * (C) ** 2
             + 18 * (A * B * C) * (D ** 3)
             + 16 * (A * E) * (C ** 4)
             - 4 * A * (C ** 3) * (D ** 2)
             - 27 * (E ** 2) * (B ** 4)
             + 18 * (C * D * E) * (B ** 3)
             - 4 * (B * D) ** 3
             - 4 * E * (B ** 2) * (C ** 3)
             + (B ** 2) * (C ** 2) * (D ** 2)
             )

    P = 8 * A * C - 3 * (B ** 2)

    R_d = B ** 3 + 8 * D * (A ** 2) - 4 * A * B * C

    delta_0 = C ** 2 - 3 * B * D + 12 * A * E

    D_d = 64 * (A ** 3) * E - 16 * (A ** 2) * (C ** 2) + 16 * A * (B ** 2) * C - 16 * (A ** 2) * B * D - 3 * (B ** 4)

    if delta < 0:
        overlap = True
    elif delta > 0 and (P < 0 and D_d < 0):
        overlap = True
    elif delta == 0:
        if P < 0 and D_d < 0 and delta_0 != 0:
            overlap = True
        elif D_d > 0 or (P > 0 and (D_d != 0 or R_d != 0)):
            overlap = False
        elif delta_0 == 0 and D_d != 0:
            overlap = True
        elif D_d == 0 and P < 0:
            overlap = True
        elif D_d == 0 and P > 0 and R_d == 0:
            overlap = False
        else:
            overlap = False

    print("----------------------------------")
    print(" ")
    print("Inputs:")
    print("Boundary radius", R)
    print("Ellipse center x-position:", xc)
    print("Ellipse center y-position:", yc)
    print("Ellipse rotation angle:", theta)
    print("Major axis:", longAxis)
    print("Minor axis:", shortAxis)
    print(" ")
    print("Outputs:")
    print("A:", A)
    print("B:", B)
    print("C:", C)
    print("D:", D)
    print("E:", E)
    print("P:", P)
    print("R_d:", R_d)
    print("D_d:", D_d)
    print("delta:", delta)
    print("delta_0:", delta_0)
    print("----------------------------------")

    return overlap


def init_Circ_H_Gr(n, a, b, R):
    v = np.ceil(R / a)
    h = np.ceil(R / b)

    N = int(v * h)

    xgrid = np.linspace(b - R, (b - R) + (int(h) - 1) * (2 * b), int(h))
    ygrid = np.linspace(a - R, (a - R) + (int(v) - 1) * (2 * a), int(v))
    thetaArray = np.full((N, 1), 0)

    X, Y = np.meshgrid(xgrid, ygrid)

    posArr = np.array([X.flatten(), Y.flatten()]).T

    posArr = np.c_[posArr, thetaArray]
    nDel = 0
    indexDel = []
    incorrect = []

    for p in range(len(posArr)):

        if posArr[p, 0] ** 2 + posArr[p, 1] ** 2 > R ** 2 or (
                HardBoundaryCircle_Disc(R, a, b, posArr[p, 0], posArr[p, 1], posArr[p, 2]) == True):
            indexDel.append(p)
            incorrect.append(posArr[p, :])
            nDel += 1
        else:
            pass

    redE = np.array(incorrect)

    c = 0

    for i in indexDel:
        newArr = np.delete(posArr, i - c, 0)
        posArr = newArr
        c += 1

    print(len(posArr))
    print(len(redE))
    newLength = len(posArr)

    if n >= newLength:
        pass
    elif n < newLength:
        d = 0

        while d != newLength - n:
            k = rd.randint(0, newLength - n - d)
            newArr = np.delete(posArr, k, 0)
            posArr = newArr
            d += 1

    print(len(posArr))

    return [posArr, redE]


def init_Circ_H_Rd(n, l, a, b):
    init_pos = np.zeros((n, 3))

    # Correct overlapping positions
    sortedInitPos = [[i, init_pos[i, 0], init_pos[i, 1], init_pos[i, 2], ] for i in range(n)]

    for x in range(n):

        valid = False
        overlapVar = False

        while valid == False:

            radius = rd.uniform(0, l)
            angle = rd.uniform(0, 2 * np.pi)
            init_pos[x, 0] = radius * np.cos(angle)
            init_pos[x, 1] = radius * np.sin(angle)
            init_pos[x, 2] = rd.uniform(0, 2 * np.pi)

            if (init_pos[x, 0] ** 2 + init_pos[x, 1] ** 2) > ((l - 2 * b) ** 2):
                if (init_pos[x, 0] ** 2 + init_pos[x, 1] ** 2) > (l ** 2) or (
                        init_pos[x, 0] ** 2 + init_pos[x, 1] ** 2) > (l - a) ** 2 or (
                        HardBoundaryCircle_Disc(l, a, b, init_pos[x, 0], init_pos[x, 1], init_pos[x, 2]) == True):
                    print('init - border overlap')
                    continue
                else:

                    for j in range(n):
                        if j != x and (overlap_Ellipse(init_pos[x, 0], init_pos[x, 1], init_pos[x, 2], init_pos[j, 0],
                                                       init_pos[j, 1], init_pos[j, 2], a, b) == True):
                            overlapVar = True
                            print('iter' + str(x) + ' init - particle overlap')
                            break
                        else:
                            print('iter' + str(x) + 'init - no particle overlap')
                            overlapVar = False

                    if overlapVar == True:
                        continue
                    else:
                        valid = True


            else:
                for j in range(n):
                    if j != x and (overlap_Ellipse(init_pos[x, 0], init_pos[x, 1], init_pos[x, 2], init_pos[j, 0],
                                                   init_pos[j, 1], init_pos[j, 2], a, b) == True):
                        overlapVar = True
                        print('iter' + str(x) + ' init - particle overlap')
                        break
                    else:
                        print('iter' + str(x) + 'init - no particle overlap')
                        overlapVar = False

                if overlapVar == True:
                    continue
                else:
                    valid = True

    return init_pos


#@jit(cache=True)
def MC_Circ_Hard(PosArray, d_pos, d_ang, steps, n, l, a, b):
    global moves
    global accepted_moves

    moves = 0
    accepted_moves = 0
    imageName = "step_figure.png"
    numE = len(PosArray)
    kVal = b / a
    file_name = "MonteCarlo_Circle_SimNotes.txt"

    ogdir = os.getcwd()

    main_folder_name = os.path.join(ogdir, f'circle_R{l}_n_{numE}_k_{kVal}_HardBC')
    folder_name = os.path.join(main_folder_name, 'instanceRun')

    if os.path.exists(main_folder_name):
        pass
    else:
        os.makedirs(main_folder_name)

    if os.path.exists(folder_name):
        expand = 0
        while True:
            expand += 1
            new_folder_name = folder_name + str(expand)
            if os.path.exists(new_folder_name):
                continue
            else:
                folder_name = new_folder_name
                break

    os.makedirs(folder_name)

    ##### Save Initial State

    fileNameArray = 'PosArray.csv'
    complete_name = os.path.join(folder_name, fileNameArray)
    np.savetxt(complete_name, PosArray, delimiter=',')

    stepsArr = np.linspace(0, steps, steps + 1)

    n = len(PosArray)

    fullImageName = os.path.join(folder_name, imageName)

    ###Plot initial state

    fileNameArray = 'PosArray.csv'
    complete_name = os.path.join(folder_name, fileNameArray)
    np.savetxt(complete_name, PosArray, delimiter=',')

    fixCount = 0
    plotCount = 0

    for u in range(steps):

        for w in range(n):

            if fixCount == (np.ceil(steps / 50)):
                fixCount = 0
                if accepted_moves / moves < 0.47:
                    d_ang = 0.9 * d_ang
                    d_pos = 0.9 * d_ang
                elif accepted_moves / moves < 0.37:
                    d_ang = 0.75 * d_ang
                    d_pos = 0.75 * d_ang
                elif accepted_moves / moves < 0.27:
                    d_ang = 0.6 * d_ang
                    d_pos = 0.6 * d_ang
                elif accepted_moves / moves < 0.17:
                    d_ang = 0.35 * d_ang
                    d_pos = 0.35 * d_ang
                elif accepted_moves / moves < 0.07:
                    d_ang = 0.2 * d_ang
                    d_pos = 0.2 * d_ang
                elif accepted_moves / moves > 0.57:
                    d_ang = 1.1 * d_ang
                    d_pos = 1.1 * d_ang
                elif accepted_moves / moves > 0.67:
                    d_ang = 1.35 * d_ang
                    d_pos = 1.35 * d_ang
                elif accepted_moves / moves > 0.77:
                    d_ang = 1.5 * d_ang
                    d_pos = 1.5 * d_ang
                elif accepted_moves / moves > 0.87:
                    d_ang = 1.65 * d_ang
                    d_pos = 1.65 * d_ang
                elif accepted_moves / moves > 0.97:
                    d_ang = 1.8 * d_ang
                    d_pos = 1.8 * d_ang

            x = d_pos * rd.uniform(-1, 1)
            y = d_pos * rd.uniform(-1, 1)
            t = d_ang * rd.uniform(-1, 1)

            # print(x,y)

            testX = PosArray[w, 0] + x
            testY = PosArray[w, 1] + y
            testT = PosArray[w, 2] + t

            if (testX ** 2 + testY ** 2) > ((l - 2 * b) ** 2):
                if (testX ** 2 + testY ** 2) > ((l - a) ** 2):
                    overlapVar = True
                    # print("border overlap")
                elif (HardBoundaryCircle_Disc(l, a, b, testX, testY, testT) == True):
                    overlapVar = True
                    # print("border overlap")
                else:
                    for j in range(n):

                        rij = np.sqrt((testX - PosArray[j, 0]) ** 2 + (testY - PosArray[j, 1]) ** 2)

                        if rij < (2 * b):
                            if j != w and (
                                    overlap_Ellipse(testX, testY, testT, PosArray[j, 0], PosArray[j, 1], PosArray[j, 2],
                                                    a, b) == True):
                                overlapVar = True
                                # print("particle overlap")
                                break
                            else:
                                overlapVar = False
                        else:
                            overlapVar = False
            else:
                for j in range(n):

                    rij = np.sqrt((testX - PosArray[j, 0]) ** 2 + (testY - PosArray[j, 1]) ** 2)

                    if rij < (2 * b):
                        if j != w and (
                                overlap_Ellipse(testX, testY, testT, PosArray[j, 0], PosArray[j, 1], PosArray[j, 2], a,
                                                b) == True):
                            overlapVar = True
                            # print("particle overlap")
                            break
                        else:
                            overlapVar = False
                    else:
                        overlapVar = False

            if overlapVar == True:
                pass
            elif overlapVar == False:
                accepted_moves += 1
                PosArray[w, 0] = testX
                PosArray[w, 1] = testY
                PosArray[w, 2] = testT

            moves += 1

            fixCount += 1

            x = 0
            y = 0

        plotCount += 1

        # if plotCount == (np.ceil(steps / 100)):
        #     plotCount = 0

        #     ######################################

        #     #### Periodically Save Snapshots #####

        #     ######################################

        #     fileNameArray = 'PosArray.csv'
        #     new_fileName = fileNameArray.split(".csv")[0] + str(u) + ".csv"
        #     complete_name = os.path.join(folder_name, new_fileName)
        #     np.savetxt(complete_name, PosArray, delimiter=',')

    reduced_density = n * np.pi * a * b / (np.pi * l ** 2)

    print(accepted_moves)
    print(moves)
    print(accepted_moves / moves)

    complete_name = os.path.join(folder_name, file_name)

    text_file = open(complete_name, 'w+')
    text_file.write("Parameters" + "\r\n")
    text_file.write("- - - - -" + "\r\n")
    text_file.write("Monte Carlo steps: " + str(steps) + "\r\n")
    text_file.write("R: " + str(l) + "\r\n")
    text_file.write("d_pos / step size: " + str(d_pos) + "\r\n")
    text_file.write("d_ang / step size: " + str(d_ang) + "\r\n")
    text_file.write("# of Ellipse: " + str(n) + "\r\n")
    text_file.write("reduced density: " + str(reduced_density) + "\r\n")
    text_file.write("Semi Minor Axis: " + str(a) + "\r\n")
    text_file.write("Semi Major Axis: " + str(b) + "\r\n")
    text_file.write("Accepted Moves: " + str(accepted_moves) + "\r\n")
    text_file.write("Total Moves: " + str(moves) + "\r\n")
    text_file.write("Acceptance Rate: " + str(100 * (accepted_moves / moves)) + " %" + "\r\n")

    text_file.close()

    fileNameArray = 'FinalPosArray.csv'
    complete_name = os.path.join(folder_name, fileNameArray)
    np.savetxt(complete_name, PosArray, delimiter=',')


def init_Ann_H_Gr(n, a, b, R, r2):
    v = np.ceil(R / a)
    h = np.ceil(R / b)

    N = int(v * h)

    xgrid = np.linspace(b - R, (b - R) + (int(h) - 1) * (2 * b), int(h))
    ygrid = np.linspace(a - R, (a - R) + (int(v) - 1) * (2 * a), int(v))
    thetaArray = np.full((N, 1), 0)

    X, Y = np.meshgrid(xgrid, ygrid)

    posArr = np.array([X.flatten(), Y.flatten()]).T

    posArr = np.c_[posArr, thetaArray]
    nDel = 0
    indexDel = []
    incorrect = []

    for p in range(len(posArr)):

        if posArr[p, 0] ** 2 + posArr[p, 1] ** 2 > R ** 2 or posArr[p, 0] ** 2 + posArr[p, 1] ** 2 < r2 ** 2 or (
                HardBoundaryCircle_Disc(R, a, b, posArr[p, 0], posArr[p, 1], posArr[p, 2]) == True) or (
                HardBoundaryCircle_Disc(r2, a, b, posArr[p, 0], posArr[p, 1], posArr[p, 2]) == True):
            indexDel.append(p)
            incorrect.append(posArr[p, :])
            nDel += 1
        else:
            pass

    redE = np.array(incorrect)

    c = 0

    for i in indexDel:
        newArr = np.delete(posArr, i - c, 0)
        posArr = newArr
        c += 1

    print(len(posArr))
    print(len(redE))
    newLength = len(posArr)

    if n >= newLength:
        pass
    elif n < newLength:
        d = 0

        while d != newLength - n:
            k = rd.randint(0, newLength - n - d)
            newArr = np.delete(posArr, k, 0)
            posArr = newArr
            d += 1

    print(len(posArr))

    return [posArr, redE]


def init_Ann_H_GrC(n, a, b, R, r2, d, e):
    ## When an exact contact distance method is reliable use that instead of xStar

    a_tilde = a * (1 + (d / 100))
    b_tilde = b + a * ((e / 100))

    Nv1 = int(np.floor(r2 / a_tilde))

    if Nv1 % 2 == 0:
        ygrid1 = np.linspace(-(Nv1 * a_tilde), Nv1 * a_tilde, Nv1)

    else:
        ygrid1 = np.linspace(-a_tilde * (Nv1 - 1), (Nv1 - 1) * a_tilde, Nv1)

    diffArray = []
    totalN = 0
    for k in range(Nv1):
        xStar = np.sqrt(r2 ** 2 - (ygrid1[k]) ** 2)
        diff = R - xStar
        Nh1 = 2 * int(np.floor(diff / (2 * b_tilde)))
        diffArray.append([xStar, Nh1])
        totalN += Nh1

    initalizedArray = np.zeros((1, 3))
    corrected = False
    countVar = 0
    for w in range(len(diffArray)):

        if diffArray[w][1] != 0:

            if diffArray[w][1] > 1:
                right = np.array(
                    [((2 * (s + 1) - 1) * b_tilde + diffArray[w][0]) for s in range(int(diffArray[w][1] / 2))])
                left = -1 * right
                xVal = np.concatenate((right, left), axis=0)
            else:
                xVal = np.array((diffArray[w][0] + b_tilde))

            # xVal=np.linspace(-chordArray[w][0]+b_tilde,chordArray[w][0]-b_tilde,chordArray[w][1])
            for vals in xVal:
                interimArray = np.append(initalizedArray, [[vals, ygrid1[w], 0]], axis=0)
                initalizedArray = interimArray
        else:
            pass

    Nv2_half = int(np.floor((R - r2) / (2 * a_tilde)))

    if Nv2_half != 0:

        up = np.linspace(r2 + 2 * a_tilde, r2 + (2 * Nv2_half) * a_tilde, Nv2_half)
        down = -1 * up
        ygrid2 = np.concatenate((up, down), axis=0)

        chordArray = []
        for k in range(2 * Nv2_half):
            chord = np.sqrt(R ** 2 - (ygrid2[k]) ** 2)
            NC = int(np.floor(2 * chord / (2 * b_tilde)))
            chordArray.append([chord, NC])
            totalN += NC

        corrected = False
        countVar = 0
        for w in range(len(chordArray)):

            if chordArray[w][1] != 0:

                if chordArray[w][1] > 1:
                    if chordArray[w][1] % 2 == 0:
                        right = np.array([((2 * (s + 1) - 1) * b_tilde) for s in range(int(chordArray[w][1] / 2))])
                        left = -1 * right
                        xVal = np.concatenate((right, left), axis=0)
                    else:
                        right = np.array([((2 * (s + 1)) * b_tilde) for s in range(int(chordArray[w][1] / 2))])
                        left = -1 * right
                        xVal = np.concatenate((right, left), axis=0)
                        interimX = np.append(xVal, [0], axis=0)
                        xVal = interimX
                else:
                    xVal = np.array([0])

                # xVal=np.linspace(-chordArray[w][0]+b_tilde,chordArray[w][0]-b_tilde,chordArray[w][1])
                for vals in xVal:
                    interimArray = np.append(initalizedArray, [[vals, ygrid2[w], 0]], axis=0)
                    initalizedArray = interimArray
            else:
                pass

    posArr = np.delete(initalizedArray, 0, 0)

    nDel = 0
    indexDel = []
    incorrect = []
    """
    for p in range(len(posArr)):
        
        if posArr[p,0]**2 + posArr[p,1]**2 > R**2 or (HardBoundaryCircle_Disc(R,a,b,posArr[p,0],posArr[p,1],posArr[p,2])==True) :
            indexDel.append(p)
            incorrect.append(posArr[p,:])
            nDel+=1
        else:       
            pass
    
    """
    redE = np.array(incorrect)

    c = 0

    for i in indexDel:
        newArr = np.delete(posArr, i - c, 0)
        posArr = newArr
        c += 1

    print(len(posArr))
    print(len(redE))
    newLength = len(posArr)

    if n >= newLength:
        pass
    elif n < newLength:
        de = 0

        while de != newLength - n:
            k = rd.randint(0, newLength - n - de)
            newArr = np.delete(posArr, k, 0)
            posArr = newArr
            de += 1

    print(len(posArr))

    return [posArr, redE]


def init_Ann_H_Rd(n, l, r2, a, b):
    # FIX 0.2

    init_pos = np.zeros((n, 3))

    # Correct overlapping positions
    sortedInitPos = [[i, init_pos[i, 0], init_pos[i, 1], init_pos[i, 2], ] for i in range(n)]

    for x in range(n):

        valid = False
        overlapVar = False

        while valid == False:

            radius = rd.uniform(r2, l)
            angle = rd.uniform(0, 2 * np.pi)
            init_pos[x, 0] = radius * np.cos(angle)
            init_pos[x, 1] = radius * np.sin(angle)
            init_pos[x, 2] = rd.uniform(0, 2 * np.pi)

            if (init_pos[x, 0] ** 2 + init_pos[x, 1] ** 2) > ((l - b) ** 2) or (
                    init_pos[x, 0] ** 2 + init_pos[x, 1] ** 2) < ((r2 + b) ** 2):

                if (init_pos[x, 0] ** 2 + init_pos[x, 1] ** 2) > (l - a) ** 2 or (
                        HardBoundaryCircle_Disc(l - 0.2, a, b, init_pos[x, 0], init_pos[x, 1], init_pos[x, 2]) == True):
                    # print('init - border overlap')
                    continue
                elif (init_pos[x, 0] ** 2 + init_pos[x, 1] ** 2) < ((r2 + a) ** 2) or (
                        GeometricPotential(init_pos[x, 0], init_pos[x, 1], init_pos[x, 2], 0, 0, 0, b, a, r2,
                                           r2) > np.sqrt(init_pos[x, 0] ** 2 + init_pos[x, 1] ** 2)):
                    # print('init - border overlap')
                    continue
                else:
                    pass

                for j in range(n):
                    if j != x and (overlap_Ellipse(init_pos[x, 0], init_pos[x, 1], init_pos[x, 2], init_pos[j, 0],
                                                   init_pos[j, 1], init_pos[j, 2], a, b) == True):
                        overlapVar = True
                        # print('iter' + str(x)+' init - particle overlap')
                        break
                    else:
                        # print('iter' + str(x)+'init - no particle overlap')
                        overlapVar = False

                if overlapVar == True:
                    continue
                else:
                    valid = True


            else:
                for j in range(n):
                    if j != x and (overlap_Ellipse(init_pos[x, 0], init_pos[x, 1], init_pos[x, 2], init_pos[j, 0],
                                                   init_pos[j, 1], init_pos[j, 2], a, b) == True):
                        overlapVar = True
                        # print('iter' + str(x)+' init - particle overlap')
                        break
                    else:
                        # print('iter' + str(x)+'init - no particle overlap')
                        overlapVar = False

                if overlapVar == True:
                    continue
                else:
                    valid = True

    return init_pos


## -- Grid state based on Chord instead of lattice grid

def init_Circ_H_GrC(n, a, b, R, d, e):
    a_tilde = a * (1 + (d / 100))
    b_tilde = b + a * ((e / 100))

    v = int(np.floor(R / a_tilde))

    if v % 2 == 0:
        ygrid = np.linspace(-(v * a_tilde), v * a_tilde, v)

    else:
        ygrid = np.linspace(-a_tilde * (v - 1), (v - 1) * a_tilde, v)

    chordArray = []
    totalN = 0
    for k in range(v):
        chord = np.sqrt(R ** 2 - (ygrid[k]) ** 2)
        NC = int(np.floor(2 * chord / (2 * b_tilde)))
        chordArray.append([chord, NC])
        totalN += NC

    initalizedArray = np.zeros((1, 3))
    corrected = False
    countVar = 0
    for w in range(len(chordArray)):

        if chordArray[w][1] != 0:

            if chordArray[w][1] > 1:
                if chordArray[w][1] % 2 == 0:
                    right = np.array([((2 * (s + 1) - 1) * b_tilde) for s in range(int(chordArray[w][1] / 2))])
                    left = -1 * right
                    xVal = np.concatenate((right, left), axis=0)
                else:
                    right = np.array([((2 * (s + 1)) * b_tilde) for s in range(int(chordArray[w][1] / 2))])
                    left = -1 * right
                    xVal = np.concatenate((right, left), axis=0)
                    interimX = np.append(xVal, [0], axis=0)
                    xVal = interimX
            else:
                xVal = np.array([0])

            # xVal=np.linspace(-chordArray[w][0]+b_tilde,chordArray[w][0]-b_tilde,chordArray[w][1])
            for vals in xVal:
                interimArray = np.append(initalizedArray, [[vals, ygrid[w], 0]], axis=0)
                initalizedArray = interimArray
        else:
            pass

    posArr = np.delete(initalizedArray, 0, 0)

    nDel = 0
    indexDel = []
    incorrect = []
    """
    for p in range(len(posArr)):
        
        if posArr[p,0]**2 + posArr[p,1]**2 > R**2 or (HardBoundaryCircle_Disc(R,a,b,posArr[p,0],posArr[p,1],posArr[p,2])==True) :
            indexDel.append(p)
            incorrect.append(posArr[p,:])
            nDel+=1
        else:       
            pass
    
    """
    redE = np.array(incorrect)

    c = 0

    for i in indexDel:
        newArr = np.delete(posArr, i - c, 0)
        posArr = newArr
        c += 1

    print(len(posArr))
    print(len(redE))
    newLength = len(posArr)

    if n >= newLength:
        pass
    elif n < newLength:
        de = 0

        while de != newLength - n:
            k = rd.randint(0, newLength - n - de)
            newArr = np.delete(posArr, k, 0)
            posArr = newArr
            de += 1

    print(len(posArr))

    return [posArr, redE]


#@jit()
def MC_Ann_Hard(PosArray, d_pos, d_ang, steps, n, l, r2, a, b):
    # global moves
    # global accepted_moves

    moves = 0
    accepted_moves = 0
    imageName = "step_figure.png"
    numE = len(PosArray)
    kVal = b / a
    file_name = "MonteCarlo_Annulus_SimNotes.txt"

    ogdir = os.getcwd()

    main_folder_name = os.path.join(ogdir, 'annulus_R{}_r{}_n_{}_k_{}_HardBC'.format(l, r2, numE, kVal))
    folder_name = os.path.join(main_folder_name, 'instanceRun')

    if os.path.exists(main_folder_name):
        pass
    else:
        os.makedirs(main_folder_name)

    if os.path.exists(folder_name):
        expand = 0
        while True:
            expand += 1
            new_folder_name = folder_name + str(expand)
            if os.path.exists(new_folder_name):
                continue
            else:
                folder_name = new_folder_name
                break

    os.makedirs(folder_name)

    ##### Save Initial State

    fileNameArray = 'PosArray.csv'
    complete_name = os.path.join(folder_name, fileNameArray)
    np.savetxt(complete_name, PosArray, delimiter=',')

    stepsArr = np.linspace(0, steps, steps + 1)

    n = len(PosArray)

    fullImageName = os.path.join(folder_name, imageName)

    ###Save Initial State

    fileNameArray = 'PosArray.csv'
    complete_name = os.path.join(folder_name, fileNameArray)
    np.savetxt(complete_name, PosArray, delimiter=',')

    minPos = .05 * a
    minAng = .025

    fixCount = 0
    plotCount = 0

    for u in range(steps):

        for w in range(n):

            if fixCount == (np.ceil(steps / 100)):
                fixCount = 0

                if d_ang > minAng and d_pos > minPos:
                    if accepted_moves / moves < 0.47:
                        d_ang = 0.9 * d_ang
                        d_pos = 0.9 * d_ang
                    elif accepted_moves / moves < 0.37:
                        d_ang = 0.75 * d_ang
                        d_pos = 0.75 * d_ang
                    elif accepted_moves / moves < 0.27:
                        d_ang = 0.6 * d_ang
                        d_pos = 0.6 * d_ang
                    elif accepted_moves / moves < 0.17:
                        d_ang = 0.35 * d_ang
                        d_pos = 0.35 * d_ang
                    elif accepted_moves / moves < 0.07:
                        d_ang = 0.2 * d_ang
                        d_pos = 0.2 * d_ang
                    elif accepted_moves / moves > 0.57:
                        d_ang = 1.1 * d_ang
                        d_pos = 1.1 * d_ang
                    elif accepted_moves / moves > 0.67:
                        d_ang = 1.35 * d_ang
                        d_pos = 1.35 * d_ang
                    elif accepted_moves / moves > 0.77:
                        d_ang = 1.5 * d_ang
                        d_pos = 1.5 * d_ang
                    elif accepted_moves / moves > 0.87:
                        d_ang = 1.65 * d_ang
                        d_pos = 1.65 * d_ang
                    elif accepted_moves / moves > 0.97:
                        d_ang = 1.8 * d_ang
                        d_pos = 1.8 * d_ang
                else:

                    d_ang = minAng
                    d_pos = minPos

            x = d_pos * rd.uniform(-1, 1)
            y = d_pos * rd.uniform(-1, 1)
            t = d_ang * rd.uniform(-1, 1)

            # print(x,y)

            testX = PosArray[w, 0] + x
            testY = PosArray[w, 1] + y
            testT = PosArray[w, 2] + t

            if (testX ** 2 + testY ** 2) > ((l - b) ** 2) or (testX ** 2 + testY ** 2) < ((r2 + 2 * b) ** 2):
                if (testX ** 2 + testY ** 2) > ((l - b) ** 2) and (testX ** 2 + testY ** 2) < ((r2 + 2 * b) ** 2):
                    if (testX ** 2 + testY ** 2) > ((l - a) ** 2):
                        overlapVar = True
                        # print("border overlap")
                    elif (testX ** 2 + testY ** 2) < ((r2 + a) ** 2):
                        overlapVar = True
                        # print("border overlap")
                    elif HardBoundaryCircle_Disc(r2 * 1.02, a, b, testX, testY, testT) == True:
                        overlapVar = True
                        # print("border overlap")
                    elif HardBoundaryCircle_Disc(l, a, b, testX, testY, testT) == True:
                        overlapVar = True
                        # print("border overlap")
                    else:
                        for j in range(n):

                            rij = np.sqrt((testX - PosArray[j, 0]) ** 2 + (testY - PosArray[j, 1]) ** 2)

                            if rij < (2 * b):
                                if j != w and (overlap_Ellipse(testX, testY, testT, PosArray[j, 0], PosArray[j, 1],
                                                               PosArray[j, 2], a, b) == True):
                                    overlapVar = True
                                    # print("particle overlap")
                                    break
                                else:
                                    overlapVar = False
                            else:
                                overlapVar = False
                elif (testX ** 2 + testY ** 2) > ((l - b) ** 2):
                    if (testX ** 2 + testY ** 2) > ((l - a) ** 2):
                        overlapVar = True
                        # print("border overlap")
                    elif (HardBoundaryCircle_Disc(l, a, b, testX, testY, testT) == True):
                        overlapVar = True
                        # print("border overlap")
                    else:
                        for j in range(n):

                            rij = np.sqrt((testX - PosArray[j, 0]) ** 2 + (testY - PosArray[j, 1]) ** 2)

                            if rij < (2 * b):
                                if j != w and (overlap_Ellipse(testX, testY, testT, PosArray[j, 0], PosArray[j, 1],
                                                               PosArray[j, 2], a, b) == True):
                                    overlapVar = True
                                    # print("particle overlap")
                                    break
                                else:
                                    overlapVar = False
                            else:
                                overlapVar = False
                else:

                    if (testX ** 2 + testY ** 2) < ((r2 + a) ** 2):
                        overlapVar = True
                        # print("border overlap")
                    elif (HardBoundaryCircle_Disc(r2 * 1.02, a, b, testX, testY, testT) == True):
                        overlapVar = True
                        # print("border overlap")
                    else:
                        for j in range(n):

                            rij = np.sqrt((testX - PosArray[j, 0]) ** 2 + (testY - PosArray[j, 1]) ** 2)

                            if rij < (2 * b):
                                if j != w and (overlap_Ellipse(testX, testY, testT, PosArray[j, 0], PosArray[j, 1],
                                                               PosArray[j, 2], a, b) == True):
                                    overlapVar = True
                                    # print("particle overlap")
                                    break
                                else:
                                    overlapVar = False
                            else:
                                overlapVar = False
            else:
                for j in range(n):

                    rij = np.sqrt((testX - PosArray[j, 0]) ** 2 + (testY - PosArray[j, 1]) ** 2)

                    if rij < (2 * b):
                        if j != w and (
                                overlap_Ellipse(testX, testY, testT, PosArray[j, 0], PosArray[j, 1], PosArray[j, 2], a,
                                                b) == True):
                            overlapVar = True
                            # print("particle overlap")
                            break
                        else:
                            overlapVar = False
                    else:
                        overlapVar = False

            if overlapVar == True:
                pass
            elif overlapVar == False:
                accepted_moves += 1
                PosArray[w, 0] = testX
                PosArray[w, 1] = testY
                PosArray[w, 2] = testT

            moves += 1

            fixCount += 1

            x = 0
            y = 0

        plotCount += 1

        # if plotCount == (np.ceil(steps / 100)):
        #     plotCount = 0

        #     ######################################

        #     #### Periodically Save Data Snapshots #####

        #     fileNameArray = 'PosArray.csv'
        #     new_fileName = fileNameArray.split(".csv")[0] + str(u) + ".csv"
        #     complete_name = os.path.join(folder_name, new_fileName)
        #     np.savetxt(complete_name, PosArray, delimiter=',')

    reduced_density = n * np.pi * a * b / ((np.pi * l ** 2) - (np.pi * r2 ** 2))

    print(accepted_moves)
    print(moves)
    print(accepted_moves / moves)

    complete_name = os.path.join(folder_name, file_name)

    text_file = open(complete_name, 'w+')
    text_file.write("Parameters" + "\r\n")
    text_file.write("- - - - -" + "\r\n")
    text_file.write("Monte Carlo steps: " + str(steps) + "\r\n")
    text_file.write("R: " + str(l) + "\r\n")
    text_file.write("r: " + str(r2) + "\r\n")
    text_file.write("d_pos / step size: " + str(d_pos) + "\r\n")
    text_file.write("d_ang / step size: " + str(d_ang) + "\r\n")
    text_file.write("# of Ellipse: " + str(n) + "\r\n")
    text_file.write("reduced density: " + str(reduced_density) + "\r\n")
    text_file.write("Semi Minor Axis: " + str(a) + "\r\n")
    text_file.write("Semi Major Axis: " + str(b) + "\r\n")
    text_file.write("Accepted Moves: " + str(accepted_moves) + "\r\n")
    text_file.write("Total Moves: " + str(moves) + "\r\n")
    text_file.write("Acceptance Rate: " + str(100 * (accepted_moves / moves)) + " %" + "\r\n")

    text_file.close()

    fileNameArray = 'FinalPosArray.csv'
    complete_name = os.path.join(folder_name, fileNameArray)
    np.savetxt(complete_name, PosArray, delimiter=',')



R = 25
r = 0
A = 3
B = 0.5
MCSteps = 10000
stepXY = 0.5 * R
stepTh = np.pi / 2
dely = 0  # % of a
delx = 0  # % of a

N = 200

t0 = time.perf_counter()

initialState = init_Circ_H_GrC(N, A, B, R, dely, delx)[0]

MC_Circ_Hard(initialState, stepXY, stepTh, MCSteps, N, R, A, B)

print(f"Execution time: {time.perf_counter() - t0} seconds.")