import numpy as np
import os

from matplotlib.patches import Ellipse, Rectangle, Circle, Wedge
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib import cm

from joblib import Parallel, delayed
    

def PlotterAnnulus(array, a, b, L, r2, Mode, fileName):
    finalC = len(array)

    fig, ax = plt.subplots()

    dumArr1 = np.zeros((1000,))
    dumArr2 = np.full((1000,), L)

    plt.plot(dumArr1, dumArr2)

    ells = [Ellipse(xy=[array[w, 0], array[w, 1]], width=2 * b, height=2 * a, angle=array[w, 2] * 360 / (2 * np.pi),
                    color='black', fill=False, label=str(w + 1)) for w in range(finalC)]

    for count in range(finalC):
        ax.add_patch(ells[count])
        # plt.text(array[count,0],array[count,1],str(count+1))
        # print(count+1)
        # print(array[count,])

    Boundary = Circle(xy=[0, 0], radius=L, fill=False)

    ax.add_patch(Boundary)

    Boundary = Circle(xy=[0, 0], radius=r2, fill=False)

    ax.add_patch(Boundary)

    ax.set_xlim(-(1.25 * L), L + (0.25 * L))
    ax.set_ylim(-(1.25 * L), L + (0.25 * L))
    ax.set_aspect(aspect='equal')

    plt.title('N=' + str(finalC) + ' b=' + str(b) + ' k=' + str(b / a) + ' R=' + str(L) + ' r=' + str(
        r2) + ' \n' + r'$\eta =$' + str((finalC * a * b) / (L ** 2 - r2 ** 2)))

    if Mode == "D":
        plt.show()
    elif Mode == "S":
        plt.savefig(fileName, dpi=150)
        plt.clf()
        plt.close('all')


def LocalOrderParameter(array, fileName, a, b, R, r2):
    # Lines (x&y) -25,-20,-15,-10,-5,0,5,10,15,20,25

    # Sections
    # 1: -25<x=<-20 & -25<y<= -20
    # 2: -25<x=<-20 & -20<y<= -15
    # 3: -25<x=<-20 & -15<y<= -10
    # 4: -25<x=<-20 & -10<y<= -5
    # 5: -25<x=<-20 & -5<y<= 0
    # 6: -25<x=<-20 & 0<y<= 5
    # 7: -25<x=<-20 & 5<y<= 10
    # 8: -25<x=<-20 & 10<y<= 15
    # 9: -25<x=<-20 & 15<y<= 20
    # 10: -25<x=<-20 & 20<y<= 25

    orderParameter = np.zeros((100, 3))

    # Entry 1 is sin(2*theta), Entry 2 is cos(2*theta), Entry 3 is N of Particles in region

    for element in array:

        if element[0] > -25 and element[0] <= -20:

            if element[1] > -25 and element[1] <= -20:

                orderParameter[0, 0] += np.sin(2 * element[2])
                orderParameter[0, 1] += np.cos(2 * element[2])
                orderParameter[0, 2] += 1

            elif element[1] > -20 and element[1] <= -15:

                orderParameter[1, 0] += np.sin(2 * element[2])
                orderParameter[1, 1] += np.cos(2 * element[2])
                orderParameter[1, 2] += 1

            elif element[1] > -15 and element[1] <= -10:

                orderParameter[2, 0] += np.sin(2 * element[2])
                orderParameter[2, 1] += np.cos(2 * element[2])
                orderParameter[2, 2] += 1

            elif element[1] > -10 and element[1] <= -5:

                orderParameter[3, 0] += np.sin(2 * element[2])
                orderParameter[3, 1] += np.cos(2 * element[2])
                orderParameter[3, 2] += 1

            elif element[1] > -5 and element[1] <= 0:

                orderParameter[4, 0] += np.sin(2 * element[2])
                orderParameter[4, 1] += np.cos(2 * element[2])
                orderParameter[4, 2] += 1

            elif element[1] > 0 and element[1] <= 5:

                orderParameter[5, 0] += np.sin(2 * element[2])
                orderParameter[5, 1] += np.cos(2 * element[2])
                orderParameter[5, 2] += 1

            elif element[1] > 5 and element[1] <= 10:

                orderParameter[6, 0] += np.sin(2 * element[2])
                orderParameter[6, 1] += np.cos(2 * element[2])
                orderParameter[6, 2] += 1

            elif element[1] > 10 and element[1] <= 15:

                orderParameter[7, 0] += np.sin(2 * element[2])
                orderParameter[7, 1] += np.cos(2 * element[2])
                orderParameter[7, 2] += 1

            elif element[1] > 15 and element[1] <= 20:

                orderParameter[8, 0] += np.sin(2 * element[2])
                orderParameter[8, 1] += np.cos(2 * element[2])
                orderParameter[8, 2] += 1

            else:

                orderParameter[9, 0] += np.sin(2 * element[2])
                orderParameter[9, 1] += np.cos(2 * element[2])
                orderParameter[9, 2] += 1

        elif element[0] > -20 and element[0] <= -15:

            if element[1] > -25 and element[1] <= -20:

                orderParameter[10, 0] += np.sin(2 * element[2])
                orderParameter[10, 1] += np.cos(2 * element[2])
                orderParameter[10, 2] += 1

            elif element[1] > -20 and element[1] <= -15:

                orderParameter[11, 0] += np.sin(2 * element[2])
                orderParameter[11, 1] += np.cos(2 * element[2])
                orderParameter[11, 2] += 1

            elif element[1] > -15 and element[1] <= -10:

                orderParameter[12, 0] += np.sin(2 * element[2])
                orderParameter[12, 1] += np.cos(2 * element[2])
                orderParameter[12, 2] += 1

            elif element[1] > -10 and element[1] <= -5:

                orderParameter[13, 0] += np.sin(2 * element[2])
                orderParameter[13, 1] += np.cos(2 * element[2])
                orderParameter[13, 2] += 1

            elif element[1] > -5 and element[1] <= 0:

                orderParameter[14, 0] += np.sin(2 * element[2])
                orderParameter[14, 1] += np.cos(2 * element[2])
                orderParameter[14, 2] += 1

            elif element[1] > 0 and element[1] <= 5:

                orderParameter[15, 0] += np.sin(2 * element[2])
                orderParameter[15, 1] += np.cos(2 * element[2])
                orderParameter[15, 2] += 1

            elif element[1] > 5 and element[1] <= 10:

                orderParameter[16, 0] += np.sin(2 * element[2])
                orderParameter[16, 1] += np.cos(2 * element[2])
                orderParameter[16, 2] += 1

            elif element[1] > 10 and element[1] <= 15:

                orderParameter[17, 0] += np.sin(2 * element[2])
                orderParameter[17, 1] += np.cos(2 * element[2])
                orderParameter[17, 2] += 1

            elif element[1] > 15 and element[1] <= 20:

                orderParameter[18, 0] += np.sin(2 * element[2])
                orderParameter[18, 1] += np.cos(2 * element[2])
                orderParameter[18, 2] += 1

            else:

                orderParameter[19, 0] += np.sin(2 * element[2])
                orderParameter[19, 1] += np.cos(2 * element[2])
                orderParameter[19, 2] += 1

        elif element[0] > -15 and element[0] <= -10:

            if element[1] > -25 and element[1] <= -20:

                orderParameter[20, 0] += np.sin(2 * element[2])
                orderParameter[20, 1] += np.cos(2 * element[2])
                orderParameter[20, 2] += 1

            elif element[1] > -20 and element[1] <= -15:

                orderParameter[21, 0] += np.sin(2 * element[2])
                orderParameter[21, 1] += np.cos(2 * element[2])
                orderParameter[21, 2] += 1

            elif element[1] > -15 and element[1] <= -10:

                orderParameter[22, 0] += np.sin(2 * element[2])
                orderParameter[22, 1] += np.cos(2 * element[2])
                orderParameter[22, 2] += 1

            elif element[1] > -10 and element[1] <= -5:

                orderParameter[23, 0] += np.sin(2 * element[2])
                orderParameter[23, 1] += np.cos(2 * element[2])
                orderParameter[23, 2] += 1

            elif element[1] > -5 and element[1] <= 0:

                orderParameter[24, 0] += np.sin(2 * element[2])
                orderParameter[24, 1] += np.cos(2 * element[2])
                orderParameter[24, 2] += 1

            elif element[1] > 0 and element[1] <= 5:

                orderParameter[25, 0] += np.sin(2 * element[2])
                orderParameter[25, 1] += np.cos(2 * element[2])
                orderParameter[25, 2] += 1

            elif element[1] > 5 and element[1] <= 10:

                orderParameter[26, 0] += np.sin(2 * element[2])
                orderParameter[26, 1] += np.cos(2 * element[2])
                orderParameter[26, 2] += 1

            elif element[1] > 10 and element[1] <= 15:

                orderParameter[27, 0] += np.sin(2 * element[2])
                orderParameter[27, 1] += np.cos(2 * element[2])
                orderParameter[27, 2] += 1

            elif element[1] > 15 and element[1] <= 20:

                orderParameter[28, 0] += np.sin(2 * element[2])
                orderParameter[28, 1] += np.cos(2 * element[2])
                orderParameter[28, 2] += 1

            else:

                orderParameter[29, 0] += np.sin(2 * element[2])
                orderParameter[29, 1] += np.cos(2 * element[2])
                orderParameter[29, 2] += 1

        elif element[0] > -10 and element[0] <= -5:

            if element[1] > -25 and element[1] <= -20:

                orderParameter[30, 0] += np.sin(2 * element[2])
                orderParameter[30, 1] += np.cos(2 * element[2])
                orderParameter[30, 2] += 1

            elif element[1] > -20 and element[1] <= -15:

                orderParameter[31, 0] += np.sin(2 * element[2])
                orderParameter[31, 1] += np.cos(2 * element[2])
                orderParameter[31, 2] += 1

            elif element[1] > -15 and element[1] <= -10:

                orderParameter[32, 0] += np.sin(2 * element[2])
                orderParameter[32, 1] += np.cos(2 * element[2])
                orderParameter[32, 2] += 1

            elif element[1] > -10 and element[1] <= -5:

                orderParameter[33, 0] += np.sin(2 * element[2])
                orderParameter[33, 1] += np.cos(2 * element[2])
                orderParameter[33, 2] += 1

            elif element[1] > -5 and element[1] <= 0:

                orderParameter[34, 0] += np.sin(2 * element[2])
                orderParameter[34, 1] += np.cos(2 * element[2])
                orderParameter[34, 2] += 1

            elif element[1] > 0 and element[1] <= 5:

                orderParameter[35, 0] += np.sin(2 * element[2])
                orderParameter[35, 1] += np.cos(2 * element[2])
                orderParameter[35, 2] += 1

            elif element[1] > 5 and element[1] <= 10:

                orderParameter[36, 0] += np.sin(2 * element[2])
                orderParameter[36, 1] += np.cos(2 * element[2])
                orderParameter[36, 2] += 1

            elif element[1] > 10 and element[1] <= 15:

                orderParameter[37, 0] += np.sin(2 * element[2])
                orderParameter[37, 1] += np.cos(2 * element[2])
                orderParameter[37, 2] += 1

            elif element[1] > 15 and element[1] <= 20:

                orderParameter[38, 0] += np.sin(2 * element[2])
                orderParameter[38, 1] += np.cos(2 * element[2])
                orderParameter[38, 2] += 1

            else:

                orderParameter[39, 0] += np.sin(2 * element[2])
                orderParameter[39, 1] += np.cos(2 * element[2])
                orderParameter[39, 2] += 1

        elif element[0] > -5 and element[0] <= 0:

            if element[1] > -25 and element[1] <= -20:

                orderParameter[40, 0] += np.sin(2 * element[2])
                orderParameter[40, 1] += np.cos(2 * element[2])
                orderParameter[40, 2] += 1

            elif element[1] > -20 and element[1] <= -15:

                orderParameter[41, 0] += np.sin(2 * element[2])
                orderParameter[41, 1] += np.cos(2 * element[2])
                orderParameter[41, 2] += 1

            elif element[1] > -15 and element[1] <= -10:

                orderParameter[42, 0] += np.sin(2 * element[2])
                orderParameter[42, 1] += np.cos(2 * element[2])
                orderParameter[42, 2] += 1

            elif element[1] > -10 and element[1] <= -5:

                orderParameter[43, 0] += np.sin(2 * element[2])
                orderParameter[43, 1] += np.cos(2 * element[2])
                orderParameter[43, 2] += 1

            elif element[1] > -5 and element[1] <= 0:

                orderParameter[44, 0] += np.sin(2 * element[2])
                orderParameter[44, 1] += np.cos(2 * element[2])
                orderParameter[44, 2] += 1

            elif element[1] > 0 and element[1] <= 5:

                orderParameter[45, 0] += np.sin(2 * element[2])
                orderParameter[45, 1] += np.cos(2 * element[2])
                orderParameter[45, 2] += 1

            elif element[1] > 5 and element[1] <= 10:

                orderParameter[46, 0] += np.sin(2 * element[2])
                orderParameter[46, 1] += np.cos(2 * element[2])
                orderParameter[46, 2] += 1

            elif element[1] > 10 and element[1] <= 15:

                orderParameter[47, 0] += np.sin(2 * element[2])
                orderParameter[47, 1] += np.cos(2 * element[2])
                orderParameter[47, 2] += 1

            elif element[1] > 15 and element[1] <= 20:

                orderParameter[48, 0] += np.sin(2 * element[2])
                orderParameter[48, 1] += np.cos(2 * element[2])
                orderParameter[48, 2] += 1

            else:

                orderParameter[49, 0] += np.sin(2 * element[2])
                orderParameter[49, 1] += np.cos(2 * element[2])
                orderParameter[49, 2] += 1

        elif element[0] > 0 and element[0] <= 5:

            if element[1] > -25 and element[1] <= -20:

                orderParameter[50, 0] += np.sin(2 * element[2])
                orderParameter[50, 1] += np.cos(2 * element[2])
                orderParameter[50, 2] += 1

            elif element[1] > -20 and element[1] <= -15:

                orderParameter[51, 0] += np.sin(2 * element[2])
                orderParameter[51, 1] += np.cos(2 * element[2])
                orderParameter[51, 2] += 1

            elif element[1] > -15 and element[1] <= -10:

                orderParameter[52, 0] += np.sin(2 * element[2])
                orderParameter[52, 1] += np.cos(2 * element[2])
                orderParameter[52, 2] += 1

            elif element[1] > -10 and element[1] <= -5:

                orderParameter[53, 0] += np.sin(2 * element[2])
                orderParameter[53, 1] += np.cos(2 * element[2])
                orderParameter[53, 2] += 1

            elif element[1] > -5 and element[1] <= 0:

                orderParameter[54, 0] += np.sin(2 * element[2])
                orderParameter[54, 1] += np.cos(2 * element[2])
                orderParameter[54, 2] += 1

            elif element[1] > 0 and element[1] <= 5:

                orderParameter[55, 0] += np.sin(2 * element[2])
                orderParameter[55, 1] += np.cos(2 * element[2])
                orderParameter[55, 2] += 1

            elif element[1] > 5 and element[1] <= 10:

                orderParameter[56, 0] += np.sin(2 * element[2])
                orderParameter[56, 1] += np.cos(2 * element[2])
                orderParameter[56, 2] += 1

            elif element[1] > 10 and element[1] <= 15:

                orderParameter[57, 0] += np.sin(2 * element[2])
                orderParameter[57, 1] += np.cos(2 * element[2])
                orderParameter[57, 2] += 1

            elif element[1] > 15 and element[1] <= 20:

                orderParameter[58, 0] += np.sin(2 * element[2])
                orderParameter[58, 1] += np.cos(2 * element[2])
                orderParameter[58, 2] += 1

            else:

                orderParameter[59, 0] += np.sin(2 * element[2])
                orderParameter[59, 1] += np.cos(2 * element[2])
                orderParameter[59, 2] += 1

        elif element[0] > 5 and element[0] <= 10:

            if element[1] > -25 and element[1] <= -20:

                orderParameter[60, 0] += np.sin(2 * element[2])
                orderParameter[60, 1] += np.cos(2 * element[2])
                orderParameter[60, 2] += 1

            elif element[1] > -20 and element[1] <= -15:

                orderParameter[61, 0] += np.sin(2 * element[2])
                orderParameter[61, 1] += np.cos(2 * element[2])
                orderParameter[61, 2] += 1

            elif element[1] > -15 and element[1] <= -10:

                orderParameter[62, 0] += np.sin(2 * element[2])
                orderParameter[62, 1] += np.cos(2 * element[2])
                orderParameter[62, 2] += 1

            elif element[1] > -10 and element[1] <= -5:

                orderParameter[63, 0] += np.sin(2 * element[2])
                orderParameter[63, 1] += np.cos(2 * element[2])
                orderParameter[63, 2] += 1

            elif element[1] > -5 and element[1] <= 0:

                orderParameter[64, 0] += np.sin(2 * element[2])
                orderParameter[64, 1] += np.cos(2 * element[2])
                orderParameter[64, 2] += 1

            elif element[1] > 0 and element[1] <= 5:

                orderParameter[65, 0] += np.sin(2 * element[2])
                orderParameter[65, 1] += np.cos(2 * element[2])
                orderParameter[65, 2] += 1

            elif element[1] > 5 and element[1] <= 10:

                orderParameter[66, 0] += np.sin(2 * element[2])
                orderParameter[66, 1] += np.cos(2 * element[2])
                orderParameter[66, 2] += 1

            elif element[1] > 10 and element[1] <= 15:

                orderParameter[67, 0] += np.sin(2 * element[2])
                orderParameter[67, 1] += np.cos(2 * element[2])
                orderParameter[67, 2] += 1

            elif element[1] > 15 and element[1] <= 20:

                orderParameter[68, 0] += np.sin(2 * element[2])
                orderParameter[68, 1] += np.cos(2 * element[2])
                orderParameter[68, 2] += 1

            else:

                orderParameter[69, 0] += np.sin(2 * element[2])
                orderParameter[69, 1] += np.cos(2 * element[2])
                orderParameter[69, 2] += 1

        elif element[0] > 10 and element[0] <= 15:

            if element[1] > -25 and element[1] <= -20:

                orderParameter[70, 0] += np.sin(2 * element[2])
                orderParameter[70, 1] += np.cos(2 * element[2])
                orderParameter[70, 2] += 1

            elif element[1] > -20 and element[1] <= -15:

                orderParameter[71, 0] += np.sin(2 * element[2])
                orderParameter[71, 1] += np.cos(2 * element[2])
                orderParameter[71, 2] += 1

            elif element[1] > -15 and element[1] <= -10:

                orderParameter[72, 0] += np.sin(2 * element[2])
                orderParameter[72, 1] += np.cos(2 * element[2])
                orderParameter[72, 2] += 1

            elif element[1] > -10 and element[1] <= -5:

                orderParameter[73, 0] += np.sin(2 * element[2])
                orderParameter[73, 1] += np.cos(2 * element[2])
                orderParameter[73, 2] += 1

            elif element[1] > -5 and element[1] <= 0:

                orderParameter[74, 0] += np.sin(2 * element[2])
                orderParameter[74, 1] += np.cos(2 * element[2])
                orderParameter[74, 2] += 1

            elif element[1] > 0 and element[1] <= 5:

                orderParameter[75, 0] += np.sin(2 * element[2])
                orderParameter[75, 1] += np.cos(2 * element[2])
                orderParameter[75, 2] += 1

            elif element[1] > 5 and element[1] <= 10:

                orderParameter[76, 0] += np.sin(2 * element[2])
                orderParameter[76, 1] += np.cos(2 * element[2])
                orderParameter[76, 2] += 1

            elif element[1] > 10 and element[1] <= 15:

                orderParameter[77, 0] += np.sin(2 * element[2])
                orderParameter[77, 1] += np.cos(2 * element[2])
                orderParameter[77, 2] += 1

            elif element[1] > 15 and element[1] <= 20:

                orderParameter[78, 0] += np.sin(2 * element[2])
                orderParameter[78, 1] += np.cos(2 * element[2])
                orderParameter[78, 2] += 1

            else:

                orderParameter[79, 0] += np.sin(2 * element[2])
                orderParameter[79, 1] += np.cos(2 * element[2])
                orderParameter[79, 2] += 1

        elif element[0] > 15 and element[0] <= 20:

            if element[1] > -25 and element[1] <= -20:

                orderParameter[80, 0] += np.sin(2 * element[2])
                orderParameter[80, 1] += np.cos(2 * element[2])
                orderParameter[80, 2] += 1

            elif element[1] > -20 and element[1] <= -15:

                orderParameter[81, 0] += np.sin(2 * element[2])
                orderParameter[81, 1] += np.cos(2 * element[2])
                orderParameter[81, 2] += 1

            elif element[1] > -15 and element[1] <= -10:

                orderParameter[82, 0] += np.sin(2 * element[2])
                orderParameter[82, 1] += np.cos(2 * element[2])
                orderParameter[82, 2] += 1

            elif element[1] > -10 and element[1] <= -5:

                orderParameter[83, 0] += np.sin(2 * element[2])
                orderParameter[83, 1] += np.cos(2 * element[2])
                orderParameter[83, 2] += 1

            elif element[1] > -5 and element[1] <= 0:

                orderParameter[84, 0] += np.sin(2 * element[2])
                orderParameter[84, 1] += np.cos(2 * element[2])
                orderParameter[84, 2] += 1

            elif element[1] > 0 and element[1] <= 5:

                orderParameter[85, 0] += np.sin(2 * element[2])
                orderParameter[85, 1] += np.cos(2 * element[2])
                orderParameter[85, 2] += 1

            elif element[1] > 5 and element[1] <= 10:

                orderParameter[86, 0] += np.sin(2 * element[2])
                orderParameter[86, 1] += np.cos(2 * element[2])
                orderParameter[86, 2] += 1

            elif element[1] > 10 and element[1] <= 15:

                orderParameter[87, 0] += np.sin(2 * element[2])
                orderParameter[87, 1] += np.cos(2 * element[2])
                orderParameter[87, 2] += 1

            elif element[1] > 15 and element[1] <= 20:

                orderParameter[88, 0] += np.sin(2 * element[2])
                orderParameter[88, 1] += np.cos(2 * element[2])
                orderParameter[88, 2] += 1

            else:

                orderParameter[89, 0] += np.sin(2 * element[2])
                orderParameter[89, 1] += np.cos(2 * element[2])
                orderParameter[89, 2] += 1

        else:

            if element[1] > -25 and element[1] <= -20:

                orderParameter[90, 0] += np.sin(2 * element[2])
                orderParameter[90, 1] += np.cos(2 * element[2])
                orderParameter[90, 2] += 1

            elif element[1] > -20 and element[1] <= -15:

                orderParameter[91, 0] += np.sin(2 * element[2])
                orderParameter[91, 1] += np.cos(2 * element[2])
                orderParameter[91, 2] += 1

            elif element[1] > -15 and element[1] <= -10:

                orderParameter[92, 0] += np.sin(2 * element[2])
                orderParameter[92, 1] += np.cos(2 * element[2])
                orderParameter[92, 2] += 1

            elif element[1] > -10 and element[1] <= -5:

                orderParameter[93, 0] += np.sin(2 * element[2])
                orderParameter[93, 1] += np.cos(2 * element[2])
                orderParameter[93, 2] += 1

            elif element[1] > -5 and element[1] <= 0:

                orderParameter[94, 0] += np.sin(2 * element[2])
                orderParameter[94, 1] += np.cos(2 * element[2])
                orderParameter[94, 2] += 1

            elif element[1] > 0 and element[1] <= 5:

                orderParameter[95, 0] += np.sin(2 * element[2])
                orderParameter[95, 1] += np.cos(2 * element[2])
                orderParameter[95, 2] += 1

            elif element[1] > 5 and element[1] <= 10:

                orderParameter[96, 0] += np.sin(2 * element[2])
                orderParameter[96, 1] += np.cos(2 * element[2])
                orderParameter[96, 2] += 1

            elif element[1] > 10 and element[1] <= 15:

                orderParameter[97, 0] += np.sin(2 * element[2])
                orderParameter[97, 1] += np.cos(2 * element[2])
                orderParameter[97, 2] += 1

            elif element[1] > 15 and element[1] <= 20:

                orderParameter[98, 0] += np.sin(2 * element[2])
                orderParameter[98, 1] += np.cos(2 * element[2])
                orderParameter[98, 2] += 1

            else:

                orderParameter[99, 0] += np.sin(2 * element[2])
                orderParameter[99, 1] += np.cos(2 * element[2])
                orderParameter[99, 2] += 1

    bounds = [-25, -20, -15, -10, -5, 0, 5, 10, 15, 20]

    patches = [Rectangle([i, j], 5, 5, angle=0.0) for i in bounds for j in bounds]

    lambdaOP = np.zeros((100,))

    for k in range(100):

        if orderParameter[k, 2] != 0:
            lambdaOP[k] = np.sqrt(
                (orderParameter[k, 0] / orderParameter[k, 2]) ** 2 + (orderParameter[k, 1] / orderParameter[k, 2]) ** 2)
        else:
            lambdaOP[k] = 0

    bulkLambda = np.sqrt(
        (np.sum(orderParameter[:, 0]) / (len(array))) ** 2 + (np.sum(orderParameter[:, 1]) / (len(array))) ** 2)

    finalC = len(array)

    fig, ax = plt.subplots()

    L = R

    ells = [Ellipse(xy=[array[w, 0], array[w, 1]], width=2 * b, height=2 * a, angle=array[w, 2] * 360 / (2 * np.pi),
                    color='black', fill=False) for w in range(finalC)]

    for count in range(finalC):
        ax.add_patch(ells[count])
        # plt.text(array[count,0],array[count,1],str(count+1))
        print(count + 1)
        print(array[count,])

    Boundary = Circle(xy=[0, 0], radius=L, fill=False)

    ax.add_patch(Boundary)

    Boundary = Circle(xy=[0, 0], radius=r2, fill=False)

    ax.add_patch(Boundary)

    p = PatchCollection(patches)
    p.set(array=lambdaOP, cmap='Reds', alpha=0.7)
    ax.add_collection(p)
    ax.set_xlim(-30, 30)
    ax.set_ylim(-30, 30)
    ax.axis('equal')
    fig.colorbar(p)

    # textstr = 'NSE=%.2f\nRMSE=%.2f\n'%(1, 2)

    textstr = 'Bulk Order Parameter ' + r'$\overline{\Lambda}=$' + str(bulkLambda)

    plt.text(0, -35, textstr, ha='center', fontsize=8)
    plt.title('Local Order Parameter ' + r'$\Lambda$' + ' Contour Map', fontsize=10)
    plt.subplots_adjust(bottom=0.15)
    plt.savefig(fileName, dpi=300)
    plt.clf()
    plt.close('all')


def PlotFile(location, R, r2, a, b):
    files = os.listdir(location)

    PosFiles = []

    for fileNames in files:

        if fileNames.startswith('PosArray'):
            PosFiles.append(fileNames)
        else:
            pass

    imagesFolder = os.path.join(location, 'Plots')

    if os.path.exists(imagesFolder):
        pass
    else:
        os.makedirs(imagesFolder)

    imagesFolder2 = os.path.join(location, 'Lambda Plots')

    if os.path.exists(imagesFolder2):
        pass
    else:
        os.makedirs(imagesFolder2)

    def PlotIter(currentFile, iF1, iF2):

        Loc = os.path.join(location, currentFile)
        PosArray = np.loadtxt(Loc, dtype=float, delimiter=',')
        imageName = os.path.join(iF1, 'step_figure' + currentFile[8:len(currentFile) - 4] + '.png')
        PlotterAnnulus(PosArray, a, b, R, r2, 'S', imageName)

        imageName = os.path.join(iF2, 'lambda_orderParameter' + currentFile[8:len(currentFile) - 4] + '.png')
        LocalOrderParameter(PosArray, imageName, a, b, R, r2)

    results = Parallel(n_jobs=2)(delayed(PlotIter)(p, imagesFolder, imagesFolder2) for p in PosFiles)

    print('Images plotted #########')


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def PolarSystem(n_r, n_theta, a, b, posArray, R, r, fileName):
    #### Transform to Polar coords

    shapePos = np.shape(posArray)
    newPosArray = np.zeros(shapePos)

    for k in range(len(posArray)):
        newPosArray[k, 0] = np.sqrt(posArray[k, 0] ** 2 + posArray[k, 1] ** 2)

        newPosArray[k, 1] = np.arctan2(posArray[k, 1], posArray[k, 0])

        newPosArray[k, 2] = posArray[k, 2]

    segments = np.zeros(shape=(n_r, n_theta, 7), dtype="complex_")

    deltaR = R / n_r
    deltaTheta = 2 * np.pi / n_theta

    r_discretized = [(h + 1) * deltaR for h in range(n_r)]
    theta_discretized = [(h + 1) * deltaTheta for h in range(n_theta)]

    # print(r_discretized)
    # print(theta_discretized)

    for i in range(n_theta):
        for c in range(n_r):
            segments[c, i, 0] = r_discretized[c]
            segments[c, i, 1] = theta_discretized[i]

    ## SORTING SEGMENTS IS ONLY THING NOT WORKING PROPERLY

    for p in newPosArray:

        r_i = p[0]

        if p[1] > 2 * np.pi:
            theta_i = p[1] % (2 * np.pi)
        elif p[1] < 0 and abs(p[1]) > 2 * np.pi:
            theta_i = ((p[1]) % (2 * np.pi)) + 2 * np.pi
        elif p[1] < 0 and abs(p[1]) <= 2 * np.pi:
            theta_i = (p[1]) + 2 * np.pi
        else:
            theta_i = p[1]

        """
        if r_i <= r_discretized[0]:
            segment_r=0
        elif r_i <= r_discretized[n_r-1] and r_i > r_discretized[n_r-2]:
            segment_r=n_r-1
        else:
            for f in range(len(r_discretized[1:n_r-1])):

                if r_i > r_discretized[f+1] and r_i <= r_discretized[f+2]:
                    segment_r= f+1
                else:
                    pass

        if theta_i > theta_discretized[n_theta-1] and theta_i <= theta_discretized[0]:
            segment_theta=0
        elif    
        """

        trialIndex = find_nearest(r_discretized, r_i)

        if trialIndex == n_r - 1:
            segment_r = trialIndex
        elif trialIndex == 0:
            if r_i > r_discretized[trialIndex]:
                segment_r = 1
            else:
                segment_r = 0
        else:
            if r_i > r_discretized[trialIndex]:
                segment_r = trialIndex + 1
            else:
                segment_r = trialIndex

        trialIndex = find_nearest(theta_discretized, theta_i)

        if trialIndex == n_theta - 1:
            if theta_i > theta_discretized[trialIndex]:
                segment_theta = 0
            else:
                segment_theta = trialIndex

        elif trialIndex == 0:
            if theta_i > theta_discretized[trialIndex]:
                segment_theta = 1
            else:
                segment_theta = 0

        else:
            if theta_i > theta_discretized[trialIndex]:
                segment_theta = trialIndex + 1
            else:
                segment_theta = trialIndex

        segments[segment_r, segment_theta, 2] += 1

        ## Scalar Orientational Order Parameter Lambda summation
        ## Index 3 is Sin2theta Index 4 is Cos2theta

        segments[segment_r, segment_theta, 3] += np.sin(2 * p[2])
        segments[segment_r, segment_theta, 4] += np.cos(2 * p[2])

        ## Positional Order Parameter for x gamma_x summation
        k_x = 2 * np.pi / R
        segments[segment_r, segment_theta, 5] += np.exp(1j * (p[0] * np.cos(p[1]) * k_x))

        ## Positional Order Parameter for y gamma_y summation
        k_y = 2 * np.pi / R
        segments[segment_r, segment_theta, 6] += np.exp(1j * (p[0] * np.sin(p[1]) * k_y))

    orderParameters = np.zeros((n_r * n_theta, 3))
    ## 0 index is lambda , 1 index is gamma x, 2 index is gamma y

    print(segments[:, :, 2])

    print(segments[:, :, 3])

    print(segments[:, :, 4])

    print(segments[:, :, 5])

    print(segments[:, :, 6])

    i = 0
    c = 0

    counterVar = 0

    for i in range(n_theta):
        for c in range(n_r):
            if segments[c, i, 2] != 0:
                orderParameters[counterVar, 0] = np.sqrt(np.real(
                    (segments[c, i, 3] / segments[c, i, 2]) ** 2 + (segments[c, i, 4] / segments[c, i, 2]) ** 2))

                orderParameters[counterVar, 1] = np.abs(segments[c, i, 5] / segments[c, i, 2])

                orderParameters[counterVar, 2] = np.abs(segments[c, i, 6] / segments[c, i, 2])
            else:
                orderParameters[counterVar, 0] = 0

                orderParameters[counterVar, 1] = 0

                orderParameters[counterVar, 2] = 0

            counterVar += 1

    print(orderParameters[:, 0])

    bulkLambda = np.sqrt(
        np.real((np.sum(segments[:, :, 3]) / len(posArray)) ** 2 + (np.sum(segments[:, :, 4]) / len(posArray)) ** 2))
    bulkGammaX = np.abs(np.sum(segments[:, :, 5]) / len(posArray))
    bulkGammaY = np.abs(np.sum(segments[:, :, 6]) / len(posArray))

    WedgeArray = [Wedge((0, 0), (h + 1) * (((R - r)) / n_r) + r, theta1=(k + 1) * (1 / n_theta) * 360,
                        theta2=(k + 2) * (1 / n_theta) * 360, width=(R - r) / n_r) for k in range(n_theta) for h in
                  range(n_r)]

    fig = plt.figure()
    ax = fig.add_subplot(aspect='equal')

    ### Lambda OP Plot

    array = posArray

    ells = [Ellipse(xy=[array[w, 0], array[w, 1]], width=2 * b, height=2 * a, angle=array[w, 2] * 360 / (2 * np.pi),
                    color='black', fill=False) for w in range(len(posArray))]

    p = PatchCollection(WedgeArray)
    p.set(array=orderParameters[:, 0], cmap='Reds', alpha=0.7)
    ax.add_collection(p)

    ax.axis('equal')
    fig.colorbar(p)

    # textstr = 'NSE=%.2f\nRMSE=%.2f\n'%(1, 2)

    textstr = 'Bulk Order Parameter ' + r'$\overline{\Lambda}=$' + str(bulkLambda)

    imageName = fileName + '_lambda.png'

    plt.text(0, -35, textstr, ha='center', fontsize=8)
    plt.title('Local Order Parameter ' + r'$\Lambda$' + ' Contour Map', fontsize=10)
    plt.subplots_adjust(bottom=0.15)
    plt.savefig(imageName, dpi=150)
    plt.clf()
    plt.close('all')

    fig = plt.figure()
    ax = fig.add_subplot(aspect='equal')

    ### Gamma X OP Plot

    p = PatchCollection(WedgeArray)
    p.set(array=orderParameters[:, 1], cmap='Blues', alpha=0.7)
    ax.add_collection(p)

    ax.axis('equal')
    fig.colorbar(p)

    # textstr = 'NSE=%.2f\nRMSE=%.2f\n'%(1, 2)

    textstr = 'Bulk Order Parameter ' + r'$\overline{\gamma_X}=$' + str(bulkGammaX)

    imageName = fileName + '_gammaX.png'

    plt.text(0, -35, textstr, ha='center', fontsize=8)
    plt.title('Local Order Parameter ' + r'$\gamma_X$' + ' Contour Map', fontsize=10)
    plt.subplots_adjust(bottom=0.15)
    plt.savefig(imageName, dpi=150)
    plt.clf()
    plt.close('all')

    fig = plt.figure()
    ax = fig.add_subplot(aspect='equal')

    ### Gamma X OP Plot

    p = PatchCollection(WedgeArray)
    p.set(array=orderParameters[:, 2], cmap='Greens', alpha=0.7)
    ax.add_collection(p)

    ax.axis('equal')
    fig.colorbar(p)

    # textstr = 'NSE=%.2f\nRMSE=%.2f\n'%(1, 2)

    textstr = 'Bulk Order Parameter ' + r'$\overline{\gamma_Y}=$' + str(bulkGammaY)

    imageName = fileName + '_gammaY.png'

    plt.text(0, -35, textstr, ha='center', fontsize=8)
    plt.title('Local Order Parameter ' + r'$\gamma_Y$' + ' Contour Map', fontsize=10)
    plt.subplots_adjust(bottom=0.15)
    plt.savefig(imageName, dpi=150)
    plt.clf()
    plt.close('all')

    # eta = len(posArray)*(a*b)/(R**2 - r**2)

    return [bulkLambda, bulkGammaX, bulkGammaY]


# for CountN in [90,100,120]:


def NewPlotFile(location, R, r2, a, b, rDiv, angDiv, OP):
    files = os.listdir(location)

    PosFiles = []

    for fileNames in files:

        if fileNames.startswith('PosArray'):
            PosFiles.append(fileNames)
        else:
            pass

    imagesFolder = os.path.join(location, 'Plots')

    if os.path.exists(imagesFolder):
        pass
    else:
        os.makedirs(imagesFolder)

    imagesFolder2 = os.path.join(location, 'Order Parameter Plots')

    if os.path.exists(imagesFolder2):
        pass
    else:
        os.makedirs(imagesFolder2)

    # def PlotIter(currentFile,iF1,iF2):

    iF1 = imagesFolder
    iF2 = imagesFolder2

    initLoc = os.path.join(location, PosFiles[0])
    initialArray = np.loadtxt(initLoc, dtype=float, delimiter=',')

    eta = len(initialArray) * (a * b) / (R ** 2 - r2 ** 2)

    totalResults = []

    for currentFile in PosFiles:

        Loc = os.path.join(location, currentFile)
        PosArray = np.loadtxt(Loc, dtype=float, delimiter=',')
        imageName = os.path.join(iF1, 'step_figure' + currentFile[8:len(currentFile) - 4] + '.png')
        PlotterAnnulus(PosArray, a, b, R, r2, 'S', imageName)

        imageName = os.path.join(iF2, currentFile[8:len(currentFile) - 4])
        # LocalOrderParameter(PosArray,imageName,a,b,R,r2)

        if OP == True:
            res = PolarSystem(rDiv, angDiv, a, b, PosArray, R, r2, imageName)

            totalResults.append(res)
        else:
            pass

    stepArray = np.linspace(0, len(PosFiles), len(PosFiles))
    if OP == True:

        fig = plt.figure()
        ax = fig.add_subplot(aspect='equal')
        plt.plot(stepArray, np.array(totalResults)[:, 0])
        plt.title('Bulk Order Parameter ' + r'$\overline{\Lambda}$' + ' Over Whole Run', fontsize=10)
        imageName = os.path.join(location, 'bulkLambdaRun.png')
        plt.savefig(imageName, dpi=150)
        plt.clf()
        plt.close('all')

        fig = plt.figure()
        ax = fig.add_subplot(aspect='equal')
        plt.plot(stepArray, np.array(totalResults)[:, 1])
        plt.title('Bulk Order Parameter ' + r'$\overline{\gamma_X}$' + ' Over Whole Run', fontsize=10)
        imageName = os.path.join(location, 'bulkGammaXRun.png')
        plt.savefig(imageName, dpi=150)
        plt.clf()
        plt.close('all')

        fig = plt.figure()
        ax = fig.add_subplot(aspect='equal')
        plt.plot(stepArray, np.array(totalResults)[:, 2])
        plt.title('Bulk Order Parameter ' + r'$\overline{\gamma_Y}$' + ' Over Whole Run', fontsize=10)
        imageName = os.path.join(location, 'bulkGammaYRun.png')
        plt.savefig(imageName, dpi=150)
        plt.clf()
        plt.close('all')

    else:
        pass

    print('Images plotted #########')

    
if __name__ == "__main__":

    NewPlotFile(testLoc, 25, 8, 0.25, 5, 5, 8, False)
