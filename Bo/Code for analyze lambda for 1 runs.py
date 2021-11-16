##### D:\Cloud Drive\Research\Liquid Crystals Project\Computational Material\annulus_R25_r4_n_40_k_20.0_HardBC\instanceRun
from typing import overload
import numpy as np
import os
import random as rd
import math
import numba

from matplotlib.patches import Ellipse, Rectangle, Circle
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib import cm, transforms
import matplotlib as mpl

def LambdaTimePlot(location):
        # open loc
    files = os.listdir(location)

    PosFiles = []
    StepNum = []
    for fileNames in files:
        
        if fileNames.startswith('PosArray'): #PosArray is 8 characters, .csv is 4 characters
            PosFiles.append(fileNames)
            if fileNames[8]=='.':
                StepNum.append(0)
            else:
                StepNum.append(int(fileNames[8:len(fileNames)-4]))

        else:
            pass
    LambdaArray=[]

    for file in PosFiles:

        fName = os.path.join(location,file)

        PosArray = np.loadtxt(fName,dtype=float,delimiter=',')

        Cos= np.mean(np.cos(2*PosArray[:,2]))
        Sin = np.mean(np.sin(2*PosArray[:,2]))

        LambdaArray.append(np.sqrt(Cos**2 + Sin**2))

    plt.plot(StepNum,LambdaArray,'x')
    plt.savefig(os.path.join(location,'Figure'),dpi=150)
    plt.clf()
    plt.close('all')
    
    #exclude first 10 points：
    lst=np.array([StepNum,LambdaArray])
    lstofpoints=lst.T
    lstofpoints_sorted=lstofpoints[np.lexsort(lstofpoints[:,::-1].T)]
    lstofpoints_selected=lstofpoints_sorted[11:]
    StepNum_selected=lstofpoints_selected.T[0]
    LambdaArray_selected=lstofpoints_selected.T[1]
    eqLambda=np.mean(LambdaArray_selected)
    text_file = open(os.path.join(location,'eqLambda.txt'),'w+')
    text_file.write("eqLambda= " + str(eqLambda))
    text_file.close()


'''
location=r'/Users/bozhang/Desktop/PHYS437A/code/rods_box_L8_n_400_b_0.5_HardBC'
files = os.listdir(location)
Posfiles=[]
Stepnum=[]
for filenames in files:
    if filenames.startswith('InstanceRun'): 
        Posfiles.append(filenames)
        if filenames[11]=='':
                Stepnum.append(0)
            else:
                Stepnum.append(int(filenames[11:len(filenames)-1]))
    else:
        pass

lstofEqLambda=[]
for i in Stepnum:
los=location+Stepnum
lstofEqLambda.append(LambdaTimePlot(los))
'''





runs=10
l=1.6
for i in range(runs+1):
    if i !=0:
        loc= r'/Users/bozhang/Desktop/PHYS437A/code/rods_box_L{0}_n_50_b_0.5_HardBC/instanceRun{1}'.format(l,i)
        LambdaTimePlot(loc)


#### Getting One equil. Lambdap

        ### 1 Run --> Lambda vs Time Plot


        ###  Ignore 10-15% of the run and find equilibrium Lambda -> 'equilibrium' Lambda 



### Run this again (identical same number of MC Steps and same % ignored)  say for X runs (could be 10 could be 5) --> each one gets an equilibrium lambda 

##### mean equilibrium Lambda is (Sum ( eqLambda_i) i =1 to X)/X 

###### std dev or variance  var(eqLambda) = 
