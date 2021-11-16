from math import e
import os
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
#Gain data from different runs
#Save it in a list and caulculate mean and std
#graph Normal distribution
loc= r'/Users/bozhang/Desktop/PHYS437A/code'
l=1.6
runs=10
numE=50
b=0.5
main_folder_name = os.path.join(loc,'rods_box_L{}_n_{}_b_{}_HardBC'.format(l,numE,b))
eqLambda=[0]*runs
for i in range(runs+1):
    if i !=0:
        folder_name = os.path.join(main_folder_name,'instanceRun{}'.format(i)) #i is the number after instanceRun
        if os.path.exists(folder_name):
            file_name = "eqLambda.txt"
            eqLambdaposition=os.path.join(folder_name,file_name)
            eqLambdatxt=open(eqLambdaposition)
            eqLambda[i-1]=float((eqLambdatxt.read())[10:])
print(eqLambda)
meaneqLambda=np.mean(eqLambda)
stdeqLambda=np.std(eqLambda, ddof=1)
runs=len(eqLambda)
print(meaneqLambda)
print(stdeqLambda)

# *****density can't be larger than 1 for 2D



def demoNormalDistri():
    x = np.linspace(meaneqLambda - 5 * stdeqLambda, meaneqLambda + 5 * stdeqLambda, 10000)
    y = np.exp(-(x - meaneqLambda) ** 2 / (2 * stdeqLambda ** 2)) / (np.sqrt(2 * np.pi) * stdeqLambda)
    plt.plot(x, y)
    plt.grid(True)
    plt.xlabel(r'$\bar{\Lambda}$',fontsize=10)
    plt.ylabel('Probability Density',fontsize=10)
    plt.title('N= ' +str(numE) + '  Rod length= ' + str(int(2*b)) + "  MCsteps: " + '1,500,000'+'\n'+'Box length= '+str(l) + r'  $\rho =$' +str(numE*((2*b)**2)/(l**2))[:5] \
        +'  Runs= '+str(runs)+ '\n' + 'MeanEqLambda= ' +str(meaneqLambda)[:5] + '  StdEqLambda= ' + str(stdeqLambda)[:5], fontsize=10)
    plt.savefig(os.path.join(loc,main_folder_name,'Normal Distribution of EqLambda'),dpi=150)
    plt.clf()
    plt.close('all')
demoNormalDistri()

def demoDistri():
    x=list(range(10))
    y=eqLambda
    plt.plot(x,y,'x')
    plt.xticks([])
    plt.yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9, 1])
    plt.xlabel('Runs',fontsize=10)
    plt.ylabel(r'$\bar{\Lambda}$',fontsize=10)
    plt.title('N= ' +str(numE) + '  Rod length= ' + str(int(2*b)) + "  MCsteps: " + '1,500,000'+'\n'+'Box length= '+str(l) + r'  $\rho =$' +str(numE*((2*b)**2)/(l**2))[:5] \
        + '  MeanEqLambda= ' +str(meaneqLambda)[:5] + '  StdEqLambda= ' + str(stdeqLambda)[:5], fontsize=10)
    plt.savefig(os.path.join(loc,main_folder_name,'Distribution of EqLambda'),dpi=150)
demoDistri()

# confidence range:
confirange=[meaneqLambda-1.96*stdeqLambda/np.sqrt(np.pi),meaneqLambda+1.96*stdeqLambda/np.sqrt(np.pi)]

text_file = open(os.path.join(main_folder_name,'eqLambdadata.txt'),'w+')
text_file.write("'equalibrium' lambda" + "\r\n")
text_file.write("- - - - -" + "\r\n")
text_file.write("list of eqLambda: " + str(eqLambda) + "\r\n")
text_file.write("mean of eqLambda:" + str(meaneqLambda) + "\r\n")
text_file.write("std of eqLambda: " + str(stdeqLambda) + "\r\n")
text_file.write("var of eqLambda: " + str(stdeqLambda**2) + "\r\n")
text_file.write("95% confidence interval: " + str(confirange) + "\r\n")
text_file.write("Monte Carlo steps: " + '1,500,000' + "\r\n")
text_file.write("L: " + str(l) + "\r\n")
text_file.write("Rod Length: " + str(2*b) + "\r\n")
text_file.write("# of Rods: " + str(numE) + "\r\n")
text_file.write("density: " + str(numE*((2*b)**2)/(l**2)) + "\r\n")


text_file.close()