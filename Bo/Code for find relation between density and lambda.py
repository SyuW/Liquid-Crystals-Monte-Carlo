from math import e
import os
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import latexify
# plot graph on density and lambda
# 1.read data from  L=2,2.5,3,3.5,4...
numE=50
b=0.5
loc=r'/Users/bozhang/Desktop/PHYS437A/code'
ListofL=[1.5, 1.6, 1.9, 1.75, 2, 2.5, 2.25, 2.75, 3, 3.5, 3.25, 3.75, 4]
Listofdensity=list(map(lambda x : numE/(x**2), ListofL))
stdeqLambda=[]
vareqLambda=[]
MeaneqLambda=[]
for i in range(len(ListofL)):
    main_folder_name = os.path.join(loc,'rods_box_L{}_n_{}_b_{}_HardBC'.format(ListofL[i],numE,b))
    file_name = os.path.join(main_folder_name,'eqLambdadata.txt')
    if os.path.exists(file_name):
        eqLambdatxt=open(file_name)
        stdeqLambdainfo=eqLambdatxt.readlines()[4]
        stdeqLambda.append(float((stdeqLambdainfo.strip().split(": "))[1]))
for i in range(len(ListofL)):
    main_folder_name = os.path.join(loc,'rods_box_L{}_n_{}_b_{}_HardBC'.format(ListofL[i],numE,b))
    file_name = os.path.join(main_folder_name,'eqLambdadata.txt')
    if os.path.exists(file_name):
        eqLambdatxt=open(file_name)
        vareqLambdainfo=eqLambdatxt.readlines()[5]
        vareqLambda.append(float((vareqLambdainfo.strip().split(": "))[1]))
        vareqLambda1=list(map(lambda x: x*1000,vareqLambda))  
for i in range(len(ListofL)):
    main_folder_name = os.path.join(loc,'rods_box_L{}_n_{}_b_{}_HardBC'.format(ListofL[i],numE,b))
    file_name = os.path.join(main_folder_name,'eqLambdadata.txt')
    if os.path.exists(file_name):
        eqLambdatxt=open(file_name)
        MeaneqLambdainfo=eqLambdatxt.readlines()[3]
        MeaneqLambda.append(float((MeaneqLambdainfo.strip().split(":"))[1]))         
print(stdeqLambda)
print(vareqLambda)
print(MeaneqLambda)

def ρ(N,L,a):
    return N*(L**2)/(a**2)
plt.plot(Listofdensity,vareqLambda1,'x')
plt.xlabel(r'$\rho$', fontsize=10)
plt.ylabel(r'$\sigma^{2}*10^{3}$', fontsize=10)
plt.title('N= ' +str(numE) + '  Rod length(L)= ' + str(int(2*b)) + "  Monte Carlo steps: " + '1,500,000',fontsize=10 )
plt.savefig(os.path.join(loc,'FigureofVarEqLambda'),dpi=150)
plt.clf()
plt.close('all')
plt.errorbar(Listofdensity,MeaneqLambda,stdeqLambda,fmt='o',ecolor='r',color='b',elinewidth=2,capsize=4)
plt.xlabel('density  '+r'$\mathrm{\rho}(N, L, a)=\frac{NL^{2}}{a^{2}}$'+'\n'+' '+'\n'+'', fontsize=10)
plt.ylabel(r'$\bar{\Lambda}$', fontsize=10)
plt.title('N= ' +str(numE) + '   Rod length(L)= ' + str(int(2*b))+ '\n' + "  MC steps= " + '1,500,000' +  '  Runs = 10 for each '+r'$\rho$',fontsize=10 )
plt.savefig(os.path.join(loc,'FigureofmeanEqLambda'),dpi=150)
plt.clf()
plt.close('all')
