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


def cross(p1,p2,p3): # cross product
    x1=p2[0]-p1[0]
    y1=p2[1]-p1[1]
    x2=p3[0]-p1[0]
    y2=p3[1]-p1[1]
    return x1*y2-x2*y1    

    
def atan2(x,y):

    if x>0:
        val = np.arctan(y/x)
    elif x<0 and y>=0:
        val = np.arctan(y/x) + np.pi
    elif x<0 and y<0:
        val = np.arctan(y/x) - np.pi
    elif x==0 and y>0:
        val = np.pi/2
    elif x==0 and y<0:
        val = np.pi/2
    
    return val

def dist(x1,y1,x2,y2):
    return np.sqrt((x1-x2)**2 + (y1-y2)**2)

def LSIntersection(L1,L2): #dcide whether intersect for 2 lines segement[p1,p2] and [p3,p4]
    p1=L1[0]
    p2=L1[1]
    p3=L2[0]
    p4=L2[1]
    # use p1p2 and p3p4 as the diagonal line of two rectangle, if two rectangle has intersect:
    if(max(p1[0],p2[0])>=min(p3[0],p4[0]) and max(p3[0],p4[0])>=min(p1[0],p2[0]) and max(p1[1],p2[1])>=min(p3[1],p4[1]) and max(p3[1],p4[1])>=min(p1[1],p2[1])): 
        # if P1 and P2 lie on both sides of line segment p3p4
        # if P3 and P4 lie on both sides of line segment p1p2
        if(cross(p1,p2,p3)*cross(p1,p2,p4)<=0 and cross(p3,p4,p1)*cross(p3,p4,p2)<=0):
            overlap=1
        else:
            overlap=0
    else:
        overlap=0
    return overlap


def conv(seg,L):
        #convert segment(x,y,theta) to 2 end points[pminus,pplus]
            x=seg[0]
            y=seg[1]
            theta=seg[2]
            pminus=[x-math.cos(theta)*L*0.5,y-math.sin(theta)*L*0.5]
            pplus=[x+math.cos(theta)*L*0.5,y+math.sin(theta)*L*0.5]
            ni=[pminus,pplus]
            return ni

def RodOverlap(seg1,L1,seg2,L2): #dcide whether intersect for 2 lines segement[p1,p2] and [p3,p4]
        
        p1,p2=conv(seg1,L1)
        p3,p4=conv(seg2,L2)
        # use p1p2 and p3p4 as the diagonal line of two rectangle, if two rectangle has intersect:
        if(max(p1[0],p2[0])>=min(p3[0],p4[0]) and max(p3[0],p4[0])>=min(p1[0],p2[0]) and max(p1[1],p2[1])>=min(p3[1],p4[1]) and max(p3[1],p4[1])>=min(p1[1],p2[1])): 
            # if P1 and P2 lie on both sides of line segment p3p4
            # if P3 and P4 lie on both sides of line segment p1p2
            if(cross(p1,p2,p3)*cross(p1,p2,p4)<=0 and cross(p3,p4,p1)*cross(p3,p4,p2)<=0):
                overlap=True
            else:
                overlap=False
        else:
            overlap=False
        return overlap # 1 if intersect, 0 if not intersect






def RodCircleOverlap(xyt,L,R,InOut):

    if InOut=='out':
        p1=conv(xyt,L)[0]
        p2=conv(xyt,L)[1]

        if dist(p1[0],p1[1],0,0)>R or dist(p2[0],p2[1],0,0)>R:
            overlap=True
        else:
            overlap=False
    elif InOut=='in':
        p1=conv(xyt,L)[0]
        p2=conv(xyt,L)[1]

        if dist(p1[0],p1[1],0,0)<R or dist(p2[0],p2[1],0,0)<R:
            overlap=True
        else:
            overlap=False

    return overlap




# l is box size or circle radius depending on choice
# b is half of rod length
# n is number of particles


def init_Rod_Circle_Open(n,l,b):


    init_pos = np.zeros((n,3))

    #Correct overlapping positions
    sortedInitPos = [[i,init_pos[i,0],init_pos[i,1],init_pos[i,2],] for i in range(n)]

    for x in range(n):

        valid = False
        overlapVar = False

        while valid == False:

            radius = rd.uniform(0,l)
            angle = rd.uniform(0,2*np.pi)
            init_pos[x,0]= radius*np.cos(angle)
            init_pos[x,1]= radius*np.sin(angle)
            init_pos[x,2]=rd.uniform(0,2*np.pi)
            

            

            if (init_pos[x,0]**2 + init_pos[x,1]**2)>(l**2):
                    print('init - border overlap')
                    continue
            

            else:
                for j in range(n):

                    if j!=x:
                        if (RodOverlap(init_pos[x,:],2*b,init_pos[j,:],2*b)==True):
                            overlapVar = True
                            print('iter' + str(x)+' init - particle overlap')
                            break
                        else:
                            pass
                    else:
                        print('iter' + str(x)+'init - no particle overlap')
                        overlapVar = False

                if overlapVar==True:
                    continue
                else:
                    valid=True


               
    
    return init_pos

def init_Rod_Circle_Hard(n,l,b):


    init_pos = np.zeros((n,3))

    #Correct overlapping positions
    sortedInitPos = [[i,init_pos[i,0],init_pos[i,1],init_pos[i,2],] for i in range(n)]

    for x in range(n):

        valid = False
        overlapVar = False

        while valid == False:

            radius = rd.uniform(0,l)
            angle = rd.uniform(0,2*np.pi)
            init_pos[x,0]= radius*np.cos(angle)
            init_pos[x,1]= radius*np.sin(angle)
            init_pos[x,2]=rd.uniform(0,2*np.pi)
            

            

            if (init_pos[x,0]**2 + init_pos[x,1]**2)>(l**2) or RodCircleOverlap(init_pos[x,:],2*b,l,'out')==True:
                    print('init - border overlap')
                    continue
            

            else:
                for j in range(n):

                    if j!=x:
                        if (RodOverlap(init_pos[x,:],2*b,init_pos[j,:],2*b)==True):
                            overlapVar = True
                            print('iter' + str(x)+' init - particle overlap')
                            break
                        else:
                            pass
                    else:
                        print('iter' + str(x)+'init - no particle overlap')
                        overlapVar = False

                if overlapVar==True:
                    continue
                else:
                    valid=True


               
    
    return init_pos


def init_Rod_Box_Open(n,l,b):


    init_pos = np.zeros((n,3))

    #Correct overlapping positions
    sortedInitPos = [[i,init_pos[i,0],init_pos[i,1],init_pos[i,2],] for i in range(n)]

    for x in range(n):

        valid = False
        overlapVar = False

        while valid == False:

            init_pos[x,0]= rd.uniform(-l/2,l/2)
            init_pos[x,1]= rd.uniform(-l/2,l/2)
            init_pos[x,2]= rd.uniform(0,2*np.pi)
            

            

            if abs(init_pos[x,0])>l/2 or abs(init_pos[x,1])>l/2 :
                    print('init - border overlap')
                    continue
            

            else:
                for j in range(n):

                    if j!=x:
                        if (RodOverlap(init_pos[x,:],2*b,init_pos[j,:],2*b)==True):
                            overlapVar = True
                            print('iter' + str(x)+' init - particle overlap')
                            break
                        else:
                            pass
                    else:
                        print('iter' + str(x)+'init - no particle overlap')
                        overlapVar = False

                if overlapVar==True:
                    continue
                else:
                    valid=True


               
    
    return init_pos

def init_Rod_Box_Hard(n,l,b):


    init_pos = np.zeros((n,3))

    #Correct overlapping positions
    sortedInitPos = [[i,init_pos[i,0],init_pos[i,1],init_pos[i,2],] for i in range(n)]

    for x in range(n):

        valid = False
        overlapVar = False

        while valid == False:

             
            init_pos[x,0]= rd.uniform(-l/2,l/2)
            init_pos[x,1]= rd.uniform(-l/2,l/2)
            init_pos[x,2]= rd.uniform(0,2*np.pi)
            

            

            #if abs(init_pos[x,0])>l/2 or abs(init_pos[x,1])>l/2 or (RodOverlap(init_pos[x,:],2*b,[init_pos[x,0],-l/2,0],l)==True) or (RodOverlap(init_pos[x,:],2*b,[init_pos[x,0],l/2,0],l)==True) or (RodOverlap(init_pos[x,:],2*b,[-l/2,init_pos[x,1],np.pi/2],l)==True) or (RodOverlap(init_pos[x,:],2*b,[l/2,init_pos[x,1],np.pi/2],l)==True):
            if abs(init_pos[x,0])>l/2 or abs(init_pos[x,1])>l/2 or (RodOverlap(init_pos[x,:],2*b,[0,-l/2,0],l)==True) or (RodOverlap(init_pos[x,:],2*b,[0,l/2,0],l)==True) or (RodOverlap(init_pos[x,:],2*b,[-l/2,0,np.pi/2],l)==True) or (RodOverlap(init_pos[x,:],2*b,[l/2,0,np.pi/2],l)==True):

                    print('init - border overlap')
                    continue
            else:
                for j in range(n):

                    if j!=x:
                        if (RodOverlap(init_pos[x,:],2*b,init_pos[j,:],2*b)==True):
                            overlapVar = True
                            print('iter' + str(x)+' init - particle overlap')
                            break
                        else:
                            pass
                    else:
                        print('iter' + str(x)+'init - no particle overlap')
                        overlapVar = False

                if overlapVar==True:
                    continue
                else:
                    valid=True


               
    
    return init_pos

def init_aligned_Rod_Box_Hard(n,l,b):
    X=np.linspace(start=b,stop=l-b,num=n)
    Y=n*[0.5*l]
    THETA=n*[0.5*math.pi]
    init_pos= np.array(list(map(list, zip(*[X,Y,THETA]))))
    return init_pos
# mode 'S' to save image with fileName string
# mode 'D' to display image
# Need to set some small a (since you are plotting a rectangle)



def PlotterCircle_Rod(array,a,b,L,Mode,fileName):

    finalC = len(array)

    f,(ax) = plt.subplots(1,1)
    f.subplots_adjust(hspace=0,wspace=0)
    

    diagL = np.sqrt(a**2 + b**2)



    ells = [Rectangle(xy=[array[w,0]-diagL*np.cos(array[w,2]),array[w,1]-diagL*np.sin(array[w,2])],width=2*b, height=2*a,angle= (-np.arctan(a/b)+(array[w,2]))*360/(2*np.pi),color='black', fill=False,label=str(w+1)) for w in range(finalC)]




        



    for count in range(finalC):

        ax.add_patch(ells[count])
        #plt.text(array[count,0],array[count,1],str(count+1))
        print(count+1)
        print(array[count,])
        

    Boundary = Circle(xy=[0,0],radius=L,fill=False)

    ax.add_patch(Boundary)


    ax.set_xlim(-(1.25*L),L+(0.25*L))
    ax.set_ylim(-(1.25*L),L+(0.25*L))
    ax.set_aspect(aspect='equal')

    plt.title('N=' +str(finalC) + ' b=' + str(b) + ' k=' + str(b/a)+ ' R='+str(L)  + ' \n' + r'$\rho =$' +str((finalC)/(np.pi*(L**2))) )

    if Mode=="D":
        plt.show()
    elif Mode=="S":
        plt.savefig(fileName,dpi=150)
        plt.clf()
        plt.close('all')


def PlotterBox_Rod(array,a,b,L,Mode,fileName):

    finalC = len(array)

    f,(ax) = plt.subplots(1,1)
    f.subplots_adjust(hspace=0,wspace=0)
    

    diagL = np.sqrt(a**2 + b**2)



    ells = [Rectangle(xy=[array[w,0]-diagL*np.cos(array[w,2]),array[w,1]-diagL*np.sin(array[w,2])],width=2*b, height=2*a,angle= (-np.arctan(a/b)+(array[w,2]))*360/(2*np.pi),color='black', fill=True,label=str(w+1)) for w in range(finalC)]




        



    for count in range(finalC):

        ax.add_patch(ells[count])
        #plt.text(array[count,0],array[count,1],str(count+1))
        print(count+1)
        print(array[count,])
        

    Boundary = Rectangle(xy=[-L/2,-L/2],width=L,height=L,fill=False)

    ax.add_patch(Boundary)

    

    ax.set_xlim(-(1.25*L),L+(0.25*L))
    ax.set_ylim(-(1.25*L),L+(0.25*L))
    ax.set_aspect(aspect='equal')

    plt.title('N=' +str(finalC) + ' b=' + str(b) + ' k=' + str(b/a)+ ' L='+str(L) + ' \n' + r'$\rho =$' +str((finalC)/(np.pi*(L**2))) )



def MC_Ann_Hard(PosArray,d_pos,d_ang,steps,n,l,r2,b):
    
    #global moves
    #global accepted_moves

    moves = 0
    accepted_moves = 0
    imageName = "step_figure.png"
    numE=len(PosArray)
    file_name = "MonteCarlo_Annulus_SimNotes.txt"

    ogdir = os.getcwd()


    main_folder_name = os.path.join(ogdir,'rods_annulus_R{}_r{}_n_{}_b_{}_HardBC'.format(l,r2,numE,b))
    folder_name = os.path.join(main_folder_name,'instanceRun')


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
    complete_name = os.path.join(folder_name,fileNameArray)
    np.savetxt(complete_name,PosArray,delimiter=',')


    stepsArr = np.linspace(0,steps,steps+1)
    
   
    n=len(PosArray)
    


    
    fullImageName=os.path.join(folder_name,imageName)    

    ###Save Initial State

    fileNameArray = 'PosArray.csv'
    complete_name = os.path.join(folder_name,fileNameArray)
    np.savetxt(complete_name,PosArray,delimiter=',')

    minPos = .025*(b)
    minAng = .025


    fixCount = 0
    plotCount =0

    for u in range(steps):
    
    
        for w in range(n):

            if fixCount==(np.ceil(steps/10)):
                fixCount=0
                
                if d_ang>minAng and d_pos>minPos:
                    if accepted_moves/moves < 0.07:
                        d_ang= 0.2*d_ang
                        d_pos=0.2*d_pos
                    elif accepted_moves/moves < 0.17:
                        d_ang= 0.35*d_ang
                        d_pos=0.35*d_pos
                    elif accepted_moves/moves < 0.27:
                        d_ang= 0.6*d_ang
                        d_pos=0.6*d_pos
                    elif accepted_moves/moves < 0.37:
                        d_ang= 0.75*d_ang
                        d_pos=0.75*d_pos
                    elif accepted_moves/moves < 0.47:
                        d_ang= 0.9*d_ang
                        d_pos=0.9*d_pos
                    elif accepted_moves/moves > 0.97:
                        d_ang= 1.8*d_ang
                        d_pos=1.8*d_pos
                    elif accepted_moves/moves > 0.87:
                        d_ang= 1.65*d_ang
                        d_pos=1.65*d_pos
                    elif accepted_moves/moves > 0.77:
                        d_ang= 1.5*d_ang
                        d_pos=1.5*d_pos
                    elif accepted_moves/moves > 0.67:
                        d_ang= 1.35*d_ang
                        d_pos=1.35*d_pos
                    elif accepted_moves/moves > 0.57:
                        d_ang= 1.1*d_ang
                        d_pos=1.1*d_pos
                      
                else:

                    d_ang=minAng
                    d_pos=minPos




            x = d_pos*rd.uniform(-1,1)
            y = d_pos*rd.uniform(-1,1)
            t = d_ang*rd.uniform(-1,1)

            #print(x,y)

            testX = PosArray[w,0] + x
            testY = PosArray[w,1] + y
            testT = PosArray[w,2] + t


            if RodCircleOverlap([testX,testY,testT],2*b,l,'out')==True or RodCircleOverlap([testX,testY,testT],2*b,r2,'in')==True or dist(testX,testY,0,0)>l or dist(testX,testY,0,0)<(r2+(.025*(2*b))):
                
                overlapVar=True
            else:
                for j in range(n):
                
                    rij=np.sqrt((testX-PosArray[j,0])**2 + (testY-PosArray[j,1])**2)
                    
                    if rij<(2*b):
                        if j!=w and (RodOverlap([testX,testY,testT],2*b,PosArray[j,:],2*b)==True):
                            overlapVar = True
                            #print("particle overlap")
                            break
                        else:
                            overlapVar=False
                    else:
                        overlapVar=False
        
                    

            if overlapVar==True:
                pass
            elif overlapVar==False:
                accepted_moves +=1
                PosArray[w,0] = testX
                PosArray[w,1] = testY
                PosArray[w,2] = testT
            
            

            
            
            
            moves +=1

            fixCount+=1
            

            x=0
            y=0
            
        plotCount+=1

            


        


        if plotCount==(np.ceil(steps/10)):

            plotCount=0

            ######################################

            #### Periodically Save Data Snapshots #####
    
            
            

            fileNameArray = 'PosArray.csv'
            new_fileName=fileNameArray.split(".csv")[0] +str(u) + ".csv"
            complete_name = os.path.join(folder_name,new_fileName)
            np.savetxt(complete_name,PosArray,delimiter=',')

    density = n/((np.pi*l**2) - (np.pi*r2**2))

    
    print(accepted_moves)
    print(moves)
    print(accepted_moves/moves)


    complete_name = os.path.join(folder_name,file_name)

    text_file = open(complete_name,'w+')
    text_file.write("Parameters" + "\r\n")
    text_file.write("- - - - -" + "\r\n")
    text_file.write("Monte Carlo steps: " + str(steps) + "\r\n")
    text_file.write("R: " + str(l) + "\r\n")
    text_file.write("r: " + str(r2) + "\r\n")
    text_file.write("d_pos / step size: " + str(d_pos) + "\r\n")
    text_file.write("d_ang / step size: " + str(d_ang) + "\r\n")
    text_file.write("# of Rods: " + str(n) + "\r\n")
    text_file.write("density: " + str(density) + "\r\n")
    text_file.write("Half Rod Length: " + str(b) + "\r\n")
    text_file.write("Accepted Moves: " + str(accepted_moves) + "\r\n")
    text_file.write("Total Moves: " + str(moves) + "\r\n")
    text_file.write("Acceptance Rate: " + str(100*(accepted_moves/moves)) + " %" +"\r\n")

    text_file.close()


    fileNameArray = 'FinalPosArray.csv'
    complete_name = os.path.join(folder_name,fileNameArray)
    np.savetxt(complete_name,PosArray,delimiter=',')



def MC_Box_Hard(PosArray,d_pos,d_ang,steps,n,l,b):
    
    #global moves
    #global accepted_moves

    moves = 0
    accepted_moves = 0
    imageName = "step_figure.png"
    numE=len(PosArray)
    file_name = "MonteCarlo_Box_SimNotes.txt"

    ogdir = os.getcwd()


    main_folder_name = os.path.join(ogdir,'rods_box_L{}_n_{}_b_{}_HardBC'.format(l,numE,b))
    folder_name = os.path.join(main_folder_name,'instanceRun')


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
    complete_name = os.path.join(folder_name,fileNameArray)
    np.savetxt(complete_name,PosArray,delimiter=',')


    stepsArr = np.linspace(0,steps,steps+1)
    
   
    n=len(PosArray)
    


    
    fullImageName=os.path.join(folder_name,imageName)    

    ###Save Initial State

    fileNameArray = 'PosArray.csv'
    complete_name = os.path.join(folder_name,fileNameArray)
    np.savetxt(complete_name,PosArray,delimiter=',')

    minPos = .000025*(b)
    minAng = .000025


    fixCount = 0
    plotCount =0

    for u in range(steps):
    
    
        for w in range(n):

            if fixCount==(np.ceil(steps/100)):
                fixCount=0
                
                if d_ang>minAng and d_pos>minPos:
                    if accepted_moves/moves < 0.47:
                        d_ang= 0.9*d_ang
                        d_pos=0.9*d_ang
                    elif accepted_moves/moves < 0.37:
                        d_ang= 0.75*d_ang
                        d_pos=0.75*d_ang
                    elif accepted_moves/moves < 0.27:
                        d_ang= 0.6*d_ang
                        d_pos=0.6*d_ang
                    elif accepted_moves/moves < 0.17:
                        d_ang= 0.35*d_ang
                        d_pos=0.35*d_ang
                    elif accepted_moves/moves < 0.07:
                        d_ang= 0.2*d_ang
                        d_pos=0.2*d_ang
                    elif accepted_moves/moves > 0.57:
                        d_ang= 1.1*d_ang
                        d_pos=1.1*d_ang
                    elif accepted_moves/moves > 0.67:
                        d_ang= 1.35*d_ang
                        d_pos=1.35*d_ang
                    elif accepted_moves/moves > 0.77:
                        d_ang= 1.5*d_ang
                        d_pos=1.5*d_ang
                    elif accepted_moves/moves > 0.87:
                        d_ang= 1.65*d_ang
                        d_pos=1.65*d_ang
                    elif accepted_moves/moves > 0.97:
                        d_ang= 1.8*d_ang
                        d_pos=1.8*d_ang
                else:

                    d_ang=minAng
                    d_pos=minPos





            x = d_pos*rd.uniform(-1,1)
            y = d_pos*rd.uniform(-1,1)
            t = d_ang*rd.uniform(-1,1)

            #print(x,y)

            testX = PosArray[w,0] + x
            testY = PosArray[w,1] + y
            testT = PosArray[w,2] + t


            if abs(testX)>l/2 or abs(testY)>l/2 or (RodOverlap([testX,testY,testT],2*b,[0,-l/2,0],l)==True) or (RodOverlap([testX,testY,testT],2*b,[0,l/2,0],l)==True) or (RodOverlap([testX,testY,testT],2*b,[-l/2,0,np.pi/2],l)==True) or (RodOverlap([testX,testY,testT],2*b,[l/2,0,np.pi/2],l)==True):
                
                overlapVar=True
            else:
                for j in range(n):
                
                    rij=np.sqrt((testX-PosArray[j,0])**2 + (testY-PosArray[j,1])**2)
                    
                    if rij<(2*b):
                        if j!=w and (RodOverlap([testX,testY,testT],2*b,PosArray[j,:],2*b)==True):
                            overlapVar = True
                            #print("particle overlap")
                            break
                        else:
                            overlapVar=False
                    else:
                        overlapVar=False
        
                    

            if overlapVar==True:
                pass
            elif overlapVar==False:
                accepted_moves +=1
                PosArray[w,0] = testX
                PosArray[w,1] = testY
                PosArray[w,2] = testT
            
            

            
            
            
            moves +=1

            fixCount+=1
            

            x=0
            y=0
            
        plotCount+=1

            


        


        if plotCount==(np.ceil(steps/100)):

            plotCount=0

            ######################################

            #### Periodically Save Data Snapshots #####
    
            
            

            fileNameArray = 'PosArray.csv'
            new_fileName=fileNameArray.split(".csv")[0] +str(u) + ".csv"
            complete_name = os.path.join(folder_name,new_fileName)
            np.savetxt(complete_name,PosArray,delimiter=',')

    density = n/(l**2)

    
    print(accepted_moves)
    print(moves)
    print(accepted_moves/moves)


    complete_name = os.path.join(folder_name,file_name)

    text_file = open(complete_name,'w+')
    text_file.write("Parameters" + "\r\n")
    text_file.write("- - - - -" + "\r\n")
    text_file.write("Monte Carlo steps: " + str(steps) + "\r\n")
    text_file.write("L: " + str(l) + "\r\n")
    text_file.write("d_pos / step size: " + str(d_pos) + "\r\n")
    text_file.write("d_ang / step size: " + str(d_ang) + "\r\n")
    text_file.write("# of Rods: " + str(n) + "\r\n")
    text_file.write(" density: " + str(density) + "\r\n")
    text_file.write("Half Rod Length: " + str(b) + "\r\n")
    text_file.write("Accepted Moves: " + str(accepted_moves) + "\r\n")
    text_file.write("Total Moves: " + str(moves) + "\r\n")
    text_file.write("Acceptance Rate: " + str(100*(accepted_moves/moves)) + " %" +"\r\n")

    text_file.close()


    fileNameArray = 'FinalPosArray.csv'
    complete_name = os.path.join(folder_name,fileNameArray)
    np.savetxt(complete_name,PosArray,delimiter=',')



N=50
B=0.5 #(rod length will be 4)
R=3.5  #(Radius of circle outer boundary or box side length)
r=0
stepXY= .15*R  ##Usually start at 15-25% of R
stepTh= np.pi/3

MCSteps = 1500000

initialState= init_Rod_Box_Hard(N,R,B)
# initialState= init_aligned_Rod_Box_Hard(N,R,B)
#initialState= init_Rod_Circle_Hard(N,R,B)

MC_Box_Hard(initialState,stepXY,stepTh,MCSteps,N,R,B)
#MC_Ann_Hard(initialState,stepXY,stepTh,MCSteps,N,R,r,B)





