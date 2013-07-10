#!/usr/bin/python
#/opt/local/bin/python2.7
import numpy
import sys
sys.path.append("/users/o/m/omyers/puthere/lib/python2.7/site-packages")
sys.path.append("/users/o/m/omyers/datasphere/ECproject")
sys.path.append("/users/o/m/omyers/puthere/matplotlib/build/lib.linux-x86_64-2.7")
import ECclass as ec
import pylab as pl
import os
from scipy.integrate import odeint
import shutil
import time as thetime
import argparse
#import gc
#import weakref

def main():

    
    # before anything get us into the NormAll directory. this is the directory that will hold the
    # directories with the different data sets. We need to start keeping track of phase diagrams and
    # PC sections

    parser = argparse.ArgumentParser()

    parser.add_argument('--dir', action = 'store', dest = "dir",type = str,required = True)

    parser.add_argument("--surf", action = "store", dest =  "surf", type = float,required = True) 
    parser.add_argument("--coef", action = "store", dest =  "coef", type = float,required = True) 
    parser.add_argument("--k"   , action = "store", dest =  "k"   , type = float,required = True) 
    parser.add_argument("--w"   , action = "store", dest =  "w"   , type = float,required = True)  
    parser.add_argument("--damp", action = "store", dest =  "damp", type = float,required = True) 
    parser.add_argument("--g"   , action = "store", dest =  "g"   , type = float,required = True)  
    parser.add_argument("--dt"  , action = "store", dest =  "dt"   , type = float,required = True)  
    parser.add_argument("--bnum"   , action = "store", dest =  "bnum"   , type = int,required = True)  

    inargs = parser.parse_args()

    # define the lower left corner of block
    initvx = -4.0
    initx = 0

    # define dimensions of block
    xby = 6.0
    vxby = 8.0

    # define number of points in each direction
    numx =  30.0
    numvx = 30.0

    #numx =  80.0
    #numvx = 80.0

    # distance between points
    incx = xby/numx
    incvx = vxby/numvx


    dt = inargs.dt
    # total number of iterations to perform
    #totIter = 1000
    totIter = 3
    totTime = totIter*dt
    time = pl.arange(0.0,totTime,dt)


    os.chdir(os.path.expanduser("~/Data/EC/2DBlock/"+inargs.dir))

    # parameters for what should be chaotic but orderd trajectories
    surface = inargs.surf
    coefficient = inargs.coef
    wave_num    = inargs.k   
    omega    = inargs.w   
    damping = inargs.damp
    grav    = inargs.g   

    # how many cells is till periodicity use x = n*pi/k (n must be even #)
    modNum = 2*pl.pi/wave_num
    
    # make ec object
    elc = ec.electricCurtain()
    
    # define the lower left corner of block
    initialvx = initvx
    initialvy = 0.0
    initialx  = initx 
    initialy  = 1.0

    # define number of points in each direction
    num_of_x  = numx  
    num_of_vx = numvx 

    # distance between points
    increment_x  = incx  
    increment_vx = incvx 

    # initial conditions vector
    # set up: [xdot,ydot,x,y]
    x0 = pl.array([initialvx,initialvy,initialx,initialy])

    # a simple count variable to be used as file names. (was having repeat issues with file names
    # whun i was uing the sum of the indicies. silly owen.
    curp = 0 

    filearr = pl.array([])

    for kapa in range (int(num_of_vx)):
        for alpha in range(int(num_of_x)):

            # make a file for the curent particles solution
            curp += 1
            curpstr = str(inargs.bnum) +"_"+ str(curp)

            # keeptrack of the fiels for later
            filearr = pl.append(filearr,curpstr)
            
            curpdatfile = open(curpstr,"a")
            
            apx = ec.surfCentreLineApx(coefficient,wave_num,omega,damping)

            # itial conditions to next point
            x0 = pl.array([initialvx+kapa*increment_vx,initialvy,initialx+alpha*increment_x,initialy])

            sol = odeint(apx.f,x0,time)
            
            for a in range(len(sol[:,0])):
                sol[a,2] = sol[a,2]%modNum

            for i in range(len(sol)):
                # add the first particles solution to the data file
                toadd = "%15.6f %15.6f %15.6f %15.6f"%(sol[i,0],sol[i,1],sol[i,2],sol[i,3])
                toadd += "\n"
                curpdatfile.write(toadd)

            curpdatfile.close()
    
if __name__ == '__main__':
    main()
