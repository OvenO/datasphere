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
#import gc
#import weakref

def main():

    
    # before anything get us into the NormAll directory. this is the directory that will hold the
    # directories with the different data sets. We need to start keeping track of phase diagrams and
    # PC sections
    os.chdir(os.path.expanduser("~/Data/EC/2DBlock"))

    dt = .1
    # total number of iterations to perform
    totIter = 1200
    totTime = totIter*dt
    time = pl.arange(0.0,totTime,dt)
    
    surf = 1.0
    coef = .32
    k = 1.0
    w = 1.0
    damp = .1
    g = .1

    # how many cells is till periodicity use x = n*pi/k (n must be even #)
    modNum = 2*pl.pi/k
    
    # make ec object
    elc = ec.electricCurtain()
    
    # define the lower left corner of block
    initvx = -4.0
    initx = 0

    # define dimensions of block
    xby = 6.0
    vxby = 8.0

    # define number of points in each direction
    numx =  200.0
    numvx = 200.0

    # distance between points
    incx = xby/numx
    incvx = vxby/numvx

    # because we don't need the WHOLE solution to watch what happens and the block goes forward in
    # time lets define a variable that determines how many points we skip before we keep one.
    # lets try this. can only implement it after the poincare section points are taken.
    # SKIP CAN NOT BE LESS THAN 1
    skip = 2

    # initial conditions vector
    # set up: [xdot,ydot,x,y]
    x0 = pl.array([initvx,0.0,initx,1.0])

    # use time to name directories
    timestr = thetime.asctime()
    newtimestr = "Block_"
    for a,b in enumerate(timestr.split()):
        newtimestr+=b+"_"

    newtimestr = newtimestr.replace(":","") 

    os.mkdir(newtimestr)
    os.chdir(newtimestr)

    # make text file with all extra information
    outFile = open("info.txt","w")
    outFile.write("Info \n coefficient: " + str(coef) \
            + "\nwave number: " +str(k)\
            + "\nomega: " + str(w)\
            + "\ndamping: " + str(damp)\
            + "\ng: " + str(g)\
            + "\ntime step: " + str(dt)\
            + "\ntotal time: " + str(dt*totIter)\
            + "\ntotal iterations: " + str(totIter)\
            + "\n\nInitial Conditions Info: "\
            + "\nLower left corner of block: " \
            + "\nx: " +str(initx) \
            + "\nvx: " +str(initvx)\
            + "\nDimensions of block:" \
            + "\nxby: " + str(xby) \
            + "\nvxby: " + str(vxby) \
            + "\nNumber of points in respective directions:" \
            + "\nnumx: " + str(numx) \
            + "\nnumvx: " + str(numvx))

    outFile.close()

    # make the files that will hold the data
    datfile = open("data.txt","a")
   
    datfile.write('''TITLE = "2DBLOCK"''')
    datfile.write("\n")
    datfile.write("VARIABLES = vx, x")
    datfile.write("\n")

    # same thing for poindat.txt
    poinfile = open("poindat.txt","a")
   
    poinfile.write('''TITLE = "2DBLOCK"''')
    poinfile.write("\n")
    poinfile.write("VARIABLES = vx, x")
    poinfile.write("\n")


    filearr = pl.array([])

    # a simple count variable to be used as file names. (was having repeat issues with file names
    # whun i was uing the sum of the indicies. silly owen.
    curp = 0 

    for kapa in range (int(numvx)):
        for alpha in range(int(numx)):
            
            # make a file for the curent particles solution
            curp += 1
            curpstr = str(curp)

            # keeptrack of the fiels for later
            filearr = pl.append(filearr,curpstr)
            
            curpdatfile = open(curpstr,"a")

            apx = ec.surfCentreLineApx(coef,k,w,damp)

            # itial conditions to next point
            x0 = pl.array([initvx+kapa*incvx,0.0,initx+alpha*incx,1.0])

            sol = odeint(apx.f,x0,time)
            
            for a in range(len(sol[:,0])):
                sol[a,2] = sol[a,2]%modNum

            for i in range(len(sol)):
                # add the first particles solution to the data file
                toadd = "%15.6f %15.6f %15.6f %15.6f"%(sol[i,0],sol[i,1],sol[i,2],sol[i,3])
                toadd += "\n"
                curpdatfile.write(toadd)

            curpdatfile.close()
    
    checktime = 2*pl.pi/w
    checknext = 0.0
    checkpoint = 0

    for i in range(totIter):
        # DATAPACKING=POINT should mean that the format is as such:
        # vx   vy   x    y
        # #'s  #'s  #'s  #'s
        # #'s  #'s  #'s  #'s
        # #'s  #'s  #'s  #'s
        # #'s  #'s  #'s  #'s
        datfile.write("ZONE   I="+str(len(filearr))+" DATAPACKING=POINT")
        datfile.write("\n")

        if (i==checkpoint):
            poinfile.write("ZONE   I="+str(len(filearr))+" DATAPACKING=POINT")
            poinfile.write("\n")


        for a,b in enumerate(filearr):
            curfile = open(b,"r")
            lines = curfile.readlines()
            tofile = lines[i]
            datfile.write(tofile)

            if (i==checkpoint):
                poinfile.write(tofile)

            curfile.close()

        if (i==checkpoint):
            checknext += checktime
            checkpoint = int(checknext/dt+.5)
    

    datfile.close()
    poinfile.close()
        
    os.system("rm 1* 2* 3* 4* 5* 6* 7* 8* 9*")
        
    os.chdir("..")

if __name__ == '__main__':
    main()
