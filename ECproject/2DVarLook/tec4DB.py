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
    os.chdir(os.path.expanduser("~/Data/EC/4DBlock"))

    dt = .05 
    # total number of iterations to perform
    totIter = 500
    totTime = totIter*dt
    time = pl.arange(0.0,totTime,dt)
    
    surf = 1.0
    coef = 5.0
    k = 1.0
    w = 2.0
    damp = .07
    g = .1

    # how many cells is till periodicity use x = n*pi/k (n must be even #)
    modNum = 2*pl.pi/k
    
    # make ec object
    elc = ec.electricCurtain()
    
    # define the lower left corner of block
    initvx = 0.0
    initvy = 0.0
    initx = 0.0
    inity = 0.0

    # define dimensions of block
    xby = 2.0*pl.pi
    yby = 10.0
    vxby = 0.0
    vyby = 0.0

    # define number of points in each direction
    numx =  100.0
    numy =  100.0
    numvx = 0.0
    numvy = 0.0

    # distance between points
    incx = xby/numx
    incy = yby/numy
    incvx = 0.0 
    incvy = 0.0 
    #incvx = vxby/numvx
    #incvy = vyby/numvy

    # because we don't need the WHOLE solution to watch what happens and the block goes forward in
    # time lets define a variable that determines how many points we skip before we keep one.
    # lets try this. can only implement it after the poincare section points are taken.
    # SKIP CAN NOT BE LESS THAN 1
    skip = 2
    

    # initial conditions vector
    # set up: [xdot,ydot,x,y]
    x0 = pl.array([initvx,initvy,initx,inity])

    # use time to name directories
    timestr = thetime.asctime()
    newtimestr = "Block_"
    for a,b in enumerate(timestr.split()):
        newtimestr+=b+"_"

    newtimestr = newtimestr.replace(":","") 

    os.mkdir(newtimestr)
    os.chdir(newtimestr)

    # make text file with all extra information
    outFile = open("info.dat","w")
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
            + "\ny: " +str(inity) \
            + "\nvx: " +str(initvx)\
            + "\nvy: " +str(initvy) \
            + "\nDimensions of block:" \
            + "\nxby: " + str(xby) \
            + "\nyby: " + str(yby) \
            + "\nvxby: " + str(vxby) \
            + "\nvyby: " + str(vyby) \
            + "\nNumber of points in respective directions:" \
            + "\nnumx: " + str(numx) \
            + "\nnumy: " + str(numy) \
            + "\nnumvx: " + str(numvx) \
            + "\nnumvy: " + str(numvy) )

    outFile.close()

    # make the files that will hold the data
    datfile = open("data.txt","a")
   
    datfile.write('''TITLE = "4DBLOCK"''')
    datfile.write("\n")
    datfile.write("VARIABLES = vx, vy, x, y")
    datfile.write("\n")

    filearr = pl.array([])

    for alpha in range(int(numx)):
        for beta in range(int(numy)):
            for kapa in range (int(numvx)):
                for gama in range(int(numvy)):

                    # make a file for the curent particles solution
                    curpstr =  str(alpha)+str(beta)+str(kapa)+str(gama)

                    # keeptrack of the fiels for later
                    filearr = pl.append(filearr,curpstr)
                    
                    curpdatfile = open(curpstr,"a")
                    

                    apx = ec.CentreLineApx(coef,k,w,damp,surf,g)

                    # itial conditions to next point
                    x0 = pl.array([initvx+alpha*incvx,initvy+beta*incvy,initx+kapa*incx,inity+gama*incy])

                    sol = odeint(apx.f,x0,time)
                    
                    for a in range(len(sol[:,0])):
                        sol[a,2] = sol[a,2]%modNum


                    for i in range(len(sol)):
                        # add the first particles solution to the data file
                        toadd = "%15.6f %15.6f %15.6f %15.6f"%(sol[i,0],sol[i,1],sol[i,2],sol[i,3])
                        toadd += "\n"
                        curpdatfile.write(toadd)


                    curpdatfile.close()

    
    for i in range(int(totIter/skip)-totIter%skip):
        # DATAPACKING=POINT should mean that the format is as such:
        # vx   vy   x    y
        # #'s  #'s  #'s  #'s
        # #'s  #'s  #'s  #'s
        # #'s  #'s  #'s  #'s
        # #'s  #'s  #'s  #'s
        datfile.write("ZONE   I="+str(len(filearr))+" DATAPACKING=POINT")
        datfile.write("\n")

        for a,b in enumerate(filearr):
            curfile = open(b,"r")
            lines = curfile.readlines()
            tofile = lines[i*skip]
            datfile.write(tofile)

            curfile.close()

    datfile.close()
        
    for apple,bannana in enumerate(filearr):
        os.remove(bannana)
        
    os.chdir("..")

if __name__ == '__main__':
    main()
