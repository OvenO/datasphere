#!/usr/bin/python
#/opt/local/bin/python2.7
import sys
import numpy
sys.path.append("/users/o/m/omyers/datasphere/PeriodicPotential")
import PPclass as pp
sys.path.append("/users/o/m/omyers/puthere/lib/python2.7/site-packages")
sys.path.append("/users/o/m/omyers/puthere/matplotlib/build/lib.linux-x86_64-2.7")
import pylab as pl
import os
from scipy.integrate import odeint

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*--
# The output of a run will be fed in as the initial conditions for the next run
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-


def main():
    bifurcate = True

    dt = .005 
    # total number of iterations to perform
    totIter = 1000000
    totTime = totIter*dt
    time = pl.arange(0.0,totTime,dt)

    coef = 0.0
    k = 1.0
    w = 1.0
    damp = .1

    numParamChecks = 2000

   # how many cells is till periodicity use x = n*pi/k (n must be even #)
    modNum = 2*pl.pi/k
    
    # increase coeff by:
    incCf = .003


    # here is the bifurcation diagram (i hope)
    if (bifurcate):
        fig8 = pl.figure()
        ax8 = fig8.add_subplot(111)
        
    x0 = pl.array([0.0,1.8])

    for j in range(numParamChecks):
        print(j)
        # initial conditions vector
        # set up: [xdot,x]
        coef += incCf
        apx = pp.periodicPotential(coef,damp)
        sol = odeint(apx.f,x0,time)
        for a in range(len(sol[:,0])):
            sol[a,1] = sol[a,1]%modNum

        
        checkT = 0.0
        poinCarSx= pl.array([])
        poinCarSxdot = pl.array([])
        intst = 1
        checkPoint = 0
        while (checkPoint < totIter):
        #    print(checkPoint)
            poinCarSx = pl.append(poinCarSx,sol[checkPoint,1])
            poinCarSxdot = pl.append(poinCarSxdot,sol[checkPoint,0])
            
            checkT = intst*2*pl.pi/(w*dt)
            # the +.5 efectivly rounds the number apropriatly
            checkPoint = int(checkT +.5)
            intst += 1

        ax8.scatter(pl.zeros(50)+coef,poinCarSx[-50:],s=.001)

        x0 = sol[-1,:]
        

    ax8.set_title("Bifurcation Diagram")
    ax8.set_ylabel("Position $x$")
    ax8.set_xlabel("Coefficient $j$")
    fig8.savefig("trybib.pdf",dpi = 200)
    #pl.show()    

    # make text file with all extra information
    outFile = open("info.dat","w")
    outFile.write("Info \n coefficient: " + str(coef) \
            + "\nwave number: " +str(k)\
            + "\nomega: " + str(w)\
            + "\ndamping: " + str(damp)\
            + "\ntime step: " + str(dt)\
            + "\ntotal time: " + str(dt*totIter)\
            + "\ntotal iterations: " + str(totIter) \
            + "\nInitial conditions: " +str(x0))
    outFile.close()


if __name__ == '__main__':
    main()
