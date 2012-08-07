#!/usr/bin/python
#/opt/local/bin/python2.7
#import pylab as pl
import sys
sys.path.append("/users/o/m/omyers/puthere/lib/python2.7/site-packages")
import scipy as pl
import ECclass as ec
import os
from scipy.integrate import odeint

def main():
    # NEED TO HAVE EMPTY DIRECTORY "Data0"
    
    # before anything get us into the NormAll directory. this is the directory that will hold the
    # directories with the different data sets. We need to start keeping track of phase diagrams and
    # PC sections
    os.chdir(os.path.expanduser("~/Data/EC/NormAll"))

    bifurcate = True

    dt = .01 
    # total number of iterations to perform
    totIter = 100000
    totTime = totIter*dt
    time = pl.arange(0.0,totTime,dt)

    surf = 1.0
    coef = 0.205
    k = 1.0
    w = 1.0
    damp = .1

    numParamChecks = 1000

   # how many cells is till periodicity use x = n*pi/k (n must be even #)
    modNum = 2*pl.pi/k
    
    # increase coeff by:
    incCf = .00001
    # make ec object
    elc = ec.electricCurtain()

    # lets name all the new folders just Data#/ 
    # don't want to write over any folder so check and see what number we are already on
    allDir = os.listdir(".")
    numdir = 0
    for l,z in enumerate(allDir):
        numdir = l        
    numdir += 2
    
    os.mkdir("Data" + str(numdir))
    os.chdir("Data" + str(numdir))

    


    # here is the bifurcation diagram (i hope)
    if (bifurcate):
        fig8 = pl.figure()
        ax8 = fig8.add_subplot(111)

    for j in range(numParamChecks):
        print(j)
        # initial conditions vector
        # set up: [xdot,ydot,x,y]
        x0 = pl.array([0.0,0.0,2.5,surf])
        coef += incCf
        apx = ec.surfCentreLineApx(coef,k,w,damp)
        sol = odeint(apx.f,x0,time)
        for a in range(len(sol[:,0])):
            sol[a,2] = sol[a,2]%modNum

        
        checkT = 0.0
        poinCarSx= pl.array([])
        poinCarSxdot = pl.array([])
        intst = 1
        checkPoint = 0
        while (checkPoint < totIter):
        #    print(checkPoint)
            poinCarSx = pl.append(poinCarSx,sol[checkPoint,2])
            poinCarSxdot = pl.append(poinCarSxdot,sol[checkPoint,0])
            checkT = intst*2*pl.pi/(w*dt)
            # the +.5 efectivly rounds the number apropriatly
            checkPoint = int(checkT +.5)
            intst += 1

        ax8.scatter(pl.zeros(50)+coef,poinCarSx[-50:],s=.001)

    ax8.set_title("Bifurcation Diagram")
    ax8.set_ylabel("Position $x$")
    ax8.set_xlabel("Coefficient $\jmath")
    fig8.savefig("bifurcation.pdf",dpi = 400)
        
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

    os.chdir("..")
    

if __name__ == '__main__':
    main()
