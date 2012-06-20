#!/usr/bin/python
#/opt/local/bin/python2.7
import pylab as pl
import ECclass as ec
import os
from scipy.integrate import odeint

def main():
    bifurcate = True

    dt = .01 
    # total number of iterations to perform
    totIter = 100000
    totTime = totIter*dt
    time = pl.arange(0.0,totTime,dt)

    coef = 0.2
    k = 1.0
    w = 1.0
    damp = .1

    numParamChecks = 2000

   # how many cells is till periodicity use x = n*pi/k (n must be even #)
    modNum = 2*pl.pi/k
    
    # increase coeff by:
    incCf = .00075
    # make ec object
    elc = ec.electricCurtain()


    # here is the bifurcation diagram (i hope)
    if (bifurcate):
        fig8 = pl.figure()
        ax8 = fig8.add_subplot(111)

    for j in range(numParamChecks):
        print(j)
        # initial conditions vector
        # set up: [xdot,ydot,x,y]
        x0 = pl.array([0.0,0.0,2.5,1.0])
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
    ax8.set_xlabel("Coefficient $j$")
    fig8.savefig("trybib.pdf",dpi = 400)
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
