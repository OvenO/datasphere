#!/usr/bin/python
#/opt/local/bin/python2.7
#import pylab as pl
import os
import sys
sys.path.append("/users/o/m/omyers/puthere/lib/python2.7/site-packages")
sys.path.append("/users/o/m/omyers/puthere/matplotlib/build/lib.linux-x86_64-2.7")
sys.path.append("/users/o/m/omyers/datasphere/ECproject")
import pylab as pl
import ECclass as ec
from scipy.integrate import odeint
import matplotlib as mpl

def main():
    datfile = open("bifdat.txt","w")
    bifurcate = True

    dt = .05 
    # total number of iterations to perform
    totIter = 50000
    totTime = totIter*dt
    time = pl.arange(0.0,totTime,dt)
    origtime = time

    coef = .255
    #coef = .26
    k = 1.0
    # NEVER CHANGE FROM 1 !!!!!! will chage modulus from feedback info
    w = 1.0
    damp = .10

    numParamChecks = 500
    #numParamChecks = 1000

   # how many cells is till periodicity use x = n*pi/k (n must be even #)
    modNum = 2*pl.pi/k
    
    # increase coeff by:
    #incCf = .0009
    incCf = .00004
    # make ec object
    elc = ec.electricCurtain()

    # font stuff 
    font = {'family' : 'normal', 'weight' : 'bold', 'size' : 22}
    mpl.rc('font',**font) 

    # here is the bifurcation diagram (i hope)
    fig8 = pl.figure()
    ax8 = fig8.add_subplot(111)

    # initial conditions vector
    # set up: [xdot,ydot,x,y]
    x0 = pl.array([0.0,0.0,3.5,1.0])
    
    # this is to keep track of the time and look for intervals of 2*pi in our new way 
    l = 0
    for j in range(numParamChecks):
        print(j)
        coef += incCf
        apx = ec.surfCentreLineApx(coef,k,w,damp)
        sol = odeint(apx.f,x0,time)
        for a in range(len(sol[:,0])):
            sol[a,2] = sol[a,2]%modNum
            datfile.write(str(sol[a,2])+"\n")

        
        poinCarSx= pl.array([])
        poinCarSxdot = pl.array([])
        
        # this is just to go through the array corectly
        othercount = 0
        while (othercount < totIter):
            #print("l is: " + str(l))        
            #print("l%tot is: " + str(l%totIter))
            #print("othercount is: " + str(othercount))

            # new way of looking for intervals of 2*pi
            if ((l*dt)%(2*pl.pi)<dt):
                #print("made it: " + str(l))
                poinCarSx = pl.append(poinCarSx,sol[l%totIter,2])
                poinCarSxdot = pl.append(poinCarSxdot,sol[l%totIter,0])

            othercount += 1
            l+=1
            
        #print("len poin: " +str(len(poinCarSx)))
        ax8.scatter(pl.zeros(50)+coef,poinCarSx[-50:],s=.1)
        x0 = sol[-1,:]
        #time = pl.arange(time[-1]+dt,time[-1]+totTime+dt,dt)
        time = pl.arange(time[-1]+dt,time[-1]+totTime+dt,dt)
         
    
    ax8.set_title("Bifurcation Diagram")
    ax8.set_ylabel("Position $x$")
    ax8.set_xlabel("Coefficient $A'$")
    #fig8.savefig("trybif.eps")
    fig8.savefig("010313trybif14.pdf",dpi = 300)
    #pl.show()    

    # make text file with all extra information
    outFile = open("info.txt","w")
    outFile.write("Info \n coefficient: " + str(coef) \
            + "\nwave number: " +str(k)\
            + "\nomega: " + str(w)\
            + "\ndamping: " + str(damp)\
            + "\ntime step: " + str(dt)\
            + "\ntotal time: " + str(dt*totIter)\
            + "\ntotal iterations: " + str(totIter) \
            + "\nInitial conditions: " +str(x0))
    outFile.close()
    
    datfile.close() 

if __name__ == '__main__':
    main()
