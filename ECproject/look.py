#!/usr/bin/python
#/opt/local/bin/python2.7
#import pylab as pl
import ECclass as ec
import os
import sys
sys.path.append("/users/o/m/omyers/puthere/lib/python2.7/site-packages")
from scipy.integrate import odeint
import scipy as pl

def main():
    # NEED TO HAVE EMPTY DIRECTORY "Data0"
    
    # before anything get us into the NormAll directory. this is the directory that will hold the
    # directories with the different data sets. We need to start keeping track of phase diagrams and
    # PC sections
    os.chdir(os.path.expanduser("~/Data/EC/NormAll"))
    # lets name all the new folders just Data#/ 
    # don't want to write over any folder so check and see what number we are already on
    allDir = os.listdir(".")
    numdir = 0
    for l,z in enumerate(allDir):
        numdir = l        
    numdir += 1
    
    os.mkdir("Data" + str(numdir))
    os.chdir("Data" + str(numdir))


    bifurcate = False

    dt = .005 
    # total number of iterations to perform
    totIter = 1000000
    totTime = totIter*dt
    time = pl.arange(0.0,totTime,dt)

    coef = 0.205
    k = 1.0
    w = 1.0
    damp = .1

    numParamChecks = 1

   # how many cells is till periodicity use x = n*pi/k (n must be even #)
    modNum = 2*pl.pi/k
    
    # increase coeff by:
    incCf = .0001
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

        # put in the bifurcation points for this run of the parameter only if we specify that we
        # want to get the bifurcation diagram
#        if (bifurcate):
#            ax8.scatter(pl.zeros(500)+coef,poinCarSx[-500:],s=.3)
#
#        ax8.set_title("Bifurcation Diagram")
#        ax8.set_ylabel("Position $x$")
#        ax8.set_xlabel("Coefficient")
#        fig8.savefig("bifurcation.png")
#
#        fig1 = pl.figure()
#        ax1 = fig1.add_subplot(111)
#        ax1.scatter(poinCarSx,poinCarSxdot,s=0.1, label="Time \
#                Slices",marker="o",linewidths=None,edgecolors= None)
#        ax1.scatter([0.0,pl.pi,2*pl.pi],[0.0,0.0,0.0],color = "Red",label="Electrodes",marker="o",)
#        ax1.set_title("Time Slice")
#        ax1.legend(loc = "best")
#        ax1.set_xlabel("$x$")
#        ax1.set_ylabel("Velocity $x$")
#        ax1.axis([4.75,5.25,.01,.02])
#        fig1.savefig("timeSL.pdf",dpi = 400)
#
#
#        fig3 = pl.figure()
#        ax3 = fig3.add_subplot(111)
#        # just plot the last 10,000 pts
#        ax3.scatter(sol[(totIter-10000):,2],sol[(totIter-10000):,0],label="Last Part Phase Space",s=.2)
#        ax3.scatter([0.0,pl.pi,2*pl.pi],[0.0,0.0,0.0],color = "Red", marker="o",label="Electrodes")
#        ax3.set_title("Phase Space")
#        ax3.legend(loc = "best")
#        ax3.set_xlabel("$x$")
#        ax3.set_ylabel("Velocity $x$")
#        fig3.savefig("phaseSpace.png")
#
#        fig4 = pl.figure()
#        ax4 = fig4.add_subplot(111)
#        # just plot the last 10,000 pts
#        ax4.plot(pl.arange(0,dt*10000,dt),sol[(totIter-10000):,2],label="Last Part $x(t)$")
#        ax4.scatter([0.0,0.0,0.0],[0.0,pl.pi,2*pl.pi],color = "Red", marker="o",label="Electrodes")
#        ax4.set_title("Position of Time")
#        ax4.legend(loc = "best")
#        ax4.set_xlabel("$x$")
#        ax4.set_ylabel("$\omega t$")
#        fig4.savefig("xoft.png")
#
#        # lets make a couple of more figures that are just the last 2000 pts
#        fig5 = pl.figure()
#        ax5 = fig5.add_subplot(111)
#        # just plot the last 2,000 pts
#        ax5.plot(sol[(totIter-2000):,2],sol[(totIter-2000):,0],label="Last Part Phase Space")
#        ax5.scatter([0.0,pl.pi,2*pl.pi],[0.0,0.0,0.0],color = "Red", marker="o",label="Electrodes")
#        ax5.set_title("Phase Space")
#        ax5.legend(loc = "best")
#        ax5.set_xlabel("$x$")
#        ax5.set_ylabel("Velocity $x$")
#        fig5.savefig("restrictPhaseSpace.png")
#
#        fig6 = pl.figure()
#        ax6 = fig6.add_subplot(111)
#        # just plot the last 2,000 pts
#        ax6.plot(pl.arange(0,dt*2000,dt),sol[(totIter-2000):,2],label="Last Part $x(t)$")
#        ax6.scatter([0.0,0.0,0.0],[0.0,pl.pi,2*pl.pi],color = "Red", marker="o",label="Electrodes")
#        ax6.set_title("Position of Time")
#        ax6.legend(loc = "best")
#        ax6.set_xlabel("$x$")
#        ax6.set_ylabel("$\omega$ $t$")
#        fig6.savefig("restrictxoft.png")
#
#        fig7 = pl.figure()
#        ax7 = fig7.add_subplot(111)
#        ax7.scatter(poinCarSx[-500:],poinCarSxdot[-500:],s=.2, label="Last 500 Time Slices")
#        ax7.scatter([0.0,pl.pi,2*pl.pi],[0.0,0.0,0.0],color = "Red", \
#                label="Electrodes",s=.1)
#        ax7.set_title("Last Time Slices")
#        ax7.legend(loc = "best")
#        ax7.set_xlabel("$x$")
#        ax7.set_ylabel("Velocity $x$")
#        fig7.savefig("restricttimeSL.png")


        
        
    #make data file with x
    dataOut = open("data.txt","w")
    for alpha,beta in enumerate(sol[(totIter-10000):,2]):
        dataOut.write(str(beta)+"\n")
    dataOut.close()
    
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
