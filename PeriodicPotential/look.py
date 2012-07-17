#!/usr/bin/python
#/opt/local/bin/python2.7
import pylab as pl
import PPclass as pp
import os
from scipy.integrate import odeint

def main():
    # NEED TO HAVE EMPTY DIRECTORY "Data0"
    
    # before anything get us into the NormAll directory. this is the directory that will hold the
    # directories with the different data sets. We need to start keeping track of phase diagrams and
    # PC sections
    os.chdir(os.path.expanduser("~/Data/PP/NormAll"))
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

    dt = .05 
    # total number of iterations to perform
    totIter = 10000
    totTime = totIter*dt
    time = pl.arange(0.0,totTime,dt)

    coef = 1.0
    damp = .1


    numParamChecks = 1

    inccoef = .05
    
    # Do we want to wrap tragectories into one or two "unit cells"?
    wrapeUp = True

    # and if we DO want to wrap it all into a box how many cells is till periodicity use x = n*pi (n must be even #)
    modNum = 2*pl.pi
    

    for j in range(numParamChecks):
        # initial conditions vector
        # set up: [xdot,x]
        x0 = pl.array([0.0,2.5])
        apx = pp.periodicPotential(coef,damp)
        sol = odeint(apx.f,x0,time)

        # if we want to wrap it all up here is where that happens
        if(wrapup):
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
            checkT = intst*2*pl.pi/dt
            # the +.5 efectivly rounds the number apropriatly
            checkPoint = int(checkT +.5)
            intst += 1


        print(poinCarSx)
        
        fig1 = pl.figure()
        ax1 = fig1.add_subplot(111)
        ax1.scatter(poinCarSx,poinCarSxdot,s=0.1, label="Time \
                Slices",marker="o",linewidths=None)
        ax1.scatter([0.0,pl.pi,2*pl.pi],[0.0,0.0,0.0],color = "Red",label="Electrodes",marker="o",)
        ax1.set_title("Time Slice")
        ax1.legend(loc = "best")
        ax1.set_xlabel("$x$")
        ax1.set_ylabel("Velocity $x$")
        #ax1.axis([4.75,5.25,.01,.02])
        fig1.savefig("timeSL.pdf",dpi = 300)


        fig3 = pl.figure()
        ax3 = fig3.add_subplot(111)
        # just plot the last 10,000 pts
        ax3.scatter(sol[(totIter-10000):,1],sol[(totIter-10000):,0],label="Last Part Phase Space",s=.2)
        ax3.scatter([0.0,pl.pi,2*pl.pi],[0.0,0.0,0.0],color = "Red", marker="o",label="Electrodes")
        ax3.set_title("Phase Space")
        ax3.legend(loc = "best")
        ax3.set_xlabel("$x$")
        ax3.set_ylabel("Velocity $x$")
        fig3.savefig("phaseSpace.png")

        fig4 = pl.figure()
        ax4 = fig4.add_subplot(111)
        # just plot the last 10,000 pts
        ax4.plot(pl.arange(0,dt*10000,dt),sol[(totIter-10000):,1],label="Last Part $x(t)$")
        ax4.scatter([0.0,0.0,0.0],[0.0,pl.pi,2*pl.pi],color = "Red", marker="o",label="Electrodes")
        ax4.set_title("Position of Time")
        ax4.legend(loc = "best")
        ax4.set_xlabel("$x$")
        ax4.set_ylabel("$\omega t$")
        fig4.savefig("xoft.png")

        # lets make a couple of more figures that are just the last 2000 pts
        fig5 = pl.figure()
        ax5 = fig5.add_subplot(111)
        # just plot the last 2,000 pts
        ax5.plot(sol[(totIter-2000):,1],sol[(totIter-2000):,0],label="Last Part Phase Space")
        ax5.scatter([0.0,pl.pi,2*pl.pi],[0.0,0.0,0.0],color = "Red", marker="o",label="Electrodes")
        ax5.set_title("Phase Space")
        ax5.legend(loc = "best")
        ax5.set_xlabel("$x$")
        ax5.set_ylabel("Velocity $x$")
        fig5.savefig("restrictPhaseSpace.png")

        fig6 = pl.figure()
        ax6 = fig6.add_subplot(111)
        # just plot the last 2,000 pts
        ax6.plot(pl.arange(0,dt*2000,dt),sol[(totIter-2000):,1],label="Last Part $x(t)$")
        ax6.scatter([0.0,0.0,0.0],[0.0,pl.pi,2*pl.pi],color = "Red", marker="o",label="Electrodes")
        ax6.set_title("Position of Time")
        ax6.legend(loc = "best")
        ax6.set_xlabel("$x$")
        ax6.set_ylabel("$\omega$ $t$")
        fig6.savefig("restrictxoft.png")

        fig7 = pl.figure()
        ax7 = fig7.add_subplot(111)
        ax7.scatter(poinCarSx[-500:],poinCarSxdot[-500:],s=.2, label="Last 500 Time Slices")
        ax7.scatter([0.0,pl.pi,2*pl.pi],[0.0,0.0,0.0],color = "Red", \
                label="Electrodes",s=.1)
        ax7.set_title("Last Time Slices")
        ax7.legend(loc = "best")
        ax7.set_xlabel("$x$")
        ax7.set_ylabel("Velocity $x$")
        fig7.savefig("restricttimeSL.png")
       
        coef += inccoef


        
        
    #make data file with x
    #dataOut = open("data.txt","w")
    #for alpha,beta in enumerate(sol[(totIter-10000):,1]):
    #    dataOut.write(str(beta)+"\n")
    #dataOut.close()
    
    # make text file with all extra information
    outFile = open("info.dat","w")
    outFile.write("Info \n coefficient: " + str(coef) \
            + "\ndamping: " + str(damp)\
            + "\ntime step: " + str(dt)\
            + "\ntotal time: " + str(dt*totIter)\
            + "\ntotal iterations: " + str(totIter) \
            + "\nInitial conditions: " +str(x0))
    outFile.close()
    
    os.chdir("..")
    

if __name__ == '__main__':
    main()
