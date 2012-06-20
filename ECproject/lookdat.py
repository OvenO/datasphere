#!/usr/bin/python
#/opt/local/bin/python2.7

#===============================================================================================
# The purpose of this program is to generate data files for a veriaty of initial conditions. We will
# loop this program in such a way to make a block of innitial conditions.
#===============================================================================================


import pylab as pl
import ECclass as ec
import os
from scipy.integrate import odeint

def main():
    # NEED TO HAVE EMPTY DIRECTORY "Data0"
    
    # before anything get us into the NormAllDat directory. this is the directory that will hold the
    # directories with the different data sets. We need to start keeping track of phase diagrams and
    # PC sections
    os.chdir("NormAllDat")
    # lets name all the new folders just Data#/ 
    # don't want to write over any folder so check and see what number we are already on
    allDir = os.listdir(".")
    numdir = 0
    for l,z in enumerate(allDir):
        numdir = l        
    numdir += 2
    
    os.mkdir("Data" + str(numdir))
    os.chdir("Data" + str(numdir))


    bifurcate = False

    dt = .005 
    # total number of iterations to perform
    totIter = 10000
    totTime = totIter*dt
    time = pl.arange(0.0,totTime,dt)

    coef = 0.275
    k = 1.0
    w = 1.0
    damp = .1

    #lower left corner of block
    initialx = 1.0
    initialvx = 0.0

    numptsx = 20
    numptsvx = 20

   # how many cells is till periodicity use x = n*pi/k (n must be even #)
    modNum = 2*pl.pi/k
    
    # increase x  by:
    incx = .1

    # increase vx  by:
    incvx = .1


    # make ec object
    elc = ec.electricCurtain()
    
    for alpha in range(numptsx):
        for beta in range(numptsvx):
            print("alpha: " + str(alpha))
            print("beta: " + str(beta))
            
            # initial conditions vector
            # set up: [xdot,ydot,x,y]
            x0 = pl.array([initialvx+beta*incvx,0.0,initialx+alpha*incx,1.0])
            # coef += incCf
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

            #make data file with x, vx 
            dataOut = open(str(alpha)+"_"+str(beta)+"data.txt","w")
            #first make line for designating coulumns
            dataOut.write("{0:s},{1:15s},{2:15s}".format("x position", "x velocity","time"))
            dataOut.write("\n")
            for line in range(totIter):
                dataOut.write("{0:.5f},{1:15.5f},{2:15.3f}".format(sol[line,2],sol[line,0],line*dt))
                dataOut.write("\n")
            dataOut.close()

            #make data file with poincare data
            poinDataOut = open(str(alpha)+"_"+str(beta)+"poincaredata.txt","w")
            #first make line for designating coulumns
            poinDataOut.write("{0:s},{1:15s},{2:15s}".format("t sliced x","t sliced vx", "time"))
            poinDataOut.write("\n")
            for line in range(len(poinCarSx)):
                poinDataOut.write("{0:.5f},{1:15.5f},{2:15.3f}".format(poinCarSx[line],poinCarSxdot[line],line*dt))
                poinDataOut.write("\n")
            poinDataOut.close()






        
        
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
