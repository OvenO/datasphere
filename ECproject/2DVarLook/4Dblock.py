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
#import gc
#import weakref

def main():

    # gc.set_debug(gc.DEBUG_LEAK)
    # NEED TO HAVE EMPTY DIRECTORY "Data0"
    
    # before anything get us into the NormAll directory. this is the directory that will hold the
    # directories with the different data sets. We need to start keeping track of phase diagrams and
    # PC sections
    os.chdir(os.path.expanduser("~/Data/EC/4DBlock"))

    dt = .005 
    # total number of iterations to perform
    totIter = 100000
    totTime = totIter*dt
    time = pl.arange(0.0,totTime,dt)
    
    surf = 1.0
    coef = 1.5
    k = 1.0
    w = 1.0
    damp = .1
    g = .13

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
    yby = 3.0
    vxby = 10.0
    vyby = 10.0
    # define number of points in each direction
    numx =  2
    numy =  2
    numvx = 2 
    numvy = 2 
    # distance between points
    incx = xby/numx
    incy = yby/numy
    incvx = vxby/numvx
    incvy = vyby/numvy
    
    


    # initial conditions vector
    # set up: [xdot,ydot,x,y]
    x0 = pl.array([initvx,initvy,initx,inity])

    # lets name all the new folders just Data#/ 
    # don't want to write over any folder so check and see what number we are already on
    allDir = os.listdir(".")
    numdir = 0
    for l,z in enumerate(allDir):
        numdir = l        
    numdir += 1
    
    os.mkdir("Data" + str(numdir))
    os.chdir("Data" + str(numdir))

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
    datfile = open("data.txt","w")
    datfile.close()
    poinfile = open("poincar.txt","w")
    poinfile.close()

    # j will just keep track of how many times we run a particle. really jsut for the the labeling
    # in the data files
    j = 0
    
    checkT = 2*pl.pi/w

    for alpha in range(numx):
        for beta in range(numy):
            for kapa in range (numvx):
                for gama in range(numvy):
                    apx = ec.CentreLineApx(coef,k,w,damp,surf,g)

                    # itial conditions to next point
                    x0 = pl.array([initvx+alpha*incvx,initvy+beta*incvy,initx+kapa*incx,inity+gama*incy])

                    sol = odeint(apx.f,x0,time)
                    
                    for a in range(len(sol[:,0])):
                        sol[a,2] = sol[a,2]%modNum

                    # define the time slice information
                    checknextT = 0.0
                    checkPoint = 0

                    # coppy the data file to temp and read from the temp in order to rewrite the
                    # data file. "juggeling to avoid having to store huge arrays.
                    shutil.copyfile("data.txt","temp.txt")
                    readdat = open("temp.txt","r")
                    writedat = open("data.txt","w")

                    # do the same stuff for the poincare data
                    shutil.copyfile("poincar.txt","pointemp.txt")
                    readpoin = open("pointemp.txt","r")
                    writepoin = open("poincar.txt","w")
                    # this variarable (pstr) is only there to make the "lables" line more readable
                    pstr = "p"+str(j)
                    newlable = "%15s %15s %15s %15s"%(pstr+"_vx",pstr+"_vy",pstr+"_x",pstr+"_y")

                    # for the poincare section data we will just change the lable by adding a "TS"
                    # at the begining denoting "Time Slice"
                    newpoinlable = "%15s %15s %15s %15s"%("TS"+pstr+"_vx","TS"+pstr+"_vy","TS"+pstr+"_x","TS"+pstr+"_y")

                    # if the folder is empty need to just write first line
                    if alpha==beta==gama==kapa==0:
                        writedat.write(newlable + "\n")
                        writepoin.write(newpoinlable + "\n")
                        for i in range(len(sol)):
                            # add the first particles solution to the data file
                            toadd = "%15.6f %15.6f %15.6f %15.6f"%(sol[i,0],sol[i,1],sol[i,2],sol[i,3])
                            toadd += "\n"
                            writedat.write(toadd)

                            # this if statment cecks to see if we are at the right point to grab the
                            # timesliced data from the soulution. checkPoint is just the integer
                            # amount of "time" of a driving cycle or whatever
                            if i==checkPoint:
                                # add first poincare section to poincare data file
                                writepoin.write(toadd)

                                # For aditional coments (and outrage) please see 2Dvarlookcoef.py
                                # coments. (purpose -> stop propigation error)
                                checknextT += checkT

                                # the +.5 efectivly rounds the number apropriatly
                                checkPoint = int(checknextT/dt +.5)
                    else:
                        oldlable = readdat.readline()
                        # pop off the "\n" so we can append the new lable
                        oldlable = oldlable[:-1]
                        # append the new lable and write it to the file
                        writedat.write(oldlable + newlable + "\n")

                        # same as above just poincare data
                        oldpoinlable = readpoin.readline()
                        oldpoinlable = oldpoinlable[:-1]
                        writepoin.write(oldpoinlable + newpoinlable + "\n")

                        for i in range(len(sol)):
                            toadd = "%15.6f %15.6f %15.6f %15.6f"%(sol[i,0],sol[i,1],sol[i,2],sol[i,3])
                            curline = readdat.readline()
                            curline = curline[:-1]
                            newline = curline + toadd + "\n"
                            writedat.write(newline)

                            if i==checkPoint:

                                writepoin.write(newline)

                                # For aditional coments (and outrage) please see 2Dvarlookcoef.py
                                # coments. (purpose -> stop propigation error)
                                checknextT += checkT

                                # the +.5 efectivly rounds the number apropriatly
                                checkPoint = int(checknextT/dt +.5)

                    
                    # close the data file
                    readdat.close()
                    writedat.close()

                    # close the poincare file
                    readpoin.close()
                    writepoin.close()

                    # increase the number of particles we have run
                    j += 1
                    

    os.remove("temp.txt")
    os.chdir("..")

        #print("last phase space info for feedback \nin form: [xdot,ydot,x,y]")
        #print("this is the time:")
        #print(time[-1])


        #fig1 = pl.figure()
        #ax1 = fig1.add_subplot(111)
        #ax1.scatter(poinCarSx,poinCarSxdot,s=.1, label="Time Slices")
        #ax1.scatter([0.0,pl.pi,2*pl.pi],[0.0,0.0,0.0],color = "Red", marker="o",label="Electrodes")
        #ax1.set_title("Time Slice")
        #ax1.legend(loc = "best")
        #ax1.set_xlabel("$x$")
        #ax1.set_ylabel("Velocity $x$")
        #fig1.savefig("timeSL_x.png")
        #
        #fig13 = pl.figure()
        #ax13 = fig13.add_subplot(111)
        #ax13.scatter(poinCarSy,poinCarSydot,s=.1, label="Time Slices")
        #ax13.scatter([0.0,0.0,0.0],[0.0,0.0,0.0],color = "Red", marker="o",label="Electrodes")
        #ax13.set_title("Time Slice")
        #ax13.legend(loc = "best")
        #ax13.set_xlabel("$y$")
        #ax13.set_ylabel("Velocity $y$")
        #fig13.savefig("timeSL_y.png")

        #fig3 = pl.figure()
        #ax3 = fig3.add_subplot(111)
        ## just plot the last 10,000 pts
        #ax3.scatter(sol[(totIter-18000):,2],sol[(totIter-18000):,0],label="Last Part Phase \
        #        Space",s=.2)
        #ax3.scatter([0.0,pl.pi,2*pl.pi],[0.0,0.0,0.0],color = "Red", marker="o",label="Electrodes")
        #ax3.set_title("Phase Space")
        #ax3.legend(loc = "best")
        #ax3.set_xlabel("$x$")
        #ax3.set_ylabel("Velocity $x$")
        #fig3.savefig("phaseSpace_x.png")

        #fig12 = pl.figure()
        #ax12 = fig12.add_subplot(111)
        ## just plot the last 10,000 pts
        #ax12.plot(sol[(totIter-18000):,3],sol[(totIter-18000):,1],label="Last Part Phase Space")
        #ax12.scatter([0.0,0.0,0.0],[0.0,0.0,0.0],color = "Red", marker="o",label="Electrodes")
        #ax12.set_title("Phase Space")
        #ax12.legend(loc = "best")
        #ax12.set_xlabel("$y$")
        #ax12.set_ylabel("Velocity $y$")
        #fig12.savefig("phaseSpace_y.png")

        #fig14 = pl.figure()
        #ax14 = fig14.add_subplot(111)
        ## just plot the last 10,000 pts
        #ax14.plot(sol[(totIter-18000):,2],sol[(totIter-18000):,3],label="Last Part Real Space")
        #ax14.scatter([0.0,pl.pi,2*pl.pi],[0.0,0.0,0.0],color = "Red", marker="o",label="Electrodes")
        #ax14.set_title("Real Space")
        #ax14.legend(loc = "best")
        #ax14.set_xlabel("$x$")
        #ax14.set_ylabel("$y$")
        #fig14.savefig("realSpace.png")



        #fig4 = pl.figure()
        #ax4 = fig4.add_subplot(111)
        ## just plot the last 10,000 pts
        #ax4.plot(pl.arange(0,dt*18000,dt),sol[(totIter-18000):,2],label="Last Part $x(t)$")
        #ax4.scatter([0.0,0.0,0.0],[0.0,pl.pi,2*pl.pi],color = "Red", marker="o",label="Electrodes")
        #ax4.set_title("Position of Time")
        #ax4.legend(loc = "best")
        #ax4.set_xlabel("$\omega t$")
        #ax4.set_ylabel("$x$")
        #fig4.savefig("xoft.png")

        #fig11 = pl.figure()
        #ax11 = fig11.add_subplot(111)
        ## just plot the last 10,000 pts
        #ax11.plot(pl.arange(0,dt*18000,dt),sol[(totIter-18000):,3],label="Last Part $y(t)$")
        #ax11.scatter([0.0,0.0,0.0],[0.0,0.0,0.0],color = "Red", marker="o",label="Electrodes")
        #ax11.set_title("Position of Time")
        #ax11.legend(loc = "best")
        #ax11.set_xlabel("$\omega t$")
        #ax11.set_ylabel("$y$")
        #fig11.savefig("yoft.png")


        ### lets make a couple of more figures that are just the last 2000 pts
        ##fig5 = pl.figure()
        ##ax5 = fig5.add_subplot(111)
        ### just plot the last 2,000 pts
        ##ax5.plot(sol[(totIter-10000):,2],sol[(totIter-10000):,0],label="Last Part Phase Space")
        ##ax5.scatter([0.0,pl.pi,2*pl.pi],[0.0,0.0,0.0],color = "Red", marker="o",label="Electrodes")
        ##ax5.set_title("Phase Space")
        ##ax5.legend(loc = "best")
        ##ax5.set_xlabel("$x$")
        ##ax5.set_ylabel("Velocity $x$")
        ##fig5.savefig("restrictPhaseSpace_x.png")

        ### lets make a couple of more figures that are just the last 2000 pts
        ##fig10 = pl.figure()
        ##ax10 = fig10.add_subplot(111)
        ### just plot the last 2,000 pts
        ##ax10.plot(sol[(totIter-10000):,3],sol[(totIter-10000):,1],label="Last Part Phase Space")
        ##ax10.scatter([0.0,0.0,0.0],[0.0,0.0,0.0],color = "Red", marker="o",label="Electrodes")
        ##ax10.set_title("Phase Space")
        ##ax10.legend(loc = "best")
        ##ax10.set_xlabel("$y$")
        ##ax10.set_ylabel("Velocity $y$")
        ##fig10.savefig("restrictPhaseSpace_y.png")


        ##fig6 = pl.figure()
        ##ax6 = fig6.add_subplot(111)
        ### just plot the last 2,000 pts
        ##ax6.plot(pl.arange(0,dt*10000,dt),sol[(totIter-10000):,2],label="Last Part $x(t)$")
        ##ax6.scatter([0.0,0.0,0.0],[0.0,pl.pi,2*pl.pi],color = "Red", marker="o",label="Electrodes")
        ##ax6.set_title("Position of Time")
        ##ax6.legend(loc = "best")
        ##ax6.set_xlabel("$x$")
        ##ax6.set_ylabel("$\omega t$")
        ##fig6.savefig("restrictxoft.png")

        ##fig9 = pl.figure()
        ##ax9 = fig9.add_subplot(111)
        ### just plot the last 2,000 pts
        ##ax9.plot(pl.arange(0,dt*10000,dt),sol[(totIter-10000):,3],label="Last Part $y(t)$")
        ##ax9.scatter([0.0,0.0,0.0],[0.0,0.0,0.0],color = "Red", marker="o",label="Electrodes")
        ##ax9.set_title("Position of Time")
        ##ax9.legend(loc = "best")
        ##ax9.set_xlabel("$y$")
        ##ax9.set_ylabel("$\omega t$")
        ##fig9.savefig("restrictyoft.png")

        #fig7 = pl.figure()
        #ax7 = fig7.add_subplot(111)
        #ax7.scatter(poinCarSx[-500:],poinCarSxdot[-500:],s=.2, label="Last 500 Time Slices")
        #ax7.scatter([0.0,pl.pi,2*pl.pi],[0.0,0.0,0.0],color = "Red", \
        #        marker="o",label="Electrodes",s=.1)
        #ax7.set_title("Last Time Slices")
        #ax7.legend(loc = "best")
        #ax7.set_xlabel("$x$")
        #ax7.set_ylabel("Velocity $x$")
        #fig7.savefig("restricttimeSL_x.png")

        #fig8 = pl.figure()
        #ax8 = fig8.add_subplot(111)
        #ax8.scatter(poinCarSy[-500:],poinCarSydot[-500:],s=.2, label="Last 500 Time Slices")
        #ax8.scatter([0.0,0.0,0.0],[0.0,0.0,0.0],color = "Red", \
        #        marker="o",label="Electrodes",s=.1)
        #ax8.set_title("Last Time Slices")
        #ax8.legend(loc = "best")
        #ax8.set_xlabel("$y$")
        #ax8.set_ylabel("Velocity $y$")
        #fig8.savefig("restricttimeSL_y.png")


        
        
        #gc.collect()
        #print(weakref.getweakrefcount(time))
        #print(weakref.getweakrefcount(fig1))
        #print(weakref.getweakrefcount(ax1))

        # try to garbage collect by using the gc module

        # try to garbage colect (manualy) some of the matplotlib objects inorder to fix the memory
        # error
        #del fig1
        #del fig13
        #del fig3
        #del fig12
        #del fig14
        #del fig4
        #del fig11
        #del fig7
        #del fig8

        #del ax1
        #del ax13
        #del ax3
        #del ax12
        #del ax14
        #del ax4
        #del ax11
        #del ax7
        #del ax8

        #del sol
        #del poinCarSx
        #del poinCarSxdot
        #del poinCarSy


        #fig1 = None
        #fig13 = None
        #fig3 = None
        #fig12 = None
        #fig14 = None
        #fig4 = None
        #fig11 = None
        #fig7 = None
        #fig8 = None

        #ax1 = None
        #ax13 = None
        #ax3 = None
        #ax12 = None
        #ax14 = None
        #ax4 = None
        #ax11 = None
        #ax7 = None
        #ax8 = None

        #del sol
        #del poinCarSx
        #del poinCarSxdot
        #del poinCarSy
        #del poinCarSydot
       #del poinCarSydot


    

if __name__ == '__main__':
    main()
