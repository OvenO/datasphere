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
import argparse
import time as thetime
#import gc
#import weakref

def main():
    # gc.set_debug(gc.DEBUG_LEAK)
    # NEED TO HAVE EMPTY DIRECTORY "Data0"
    
    # before anything get us into the NormAll directory. this is the directory that will hold the
    # directories with the different data sets. We need to start keeping track of phase diagrams and
    # PC sections
    os.chdir(os.path.expanduser("~/Data/EC/2DNormAllVar/VarCoef"))
    
    # YAY ARGPARSE
    parser = argparse.ArgumentParser()
    # this is what we want g to be
    parser.add_argument("-g",action="store",dest="g",type=float,default=.02) 
    # the damping coefficient...
    parser.add_argument("-b",action="store",dest="b",type=float,default=.05)
            # stable for .1 (right now 11/8/12) 

    # starting coefficient...
    parser.add_argument("-c",action="store",dest="c",type=float,default=2.5)
    # number of parameter checks...
    parser.add_argument("-n",action="store",dest="n",type=int,default=20)
    # increase in coefficent per param check...
    parser.add_argument("-i",action="store",dest="i",type=float,default= .02)
    # pass in directory name also :/
    timestring = thetime.asctime()
    newtimestr = ""
    for l,m in enumerate(timestring.split()):
        newtimestr += m + "_"
    newtimestr = newtimestr.replace(":","")
    parser.add_argument("-d",action="store",dest="d",type=str,default=newtimestr)

    inargs = parser.parse_args()

    dt = .01 
    # total number of iterations to perform
    totIter = 80000
    totTime = totIter*dt
    time = pl.arange(0.0,totTime,dt)
    
    surf = 1.0
    coef = inargs.c
    k = 1.0
    w = 2.0
    damp = inargs.b
    g = inargs.g

    numParamChecks = inargs.n

    # how many cells is till periodicity use x = n*pi/k (n must be even #)
    modNum = 2*pl.pi/k
    
    # increase coeff by:
    incCf = inargs.i
    # make ec object
    elc = ec.electricCurtain()
    
    # initial conditions
    initvx = 0.0
    initvy = 0.0
    initx = 1.5
    inity = 2.2

    # initial conditions vector
    # set up: [xdot,ydot,x,y]
    x0 = pl.array([initvx,initvy,initx,inity])

    # NEW DATA FILE MANAGMENT SYSTEM
    # Date and time... thats all for now
    curtime = inargs.d
    os.mkdir(curtime)
    os.chdir(curtime)

    # make text file with all extra information
    outFile = open("info.txt","w")
    outFile.write("Info \n initial coefficient: " + str(coef) \
            + "\nincreace in coefficient: " + str(incCf)\
            + "\nwave number: " +str(k)\
            + "\nomega: " + str(w)\
            + "\ndamping: " + str(damp)\
            + "\ng: " + str(g)\
            + "\ntime step: " + str(dt)\
            + "\ntotal time: " + str(dt*totIter)\
            + "\ntotal iterations: " + str(totIter)\
            + "\nInitial Conditions: \n" +
            "initial x: " +str(initx) \
            +"\ninitial y: " +str(inity) \
            +"\ninitial vx: " +str(initvx)\
            +"\ninitial vy: " +str(initvy) )
    outFile.close()

    # make the files that will hold the data
    datfile = open("data.txt","w")
    datfile.close()
    poinfile = open("poincar.txt","w")
    poinfile.close()


    #checkT = 2*pl.pi/w
    #checkPoint = 0
    l = 0
    for j in range(numParamChecks):
        apx = ec.CentreLineApx(coef,k,w,damp,surf,g)
        sol = odeint(apx.f,x0,time)
        
        for a in range(len(sol[:,0])):
            sol[a,2] = sol[a,2]%modNum

        #coef += incCf
        damp += incCf

        #newstart = checkPoint%totIter

        #checknextT = 0.0
        #checkPoint = 0
        

        # copy the data file to temp and read from the temp inorder to rewrite the data file.
        # "Juggeling" to avoid having to store huge arrays.
        shutil.copyfile("data.txt","temp.txt")
        readdat = open("temp.txt","r")
        writedat = open("data.txt","w")
    
        # do the same stuff for the poincare data
        shutil.copyfile("poincar.txt","pointemp.txt")
        readpoin= open("pointemp.txt","r")
        writepoin = open("poincar.txt","w")
        # this variarable (pstr) is only there to make the "lables" line more readable
        pstr = "p"+str(j)
        newlable = "%15s %15s %15s %15s"%(pstr+"_vx",pstr+"_vy",pstr+"_x",pstr+"_y")

        # for the poincare section data we will just change the lable by adding a "TS" at the
        # begining detoting "Time Slice"
        newpoinlable = "%15s %15s %15s %15s"%("TS"+pstr+"_vx","TS"+pstr+"_vy","TS"+pstr+"_x","TS"+pstr+"_y")

        # if the folder is empty need to just write first line
        if j==0:
            writedat.write(newlable + "\n")
            writepoin.write(newpoinlable + "\n")
            for i in range(len(sol)):
                # add the first particles solution to the data file
                toadd = "%15.6f %15.6f %15.6f %15.6f"%(sol[i,0],sol[i,1],sol[i,2],sol[i,3])
                toadd += "\n"
                writedat.write(toadd)
                
                # this if statements cecks to see if we are at the right point to grab the
                # timesliced data from the solution. 
                # IN ORDER TO AVOID PROPAGATING ERROR THAT IS F-ING UP THE POINCARE SECTIONS WE NEED
                # TO .....
                # GOD DAMN IT OWEN

                # this should work better for poincare sections
                if ((l*dt)%(2*pl.pi/w)<dt):
                    # add first poincare section to poincare data file
                    writepoin.write(toadd)
                    l+=1
                    #checknextT += checkT 
                    # the +.5 efectivly rounds the number apropriatly
                    #checkPoint = int(checknextT/dt+.5)
                l+=1
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

                #if i==(checkPoint+newstart):
                # this should work better for poincare sections
                if ((l*dt)%(2*pl.pi/w)<dt):
                    print((l*dt)%(2*pl.pi/w))
                    writepoin.write(newline) 
                    #checknextT += checkT
                    # the +.5 efectivly rounds the number apropriatly
                    #checkPoint = int(checknextT/dt+.5)
                     
                l+=1
        # close the data file
        readdat.close()
        writedat.close()
        
        # close the poincare file
        readpoin.close()
        writepoin.close()


        # this is the feedback information to make the transitions smooth
        x0 = sol[-1,:]
        time = pl.arange(time[-1]+dt,time[-1]+totTime+dt,dt)

    os.remove("temp.txt")
    os.remove("pointemp.txt")
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


    

if __name__ == '__main__':
    main()
