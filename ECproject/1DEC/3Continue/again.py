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
import argparse
#import gc
#import weakref

def get_init_arr(lines):
    arr = pl.array([])
    count = -1
    ln_splited = lines[count].split()
    while ("ZONE" not in ln_splited[0]):
        for i,j in enumerate(ln_splited):
            arr = pl.append(arr,float(j))
        count -=1
        ln_splited = lines[count].split()

    arr = arr.reshape(abs(count)-1,4)

    return arr, (abs(count)-1)

def main():
    
    # before anything get us into the NormAll directory. this is the directory that will hold the
    # directories with the different data sets. We need to start keeping track of phase diagrams and
    # PC sections

    parser = argparse.ArgumentParser()

    parser.add_argument('--dir', action = 'store', dest = "dir",type = str,required = True)
    parser.add_argument('--file', action = 'store', dest = "file",type = str,required = True)

    parser.add_argument("--surf", action = "store", dest =  "surf", type = float,required = True) 
    parser.add_argument("--k"   , action = "store", dest =  "k"   , type = float,required = True) 
    parser.add_argument("--w"   , action = "store", dest =  "w"   , type = float,required = True)  
    parser.add_argument("--damp", action = "store", dest =  "damp", type = float,required = True) 
    parser.add_argument("--g"   , action = "store", dest =  "g"   , type = float,required = True)  
    parser.add_argument("--dt"  , action = "store", dest =  "dt"   , type = float,required = True)  

    inargs = parser.parse_args()

    dt = inargs.dt
    # total number of iterations to perform
    totIter = 40000
    totTime = totIter*dt
    time = pl.arange(0.0,totTime,dt)


    os.chdir(os.path.expanduser("~/Data/EC/2DBlock/Old/"+inargs.dir))

    # parameters for what should be chaotic but orderd trajectories
    surface = inargs.surf
    wave_num    = inargs.k   
    omega    = inargs.w   
    damping = inargs.damp
    grav    = inargs.g   


    # how many cells is till periodicity use x = n*pi/k (n must be even #)
    modNum = 2*pl.pi/wave_num
    
    # make ec object
    elc = ec.electricCurtain()
    
    
    # get all the initial conditions we age going to do
    cur_file = open(inargs.file,"r")
    cur_lines = cur_file.readlines()
    cur_file.close()

    cur_file = open(inargs.file,"a")

    # get coef from file
    coefficient = float(cur_lines[0].split()[-1])
    
    # get the new initial conditions as an array in the form
    # arr[i,j]
    # i denotes a particular particle
    # j denotes x,y,vx,vy 
    # j=0 ---> vx
    # j=1 ---> vy
    # j=2 ---> x
    # j=3 ---> y
    # thus function also returns the number of particles so we can loop over it
    init,p_num = get_init_arr(cur_lines)

    
    all_poin = pl.array([])

    # count the number of poin sections. Needd for reshaping all_poin array
    num_ts = 0
    print("pnum is : " + str(p_num))
    print("len(init) : " + str(len(init[:,0])))
    for i in xrange(p_num):

        apx = ec.surfCentreLineApx(coefficient,wave_num,omega,damping)

        # itial conditions to next point
        x0 = init[i,:]

        sol = odeint(apx.f,x0,time)

        
        for a in range(len(sol[:,0])):
            sol[a,2] = sol[a,2]%modNum
            if(((a*dt)%(2.0*pl.pi/omega))<dt):
                print(((a*dt)%(2.0*pl.pi/omega))-dt)
                all_poin = pl.append(all_poin,sol[a,:])
                
                if(i==0):
                    num_ts += 1

    # Now reshape all_poin and put it in the file corectly 
    print("p_num is: " +str(p_num))
    print("num_ts is: "+ str(num_ts))
    all_poin = all_poin.reshape(p_num,num_ts,4)

    # add the poin sections back to file
    for a in range(num_ts):

        cur_file.write("ZONE   I="+str(len(["used","to","be","filearr"]))+" DATAPACKING=POINT")
        cur_file.write("\n") 

        for b in range(p_num):
            # add the first particles solution to the data file
            toadd = "%15.6f %15.6f %15.6f %15.6f"%(all_poin[b,a,0],all_poin[b,a,1],all_poin[b,a,2],all_poin[b,a,3])
            toadd += "\n"
            cur_file.write(toadd)

           
    cur_file.close()

        ## make a file for the curent particles solution
        #curp += 1
        #curpstr = str(inargs.bnum) +"_"+ str(curp)

        ## keeptrack of the fiels for later
        #filearr = pl.append(filearr,curpstr)
        #
        #curpdatfile = open(curpstr,"a")
        #
        #apx = ec.surfCentreLineApx(coefficient,wave_num,omega,damping)

        ## itial conditions to next point
        #x0 = pl.array([initialvx+kapa*increment_vx,initialvy,initialx+alpha*increment_x,initialy])

        #sol = odeint(apx.f,x0,time)
        #
        #for a in range(len(sol[:,0])):
        #    sol[a,2] = sol[a,2]%modNum

        #for i in range(len(sol)):
        #    # add the first particles solution to the data file
        #    toadd = "%15.6f %15.6f %15.6f %15.6f"%(sol[i,0],sol[i,1],sol[i,2],sol[i,3])
        #    toadd += "\n"
        #    curpdatfile.write(toadd)

        #curpdatfile.close()
    
if __name__ == '__main__':
    main()
