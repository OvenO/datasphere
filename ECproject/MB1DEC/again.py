#!/usr/bin/python
#/opt/local/bin/python2.7

import numpy
import sys
sys.path.append("/users/o/m/omyers/datasphere/ECproject")
import ECclass as ec
import scipy as pl
import os
from scipy.integrate import odeint
import shutil
import time as thetime
import argparse
from datetime import datetime
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
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', action = 'store', dest = "dir",type = str,required = True)
    parser.add_argument('--file', action = 'store', dest = "file",type = int,required = True)
    parser.add_argument('--totiter', action = 'store', dest = "totIter",type = str, required = True)
    parser.add_argument('--sliced', action = 'store', dest = "sliced",type = str,required = True)
    inargs = parser.parse_args()
    dir = inargs.dir
    file = str(inargs.file)+'poindat.txt'
    if inargs.sliced == 'True':
        sliced = True
    elif inargs.sliced == 'False':
        sliced = False

    print('inargs.totItier: ' +str(inargs.totIter))
    print('inargs.dir: ' +str(inargs.dir))
    print('inargs.file: ' + str(inargs.file))
    print('inargs.sliced: '+str(inargs.sliced))
    totIter = float(inargs.totIter)

    # The reason we add in the inargs.file (which is just an iteger) is becase we need the directory
    # in the /tmp file to be compleately unique for each file. Becasue some files might be handeld
    # by the same node we need to distiguesh. I also want the file being delt with in a directory
    # becasue we also need an info file that goes with it -> if we have two different systems
    # running and the node is handeling both of them then they would end up with the same info file.
    # The same thing could happen with the poindat.txt files too but that is much less likely.
    print('ls /tmp: ' + str(os.listdir('/tmp')))
    print('ls /tmp/dir+inargs.file/: ' + str(os.listdir('/tmp/'+dir+str(inargs.file))))

    #info_file = open("/users/o/m/omyers/Data/EC/2DBlock/Old/"+dir+"/info.txt","r")
    # trying to use /tmp folder to speed things up
    info_file = open("/tmp/"+dir+str(inargs.file)+"/info.txt","r")
    lns = info_file.readlines()

    # Number of particles
    N = float(lns[3].split()[-1])
    # because we want to be able to sweep over many different variables we need to just reference
    # the general variable NOT being changed -> the one that is being changed is documented in the
    # first line of the actual data files.
    # we have aslo made it so the variable name and the number are on the same line. the split()[-1]
    dt = float(lns[2].split()[-1])
    # particle particle interaction stregth
    qq = float(lns[4].split()[-1])
    # number of unit cells till periodicity
    num_cell = float(lns[6].split()[-1])
    # how many cells is till periodicity use x = n*pi/k (n must be even #)
    modNum = 2*pl.pi*num_cell
    totTime = totIter*dt
    print("number of unit cells is: "+str(num_cell))
    print("dt is: " + str(dt))
    print("totTime is: "+str(totTime))

    time = pl.arange(0.0,totTime,dt)


    # try /tmp/ folder
    os.chdir("/tmp/"+dir+str(inargs.file))
    #os.chdir(os.path.expanduser("~/Data/EC/2DBlock/Old/"+dir))
    
    # open the file we are working with
    cur_file = open(file,"r")
    # get the line with the variable value we are doing the run with
    var_line = cur_file.readline()
    # Here we need to determine which it is so we can pass the right variables into the function
    if lns[1].split()[0] == "A:":
        As = pl.zerros(N) + float(lns[1].split()[-1])
        beta = float(var_lines[0].split()[-1])
    if lns[1].split()[0] == "beta:":
        beta = float(lns[1].split()[-1])
        As = pl.zerros(N) + float(var_lines[0].split()[-1])
    # get all the initial conditions we age going to do we will be using the formating that is:
    # 0-N -> velx  ex) with 10 particles arr[3] would be the 3rd particles velocity
    # N-2*N -> x
    init = pl.genfromtxt(file) 
    cur_file.close()

    # re open to append solution, we don't need to append, we just eed to write the whole soluton
    cur_file = open(file,"w")
    
    apx = ec.ExactPeriodic1D(qq,As,beta,num_cell)
    
    # this line makes it so that if we want to run again after we have already run then we can. This
    # grabs the last line to use as the first initial conditions
    x0 = init[-1,:]
    
    sol = odeint(apx.f,x0,time)

    # modulus the x positions so they are in the right spot
    sol[:,N:2*N] = sol[:,N:2*N]%modNum
    
    if sliced:
        for a in range(len(sol[:,0])):
            if(((a*dt)%(2.0*pl.pi))<dt):
                # str(sol)[1,-1] gets rid of the '[' and ']' at the begining and end of the array
                # when it is turned into a string.
                cur_file.write(str(sol[a,:])[1,-1])
    else:
        for a in range(len(sol[:,0])):
            cur_file.write(str(sol[a,:])[1,-1])


        
if __name__ == '__main__':
    main()
