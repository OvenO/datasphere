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

    #start_time=datetime.now() 

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
    print('info lines: '+str(lns))
    for i,j in enumerate(lns):
        if 'sweep var' in j:
            print('found sweep variable in info file')
            var_str = j.split()[0]
            print('var_str is: '+str(var_str))
    info_file.close()

    # open the file we are working with
    cur_file = open("/tmp/"+dir+str(inargs.file)+"/"+file,"r")
    # get the line with the variable value we are doing the run with
    var_line = cur_file.readline()
    print('var_line is: ' + var_line)

    # because we want to be able to sweep over many different variables we need to just reference
    # the general variable NOT being changed -> the one that is being changed is documented in the
    # first line of the actual data files.
    if var_str!='dt':
        dt= float(lns[1].split()[-1])
        print('dt: '+str(dt))
    else:
        dt = float(var_line.split()[-1])

    if var_str!='cycles':
        cycles = float(lns[2].split()[-1])
        print('cycles: '+str(cycles))
    else:
        cycles = float(var_line.split()[-1])

    if var_str!='N':
        N = int(lns[3].split()[-1])
        print('N: '+str(N))
    else:
        N = int(var_line.split()[-1])

    # particle particle interaction stregth
    if var_str!='B':
        B = float(lns[4].split()[-1])
        print('B: '+str(B))
    else:
        B = float(var_line.split()[-1])

    if var_str!='O_mega':
        O_mega = float(lns[5].split()[-1])
        print('O_mega: ' + str(O_mega))
    else:
        O_mega = float(var_line.split()[-1])

    if var_str!='A':
        A = float(lns[6].split()[-1])
        print('A: '+str(A))
    else:
        A = float(var_line.split()[-1])

    totTime = totIter*dt
    print("dt is: " + str(dt))
    print("totTime is: "+str(totTime))

    time = pl.arange(0.0,totTime,dt)

    # try /tmp/ folder
    os.chdir("/tmp/"+dir+str(inargs.file))
    #os.chdir(os.path.expanduser("~/Data/EC/2DBlock/Old/"+dir))
    
    # get all the initial conditions we age going to do we will be using the formating that is:
    # 0-N -> velx  ex) with 10 particles arr[3] would be the 3rd particles velocity
    # N-2*N -> x
    init = pl.genfromtxt(cur_file) 
    cur_file.close()

    # re open to append solution, we don't need to append, we just eed to write the whole soluton
    cur_file = open(file,"w")
    cur_file.write(var_line)
    
    apx = ec.OurHMF(N,O_mega,A,B)
    
    # this statment makes it so that if we want to run again after we have already run then we can. This
    # grabs the last line to use as the first initial conditions
    if len(pl.shape(init))>1: 
        x0 = init[-1,:]
    else: x0 = init
    
    print('about to run odeint')
    sol = odeint(apx.f,x0,time)
    print('ran. shape(sol): ' + str(pl.shape(sol)))
    print('first 10 lines of sol:')
    for i in range(10):
        print(str(sol[i,:]))

    if sliced:
        for a in range(len(sol[:,0])):
            if(((a*dt)%(2.0*pl.pi))<dt):
                # str(sol)[1,-1] gets rid of the '[' and ']' at the begining and end of the array
                # when it is turned into a string.
                cur_file.write(str(sol[a,:])[1:-1].replace('\n',''))
                cur_file.write('\n')
    else:
        for a in range(len(sol[:,0])):
            cur_file.write(str(sol[a,:])[1:-1].replace('\n',''))
            cur_file.write('\n')

    cur_file.close()


    #print (datetime.now() - start_time)

    os.system('gzip /tmp/'+dir+str(inargs.file)+'/'+file)
    os.system('cp /tmp/'+dir+str(inargs.file)+'/'+file+'.gz /users/o/m/omyers/Data/EC/2DBlock/Old/'+dir+'/'+file+'.gz')
    os.system('rm -r /tmp/'+dir+str(inargs.file))
    
        
if __name__ == '__main__':
    main()
