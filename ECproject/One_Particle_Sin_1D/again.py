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
    # keep track of runtime
    start_time = datetime.now()
    
    # before anything get us into the NormAll directory. this is the directory that will hold the
    # directories with the different data sets. We need to start keeping track of phase diagrams and
    # PC sections

    parser = argparse.ArgumentParser()

    parser.add_argument('--dir', action = 'store', dest = "dir",type = str,required = True)
    parser.add_argument('--file', action = 'store', dest = "file",type = int,required = True)
    parser.add_argument('--totiter', action = 'store', dest = "totIter",type = str, required = True)
    parser.add_argument('--sliced', action = 'store', dest = "sliced",type = str,required = True)
    # if this is true we are going to just save data for the last cycle of the run
    parser.add_argument('--fulast', action = 'store', dest = "fulast",type = bool,required = True)
    inargs = parser.parse_args()

    full_last = inargs.fulast
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

    for i,j in enumerate(lns):
        if 'dt' in j:
            dt = float(j.split()[-1])
        if 'beta' in j:
            beta = float(j.split()[-1])

    totTime = totIter*dt
    print("dt is: " + str(dt))
    print("totTime is: "+str(totTime))
    time = pl.arange(0.0,totTime,dt)


    # try /tmp/ folder
    os.chdir("/tmp/"+dir+str(inargs.file))
    #os.chdir(os.path.expanduser("~/Data/EC/2DBlock/Old/"+dir))

    # how many cells is till periodicity use x = n*pi/k (n must be even #)
    modNum = 2*pl.pi
    
    # get all the initial conditions we age going to do
    cur_file = open(file,"r")
    cur_lines = cur_file.readlines()
    cur_file.close()

    cur_file = open(file,"a")

    # get coef from file
    A = float(cur_lines[0].split()[-1])
    
    init,p_num = get_init_arr(cur_lines)

    all_poin = pl.array([])

    # count the number of poin sections. Needd for reshaping all_poin array
    num_ts = 0
    #print("pnum is : " + str(p_num))
    #print("len(init) : " + str(len(init[:,0])))
    for i in xrange(p_num):

        apx = ec.One_Particle_Ensble_Sin1D(A,beta)
        # itial conditions to next point
        x0 = init[i,:]

        sol = odeint(apx.f,x0,time)

        if sliced:
            for a in range(len(sol[:,0])):
                sol[a,2] = sol[a,2]%modNum
                if(((a*dt)%(2.0*pl.pi))<dt):
                    all_poin = pl.append(all_poin,sol[a,:]) 
                    if(i==0): 
                        num_ts += 1 
        else: 
            all_poin = pl.append(all_poin,sol)
            num_ts = len(time)

    # Now reshape all_poin and put it in the file corectly 
    print("p_num is: " +str(p_num))
    print("num_ts is: "+ str(num_ts))
    all_poin = all_poin.reshape(p_num,-1,4)

    # add the poin sections back to file. Starting at 1 (if not full_last) because we dont need to
    # repeat the PC section that is already there. if full_last is true then Sliced must be false
    # and we only want the last 2*pl.pi time.
    if full_last:
        begin_write = num_ts - int(2.0*pl.pi/dt)
    else:
        begin_write = 1

    for a in range(begin_write,num_ts):

        cur_file.write("ZONE   I="+str(len(["used","to","be","filearr"]))+" DATAPACKING=POINT")
        cur_file.write("\n") 
        
        # I think we are reversing the order of the points in a PC section so I'm going to try to
        # write the file backwards (this the strange indexing and the -b).
        for b in range(1,p_num+1):
            # add the first particles solution to the data file
            toadd = "%15.6f %15.6f %15.6f %15.6f"%(all_poin[-b,a,0],all_poin[-b,a,1],all_poin[-b,a,2],all_poin[-b,a,3])
            toadd += "\n"
            cur_file.write(toadd)
           
    cur_file.close()

    print (datetime.now() - start_time)
    
    os.system('gzip /tmp/'+dir+str(inargs.file)+'/'+file)
    os.system('cp /tmp/'+dir+str(inargs.file)+'/'+file+'.gz /users/o/m/omyers/Data/EC/2DBlock/Old/'+dir+'/'+file+'.gz')
    os.system('rm -r /tmp/'+dir+str(inargs.file))
    
    #os.system('cp -r /tmp/'+dir +' ' + '/users/o/m/omyers/Data/EC/2DBlock/Old/'+dir)
    #os.system('rm -r /tmp/'+dir)
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
