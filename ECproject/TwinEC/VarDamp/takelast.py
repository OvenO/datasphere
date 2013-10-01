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

def get_b_num():
    f = open("0poindat.txt")
    lines = f.readlines()
    count = 0
    for i,j in enumerate(lines[2:]):
        if "ZONE" in j:
            break
        count+=1
    return count
def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', action = 'store', dest = "dir",type = str,required = True)
    inargs = parser.parse_args()

    os.chdir(os.path.expanduser("~/Data/EC/4DBlock/Old/"+inargs.dir))

    os.mkdir("Last")
 
    b_num = get_b_num()

    os.system("cp info.txt Last/")
    
    all = os.listdir(".")
    for i,j in enumerate(all):
        if ("poindat" in j):
            curfile = open(j,"r")
            to_file = open("Last/"+j,"w")
            to_file.write(curfile.readline())
            to_file.write(curfile.readline())
            cur_dat = pl.genfromtxt(curfile,comments="Z")
            for l in range(b_num):
                to_file.write(str(cur_dat[-b_num+l,0])+" "+str(cur_dat[-b_num+l,1])+" "+str(cur_dat[-b_num+l,2])+" "+str(cur_dat[-b_num+l,3]))
                to_file.write("\n")
            to_file.close()
            curfile.close()
   
if __name__ == '__main__':
    main()
