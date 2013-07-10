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


def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', action = 'store', dest = "dir",type = str,required = True)
    parser.add_argument('-n', action = 'store', dest = "n",type = int,required = True)
    inargs = parser.parse_args()

    os.chdir(os.path.expanduser("~/Data/EC/2DBlock/Old/"+inargs.dir))
    os.mkdir("Last")

    os.system("cp info.txt Last/")
    
    all = os.listdir(".")
    for i,j in enumerate(all):
        if ("poindat" in j):
            curfile = open(j,"r")
            to_file = open("Last/"+j,"w")
            lines = curfile.readlines()
            to_file.write(lines[0])
            to_file.write(lines[1])
            count = 0
            for a in range(inargs.n):
                count-=1
                l_split = lines[count].split() 
                while("ZONE" not in l_split[0]):
                    to_file.write(lines[count])
                    count-=1
                    l_split = lines[count].split()
            to_file.close()
            curfile.close()
   
if __name__ == '__main__':
    main()
