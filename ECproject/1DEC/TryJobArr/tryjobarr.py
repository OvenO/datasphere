import numpy
import sys
sys.path.append("/users/o/m/omyers/puthere/lib/python2.7/site-packages")
sys.path.append("/users/o/m/omyers/datasphere/ECproject")
sys.path.append("/users/o/m/omyers/puthere/matplotlib/build/lib.linux-x86_64-2.7")
import ECclass as ec
import pylab as pl
import os
import time as thetime
import argparse

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', action = 'store', dest = "dir",type = str,required = True)

    inargs = parser.parse_args()

    dir = inargs.dir

    
    file = open(dir,"w")
    file.write("Hi! I'm a File")
    file.close()



if __name__ == "__main__":
    main()
            
 

