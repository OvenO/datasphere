import numpy
import sys
sys.path.append("/users/o/m/omyers/puthere/lib/python2.7/site-packages")
sys.path.append("/users/o/m/omyers/datasphere/ECproject")
sys.path.append("/users/o/m/omyers/puthere/matplotlib/build/lib.linux-x86_64-2.7")
import ECclass as ec
import pylab as pl
import os
import time as thetime

# Here we are going to just make a Directory with an info file and 1 Block of initial contitions so
# we can run continue on it. Should be faster than hm_bif method

#******************************************************************* 
#******************************************************************* 
def get_make_dir():

    # use time to name directories
    timestr = thetime.asctime()
    newtimestr = "Block_"
    for a,b in enumerate(timestr.split()):
        newtimestr+=b+"_"

    newtimestr = newtimestr.replace(":","") 
    
    print(newtimestr)
    os.system("mkdir /users/o/m/omyers/Data/EC/2DBlock/"+newtimestr)
    
    return newtimestr

#******************************************************************* 
#******************************************************************* 
def main():

    how_many = 1
    # increase coef each time by this amount
    inc_coef = .001
    
    startnum = 0

    #EC parameters
    surf =  1.0 
    # good place to start
    #coef =  .2004
    coef = .2500
    k    =  1.0 
    w    =  1.0 
    damp =  .1
    g    =  .1
    dt   =  .05

    # define the lower left corner of block
    initialvx = -1.5
    initialvy = 0.0
    initialx  = 0.0
    initialy  = 1.0
    # define dimensions of block
    xby = 6.28
    vxby = 3.0
    # define number of points in each direction
    num_of_x  = 120.0
    num_of_vx = 60.0
    #num_of_x  = 60.0
    #num_of_vx = 30.0
    # distance between points
    increment_x  = xby/num_of_x
    increment_vx = vxby/num_of_vx

    dir = get_make_dir()

    # make info file
    info_file = open("/users/o/m/omyers/Data/EC/2DBlock/"+dir+"/info.txt","w")
    info_file.write("--dir\n"+str(dir)+"\n--surf\n"+str(surf)+"\n--k\n"+str(k)+"\n--w\n"+str(w)+\
            "\n--damp\n"+str(damp)+"\n--g\n"+str(g)+"\n--dt\n"+str(dt)+" \n")
    info_file.close()

    for l in range(how_many):
        # make the file
        file = open("/users/o/m/omyers/Data/EC/2DBlock/"+dir+"/" +str(startnum+l)+"poindat.txt","w")

        # write the first few lines of the file
        file.write("coef --> " +str(coef)+"\n")
        file.write("ZONE I=something DATAPACKING=POINT\n")

        for kapa in range (int(num_of_vx)):
            for alpha in range(int(num_of_x)):
                # write vx 
                file.write(str(initialvx + kapa*increment_vx) + "  ")
                # write vy
                file.write(str(initialvy) + "  ")
                # write x
                file.write(str(initialx + alpha*increment_x) + "  ")
                # write y
                file.write(str(initialy)) 
                # write a enter
                file.write("\n")

        
        # now if we want we can coppy this file to make 1pointdat.txt, 2poindat.txt, ect 
        # but we just have to make sure we put the right coeficient as the first line
        file.close()
        coef += inc_coef

if __name__ == "__main__":
    main()
            
 

