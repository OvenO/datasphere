import sys
sys.path.append("/users/o/m/omyers/puthere/lib/python2.7/site-packages")
sys.path.append("/users/o/m/omyers/puthere/matplotlib/build/lib.linux-x86_64-2.7")
import pylab as pl
import argparse
import os


def main(): 
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', action = 'store', dest = "dir",type = str,required = True)
    parser.add_argument('--dt', action = 'store', dest = "dt",type = float,required = True)
    parser.add_argument('-w', action = 'store', dest = 'w', type = float, required = True)

    inargs = parser.parse_args()

    dt = inargs.dt
    w = inargs.w

    os.chdir("/users/o/m/omyers/Data/EC/4DBlock/"+inargs.dir)

    # find the number of total iterationos
    # open a file and seeing how long it is
    temp_file = open("0_1","r")
    temp_lines = temp_file.readlines()
    tot_iter = len(temp_lines)
    temp_file.close()
    
    # want the file array before we make the data files
    filearr = os.listdir(".")
    
    datafile = open("data.txt","w")
    
    poinfile = open("poindat.txt","w")

#*******************************************************************************************
#*******************************************************************************************
    # Im going to comment out all the data file stuff to see how much things speed up when we are
    # just doing the poincare sections.
#*******************************************************************************************
#*******************************************************************************************

    for i in xrange(tot_iter):
        # DATAPACKING=POINT should mean that the format is as such:
        # vx   vy   x    y
        # #'s  #'s  #'s  #'s
        # #'s  #'s  #'s  #'s
        # #'s  #'s  #'s  #'s
        # #'s  #'s  #'s  #'s
#        datafile.write("ZONE   I="+str(len(filearr))+" DATAPACKING=POINT")
#        datafile.write("\n")

        if (((i*dt)%(2*pl.pi/w))<dt):
            poinfile.write("ZONE   I="+str(len(filearr))+" DATAPACKING=POINT")
            poinfile.write("\n")


        for a,b in enumerate(filearr):
            curfile = open(b,"r")
            lines = curfile.readlines()
            tofile = lines[i]
#            datafile.write(tofile)

            if (((i*dt)%(2*pl.pi/w))<dt):

                poinfile.write(tofile)

            curfile.close()


    datafile.close()
    poinfile.close()
        
    os.system("rm 0* 1* 2* 3* 4* 5* 6* 7* 8* 9*")

if __name__ == "__main__":
    main()
