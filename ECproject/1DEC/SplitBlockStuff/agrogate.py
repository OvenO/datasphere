import sys
sys.path.append("/users/o/m/omyers/puthere/lib/python2.7/site-packages")
sys.path.append("/users/o/m/omyers/puthere/matplotlib/build/lib.linux-x86_64-2.7")
import pylab as pl
import argparse
import os


def main(): 
    # this is so we can make a progreesion of the right file names
    add_num_file_name = 0
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', action = 'store', dest = "dir",type = str,required = True)
    parser.add_argument('--dt', action = 'store', dest = "dt",type = float,required = True)
    parser.add_argument('-w', action = 'store', dest = 'w', type = float, required = True)
    parser.add_argument('--bnum', action = 'store', dest = 'bnum', type = int, required = True)
    parser.add_argument('--pts', action = 'store', dest = 'pts', type = int, required = True)
    parser.add_argument('--stimes', action = 'store', dest = 'stimes', type = int, required = True)
    parser.add_argument('--icount', action = 'store', dest = 'icount', type = int, required = True)
    parser.add_argument('--coef', action = 'store', dest = 'coef', type = float, required = True)

    inargs = parser.parse_args()

    dt = inargs.dt
    w = inargs.w
    bnum = inargs.bnum
    stimes = inargs.stimes
    init_count = inargs.icount
    dir = inargs.dir
    coef = inargs.coef

    os.chdir("/users/o/m/omyers/Data/EC/2DBlock/"+dir)

    # find the number of total iterationos
    # open a file and seeing how long it is
    temp_file = open(str(init_count*stimes) + "_1","r")
    temp_lines = temp_file.readlines()
    tot_iter = len(temp_lines)
    temp_file.close()
    
    # init_count comes from the number of times we are increasing the coefiissient for a giving
    # subission. Because we are limmiting our sumbisison this number will only by areound 10 or 20.
    # the purpose of the set up of this loop is to only agrogate the fiels that were most recently
    # compleated so we arent doing all of them all at once. 
    print(range(init_count*stimes,bnum))
    print("bnum is: " + str(bnum))
    print("init_count is: " + str(init_count))
    print("stimes is: " + str(stimes))
    print("tot_iter is: " + str(tot_iter))
    for alpha in range(init_count*stimes,bnum):
        print(alpha)
    
        poinfile = open(str(alpha+add_num_file_name)+"poindat.txt","w")
        poinfile.write("coef --> " + str(coef) + "\n")

        for i in xrange(tot_iter):
            # DATAPACKING=POINT should mean that the format is as such:
            # vx   vy   x    y
            # #'s  #'s  #'s  #'s
            # #'s  #'s  #'s  #'s
            # #'s  #'s  #'s  #'s
            # #'s  #'s  #'s  #'s

            if (((i*dt)%(2*pl.pi/w))<dt):
                poinfile.write("ZONE   I="+str(len(["used","to","be","filearr"]))+" DATAPACKING=POINT")
                poinfile.write("\n")
            
            print ("inargs.pts is: " + str(inargs.pts))

            for b in xrange(inargs.pts-1):
                file_str = str(alpha) + "_" + str(b+1)

                curfile = open(file_str,"r")
                lines = curfile.readlines()
                tofile = lines[i]

                if (((i*dt)%(2*pl.pi/w))<dt):

                    poinfile.write(tofile)

                    print("writing to poinfile")

                curfile.close()

        poinfile.close()
    

    for alpha in range(init_count*stimes,bnum):
        # The _* is nessasary so we don't accedentaly delete our poindat.txt files that have a #
        # heading of alpha as well
        os.system("rm " + str(alpha) + "_*")

if __name__ == "__main__":
    main()
