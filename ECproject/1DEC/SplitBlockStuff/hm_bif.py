import numpy
import sys
sys.path.append("/users/o/m/omyers/puthere/lib/python2.7/site-packages")
sys.path.append("/users/o/m/omyers/datasphere/ECproject")
sys.path.append("/users/o/m/omyers/puthere/matplotlib/build/lib.linux-x86_64-2.7")
import ECclass as ec
import pylab as pl
import os
import time as thetime

#******************************************************************* 
#******************************************************************* 
def agrogate(directory,dt,w,num_small_blocks,pts,times_subed,num_inc,coef):
    

    arg_str = "--dir " + directory + " -w " + str(w) + " --dt " + str(dt) + " --bnum " + \
            str(num_small_blocks) + " --pts " + str(pts) + " --stimes " + str(times_subed) + \
            " --icount " + str(num_inc) + " --coef " + str(coef)
    
    # write the write script file so the arguments are in
    # the python line of the script file. then this line
    os.system("cp agrogate.script ag_temp.script")
    edit_file = open("ag_temp.script","a")
    edit_file.write("python agrogate.py "+arg_str)
    print("dirictory is: "+ directory)
    edit_file.close()
    os.system("qsub ag_temp.script")
    os.remove("ag_temp.script") 
    

#******************************************************************* 
#******************************************************************* 
# this function needs the number of small blocks passed into it becuase
# it is going to see see if the number of files we get back from the
# qsub is the same as the number we submiter ... which is num_small_blocks
def are_we_done(num_small_blocks,directory):
    all_files = os.listdir(".")
    
    count = 0

    for a,b in enumerate(all_files):
        if (directory in b):
            count+=1
    
    if (count == num_small_blocks):
        return True
    else:
        return False

#******************************************************************* 
#******************************************************************* 
def are_fin_done():
    all_files = os.listdir(".")
    
    count = 0

    for a,b in enumerate(all_files):
        if ("finup" in b):
            count+=1
    
    if (count == 1):
        return True
    else:
        return False


#******************************************************************* 
#******************************************************************* 
def get_make_dir():

    # use time to name directories
    timestr = thetime.asctime()
    newtimestr = "Block_"
    for a,b in enumerate(timestr.split()):
        newtimestr+=b+"_"

    newtimestr = newtimestr.replace(":","") 
    
    os.system("mkdir /users/o/m/omyers/Data/EC/2DBlock/"+newtimestr)
    
    return newtimestr

#******************************************************************* 
#******************************************************************* 
def run_the_block(arg_str,count,coef,directory):

        arg_str += "--coef " + str(coef) + " "
        # we also need to pass in a number that uniquely defineds each little block we are
        # submiting. this is going to head the file names that get put in a single directory
        arg_str += "--bnum " +str(count)+" "

        # clean up the white space
        arg_str_arr = arg_str.split()
        arg_str = ""
        for alpha,beta in enumerate(arg_str_arr):
            arg_str += beta
            arg_str += " "
        
        # write the write script file so the arguments are in
        # the python line of the script file. then this line
        
        # origonal file
        of = open("tec2DB_splitblock.script","r")
        o_lines = of.readlines()
        # new file
        nf = open("temp.script","w")
        for i,j in enumerate(o_lines):
            if ("ECtec2DB_split" in j):
                nf.write("#PBS -N ")
                nf.write(directory)
                nf.write("\n")
            else:
                nf.write(j)
                nf.write("\n")


        nf.write("python tec2DB_splitblock.py "+arg_str)
        nf.close()
        os.system("qsub temp.script")
        os.remove("temp.script") 

#******************************************************************* 
#******************************************************************* 
#******************************************************************* 
def main():

    #EC parameters
    surf =  1.0 

    # good place to start
    #coef =  .2004
    coef = .000 +.0002*000
    k    =  1.0 
    w    =  1.0 
    damp =  .1 
    g    =  .1
    dt   =  .05

    # increase coef each time by this amount
    inc_coef = .0005

    # increase the coeficient this many times
    # should keep less than 20 so only 20 jobs are running at a time
    #num_inc = 10

    num_inc = 10

    # number of time we are going to submit num_inc so we can do may runs but only num_inc at a time
    sub_times = 20

    # need a folder for ALLLLLL of the data
    directory = get_make_dir()
    print("directory is: " + directory)
    
    # keep track of the total number of times we increase out parameter
    count = 0

    arg_str = "--dir " + directory + " "

    arg_str += "--surf " + str(surf) \
        + " --k " + str(k) \
        + " --w " + str(w) \
        + " --damp " + str(damp) \
        + " --g " + str(g) \
        + " --dt " + str(dt) + " "

    for beta in range(sub_times):
        for b in range(num_inc):

            run_the_block(arg_str,count,coef,directory)

            coef += inc_coef
            count += 1

                
        # now we are going to wait for all of the little block to be done and then agrogate the
        # files into data.txt and poin.txt
        done = False
        
        cnt = 0
        while (not done):
            cnt += 1
            thetime.sleep(10) 

            # are_we_done function checks to see if all the return files are back from q sub. if they
            # are all there... We are done :)
            done = are_we_done(num_inc,directory)

            # kill it if we have been going for too long (means something is probably wrong) 
            if cnt>100000:
                done = True
            

        # for the are_we_done fuction to work again the next time through the loop we need to
        # remove the return fiels
        os.system("rm Block_*")
        
        # agrogate takes all the individual files an puts them into data.txt and poindat.txt
        # Pass is the directory (string) that we are putting everything in so we can operate on the
        # files there.
        # Also it needs dt and w to make poindat.txt file corectly
        agrogate(directory,dt,w,count,900,beta,num_inc,coef)

       
        # This code here will keep the program from moving on until agrogate is compleately done
        #fin_done = False
        #while (not fin_done):
        #    thetime.sleep(30)
        #    fin_done = are_fin_done()

    # make an info file
    info_file = open("/users/o/m/omyers/Data/EC/2DBlock/"+ directory+"/info.txt","w")
    splited_arg = arg_str.split()
    for l,k in enumerate(splited_arg):
        info_file.write(k)    
        info_file.write(" \n ")
    info_file.close()

if __name__ == "__main__":
    main()
