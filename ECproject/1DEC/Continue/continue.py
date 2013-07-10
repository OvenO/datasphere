import os
import argparse
import time as thetime


# this function needs the number of small blocks passed into it becuase
# it is going to see see if the number of files we get back from the
# qsub is the same as the number we submiter ... which is num_small_blocks
def are_we_done(num_small_blocks):

    all_files = os.listdir(".")
    
    count = 0

    for a,b in enumerate(all_files):
        if ("cont_wh" in b):
            count+=1
    
    if (count == num_small_blocks):
        return True
    else:
        return False

#******************************************************************* 
#******************************************************************* 
def info_to_arg_str(dir):
    print("anyting")    
    arg_str = " "
    info_file = open("/users/o/m/omyers/Data/EC/2DBlock/Old/"+dir+"/info.txt","r")
    lns = info_file.readlines()

    for i,j in enumerate(lns):
        arg_str += j[:-1] + " "

    return arg_str

def main():
    count_the_time = 0.0
    print("in main")

    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', action = 'store', dest = "dir",type = str,required = True)


    # star num and end num alow us to just work on certian files if most of them look pretty good.
    # or you can do all of them at once. you have some choice.
    # snum --> start number
    # enum --> end number
    parser.add_argument('--snum', action = 'store', dest = "snum",type = int,default = 0)
    parser.add_argument('--enum', action = 'store', dest = "enum",type = int,default = 0)

    # new way... now we are set up so it will just loop over the whole file in sections. The number
    # of sections is sub_times and this number will come as an argument

    parser.add_argument('--stimes', action = 'store', dest = "stimes",type = int,default=1)
    parser.add_argument('--aonce', action = 'store', dest = "aonce",type = int,default=0)
    # incase sections start at non zero number like 200poindat.txt
    parser.add_argument('--start', action = 'store', dest = 'start',type = int,default=0)
    inargs = parser.parse_args()

    dir = inargs.dir
    snum = inargs.snum
    enum = inargs.enum
    
    start = inargs.start
    
    # number of times to submit 
    sub_times = inargs.stimes
    # how many are done at a submission
    aonce = inargs.aonce
    
    # keep program flexable so we can use --snum --enum still if we want
    if (enum == 0):
        enum = aonce
            
    print("what")
    arg_str = info_to_arg_str(dir)
    print(arg_str)


    for alpha in range(sub_times):
        print(snum+alpha*aonce,enum+alpha*aonce)
        for i in range(snum+alpha*aonce+start,enum+alpha*aonce+start):
        
            # Tell again.py what file to work on 
            arg_str += " --file " + str(i) + "poindat.txt"
        
            # write the write script file so the arguments are in
            # the python line of the script file. then this line
            os.system("cp again.script again_temp.script")
            edit_file = open("again_temp.script","a")
            edit_file.write("python again.py "+arg_str)
            edit_file.close()
            os.system("qsub again_temp.script")
            os.remove("again_temp.script") 
        
    
        # now we are going to wait for all of the little block to be done and then agrogate the
        # files into data.txt and poin.txt
        done = False
        
        cnt = 0
        while (not done):
            cnt += 1
            thetime.sleep(10) 

            count_the_time += 10.0
            print(count_the_time)

            # are_we_done function checks to see if all the return files are back from q sub. if they
            # are all there... We are done :)
            done = are_we_done(enum-snum)
    
            # kill it if we have been going for too long (means something is probably wrong) 
            if cnt>100000:
                done = True
            
    
        # for the are_we_done fuction to work again the next time through the loop we need to
        # remove the return fiels
        os.system("rm cont_wh*")


if __name__ == "__main__":
    main()
