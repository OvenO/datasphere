import time as thetime
import os

# this function needs the number of small blocks passed into it becuase
# it is going to see see if the number of files we get back from the
# qsub is the same as the number we submiter ... which is num_small_blocks
def are_we_done(num):

    all_files = os.listdir(".")
    
    count = 0

    for a,b in enumerate(all_files):
        if ("again.o" in b):
            count+=1
    
    if (count == num):
        return True
    else:
        return False

        
def main():    
    # now we are going to wait for all of the little block to be done and then agrogate the
    # files into data.txt and poin.txt
    count_the_time=0.0
    done = False
    
    cnt = 0
    while (not done):
        cnt += 1
        thetime.sleep(10) 
    
        count_the_time += 10.0
        print(count_the_time)
    
        # are_we_done function checks to see if all the return files are back from q sub. if they
        # are all there... We are done :)
        done = are_we_done(500)
    
        # kill it if we have been going for too long (means something is probably wrong) 
        if cnt>100000:
            done = True
        
    
    # for the are_we_done fuction to work again the next time through the loop we need to
    # remove the return fiels
    #os.system("rm again*")
    print("this worked")
if __name__ == "__main__":
    main()
