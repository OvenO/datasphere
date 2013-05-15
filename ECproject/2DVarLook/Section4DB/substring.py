import os
import shutil
import time

def main():

# define lower left of block as corner # 1
c1 = [,]





# Old Old Old Old Old Old Old Old Old Old Old Old Old Old Old Old Old Old Old Old Old Old Old Old Old Old Old Old Old Old Old Old
    #initial_g = .1
    #increace_g = .05
    #num_g = 5

    initial_b = .1
    increace_b = .05
    num = 4

    
    
    for i in range(num):

        # fix the time string so spaces are underscors
        thetime = time.asctime()
        newtime =""
        for l,m in enumerate(thetime.split()):
            newtime += m + "_"
        newtime = newtime[:-1]
        # what will the new g value be?
        argstr = "-b "+ str(initial_b + i*increace_b)+ " -d "+newtime
        shutil.copyfile("substring.script","temp.script")
        editfile = open("temp.script","a")
        #shutil.copyfile("substring.script","varcoef_g"+str(g)+"_b"+str(b)+".script")
        #editfile = open("substring.script","varcoef_g"+str(g)+"_b"+str(b)+".script","a")
        editfile.write("python 2Dlookvarcoef.py "+argstr)
        editfile.close()

        os.system("qsub temp.script")
    
        os.remove("temp.script")
        
        time.sleep(5)

if __name__=="__main__":
    main()
