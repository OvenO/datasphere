import pylab as pl
import os

def main():
    #in NormAll directroy which file are the data files in
    dataFile="Data3"

    #how many iterations into future
    #for poincar version
    numpoincar = 6
    #for normal version
    jump=20

    #number of plots
    numplot=4



    os.chdir(os.path.expanduser("~/Data/PP/NormAllDat/"+dataFile))
    
    files = os.listdir(".")

    # make different arrays for different data
    dat=pl.array([])
    poincardat=pl.array([])

    # seperate poincare section files and full trajectory files
    for i,j in enumerate(files):
        #print(j[-16:])
        if(j[-16:]=="poincaredata.txt"):
            poincardat=pl.append(poincardat,j)
        elif(j!="info.dat"):
                dat=pl.append(dat,j)
    
    
    for a in range(int(10000/20)-1):
        # Here we will just try to plot a block runing through ONE poincare section
        #make arrays to story thee block points to be ploted
        firstBptsx=pl.array([])
        firstBptsvx=pl.array([])
        secondBptsx=pl.array([])
        secondBptsvx=pl.array([])


        for i,j in enumerate(dat):

            print(j)
            curfile = open(j,"r")
            curlines = curfile.readlines()


            #spilit the position velocity and time from the line
            firstline=curlines[1]
            nextline=curlines[(a+1)*jump]

            #split the lines to parse the information
            firstx = float(firstline.split(",")[0])
            firstvx = float(firstline.split(",")[1])

            nextx = float(nextline.split(",")[0])
            nextvx = float(nextline.split(",")[1])

            firstBptsx = pl.append(firstBptsx,firstx)
            firstBptsvx = pl.append(firstBptsvx,firstvx)

            secondBptsx = pl.append(secondBptsx,nextx)
            secondBptsvx  = pl.append(secondBptsvx,nextvx)

            
        fig = pl.figure()
        ax = fig.add_subplot(111)
        ax.scatter(firstBptsx,firstBptsvx)
        ax.scatter(secondBptsx,secondBptsvx)
        fig.savefig(str(a)+"blocks.png")



#    # Here we will just try to plot a block runing through ONE poincare section
#    #make arrays to story thee block points to be ploted
#    firstBptsx=pl.array([])
#    firstBptsvx=pl.array([])
#    secondBptsx=pl.array([])
#    secondBptsvx=pl.array([])
#
#    
#    for i,j in enumerate(poincardat):
#        curfile = open(j,"r")
#        curlines = curfile.readlines()
#
#        print(j)
#
#        #spilit the position velocity and time from the line
#        firstline=curlines[1]
#        nextline=curlines[numpoincar]
#
#        #split the lines to parse the information
#        firstx = float(firstline.split(",")[0])
#        firstvx = float(firstline.split(",")[1])
#
#        nextx = float(nextline.split(",")[0])
#        nextvx = float(nextline.split(",")[1])
#
#        firstBptsx = pl.append(firstBptsx,firstx)
#        firstBptsvx = pl.append(firstBptsvx,firstvx)
#
#        secondBptsx = pl.append(secondBptsx,nextx)
#        secondBptsvx  = pl.append(secondBptsvx,nextvx)

        
    

    
    os.chdir("../..")
if __name__=="__main__":
    main()
