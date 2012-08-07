#import pylab as pl
import sys
sys.path.append("/users/o/m/omyers/puthere/lib/python2.7/site-packages")
import scipy as pl
import os

def whatismaxdist(compx,compvx):
    # array of distances
    distarr = pl.array([])
    
    for a in range(4):
        for i in range(a+1,5)
        distarr = pl.append(distarr,pl.sqrt((compx[a]-compx[i])**2+(compvx[a]-compvx[i])**2))
    
    return distarr

def main():
    #in NormAll directroy which file are the data files in
    dataFile="Data3"



    os.chdir(os.path.expanduser("~/Data/EC/NormAllDat/"+dataFile))
    
    files = os.listdir(".")

    # make different arrays for different data
    dat=pl.array([])
    poincardat=pl.array([])

    #---------------------------------------------------------------------------------------------
    # FOR NOW: lets just see if we can get this to work for the manifold and point atractor. Because
    # we are only looking for one point atractro (that is period one) we only need one variable to
    # store the position of the atractor. For the sake of pure rigor and learning we will be looking
    # for its full position in phase space. ie the position and velocity of the atractor. We could
    # probably get away with just looking at position but I think this is bad fourm.
    atractor = pl.array([0,0])
    # This is in the form [vel,pos] as the fuctions are in the PPclass.
    #---------------------------------------------------------------------------------------------



    # seperate poincare section files and full trajectory files
    for i,j in enumerate(files):
        #print(j[-16:])
        if(j[-16:]=="poincaredata.txt"):
            poincardat=pl.append(poincardat,j)
            #this = open(j,"r"
        elif(j!="info.dat"):
                dat=pl.append(dat,j)

    
    
    #for a in range(int(10000/20)-1):

    # Here we will just try to plot a block runing through ONE poincare section
    #make arrays to story the block points to be ploted
    firstBptsx=pl.array([])
    firstBptsvx=pl.array([])

    #initialize the figure first
    fig = pl.figure()
    ax = fig.add_subplot(111)


    for i,j in enumerate(poincardat):

        print(j)
        curfile = open(j,"r")
        curlines = curfile.readlines()


        #spilit the position velocity and time from the line
        firstline=curlines[1]

        #split the lines to parse the information
        firstx = float(firstline.split(",")[0])
        firstvx = float(firstline.split(",")[1])
        
        #firstBptsx = pl.append(firstBptsx,firstx)
        #firstBptsvx = pl.append(firstBptsvx,firstvx)

        # compair the last 5 poincare points to see if its the atractor or chaotic manifold
        # keep all the numbers in an array to compair in a loop? yes, thats what ill try to do
        compairx = pl.array([])
        compairvx= pl.array([])
        for i in range(5):
            nextline=curlines[i-5]


            #split the lines to parse the information
            nextx = float(nextline.split(",")[0])
            nextvx = float(nextline.split(",")[1])

            compairx = pl.append(secondBptsx,nextx)
            compairvx = pl.append(secondBptsvx,nextvx)

        # find maximum distance between last 5 pts
        # this function retuns an array of numbers that are the distances between the points
        # The input is just the compairx and compairvx arrays
        distances = whataredist(compairx,compairvx)
        maxdist = max(list(distances))

        ptcolor = 'r'
        if(maxdist<.2):
            ptcolor = 'b'
        ax.scatter(firstx,firstvx,color=ptcolor)
        #ax.scatter(secondBptsx,secondBptsvx)


        

            
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

