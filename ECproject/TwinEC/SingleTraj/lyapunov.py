import numpy
import sys
sys.path.append("/users/o/m/omyers/puthere/lib/python2.7/site-packages")
sys.path.append("/users/o/m/omyers/datasphere/ECproject")
sys.path.append("/users/o/m/omyers/puthere/matplotlib/build/lib.linux-x86_64-2.7")
import ECclass as ec
import pylab as pl
import os
from scipy.integrate import odeint
from scipy.integrate import ode
import shutil
import time as thetime
import argparse
from datetime import datetime
#from mpl_toolkits.mplot3d import Axes3D
import random

#**********************************************************************************************
#**********************************************************************************************
# Here we will find a second point that is some small distace (epsilon) from the x0 point.
# We need this function because we want the distace to be epsilon but the position about x0 to be
# random
def get_first_init(x0,epsilon):
    # we will use a change of coordinates to get the location of the particle relative to x0. First
    # we just find some random point a distace epsilon from the origin.

    # This  is a little more tricky in 4 dimensions
    theta = random.random()*2.0*pl.pi
    phi = random.random()*2.0*pl.pi
    psi = random.random()*2.0*pl.pi

    x =  epsilon*pl.cos(psi)
    vx = epsilon*pl.sin(psi)*pl.cos(phi)
    y =  epsilon*pl.sin(psi)*pl.sin(phi)*pl.sin(theta)
    vy = epsilon*pl.sin(psi)*pl.sin(phi)*pl.cos(theta)
    
    # change of coordinates
    return x0 + pl.array([vx,vy,x,y])

#**********************************************************************************************
#**********************************************************************************************
# function just finds the magnitude of the distance of the two input points
def distance(x1,x2):
    x_sqrd  = (x1[0]-x2[0])**2
    vx_sqrd = (x1[1]-x2[1])**2
    y_sqrd  = (x1[2]-x2[2])**2
    vy_sqrd = (x1[3]-x2[3])**2
    return pl.sqrt(x_sqrd + y_sqrd + vx_sqrd + vy_sqrd)

#**********************************************************************************************
#**********************************************************************************************
# this function needs both the final positions of the two trajectories and the variable eplilon wich
# is the starting/reset DISTANCE of the two poinnts. jsee lab book #2 pg 59 for a
# more detailed explination. "renomalize" will retun a point that is on the line contecting the
# final points of the trajectorys and will this new point will be a distance epsilon from the
# unpurtubed trjectory (x_unpurt). x_putub is final position of purtubed trajectorie
def renormalize(x_unpurt,x_puturb,epsilon):
    final_dist = distance(x_unpurt,x_puturb)
    
    # just initialize xnew
    xnew = pl.array([0.0,0.0,0.0,0.0])

    # the new renormalized vx (see lab book #2 pg 61)
    xnew[0] = x_unpurt[0]+(epsilon/final_dist)*(x_puturb[0]-x_unpurt[0])
    xnew[1] = x_unpurt[1]+(epsilon/final_dist)*(x_puturb[1]-x_unpurt[1])
    xnew[2] = x_unpurt[2]+(epsilon/final_dist)*(x_puturb[2]-x_unpurt[2])
    xnew[3] = x_unpurt[3]+(epsilon/final_dist)*(x_puturb[3]-x_unpurt[3])

    return xnew
#**********************************************************************************************
#**********************************************************************************************
def get_poin(sol,dt):
    poin = pl.array([])
    for i in range(len(sol)):
        if(((i*dt)%(2.0*pl.pi))<=dt):
            poin = pl.append(poin,sol[i,:])
    poin = poin.reshape(-1,4)
    # flip order of array
    #poin = poin[::-1,:]
    return poin
#**********************************************************************************************
#**********************************************************************************************
def plot_first_sol(sol,dt):
    
    poin = get_poin(sol,dt)
    poin[:,2]=poin[:,2]%(2*pl.pi)

    first_fig = pl.figure()
    first_ax = first_fig.add_subplot(111)
    first_ax.scatter(poin[:,2],poin[:,0],color="Black",s=.1)
    first_ax.set_xlabel("$x_1$",fontsize=25)
    first_ax.set_ylabel("$x_2$",fontsize=25)
    #first_ax.set_xlim([0,2*pl.pi])
    #first_ax.set_ylim([-1.3,1.3])
    first_fig.tight_layout()
    first_fig.savefig("first_plot.png")
    os.system("open first_plot.png")
    
#**********************************************************************************************
#**********************************************************************************************
def main():
    
    # initial purturbation size
    epsilon = 1.0e-7
    print('epsilon is: '+str(epsilon))



    # period variable should probably just be kept at the actual period (2*pl.pi) but the p
    # this variable is to alow us to changhe the lenght of time we wate before we colect th
    # distance information and renormalize. This is also nesssasary in the final calculatio
    # LE becae LE = 1/period * ln(rm/r0).
    period = 2.0*pl.pi
    print('period is: '+str(period))

    # first throw away (to get initial conditions in atractor)
    throw_away_1 = 250
    print('throw_away_1 is: '+str(throw_away_1))
    #throw_away_1 = 200
    # second throw away (to make sure puturbed trajectory is aligned)
    throw_away_2 = 100
    print('throw_away_2 is: '+str(throw_away_2))
    #throw_away_2 = 200

    # number of cycle-sets to go through
    # try 
    num = 55000
    print('number of cycles is: '+str(num))
    # works
    #num = 60000
    #num = 16000
    #num = 30000
    #num = 8000
    
    dt = .001 
    print('dt is: '+str(dt))
    # total number of sterations to perform inorder to get cycles rioht
    totIter = period/dt
    totTime = period
    time = pl.arange(0.0,totTime,dt)
    
    #time array for the first run to get initial condition in strange atractor
    first_time = pl.arange(0.0,throw_away_1*2.0*pl.pi,dt)
    
    surf = 1.0
    coef = 9.0
    k = 1.0
    w = 1.0
    damp = .05
    g = .1

    # how many cells is till periodicity use x = n*pi/k (n must be even #) modNum = 2*pl.pi/k
    modNum = 2.0*pl.pi
    
    # some random initial conditions
    initx = 3.649
    inity = 1.3
    initvx = .1237
    initvy = 0.0

    # now find initial condition in stange atractor by running far in time. print this to use it for
    # later
    x0 = pl.array([initvx,initvy,initx,inity])
    
    apx = ec.CentreLineApx(coef,k,w,damp,surf,g,dt)
    first_sol = odeint(apx.f,x0,first_time)
    #plot_first_sol(first_sol,dt)
    # now reset x0
    first_sol[:,2]=first_sol[:,2]%(2*pl.pi)
    x0 = first_sol[-1,:]
    
    x_other = get_first_init(x0,epsilon)
    
    print("x0 from first_sol is:")
    print(x0)
    print("x_other (purturbed is:")
    print x_other

    # two arrays to store solution data
    arr_orig = pl.array([])
    arr_othr = pl.array([])

    full_orig = pl.array([])
    full_othr = pl.array([])
    
    # array to keep distance after driving cycle 
    darr = pl.array([]) 
    watch = pl.array([])
    watch_le = 0.0

    for i in range(num + throw_away_2):
        sol = odeint(apx.f,x0,time)
        sol_other = odeint(apx.f,x_other,time)

        sol[:,2]=sol[:,2]%(2*pl.pi)
        sol_other[:,2]=sol_other[:,2]%(2*pl.pi)

        #arr_orig = pl.append(arr_orig,sol[-1,:])
        #arr_othr = pl.append(arr_othr,sol_other[-1,:])
    
        if i> throw_away_2: 
            #get the new distance
            darr = pl.append(darr,distance(sol[-1,:],sol_other[-1,:]))
            
            watch_le += pl.log(abs(darr[-1]/epsilon))
            cur_avg = watch_le/(i-throw_away_2+1)/period
            watch = pl.append(watch,cur_avg)


        #full_orig = pl.append(full_orig,sol)
        #full_othr = pl.append(full_othr,sol_other)

#        poin = get_poin(sol,dt)
#        poin_other = get_poin(sol_other,dt)

        x0 = sol[-1,:]
        x_other = renormalize(x0,sol_other[-1,:],epsilon)
    
    #arr_orig = arr_orig.reshape(-1,4)
    #arr_othr = arr_othr.reshape(-1,4)
    
    #full_orig = full_orig.reshape(-1,4)
    #full_othr = full_othr.reshape(-1,4)

    eps_arr = pl.zeros(len(darr))+epsilon
    le = pl.zeros(len(darr))+pl.log(abs(darr/epsilon))/period

    le_avg = 0.0
    for i,j in enumerate(le):
        le_avg += j
    le_avg = le_avg/len(le)

    print("le is")
    print(le_avg)

    fig = pl.figure()
    ax = fig.add_subplot(111)
    ax.scatter(pl.arange(len(watch)),watch,s=.1)
    #ax.set_xlabel("$x_1$",fontsize=25)
    #ax.set_ylabel("$x_2$",fontsize=25)
    #ax.set_xlim([0,2*pl.pi])
    #ax.set_ylim([-1.3,1.3])
    #fig.tight_layout()
    fig.savefig("convergence.png")
    os.system("open convergence.png")
 

    ## plot of PC sections (red Blue)
    #fig = pl.figure()
    #ax = fig.add_subplot(111)
    ##ax.scatter([0.0,pl.pi,2.0*pl.pi],[0.0,0.0,0.0],color="Red")
    ##ax.scatter(arr_orig[:,2],arr_orig[:,0],color="Red",s=.1)
    ##ax.scatter(arr_othr[:,2],arr_othr[:,0],color="Blue",s=.1)
    #ax.scatter(full_orig[:,2],full_orig[:,0],color="Red",s=.1)
    #ax.scatter(full_othr[:,2],full_othr[:,0],color="Blue",s=.1)
    #ax.set_xlabel("$x_1$",fontsize=25)
    #ax.set_ylabel("$x_2$",fontsize=25)
    ##ax.set_xlim([0,2*pl.pi])
    ##ax.set_ylim([-1.3,1.3])
    #fig.tight_layout()
    #fig.savefig("sep.png")
    #os.system("open sep.png")
    
if __name__ == '__main__':
    main()
