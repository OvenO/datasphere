import numpy
import sys
import os
sys.path.append(os.path.expanduser('~')+"/datasphere/ECproject")
import ECclass as ec
import pylab as pl
from scipy.integrate import odeint
from scipy.integrate import ode
import shutil
import time as thetime
import argparse
from datetime import datetime
from mpl_toolkits.mplot3d import Axes3D
import random


class One_Particle_Ensble_Sin1D(object):
    def __init__(self,A,beta,dt):
        self.A = A
        print('self.A is: ' + str(self.A))
        self.beta = beta
        print('self.beta is: '+str(self.beta))
        self.num_cell = 1.0
        # d here is the length of the system
        self.d = self.num_cell*2.0*pl.pi
        print('slef.d: ' +str(self.d))
        # right now As is only different becasue of different particle "densities". The reason I
        # have stated it like this is because particles with different chages would then need the qq
        # factor to actualy be q[i]*q[j]. or something like that. Lets just see if we can achive the
        # particle separation with the As method
        self.As = A
        self.dt = dt
        self.sol = pl.array([]) 
   
    # difine a function that grabs the matrix elements of the jacobiani, set_sol must have already
    # been done for this to work
    def set_sol(self,sol):
        self.sol=sol
        # to get the solution at a particular time we need the index that is assosiated with that
        # time. we get this by taking the time value wanted and deviding by dt. In order for this to
        # work with single (non array values) of time we need a self.dt to be defined.

    def J(self,which_M,t):
        x1 = self.sol[int(t/self.dt)-1,2]
        y =  self.sol[int(t/self.dt)-1,3]
        # define the matrix elements of the time dependent jacobian
        M11 = 0.0
        M12 = 1.0
        M21 = self.A*pl.cos(x1)*pl.cos(t)
        M22 = -self.beta


        if (which_M == "M11"):
            return M11
        if (which_M == "M12"):
            return M12
        if (which_M == "M21"):
            return M21
        if (which_M == "M22"):
            return M22

    def mw(self,warr,t):
        # to get the solution at a particular time we need the index that is assosiated with that
        # time. we get this by taking the time value wanted and deviding by dt. In order for this to
        # work with single (non array values) of time we need a self.dt to be defined.
        
        dotW11 = warr[2]
        dotW12 = warr[3]
        dotW21 = warr[0]*self.J("M21",t)+warr[2]*(-self.beta)
        dotW22 = warr[1]*self.J("M21",t)+warr[3]*(-self.beta)
        return [dotW11,dotW12,dotW21,dotW22]

    # just make normal functions to try to pass into odeint function. Should be much faster
    def f(self,xarr,t):
        x0dot = self.A*pl.sin(xarr[2])*pl.cos(t) - self.beta*xarr[0]
        x1dot = 0.0
        x2dot = xarr[0]
        x3dot = 0.0
        return [x0dot,x1dot,x2dot,x3dot]

# functions looks to see weather or not the curent point is in the threshold radius of the first
# point
# returns True if NOT in threshhold radius
# returns False if we found our guy
def not_close(first_pnt,curnt_pnt,thresh):
    rf = pl.array([first_pnt[0] , first_pnt[2]])
    rs = pl.array([curnt_pnt[0] , curnt_pnt[2]])
    diff = rf-rs
    r = pl.sqrt(diff[0]**2+diff[1]**2)
    #print("r is: "+str(r))

    if (r>thresh):
        return True
    else:
        return False


# find a single loop of the limit cycle. Might be periodoc over more than one cycle
# returns the solution of just that loop AND the periodicity of the loop
# takes a threshhold number. If it cant find a solution where the begining and end of the
# trajectroy lye within this threshold value than it quits and prints an error
#thresh is a distance in the phase plane
def find_one_full_closed(sol,thresh,dt):
    not_found = False
    # work our way backwards from last time value to find last period

    # first find last %2*pi position
    loc = len(sol[:,2])
    while ((loc*dt)%(2*pl.pi)>dt):
        loc-=1
    first_loc = loc 
    first_pnt = sol[first_loc,:]
    loc-=1
    # now find the next point where the orbit closes (going backward) 
    # orbits should have trajectories in multiples of 2*pi so only check those
    while ((loc*dt)%(2*pl.pi)>dt):
        loc-=1

    curnt_pnt = sol[loc,:]

    # for "slow" trajectories the point after the first may be within the threshold value. This is
    # not bad as it means the time step is definetely small enough but is messes up the next loop.
    # To fix this problem we will subtract more than one from loc. Not to much though otherwise we
    # risk crossing some 2*pi barier...though probably not. (left origonal loc-=1 for comparison).
    #loc -= 1
    # increas by pi/4
    loc -= int(pl.pi/4.0/dt)
    while (not_close(first_pnt,curnt_pnt,thresh)):
        if (loc == 0):
            print("Point in threshold not found!!")
            not_found = True
            #raise Exception("Point in threshold not found!!")
            break
        while ((loc*dt)%(2*pl.pi)>dt):
            loc-=1
        curnt_pnt = sol[loc,:]
        secnd_loc = loc
        loc-=1
    

    secnd_pnt = curnt_pnt

    if not_found:
        final = find_one_full_closed(sol,thresh*2,dt)
    else:
        final = sol[secnd_loc:first_loc+1,:]
    
    return final
def main():
    # do we want an image of the loops in 3D for different A?
    loops = True
    # do we want to make a movie of the stability multipliers in the complex plane?
    mk_stab_mov = True

    if loops:
        d3fig = pl.figure()
        d3ax = d3fig.add_subplot(111,projection='3d')

    os.mkdir('ScndLoopImgs')
    # make a directory for the stability multiplier images --> this will be a movie as a function of
    # A
    if mk_stab_mov: os.mkdir('ScndStabMovie')

    # this variable just exsits so we dont print the A value of the bifurcation point more than
    # once.
    found_bif = False

    # make file to store q (periodicity)
    q_file = open("scnd_qdata.txt","w")
    # make file to store stability multipliers
    eig_file  = open("scnd_data.txt","w")
    eig_file.write("eig1   eig2   A\n")

    dt = .0002 
    # total number of iterations to perform
    totIter = 5000000
    totTime = totIter*dt
    time = pl.arange(0.0,totTime,dt)
    
    beta = .6

    # how many cells is till periodicity use x = n*pi/k (n must be even #) modNum = 2*pl.pi/k
    modNum = 2.0*pl.pi
    
    # initial conditions for the periodobling cascade
    initx = pl.pi +.1
    inity = 0.0
    initvx = 0.0
    initvy = 0.0
    
    A_start = 0.956
    A = A_start
    A_max = .959
    # good for stability multipliers
    #A_step = .0005
    A_step = .000005
    
    count = 0

    # make arrays to keep eigen values. There willl be two eigen values so lets hve two seperate
    # arrays for them
    eigs1 = pl.array([])
    eigs2 = pl.array([])

    # file to write final poition of particle to
    final = open("scnd_final_position.txt","w")
    final.write("Last position of orbit,   A\n")
    x0 = pl.array([initvx,initvy,initx,inity])
    
    previous_q = 0.0
    while A < A_max:
        # initial conditions vector
        # set up: [xdot,ydot,x,y]
        apx = One_Particle_Ensble_Sin1D(A,beta,dt)
        sol = odeint(apx.f,x0,time)
        print("x0")
        print(x0)
        
        #sol[:,2]=sol[:,2]%(2*pl.pi)
        # find a single loop of the limit cycle. Might be periodoc over more than one cycle
        # returns the solution of just that loop AND the periodicity of the loop
        # takes a threshhold number. If it cant find a solution where the begining and end of the
        # trajectroy lye within this threshold value than it quits and prints an error
        #thresh is distance in the phase place
        #thresh = .01
        thresh = .0003
        
        # change this back for bifrucation other than FIRST PI BIF
        loop = find_one_full_closed(sol,thresh,dt)
        #loop = pl.zeros([int(2.0*pl.pi/dt),4])
        #loop[:,2]+=pl.pi


        if "stop" in loop:
            break

        loop_t = pl.arange(0.0,(len(loop))*dt,dt)
        
        if loops :
            d3ax.plot(pl.zeros(len(loop))+A,loop[:,2],loop[:,0],color="Black")


        #fig = pl.figure()
        #ax = fig.add_subplot(111)
        ##ax.scatter([0.0,pl.pi,2.0*pl.pi],[0.0,0.0,0.0],color="Red")
        ##ax.plot(loop[:,2],loop[:,0],":",color="Black")
        #ax.plot(loop[:,2],loop[:,0],color="Black")
        #ax.set_xlabel("$x_1$",fontsize=25)
        #ax.set_ylabel("$x_2$",fontsize=25)
        ##ax.set_xlim([pl.pi-pl.pi/3.0,pl.pi+pl.pi/3.0])
        ##ax.set_ylim([-.3,.3])
        #fig.tight_layout()
        #fig.savefig("ScndLoopImgs/"+str(A)+".png",dpi = 300,transparent=True)
        ##os.system("open LoopImgs/" +str(A)+".png")
        #pl.close(fig)
       
        apx.set_sol(loop)

        apx.set_sol(loop)
    
        w0 = pl.array([1.0,0.0,0.0,1.0])
        w_of_t = odeint(apx.mw,w0,loop_t,hmax=dt,hmin=dt)
        #w_of_t = odeint(apx.mw,w0,loop_t)
        
        current_q = loop_t[-1]/(2.0*pl.pi)
        # print the period of the orbit we are working on
        print("q: " + str(current_q))
        q_file.write(str(loop_t[-1]/(2.0*pl.pi))+" "+str(A)+"\n")

        if current_q > (previous_q+1.0):
            print("bifurcation point. A = " +str(A))

        previous_q=current_q
    
        # make the matrix form of w_of_t
        matrix = w_of_t[-1,:].reshape(2,2)
        
        # use linalg to get the eigen values of the W(t=q) where q is the period time of the orbit
        vals,vect = numpy.linalg.eig(matrix) 
        
        if((abs(vals[0])<=1.0) and (not found_bif)):
            print("this is the bifurcation point (l1)")
            print(A)
            found_bif = True
        if(abs(vals[1])<=1.0 and (not found_bif)):
            print("this is the bifurcation point (l2)")
            print(A)
            found_bif = True

        eigs1 = pl.append(eigs1,vals[0])
        eigs2 = pl.append(eigs2,vals[1])

        eig_file.write(str(vals[0])+" "+str(vals[1])+" "+str(A)+"\n")

        count+=1
        x0 = loop[-1,:]
        final.write(str(x0)[1:-1]+" "+str(A) +"\n")
        A += A_step
        print("A: "+str(A))

    theta = pl.arange(0,10,.05)
    A_arr = pl.arange(A_start,A,A_step)

    print('we are above')
    while len(A_arr)>len([k.real for k in eigs1]):
        A_arr = A_arr[:-1]
    while len(A_arr)<len([k.real for k in eigs1]):
        A_arr = pl.append(A_arr,A_arr[-1]+A_step)
    print('we are below')

    fig1 = pl.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(pl.cos(theta),pl.sin(theta))
    ax1.plot([k.real for k in eigs1],[l.imag for l in eigs1])
    ax1.set_xlabel("Re[$\lambda_1$]",fontsize=25)
    ax1.set_ylabel("Im[$\lambda_1$]",fontsize=25)
    fig1.tight_layout()
    fig1.savefig("scnd_eig1.png")
    os.system("open scnd_eig1.png")

    fig2 = pl.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot(pl.cos(theta),pl.sin(theta))
    ax2.plot([k.real for k in eigs2],[l.imag for l in eigs2])
    ax2.set_xlabel("Re[$\lambda_2$]",fontsize=25)
    ax2.set_ylabel("Im[$\lambda_2$]",fontsize=25)
    fig2.tight_layout()
    fig2.savefig("scnd_eig2.png")
    os.system("open scnd_eig2.png")

    fig3, ax3 = pl.subplots(2,sharex=True)
    ax3[0].plot(A_arr,[k.real for k in eigs1],color='k')
    ax3[1].plot(A_arr,[k.imag for k in eigs1],color='k')
    ax3[0].set_ylabel("Re[$\lambda_1$]",fontsize = 25)
    ax3[1].set_ylabel("Im[$\lambda_1$]",fontsize = 25)
    ax3[1].set_xlabel("$A$",fontsize = 25)
    fig3.tight_layout()
    fig3.savefig("scnd_A_vs_eig1.png")
    os.system("open scnd_A_vs_eig1.png")

    fig4, ax4 = pl.subplots(2,sharex=True)
    ax4[0].plot(A_arr,[k.real for k in eigs2], color = 'k')
    ax4[1].plot(A_arr,[k.imag for k in eigs2], color = 'k')
    ax4[0].set_ylabel("Re[$\lambda_2$]",fontsize = 25)
    ax4[1].set_ylabel("Im[$\lambda_2$]",fontsize = 25)
    ax4[1].set_xlabel("$A$",fontsize = 25)
    fig4.tight_layout()
    fig4.savefig("scnd_A_vs_eig2.png")
    os.system("open scnd_A_vs_eig2.png")

    eig_file.close()

    final.close()    
        
    ## make text file with all extra information
    #outFile = open("info.dat","w")
    #outFile.write("Info \n coefficient: " + str(coef) \
    #        + "\nwave number: " +str(k)\
    #        + "\nomega: " + str(w)\
    #        + "\ndamping: " + str(damp)\
    #        + "\ng: " + str(g)\
    #        + "\ntime step: " + str(dt)\
    #        + "\ntotal time: " + str(dt*totIter)\
    #        + "\ntotal iterations: " + str(totIter)\
    #        + "\nInitial Conditions: \n" +
    #        "initial x: " +str(initx) \
    #        +"\ninitial y: " +str(inity) \
    #        +"\ninitial vx: " +str(initvx)\
    #        +"\ninitial vy: " +str(initvy) )
    #outFile.close()
    
    # line for stable static point

    if loops:
        line = pl.arange(A_start-.01,A_start,A_step)
        pi_line = pl.zeros(len(line))+pl.pi
        z_line = pl.zeros(len(line))
        d3ax.plot(line,pi_line,z_line,color="Black")
        d3ax.set_xlabel("$A$",fontsize=25)
        d3ax.set_ylabel("$x_1$",fontsize=25)
        d3ax.set_zlabel("$x_2$",fontsize=25)
        d3fig.tight_layout()
        d3fig.savefig("scnd_loops.png",dpi=300)

    if mk_stab_mov:
        for i in range(len(eigs1)):
            s_fig = pl.figure()
            s_ax = s_fig.add_subplot(111)
            s_ax.plot(pl.cos(theta),pl.sin(theta))
            s_ax.scatter(eigs1[i].real,eigs1[i].imag,c='r',s=20)
            s_ax.scatter(eigs2[i].real,eigs2[i].imag,c='b',s=20)
            s_ax.set_xlabel("Re[$\lambda$]",fontsize=25)
            s_ax.set_ylabel("Im[$\lambda$]",fontsize=25)
            s_fig.tight_layout()
            s_fig.savefig("%(num)0.5d_stbl.png"%{"num":i})
            pl.close(s_fig)



if __name__ == '__main__':
    main()
