import numpy
import sys
import os
sys.path.append(os.path.expanduser('~')+"/datasphere/ECproject")
import ECclass as ec
import pylab as pl
from scipy.integrate import odeint
from scipy.integrate import ode
from scipy.special import polygamma
import shutil
import time as thetime
import argparse
from datetime import datetime
#from mpl_toolkits.mplot3d import Axes3D
import random


class Two_Particle_Sin1D(object):
    def __init__(self,A,beta,qq,dt):
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
        self.qq = qq
        self.sol = pl.array([]) 
   
    # difine a function that grabs the matrix elements of the jacobiani, set_sol must have already
    # been done for this to work
    def set_sol(self,sol):
        self.sol=sol
        # to get the solution at a particular time we need the index that is assosiated with that
        # time. we get this by taking the time value wanted and deviding by dt. In order for this to
        # work with single (non array values) of time we need a self.dt to be defined.


    def J(self,t):
        # the -1 in the lines below is for the right rounding with int()
        # x1 is position of p1 (particle 1) 
        x1 = self.sol[int(t/self.dt)-1,2]
        # x2 is velocity of p1
        x2 = self.sol[int(t/self.dt)-1,0]
        # x3 is position of p2
        x3 = self.sol[int(t/self.dt)-1,3]
        # x4 is velocity of p2
        x4 = self.sol[int(t/self.dt)-1,1]

        # These are the differentials of the forces of the particles. Writen like this to make the
        # matrix below easier to read f14 is force of p2 on p1 _dx1 is derivitive with respect to x1
        # Note on 1/r2 part -> goes to cubic so it will always retain its sign.
        df13_dx1 = -2.0/(x1-x3)**3 + (polygamma(2,1.0+(x1-x3)/self.d)+polygamma(2,1.0-(x1-x3)/self.d))/self.d**3
        # the final deriviative of -x3 just gives you the negative of everything above
        df13_dx3 = -df13_dx1
        df31_dx1 = 2.0/(x3-x1)**3 - (polygamma(2,1.0-(x3-x1)/self.d)+polygamma(2,1.0+(x3-x1)/self.d))/self.d**3
        df31_dx3 = -df31_dx1


        # define the matrix elements of the time dependent jacobian
        jacobian = pl.array([ \
        [0.0                                   , 1.0       , 0.0                                   , 0.0],
        [self.A*pl.cos(x1)*pl.cos(t)+df13_dx1, -self.beta, df13_dx3                              , 0.0],
        [0.0                                   , 0.0       , 0.0                                   , 1.0],
        [df31_dx1                              , 0.0       , self.A*pl.cos(x3)*pl.cos(t)+df31_dx3, -self.beta]\
        ])

        return jacobian
        # Keep this incase we can't get the actual array version to work and we switch back to wanting
        # to do it with the M## variables
        #jacobian = pl.array([ \
        #[M11=0.0                                       M12=1.0      M13=0.0                             M14=0.0                                     M15=0.0]
        #[M21=-self.A*pl.cos(x1)*pl.cos(x3)+df14_dx1    M22=0.0      M23=self.A*pl.sin(x1)*pl.sin(x3)    M24=df14_dx4                                M25=0.0]
        #[M31=0.0                                       M32=0.0      M33=0.0                             M34=0.0                                     M35=0.0]
        #[M41=0.0                                       M42=0.0      M43=0.0                             M44=0.0                                     M45=1.0]
        #[M51=df41_dx1                                  M52=0.0      M53=self.A*pl.sin(x4)*pl.sin(x3)    M54=-self.A*pl.cos(x4)*pl.cos(x3)+df41_dx4  M55=0.0]\
        #]

        #if (which_M == "M11"):
        #    return M11
        #if (which_M == "M12"):
        #    return M12
        #if (which_M == "M21"):
        #    return M21
        #if (which_M == "M22"):
        #    return M22

    def mw(self,warr,t):
        jacobian = self.J(t)

        # odeint needs a 1d array but it will be easier for us to work with a matrix. This reshpae
        # just puts it into its indended form
        W = warr.reshape(4,4)
        
        # the matrix product of the jacobian and the solution matrix (W) is what we need to solve
        dot_W_matrix = pl.dot(jacobian,W)
        
        # get it back into 1d form for odeint
        return dot_W_matrix.reshape(-1)

        # Origonal for N=1. Keeping for reference for now
        #dotW11 = warr[2]
        #dotW12 = warr[3]
        #dotW21 = warr[0]*self.J("M21",t)+warr[2]*(-self.beta)
        #dotW22 = warr[1]*self.J("M21",t)+warr[3]*(-self.beta)
        #return [dotW11,dotW12,dotW21,dotW22]

    def f(self,x,t):
        # This one line is different than in ECclass.py
        N = 2
        xdot = pl.array([])

        # modulus the x for periodicity.
        x[N:2*N]= x[N:2*N]%self.d
        # HERE ---->> 1Dify
        for i in range(N):
            temp = 0.0
            for j in range(N):
                if i == j:
                    continue
                #repulsive x interparticle force of j on i
                temp += self.qq*(x[N+i]-x[N+j])/(pl.sqrt((x[N+i]-x[N+j])**2)**3)
                # All of the forces coming from the 'same' paricle but from other 'cells' due to the
                # periodic contrains can be wraped up in a sum that converges to an aswer that can
                # be expressed in terms of polygamma functions (se pg 92 of notebook).
                # Note on the sign (xi-xj or xj-xi). Changing the sign of the xi-xj term (i.e. which
                # particle are we considering forces on) changes the direction of the force
                # apropriately.
                temp += self.qq*(polygamma(1,(self.d+x[N+i]-x[N+j])/self.d)-polygamma(1,1.0-((x[N+i]-x[N+j])/self.d)))/(self.d**2)
            # periodic force on particle i
            temp += self.A*pl.sin(x[N+i])*pl.cos(t)
            temp -= self.beta*x[i]
            xdot = pl.append(xdot,temp)
        for i in range(N):
            xdot = pl.append(xdot,x[i])
        return xdot

#*********************************************************************************************************************
#*********************************************************************************************************************
# THESE TWO FUNCTIONS ARE CURRENTLY NOT SET UP TO RUN FOR THE 2 PARTICLE CASE
# They are not nesassary for the first bifurcatoin anyway so lets just get that working first.
#*********************************************************************************************************************
#*********************************************************************************************************************
## functions looks to see weather or not the curent point is in the threshold radius of the first
## point
## returns True if NOT in threshhold radius
## returns False if we found our guy
#def not_close(first_pnt,curnt_pnt,thresh):
#    rf = pl.array([first_pnt[0] , first_pnt[2]])
#    rs = pl.array([curnt_pnt[0] , curnt_pnt[2]])
#    diff = rf-rs
#    r = pl.sqrt(diff[0]**2+diff[1]**2)
#    #print("r is: "+str(r))
#
#    if (r>thresh):
#        return True
#    else:
#        return False


## find a single loop of the limit cycle. Might be periodoc over more than one cycle
## returns the solution of just that loop AND the periodicity of the loop
## takes a threshhold number. If it cant find a solution where the begining and end of the
## trajectroy lye within this threshold value than it quits and prints an error
##thresh is a distance in the phase plane
#def find_one_full_closed(sol,thresh,dt):
#    not_found = False
#    # work our way backwards from last time value to find last period
#
#    # first find last %2*pi position
#    loc = len(sol[:,2])
#    while ((loc*dt)%(2*pl.pi)>dt):
#        loc-=1
#    first_loc = loc 
#    first_pnt = sol[first_loc,:]
#    loc-=1
#    # now find the next point where the orbit closes (going backward) 
#    # orbits should have trajectories in multiples of 2*pi so only check those
#    while ((loc*dt)%(2*pl.pi)>dt):
#        loc-=1
#
#    curnt_pnt = sol[loc,:]
#
#    # for "slow" trajectories the point after the first may be within the threshold value. This is
#    # not bad as it means the time step is definetely small enough but is messes up the next loop.
#    # To fix this problem we will subtract more than one from loc. Not to much though otherwise we
#    # risk crossing some 2*pi barier...though probably not. (left origonal loc-=1 for comparison).
#    #loc -= 1
#    # increas by pi/4
#    loc -= int(pl.pi/4.0/dt)
#    while (not_close(first_pnt,curnt_pnt,thresh)):
#        if (loc == 0):
#            print("Point in threshold not found!!")
#            not_found = True
#            #raise Exception("Point in threshold not found!!")
#            break
#        while ((loc*dt)%(2*pl.pi)>dt):
#            loc-=1
#        curnt_pnt = sol[loc,:]
#        secnd_loc = loc
#        loc-=1
#    
#
#    secnd_pnt = curnt_pnt
#
#    if not_found:
#        final = find_one_full_closed(sol,thresh*2,dt)
#    else:
#        final = sol[secnd_loc:first_loc+1,:]
#    
#    return final

def main():
    # do we want an image of the loops in 3D for different A?
    loops = False
    # do we want to make a movie of the stability multipliers in the complex plane?
    mk_stab_mov = True

    if loops:
        fig3d = pl.figure
        d3ax = fig.add_subplot(111,projection='3d')

    #os.mkdir('LoopImgs')

    # make a directory for the stability multiplier images --> this will be a movie as a function of
    # A
    if mk_stab_mov: os.mkdir('StabMovie')

    # this variable just exsits so we dont print the A value of the bifurcation point more than
    # once.
    found_bif = False

    # make file to store q (periodicity)
    q_file = open("qdata.txt","w")

    # make file to store stability multipliers
    eig_file  = open("data.txt","w")
    eig_file.write("eig1   eig2   A\n")

    dt = .001 
    # total number of iterations to perform
    # totIter = 10000000
    totIter = 50000
    totTime = totIter*dt
    time = pl.arange(0.0,totTime,dt)
    
    beta = .6
    qq = 1.0

    # how many cells is till periodicity use x = n*pi/k (n must be even #) modNum = 2*pl.pi/k
    modNum = 2.0*pl.pi
    
    # initial conditions for p1 and p2
    p1_init_x = pl.pi
    p1_init_vx = 0.0

    p2_init_x = 0.0
    p2_init_vx = 0.0

    A_start = 0.2
    A = A_start
    A_max = .8
    A_step = .01

    
    count = 0

    # make arrays to keep eigen values. There willl be four eigen values for N=2 so lets hve two seperate
    # arrays for them
    eigs1 = pl.array([])
    eigs2 = pl.array([])
    eigs3 = pl.array([])
    eigs4 = pl.array([])


    # file to write final poition of particle to
    final = open("final_position.txt","w")
    final.write("Last position of orbit,   A\n")
    x0 = pl.array([p1_init_vx,p2_init_vx,p1_init_x,p2_init_x])
    
    previous_q = 0.0
    while A < A_max:
        # initial conditions vector
        # set up: [xdot,ydot,x,y]
        apx = Two_Particle_Sin1D(A,beta,qq,dt)
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
        thresh = .00005
        
        # change this back for bifrucation other than FIRST PI BIF
        #loop = find_one_full_closed(sol,thresh,dt)
        loop = pl.zeros([int(2.0*pl.pi/dt),4])
        loop[:,2]+=pl.pi
        # the other particle needs to be at zero. Already is

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
        #fig.savefig("LoopImgs/"+str(A)+".png",dpi = 300,transparent=True)
        ##os.system("open LoopImgs/" +str(A)+".png")
        #pl.close(fig)
       
        apx.set_sol(loop)

    
        # solution matrix at t=0 is identity matrix -> we need it in 1d arrar form for odeint though
        # -> reshape takes care of the the form
        w0 = pl.identity(4).reshape(-1)
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
        matrix = w_of_t[-1,:].reshape(4,4)

        print('solution matrix is:')
        for i in range(len(matrix)):
            print(str(matrix[i,:]))

        # print('solution matrix times the initial conditions vector:')
        # This is wrong because the x0 needs to be re-ordered to b e consistant with the way the
        # Jacobian is set up.
        # print(pl.dot(matrix,x0))
        
        
        # use linalg to get the eigen values of the W(t=q) where q is the period time of the orbit
        vals,vect = numpy.linalg.eig(matrix) 
        print('the trace is' + str(matrix.trace()))
        print('determinant is: '+ str(pl.det(matrix)))
        
        if((abs(vals[0])<=1.0) and (not found_bif)):
            print("this is the bifurcation point (l1)")
            print(A)
            found_bif = True
        if(abs(vals[1])<=1.0 and (not found_bif)):
            print("this is the bifurcation point (l2)")
            print(A)
            found_bif = True
        
        print('number of eigen values is: ' + str(len(vals)))
        eigs1 = pl.append(eigs1,vals[0])
        eigs2 = pl.append(eigs2,vals[1])
        eigs3 = pl.append(eigs3,vals[2])
        eigs4 = pl.append(eigs4,vals[3])


        eig_file.write(str(vals[0])+" "+str(vals[1])+" "+str(vals[2])+" "+str(vals[3])+" "+str(A)+"\n")

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
    fig1.savefig("eig1.png")
    os.system("open eig1.png")

    fig2 = pl.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot(pl.cos(theta),pl.sin(theta))
    ax2.plot([k.real for k in eigs2],[l.imag for l in eigs2])
    ax2.set_xlabel("Re[$\lambda_2$]",fontsize=25)
    ax2.set_ylabel("Im[$\lambda_2$]",fontsize=25)
    fig2.tight_layout()
    fig2.savefig("eig2.png")
    os.system("open eig2.png")

    fig3, ax3 = pl.subplots(2,sharex=True)
    ax3[0].plot(A_arr,[k.real for k in eigs1],color='k')
    ax3[1].plot(A_arr,[k.imag for k in eigs1],color='k')
    ax3[0].set_ylabel("Re[$\lambda_1$]",fontsize = 25)
    ax3[1].set_ylabel("Im[$\lambda_1$]",fontsize = 25)
    ax3[1].set_xlabel("$A$",fontsize = 25)
    fig3.tight_layout()
    fig3.savefig("A_vs_eig1.png")
    os.system("open A_vs_eig1.png")

    fig4, ax4 = pl.subplots(2,sharex=True)
    ax4[0].plot(A_arr,[k.real for k in eigs2], color = 'k')
    ax4[1].plot(A_arr,[k.imag for k in eigs2], color = 'k')
    ax4[0].set_ylabel("Re[$\lambda_2$]",fontsize = 25)
    ax4[1].set_ylabel("Im[$\lambda_2$]",fontsize = 25)
    ax4[1].set_xlabel("$A$",fontsize = 25)
    fig4.tight_layout()
    fig4.savefig("A_vs_eig2.png")
    os.system("open A_vs_eig2.png")

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
        line = pl.arange(start_A-.01,start_A,A_step)
        pi_line = pl.zeros(len(line))+pl.pi
        z_line = pl.zeros(len(line))
        d3ax.plot(line,pi_line,z_line,color="Black")
        d3ax.set_xlabel("$A$",fontsize=25)
        d3ax.set_ylabel("$x_1$",fontsize=25)
        d3ax.set_zlabel("$x_2$",fontsize=25)
        d3fig.tight_layout()
        d3fig.savefig("loops.png",dpi=300)

    if mk_stab_mov:
        for i in range(len(eigs1)):
            s_fig = pl.figure()
            s_ax = s_fig.add_subplot(111)
            s_ax.plot(pl.cos(theta),pl.sin(theta))
            s_ax.scatter(eigs1[i].real,eigs1[i].imag,c='r',s=20)
            s_ax.scatter(eigs2[i].real,eigs2[i].imag,c='b',s=20)
            s_ax.scatter(eigs3[i].real,eigs3[i].imag,c='b',s=20)
            s_ax.scatter(eigs4[i].real,eigs4[i].imag,c='b',s=20)
            s_ax.set_xlabel("Re[$\lambda$]",fontsize=25)
            s_ax.set_ylabel("Im[$\lambda$]",fontsize=25)
            s_fig.tight_layout()
            s_fig.savefig("StabMovie/%(num)0.5d_stbl.png"%{"num":i})
            pl.close(s_fig)



if __name__ == '__main__':
    main()
