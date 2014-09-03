from datetime import datetime
import sys
sys.path.append("/users/o/m/omyers/puthere/lib/python2.7/site-packages")
import scipy as pl
from scipy.special import polygamma
import os

# List of classes and 1-2 word descriptions:
# class OurHMF(object): our special Hamiltonian Mean field theory.
# class Test(object): see if we can reproduce Molucular Dynamics results from computational book
# class SinSin2D(object):
# class Sin1D(object):
# class SecSin1D(object): ancored charge boundary conditions
# class One_Particle_Ensble_Sin1D(object):
# class HardCoreSin1D(object):
# class CentreLineApx(object):
# class TwinECApx(object):
# class MultiTwinECApxQuad(object):
# class MultiTwinECApxQuadPeriodic(object):
# class ExactPeriodic1D(object):
# class Flat2D(object):
# class travPlaneChrgApx(object):
# class surfTravPlaneChrgApx(object):
# class surfCentreLineApx(object):

# using the scipy odeint functions as a way of solving first order differential equations can be a
# little confusing so here are some notes:
# first of all the new array structure being passed back and forth for a full location in phase
# space is x0dot = x doubledot (or vx dot)
#          x1dot = y doubledot (or vy dot)
#          x2dot = x
#          x3dot = y                                
# this is the form of the array "xarr" that is being passed into the the f(self,xarr,t) function.
# the "t" in this case is the time array that we want our solution for. Find the youtube guy if you
# forget what goes where. (nonlinear pendulub ex)

#**************************************************************************************************
#**************************************************************************************************
class OurHMF(object):
   
    def __init__(self,N,I):
        # interaction stregth (dimless)
        self.I = I
        # number of particles
        self.N = N

    def f(self,x,t):
        # x[N+i] is the position of i^th particle
        # x[i] is the velocity of i^th particle

        # This initilazation should be slightly faster than re-doing "len(x)"
        xdot = pl.array([])

        for i in range(self.N):
            # temp is to have a variable to keep adding terms of the "field" and the interactions to
            # so we can do things peicewise.  start with the field contribution
            temp = 0.0
            #print('1st temp:' + str(temp))
            for j in range(self.N):
                if i == j:
                    continue
                # find phi_i and phi_j needed in the calculation of the interaction
                #phi_i = (pl.pi/self.N)*(i*2+x[self.N+i])
                #phi_j = (pl.pi/self.N)*(j*2+x[self.N+j])
                phi_i = (pl.pi/self.N)*i*2+x[self.N+i]
                phi_j = (pl.pi/self.N)*j*2+x[self.N+j]
                # interaction force of j on i
                #temp += -(pl.pi * self.I/(2.0*self.N))*pl.sin(phi_i-phi_j)*pl.cos(x[self.N+i])
                # THis is right but try rescaling to compair with HMF model
                #temp += -(pl.pi * self.I/(2.0*self.N))*pl.sin(phi_i-phi_j)
                # rescaled one
                temp += -(self.I/(2.0*self.N))*pl.sin(phi_i-phi_j)
            
            #print('2nd temp:' + str(temp))
            # What we have as temp is the second derivative of the position i.e v_dot -> so these
            # are the first in the array
            xdot = pl.append(xdot,temp)
            #print('xdot: ' + str(xdot))
        # The second in the array are just the first deriviative of the positions i.e v wich is what
        # were the first values in the input array.
        xdot = pl.append(xdot,x[:self.N])
        
        return xdot

#class OurHMF(object):
#   
#    def __init__(self,N,O_mega,A,B):
#        # O_mega is the natural frequency omega_0 over the driving frequency gamma?
#        self.O_mega = O_mega
#        # dimensionless driving amplituede 
#        self.A = A
#        # interaction stregth (dimless)
#        self.B = B
#        # number of particles
#        self.N = N
#
#    def f(self,x,t):
#        # x[N+i] is the position of i^th particle
#        # x[i] is the velocity of i^th particle
#
#        # This initilazation should be slightly faster than re-doing "len(x)"
#        xdot = pl.array([])
#
#        for i in range(self.N):
#            # temp is to have a variable to keep adding terms of the "field" and the interactions to
#            # so we can do things peicewise.  start with the field contribution
#            temp = -(self.O_mega+self.A*pl.cos(t))*pl.sin(x[self.N+i])
#            #print('1st temp:' + str(temp))
#            for j in range(self.N):
#                if i == j:
#                    continue
#                # find phi_i and phi_j needed in the calculation of the interaction
#                phi_i = (pl.pi/self.N)*(i*2+pl.sin(x[self.N+i]))
#                phi_j = (pl.pi/self.N)*(j*2+pl.sin(x[self.N+j]))
#                # interaction force of j on i
#                temp += -self.B*pl.sin(phi_i-phi_j)*pl.cos(x[self.N+i])
#            
#            #print('2nd temp:' + str(temp))
#            # What we have as temp is the second derivative of the position i.e v_dot -> so these
#            # are the first in the array
#            xdot = pl.append(xdot,temp)
#            #print('xdot: ' + str(xdot))
#        # The second in the array are just the first deriviative of the positions i.e v wich is what
#        # were the first values in the input array.
#        xdot = pl.append(xdot,x[:self.N])
#        
#        return xdot

#**************************************************************************************************
#**************************************************************************************************

# Here we are going to ancor a couple of charges at the boundary instead of having periodic boundary
# conditions. SecSin1D --> Section Sin 1D
class SecSin1D(object):
    def __init__(self,qq,As,beta,num_cell):
        self.qq = qq
        self.beta = beta
        self.num_cell = num_cell
        # d here is the length of the system
        self.d = num_cell*2.0*pl.pi
        print('slef.d: ' +str(self.d))
        # right now As is only different becasue of different particle "densities". The reason I
        # have stated it like this is because particles with different chages would then need the qq
        # factor to actualy be q[i]*q[j]. or something like that. Lets just see if we can achive the
        # particle separation with the As method
        self.As = As

    def f(self,x,t):
        N = len(x)/2
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

                # Add the forces from the two ancored charges
                # First one at x=0
                temp += self.qq*(x[N+i]-0.0)/(pl.sqrt((x[N+i]-0.0)**2)**3)
                # Second one at the other boundary i.e. self.d
                temp += self.qq*(x[N+i]-self.d)/(pl.sqrt((x[N+i]-self.d)**2)**3)

            # periodic force on particle i
            temp += self.As[i]*pl.sin(x[N+i])*pl.cos(t)
            temp -= self.beta*x[i]
            xdot = pl.append(xdot,temp)
        for i in range(N):
            xdot = pl.append(xdot,x[i])
        return xdot


class O_func(object):
    '''O_func is just a way to acess comon functions across different classes. Or you can import
    this to use them in other ways as well. These function may either take a list or a float or
    pl.ndarray.'''
    def square_wave(self,input):
        if type(input)==float:
            if input%(2.0*pl.pi)<= pl.pi:
                output = 1.0
            if input%(2.0*pl.pi)> pl.pi:
                output = -1.0
        if (type(input)==pl.ndarray)or(type(input)==list):
            output = pl.array([])
            for i,j in enumerate(input):
                if j%(2.0*pl.pi)<= pl.pi:
                    output = pl.append(output,1.0)
                if j%(2.0*pl.pi)> pl.pi:
                    output = pl.append(output,-1.0)
            print('square wave output: ' + str(output)) 
        return output

class Test(object):
    
    def __init__(self,xd,yd):
        self.sigma = 1.0
        self.epsilon = 1.0
        self.order = 2
        self.xd = xd
        self.yd = xd
        # if the distance between two particls is greater than this number than we don't bother with
        # it
        self.cutoff = 3
        # find initial conditions. Put on grid and then displace by half the lattice constant.
        

    def f(self,x,t):
        # for now masses just = 1.0
        # the 4.0 only works for 2D
        N = len(x)/4
        xdot = pl.array([])
        # modulus the y component to keep periodicity right.
        x[3*N:4*N]= x[3*N:4*N]%self.yd
        # x too
        x[2*N:3*N]= x[2*N:3*N]%self.xd

        for i in range(N):
            temp = 0.0
            for j in range(N):
                if i == j:
                    continue
                #repulsive x interparticle force of j on i
                r_temp = pl.sqrt((x[2*N+i]-x[2*N+j])**2+(x[3*N+i]-x[3*N+j])**2)
                if r_temp < self.cutoff:
                    temp += (x[2*N+i]-x[2*N+j])/(r_temp)*24.0*(2.0/(r_temp**13) - 1.0/(r_temp**7))
                # if periodic in x we are going to need to include all the wrap around forces
                for gama in range(self.order):
                    r_temp = pl.sqrt((gama*self.xd-(x[2*N+i]-x[2*N+j]))**2+(x[3*N+i]-x[3*N+j])**2)
                    if r_temp < self.cutoff:
                        temp += -(gama*self.xd-(x[2*N+i]-x[2*N+j]))/(r_temp)*24.0*(2.0/(r_temp**13) - 1.0/(r_temp**7))

                    r_temp = pl.sqrt((gama*self.xd+(x[2*N+i]-x[2*N+j]))**2+(x[3*N+i]-x[3*N+j])**2)
                    if r_temp < self.cutoff:
                        temp += (gama*self.xd+(x[2*N+i]-x[2*N+j]))/(r_temp)*24.0*(2.0/(r_temp**13) - 1.0/(r_temp**7))
            #temp -= self.beta*x[i]
            xdot = pl.append(xdot,temp)
        for i in range(N):
            temp = 0.0
            for j in range(N):
                if i == j:
                    continue
                #repulsive y interparticle force of j on i
                r_temp = pl.sqrt((x[2*N+i]-x[2*N+j])**2+(x[3*N+i]-x[3*N+j])**2)
                if r_temp < self.cutoff:
                    temp += (x[3*N+i]-x[3*N+j])/r_temp*24.0*(2.0/(r_temp**13)-1.0/(r_temp**7)) 
                # force from same particle but in other direction becasue of periodic boundar
                # conditions
                for gama in range(self.order):
                    r_temp = pl.sqrt((x[2*N+i]-x[2*N+j])**2+(gama*self.yd-(x[3*N+i]-x[3*N+j]))**2)
                    if r_temp < self.cutoff:
                        temp += -(gama*self.yd-(x[3*N+i]-x[3*N+j]))/r_temp*24.0*(2.0/(r_temp**13)-1.0/(r_temp**7)) 

                    r_temp = pl.sqrt((x[2*N+i]-x[2*N+j])**2+(gama*self.yd+(x[3*N+i]-x[3*N+j]))**2)
                    if r_temp < self.cutoff:
                        temp += (gama*self.yd+(x[3*N+i]-x[3*N+j]))/r_temp*24.0*(2.0/(r_temp**13)-1.0/(r_temp**7)) 

                    #temp += -self.qq*(gama*self.yd-(x[3*N+i]-x[3*N+j]))/(pl.sqrt((gama*self.yd-(x[3*N+i]-x[3*N+j]))**2+(x[2*N+i]-x[2*N+j])**2)**3)
                    #temp += self.qq* (gama*self.yd+(x[3*N+i]-x[3*N+j]))/(pl.sqrt((gama*self.yd+(x[3*N+i]-x[3*N+j]))**2+(x[2*N+i]-x[2*N+j])**2)**3)
            #temp -= self.beta*x[N+i]
            xdot = pl.append(xdot,temp)
        for i in range(N):
            xdot = pl.append(xdot,x[i])
        for i in range(N):
            xdot = pl.append(xdot,x[N+i])

        return xdot


class SinSin2D(object):
    def __init__(self,qq,As,beta,x_num_cell,y_num_cell,x_periodic,order):
        self.qq = qq
        self.beta = beta
        self.x_num_cell = x_num_cell
        self.y_num_cell = y_num_cell
        self.yd = self.y_num_cell * 2.0*pl.pi
        self.xd = self.x_num_cell * 2.0*pl.pi
        # right now As is only different becasue of different particle "densities". The reason I
        # have stated it like this is because particles with different chages would then need the qq
        # factor to actualy be q[i]*q[j]. or something like that. Lets just see if we can achive the
        # particle separation with the As method
        self.As = As
        # do we want periodicity in x? This is a boolean.
        self.x_periodic = x_periodic
        # the order we want to go to as far as calculating periodic forces
        self.order = order
        self.square_wave = O_func.square_wave

    def f(self,x,t):
        # for now masses just = 1.0
        # the 4.0 only works for 2D
        N = len(x)/4
        xdot = pl.array([])
        # modulus the y component to keep periodicity right.
        x[3*N:4*N]= x[3*N:4*N]%self.yd
        # x too
        if self.x_periodic:
            x[2*N:3*N]= x[2*N:3*N]%self.xd
        for i in range(N):
            temp = 0.0
            for j in range(N):
                if i == j:
                    continue
                #repulsive x interparticle force of j on i
                temp += self.qq*(x[2*N+i]-x[2*N+j])/(pl.sqrt((x[2*N+i]-x[2*N+j])**2+(x[3*N+i]-x[3*N+j])**2)**3)
                # if periodic in x we are going to need to include all the wrap around forces
                if self.x_periodic:
                    for gama in range(self.order):
                        temp += -self.qq*(gama*self.xd-(x[2*N+i]-x[2*N+j]))/(pl.sqrt((gama*self.xd-(x[2*N+i]-x[2*N+j]))**2+(x[3*N+i]-x[3*N+j])**2)**3)
                        temp += self.qq* (gama*self.xd+(x[2*N+i]-x[2*N+j]))/(pl.sqrt((gama*self.xd+(x[2*N+i]-x[2*N+j]))**2+(x[3*N+i]-x[3*N+j])**2)**3)
            # sin(x)cos(t) force on particle i
            temp+=self.As[i]*pl.sin(x[2*N+i])*pl.cos(t)
            temp -= self.beta*x[i]
            xdot = pl.append(xdot,temp)
        for i in range(N):
            temp = 0.0
            for j in range(N):
                if i == j:
                    continue
                #repulsive y interparticle force of j on i
                temp += self.qq*(x[3*N+i]-x[3*N+j])/(pl.sqrt((x[2*N+i]-x[2*N+j])**2+(x[3*N+i]-x[3*N+j])**2)**3)
                # force from same particle but in other direction becasue of periodic boundar
                # conditions
                for gama in range(self.order):
                    temp += -self.qq*(gama*self.yd-(x[3*N+i]-x[3*N+j]))/(pl.sqrt((gama*self.yd-(x[3*N+i]-x[3*N+j]))**2+(x[2*N+i]-x[2*N+j])**2)**3)
                    temp += self.qq* (gama*self.yd+(x[3*N+i]-x[3*N+j]))/(pl.sqrt((gama*self.yd+(x[3*N+i]-x[3*N+j]))**2+(x[2*N+i]-x[2*N+j])**2)**3)
            # sin(y)cos(t) force on particle i
            temp+=self.As[i]*pl.sin(x[3*N+i])*pl.cos(t)
            temp -= self.beta*x[N+i]
            xdot = pl.append(xdot,temp)
        for i in range(N):
            xdot = pl.append(xdot,x[i])
        for i in range(N):
            xdot = pl.append(xdot,x[N+i])

        return xdot

class Sin1D(object):
   
    def __init__(self,qq,As,beta,num_cell):
        self.qq = qq
        self.beta = beta
        self.num_cell = num_cell
        # d here is the length of the system
        self.d = num_cell*2.0*pl.pi
        print('slef.d: ' +str(self.d))
        # right now As is only different becasue of different particle "densities". The reason I
        # have stated it like this is because particles with different chages would then need the qq
        # factor to actualy be q[i]*q[j]. or something like that. Lets just see if we can achive the
        # particle separation with the As method
        self.As = As

    def f(self,x,t):
        N = len(x)/2
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
            temp += self.As[i]*pl.sin(x[N+i])*pl.cos(t)
            temp -= self.beta*x[i]
            xdot = pl.append(xdot,temp)
        for i in range(N):
            xdot = pl.append(xdot,x[i])
        return xdot

class HardCoreSin1D(object):
   
    def __init__(self,rad,As,beta,num_cell,spring):
        # particle radius
        self.rad = rad
        self.beta = beta
        self.num_cell = num_cell
        # d here is the length of the system
        self.d = num_cell*2.0*pl.pi
        print('slef.d: ' +str(self.d))
        # right now As is only different becasue of different particle "densities". The reason I
        # have stated it like this is because particles with different chages would then need the qq
        # factor to actualy be q[i]*q[j]. or something like that. Lets just see if we can achive the
        # particle separation with the As method
        self.As = As
        # spring constant of particles squishieness
        self.spring = spring

    def f(self,x,t):
        N = len(x)/2
        xdot = pl.array([])

        # modulus the x for periodicity.
        x[N:2*N]= x[N:2*N]%self.d
        # HERE ---->> 1Dify
        for i in range(N):
            temp = 0.0
            for j in range(N):
                if i==j: continue
                if abs(x[N+i]-x[N+j])<2.0*self.rad:
                    temp += (x[N+i]-x[N+j])*self.spring

            temp += self.As[i]*pl.sin(x[N+i])*pl.cos(t)
            temp -= self.beta*x[i]
            xdot = pl.append(xdot,temp)
        for i in range(N):
            xdot = pl.append(xdot,x[i])
        return xdot


class CentreLineApx(object):
    def __init__(self,coef,k,w,drgCoef,surf,g,dt):
        self.dt = dt
        self.coef = coef
        self.k = k
        self.w = w
        self.drg = drgCoef
        self.surf = surf
        self.g = g
        self.sol = pl.array([])
        self.square_wave = O_func.square_wave

    def set_sol(self,sol):
        self.sol=sol

   
    # just make normal functions to try to pass into odeint function. Should be much faster
    def f(self,xarr,t):
        temp1 = 0.0
        temp2 = 0.0

        # RIGHT NOW WE ARE LOOKING AT THE SQUARE WAVE VERSION. TO GO BACK TO NORMAL JUST REPLACE
        # SELF.SQUAREWAVE WITH PL.COS. !!!!!!!!!!!!!!!!!!!!!!!!!
    
        for i in range(2):
            temp1+=pl.sin(self.k*xarr[2]-i*pl.pi)*self.square_wave(self.w*t-i*pl.pi)/(pl.cosh(self.k*(self.surf+abs(xarr[3]-self.surf)))-pl.cos(self.k*xarr[2]-i*pl.pi)) 
        temp1 = temp1*self.coef
        temp1 -= self.drg*xarr[0]
        x0dot = temp1
        for i in range(2):
            temp2+=pl.sign(xarr[3]-self.surf)*pl.sinh(self.k*(self.surf+abs(xarr[3]-self.surf)))*self.square_wave(self.w*t-i*pl.pi)/(pl.cosh(self.k*(self.surf+abs(xarr[3]-self.surf)))-pl.cos(self.k*xarr[2]-i*pl.pi)) 
        temp2 = temp2*self.coef
        temp2 -= self.drg*xarr[1]
        temp2 -= pl.sign(xarr[3]-self.surf)*self.g
        x1dot = temp2
        x2dot = xarr[0]
        x3dot = xarr[1]
        return [x0dot,x1dot,x2dot,x3dot]

# A is the dimensionless interaction amplitude
# beta is the dimesionless damping
class TwinECApx(object):
    def __init__(self,interaction_amp,dimless_beta,dimless_dist_between_EC,dt):
        self.dt = dt
        self.A = interaction_amp
        self.beta = dimless_beta
        # distance between the two ECs for now just set to 2pi
        # self.d = dimless_dist_between_EC
        self.d = 2*pl.pi
        self.sol = pl.array([])

    def set_sol(self,sol):
        self.sol=sol
    
    # just make normal functions to try to pass into odeint function. Should be much faster
    def f(self,xarr,t):
        x0dot = self.A*(pl.cos(t)+0.2)*(-2*(pl.cosh(self.d-xarr[3]) + pl.cosh(xarr[3]))*(pl.cos(xarr[2])**2 - pl.cosh(self.d-xarr[3])*pl.cosh(xarr[3]))*pl.sin(xarr[2]))/((pl.cos(xarr[2]) - pl.cosh(self.d-xarr[3]))*(pl.cos(xarr[2]) + pl.cosh(self.d - xarr[3]))*(pl.cos(xarr[2]) - pl.cosh(xarr[3]))*(pl.cos(xarr[2]) + pl.cosh(xarr[3])))- self.beta*xarr[0]
        x1dot = self.A*(pl.cos(t)+0.2)*(2*pl.cos(xarr[2])*(pl.sinh(self.d - xarr[3]) - pl.sinh(xarr[3]))*(-pl.sin(xarr[2])**2 + pl.sinh(self.d - xarr[3])*pl.sinh(xarr[3])))/((pl.cos(xarr[2]) - pl.cosh(self.d - xarr[3]))*(pl.cos(xarr[2]) + pl.cosh(self.d - xarr[3]))*(pl.cos(xarr[2]) - pl.cosh(xarr[3]))*(pl.cos(xarr[2]) + pl.cosh(xarr[3])))- self.beta*xarr[1]
        x2dot = xarr[0]
        x3dot = xarr[1]
        return [x0dot,x1dot,x2dot,x3dot]

class MultiTwinECApxQuad(object):
    def __init__(self,qq,A,beta,d,quadA):
        self.qq = qq
        self.A = A
        self.quadA = quadA
        self.d = d
        self.beta = beta

    def f(self,x,t):
    
        # the 4.0 only works for 2D
        N = len(x)/4
        xdot = pl.array([])

        print('x is: '+str(x))
    
        for i in range(N):
            temp = 0.0
            for j in range(N):
                if i == j:
                    continue
                #repulsive x interparticle force of j on i
                temp += self.qq*(x[2*N+i]-x[2*N+j])/(pl.sqrt((x[2*N+i]-x[2*N+j])**2+(x[3*N+i]-x[3*N+j])**2)**3)
            # EC x force on particle i
            temp += self.A*(pl.cos(t))*(-2*(pl.cosh(self.d-x[3*N+i]) + pl.cosh(x[3*N+i]))*(pl.cos(x[2*N+i])**2 - pl.cosh(self.d-x[3*N+i])*pl.cosh(x[3*N+i]))*pl.sin(x[2*N+i]))/((pl.cos(x[2*N+i]) - pl.cosh(self.d-x[3*N+i]))*(pl.cos(x[2*N+i]) + pl.cosh(self.d - x[3*N+i]))*(pl.cos(x[2*N+i]) - pl.cosh(x[3*N+i]))*(pl.cos(x[2*N+i]) + pl.cosh(x[3*N+i]))) -self.beta*x[i]
            xdot = pl.append(xdot,temp)
        for i in range(N):
            temp = 0.0
            for j in range(N):
                if i == j:
                    continue
                #repulsive y interparticle force of j on i
                temp += self.qq*(x[3*N+i]-x[3*N+j])/(pl.sqrt((x[2*N+i]-x[2*N+j])**2+(x[3*N+i]-x[3*N+j])**2)**3)
            # EC y force on particle i
            temp += self.A*(pl.cos(t))*(2*pl.cos(x[2*N+i])*(pl.sinh(self.d - x[3*N+i]) - pl.sinh(x[3*N+i]))*(-pl.sin(x[2*N+i])**2 + pl.sinh(self.d - x[3*N+i])*pl.sinh(x[3*N+i])))/((pl.cos(x[2*N+i]) - pl.cosh(self.d - x[3*N+i]))*(pl.cos(x[2*N+i]) + pl.cosh(self.d - x[3*N+i]))*(pl.cos(x[2*N+i]) - pl.cosh(x[3*N+i]))*(pl.cos(x[2*N+i]) + pl.cosh(x[3*N+i])))-self.beta*x[N+i]
            # quadropole restoring force to center of electric curtains
            temp += -self.quadA*(pl.cos(50.0*t)+1.0)*(x[3*N+i]-self.d/2.0)
            xdot = pl.append(xdot,temp)
        for i in range(N):
            xdot = pl.append(xdot,x[i]) 
        for i in range(N):
            xdot = pl.append(xdot,x[N+i])
    
        print('len xdot is: '+str(len(xdot)))
        print('xdot is: '+str(xdot))
        return xdot

class MultiTwinECApxQuadPeriodic(object):
    def __init__(self,qq,A,beta,d,quadA,num_cell):
        self.qq = qq
        self.A = A
        self.quadA = quadA
        self.d = d
        self.beta = beta
        self.num_cell = num_cell
        self.diameter = num_cell*2.0*pl.pi
        self.square_wave = O_func.square_wave

    def f(self,x,t):
        # the 4.0 only works for 2D
        N = len(x)/4
        xdot = pl.array([])
        # modulus the x component to keep periodicity right.
        x[2*N:3*N]= x[2*N:3*N]%self.diameter

        for i in range(N):
            temp = 0.0
            for j in range(N):
                if i == j:
                    continue
                #repulsive x interparticle force of j on i
                temp += self.qq*(x[2*N+i]-x[2*N+j])/(pl.sqrt((x[2*N+i]-x[2*N+j])**2+(x[3*N+i]-x[3*N+j])**2)**3)
                # force from same particle but in other direction becasue of periodic boundar
                # conditions
                temp += -self.qq*(self.diameter-(x[2*N+i]-x[2*N+j]))/(pl.sqrt((self.diameter-(x[2*N+i]-x[2*N+j]))**2+(x[3*N+i]-x[3*N+j])**2)**3)
                # could do this forever but for now just include one more of force going around in
                # same direction as first. draw a diagam if you need to
                temp += self.qq*(self.diameter+(x[2*N+i]-x[2*N+j]))/(pl.sqrt((self.diameter+(x[2*N+i]-x[2*N+j]))**2+(x[3*N+i]-x[3*N+j])**2)**3)
            # EC x force on particle i
            temp += self.A*(self.square_wave(t))*(-2*(pl.cosh(self.d-x[3*N+i]) + \
                pl.cosh(x[3*N+i]))*(pl.cos(x[2*N+i])**2 - \
                    pl.cosh(self.d-x[3*N+i])*pl.cosh(x[3*N+i]))*pl.sin(x[2*N+i]))/((pl.cos(x[2*N+i])\
                        - pl.cosh(self.d-x[3*N+i]))*(pl.cos(x[2*N+i]) + pl.cosh(self.d - \
                            x[3*N+i]))*(pl.cos(x[2*N+i]) - pl.cosh(x[3*N+i]))*(pl.cos(x[2*N+i]) + \
                                pl.cosh(x[3*N+i]))) - self.beta*x[i]
            xdot = pl.append(xdot,temp)
        for i in range(N):
            temp = 0.0
            for j in range(N):
                if i == j:
                    continue
                #repulsive y interparticle force of j on i
                temp += self.qq*(x[3*N+i]-x[3*N+j])/(pl.sqrt((x[2*N+i]-x[2*N+j])**2+(x[3*N+i]-x[3*N+j])**2)**3)
                # force from same particle but in other direction becasue of periodic boundar
                # conditions
                temp += self.qq*(x[3*N+i]-x[3*N+j])/(pl.sqrt((self.diameter-(x[2*N+i]-x[2*N+j]))**2+(x[3*N+i]-x[3*N+j])**2)**3)
            # EC y force on particle i
            temp += self.A*(self.square_wave(t))*(2*pl.cos(x[2*N+i])*(pl.sinh(self.d - x[3*N+i]) - \
                pl.sinh(x[3*N+i]))*(-pl.sin(x[2*N+i])**2 + pl.sinh(self.d - \
                    x[3*N+i])*pl.sinh(x[3*N+i])))/((pl.cos(x[2*N+i]) - pl.cosh(self.d - \
                        x[3*N+i]))*(pl.cos(x[2*N+i]) + pl.cosh(self.d - x[3*N+i]))*(pl.cos(x[2*N+i]) \
                            - pl.cosh(x[3*N+i]))*(pl.cos(x[2*N+i]) + \
                                pl.cosh(x[3*N+i])))-self.beta*x[N+i] 
            # quadropole restoring force to center of electric curtains
            temp += -self.quadA*(pl.cos(200.0*t)+0.1)*(x[3*N+i]-self.d/2.0)
            xdot = pl.append(xdot,temp)
        for i in range(N):
            xdot = pl.append(xdot,x[i]) 
        for i in range(N):
            xdot = pl.append(xdot,x[N+i])
        return xdot

# This only works for 1D
class ExactPeriodic1D(object):
   
    def __init__(self,qq,As,beta,num_cell):
        self.qq = qq
        self.beta = beta
        self.num_cell = num_cell
        # d here is the length of the system
        self.d = num_cell*2.0*pl.pi
        # right now As is only different becasue of different particle "densities". The reason I
        # have stated it like this is because particles with different chages would then need the qq
        # factor to actualy be q[i]*q[j]. or something like that. Lets just see if we can achive the
        # particle separation with the As method
        self.As = As
        self.square_wave = O_func.square_wave

    def f(self,x,t):
        N = len(x)/2
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
                temp += self.qq*(polygamma(1,(self.d+x[N+i]-x[N+j])/self.d)-polygamma(1,1.0-((x[N+i]-x[N+j])/self.d)))/(self.d**2)
            # EC x force on particle i
            for a in range(2):
                temp+=self.As[i]*pl.sin(x[N+i]-a*pl.pi)*pl.cos(t-a*pl.pi)/(pl.cosh(1.0)-pl.cos(x[N+i]-a*pl.pi)) 
            temp -= self.beta*x[i]
            xdot = pl.append(xdot,temp)

        for i in range(N):
            xdot = pl.append(xdot,x[i])

    
        return xdot

class Flat2D(object):
    def __init__(self,qq,As,beta,yd,num_cell,x_periodic,order):
        self.qq = qq
        self.beta = beta
        self.num_cell = num_cell
        self.yd = yd
        self.xd = num_cell * 2.0*pl.pi
        # right now As is only different becasue of different particle "densities". The reason I
        # have stated it like this is because particles with different chages would then need the qq
        # factor to actualy be q[i]*q[j]. or something like that. Lets just see if we can achive the
        # particle separation with the As method
        self.As = As
        # do we want periodicity in x? This is a boolean.
        self.x_periodic = x_periodic
        # the order we want to go to as far as calculating periodic forces
        self.order = order
        self.square_wave = O_func.square_wave

    def f(self,x,t):
        # for now masses just = 1.0
        # the 4.0 only works for 2D
        N = len(x)/4
        xdot = pl.array([])
        # modulus the y component to keep periodicity right.
        x[3*N:4*N]= x[3*N:4*N]%self.yd
        # x too
        if self.x_periodic:
            x[2*N:3*N]= x[2*N:3*N]%self.xd
        for i in range(N):
            temp = 0.0
            for j in range(N):
                if i == j:
                    continue
                #repulsive x interparticle force of j on i
                temp += self.qq*(x[2*N+i]-x[2*N+j])/(pl.sqrt((x[2*N+i]-x[2*N+j])**2+(x[3*N+i]-x[3*N+j])**2)**3)
                # if periodic in x we are going to need to include all the wrap around forces
                if self.x_periodic:
                    #repulsive y interparticle force of j on i
                    temp += self.qq*(x[2*N+i]-x[2*N+j])/(pl.sqrt((x[2*N+i]-x[2*N+j])**2+(x[3*N+i]-x[3*N+j])**2)**3)
                    for gama in range(self.order):
                        temp += -self.qq*(gama*self.xd-(x[2*N+i]-x[2*N+j]))/(pl.sqrt((gama*self.xd-(x[2*N+i]-x[2*N+j]))**2+(x[3*N+i]-x[3*N+j])**2)**3)
                        temp += self.qq* (gama*self.xd+(x[2*N+i]-x[2*N+j]))/(pl.sqrt((gama*self.xd+(x[2*N+i]-x[2*N+j]))**2+(x[3*N+i]-x[3*N+j])**2)**3)
            # EC x force on particle i
            # surface set to 1.0
            for a in range(2):
                temp+=self.As[i]*pl.sin(x[2*N+i]-a*pl.pi)*pl.cos(t-a*pl.pi)/(pl.cosh(1.0)-pl.cos(x[2*N+i]-a*pl.pi)) 
            temp -= self.beta*x[i]
            xdot = pl.append(xdot,temp)
        for i in range(N):
            temp = 0.0
            for j in range(N):
                if i == j:
                    continue
                #repulsive y interparticle force of j on i
                temp += self.qq*(x[3*N+i]-x[3*N+j])/(pl.sqrt((x[2*N+i]-x[2*N+j])**2+(x[3*N+i]-x[3*N+j])**2)**3)
                # force from same particle but in other direction becasue of periodic boundar
                # conditions
                for gama in range(self.order):
                    temp += -self.qq*(gama*self.yd-(x[3*N+i]-x[3*N+j]))/(pl.sqrt((gama*self.yd-(x[3*N+i]-x[3*N+j]))**2+(x[2*N+i]-x[2*N+j])**2)**3)
                    temp += self.qq* (gama*self.yd+(x[3*N+i]-x[3*N+j]))/(pl.sqrt((gama*self.yd+(x[3*N+i]-x[3*N+j]))**2+(x[2*N+i]-x[2*N+j])**2)**3)
            temp -= self.beta*x[N+i]
            xdot = pl.append(xdot,temp)
        for i in range(N):
            xdot = pl.append(xdot,x[i])
        for i in range(N):
            xdot = pl.append(xdot,x[N+i])

    
        return xdot

# for now I am going to include the "new type" of reflecting surface...This means I will switch the
# sign of the y interactions but let the particle go below the plain of the surface. This is an
# atempt at dealing with the discontinuity in the y velocity upon a reflection.

class travPlaneChrgApx(object):
    # k should be less than 4189(inverse meters)
    # coef sould be less than 333
    def __init__(self,coef,k,w,drgCoef):
        self.coef = coef
        self.k = k
        self.w = w
        self.drg = drgCoef
        # for now keep the surface at
    def f(self,xarr,t):
        x0dot = -self.coef*pl.exp(-self.k*(1.0+abs(xarr[3]-1.0)))*pl.sin(self.w*t-self.k*xarr[2])
        x1dot = pl.sign(xarr[3]-1.0)*self.coef*pl.exp(-self.k*(1.0+abs(xarr[3]-1.0)))*pl.cos(self.w*t-self.k*xarr[2]) -\
                pl.sign(xarr[3]-1.0)*9.8
        x2dot = xarr[0]
        x3dot = xarr[1]
        return [x0dot,x1dot,x2dot,x3dot]

class surfTravPlaneChrgApx(object):
    # k should be less than 4189(inverse meters)
    # coef sould be less than 333
    def __init__(self,coef,k,w,drgCoef):
        self.coef = coef
        self.k = k
        self.w = w
        self.drg = drgCoef
    def f(self,xarr,t):
        x0dot =  -self.coef*pl.exp(-self.k*abs(xarr[3]))*pl.sin(self.w*t-self.k*xarr[2])
        x1dot = 0
        x2dot = xarr[0]
        x3dot = 0
        return [x0dot,x1dot,x2dot,x3dot]

class surfCentreLineApx(object):
    def __init__(self,coef,k,w,drgCoef):
        self.coef = coef
        self.k = k
        self.w = w
        self.drg = drgCoef
    
    # just make normal functions to try to pass into odeint function. Should be much faster
    def f(self,xarr,t):
        temp = 0.0
    
        for i in range(2):
            temp+=pl.sin(self.k*xarr[2]-i*pl.pi)*pl.cos(self.w*t-i*pl.pi)/(pl.cosh(self.k*xarr[3])-pl.cos(self.k*xarr[2]-i*pl.pi)) 
        temp = temp*self.coef
        temp -= self.drg*xarr[0]
        x0dot = temp
        x1dot = 0.0
        x2dot = xarr[0]
        x3dot = 0.0
        return [x0dot,x1dot,x2dot,x3dot]
class One_Particle_Ensble_Sin1D(object):
    def __init__(self,A,beta):
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
   
    # just make normal functions to try to pass into odeint function. Should be much faster
    def f(self,xarr,t):
        x0dot = self.A*pl.sin(xarr[2])*pl.cos(t) - self.beta*xarr[0]
        x1dot = 0.0
        x2dot = xarr[0]
        x3dot = 0.0
        return [x0dot,x1dot,x2dot,x3dot]

class One_Particle_Ensble_Sin2D(object):
    def __init__(self,A,beta):
        self.A = A
        print('self.A is: ' + str(self.A))
        self.beta = beta
        print('self.beta is: '+str(self.beta))
        self.x_num_cell = 1.0
        self.y_num_cell = 1.0
        # d here is the length of the system
        self.dx = self.x_num_cell*2.0*pl.pi
        self.dy = self.y_num_cell*2.0*pl.pi
        print('slef.dx: ' +str(self.dx))
        print('slef.dy: ' +str(self.dy))
        # right now As is only different becasue of different particle "densities". The reason I
        # have stated it like this is because particles with different chages would then need the qq
        # factor to actualy be q[i]*q[j]. or something like that. Lets just see if we can achive the
        # particle separation with the As method
        self.As = A
   
    # just make normal functions to try to pass into odeint function. Should be much faster
    def f(self,xarr,t):
        x0dot = self.A*(pl.sin(xarr[2]))*pl.cos(t) - self.beta*xarr[0]
        x1dot = self.A*(pl.sin(xarr[3]))*pl.cos(t) - self.beta*xarr[1]
        x2dot = xarr[0]
        x3dot = xarr[1]
        return [x0dot,x1dot,x2dot,x3dot]

