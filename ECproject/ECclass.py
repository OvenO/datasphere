#import pylab as pl
from datetime import datetime
import sys
sys.path.append("/users/o/m/omyers/puthere/lib/python2.7/site-packages")
import scipy as pl
import os

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

# for now I am going to include the "new type" of reflecting surface...This means I will switch the
# sign of the y interactions but let the particle go below the plain of the surface. This is an
# atempt at dealing with the discontinuity in the y velocity upon a reflection.

# This is (at first) going to be used as a "normalized" equation of motion which in this case only means our lengths are
# normalized (specificaly so surface can just be at 1 like in out 1D case). gravity is going to be the same old 9.8 and the multiplying coeficient is therefor
# going to be whatever strange number makes things look nice hahah (Eeeee...) 
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

    def set_sol(self,sol):
        self.sol=sol

    def square_wave(input):
        output = pl.array([])
        for i,j in enumerate(input):
            if j%(2.0*pl.pi)<= pl.pi:
                output = pl.append(output,1.0)
            if j%(2.0*pl.pi)> pl.pi:
                output = pl.append(output,-1.0)
                    
        return output
    
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

