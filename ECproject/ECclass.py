#import pylab as pl
from datetime import datetime
import sys
sys.path.append("/users/o/m/omyers/puthere/lib/python2.7/site-packages")
import scipy as pl
import os


# ---------------------------------------------------------------------------------------------

# This is a remake of the ElectricCurtain class that will run faster and without the problems from
# the previous. 
#
# coeff ---> this is the coefficient multiplying the acceleration function (in the case of the clc
# this does not include k). The charge to mass ratioo and the amplitude of the electric field or
# charge are wraped up in here as well as the permitivity.
#
# k ---> wave number (2*pi/wavelength)
# w ---> omega, The driving frequency
# stokeDrg ---> stokes drag coefficient
#
# Another important feture of this class is that the gravitational constand is set to the actual
# value (9.8) but it can be changed if one want to run the electric curtain on a different planet
# (or no planet)
#
# Also we have not made parameters dimensionless and THAT MUST BE DONE STILL

class electricCurtain(object):

    # Make The functions for Defining the Block and phase space location vectors
    # first set where the lower left corner is
    def setCorner(self,x,xdot,y,ydot):
        self.corner = pl.array([float(x),float(xdot),float(y),float(ydot)])
    # next set the dimensions of the block in phase space. ie the total length of each side.
    def setDims(self,xlen,xdlen,ylen,ydlen):
        self.xlen = float(xlen)
        self.xdlen = float(xdlen)
        self.ylen = float(ylen)
        self.ydlen = float(ydlen)
    # need to set the number of points per dimension
    def setnumPts(self,x,xdot,y,ydot):
        self.numx = x
        self.numxdot = xdot
        self.numy = y
        self.numydot = ydot

        self.totNumPts = x*y*xdot*ydot
    # we now have all the information to make the block
    def makeBlock(self):
        self.Block = pl.zeros([self.totNumPts +1 ,4])

        # this array holds the spacing between points in the directions x,xdot,y,ydot respectivly
        self.stepSizeArray = pl.array([self.xlen/self.numx , self.xdlen/self.numxdot ,  \
                self.ylen/self.numy , self.ydlen/self.numydot])

        # currently this loop fills the "Block" array with vectors pointing to all of the points in
        # the block. It is worth noting that the edges furthes from the corner do not contain points
        # directly on the edge. This is an easy fix if it needs to be but as the number of points in
        # the block gets large the outside edge is infentesimaly close to the block.
        index = 0
        self.Block[index] = self.corner
        for i in range(self.numx):
            for j in range(self.numxdot):
                for l in range(self.numy):
                    for k in range(self.numydot):
                        index += 1
                        self.Block[index,0] = self.Block[0,0] + \
                                i*self.stepSizeArray[0]
                        self.Block[index,1] = self.Block[0,1] + \
                                j*self.stepSizeArray[1]
                        self.Block[index,2] = self.Block[0,2] + \
                                l*self.stepSizeArray[2]
                        self.Block[index,3] = self.Block[0,3] + \
                                k*self.stepSizeArray[3]
# using the scipy odeint functions as a way of solving first order differential equations can be a
# little confusing so here are some notes:
# first of all the new array structure being passed back and forth for a full location in phase
# space is x1dot = x doubledot (or vx dot)
#          x2dot = y doubledot (or vy dot)
#          x3dot = x
#          x4dot = y                                
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
        x1dot = -self.coef*pl.exp(-self.k*(1.0+abs(xarr[3]-1.0)))*pl.sin(self.w*t-self.k*xarr[2])
        x2dot = pl.sign(xarr[3]-1.0)*self.coef*pl.exp(-self.k*(1.0+abs(xarr[3]-1.0)))*pl.cos(self.w*t-self.k*xarr[2]) -\
                pl.sign(xarr[3]-1.0)*9.8
        x3dot = xarr[0]
        x4dot = xarr[1]
        return [x1dot,x2dot,x3dot,x4dot]

class surfTravPlaneChrgApx(object):
    # k should be less than 4189(inverse meters)
    # coef sould be less than 333
    def __init__(self,coef,k,w,drgCoef):
        self.coef = coef
        self.k = k
        self.w = w
        self.drg = drgCoef
    def f(self,xarr,t):
        x1dot =  -self.coef*pl.exp(-self.k*abs(xarr[3]))*pl.sin(self.w*t-self.k*xarr[2])
        x2dot = 0
        x3dot = xarr[0]
        x4dot = 0
        return [x1dot,x2dot,x3dot,x4dot]

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
        x1dot = temp
        x2dot = 0.0
        x3dot = xarr[0]
        x4dot = 0.0
        return [x1dot,x2dot,x3dot,x4dot]
class CentreLineApx(object):
    def __init__(self,coef,k,w,drgCoef,surf,g):
        self.coef = coef
        self.k = k
        self.w = w
        self.drg = drgCoef
        self.surf = surf
        self.g = g
    
    # just make normal functions to try to pass into odeint function. Should be much faster
    def f(self,xarr,t):
        temp1 = 0.0
        temp2 = 0.0
    
        for i in range(2):
            temp1+=pl.sin(self.k*xarr[2]-i*pl.pi)*pl.cos(self.w*t-i*pl.pi)/(pl.cosh(self.k*(self.surf+abs(xarr[3]-self.surf)))-pl.cos(self.k*xarr[2]-i*pl.pi)) 
        temp1 = temp1*self.coef
        temp1 -= self.drg*xarr[0]
        x1dot = temp1
        for i in range(2):
            temp2+=pl.sign(xarr[3]-self.surf)*pl.sinh(self.k*(self.surf+abs(xarr[3]-self.surf)))*pl.cos(self.w*t-i*pl.pi)/(pl.cosh(self.k*(self.surf+abs(xarr[3]-self.surf)))-pl.cos(self.k*xarr[2]-i*pl.pi)) 
        temp2 = temp2*self.coef
        temp2 -= self.drg*xarr[1]
        temp2 -= pl.sign(xarr[3]-self.surf)*self.g
        x2dot = temp2
        x3dot = xarr[0]
        x4dot = xarr[1]
        return [x1dot,x2dot,x3dot,x4dot]

