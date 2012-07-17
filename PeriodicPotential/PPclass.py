import pylab as pl
from datetime import datetime
import os


# ---------------------------------------------------------------------------------------------

# This is going to be our class for periodic potentials

# using the scipy odeint functions as a way of solving first order differential equations can be a
# little confusing so here are some notes:
# first of all the new array structure being passed back and forth for a full location in phase
# space is x1dot = x doubledot (or vx dot)
#          x2dot = x
# this is the form of the array "xarr" that is being passed into the the f(self,xarr,t) function.
# the "t" in this case is the time array that we want our solution for. Find the youtube guy if you
# forget what goes where. (nonlinear pendulub ex)

#   Explicitly the equations of motion are xdot = v. vdot = A*cos(omega*time)sin(k*x). We will just
#   express them in their dimensionless form though with vdot = Acos(t')sin(x'). When we include
#   code for a square and Hexagonal potiential it will be the equvelent set up.


class periodicPotential(object):
    def __init__(self,coef,drgCoef):
        self.coef = coef
        self.drg = drgCoef
    
    # just make normal functions to try to pass into odeint function. Should be much faster
    def f(self,xarr,t):
        x1dot = -self.drg*xarr[0] + self.coef*pl.cos(t)*pl.sin(xarr[1])
        x2dot = xarr[0]
        return [x1dot,x2dot]

# Start with a simpe symetric periodic potential. Just to be clear the particle is in the x,y plane
# and the potential is in the x,y plane. This is not a "hopping of the surface" type problem.

# Here we are looking at the potential (cosx + cosy)*cost
class simple2DSymetricPP(object):
    def __init__(self,coef,drgCoef): 
        self.coef = coef
        self.drg = drgCoef
    # The structure of this array with the position and velocity information will be:
    # [xdot,ydot,x,y]
    def f(self,arr,t): 
        x1dot = -self.drg*arr[0] + self.coef*pl.cos(t)*pl.sin(arr[2])
        x2dot = -self.drg*arr[1] + self.coef*pl.cos(t)*pl.sin(arr[3])
        x3dot = arr[0]
        x4dot = arr[1]
        return [x1dot,x2dot,x3dot,x4dot]

# This is more or less the same as the class above but here we look at a potential of
# (cosx*cosy)*cost
class secondSimple2DSymetricPP(object):
    def __init__(self,coef,drgCoef): 
        self.coef = coef
        self.drg = drgCoef
    # The structure of this array with the position and velocity information will be:
    # [xdot,ydot,x,y]
    def f(self,arr,t): 
        x1dot = -self.drg*arr[0] + self.coef*pl.cos(t)*pl.sin(arr[2])*pl.cos(arr[3])
        x2dot = -self.drg*arr[1] + self.coef*pl.cos(t)*pl.sin(arr[3])*pl.cos(arr[2])
        x3dot = arr[0]
        x4dot = arr[1]
        return [x1dot,x2dot,x3dot,x4dot]

