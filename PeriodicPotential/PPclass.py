#import pylab as pl
from datetime import datetime
import sys
sys.path.append("/users/o/m/omyers/puthere/lib/python2.7/site-packages")
import scipy as pl
import os


class timeDependentPotential1D(object):
    def __init__(self,coef,drgCoef):
        self.coef = coef
        self.drg = drgCoef
    
    # just make normal functions to try to pass into odeint function. Should be much faster
    def f(self,xarr,t):
    
        x1dot = xarr[1]
        x1dot=self.coef*pl.sin(xarr[1])*pl.cos(t) - self.drg*xarr[0]
        return [x1dot,x2dot]
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

