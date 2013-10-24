import pylab as pl
import random
from scipy.integrate import odeint

def f(x,t):
    # for now masses just = 1.0

    # the 4.0 only works for 2D
    N = len(x)/4
    xdot = pl.array([])

    for i in range(N):
        temp = 0.0
        for j in range(N):
            if i == j:
                continue
            temp += -(x[2*N+i]-x[2*N+j])/(pl.sqrt((x[2*N+i]-x[2*N+j])**2+(x[3*N+i]-x[3*N+j])**2)**3)
        xdot = pl.append(xdot,temp)
    for i in range(N):
        temp = 0.0
        for j in range(N):
            if i == j:
                continue
            temp += -(x[3*N+i]-x[3*N+j])/(pl.sqrt((x[2*N+i]-x[2*N+j])**2+(x[3*N+i]-x[3*N+j])**2)**3)
        xdot = pl.append(xdot,temp)
    for i in range(N):
        xdot = pl.append(xdot,x[i]) 
    for i in range(N):
        xdot = pl.append(xdot,x[N+i])

    print('len xdot is: '+str(len(xdot)))
    return xdot

def main():
    # number of particles
    N=8
    x0 = pl.zeros([4*N])

    for i,j in enumerate(x0):
        x0[i] = random.random()*2
    
    # first just try to recreate what we did with 2 particles
    #x0[0] = 0.0
    #x0[1] = 0.0
    #x0[2] = 0.0
    #x0[3] = 1.0
    #x0[4] = 0.0
    #x0[5] = .2
    #x0[6] = 0.0
    #x0[7] = 0.0

    t = pl.arange(0,10,.01)

    sol = odeint(f,x0,t)

    for i in range(len(t)):
        fig = pl.figure()
        ax = fig.add_subplot(111)
        ax.set_xlim([-1,5])
        ax.set_ylim([-1,5])
        ax.scatter(sol[i,2*N:3*N],sol[i,3*N:4*N])
        print('len scatter array is: '+str(len(sol[i,2*N:3*N])))
        #ax.scatter(sol[i,5],sol[i,7])
        fig.savefig('%(number)04d.png'%{'number':i})
        pl.close(fig)

if __name__ == '__main__':
        main()
