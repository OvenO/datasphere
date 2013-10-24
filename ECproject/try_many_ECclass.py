import pylab as pl
import random
import ECclass as ec
from scipy.integrate import odeint


def main():
    # number of particles
    N=8
    qq = 0.5
    A  = 0.4
    beta = .2
    quadA = 40.0
    # distance between electrodes
    d = pl.pi
    num_cell =2.0
    diameter = num_cell*2.0*pl.pi

    x0 = pl.zeros([4*N])

    some_colors = ['b','c','g','k','r','y','b','c','g','k','r','y','b','c','g','k','r','y']

    for i,j in enumerate(x0):
        if i in range(2*N,3*N):
            x0[i] = random.random()*diameter
            continue
        x0[i] = random.random()*d

    t = pl.arange(0,30,.1)

    elec = ec.MultiTwinECApxQuadPeriodic(qq,A,beta,d,quadA,num_cell)
    
    # first just try to recreate what we did with 2 particles
    #x0[0] = 0.0
    #x0[1] = 0.0
    #x0[2] = 0.0
    #x0[3] = 0.0
    #x0[4] = 1.0
    #x0[5] = 1.5
    #x0[6] = 1.4
    #x0[7] = 1.0

    sol = odeint(elec.f,x0,t)

    for i in range(len(t)):
        fig = pl.figure()
        ax = fig.add_subplot(111)
        ax.set_xlim([0.0,diameter])
        ax.set_ylim([-0.5,pl.pi+0.5])
        ax.scatter(sol[i,2*N:3*N]%diameter,sol[i,3*N:4*N],c=some_colors)
        ax.scatter(pl.arange(0,diameter+pl.pi,pl.pi),pl.zeros(int(diameter/pl.pi+1.5)),s=5)
        ax.scatter(pl.arange(0,diameter+pl.pi,pl.pi),pl.zeros(int(diameter/pl.pi+1.5))+pl.pi,s=5)
        print('len scatter array is: '+str(len(sol[i,2*N:3*N])))
        #ax.scatter(sol[i,5],sol[i,7])
        fig.savefig('%(number)04d.png'%{'number':i})
        pl.close(fig)

if __name__ == '__main__':
        main()
