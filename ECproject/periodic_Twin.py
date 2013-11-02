import pylab as pl
import random
import ECclass as ec
from scipy.integrate import odeint


def main():
    # number of particles
    N=6
    qq = 0.5
    A  = 1.1
    # A = 0.9 and beta = 1.0 is not bad
    beta = 2.4
    quadA = 40.0
    # distance between electrodes
    d = pl.pi
    num_cell =3.0
    diameter = num_cell*2.0*pl.pi

    x0 = pl.zeros([4*N])

    some_colors = ['b','c','g','k','r','y','b','c','g','k','r','y','b','c','g','k','r','y']

    #for i,j in enumerate(x0):
    #    if i in range(2*N,4*N):
    #        x0[i] = random.random()*diameter
    #        continue
    #    x0[i] = random.random()*d
    x0[12]=1.0*pl.pi-pl.pi/2.0
    x0[13]=1.0*pl.pi+pl.pi/2.0
    x0[14]=3.0*pl.pi-pl.pi/2.0
    x0[15]=3.0*pl.pi+pl.pi/2.0
    x0[16]=5.0*pl.pi-pl.pi/2.0
    x0[17]=5.0*pl.pi+pl.pi/2.0
    x0[18]=pl.pi/2.0
    x0[19]=pl.pi/2.0
    x0[20]=pl.pi/2.0
    x0[21]=pl.pi/2.0
    x0[22]=pl.pi/2.0
    x0[23]=pl.pi/2.0

    t = pl.arange(0,100,.1)

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
        ax.set_ylim([-1.0,pl.pi+1.0])
        ax.scatter(sol[i,2*N:3*N]%diameter,sol[i,3*N:4*N],c=some_colors)
        ax.scatter(pl.arange(0,diameter+pl.pi,pl.pi),pl.zeros(int(diameter/pl.pi+1.5)),s=50)
        ax.scatter(pl.arange(0,diameter+pl.pi,pl.pi),pl.zeros(int(diameter/pl.pi+1.5))+pl.pi,s=50)
        print('len scatter array is: '+str(len(sol[i,2*N:3*N])))
        #ax.scatter(sol[i,5],sol[i,7])
        fig.savefig('%(number)04d.png'%{'number':i})
        pl.close(fig)

if __name__ == '__main__':
        main()
