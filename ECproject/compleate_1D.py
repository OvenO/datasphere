import pylab as pl
import random
import ECclass as ec
from scipy.integrate import odeint
import os
import shutil
  

def main():
    # do we want phase space imaes
    make_phase = True
    # do we want real space images
    make_real = True
    # number of particles
    N=8
    qq = 1.0
    beta = 1.0
    num_cell =5.0
    d = num_cell * 2.0 * pl.pi

    As = pl.zeros(2*N) + .5
    #As = pl.zeros(N/2)+.8
    #As = pl.append(As,pl.zeros(N/2)+.4)

    x0 = pl.zeros([2*N])

    some_colors = ['b','c','g','k','r','y','b','c','g','k','r','y','b','c','g','k','r','y']

    for i,j in enumerate(x0):
        if i in range(N,2*N):
            print(i)
            x0[i] = random.random()*d
            continue

    t = pl.arange(0,70,.07)
    
    elec = ec.ExactPeriodic1D(qq,As,beta,num_cell)
    
    ## first just try to recreate what we did with 2 particles
    #x0[0] = 0.0
    #x0[1] = 0.0
    #x0[2] = 1.0
    #x0[3] = 2.0

    sol = odeint(elec.f,x0,t)

    os.mkdir('RunImages')
    if make_phase:
        os.mkdir('RunImages/PhaseSpace')
    if make_real:
        os.mkdir('RunImages/RealSpace')

    for i in range(len(t)):
        print(sol[i,:N])
        # scatter for different As solutions
        #ax.scatter(sol[i,2*N:(2*N+N/2)],sol[i,3*N:(3*N+N/2)]%d,c='r')
        #ax.scatter(sol[i,(2*N+N/2):3*N],sol[i,(3*N+N/2):4*N]%d,c='b')
        if make_real:
            r_fig = pl.figure()
            r_ax = r_fig.add_subplot(111)
            r_ax.set_xlim([0,d])
            r_ax.set_ylim([0,2.0])

            r_ax.scatter(sol[i,N:2*N]%d,pl.zeros(N)+1.0,c='b')
            r_fig.savefig('RunImages/RealSpace/%(number)04d.png'%{'number':i})
            pl.close(r_fig)
        if make_phase:
            p_fig = pl.figure()
            p_ax = p_fig.add_subplot(111)
            p_ax.set_xlim([0,d])
            p_ax.set_ylim([-2.0,2.0])
            p_ax.scatter(sol[i,N:2*N]%d,sol[i,:N],c='b')
            p_fig.savefig('RunImages/PhaseSpace/%(number)04d.png'%{'number':i})
            pl.close(p_fig)

        #ax.scatter(pl.arange(0,diameter+pl.pi,pl.pi),pl.zeros(int(diameter/pl.pi+1.5)),s=5)
        #ax.scatter(pl.arange(0,diameter+pl.pi,pl.pi),pl.zeros(int(diameter/pl.pi+1.5))+pl.pi,s=5)
        #print('len scatter array is: '+str(len(sol[i,2*N:3*N])))
        #ax.scatter(sol[i,5],sol[i,7])

    # only works for two particles but this verifies that the 'analitic' periodic potential is
    # working.
    #print('final distanc')
    #print(sol[-1,2]-sol[-1,3])

    # be carful with this
    #if 'RunImages' in os.listdir('.'):
    #    shutil.rmtree('RunImages')

if __name__ == '__main__':
        main()
