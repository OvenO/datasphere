import pylab as pl
import random
import ECclass as ec
from scipy.integrate import odeint
  

def main():
    # number of particles 20 works well for separation. 30 Takes a long time
    N=20
    qq = 0.3
    beta = .6
    # y periodicity
    yd = 2.0*pl.pi
    num_cell = 2.0
    xd = num_cell * 2.0*pl.pi
    x_periodic = True
    # the order we want to go to as far as calculating periodic forces
    order = 3

    #N=2
    #qq = 1.3
    #beta = .6
    ## y periodicity
    #d = 2.0*pl.pi
    #num_cell =3.0
    

    As = pl.zeros(N) + 1.0
    #As = pl.zeros(N/2)+1.0
    #As = pl.append(As,pl.zeros(N/2)+.4)

    x0 = pl.zeros([4*N])

    some_colors = ['b','c','g','k','r','y','b','c','g','k','r','y','b','c','g','k','r','y']

    for i,j in enumerate(x0):
        if i in range(2*N,3*N):
            print(i)
            x0[i] = random.random()*xd
            continue
        if i in range(3*N,4*N):
            print(i)
            x0[i] = random.random()*yd
            continue

    t = pl.arange(0,40,.15)
    
    elec = ec.Flat2D(qq,As,beta,yd,num_cell,x_periodic,order)
    ## first just try to recreate what we did with 2 particles
    #x0[0] = 0.0
    #x0[1] = 0.0
    #x0[2] = 0.0
    #x0[3] = 0.0
    #x0[4] = 1.0
    #x0[5] = 1.0
    #x0[6] = 1.5
    #x0[7] = 2.5

    sol = odeint(elec.f,x0,t)

    for i in range(len(t)):
        fig = pl.figure()
        ax = fig.add_subplot(111)
        if x_periodic:
            ax.set_xlim([-0,xd])
        ax.set_ylim([0,yd])
        print(sol[i,:N])
        if x_periodic:
            ax.scatter(sol[i,2*N:(2*N+N/2)]%xd,sol[i,3*N:(3*N+N/2)]%yd,c='r')
            ax.scatter(sol[i,(2*N+N/2):3*N]%xd,sol[i,(3*N+N/2):4*N]%yd,c='b')
        else:
            ax.scatter(sol[i,2*N:(2*N+N/2)],sol[i,3*N:(3*N+N/2)]%yd,c='r')
            ax.scatter(sol[i,(2*N+N/2):3*N],sol[i,(3*N+N/2):4*N]%yd,c='b')

        #ax.scatter(pl.arange(0,diameter+pl.pi,pl.pi),pl.zeros(int(diameter/pl.pi+1.5)),s=5)
        #ax.scatter(pl.arange(0,diameter+pl.pi,pl.pi),pl.zeros(int(diameter/pl.pi+1.5))+pl.pi,s=5)
        print('len scatter array is: '+str(len(sol[i,2*N:3*N])))
        #ax.scatter(sol[i,5],sol[i,7])
        fig.savefig('%(number)04d.png'%{'number':i})
        pl.close(fig)
    if N == 2:
        print('final distanc')
        print(pl.sqrt((sol[-1,4]-sol[-1,5])**2+(sol[-1,6]-sol[-1,7])**2))

if __name__ == '__main__':
        main()
