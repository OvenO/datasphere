import pylab as pl
from scipy.integrate import odeint

def f(x,t):

    d = pl.sqrt((x[0]-x[4])**2+(x[2]-x[6])**2)
    x1dot = -(x[0]-x[4])/d**3
    x3dot = -(x[2]-x[6])/d**3
    x5dot = -(x[4]-x[0])/d**3
    x7dot = -(x[6]-x[2])/d**3
    x0dot = x[1]
    x2dot = x[3]
    x4dot = x[5]
    x6dot = x[7]

    return [x0dot,x1dot,x2dot,x3dot,x4dot,x5dot,x6dot,x7dot]

def main():
    t = pl.arange(0,10,.05)
    x0 = [0.0,0.0,0.0,0.0,.2,0.0,0.0,1.0]

    sol = odeint(f,x0,t)

    for i in range(len(t)):
        fig = pl.figure()
        ax = fig.add_subplot(111)
        ax.set_xlim([-.4,.4])
        ax.set_ylim([-.4,.4])
        ax.scatter(sol[i,0],sol[i,2])
        ax.scatter(sol[i,4],sol[i,6])
        fig.savefig('%(number)04d.png'%{'number':i})
        pl.close(fig)

if __name__ == '__main__':
        main()
