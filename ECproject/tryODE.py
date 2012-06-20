#!/opt/local/bin/python2.7
import pylab as pl
from scipy.integrate import odeint

# try the scipy ode equation solver for pendulum
def main():
    # initial conditions [theta,omega]
    y0 = [pl.pi/4,0.0]
    def pen(y,t):
        return [y[1],-pl.sin(y[0])]

    t = pl.arange(0.0,10.0,.05)

    sol = odeint(pen,y0,t)
    print(sol[:,0])
    #pl.plot(sol[:,0],sol[:,1])
    pl.plot(sol)
    pl.show()
if __name__ == "__main__":
    main()

