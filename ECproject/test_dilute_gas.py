import pylab as pl
import random
import ECclass as ec
from scipy.integrate import odeint
import os
import shutil
  

def main():
    # skip ploting some so it is faster to watch
    skip = 5
    
    # do we want phase space imaes
    make_phase = False
    # do we want real space images
    make_real = True
    # just a last bit of the phase space plot (single plot)
    make_lphase = False
    # number of particles
    N=20
    #beta = .6
    xd = 10.0
    yd =10.0

    delta_x = xd/4.0
    delta_y = yd/5.0
    
    t = pl.arange(0,20,.02)

    x0 = pl.zeros([4*N])

    some_colors = ['b','c','g','k','r','y','b','c','g','k','r','y','b','c','g','k','r','y']
    # random positions
    # first arange in a square latice and then displace by a 4th the separatoin
    # Lets just do this for the 20 particles
    # x position
    for i in range(N):
        x0[2*N+i] = delta_x/2.0 + delta_x*(i%4) + delta_x/2.0*(random.random()-.5)
    count = 0
    for i in range(N):
        if i%4==0 and i>0: count +=1
        x0[3*N+i] = delta_y/2.0 + delta_y*(count) + delta_y/2.0*(random.random()-.5)
    for i in range(N):
        # each particle starts with vel magnitude 1 in random direction
        random_angle = 2.0*pl.pi*random.random()
        # x_vel
        x0[i]= pl.cos(random_angle)
        # y_vel
        x0[i+N]= pl.sin(random_angle)

#    pl.scatter(x0[2*N:3*N],x0[3*N:4*N])
#    pl.show()
 
    
    elec = ec.Test(xd,yd)
    
    sol = odeint(elec.f,x0,t)
    
#    cur_file = open('cur_run_data'
#    # modulus the x positions so they are in the right spot
#    sol[:,2*N:3*N] = sol[:,2*N:3*N]%x_modNum
#    # modulus the x positions so they are in the right spot
#    sol[:,3*N:4*N] = sol[:,3*N:4*N]%y_modNum
#    
#    if sliced:
#        for a in range(len(sol[:,0])):
#            if(((a*dt)%(2.0*pl.pi))<dt):
#                # str(sol)[1,-1] gets rid of the '[' and ']' at the begining and end of the array
#                # when it is turned into a string.
#                cur_file.write(str(sol[a,:])[1:-1].replace('\n',''))
#                cur_file.write('\n')
#    else:
#        for a in range(len(sol[:,0])):
#            cur_file.write(str(sol[a,:])[1:-1].replace('\n',''))
#            cur_file.write('\n')
#
#    cur_file.close()


    os.mkdir('RunImages')
    if make_phase:
        os.mkdir('RunImages/PhaseSpace')
    if make_real:
        os.mkdir('RunImages/RealSpace')

    for i in range(len(t)):
        if i%skip!=0:
            continue
        print('x components: '+str(sol[i,2*N:3*N]))
        # scatter for different As solutions
        #ax.scatter(sol[i,2*N:(2*N+N/2)],sol[i,3*N:(3*N+N/2)]%d,c='r')
        #ax.scatter(sol[i,(2*N+N/2):3*N],sol[i,(3*N+N/2):4*N]%d,c='b')
        if make_real:
            r_fig = pl.figure()
            r_ax = r_fig.add_subplot(111)
            r_ax.set_xlim([0,xd])
            r_ax.set_ylim([0,yd])

            r_ax.scatter(sol[i,2*N:3*N]%xd,sol[i,3*N:4*N]%yd,c=['b','r'])
            r_fig.savefig('RunImages/RealSpace/%(number)04d.png'%{'number':i})
            pl.close(r_fig)

#        if make_phase:
#            p_fig = pl.figure()
#            p_ax = p_fig.add_subplot(111)
#            p_ax.set_xlim([0,d])
#            p_ax.set_ylim([-2.0,2.0])
#            p_ax.scatter(sol[i,N:2*N]%d,sol[i,:N],c='b')
#            p_fig.savefig('RunImages/PhaseSpace/%(number)04d.png'%{'number':i})
#            pl.close(p_fig)


#    if make_lphase:
#        lp_fig = pl.figure()
#        lp_ax = lp_fig.add_subplot(111)
#        #lp_ax.set_xlim([0,d])
#        #lp_ax.set_ylim([-2.0,2.0])
#        for a in range(N):
#            lp_ax.plot(sol[-len(sol)/5:,N+a]%d,sol[-len(sol)/5:,a],c='k')
#        lp_fig.savefig('RunImages/last_phase.png',dpi= 300)
#        pl.close(lp_fig)


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
