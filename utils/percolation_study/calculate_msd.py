import numpy as np 

def msd(x):
  y = x * x
  cs1 = np.flipud(np.cumsum(y))
  cs2 = np.flipud(np.cumsum(np.flipud(y)))
  prod = np.correlate(x, x, mode='full')
  prod = prod[prod.size // 2:]
  ans = cs1 + cs2 - 2 * prod
  return ans / np.hstack([np.arange(x.size, 1, -1), 1])

def structure_factor(x,q,box_l):
    '''
    Input:
    x: array-like, x.shape=(N_particles,N_dim, N_time), all particle postions
    q: q-value, wave vector length

    Output:
    s_q: float structure factor at specific q-value
    '''

    N_particles, N_dim, N_time = x.shape

    #TODO better distances calculation! 
    s_factor_time = np.zeros(N_time)

    for tau in range(N_time):
      s_factor_ij = [] 
      for i in range(N_particles):
        for j in range(N_particles):
          if i<j:
            a = x[i,:,tau]
            b = x[j,:,tau]

            dist_ij = a-b
            dist_ij = dist_ij - box_l*np.rint(dist_ij/box_l)

            r=np.linalg.norm(dist_ij,axis=1)
            s_factor_ij.append(np.sin(q*r)/(q*r))

        s_factor_time[tau] = np.mean(np.array(s_factor_ij))

    s_q = np.mean(s_factor_time, axis=0)
    return s_q


def self_scattering_function(x,q):
    '''
    Input:
    x: array-like, x.shape=(N_particles,N_dim, N_time), all particle postions
    q: q-value, wave vector length

    Output:
    s_qt:  array-like s_qt.shape = (N_time), self-scattering function value at specific q-value
    '''
    N_particles, N_dim, N_time = x.shape
    s_factor_i_time = np.zeros((N_particles,N_time))

    for i in range(N_particles):
        for tau in range(N_time):
            a = x[i,:,tau:]
            b = x[i,:,:-tau]

            r=np.linalg.norm(a-b,axis=1)
            # s_factor.shape = (N_time-tau)
            s_factor = np.sin(q*r)/(q*r)
            s_factor_i_time[i,tau] = np.mean(s_factor)

    s_factor_time = np.mean(s_factor_i_time, axis=0)
    return s_qt


def intermediate_scattering_function(x,q,box_l):
    '''
    Input:
    x: array-like, x.shape=(N_particles,N_dim, N_time), all particle postions
    q: q-value, wave vector length

    Output:
    fc_qt: arraylike, fc_qt.shape(N_time), intermediate scattering function value at specific q-value
    (needs to still be normalized with S(q) afterwards)

    '''
    N_particles, N_dim, N_time = x.shape
    fc_q_time = np.zeros(N_time)

    for tau in range(N_time):
      s_factor=[]
      for i in range(N_particles):
        for j in range(N_particles):
          if i<j: 
            a = x[i,:,tau:]
            b = x[j,:,:-tau]

            dist_ij_tau = a-b
            dist_ij_tau = dist_ij_tau - box_l*np.rint(dist_ij_tau/box_l)

            r=np.linalg.norm(dist_ij_tau,axis=1)
            # s_factor.shape = (N_time-tau)
            s_factor.append(np.sin(q*r)/(q*r))

      fc_q_time[tau] = np.mean(np.array(s_factor))

    return fc_q_time


def alpha_2(x):

    '''
    Input:
    x: array-like, x.shape=(N_particles,N_dim, N_time), all particle postions

    Output:
    alpha_2: array-like, alpha_2.shape=(N_time), non-gaussian parameter alpha_2 as function of time
    '''
    dx = np.zeros((N_particles,N_dim,N_time))
    dx[1:] = x[:,:,1:] - x[:,:,:-1]
    dx = dx - box_l[0]*np.rint(dx/box_l[0])
    x_pbc = np.cumsum(dx, axis=2)

    r = np.linalg.norm(x_pbc, axis=1)
    r2 = np.power(r,2)
    alpha_2 = (3/5)*N_particles*(np.mean(r2*r2, axis=0)/np.power(np.mean(r2,axis=0),2))

    return alpha_2





if __name__ == '__main__':

    N_particles=1500
    box_all = np.fromfile("Box.bin")
    box_l = box_all[3:5]

    pos = np.fromfile("positions.bin")
    pos = np.reshape(pos, (-1,3))
    pos = pos[:,:2]

    N_time = len(pos)/N_particles

    pos = np.reshape(pos, (N_particles,N_time,2))
    msdx=np.zeros((N_particles,N_time))
    msdy=np.zeros((N_particles,N_time))

    # calc diffs 
    Dx = np.zeros((N_particles,N_time))
    Dy = np.zeros((N_Particles,N_time))

    for pid in range(N_particles):
        x=pos[pid,:,0]
        y=pos[pid,:,1]

        dx = np.zeros(N_time)
        dy = np.zeros(N_time)

        pdx = x[1:] - x[:-1]
        dx[1:] = pdx - box_l[0]*np.rint(pdx/box_l[0])

        pdy = y[1:] - y[:-1]
        dy[1:] = pdy - box_l[1]*np.rint(pdy/box_l[1])

        Dx[pid] = dx
        Dy[pid] = dy 

    # calc msd 
    for pid in range(N_particles):
        msdx[pid] = msd(np.cumsum(Dx[pid]))
        msdy[pid] = msd(np.cumsum(Dy[pid]))

    r=np.mean(msdx, axis=0) + np.mean(msdy, axis=0)
    np.savetxt('msd.dat', r, newline='\n')

    # van Hove self correlation 
    times=np.array([1,10,100,1000,10000,100000])
    bin_size=0.01
    r_stop=3
    r_start=0

    r_vals = np.linspace(r_start,r_stop,bin_size)

    van_hove = np.zeros((len(times),len(r_vals)))
    for pid in range(N_particles):
        dr1 = np.sqrt(Dx[pid] +Dy[pid])
        dr_cumsum = np.cumsum(dr1)
        for time_i, time in enumerate(times):
            dR = dr_cumsum[time]
            rval_j = int(dR // bin_size) 
            van_hove[time_i,rval_j] += 1  

    np.savetxt("vanHove.dat", van_hove, newline='\n')

