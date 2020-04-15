import numpy as np 
import scipy.signal 

def msd_scipy(trajectory: np.ndarray, 
              max_delta: int = None) -> np.ndarray:
    """ This should be the fastest implementation of MSD
    We use the same tricks as in msd_numpy, but rely on SciPy
    to provide the convolution implementation. SciPy
    automatically decides whether to use FFT or not.
    :param trajectory: a numpy array containing the sequence
    of positions ordered in time.
    :param max_delta: the maximum time delta for which to 
    compute the MSD.
    :return: the mean squared displacement sequence
    """
    if max_delta is None:
        max_delta = len(trajectory)
        
    assert max_delta <= len(trajectory), \
        'max_delta must not exceed the trajectory length'
    
    squares = trajectory * trajectory
    
    # Computing:
    # bwd_sums = np.zeros(max_delta)
    # for delta in range(1, max_delta):
    #    bwd_sums[delta] = (trajectory[:-delta] ** 2).sum()
    bwd_sums = np.flipud(np.cumsum(squares))[:max_delta]
    
    # Computing:
    # fwd_sums = np.zeros(max_delta)
    # for delta in range(1, max_delta):
    #    fwd_sums[delta] = (trajectory[delta:] ** 2).sum()
    fwd_sums = np.flipud(np.cumsum(np.flipud(squares)))[:max_delta]
    
    # Computing
    # corrs = np.zeros(max_delta)
    # for delta in range(1, max_delta):
    #    corrs[delta] = (trajectory[delta:] * trajectory[:-delta]).sum()
    from scipy.signal import convolve
    corrs = convolve(trajectory, trajectory[::-1], mode='full')
    corrs = corrs[len(trajectory) - 1:][:max_delta]
    
    totals = bwd_sums + fwd_sums - 2 * corrs
    numels = np.arange(len(trajectory), 
                       len(trajectory) - max_delta, -1)
    
    return totals / numels


def msd(x):
  y = x * x
  cs1 = np.flipud(np.cumsum(y))
  cs2 = np.flipud(np.cumsum(np.flipud(y)))
  prod = scipy.signal.convolve(x, x[::-1], mode='full')
  prod = prod[prod.size // 2:]
  ans = cs1 + cs2 - 2 * prod
  return ans / np.hstack([np.arange(x.size, 1, -1), 1])

def structure_factor(x,box_l):
    '''
    Input:
    x: array-like, x.shape=(N_particles,N_time,N_dim), all particle postions

    Output:
    sq: float structure factor
    '''
    print("box: ", box_l)

    # next steps
    # make vector out of qval, with allowed 2pi/box_vector elements
    # take real gel state and calculate on patchy02
    # fix other factors
    # calculate structure factor on patchy02 (more Ram)
    # where does the sin(x)/x come frome
    # test with actual e(-iq (delta r)) calculation


    #qval = np.arange(2,1000,2)*np.pi/box_l[0] 
    qval = np.zeros((198,2))

    qval[:,0] = np.arange(2,200)*np.pi/box_l[0] 
    qval[:,1] = np.arange(2,200)*np.pi/box_l[1] 
    
    N_particles,N_time,N_dim = x.shape

    m_dim = int((N_particles)*(N_particles-1)/2)

    #pdist = np.zeros((N_particles,N_particles,2))
    dist = np.zeros((m_dim,2))   
    times = list(range(N_time))
    s_factor = np.zeros((len(times), len(qval)))
    for ti,tau in enumerate(times): 
        x_tau = x[:,tau,:]
        for i in range(2):
            p = np.reshape(x_tau[:,i], (N_particles,1))
            p = p - p.transpose()
            p = p - box_l[i]*np.rint(p/box_l[i])
            dist[:,i] = p[np.triu_indices(N_particles, k=1)]
        qp = np.dot(qval,dist.transpose())
        #r = np.linalg.norm(pdist,axis=2)

        #r = r[np.triu_indices(N_particles, k = 1)]
        #r = np.reshape(r,(-1,1))
        #qr = np.dot(qval,r.transpose())
        sf = np.exp(-1j*qp) 
        s_factor[ti] = np.sum(sf,axis=1)

    sq=np.mean(s_factor,axis=0)/N_particles
    sq = np.reshape(sq, (-1,1))
    qval=np.reshape(np.linalg.norm(qval,axis=1),(-1,1))
    sq = np.concatenate((qval,sq), axis=1)

    return sq


def vanhove(x, box_l,tau):
  '''
  Input:
  x: array-like, x.shape=(N_particles,N_time,N_dim), all particle postions
  q: q-value, wave vector length

  Output:
  vanhove, array-like. vanhove.shape = (len(bins),2), vanhove function at time tau
  '''

  N_particles,N_time,N_dim = x.shape
  bins = np.linspace(0,5,501)
  print(bins)
  tau_vals = []
  for i in range(N_particles):
      a = x[i,tau:,:]
      b = x[i,:-tau,:]

      dist_ii_tau = a-b
      dist_ii_tau = dist_ii_tau - box_l*np.rint(dist_ii_tau/box_l)

      r=np.linalg.norm(dist_ii_tau,axis=1)
      tau_vals.extend(r)
  
  hist, bin_edges = np.histogram(tau_vals, bins=bins, density=True)
  hist = np.reshape(hist,(-1,1))
  bins = np.reshape(bins,(-1,1))
  vanhove = np.concatenate((bins[1:],hist), axis=1)

  return vanhove

def self_scattering_function(x,q,box_l,times):
    '''
    Input:
    x: array-like, x.shape=(N_particles,N_time,N_dim), all particle postions
    q: q-value, wave vector length

    Output:
    fc_qt: arraylike, fc_qt.shape(N_time), intermediate scattering function value at specific q-value
    (needs to still be normalized with S(q) afterwards)

    '''
    N_particles,N_time,N_dim = x.shape

    fs_q_time = np.zeros(len(times))

    pdist = np.zeros((N_particles,N_particles,2))
    x_0 = x[:,0,:]
    for ti,tau in enumerate(times):
      x_tau = x[:,tau,:]
      p = x_0 - x_tau
      pdist = p - box_l*np.rint(p/box_l)

      r = np.linalg.norm(pdist,axis=1)
      s_factor = np.sin(q*r)/(q*r)

      fs_q_time[ti] = np.mean(s_factor) 
 
    return fs_q_time 

# Broken, coz missing S(q) normalization 
def intermediate_scattering_function(x, N_particles, times, q,box_l):
    '''
    Input:
    x: array-like, x.shape=(N_particles,N_dim, N_time), all particle postions

    q: q-value, wave vector length

    Output:
    fc_qt: arraylike, fc_qt.shape(N_time), intermediate scattering function value at specific q-value
    (needs to still be normalized with S(q) afterwards)

    '''
    fc_q_time = np.zeros(len(times))

    pdist = np.zeros((N_particles,N_particles,2))
    xi_0 = x[:,0,:]
    for ti,tau in enumerate(times):
      xj_tau = x[:,tau,:]
      print(tau)
      for i in range(2):
          p0 = np.reshape(xi_0[:,i], (N_particles,1))
          ptau  = np.reshape(xj_tau[:,i], (N_particles,1))
          p = p0 - ptau.transpose()
          pdist[:,:,i] = p - box_l[i]*np.rint(p/box_l[i])

      r = np.linalg.norm(pdist,axis=2)
      r = r[np.triu_indices_from(r)]
      s_factor = np.sin(q*r)/(q*r)

      fc_q_time[ti] = np.mean(s_factor) 
 
    return fc_q_time 


#Note: now broken 
def alpha_2(pos):

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

import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-methods', nargs='+', choices=['msd', 'sq', 'fcqt','fsqt','vanhove'], type=str, required=True)
    args = parser.parse_args()


    N_particles=1500
    box_all = np.fromfile("Box.bin")
    box_l = box_all[3:5]

    pos = np.fromfile("positions.bin")
    pos = np.reshape(pos, (-1,3))
    pos = pos[:,:2]

    N_time = int(len(pos)/N_particles)

    x = np.zeros((N_particles,N_time,2))
    for i in range(N_particles):
      x[i] = pos[i::N_particles]

    if 'msd' in args.methods:
      msdx=np.zeros((N_particles,N_time))
      msdy=np.zeros((N_particles,N_time))

      # calc diffs 
      Dx = np.zeros((N_particles,N_time))
      Dy = np.zeros((N_particles,N_time))

      for pid in range(N_particles):

          x=pos[pid::N_particles,0]
          y=pos[pid::N_particles,1]

          dx = np.zeros(N_time)
          dy = np.zeros(N_time)

          pdx = x[1:] - x[:-1]
          dx[1:] = pdx - box_l[0]*np.rint(pdx/box_l[0])

          pdy = y[1:] - y[:-1]
          dy[1:] = pdy - box_l[1]*np.rint(pdy/box_l[1])

          Dx[pid] = dx
          Dy[pid] = dy 

      for pid in range(N_particles):
          msdx[pid] = msd_scipy(np.cumsum(Dx[pid]))
          msdy[pid] = msd_scipy(np.cumsum(Dy[pid]))

      r=np.mean(msdx, axis=0) + np.mean(msdy, axis=0)
      np.savetxt('msd.dat', r, newline='\n')

    if 'vanhove' in args.methods:
      times = [1,10,100,1000,10000,100000]
      for tau in times: 
        histogram = vanhove(x,box_l,tau)
        np.savetxt("vanhove_{}.dat".format(tau), histogram, newline='\n')

    if 'sq' in args.methods:
      sq = structure_factor(x,box_l)
      np.savetxt("structure_factor.dat", sq, newline='\n')

    if 'fsqt' in args.methods:
      q_vals = [0.1, 0.36, 1.81,3.26, 4.70, 6.15,7.60,9.04,10.49,16.28]
      times = list(range(N_time))
      wtimes = np.reshape(np.array(times), (-1,1))
      for q in q_vals:
          print(q)
          fcqt = self_scattering_function(x,q,box_l,times)
          fcqt = np.reshape(fcqt, (-1,1))
          result = np.concatenate((wtimes,fcqt), axis=1)
          np.savetxt("fst_{}.dat".format(q), result, newline='\n')

    if 'fcqt' in args.methods:
      q_vals = [0.1, 0.36, 1.81,3.26, 4.70, 6.15,7.60,9.04,10.49,16.28]
      times = list(range(N_time))[1::10]
      wtimes = np.reshape(np.array(times), (-1,1))
      for q in q_vals:
          print(q)
          fcqt = intermediate_scattering_function(x, N_particles, times, q,box_l)
          fcqt = np.reshape(fcqt, (-1,1))
          result = np.concatenate((wtimes,fcqt), axis=1)
          np.savetxt("fct_q_{}.dat".format(q), result, newline='\n')

