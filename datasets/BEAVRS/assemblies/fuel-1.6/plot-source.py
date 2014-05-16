import matplotlib.pyplot as plt
import numpy as np
import h5py

import sklearn.mixture.dpgmm as dpgmm

batches = [5,10,15,20,25]

for batch in batches:

  file = h5py.File('statepoint.' + str(batch) + '.h5')
  source = file['source_bank'][...]
  num_particles = len(source)

  sites = np.zeros((num_particles, 3))

  for particle in range(num_particles):
    sites[particle, :] = source[particle][1]

  fig = plt.figure()
  plt.scatter(sites[:,0], sites[:,1], linewidths=1)
  plt.title('Source Sites - 1.6% Enr. Assembly')
  plt.xlabel('x-coordinate')
  plt.ylabel('y-coordinate')
  plt.savefig('sources-batch-' + str(batch) + '.png')

  model = dpgmm.DPGMM(n_components=250, alpha=1000, n_iter=1000)
  model.fit(sites)
  print '# components = %d' % model.n_components
  print model.converged_
  print model.means_

  file.close()