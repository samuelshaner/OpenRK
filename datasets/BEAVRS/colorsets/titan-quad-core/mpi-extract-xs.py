import openmc
import openmc.statepoint
from datasets.energy_groups import group_structures
from infermc.process import XSTallyExtractor, MicroXSTallyExtractor
from infermc.multigroupxs import xs_types
import infermc.plotter as plotter
from mpi4py import MPI
from geometry import geometry
import numpy as np



# Get MPI metadata
comm = MPI.COMM_WORLD
num_procs = comm.Get_size()
rank = comm.Get_rank()

first_batch = 55
last_batch = 70
batch_interval = 5

groups = group_structures['CASMO']['2-group']


#######  FIND MAX AND MIN CROSS-SECTION VALUES FOR EACH NUCLIDE AND GROUP ######

openmc_geometry = openmc.get_openmc_geometry(geometry)
nuclides = openmc_geometry.getAllNuclides()
tally_nuclides = [nuclides['H-1'][0], nuclides['B-10'][0], nuclides['O-16'][0],
                  nuclides['U-235'][0], nuclides['U-238'][0], nuclides['Zr-90'][0],
                  nuclides['Fe-56'][0]]

num_nuclides = len(tally_nuclides)
num_xs_types = len(xs_types)

min_xs = 1e5 * np.ones((num_nuclides, num_xs_types, 2))
max_xs = -1. * np.ones((num_nuclides, num_xs_types, 2))

# Loop over all batches
curr_batch = first_batch + rank * batch_interval

while (curr_batch >= first_batch) and (curr_batch <= last_batch):

  print('proc {0} opening batch {1}'.format(rank, curr_batch))

  filename = '../quad-core/statepoint.{0:03d}.h5'.format(curr_batch)

  # Initialize a handle on the OpenMC statepoint file
  statepoint = openmc.statepoint.StatePoint(filename)

  # Build MicroXS objects from the Tally data
  extractor = MicroXSTallyExtractor(statepoint)
  extractor.extractAllMultiGroupXS(groups, 'distribcell')

  # Look at max and min xs for each nuclide, xs type, and group
  for i, nuclide in enumerate(tally_nuclides):
    for j, xs_type in enumerate(xs_types):
      for k, group in enumerate([1,2]):
        curr_min_xs = extractor.getMinXS(xs_type, nuclide, 'distribcell', group)
        curr_max_xs = extractor.getMaxXS(xs_type, nuclide, 'distribcell', group)
        min_xs[i,j,k] = min(min_xs[i,j,k], curr_min_xs)
        max_xs[i,j,k] = max(max_xs[i,j,k], curr_max_xs)

  openmc.reset_auto_ids()
  del extractor, statepoint

  # Update current batch
  curr_batch += num_procs * batch_interval

all_min_xs = np.zeros((num_nuclides, num_xs_types, 2))
all_max_xs = np.zeros((num_nuclides, num_xs_types, 2))

for i, nuclide in enumerate(tally_nuclides):
  for j, xs_type in enumerate(xs_types):
    comm.Allreduce([min_xs[i,j,:], MPI.DOUBLE],
                   [all_min_xs[i,j,:], MPI.DOUBLE], op=MPI.MIN)
    comm.Allreduce([max_xs[i,j,:], MPI.DOUBLE],
                   [all_max_xs[i,j,:], MPI.DOUBLE], op=MPI.MAX)



#######  FIND MAX AND MIN CROSS-SECTION VALUES FOR EACH NUCLIDE AND GROUP ######

# Loop over all batches
curr_batch = first_batch + rank * batch_interval

while (curr_batch >= first_batch) and (curr_batch <= last_batch):

  print('proc {0} opening batch {1}'.format(rank, curr_batch))

  filename = '../quad-core/statepoint.{0:03d}.h5'.format(curr_batch)

  # Initialize a handle on the OpenMC statepoint file
  statepoint = openmc.statepoint.StatePoint(filename)

  ## MICROS
  extractor = MicroXSTallyExtractor(statepoint)
  extractor.extractAllMultiGroupXS(groups, 'cell')
  extractor.extractAllMultiGroupXS(groups, 'distribcell')

  for i, nuclide in enumerate(tally_nuclides):
    for j, xs_type in enumerate(xs_types):
      plotter.scatter_micro_xs(extractor, xs_types=[xs_type], nuclides=[nuclide],
                               domain_types=['distribcell', 'cell'],
                               colors=['cell', 'cell'],
                               filename='batch-{0}'.format(curr_batch),
                               xlim=[0.9*all_min_xs[i,j,0], 1.1*all_max_xs[i,j,0]],
                               ylim=[0.9*all_min_xs[i,j,1], 1.1*all_max_xs[i,j,1]])


  openmc.reset_auto_ids()
  del extractor, statepoint

  # Update current batch
  curr_batch += num_procs * batch_interval
