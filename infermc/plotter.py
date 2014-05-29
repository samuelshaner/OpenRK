from infermc.process import *
from opencsg.checkvalue import *
from statepoint import StatePoint
import os
import matplotlib.pyplot as plt
import numpy as np
import numpy.random


## A static variable for the output directory in which to save plots
directory = "plots/"

## The number of colors to use when creating a random color map for plots
num_colors = 50

## An array of random floats that represents a random color map for plots
color_map = np.random.random_sample((num_colors,))


def plot_fluxes(geometry, statepoint, energies=[0], gridsize=250):

  global directory

  # Make directory if it does not exist
  if not os.path.exists(directory):
    os.makedirs(directory)

  if not isinstance(geometry, opencsg.Geometry):
    msg = 'Unable to plot the scalar flux since {0} is not an ' \
          'OpenCSG Geometry class object'.format(geometry)
    raise ValueError(msg)

  if not isinstance(statepoint, StatePoint):
    msg = 'Unable to plot the scalar flux since {0} is not an ' \
          'OpenMC StatePoint class object'.format(statepoint)
    raise ValueError(msg)

  if isinstance(energies, (list, tuple, np.array)):

    for energy in energies:

      if not is_integer(energy):
        msg ='Unable to plot the scalar flux since the energies list ' \
             'contains {0} which is not an integer value'.format(energy)
        raise ValueError(msg)

      elif energy < 0:
        msg = 'Unable to plot the scalar flux since the energies list ' \
              'contains {0} which is a negative value'.format(energies)
        raise ValueError(msg)

  elif is_integer(energies):

    if energies < 0:
      msg = 'Unable to plot the scalar flux since the energies argument ' \
            'is {0} which is a negative value'.format(energies)
      raise ValueError(msg)

    energies = [energies]

  else:
    msg = 'Unable to plot the scalar flux since the energies {0} is not ' \
          'an integer or Python list/tuple or NumPy array'.format(energies)
    raise ValueError(msg)

  if not is_integer(gridsize):
    msg = 'Unable to plot the scalar flux since the gridsize {0} is not ' \
          'an integer'.format(gridsize)
    raise ValueError(msg)

  if gridsize <= 0:
    msg = 'Unable to plot the scalar flux with a negative gridsize ' \
          '{0}'.format(gridsize)
    raise ValueError(msg)

  print('Plotting the scalar flux tallies...')


  # Initialize a numpy array for the groupwise scalar fluxes
  fluxes = numpy.zeros((len(energies), gridsize, gridsize))

  # Retrieve the bounding box for the geometry
  xmin = geometry.getMinX()
  xmax = geometry.getMaxX()
  ymin = geometry.getMinY()
  ymax = geometry.getMaxY()

  # Initialize numpy arrays for the grid points
  xcoords = np.linspace(xmin, xmax, gridsize)
  ycoords = np.linspace(ymin, ymax, gridsize)

  tally_extractor = XSTallyExtractor(statepoint=statepoint, geometry=geometry)
  cells_to_tallies = tally_extractor.getCellsToTallies()

  num_regions = geometry._num_regions
  region_fluxes = np.zeros((num_regions, 2))

#  volumes = geometry._region_volumes

  for region in range(num_regions):
    region_fluxes[region, :] = \
      tally_extractor.getDistribcellTallyData(region, 'flux')

  # Normalize the flux to the volume
#  region_fluxes[:, 0] /= volumes
#  region_fluxes[:, 1] /= volumes

  for i in range(gridsize):
    for j in range(gridsize):

      # Find the flat source region IDs for each grid point
      x = xcoords[i]
      y = ycoords[j]

      region = geometry.getRegionId(x=x, y=y)
      flux = region_fluxes[region, :]

#      cell_id = path[-1]
#      tally_id = cells_to_tallies[cell_id]
#      filters = [('distribcell', path)]
#      flux = statepoint.get_value(tally_id, filters, 0)

      # Get the scalar flux for each energy
      for index, energy in enumerate(energies):
        fluxes[index, j, i] = flux[index]

  # Loop over all energy group and create a plot
  for index, energy in enumerate(energies):

    # Plot a 2D color map of the flat source regions
    fig = plt.figure()
    plt.imshow(np.flipud(fluxes[index, :, :]), extent=[xmin, xmax, ymin, ymax])
    plt.colorbar()
    plt.title('Scalar Flux in Energy {0}'.format(energy))
    filename = '{0}flux-energy-{1}.png'.format(directory, energy)
    fig.savefig(filename, bbox_inches='tight')