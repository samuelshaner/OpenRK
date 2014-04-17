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
    print('Unable to plot the scalar flux since % is not an OpenCSG Geometry '
          'class object' % str(geometry))

  if not isinstance(statepoint, StatePoint):
    print('Unable to plot the scalar flux since %s is not an OpenMC StatePoint'
          'class object' % str(statepoint))

  if isinstance(energies, (list, tuple, np.array)):

    for energy in energies:

      if not is_integer(energy):
        print('Unable to plot the scalar flux since the energies list '
              'contains %s which is not an integer value' % str(energy))

      elif energy < 0:
        print('Unable to plot the scalar flux since the energies list '
              'contains %d which is a negative value' % energies)

  elif is_integer(energies):

    if energies < 0:
      print('Unable to plot the scalar flux since the energies argument '
            'is %d which is a negative value' % energies)

    energies = [energies]

  else:
    print('Unable to plot the scalar flux since the energies %s is not '
          'an integer or Python list/tuple or NumPy array' % str(energies))

  if not is_integer(gridsize):
    print('Unable to plot the scalar flux since the gridsize %s is not '
          'an integer' % str(gridsize))

  if gridsize <= 0:
    print('Unable to plot the scalar flux with a negative gridsize '
          '%d' % gridsize)

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

  tally_extractor = TallyExtractor(statepoint=statepoint, geometry=geometry)
  cells_to_tallies = tally_extractor.getCellsToTallies()

  for i in range(gridsize):
    for j in range(gridsize):

      # Find the flat source region IDs for each grid point
      x = xcoords[i]
      y = ycoords[j]

      coords = geometry.findCoords(x=x, y=y)
      path = get_path(coords)

      cell_id = path[-1]
      tally_id = cells_to_tallies[cell_id]
      filters = [('distribcell', path)]
      flux = statepoint.get_value(tally_id, filters, 0)

      # Get the scalar flux for each energy
      for index, energy in enumerate(energies):
        fluxes[index, j, i] = flux[index]

  # Loop over all energy group and create a plot
  for index, energy in enumerate(energies):

    # Plot a 2D color map of the flat source regions
    fig = plt.figure()
    plt.imshow(np.flipud(fluxes[index, :, :]), extent=[xmin, xmax, ymin, ymax])
    plt.colorbar()
    plt.title('Scalar Flux in Energy ' + str(energy))
    filename = directory + 'flux-energy-' + str(energy) + '.png'
    fig.savefig(filename, bbox_inches='tight')