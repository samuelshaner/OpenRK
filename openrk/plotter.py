__author__ = 'Sam Shaner'
__email__ = 'shaner@mit.edu'


import matplotlib
import openrk
from mesh import *

# force headless backend, or set 'backend' to 'Agg'
# in your ~/.matplotlib/matplotlibrc
matplotlib.use('Agg')

import matplotlib.pyplot as plt

# Force non-interactive mode, or set 'interactive' to False
# in your ~/.matplotlib/matplotlibrc
plt.ioff()

import numpy as np
import numpy.random
import os

## A static variable for the output directory in which to save plots
SUBDIRECTORY = "plots/"

def plot_flux(mesh, name, energy_groups=[0], gridsize=250):

  global SUBDIRECTORY

  # Make directory if it does not exist
  if not os.path.exists(SUBDIRECTORY):
    os.makedirs(SUBDIRECTORY)

  # Error checking
  if not isinstance(mesh, Mesh):
    msg = 'Unable to plot the cells since input was not ' \
          'a Mesh class object'
    raise ValueError(msg)

  if not is_integer(gridsize):
    msg = 'Unable to plot the cells since the gridsize {0} ' \
          'is not an integer'.format(gridsize)
    raise ValueError(msg)

  if gridsize <= 0:
    msg = 'Unable to plot the cells with a negative ' \
          'gridsize {0}'.format(gridsize)
    raise ValueError(msg)

  if not isinstance(energy_groups, list):
    energy_groups = [energy_groups]

  print('Plotting the fluxes...')

  # Initialize a numpy array for the groupwise scalar fluxes
  fluxes = numpy.zeros((len(energy_groups), gridsize, gridsize))

  TINY_MOVE = 1.0e-8
  bounds = mesh.getBounds()

  # Retrieve the bounding box for the geometry
  xmin = bounds[0] + TINY_MOVE
  xmax = bounds[1] - TINY_MOVE
  ymin = bounds[2] + TINY_MOVE
  ymax = bounds[3] - TINY_MOVE
  
  # Initialize numpy arrays for the grid points
  xcoords = np.linspace(xmin, xmax, gridsize)
  ycoords = np.linspace(ymin, ymax, gridsize)

  for i in range(gridsize):
    for j in range(gridsize):

      # Find the flat source region IDs for each grid point
      x = xcoords[i]
      y = ycoords[j]

      cell = mesh.findCell(x, y)

      for index, group in enumerate(energy_groups):
        fluxes[index][j][i] = mesh.getFlux(name, cell, group)

  # Loop over all energy group and create a plot
  for index, group in enumerate(energy_groups):

    # Plot a 2D color map of the flat source regions
    fig = plt.figure()
    plt.imshow(np.flipud(fluxes[index,:,:]), extent=[xmin, xmax, ymin, ymax])
    #plt.colorbar()
    plt.title('Mesh Cell Scalar Flux in Group ' + str(group))
    filename = SUBDIRECTORY + 'mesh-flux-group-' + str(group) + '.png'
    fig.savefig(filename, bbox_inches='tight')


