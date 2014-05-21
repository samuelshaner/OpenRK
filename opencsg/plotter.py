__author__ = 'Will Boyd'
__email__ = 'wboyd@mit.edu'


import matplotlib
from opencsg import *

# force headless backend, or set 'backend' to 'Agg'
# in your ~/.matplotlib/matplotlibrc
matplotlib.use('Agg')

import matplotlib.pyplot as plt

# Force non-interactive mode, or set 'interactive' to False
# in your ~/.matplotlib/matplotlibrc
plt.ioff()

import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np
import numpy.random
import os, sys


## A static variable for the output directory in which to save plots
subdirectory = "plots/"

## The number of colors to use when creating a random color map for plots
num_colors = 50

## An array of random floats that represents a random color map for plots
color_map = np.random.random_sample((num_colors,))


def plot_cells(geometry, gridsize=250):

  global subdirectory

  # Make directory if it does not exist
  if not os.path.exists(subdirectory):
    os.makedirs(subdirectory)

  # Error checking
  if not isinstance(geometry, Geometry):
    exit('Unable to plot the cells since input was ' + \
         'not a Geometry class object')

  if not is_integer(gridsize):
    exit('Unable to plot the cells since the gridsize is' + \
         'is not an integer' % str(gridsize))

  if gridsize <= 0:
    exit('Unable to plot the cells with a negative ' \
         'gridsize (%d)' % gridsize)

  print('Plotting the Cells...')

  # Initialize a NumPy array for the surface colors
  surface = numpy.zeros((gridsize, gridsize))

  # Retrieve the bounding box for the Geometry
  xmin = geometry.getMinX()
  xmax = geometry.getMaxX()
  ymin = geometry.getMinY()
  ymax = geometry.getMaxY()

  # Initialize numpy arrays for the grid points
  xcoords = np.linspace(xmin, xmax, gridsize)
  ycoords = np.linspace(ymin, ymax, gridsize)

  # Find the flat source region IDs for each grid point
  for i in range(gridsize):
    for j in range(gridsize):

      cell = geometry.findCell(x=xcoords[i], y=ycoords[j])
      surface[i][j] = color_map[cell.getId() % num_colors]

  # Plot a 2D color map of the flat source regions
  fig = plt.figure()
  plt.imshow(surface, extent=[xmin, xmax, ymin, ymax])
  plt.title('Cells')
  filename = subdirectory + 'cells.png'
  fig.savefig(filename, bbox_inches='tight')
  plt.close(fig)



def plot_materials(geometry, gridsize=250):

  global subdirectory

  # Make directory if it does not exist
  if not os.path.exists(subdirectory):
    os.makedirs(subdirectory)

  # Error checking
  if not isinstance(geometry, Geometry):
    exit('Unable to plot the materials since input was ' + \
         'not a Geometry class object')

  if not is_integer(gridsize):
    exit('Unable to plot the materials since the gridsize is' + \
         'is not an integer' % str(gridsize))

  if gridsize <= 0:
    exit('Unable to plot the materials with a negative ' \
         'gridsize (%d)' % gridsize)

  print('Plotting the Materials...')

  # Initialize a NumPy array for the surface colors
  surface = numpy.zeros((gridsize, gridsize))

  # Retrieve the bounding box for the Geometry
  xmin = geometry.getMinX()
  xmax = geometry.getMaxX()
  ymin = geometry.getMinY()
  ymax = geometry.getMaxY()

  # Initialize numpy arrays for the grid points
  xcoords = np.linspace(xmin, xmax, gridsize)
  ycoords = np.linspace(ymin, ymax, gridsize)

  # Find the flat source region IDs for each grid point
  for i in range(gridsize):
    for j in range(gridsize):

      cell = geometry.findCell(x=xcoords[i], y=ycoords[j])
      surface[i][j] = color_map[cell.getFill().getId() % num_colors]

  # Plot a 2D color map of the flat source regions
  fig = plt.figure()
  plt.imshow(surface, extent=[xmin, xmax, ymin, ymax])
  plt.title('Materials')
  filename = subdirectory + 'materials.png'
  fig.savefig(filename, bbox_inches='tight')
  plt.close(fig)


def plot_regions(geometry, gridsize=250):

  global subdirectory

  # Make directory if it does not exist
  if not os.path.exists(subdirectory):
    os.makedirs(subdirectory)

  # Error checking
  if not isinstance(geometry, Geometry):
    exit('Unable to plot the regions since input was ' + \
         'not a Geometry class object')

  if not is_integer(gridsize):
    exit('Unable to plot the regions since the gridsize is' + \
         'is not an integer' % str(gridsize))

  if gridsize <= 0:
    exit('Unable to plot the regions with a negative ' \
         'gridsize (%d)' % gridsize)

  print('Plotting the Regions...')

  # Initialize a NumPy array for the surface colors
  surface = numpy.zeros((gridsize, gridsize))

  # Retrieve the bounding box for the Geometry
  xmin = geometry.getMinX()
  xmax = geometry.getMaxX()
  ymin = geometry.getMinY()
  ymax = geometry.getMaxY()

  # Initialize numpy arrays for the grid points
  xcoords = np.linspace(xmin, xmax, gridsize)
  ycoords = np.linspace(ymin, ymax, gridsize)

  # Initialize the offsets used for computing region IDs
  geometry.initializeCellOffsets()

  # Find the flat source region IDs for each grid point
  for i in range(gridsize):
    for j in range(gridsize):

      region_id = geometry.getRegionId(x=xcoords[i], y=ycoords[j])
      surface[i][j] = color_map[region_id % num_colors]

  # Plot a 2D color map of the flat source regions
  fig = plt.figure()
  plt.imshow(surface, extent=[xmin, xmax, ymin, ymax])
  plt.title('Regions')
  filename = subdirectory + 'regions.png'
  fig.savefig(filename, bbox_inches='tight')
  plt.close(fig)
