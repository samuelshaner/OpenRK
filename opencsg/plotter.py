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


def plot_cells(geometry, plane='xy', midplane=0., gridsize=250):

  global subdirectory

  # Make directory if it does not exist
  if not os.path.exists(subdirectory):
    os.makedirs(subdirectory)

  # Error checking
  if not isinstance(geometry, Geometry):
    msg = 'Unable to plot the cells since input was not ' \
          'a Geometry class object'
    raise ValueError(msg)

  if not is_integer(gridsize):
    msg = 'Unable to plot the cells since the gridsize {0} is' \
          'is not an integer'.format(gridsize)
    raise ValueError(msg)

  if gridsize <= 0:
    msg = 'Unable to plot the cells with a negative ' \
          'gridsize {0}'.format(gridsize)
    raise ValueError(msg)

  if plane not in ['xy', 'xz', 'yz']:
    msg = 'Unable to plot the cells with an invalid ' \
        'plane {0}. Plane options xy, xz, yz'.format(plane)
    raise ValueError(msg)

  if not is_float(midplane):
    msg = 'Unable to plot the cells since the midplane {0} is' \
          'is not a float'.format(gridsize)
    raise ValueError(msg)

  print('Plotting the Cells...')

  # Initialize a NumPy array for the surface colors
  surface = numpy.zeros((gridsize, gridsize))

  # Retrieve the bounding box for the Geometry
  bounds = geometry.getBounds()
  if plane == 'xy':
    xcoords = np.linspace(bounds[0], bounds[1], gridsize)
    ycoords = np.linspace(bounds[2], bounds[3], gridsize)
    if midplane < bounds[4] or midplane > bounds[5]:
      msg = 'Unable to plot midplane at z={0} as it must lie ' \
          'between the z bounds[{1},{2}]'.format(midplane, \
                                                   bounds[4], bounds[5])
      raise ValueError(msg)
    del bounds[4:]
  elif plane == 'xz':
    xcoords = np.linspace(bounds[0], bounds[1], gridsize)
    zcoords = np.linspace(bounds[4], bounds[5], gridsize)
    if midplane < bounds[2] or midplane > bounds[3]:
      msg = 'Unable to plot midplane at y={0} as it must lie ' \
          'between the y bounds[{1},{2}]'.format(midplane, \
                                                   bounds[2], bounds[3])
      raise ValueError(msg)
    del bounds[2:4]
  else:
    ycoords = np.linspace(bounds[2], bounds[3], gridsize)
    zcoords = np.linspace(bounds[4], bounds[5], gridsize)
    if midplane < bounds[0] or midplane > bounds[1]:
      msg = 'Unable to plot midplane at x={0} as it must lie ' \
          'between the x bounds[{1},{2}]'.format(midplane, \
                                                   bounds[0], bounds[1])
      raise ValueError(msg)
    del bounds[:2]
  
  # Find the flat source region IDs for each grid point
  for i in range(gridsize):
    for j in range(gridsize):

      if plane == 'xy':
        cell = geometry.findCell(x=xcoords[i], y=ycoords[j], z=midplane)
      elif plane == 'xz':
        cell = geometry.findCell(x=xcoords[i], y=midplane, z=zcoords[j])
      else:
        cell = geometry.findCell(x=midplane, y=ycoords[i], z=zcoords[j])  

      surface[j][i] = color_map[cell._id % num_colors]

  # Plot a 2D color map of the flat source regions
  fig = plt.figure()
  surface = np.flipud(surface)
  plt.imshow(surface, extent=bounds)
  plt.title('Cells ' + plane)
  filename = subdirectory + 'cells-' + plane + '.png'
  fig.savefig(filename, bbox_inches='tight')
  plt.close(fig)



def plot_materials(geometry, plane='xy', midplane=0., gridsize=250):

  global subdirectory

  # Make directory if it does not exist
  if not os.path.exists(subdirectory):
    os.makedirs(subdirectory)

  # Error checking
  if not isinstance(geometry, Geometry):
    msg = 'Unable to plot the materials since input was not ' \
          'a Geometry class object'
    raise ValueError(msg)

  if not is_integer(gridsize):
    msg = 'Unable to plot the materials since the gridsize {0} is' \
          'is not an integer'.format(gridsize)
    raise ValueError(msg)

  if gridsize <= 0:
    msg = 'Unable to plot the materials with a negative ' \
          'gridsize {0}'.format(gridsize)
    raise ValueError(msg)

  if plane not in ['xy', 'xz', 'yz']:
    msg = 'Unable to plot the materials with an invalid ' \
        'plane {0}. Plane options xy, xz, yz'.format(plane)
    raise ValueError(msg)

  if not is_float(midplane):
    msg = 'Unable to plot the materials since the midplane {0} is' \
          'is not a float'.format(gridsize)
    raise ValueError(msg)

  print('Plotting the Materials...')

  # Initialize a NumPy array for the surface colors
  surface = numpy.zeros((gridsize, gridsize))

  # Retrieve the bounding box for the Geometry
  bounds = geometry.getBounds()
  if plane == 'xy':
    xcoords = np.linspace(bounds[0], bounds[1], gridsize)
    ycoords = np.linspace(bounds[2], bounds[3], gridsize)
    if midplane < bounds[4] or midplane > bounds[5]:
      msg = 'Unable to plot midplane at z={0} as it must lie ' \
          'between the z bounds[{1},{2}]'.format(midplane, \
                                                   bounds[4], bounds[5])
      raise ValueError(msg)
    del bounds[4:]
  elif plane == 'xz':
    xcoords = np.linspace(bounds[0], bounds[1], gridsize)
    zcoords = np.linspace(bounds[4], bounds[5], gridsize)
    if midplane < bounds[2] or midplane > bounds[3]:
      msg = 'Unable to plot midplane at y={0} as it must lie ' \
          'between the y bounds[{1},{2}]'.format(midplane, \
                                                   bounds[2], bounds[3])
      raise ValueError(msg)
    del bounds[2:4]
  else:
    ycoords = np.linspace(bounds[2], bounds[3], gridsize)
    zcoords = np.linspace(bounds[4], bounds[5], gridsize)
    if midplane < bounds[0] or midplane > bounds[1]:
      msg = 'Unable to plot midplane at x={0} as it must lie ' \
          'between the x bounds[{1},{2}]'.format(midplane, \
                                                   bounds[0], bounds[1])
      raise ValueError(msg)
    del bounds[:2]

  # Find the flat source region IDs for each grid point
  for i in range(gridsize):
    for j in range(gridsize):

      if plane == 'xy':
        cell = geometry.findCell(x=xcoords[i], y=ycoords[j], z=midplane)
      elif plane == 'xz':
        cell = geometry.findCell(x=xcoords[i], y=midplane, z=zcoords[j])
      else:
        cell = geometry.findCell(x=midplane, y=ycoords[i], z=zcoords[j])  

      surface[j][i] = color_map[cell._fill._id % num_colors]

  # Plot a 2D color map of the flat source regions
  fig = plt.figure()
  surface = np.flipud(surface)
  plt.imshow(surface, extent=bounds)
  plt.title('Materials ' + plane)
  filename = subdirectory + 'materials-' + plane + '.png'
  fig.savefig(filename, bbox_inches='tight')
  plt.close(fig)


def plot_regions(geometry, plane='xy', midplane=0., gridsize=250):

  global subdirectory

  # Make directory if it does not exist
  if not os.path.exists(subdirectory):
    os.makedirs(subdirectory)

  # Error checking
  if not isinstance(geometry, Geometry):
    msg = 'Unable to plot the regions since input was not ' \
          'a Geometry class object'
    raise ValueError(msg)

  if not is_integer(gridsize):
    msg = 'Unable to plot the regions since the gridsize {0} is' \
          'is not an integer'.format(gridsize)
    raise ValueError(msg)

  if gridsize <= 0:
    msg = 'Unable to plot the regions with a negative ' \
          'gridsize {0}'.format(gridsize)
    raise ValueError(msg)

  if plane not in ['xy', 'xz', 'yz']:
    msg = 'Unable to plot the regions with an invalid ' \
        'plane {0}. Plane options xy, xz, yz'.format(plane)
    raise ValueError(msg)

  if not is_float(midplane):
    msg = 'Unable to plot the regions since the midplane {0} is' \
          'is not a float'.format(gridsize)
    raise ValueError(msg)

  print('Plotting the Regions...')

  # Initialize a NumPy array for the surface colors
  surface = numpy.zeros((gridsize, gridsize))

  # Retrieve the bounding box for the Geometry
  bounds = geometry.getBounds()
  if plane == 'xy':
    xcoords = np.linspace(bounds[0], bounds[1], gridsize)
    ycoords = np.linspace(bounds[2], bounds[3], gridsize)
    if midplane < bounds[4] or midplane > bounds[5]:
      msg = 'Unable to plot midplane at z={0} as it must lie ' \
          'between the z bounds[{1},{2}]'.format(midplane, \
                                                   bounds[4], bounds[5])
      raise ValueError(msg)
    del bounds[4:]
  elif plane == 'xz':
    xcoords = np.linspace(bounds[0], bounds[1], gridsize)
    zcoords = np.linspace(bounds[4], bounds[5], gridsize)
    if midplane < bounds[2] or midplane > bounds[3]:
      msg = 'Unable to plot midplane at y={0} as it must lie ' \
          'between the y bounds[{1},{2}]'.format(midplane, \
                                                   bounds[2], bounds[3])
      raise ValueError(msg)
    del bounds[2:4]
  else:
    ycoords = np.linspace(bounds[2], bounds[3], gridsize)
    zcoords = np.linspace(bounds[4], bounds[5], gridsize)
    if midplane < bounds[0] or midplane > bounds[1]:
      msg = 'Unable to plot midplane at x={0} as it must lie ' \
          'between the x bounds[{1},{2}]'.format(midplane, \
                                                   bounds[0], bounds[1])
      raise ValueError(msg)
    del bounds[:2]

  # Initialize the offsets used for computing region IDs
  geometry.initializeCellOffsets()

  # Find the flat source region IDs for each grid point
  for i in range(gridsize):
    for j in range(gridsize):

      if plane == 'xy':
        region_id = geometry.getRegionId(x=xcoords[i], y=ycoords[j], z=midplane)
      elif plane == 'xz':
        region_id = geometry.getRegionId(x=xcoords[i], y=midplane, z=zcoords[j])
      else:
        region_id = geometry.getRegionId(x=midplane, y=ycoords[i], z=zcoords[j])

      surface[j][i] = color_map[region_id % num_colors]

  # Plot a 2D color map of the flat source regions
  fig = plt.figure()
  surface = np.flipud(surface)
  plt.imshow(surface, extent=bounds)
  plt.title('Regions ' + plane)
  filename = subdirectory + 'regions-' + plane + '.png'
  fig.savefig(filename, bbox_inches='tight')
  plt.close(fig)
