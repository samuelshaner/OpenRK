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

import numpy as np
import numpy.random
import os

## A static variable for the output directory in which to save plots
subdirectory = "plots/"


def plot_cells(geometry, plane='xy', offset=0., gridsize=250):

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

  if not is_float(offset):
    msg = 'Unable to plot the cells since the offset {0} is' \
          'is not a float'.format(gridsize)
    raise ValueError(msg)

  print('Plotting the Cells...')

  # Get the number of Cells filled with Materials
  cells = geometry.getAllMaterialCells()
  num_cells = len(cells)

  # Create array of equally spaced randomized floats as a color map for plots
  # Seed the NumPy random number generator to ensure reproducible color maps
  numpy.random.seed(1)
  color_map = np.linspace(0., 1., num_cells, endpoint=False)
  numpy.random.shuffle(color_map)

  # Initialize a NumPy array for the surface colors
  surface = numpy.zeros((gridsize, gridsize))

  # Retrieve the pixel coordinates
  coords = get_pixel_coords(geometry, plane, offset, gridsize)

  # Find the flat source region IDs for each grid point
  for i in range(gridsize):
    for j in range(gridsize):

      if plane == 'xy':
        cell = geometry.findCell(x=coords['x'][i], y=coords['y'][j], z=offset)
      elif plane == 'xz':
        cell = geometry.findCell(x=coords['x'][i], y=offset, z=coords['z'][j])
      else:
        cell = geometry.findCell(x=offset, y=coords['y'][i], z=coords['z'][j])  

      surface[j][i] = color_map[cell._id % num_cells]

  # Plot a 2D color map of the flat source regions
  fig = plt.figure()
  surface = np.flipud(surface)
  plt.imshow(surface, extent=coords['bounds'], interpolation="nearest")
  plt.title('Cells ' + plane)
  filename = subdirectory + 'cells-' + plane + '.png'
  fig.savefig(filename, bbox_inches='tight')
  plt.close(fig)



def plot_materials(geometry, plane='xy', offset=0., gridsize=250):

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

  if not is_float(offset):
    msg = 'Unable to plot the materials since the offset {0} is' \
          'is not a float'.format(gridsize)
    raise ValueError(msg)

  print('Plotting the Materials...')

  # Get the number of Cells filled with Materials
  materials = geometry.getAllMaterials()
  num_materials = len(materials)

  # Create array of equally spaced randomized floats as a color map for plots
  # Seed the NumPy random number generator to ensure reproducible color maps
  numpy.random.seed(1)
  color_map = np.linspace(0., 1., num_materials, endpoint=False)
  numpy.random.shuffle(color_map)

  # Initialize a NumPy array for the surface colors
  surface = numpy.zeros((gridsize, gridsize))

  # Retrieve the pixel coordinates
  coords = get_pixel_coords(geometry, plane, offset, gridsize)

  # Find the flat source region IDs for each grid point
  for i in range(gridsize):
    for j in range(gridsize):

      if plane == 'xy':
        cell = geometry.findCell(x=coords['x'][i], y=coords['y'][j], z=offset)
      elif plane == 'xz':
        cell = geometry.findCell(x=coords['x'][i], y=offset, z=coords['z'][j])
      else:
        cell = geometry.findCell(x=offset, y=coords['y'][i], z=coords['z'][j])  

      surface[j][i] = color_map[cell._fill._id % num_materials]

  # Plot a 2D color map of the flat source regions
  fig = plt.figure()
  surface = np.flipud(surface)
  plt.imshow(surface, extent=coords['bounds'], interpolation="nearest")
  plt.title('Materials ' + plane)
  filename = subdirectory + 'materials-' + plane + '.png'
  fig.savefig(filename, bbox_inches='tight')
  plt.close(fig)


def plot_regions(geometry, plane='xy', offset=0., gridsize=250):

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

  if not is_float(offset):
    msg = 'Unable to plot the regions since the offset {0} is' \
          'is not a float'.format(gridsize)
    raise ValueError(msg)

  print('Plotting the Regions...')

  # Initialize the offsets used for computing region IDs
  geometry.initializeCellOffsets()
  num_regions = geometry._num_regions

  # Create array of equally spaced randomized floats as a color map for plots
  # Seed the NumPy random number generator to ensure reproducible color maps
  numpy.random.seed(1)
  color_map = np.linspace(0., 1., num_regions, endpoint=False)
  numpy.random.shuffle(color_map)

  # Initialize a NumPy array for the surface colors
  surface = numpy.zeros((gridsize, gridsize))

  # Retrieve the pixel coordinates
  coords = get_pixel_coords(geometry, plane, offset, gridsize)

  # Find the flat source region IDs for each grid point
  for i in range(gridsize):
    for j in range(gridsize):

      if plane == 'xy':
        region_id = geometry.getRegionId(x=coords['x'][i], y=coords['y'][j], z=offset)
      elif plane == 'xz':
        region_id = geometry.getRegionId(x=coords['x'][i], y=offset, z=coords['z'][j])
      else:
        region_id = geometry.getRegionId(x=offset, y=coords['y'][i], z=coords['z'][j])

      surface[j][i] = color_map[region_id % num_regions]

  # Plot a 2D color map of the flat source regions
  fig = plt.figure()
  surface = np.flipud(surface)
  plt.imshow(surface, extent=coords['bounds'], interpolation="nearest")
  plt.title('Regions ' + plane)
  filename = subdirectory + 'regions-' + plane + '.png'
  fig.savefig(filename, bbox_inches='tight')
  plt.close(fig)


def plot_neighbor_cells(geometry, plane='xy', offset=0.,
                        gridsize=250, unique=False):

  global subdirectory

  # Make directory if it does not exist
  if not os.path.exists(subdirectory):
    os.makedirs(subdirectory)

  # Error checking
  if not isinstance(geometry, Geometry):
    msg = 'Unable to plot the neighbor cells since input was not ' \
          'a Geometry class object'
    raise ValueError(msg)

  if not is_integer(gridsize):
    msg = 'Unable to plot the neighbor cells since the gridsize {0} is' \
          'is not an integer'.format(gridsize)
    raise ValueError(msg)

  if gridsize <= 0:
    msg = 'Unable to plot the neighbor cells with a negative ' \
          'gridsize {0}'.format(gridsize)
    raise ValueError(msg)

  if plane not in ['xy', 'xz', 'yz']:
    msg = 'Unable to plot the neighbor cells with an invalid ' \
        'plane {0}. Plane options xy, xz, yz'.format(plane)
    raise ValueError(msg)

  if not is_float(offset):
    msg = 'Unable to plot the neighbor cells since the offset {0} is' \
          'is not a float'.format(gridsize)
    raise ValueError(msg)

  print('Plotting the Neighbor Cells...')

  # Initialize the offsets used for computing region IDs
  geometry.initializeCellOffsets()

  # Build the neighbor Cells/Universes
  geometry.buildNeighbors()

  if unique:
    num_neighbors = geometry._num_neighbors
  else:
    num_neighbors = geometry._num_unique_neighbors

  # Create array of equally spaced randomized floats as a color map for plots
  # Seed the NumPy random number generator to ensure reproducible color maps
  numpy.random.seed(1)
  color_map = np.linspace(0., 1., num_neighbors, endpoint=False)
  numpy.random.shuffle(color_map)

  print 'neighbor color map: {0}'.format(color_map)

  # Initialize a NumPy array for the surface colors
  surface = numpy.zeros((gridsize, gridsize))

  # Retrieve the pixel coordinates
  surf = get_pixel_coords(geometry, plane, offset, gridsize)

  # Find the flat source region IDs for each grid point
  for i in range(gridsize):
    for j in range(gridsize):

      if plane == 'xy':
        region_id = geometry.getRegionId(x=surf['x'][i], y=surf['y'][j], z=offset)
      elif plane == 'xz':
        region_id = geometry.getRegionId(x=surf['x'][i], y=offset, z=surf['z'][j])
      else:
        region_id = geometry.getRegionId(x=offset, y=surf['y'][i], z=surf['z'][j])

      if unique:
        neighbors = geometry._regions_to_unique_neighbors[region_id]
        neighbor_id = geometry._unique_neighbor_ids[neighbors]
      else:
        neighbors = geometry._regions_to_neighbors[region_id]
        neighbor_id = geometry._neighbor_ids[neighbors]

      surface[j][i] = color_map[neighbor_id % num_neighbors]


  # Plot a 2D color map of the flat source regions
  fig = plt.figure()
  surface = np.flipud(surface)
  plt.imshow(surface, extent=surf['bounds'], interpolation="nearest")

  if unique:
    plt.title('Unique Neighbor Cells ' + plane)
    filename = subdirectory + 'unique-neighbor-cells-' + plane + '.png'
  else:
    plt.title('Neighbor Cells ' + plane)
    filename = subdirectory + 'neighbor-cells-' + plane + '.png'

  fig.savefig(filename, bbox_inches='tight')
  plt.close(fig)


def get_pixel_coords(geometry, plane, offset, gridsize):

  # initialize variables to be returned
  bounds = geometry.getBounds()
  xcoords = None
  ycoords = None
  zcoords = None
  coords = dict()

  if plane == 'xy':
    xcoords = np.linspace(bounds[0], bounds[1], gridsize)
    ycoords = np.linspace(bounds[2], bounds[3], gridsize)
    if offset < bounds[4] or offset > bounds[5]:
      msg = 'Unable to plot offset at z={0} as it must lie ' \
          'between the z bounds[{1},{2}]'.format(offset, \
                                                   bounds[4], bounds[5])
      raise ValueError(msg)
    del bounds[4:]
  elif plane == 'xz':
    xcoords = np.linspace(bounds[0], bounds[1], gridsize)
    zcoords = np.linspace(bounds[4], bounds[5], gridsize)
    if offset < bounds[2] or offset > bounds[3]:
      msg = 'Unable to plot offset at y={0} as it must lie ' \
          'between the y bounds[{1},{2}]'.format(offset, \
                                                   bounds[2], bounds[3])
      raise ValueError(msg)
    del bounds[2:4]
  else:
    ycoords = np.linspace(bounds[2], bounds[3], gridsize)
    zcoords = np.linspace(bounds[4], bounds[5], gridsize)
    if offset < bounds[0] or offset > bounds[1]:
      msg = 'Unable to plot offset at x={0} as it must lie ' \
          'between the x bounds[{1},{2}]'.format(offset, \
                                                   bounds[0], bounds[1])
      raise ValueError(msg)
    del bounds[:2]
  
  # add attributes to coords dictionary
  coords['x'] = xcoords
  coords['y'] = ycoords
  coords['z'] = zcoords
  coords['bounds'] = bounds

  return coords
