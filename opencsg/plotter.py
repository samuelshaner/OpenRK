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
SUBDIRECTORY = "plots/"


def plot_cells(geometry, plane='xy', offset=0., gridsize=250):

  global SUBDIRECTORY

  # Make directory if it does not exist
  if not os.path.exists(SUBDIRECTORY):
    os.makedirs(SUBDIRECTORY)

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
          'is not a float'.format(offset)
    raise ValueError(msg)

  print('Plotting the Cells...')

  # Get the number of Cells filled with Materials
  cells = geometry.getAllMaterialCells()
  num_cells = len(cells)

  # Create array of equally spaced randomized floats as a color map for plots
  # Seed the NumPy random number generator to ensure reproducible color maps
  numpy.random.seed(1)
  colors = np.linspace(0., 1., num_cells, endpoint=False)
  numpy.random.shuffle(colors)

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

      # If we did not find a Cell for this region, use a NaN "bad" number color
      if cell is None:
        surface[j][i] = np.nan
      else:
        surface[j][i] = colors[cell._id % num_cells]

  # Make Matplotlib color "bad" numbers (ie, NaN, INF) with transparent pixels
  cmap = plt.get_cmap('spectral')
  cmap.set_bad(alpha=0.0)

  # Plot a 2D color map of the flat source regions
  fig = plt.figure()
  surface = np.flipud(surface)
  plt.imshow(surface, extent=coords['bounds'],
             interpolation='nearest', cmap=cmap)
  plt.title('Cells ' + plane)
  filename = SUBDIRECTORY + 'cells-' + plane + '.png'
  fig.savefig(filename, bbox_inches='tight')
  plt.close(fig)



def plot_materials(geometry, plane='xy', offset=0., gridsize=250):

  global SUBDIRECTORY

  # Make directory if it does not exist
  if not os.path.exists(SUBDIRECTORY):
    os.makedirs(SUBDIRECTORY)

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
          'is not a float'.format(offset)
    raise ValueError(msg)

  print('Plotting the Materials...')

  # Get the number of Cells filled with Materials
  materials = geometry.getAllMaterials()
  num_materials = len(materials)

  # Create array of equally spaced randomized floats as a color map for plots
  # Seed the NumPy random number generator to ensure reproducible color maps
  numpy.random.seed(1)
  colors = np.linspace(0., 1., num_materials, endpoint=False)
  numpy.random.shuffle(colors)

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

      # If we did not find a Cell for this region, use a NaN "bad" number color
      if cell is None:
        surface[j][i] = np.nan
      else:
        surface[j][i] = colors[cell._fill._id % num_materials]

  # Make Matplotlib color "bad" numbers (ie, NaN, INF) with transparent pixels
  cmap = plt.get_cmap('spectral')
  cmap.set_bad(alpha=0.0)

  # Plot a 2D color map of the flat source regions
  fig = plt.figure()
  surface = np.flipud(surface)
  plt.imshow(surface, extent=coords['bounds'],
             interpolation='nearest', cmap=cmap)
  plt.title('Materials ' + plane)
  filename = SUBDIRECTORY + 'materials-' + plane + '.png'
  fig.savefig(filename, bbox_inches='tight')
  plt.close(fig)


def plot_regions(geometry, plane='xy', offset=0., gridsize=250):

  global SUBDIRECTORY

  # Make directory if it does not exist
  if not os.path.exists(SUBDIRECTORY):
    os.makedirs(SUBDIRECTORY)

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
          'is not a float'.format(offset)
    raise ValueError(msg)

  print('Plotting the Regions...')

  # Initialize the offsets used for computing region IDs
  geometry.initializeCellOffsets()
  num_regions = geometry._num_regions

  # Create array of equally spaced randomized floats as a color map for plots
  # Seed the NumPy random number generator to ensure reproducible color maps
  numpy.random.seed(1)
  colors = np.linspace(0., 1., num_regions, endpoint=False)
  numpy.random.shuffle(colors)

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

      # If we did not find a region for this region, use a NaN "bad" number color
      if np.isnan(region_id):
        surface[j][i] = region_id
      else:
       surface[j][i] = colors[region_id % num_regions]

  # Make Matplotlib color "bad" numbers (ie, NaN, INF) with transparent pixels
  cmap = plt.get_cmap('spectral')
  cmap.set_bad(alpha=0.0)

  # Plot a 2D color map of the flat source regions
  fig = plt.figure()
  surface = np.flipud(surface)
  plt.imshow(surface, extent=coords['bounds'],
             interpolation='nearest', cmap=cmap)
  plt.title('Regions ' + plane)
  filename = SUBDIRECTORY + 'regions-' + plane + '.png'
  fig.savefig(filename, bbox_inches='tight')
  plt.close(fig)


def plot_neighbor_cells(geometry, plane='xy', offset=0.,
                        gridsize=250, unique=False):

  global SUBDIRECTORY

  # Make directory if it does not exist
  if not os.path.exists(SUBDIRECTORY):
    os.makedirs(SUBDIRECTORY)

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
          'is not a float'.format(offset)
    raise ValueError(msg)

  print('Plotting the Neighbor Cells...')

  # Initialize the offsets used for computing region IDs
  geometry.initializeCellOffsets()

  # Build the neighbor Cells/Universes
  geometry.buildNeighbors()

  if unique:
    num_neighbors = geometry._num_unique_neighbors
  else:
    num_neighbors = geometry._num_neighbors

  # Create array of equally spaced randomized floats as a color map for plots
  # Seed the NumPy random number generator to ensure reproducible color maps
  numpy.random.seed(1)
  colors = np.linspace(0., 1., num_neighbors, endpoint=False)
  numpy.random.shuffle(colors)

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

      # If we did not find a region for this region, use a NaN "bad" number color
      if np.isnan(region_id):
        surface[j][i] = region_id
      elif unique:
        neighbors = geometry._regions_to_unique_neighbors[region_id]
        neighbor_id = geometry._unique_neighbor_ids[neighbors]
        surface[j][i] = colors[neighbor_id % num_neighbors]
      else:
        neighbors = geometry._regions_to_neighbors[region_id]
        neighbor_id = geometry._neighbor_ids[neighbors]
        surface[j][i] = colors[neighbor_id % num_neighbors]

  # Make Matplotlib color "bad" numbers (ie, NaN, INF) with transparent pixels
  cmap = plt.get_cmap('spectral')
  cmap.set_bad(alpha=0.0)

  # Plot a 2D color map of the neighbor cells
  fig = plt.figure()
  surface = np.flipud(surface)
  plt.imshow(surface, extent=surf['bounds'],
             interpolation='nearest', cmap=cmap)

  if unique:
    plt.title('Unique Neighbor Cells ' + plane)
    filename = SUBDIRECTORY + 'unique-neighbor-cells-' + plane + '.png'
  else:
    plt.title('Neighbor Cells ' + plane)
    filename = SUBDIRECTORY + 'neighbor-cells-' + plane + '.png'

  fig.savefig(filename, bbox_inches='tight')
  plt.close(fig)


def plot_segments(rays, geometry, plane='xy', offset=0.,
                  gridsize=250, linewidths=1):

  global SUBDIRECTORY

  # Make directory if it does not exist
  if not os.path.exists(SUBDIRECTORY):
    os.makedirs(SUBDIRECTORY)

  # Error checking
  for ray in rays:
    if not isinstance(ray, Ray):
      msg = 'Unable to plot the segments since input does not ' \
            'completely consist of ray objects'
      raise ValueError(msg)

  if not isinstance(geometry, Geometry):
    msg = 'Unable to plot the segments since input was not ' \
          'a Geometry class object'
    raise ValueError(msg)

  if not is_integer(gridsize):
    msg = 'Unable to plot the segments since the gridsize {0} is' \
          'is not an integer'.format(gridsize)
    raise ValueError(msg)

  if gridsize <= 0:
    msg = 'Unable to plot the segments with a negative ' \
          'gridsize {0}'.format(gridsize)
    raise ValueError(msg)

  if not is_integer(linewidths):
    msg = 'Unable to plot the segments since the linewidths {0} is' \
          'is not an integer'.format(linewidths)
    raise ValueError(msg)

  if linewidths <= 0:
    msg = 'Unable to plot the segments with a negative ' \
          'linewidths {0}'.format(linewidths)
    raise ValueError(msg)

  if plane not in ['xy', 'xz', 'yz']:
    msg = 'Unable to plot the segments with an invalid ' \
        'plane {0}. Plane options xy, xz, yz'.format(plane)
    raise ValueError(msg)

  if not is_float(offset):
    msg = 'Unable to plot the neighbor cells since the offset {0} is' \
          'is not a float'.format(offset)
    raise ValueError(msg)

  print('Plotting the Segments...')

  # Get the number of regions
  num_regions = geometry._num_regions

  # Create array of equally spaced randomized floats as a color map for plots
  # Seed the NumPy random number generator to ensure reproducible color maps
  numpy.random.seed(1)
  color_map = list()
  for i in xrange(num_regions):
    color_map.append(np.random.rand(4))

  # Initialize a NumPy array for the segment colors
  colors = list()

  # Retrieve the pixel coordinates
  coords = get_pixel_coords(geometry, plane, offset, gridsize)

  # Generate start and end points for segments and assign color by region id
  segments = list()
  for ray in rays:
    start = Point()
    start.setCoords(ray._point._coords)
    dir = ray._direction.toPolar()
    for segment in xrange(ray._num_segments):
      colors.append(color_map[ray._segments[segment]._region_id % num_regions])
      length = ray._segments[segment]._length
      x = start._coords[0] + length*np.sin(dir[2])*np.cos(dir[1])
      y = start._coords[1] + length*np.sin(dir[2])*np.sin(dir[1])
      z = start._coords[2] + length*np.cos(dir[2])
      end = np.array([x, y, z])
      if plane == 'xy':
        segments.append([start._coords[:2], end[:2]])
        start = Point()
        start.setCoords(end + TINY_BIT*ray._direction._comps)
      elif plane == 'xz':
        segments.append([start._coords[::2], end[::2]])
        start = Point()
        start.setCoords(end + TINY_BIT*ray._direction._comps)
      elif plane == 'yz':
        segments.append([start._coords[1:], end[1:]])
        start = Point()
        start.setCoords(end + TINY_BIT*ray._direction._comps)

  colors = np.array(colors)

  # Plot a 2D color map of the segments
  lc = matplotlib.collections.LineCollection(segments, colors=colors, linewidths=linewidths)
  fig, ax = plt.subplots()
  ax.add_collection(lc)
  plt.xlim(coords['x'][0], coords['x'][-1])
  plt.ylim(coords['y'][0], coords['y'][-1])
  plt.title('Segments ' + plane)
  ax.margins(0)
  filename = SUBDIRECTORY + 'segments-' + plane + '.png'
  fig.savefig(filename, bbox_inches='tight')

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
