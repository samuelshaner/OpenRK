__author__ = 'Will Boyd'
__email__ = 'wboyd@mit.edu'


import matplotlib
from opencg import *

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


def get_unique_integers(data):
  '''Replace unique values in array with integers from monotonic sequence'''

  # Inspired by the following post on StackOverflow:
  # http://stackoverflow.com/questions/15709169/numpy-replace-groups-of-elements-with-integers-incrementally

  values, indices, inverse = np.unique(data, True, True)
  values, inverse = np.unique(indices[inverse], False, True)
  inverse = np.reshape(inverse, data.shape)
  return inverse


def plot_cells(geometry, plane='xy', offset=0., gridsize=250,
               xlim=None, ylim=None, zlim=None):

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
    msg = 'Unable to plot the cells since the gridsize {0} ' \
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
    msg = 'Unable to plot the cells since the offset {0} ' \
          'is not a float'.format(offset)
    raise ValueError(msg)

  print('Plotting the Cells...')

  # Initialize a NumPy array for the surface colors
  surface = numpy.zeros((gridsize, gridsize), dtype=np.int64)

  # Retrieve the pixel coordinates
  coords = get_pixel_coords(geometry, plane, offset, gridsize, xlim, ylim, zlim)

  # Find the flat source region IDs for each grid point
  for i in range(gridsize):
    for j in range(gridsize):

      if plane == 'xy':
        cell = geometry.findCell(x=coords['x'][i], y=coords['y'][j], z=offset)
      elif plane == 'xz':
        cell = geometry.findCell(x=coords['x'][i], y=offset, z=coords['z'][j])
      else:
        cell = geometry.findCell(x=offset, y=coords['y'][i], z=coords['z'][j])  

      # If we did not find a Cell for this region, use a -1 "bad" number color
      if cell is None:
        surface[j][i] = -1
      else:
        surface[j][i] = cell._id

  # Get the number of Material Cells in the Geometry
  material_cells = geometry.getAllMaterialCells()
  num_cells = len(material_cells)

  # Create array of all Cell IDs and randomly (but reproducibly) permute it
  cell_ids = [cell_id for cell_id in material_cells]
  numpy.random.seed(1)
  numpy.random.shuffle(cell_ids)

  # Create an array of the colors (array indices) for each value in the surface
  colors = list()
  for cell_id in surface.ravel():
    if cell_id != -1:
      colors.append(cell_ids.index(cell_id))
    else:
      colors.append(np.nan)

  colors = np.asarray(colors)
  colors.shape = surface.shape

  # Make Matplotlib color "bad" numbers (ie, NaN, INF) with transparent pixels
  cmap = plt.get_cmap('spectral')
  cmap.set_bad(alpha=0.0)

  # Plot a 2D color map of the Materials
  fig = plt.figure()
  colors = np.flipud(colors)
  plt.imshow(colors, extent=coords['bounds'],
             interpolation='nearest', cmap=cmap, vmin=0, vmax=num_cells)
  plt.title('Cells ' + plane)
  filename = SUBDIRECTORY + 'cells-' + plane + '.png'
  fig.savefig(filename, bbox_inches='tight')
  plt.close(fig)



def plot_materials(geometry, plane='xy', offset=0., gridsize=250,
                   xlim=None, ylim=None, zlim=None):

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
    msg = 'Unable to plot the materials since the gridsize {0} ' \
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
    msg = 'Unable to plot the materials since the offset {0} ' \
          'is not a float'.format(offset)
    raise ValueError(msg)

  print('Plotting the Materials...')

  # Initialize a NumPy array for the surface colors
  surface = numpy.zeros((gridsize, gridsize))

  # Retrieve the pixel coordinates
  coords = get_pixel_coords(geometry, plane, offset, gridsize, xlim, ylim, zlim)

  # Find the flat source region IDs for each grid point
  for i in range(gridsize):
    for j in range(gridsize):

      if plane == 'xy':
        cell = geometry.findCell(x=coords['x'][i], y=coords['y'][j], z=offset)
      elif plane == 'xz':
        cell = geometry.findCell(x=coords['x'][i], y=offset, z=coords['z'][j])
      else:
        cell = geometry.findCell(x=offset, y=coords['y'][i], z=coords['z'][j])  

      # If we did not find a Cell for this region, use a -1 "bad" number color
      if cell is None:
        surface[j][i] = -1
      else:
        surface[j][i] = cell._fill._id

  # Get the number of Materials in the Geometry
  materials = geometry.getAllMaterials()
  num_materials = len(materials)

  # Create array of all Material IDs and randomly (but reproducibly) permute it
  material_ids = [material_id for material_id in materials]
  numpy.random.seed(1)
  numpy.random.shuffle(material_ids)

  # Create an array of the colors (array indices) for each value in the surface
  colors = list()
  for material_id in surface.ravel():
    if material_id != -1:
      colors.append(material_ids.index(material_id))
    else:
      colors.append(np.nan)

  colors = np.asarray(colors)
  colors.shape = surface.shape

  # Make Matplotlib color "bad" numbers (ie, NaN, INF) with transparent pixels
  cmap = plt.get_cmap('spectral')
  cmap.set_bad(alpha=0.0)

  # Plot a 2D color map of the Materials
  fig = plt.figure()
  colors = np.flipud(colors)
  plt.imshow(colors, extent=coords['bounds'],
             interpolation='nearest', cmap=cmap, vmin=0, vmax=num_materials)
  plt.title('Materials ' + plane)
  filename = SUBDIRECTORY + 'materials-' + plane + '.png'
  fig.savefig(filename, bbox_inches='tight')
  plt.close(fig)


def plot_regions(geometry, plane='xy', offset=0., gridsize=250,
                   xlim=None, ylim=None, zlim=None):

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
    msg = 'Unable to plot the regions since the gridsize {0} ' \
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
    msg = 'Unable to plot the regions since the offset {0} ' \
          'is not a float'.format(offset)
    raise ValueError(msg)

  print('Plotting the Regions...')

  # Initialize a NumPy array for the surface colors
  surface = numpy.zeros((gridsize, gridsize))

  # Retrieve the pixel coordinates
  coords = get_pixel_coords(geometry, plane, offset, gridsize, xlim, ylim, zlim)

  # Find the flat source region IDs for each grid point
  for i in range(gridsize):
    for j in range(gridsize):

      if plane == 'xy':
        region_id = geometry.getRegionId(x=coords['x'][i], y=coords['y'][j], z=offset)
      elif plane == 'xz':
        region_id = geometry.getRegionId(x=coords['x'][i], y=offset, z=coords['z'][j])
      else:
        region_id = geometry.getRegionId(x=offset, y=coords['y'][i], z=coords['z'][j])

      # If we did not find a region for this region, use a -1 "bad" number color
      if np.isnan(region_id):
        surface[j][i] = -1
      else:
       surface[j][i] = region_id

  # Get the number of Regions in the Geometry
  num_regions = geometry._num_regions
  regions = np.arange(num_regions)

  # Create array of all Region IDs and randomly (but reproducibly) permute it
  region_ids = [region_id for region_id in regions]
  numpy.random.seed(1)
  numpy.random.shuffle(region_ids)

  # Create an array of the colors (array indices) for each value in the surface
  colors = list()
  for region_id in surface.ravel():
    if region_id != -1:
      colors.append(region_ids.index(region_id))
    else:
      colors.append(np.nan)

  colors = np.asarray(colors)
  colors.shape = surface.shape

  # Make Matplotlib color "bad" numbers (ie, NaN, INF) with transparent pixels
  cmap = plt.get_cmap('spectral')
  cmap.set_bad(alpha=0.0)

  # Plot a 2D color map of the regions
  fig = plt.figure()
  colors = np.flipud(colors)
  plt.imshow(colors, extent=coords['bounds'],
             interpolation='nearest', cmap=cmap, vmin=0, vmax=num_regions)
  plt.title('Regions ' + plane)
  filename = SUBDIRECTORY + 'regions-' + plane + '.png'
  fig.savefig(filename, bbox_inches='tight')
  plt.close(fig)


def plot_neighbor_cells(geometry, plane='xy', offset=0.,
                        gridsize=250, first_level=0, unique=False,
                        xlim=None, ylim=None, zlim=None):

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
    msg = 'Unable to plot the neighbor cells since the gridsize {0} ' \
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
    msg = 'Unable to plot the neighbor cells since the offset {0} ' \
          'is not a float'.format(offset)
    raise ValueError(msg)

  print('Plotting the Neighbor Cells...')

  # Initialize the offsets used for computing region IDs
  geometry.initializeCellOffsets()

  # Build the neighbor Cells/Universes
  geometry.buildNeighbors()

  # Initialize a NumPy array for the surface colors
  surface = numpy.zeros((gridsize, gridsize))

  # Retrieve the pixel coordinates
  surf = get_pixel_coords(geometry, plane, offset, gridsize, xlim, ylim, zlim)

  # Find the flat source region IDs for each grid point
  for i in range(gridsize):
    for j in range(gridsize):

      if plane == 'xy':
        region_id = geometry.getRegionId(x=surf['x'][i], y=surf['y'][j], z=offset)
      elif plane == 'xz':
        region_id = geometry.getRegionId(x=surf['x'][i], y=offset, z=surf['z'][j])
      else:
        region_id = geometry.getRegionId(x=offset, y=surf['y'][i], z=surf['z'][j])

      # If we did not find a region for this region, use a -1 "bad" number color
      if np.isnan(region_id):
        surface[j][i] = -1
      elif unique:
        neighbors_hash = geometry.getUniqueNeighborsHash(region_id, first_level)
        surface[j][i] = neighbors_hash
      else:
        neighbors_hash = geometry.getNeighborsHash(region_id, first_level)
        surface[j][i] = neighbors_hash

  # Get the number of neighbors in the plot
  num_colors = np.unique(surface).size

  # Create array of equally spaced randomized floats as a color map for plots
  # Seed the NumPy random number generator to ensure reproducible color maps
  numpy.random.seed(1)
  colors = np.arange(0., num_colors, 1, dtype=np.int64)
  numpy.random.shuffle(colors)

  # Replace neighbor IDs with monotonically increasing integers (starting at 0)
  surface = get_unique_integers(surface)
  surface = colors[surface]

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


def plot_segments(rays, geometry, plane='xy', linewidths=1):

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
  plt.title('Segments ' + plane)
  ax.margins(0)
  filename = SUBDIRECTORY + 'segments-' + plane + '.png'
  fig.savefig(filename, bbox_inches='tight')


def get_pixel_coords(geometry, plane, offset, gridsize, xlim, ylim, zlim):

  # initialize variables to be returned
  bounds = geometry.getBounds()
  xcoords = None
  ycoords = None
  zcoords = None
  coords = dict()

  if not xlim is None:
    bounds[0] = xlim[0]
    bounds[1] = xlim[1]

  if not ylim is None:
    bounds[2] = ylim[0]
    bounds[3] = ylim[1]

  if not zlim is None:
    bounds[4] = zlim[0]
    bounds[5] = zlim[1]

  if plane == 'xy':
    xcoords = np.linspace(bounds[0], bounds[1], gridsize)
    ycoords = np.linspace(bounds[2], bounds[3], gridsize)
    if offset < bounds[4] or offset > bounds[5]:
      msg = 'Unable to plot offset at z={0} as it must lie ' \
          'between the z bounds [{1},{2}]'.format(offset, \
                                                   bounds[4], bounds[5])
      raise ValueError(msg)
    del bounds[4:]
  elif plane == 'xz':
    xcoords = np.linspace(bounds[0], bounds[1], gridsize)
    zcoords = np.linspace(bounds[4], bounds[5], gridsize)
    if offset < bounds[2] or offset > bounds[3]:
      msg = 'Unable to plot offset at y={0} as it must lie ' \
          'between the y bounds [{1},{2}]'.format(offset, \
                                                   bounds[2], bounds[3])
      raise ValueError(msg)
    del bounds[2:4]
  else:
    ycoords = np.linspace(bounds[2], bounds[3], gridsize)
    zcoords = np.linspace(bounds[4], bounds[5], gridsize)
    if offset < bounds[0] or offset > bounds[1]:
      msg = 'Unable to plot offset at x={0} as it must lie ' \
          'between the x bounds [{1},{2}]'.format(offset, \
                                                   bounds[0], bounds[1])
      raise ValueError(msg)
    del bounds[:2]
  
  # add attributes to coords dictionary
  coords['x'] = xcoords
  coords['y'] = ycoords
  coords['z'] = zcoords
  coords['bounds'] = bounds

  return coords
