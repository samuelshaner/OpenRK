import matplotlib
from opencsg import *

# force headless backend, or set 'backend' to 'Agg'
# in your ~/.matplotlib/matplotlibrc
#matplotlib.use('Agg')

import matplotlib.pyplot as plt

# Force non-interactive mode, or set 'interactive' to False
# in your ~/.matplotlib/matplotlibrc
#plt.ioff()

from infermc.process import XSTallyExtractor
from infermc.multigroupxs import *
from openmc.statepoint import StatePoint
import numpy as np
import os

# Type-checking support
from typecheck import accepts, Or, Exact


## A static variable for the output directory in which to save plots
DIRECTORY = "plots/"

## List of supported file formats for saving Matplotlib figures
FILE_FORMATS = ['png', 'jpg', 'pdf', 'svg', 'eps', 'pkl']


# Build color maps
@accepts(Geometry)
def get_color_maps(geometry):

  # Initialize the offsets used for computing region IDs
  geometry.initializeCellOffsets()

  # Build the neighbor Cells/Universes
  geometry.buildNeighbors()

  materials = geometry.getAllMaterials()
  cells = geometry.getAllMaterialCells()
  universes = geometry.getAllMaterialUniverses()

  num_materials = len(materials)
  num_cells = len(cells)
  num_universes = len(universes)
  num_regions = geometry._num_regions
  num_neighbors = geometry._num_neighbors
  num_unique_neighbors = geometry._num_unique_neighbors

  # Create arrays of equally spaced randomized floats as a color map for plots
  # Seed the NumPy random number generator to ensure reproducible color maps

  # Initialize dictionary for all color maps
  # Keys    - color type string (ie, 'materials', 'cells', 'neighbors', etc.)
  # Values  - randomized NumPy array of floats in [0, 1)
  color_maps = dict()

  # Create color map for Materials
  np.random.seed(1)
  color_map = np.linspace(0., 1., num_materials, endpoint=False)
  np.random.shuffle(color_map)
  color_maps['material'] = color_map

  # Create color map for Cells
  np.random.seed(1)
  color_map = np.linspace(0., 1., num_cells, endpoint=False)
  np.random.shuffle(color_map)
  color_maps['cell'] = color_map

  # Create color map for Universes
  np.random.seed(1)
  color_map = np.linspace(0., 1., num_universes, endpoint=False)
  np.random.shuffle(color_map)
  color_maps['universe'] = color_map

  # Create color map for Regions
  np.random.seed(1)
  color_map = np.linspace(0., 1., num_regions, endpoint=False)
  np.random.shuffle(color_map)
  color_maps['region'] = color_map

  # Create color map for neighbor Cells/Universes
  np.random.seed(1)
  color_map = np.linspace(0., 1., num_neighbors, endpoint=False)
  np.random.shuffle(color_map)
  color_maps['neighbors'] = color_map

  # Create color map for neighbor Cells/Universes
  np.random.seed(1)
  color_map = np.linspace(0., 1., num_unique_neighbors, endpoint=False)
  np.random.shuffle(color_map)
  color_maps['unique neighbors'] = color_map

  return color_maps


@accepts(XSTallyExtractor, str, [str],
         (int, int), [str], str, str)
def scatter_multigroup_xs(extractor, xs_type, domain_types=['distribcell'],
                          energy_groups=(1,2), colors=['domain_type'],
                          filename='multigroup-xs', extension='png'):

  global DIRECTORY

  # Make directory if it does not exist
  if not os.path.exists(DIRECTORY):
    os.makedirs(DIRECTORY)

  # Creat a list of colors the length of the number of domains to plot
  if len(colors) == 1:
    colors = [colors[0] for i in range(len(domain_types))]

  # Replace "domain_type" color with corresponding domain type in domain_types
  for i, color in enumerate(colors):
    if color == 'domain_type':
      colors[i] = domain_types[i]

  geometry = extractor._opencsg_geometry

  # Get color maps for each domain within each domain type
  color_maps = get_color_maps(geometry)

  # Create a color for each type of domain
  domain_colors = np.linspace(0, 1, len(domain_types), endpoint=False)

  # Get lists of all Materials, Cells, Universes
  materials = geometry.getAllMaterials()
  cells = geometry.getAllMaterialCells()
  universes = geometry.getAllMaterialUniverses()

  # Get the number of Materials, Cells, Universes, etc.
  num_materials = len(materials)
  num_cells = len(cells)
  num_universes = len(universes)
  num_regions = geometry._num_regions
  num_neighbors = geometry._num_neighbors
  num_unique_neighbors = geometry._num_unique_neighbors

  # Initialize an empty dictionary of data
  data = dict()


  for i, domain_type in enumerate(domain_types):

    # Get the damn data
    if domain_type == 'distribcell':

      # Data is indexed by: (group #1, group #2, color)
      data[domain_type] = np.zeros((num_regions, 3))

      for region in range(num_regions):
        path = extractor.getPath(region)
        cell_id = path[-1]
        universe_id = path[-2]

        xs = extractor.getMultiGroupXS(xs_type, cell_id, 'distribcell')
        data[domain_type][region, 0] = xs.getXS([energy_groups[0]], [region])
        data[domain_type][region, 1] = xs.getXS([energy_groups[1]], [region])

        # Find color!!!
        if colors[i] == 'material':
          material_id = cells[cell_id]._fill._id
          color = color_maps[colors[i]][material_id % num_materials]
          data[domain_type][region, 2] = color

        elif colors[i] == 'cell':
          color = color_maps[colors[i]][cell_id % num_cells]
          data[domain_type][region, 2] = color

        elif colors[i] == 'distribcell':
          color = color_maps[colors[i]][region % num_regions]
          data[domain_type][region, 2] = color

        elif colors[i] == 'universe':
          color = color_maps[colors[i]][universe_id % num_universes]
          data[domain_type][region, 2] = color

        elif colors[i] == 'neighbors':
          neighbors = geometry._regions_to_neighbors[region]
          neighbor_id = geometry._neighbor_ids[neighbors]
          color = color_maps[colors[i]][neighbor_id % num_neighbors]
          data[domain_type][region, 2] = color

        elif colors[i] == 'unique neighbors':
          unique_neighbors = geometry._regions_to_unique_neighbors[region]
          neighbor_id = geometry._unique_neighbor_ids[unique_neighbors]
          color = color_maps[colors[i]][neighbor_id % num_unique_neighbors]
          data[domain_type][region, 2] = color

        else:
          data[domain_type][region, 2] = domain_colors[i]


    elif domain_type == 'material':

      # Data is indexed by: (group #1, group #2, color)
      data[domain_type] = np.zeros((num_materials, 3))

      for j, material_id in enumerate(materials):
        xs = extractor.getMultiGroupXS(xs_type, material_id, 'material')
        data[domain_type][j, 0] = xs.getXS([energy_groups[0]], [material_id])
        data[domain_type][j, 1] = xs.getXS([energy_groups[1]], [material_id])

        # Find color!!!
        if colors[i] == 'material':
          color = color_maps[colors[i]][material_id % num_materials]
          data[domain_type][j, 2] = color

        else:
          data[domain_type][j, 2] = domain_colors[i]


    elif domain_type == 'cell':

      # Data is indexed by: (group #1, group #2, color)
      data[domain_type] = np.zeros((num_cells, 3))

      for j, cell_id in enumerate(cells):
        xs = extractor.getMultiGroupXS(xs_type, cell_id, 'cell')
        data[domain_type][j, 0] = xs.getXS([energy_groups[0]], [cell_id])
        data[domain_type][j, 1] = xs.getXS([energy_groups[1]], [cell_id])

        # Find color!!!
        if colors[i] == 'material':
          material_id = cells[cell_id]._fill._id
          color = color_maps[colors[i]][material_id % num_materials]
          data[domain_type][j, 2] = color

        elif colors[i] == 'cell':
          color = color_maps[colors[i]][cell_id % num_cells]
          data[domain_type][j, 2] = color

        else:
          data[domain_type][j, 2] = domain_colors[i]


    elif domain_type == 'universe':

      # Data is indexed by: (group #1, group #2, color)
      data[domain_type] = np.zeros((num_universes, 3))

      for j, universe_id in enumerate(universes):
        xs = extractor.getMultiGroupXS(xs_type, universe_id, 'universe')
        data[domain_type][j, 0] = xs.getXS([energy_groups[0]], [universe_id])
        data[domain_type][j, 1] = xs.getXS([energy_groups[1]], [universe_id])

        # Find color!!!
        if colors[i] == 'universe':
          color = color_maps[colors[i]][universe_id % num_universes]
          data[domain_type][j, 2] = color

        else:
          data[domain_type][j, 2] = domain_colors[i]


  # FIXME: Deal with scatter matrix!!


  fig = plt.figure()
  ax = fig.add_subplot(111)

  if 'distribcell' in domain_types:
    plt.scatter(data['distribcell'][:,0], data['distribcell'][:,1],
                c=data['distribcell'][:,2], edgecolors='k', picker=True)

  if 'material' in domain_types:
    plt.scatter(data['material'][:,0], data['material'][:,1],
                c=data['material'][:,2], edgecolors='k', s=80, picker=True)

  if 'cell' in domain_types:
    plt.scatter(data['material'][:,0], data['material'][:,1],
                color=data['material'][:,2], edgecolors='k', s=80, picker=True)

  if 'universe' in domain_types:
    plt.scatter(data['universe'][:,0], data['universe'][:,1],
                color=data['universe'][:,2], edgecolors='k', s=80, picker=True)

  plt.xlabel('Group {0} [cm^-1]'.format(energy_groups[0]))
  plt.ylabel('Group {0} [cm^-1]'.format(energy_groups[1]))
  plt.title('{0} {1} Cross-Section'.format(domain_type.capitalize(),
                                           xs_type.capitalize()))
  plt.grid()

  filename = DIRECTORY + '/' + filename + '.' + extension

  if extension is 'pkl':
    import pickle
    ax = plt.subplot(111)
    pickle.dump(ax, file(filename, 'w'))
  else:
    plt.savefig(filename, bbox_inches='tight')

  plt.close(fig)