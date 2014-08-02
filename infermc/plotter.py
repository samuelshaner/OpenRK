import matplotlib
from opencsg import *

# force headless backend, or set 'backend' to 'Agg'
# in your ~/.matplotlib/matplotlibrc
#matplotlib.use('Agg')

import matplotlib.pyplot as plt

# Force non-interactive mode, or set 'interactive' to False
# in your ~/.matplotlib/matplotlibrc
#plt.ioff()


from infermc.process import *
from infermc.multigroupxs import *
from opencsg.checkvalue import *
from statepoint import StatePoint
import numpy as np
import numpy.random
import os


## A static variable for the output directory in which to save plots
directory = "plots/"


## List of supported file formats for saving Matplotlib figures
file_formats = ['png', 'jpg', 'pdf', 'svg', 'eps', 'pkl']


# Build color maps
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
  numpy.random.seed(1)
  color_map = np.linspace(0., 1., num_materials, endpoint=False)
  numpy.random.shuffle(color_map)
  color_maps['material'] = color_map

  # Create color map for Cells
  numpy.random.seed(1)
  color_map = np.linspace(0., 1., num_cells, endpoint=False)
  numpy.random.shuffle(color_map)
  color_maps['cell'] = color_map

  # Create color map for Universes
  numpy.random.seed(1)
  color_map = np.linspace(0., 1., num_universes, endpoint=False)
  numpy.random.shuffle(color_map)
  color_maps['universe'] = color_map

  # Create color map for Regions
  numpy.random.seed(1)
  color_map = np.linspace(0., 1., num_regions, endpoint=False)
  numpy.random.shuffle(color_map)
  color_maps['region'] = color_map

  # Create color map for neighbor Cells/Universes
  numpy.random.seed(1)
  color_map = np.linspace(0., 1., num_neighbors, endpoint=False)
  numpy.random.shuffle(color_map)
  color_maps['neighbors'] = color_map

  # Create color map for neighbor Cells/Universes
  numpy.random.seed(1)
  color_map = np.linspace(0., 1., num_unique_neighbors, endpoint=False)
  numpy.random.shuffle(color_map)
  color_maps['unique neighbors'] = color_map

  return color_maps


def scatter_multigroup_xs(extractor, xs_type, domain_types=['distribcell'],
                          energy_groups=[1,2], colors=['domain_type'],
                          filename='multigroup-xs', extension='png'):

  if not isinstance(extractor, XSTallyExtractor):
    msg = 'Unable to scatter plot cross-sections since {0} is not a ' \
          'XSTallyExtractor class object'.format(extractor)
    raise ValueError(msg)

  elif not xs_type in xs_types:
    msg = 'Unable to scatter plot cross-sections since {0} is not a ' \
          'valid cross-section type'.format(xs_type)
    raise ValueError(msg)

  elif not isinstance(domain_types, (list, tuple, np.ndarray)):
    msg = 'Unable to scatter plot cross-sections since the domain ' \
          'types is not a Python tuple/list or NumPy array'
    raise ValueError(msg)

  for domain_type in domain_types:
    if not domain_type in ['material', 'cell', 'distribcell', 'universe']:
      msg = 'Unable to scatter plot cross-sections since {0} is not a ' \
            'valid domain type'.format(domain_type)
      raise ValueError(msg)

  if not isinstance(energy_groups, (list, tuple, np.ndarray)):
    msg = 'Unable to scatter plot cross-sections since {0} is not ' \
          'a Python tuple/list or NumPy array'.format(energy_groups)
    raise ValueError(msg)

  elif len(energy_groups) != 2:
    msg = 'Unable to scatter plot cross-sections since {0} groups were ' \
          'input but only two groups is allowed'.format(len(energy_groups))
    raise ValueError(msg)

  for group in energy_groups:
    if not is_integer(group) or group <= 0:
      msg = 'Unable to scatter plot cross-sections since {0} is not a ' \
            'positive integer energy group'.format(group)
      raise ValueError(msg)

  if not isinstance(colors, (list, tuple, np.ndarray)):
    msg = 'Unable to scatter plot cross-sections since {0} is not ' \
          'a Python tuple/list of NumPy array of colors'.format(energy_groups)
    raise ValueError(msg)

  elif len(colors) != 1 and len(colors) != len(domain_types):
    msg = 'Unable to scatter plot cross-sections since the number of colors ' \
          '{0} is not 1 or the number of domain types'.format(len(colors))
    raise ValueError(msg)

  for color in colors:
    if not color in ['domain_type', 'material', 'cell', 'distribcell'
                     'universe', 'neighbors', 'unique neighbors']:
      msg = 'Unable to scatter plot cross-sections with color {0} which ' \
            'is not a supported type'.format(color)
      raise ValueError(msg)

  if not is_string(filename):
    msg = 'Unable to scatter plot cross-sections since {0} is not a ' \
          'valid filename string'.format(filename)
    raise ValueError(msg)

  elif not is_string(extension) or not extension in file_formats:
    msg = 'Unable to scatter plot cross-sections since {0} ' \
          'is not a valid string file extension'.format(format)
    raise ValueError(msg)

  global directory

  # Make directory if it does not exist
  if not os.path.exists(directory):
    os.makedirs(directory)

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
        data[domain_type][region, 0] = xs.getXS(energy_groups[0], region)
        data[domain_type][region, 1] = xs.getXS(energy_groups[1], region)

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
        data[domain_type][j, 0] = xs.getXS(energy_groups[0], material_id)
        data[domain_type][j, 1] = xs.getXS(energy_groups[1], material_id)

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
        data[domain_type][j, 0] = xs.getXS(energy_groups[0], cell_id)
        data[domain_type][j, 1] = xs.getXS(energy_groups[1], cell_id)

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
        data[domain_type][j, 0] = xs.getXS(energy_groups[0], universe_id)
        data[domain_type][j, 1] = xs.getXS(energy_groups[1], universe_id)

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

  filename = directory + '/' + filename + '.' + extension

  if extension is 'pkl':
    import pickle
    ax = plt.subplot(111)
    pickle.dump(ax, file(filename, 'w'))
  else:
    plt.savefig(filename, bbox_inches='tight')

  plt.close(fig)





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