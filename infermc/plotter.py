from infermc.process import *
from infermc.multigroupxs import *
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


def scatter_multigroup_xs(extractor, xs_type, domain_types=['distribcell'],
                          energy_groups=[1,2], colors=['domain_type'],
                          legend=False, filename='multigroup-xs'):

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
    if not color in ['domain_type', 'material', 'cell',
                     'distribcell', 'universe', 'neighbors']:
      msg = 'Unable to scatter plot cross-sections with color {0} which ' \
            'is not a supported type'.format(color)
      raise ValueError(msg)

  if not isinstance(legend, (bool, np.bool)):
    msg = 'Unable to scatter plot cross-sections with legend {0} which is ' \
          'not a boolean type'.format(legend)
    raise ValueError(msg)

  elif not is_string(filename):
    msg = 'Unable to scatter plot cross-sections since {0} is not a ' \
          'valid filename string'.format(filename)
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
  materials = geometry.getAllMaterials()
  cells = geometry.getAllMaterialCells()
  universes = geometry.getAllMaterialUniverses()

  num_regions = geometry._num_regions
  num_materials = len(materials)
  num_cells = len(cells)
  num_universes = len(universes)


  # Build neighbor Cells/Universes if we will use it for coloring
  if 'neighbors' in colors:
    geometry.buildNeighbors()


  domain_colors = np.random.rand(len(domain_types))

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
          data[domain_type][region, 2] = cells[cell_id]._fill._id
        elif colors[i] == 'cell':
          data[domain_type][region, 2] = cell_id
        elif colors[i] == 'distribcell':
          data[domain_type][region, 2] = region
        elif colors[i] == 'universe':
          data[domain_type][region, 2] = universes[universe_id]
        elif colors[i] == 'neighbors':
          neighbors = geometry.getUniqueNeighbors(region)
          data[domain_type][region, 2] = hash(tuple(neighbors))
        else:
          data[domain_type][region, 2] = domain_colors[i]


    elif domain_type == 'material':

      # Data is indexed by: (group #1, group #2, color)
      data[domain_type] = np.zeros((num_materials, 3))

      for j, material in enumerate(materials):
        xs = extractor.getMultiGroupXS(xs_type, material._id, 'material')
        data[domain_type][j, 0] = xs.getXS(energy_groups[0], material._id)
        data[domain_type][j, 1] = xs.getXS(energy_groups[1], material._id)

        # Find color!!!
        if colors[i] == 'material':
          data[domain_type][j, 2] = material._id
        else:
          data[domain_type][j, 2] = domain_colors[i]


    elif domain_type == 'cell':

      # Data is indexed by: (group #1, group #2, color)
      data[domain_type] = np.zeros((num_cells, 3))

      for j, cell in enumerate(cells):
        xs = extractor.getMultiGroupXS(xs_type, cell._id, 'cell')
        data[domain_type][j, 0] = xs.getXS(energy_groups[0], cell._id)
        data[domain_type][j, 1] = xs.getXS(energy_groups[1], cell._id)

        # Find color!!!
        if colors[i] == 'cell':
          data[domain_type][j, 2] = cell._id
        elif colors[i] == 'material':
          data[domain_type][j, 2] = cell._fill._id
        else:
          data[domain_type][j, 2] = domain_colors[i]


    elif domain_type == 'universe':

      # Data is indexed by: (group #1, group #2, color)
      data[domain_type] = np.zeros((num_universes, 3))

      for j, universe in enumerate(universes):
        xs = extractor.getMultiGroupXS(xs_type, universe._id, 'universe')
        data[domain_type][j, 0] = xs.getXS(energy_groups[0], universe._id)
        data[domain_type][j, 1] = xs.getXS(energy_groups[1], universe._id)

        # Find color!!!
        if colors[i] == 'universe':
          data[domain_type][j, 2] = universe._id
        else:
          data[domain_type][j, 2] = domain_colors[i]


  # FIXME: Deal with scatter matrix!!


  fig = plt.figure()
  ax = fig.add_subplot(111)

  if 'distribcell' in domain_types:
    plt.scatter(data['distribcell'][:,0], data['distribcell'][:,1],
                c=data['distribcell'][:,2], edgecolors='k')

  if 'material' in domain_types:
    plt.scatter(data['material'][:,0], data['material'][:,1],
                c=data['material'][:,2], edgecolors='k', s=80)

  if 'cell' in domain_types:
    plt.scatter(data['material'][:,0], data['material'][:,1],
                color=data['material'][:,2], edgecolors='k', s=80)

  if 'universe' in domain_types:
    plt.scatter(data['universe'][:,0], data['universe'][:,1],
                color=data['universe'][:,2], edgecolors='k', s=80)

  plt.xlabel('Group {0} [cm^-1]'.format(energy_groups[0]))
  plt.ylabel('Group {0} [cm^-1]'.format(energy_groups[1]))
  plt.title('{0} {1} Cross-Section'.format(domain_type.capitalize(),
                                           xs_type.capitalize()))
  plt.grid()
  filename = directory + '/' + filename + '.png'
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