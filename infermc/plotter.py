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
                          energy_groups=[1,2], filename='multigroup-xs'):

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

  elif not is_string(filename):
    msg = 'Unable to scatter plot cross-sections since {0} is not a ' \
          'valid filename string'.format(filename)
    raise ValueError(msg)

  global directory

  # Make directory if it does not exist
  if not os.path.exists(directory):
    os.makedirs(directory)

  # Get the damn data
  if 'distribcell' in domain_types:
    num_regions = extractor._opencsg_geometry._num_regions

    distribcell_data = np.zeros((num_regions, 2))

    for region in range(num_regions):
      path = extractor.getPath(region)
      cell_id = path[-1]

      xs = extractor.getMultiGroupXS(xs_type, cell_id, 'distribcell')
      distribcell_data[region, 0] = xs.getXS(energy_groups[0], region)
      distribcell_data[region, 1] = xs.getXS(energy_groups[1], region)


  # FIXME: Deal with scatter matrix!!
  # FIXME: Deal with other domain types!!

  if 'material' in domain_types:

    openmc_geometry = extractor._openmc_geometry
    openmc_materials = openmc_geometry.getAllMaterials()
    num_materials = len(openmc_materials)

    material_data = np.zeros((num_materials, 2))

    for i, material in enumerate(openmc_materials):
      xs = extractor.getMultiGroupXS(xs_type, material._id, 'material')
      material_data[i, 0] = xs.getXS(energy_groups[0], material._id)
      material_data[i, 1] = xs.getXS(energy_groups[1], material._id)


  if 'cell' in domain_types:
    exit('not yet implemented')
  if 'universe' in domain_type:
    exit('not yet implemented')


  fig = plt.figure()
  ax = fig.add_subplot(111)

  if 'distribcell' in domain_types:
    plt.scatter(distribcell_data[:,0], distribcell_data[:,1], color='b', edgecolors='k')

  if 'material' in domain_types:
    plt.scatter(material_data[:,0], material_data[:,1], color='r', edgecolors='k', s=80)

  plt.xlabel('Group {0} [cm^-1]'.format(energy_groups[0]))
  plt.ylabel('Group {0} [cm^-1]'.format(energy_groups[1]))
  plt.title('{0} {1}'.format(domain_type.capitalize(), xs_type))
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