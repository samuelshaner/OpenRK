import matplotlib

# force headless backend, or set 'backend' to 'Agg'
# in your ~/.matplotlib/matplotlibrc
matplotlib.use('Agg')

import matplotlib.pyplot as plt

# Force non-interactive mode, or set 'interactive' to False
# in your ~/.matplotlib/matplotlibrc
plt.ioff()

import openmc
import opencsg
from infermc.process import XSTallyExtractor, MicroXSTallyExtractor
import numpy as np
import os


# A static variable for the output directory in which to save plots
DIRECTORY = "plots/"

# List of supported file formats for saving Matplotlib figures
FILE_FORMATS = ['png', 'jpg', 'pdf', 'svg', 'eps', 'pkl']

# Dictionary of data point sizes for each domain type
SCATTER_SIZES = dict()
SCATTER_SIZES['distribcell'] = 20
SCATTER_SIZES['material'] = 100
SCATTER_SIZES['cell'] = 100
SCATTER_SIZES['universe'] = 100


# Build color maps
def get_color_maps(geometry):

  # Initialize the offsets used for computing region IDs
  geometry.initializeCellOffsets()

  # Build the neighbor Cells/Universes
  geometry.buildNeighbors()

  materials = geometry.getAllMaterials()
  cells = geometry.getAllMaterialCells()
  universes = geometry.getAllMaterialUniverses()

#  num_materials = len(materials)
  num_materials = max(materials.keys()) - min(materials.keys()) + 1
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
  range = dict()

  # Create color map for Materials
  np.random.seed(1)
  color_map = np.linspace(0., 1., num_materials, endpoint=False)
  np.random.shuffle(color_map)
  color_maps['material'] = color_map
  range['material'] = num_materials

  # Create color map for Cells
  np.random.seed(1)
  color_map = np.linspace(0., 1., num_cells, endpoint=False)
  np.random.shuffle(color_map)
  color_maps['cell'] = color_map
  range['cell'] = num_cells

  # Create color map for Distribcells
  np.random.seed(1)
  color_map = np.linspace(0., 1., num_cells, endpoint=False)
  np.random.shuffle(color_map)
  color_maps['distribcell'] = color_map
  range['distribcell'] = num_cells

  # Create color map for Universes
  np.random.seed(1)
  color_map = np.linspace(0., 1., num_universes, endpoint=False)
  np.random.shuffle(color_map)
  color_maps['universe'] = color_map
  range['universe'] = num_universes

  # Create color map for Regions
  np.random.seed(1)
  color_map = np.linspace(0., 1., num_regions, endpoint=False)
  np.random.shuffle(color_map)
  color_maps['region'] = color_map
  range['region'] = num_regions

  # Create color map for neighbor Cells/Universes
  np.random.seed(1)
  color_map = np.linspace(0., 1., num_neighbors, endpoint=False)
  np.random.shuffle(color_map)
  color_maps['neighbors'] = color_map
  range['neighbors'] = num_neighbors

  # Create color map for neighbor Cells/Universes
  np.random.seed(1)
  color_map = np.linspace(0., 1., num_unique_neighbors, endpoint=False)
  np.random.shuffle(color_map)
  color_maps['unique neighbors'] = color_map
  range['unique neighbors'] = num_unique_neighbors

  return color_maps, range


def scatter_multigroup_xs(extractor, xs_type, domain_types=['distribcell'],
                          energy_groups=(1,2), colors=['domain_type'],
                          filename='multigroup-xs', extension='png'):

  global DIRECTORY
  global SCATTER_SIZES

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
  color_maps, range = get_color_maps(geometry)

  # Create a color for each type of domain
  domain_colors = np.linspace(0, 1, len(domain_types), endpoint=False)

  # Get lists of all Materials, Cells, Universes
  domains = dict()
  domains['material'] = geometry.getAllMaterials()
  domains['cell'] = geometry.getAllMaterialCells()
  domains['universe'] = geometry.getAllMaterialUniverses()
  domains['distribcell'] = geometry.getAllMaterialCells()

  # Get the number of Materials, Cells, Universes, etc.
  num_domains = dict()
  num_domains['material'] = len(domains['material'])
  num_domains['cell'] = len(domains['cell'])
  num_domains['universe'] = len(domains['universe'])
  num_domains['distribcell'] = geometry._num_regions
  num_domains['neighbors'] = geometry._num_neighbors
  num_domains['unique neighbors'] = geometry._num_unique_neighbors

  # Initialize an empty dictionary of data for each type of domain
  data = dict()

  for i, domain_type in enumerate(domain_types):

    # Data is indexed by: (group #1, group #2, color)
    data[domain_type] = np.zeros((num_domains[domain_type], 3))

    # Counter for the number of subdomains within each domain
    subdomain_counter = 0

    for domain_id in domains[domain_type]:

      # Get the MultiGroupXS object for this domain
      xs = extractor.getMultiGroupXS(xs_type, domain_id, domain_type)

      # Get subdomains IDs and their corresponding offsets and convert
      # to indices for the data array
      subdomains = xs.getSubDomains()
      offsets = np.array(xs.getSubDomainOffsets()) + subdomain_counter
      subdomain_counter += len(offsets)

      # Get the cross-section data for all subdomain and store to the array
      data[domain_type][offsets, 0:2] = np.squeeze(xs.getXS(energy_groups))

      if colors[i] == 'neighbors':
        neighbors = map(geometry._regions_to_neighbors.get, subdomains)
        neighbor_ids = np.array(map(geometry._neighbor_ids.get, neighbors))
        color = color_maps[colors[i]][neighbor_ids % range[colors[i]]]
        data[domain_type][offsets, 2] = color

      elif colors[i] == 'unique neighbors':
        neighbors = map(geometry._regions_to_unique_neighbors.get, subdomains)
        neighbor_ids = np.array(map(geometry._unique_neighbor_ids.get, neighbors))
        color = color_maps[colors[i]][neighbor_ids % range[colors[i]]]
        data[domain_type][offsets, 2] = color

      elif colors[i] in ['material', 'cell', 'universe', 'distribcell']:
        color_id = xs._colors[colors[i]]
        color = color_maps[colors[i]][color_id % range[colors[i]]]
        data[domain_type][offsets, 2] = color

      else:
        data[domain_type][offsets, 2] = domain_colors[i]

  fig = plt.figure()

  for i, domain_type in enumerate(domain_types):
    plt.scatter(data[domain_type][:,0], data[domain_type][:,1],
                c=data[domain_type][:,2], edgecolors='k',
                s=SCATTER_SIZES[domain_type], picker=True)

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


def scatter_micro_xs(extractor, xs_type, nuclide, domain_types=['distribcell'],
                     energy_groups=(1,2), filename='micro-xs', extension='png'):

  global DIRECTORY
  global SCATTER_SIZES

  # Make directory if it does not exist
  if not os.path.exists(DIRECTORY):
    os.makedirs(DIRECTORY)

  multigroup_xs = extractor._multigroup_xs

  fig = plt.figure()
  legend = list()

  # Get color maps for each domain within each domain type
  color_maps, range = get_color_maps(extractor._opencsg_geometry)

  for domain_type in domain_types:

    colors = iter(matplotlib.cm.hsv(np.linspace(0, 1, range[domain_type])))

    for domain_id in multigroup_xs[domain_type].keys():

      xs = extractor._multigroup_xs[domain_type][domain_id][xs_type]

      if not xs.containsNuclide(nuclide):
        continue

      data = xs.getXS(groups=energy_groups, nuclides=[nuclide])
      plt.scatter(data[:,0,:].ravel(), data[:,1,:].ravel(), color=next(colors),
                  s=SCATTER_SIZES[domain_type], edgecolors='k', picker=True)
      legend.append('{0} {1}'.format(domain_type.capitalize(), domain_id))

  plt.xlabel('Group {0} [barns]'.format(energy_groups[0]))
  plt.ylabel('Group {0} [barns]'.format(energy_groups[1]))
  plt.title('{0} {1} Cross-Section'.format(nuclide._name, xs_type.capitalize()))
  plt.grid()
#  plt.legend(legend)

  filename = DIRECTORY + '/' + filename + '.' + extension

  if extension is 'pkl':
    import pickle
    ax = plt.subplot(111)
    pickle.dump(ax, file(filename, 'w'))
  else:
    plt.savefig(filename, bbox_inches='tight')

  plt.close(fig)
