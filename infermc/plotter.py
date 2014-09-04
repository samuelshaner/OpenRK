import matplotlib

# force headless backend, or set 'backend' to 'Agg'
# in your ~/.matplotlib/matplotlibrc
matplotlib.use('Agg')

import matplotlib.pyplot as plt

# Force non-interactive mode, or set 'interactive' to False
# in your ~/.matplotlib/matplotlibrc
plt.ioff()

import infermc
import numpy as np
import os, itertools


# A static variable for the output directory in which to save plots
DIRECTORY = 'plots'

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

  materials = geometry.getAllMaterials()
  cells = geometry.getAllMaterialCells()
  universes = geometry.getAllMaterialUniverses()

  num_materials = len(materials)
  num_cells = len(cells)
  num_universes = len(universes)

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
  color_maps['material'] = itertools.cycle(matplotlib.cm.Spectral(color_map))

  # Create color map for Cells
  np.random.seed(1)
  color_map = np.linspace(0., 1., num_cells, endpoint=False)
  np.random.shuffle(color_map)
  color_maps['cell'] = itertools.cycle(matplotlib.cm.hsv(color_map))

  # Create color map for Distribcells
  np.random.seed(1)
  color_map = np.linspace(0., 1., num_cells, endpoint=False)
  np.random.shuffle(color_map)
  color_maps['distribcell'] = itertools.cycle(matplotlib.cm.Set1(color_map))

  # Create color map for Universes
  np.random.seed(1)
  color_map = np.linspace(0., 1., num_universes, endpoint=False)
  np.random.shuffle(color_map)
  color_maps['universe'] = itertools.cycle(matplotlib.cm.Dark2(color_map))

  return color_maps


def scatter_multigroup_xs(extractor, filename, xs_types='all',
                          domain_types='all', energy_groups=(1,2),
                          colors=['domain_type'], extension='png',
                          xlim=None, ylim=None):

  global DIRECTORY
  global SCATTER_SIZES

  if domain_types == 'all':
    domain_types = infermc.domain_types

  if xs_types == 'all':
    xs_types = infermc.xs_types

  # Creat a list of colors the length of the number of domains to plot
  if len(colors) == 1:
    colors = [colors[0] for i in range(len(domain_types))]

  # Replace "domain_type" color with corresponding domain type in domain_types
  for i, color in enumerate(colors):
    if color == 'domain_type':
      colors[i] = domain_types[i]

  # Get color maps for each domain within each domain type
  geometry = extractor._opencsg_geometry
  color_maps = get_color_maps(geometry)

  # Get lists of all Materials, Cells, Universes
  domains = dict()
  domains['material'] = geometry.getAllMaterials()
  domains['cell'] = geometry.getAllMaterialCells()
  domains['universe'] = geometry.getAllMaterialUniverses()
  domains['distribcell'] = geometry.getAllMaterialCells()

  axis_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)

  for i, xs_type in enumerate(xs_types):

    if xs_type == 'scatter matrix':
      continue

    # Create a new Matplotlib figure
    fig = plt.figure()

    for j, domain_type in enumerate(domain_types):
      for domain_id in domains[domain_type]:

        # Get the color for this domain ID
        color = next(color_maps[colors[j]])

        # Get the MultiGroupXS object for this domain
        xs = extractor.getMultiGroupXS(xs_type, domain_id, domain_type)

        # Get the cross-section data for all subdomain and store to the array
        data = xs.getXS(groups=energy_groups)

        # Get the color for this domain ID
        color = next(color_maps[colors[j]])

        # Plot the data for this domain
        plt.scatter(data[:,0,...].ravel(), data[:,1,...].ravel(), c=color,
                    edgecolors='k', s=SCATTER_SIZES[domain_type])

    plt.xlabel('Group {0} [cm^-1]'.format(energy_groups[0]))
    plt.ylabel('Group {0} [cm^-1]'.format(energy_groups[1]))
    plt.title('{0} {1} Cross-Section'.format(xs_type.capitalize()))
    plt.grid()

    if not xlim is None:
      plt.xlim(xlim)

    if not ylim is None:
      plt.ylim(ylim)

    ax = plt.subplot(111)
    ax.yaxis.set_major_formatter(axis_formatter)
    ax.xaxis.set_major_formatter(axis_formatter)

    subdirectory = 'macroxs/{0}'.format(xs_type)
    full_directory = DIRECTORY + '/' + subdirectory

    if not os.path.exists(full_directory):
      os.makedirs(full_directory)

    full_filename = full_directory + '/' + filename + '.' + extension

    if extension is 'pkl':
      import pickle
      ax = plt.subplot(111)
      pickle.dump(ax, file(full_filename, 'w'))
    else:
      plt.savefig(full_filename, bbox_inches='tight')

    plt.close(fig)


def scatter_micro_xs(extractor, filename, nuclides='all', xs_types='all',
                     domain_types='all', energy_groups=(1,2),
                     colors=['domain_type'], extension='png',
                     xlim=None, ylim=None):

  global DIRECTORY
  global SCATTER_SIZES

  if nuclides == 'all':
    all_nuclides = extractor._openmc_geometry.getAllNuclides()
    nuclides = list()
    for nuclide_name, nuclide_tuple in all_nuclides.items():
      nuclides.append(nuclide_tuple[0])

  if domain_types == 'all':
    domain_types = infermc.domain_types

  if xs_types == 'all':
    xs_types = infermc.xs_types

  # Creat a list of colors the length of the number of domains to plot
  if len(colors) == 1:
    colors = [colors[0] for i in range(len(domain_types))]

  # Replace "domain_type" color with corresponding domain type in domain_types
  for i, color in enumerate(colors):
    if color == 'domain_type':
      colors[i] = domain_types[i]

  # Get color maps for each domain within each domain type
  geometry = extractor._opencsg_geometry
  color_maps = get_color_maps(geometry)

  # Create a color for each type of domain
  domain_colors = np.linspace(0, 1, len(domain_types), endpoint=False)

  # Get lists of all Materials, Cells, Universes
  domains = dict()
  domains['material'] = geometry.getAllMaterials()
  domains['cell'] = geometry.getAllMaterialCells()
  domains['universe'] = geometry.getAllMaterialUniverses()
  domains['distribcell'] = geometry.getAllMaterialCells()

  axis_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)

  for nuclide in nuclides:
    for i, xs_type in enumerate(xs_types):

      if xs_type == 'scatter matrix':
        continue

      # Create a new Matplotlib figure for this Nuclide cross-section
      fig = plt.figure()

      for j, domain_type in enumerate(domain_types):
        for domain_id in domains[domain_type]:

          # Get the color for this domain ID
          color = next(color_maps[colors[j]])

          # Get the MultiGroupXS object for this domain
          xs = extractor.getMultiGroupXS(xs_type, domain_id, domain_type)

          if not xs.containsNuclide(nuclide):
            continue

          # Get the cross-section data for all subdomain and store to the array
          data = xs.getXS(groups=energy_groups, nuclides=[nuclide])

          # Plot the data for this domain
          plt.scatter(data[:,0,...].ravel(), data[:,1,...].ravel(), c=color,
                      edgecolors='k', s=SCATTER_SIZES[domain_type])

      plt.xlabel('Group {0} [barns]'.format(energy_groups[0]))
      plt.ylabel('Group {0} [barns]'.format(energy_groups[1]))
      plt.title('{0} {1} Cross-Section'.format(nuclide._name,
                                               xs_type.capitalize()))
      plt.grid()

      if not xlim is None:
        plt.xlim(xlim)

      if not ylim is None:
        plt.ylim(ylim)

      ax = plt.subplot(111)
      ax.yaxis.set_major_formatter(axis_formatter)
      ax.xaxis.set_major_formatter(axis_formatter)

      subdirectory = 'microxs/{0}/{1}'.format(nuclide._name, xs_type)
      full_directory = DIRECTORY + '/' + subdirectory

      if not os.path.exists(full_directory):
        os.makedirs(full_directory)

      full_filename = full_directory + '/' + filename + '.' + extension

      if extension is 'pkl':
        import pickle
        ax = plt.subplot(111)
        pickle.dump(ax, file(full_filename, 'w'))
      else:
        plt.savefig(full_filename, bbox_inches='tight')

      plt.close(fig)