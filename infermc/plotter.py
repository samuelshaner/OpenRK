import matplotlib

# force headless backend, or set 'backend' to 'Agg'
# in your ~/.matplotlib/matplotlibrc
matplotlib.use('Agg')

import pylab as plt

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
SCATTER_SIZES['distribcell'] = 30
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
  num_regions = geometry._num_regions

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
  geometry = extractor._opencg_geometry
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
    all_nuclides = extractor._openmc_geometry.get_all_nuclides()
    nuclides = list()
    for nuclide_name, nuclide_tuple in all_nuclides.items():
      nuclides.append(nuclide_tuple[0])

  if domain_types == 'all':
    domain_types = infermc.domain_types

  if xs_types == 'all':
    xs_types = infermc.xs_types

  # Create a list of colors the length of the number of domains to plot
  if len(colors) == 1:
    colors = [colors[0] for i in range(len(domain_types))]

  # Replace "domain_type" color with corresponding domain type in domain_types
  for i, color in enumerate(colors):
    if color == 'domain_type':
      colors[i] = domain_types[i]

  # Get color maps for each domain within each domain type
  geometry = extractor._opencg_geometry
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
          try:
            xs = extractor.getMultiGroupXS(xs_type, domain_id, domain_type)

            if not xs.containsNuclide(nuclide):
              continue

            # Get the cross-section data for all subdomain and store to the array
            data = xs.getXS(groups=energy_groups, nuclides=[nuclide])

            # Plot the data for this domain
            plt.scatter(data[:,0,...].ravel(), data[:,1,...].ravel(), c=color,
                        edgecolors='k', s=SCATTER_SIZES[domain_type])

          except (KeyError, ValueError):
            # If the xs does not exist, continue
            pass

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


def scatter_rxn_rate_flux(multigroup_xs, filename, nuclide,
                          energy_group, uncertainties=False,
                          extension='png', xlim=None, ylim=None):

  if multigroup_xs._xs_type == 'scatter matrix':
    msg = 'Unable to make scatter plot for scattering matrices'
    raise ValueError(msg)

  elif multigroup_xs._xs_type == 'capture':
    msg = 'Unable to make scatter plot for capture cross-sections'
    raise ValueError(msg)

  elif multigroup_xs._xs_type == 'transport':
    msg = 'Unable to make scatter plot for capture cross-sections'
    raise ValueError(msg)

  global DIRECTORY
  global SCATTER_SIZES

  xs_type = multigroup_xs._xs_type
  group_index = multigroup_xs._energy_groups.getGroupIndices([energy_group])
  nuclide_index = multigroup_xs.getNuclideIndices([nuclide])

  # Extract and clean the Tally data
  tally_data, zero_indices = multigroup_xs.getTallyData()
  rxn_rate = tally_data[xs_type][0, :, group_index, nuclide_index, ...].ravel()
  flux = tally_data['flux'][0, :, group_index, ...].ravel()

  # Normalize the data so it can be viewed on a reasonable scale
  flux /= np.linalg.norm(flux)
  rxn_rate /= np.linalg.norm(rxn_rate)

  fig = plt.figure()
  ax = fig.add_subplot(111)
  axis_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)

  # Get the cross-section type and domain
  xs_type = multigroup_xs._xs_type
  domain = multigroup_xs._domain

  plt.scatter(flux, rxn_rate, s=SCATTER_SIZES['distribcell'])

  plt.xlabel('Flux')
  plt.ylabel('RXN Rate')
  plt.title('{0} {1} RXN Rate / Flux'.format(nuclide._name, xs_type.capitalize()))

  if not xlim is None:
    plt.xlim(xlim)

  if not ylim is None:
    plt.ylim(ylim)

  ax.yaxis.set_major_formatter(axis_formatter)
  ax.xaxis.set_major_formatter(axis_formatter)

  subdirectory = 'rxn-rate-flux/distribcell-{0}/{1}/{2}'.format(domain._id,
                                                               nuclide._name,
                                                               xs_type)
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


def scatter_neighbor_xs(multigroup_xs, filename, nuclide,
                        energy_groups=(1,2), uncertainties=False,
                        extension='png', xlim=None, ylim=None):

  if multigroup_xs._xs_type == 'scatter matrix':
    msg = 'Unable to make scatter plot for scattering matrices'
    raise ValueError(msg)

  global DIRECTORY
  global SCATTER_SIZES

  fig = plt.figure()
  ax = fig.add_subplot(111)
  axis_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)

  # Get the cross-section type and domain
  xs_type = multigroup_xs._xs_type
  domain = multigroup_xs._domain

  # Get the neighbor ID for each data point
  neighbors = multigroup_xs.getSubDomainNeighbors()
  unique_neighbors = np.unique(neighbors)

  # Initialize a color map for each neighbor
  np.random.seed(1)
  color_map = np.linspace(0., 1., len(unique_neighbors), endpoint=False)
  np.random.shuffle(color_map)
  color_map = itertools.cycle(matplotlib.cm.hsv(color_map))

  # Initialize lists to store plot handles and labels to create a legend
  plots = list()
  labels = list()

  for unique_neighbor in unique_neighbors:

    subdomains = multigroup_xs.getNeighborSubDomains(unique_neighbor)

    # Get the cross-section data for all subdomains
    data = multigroup_xs.getXS(groups=energy_groups,
                               subdomains=subdomains, nuclides=[nuclide])

    # If the user requested uncertainties, find 1-sigma "radius" for each point
    if uncertainties:
      x_std_dev = multigroup_xs.getXS(groups=[energy_groups[0]],
                                      subdomains=subdomains,
                                    nuclides=[nuclide], metric='std_dev')
      y_std_dev = multigroup_xs.getXS(groups=[energy_groups[1]],
                                      subdomains=subdomains,
                                    nuclides=[nuclide], metric='std_dev')
      radii = np.sqrt(x_std_dev**2 + y_std_dev**2)

    label = 'neighbor {0}'.format(unique_neighbor)
    labels.append(label)

    # Plot the data for all subdomains color-coded by neighbor ID
    if uncertainties:
      plots.append(plt.scatter(data[:,0,...].ravel(), data[:,1,...].ravel(),
                   c=next(color_map), edgecolors='k', s=radii, alpha=0.6))
    else:
      plots.append(ax.scatter(data[:,0,...].ravel(), data[:,1,...].ravel(),
                   c=next(color_map), edgecolors='k', s=SCATTER_SIZES['distribcell']))

  plt.xlabel('Group {0} [barns]'.format(energy_groups[0]))
  plt.ylabel('Group {0} [barns]'.format(energy_groups[1]))
  plt.title('{0} {1} Cross-Section'.format(nuclide._name, xs_type.capitalize()))
  plt.legend(plots, labels, scatterpoints=1, loc='best', fontsize=12)
  plt.grid()

  if not xlim is None:
    plt.xlim(xlim)

  if not ylim is None:
    plt.ylim(ylim)

  ax.yaxis.set_major_formatter(axis_formatter)
  ax.xaxis.set_major_formatter(axis_formatter)

  subdirectory = 'neighbors/distribcell-{0}/{1}/{2}'.format(domain._id,
                                                            nuclide._name,
                                                            xs_type)
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


def scatter_all_neighbors(micro_extractor, filename, energy_groups=(1,2),
                          uncertainties=False, extension='png',
                          xlim=None, ylim=None):

  distribcell_xs = micro_extractor._multigroup_xs['distribcell']

  for domain_id in distribcell_xs:
    for xs_type in distribcell_xs[domain_id]:

      if xs_type == 'scatter matrix':
        continue

      multigroup_xs = distribcell_xs[domain_id][xs_type]
      nuclides = multigroup_xs._nuclides

      for nuclide in nuclides:
        scatter_neighbor_xs(multigroup_xs, filename, nuclide, energy_groups,
                            uncertainties, extension, xlim, ylim)