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
SCATTER_SIZES['distribcell'] = 40
SCATTER_SIZES['material'] = 150
SCATTER_SIZES['cell'] = 150
SCATTER_SIZES['universe'] = 150

# Dictionary of data point markers for each domain type
MARKERS = dict()
MARKERS['distribcell'] = "o"
MARKERS['material'] = "s"
MARKERS['cell'] = "8"
MARKERS['universe'] = "D"

# Dictionary of data point linewidths for each domain type
LINEWIDTHS = dict()
LINEWIDTHS['distribcell'] = 0
LINEWIDTHS['material'] = 1
LINEWIDTHS['cell'] = 1
LINEWIDTHS['universe'] = 1


# Build color maps
def get_color_maps(geometry):

  # Initialize the offsets used for computing region IDs
  geometry.initializeCellOffsets()

  materials = geometry.getAllMaterials()
  material_cells = geometry.getAllMaterialCells()
  universes = geometry.getAllMaterialUniverses()

  num_materials = len(materials)
  num_cells = len(material_cells)
  num_universes = len(universes)

  # Initialize dictionary for all color maps
  # Keys    - color type string (ie, 'materials', 'cells', etc.)
  # Values  - randomized NumPy arrays of non-negative integers
  color_maps = dict()

  # Create color map for Materials
  material_ids = [material_id for material_id in materials]
  np.random.seed(1)
  np.random.shuffle(material_ids)
  color_maps['material'] = (material_ids, num_materials)

  # Create color map for Cells
  cell_ids = [cell_id for cell_id in material_cells]
  np.random.seed(1)
  np.random.shuffle(cell_ids)
  color_maps['cell'] = (cell_ids, num_cells)

  # Create color map for Distribcells
  cell_ids = [cell_id for cell_id in material_cells]
  np.random.seed(1)
  np.random.shuffle(cell_ids)
  color_maps['distribcell'] = (cell_ids, num_cells)

  # Create color map for Universes
  universe_ids = [universe_id for universe_id in universes]
  np.random.seed(1)
  np.random.shuffle(universe_ids)
  color_maps['universe'] = (universe_ids, num_universes)

  return color_maps


def scatter_multigroup_xs(extractor, filename, xs_types='all',
                          domain_types='all', energy_groups=(1,2),
                          colors=['domain_type'], extension='png',
                          xlim=None, ylim=None):

  global DIRECTORY
  global SCATTER_SIZES, LINEWIDTHS, MARKERS

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
  cmap = plt.get_cmap('spectral')

  # Get lists of all Materials, Cells, Universes
  domains = dict()
  domains['material'] = geometry.getAllMaterials()
  domains['cell'] = geometry.getAllMaterialCells()
  domains['distribcell'] = geometry.getAllMaterialCells()
  domains['universe'] = geometry.getAllMaterialUniverses()

  axis_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)

  for i, xs_type in enumerate(xs_types):

    if xs_type == 'scatter matrix':
      continue

    # Create a new Matplotlib figure
    fig = plt.figure()

    for j, domain_type in enumerate(domain_types):
      for domain_id in domains[domain_type]:

        # Get the MultiGroupXS object for this domain
        xs = extractor.getMultiGroupXS(xs_type, domain_id, domain_type)

        # Get the cross-section data for all subdomain and store to the array
        data = xs.getXS(groups=energy_groups)

        # Get the color for this domain ID
        color = color_maps[domain_type][0].index(domain_id)
        max_color = color_maps[domain_type][1]
        colors = np.ones(data.size/2) * color

        # Plot the data for this domain
        plt.scatter(data[:,0,...].ravel(), data[:,1,...].ravel(),
                    cmap=cmap, c=colors, vmin=0, vmax=max_color,
                    lw=LINEWIDTHS[domain_type], marker=MARKERS[domain_type],
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
  global SCATTER_SIZES, LINEWIDTHS, MARKERS

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
  cmap = plt.get_cmap('spectral')

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

          # Get the MultiGroupXS object for this domain
          try:
            xs = extractor.getMultiGroupXS(xs_type, domain_id, domain_type)

            if not xs.containsNuclide(nuclide):
              continue

            # Get the cross-section data for all subdomain and store to the array
            data = xs.getXS(groups=energy_groups, nuclides=[nuclide])

            # Get the color for this domain ID
            color = color_maps[domain_type][0].index(domain_id)
            max_color = color_maps[domain_type][1]
            colors = np.ones(data.size/2) * color

            # Plot the data for this domain
            plt.scatter(data[:,0,...].ravel(), data[:,1,...].ravel(),
                        cmap=cmap, c=colors, vmin=0, vmax=max_color,
                        lw=LINEWIDTHS[domain_type], marker=MARKERS[domain_type],
                        edgecolors='k', s=SCATTER_SIZES[domain_type])

          # If the xs does not exist, continue
          except (KeyError, ValueError):
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

  plt.scatter(flux, rxn_rate, s=SCATTER_SIZES['distribcell'],
              lw=LINEWIDTHS['distribcell'], marker=MARKERS['distribcell'])

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
                   lw=LINEWIDTHS['distribcell'], marker=MARKERS['distribcell'],
                   c=next(color_map), edgecolors='k', s=radii, alpha=0.6))
    else:
      plots.append(ax.scatter(data[:,0,...].ravel(), data[:,1,...].ravel(),
                   lw=LINEWIDTHS['distribcell'], marker=MARKERS['distribcell'],
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