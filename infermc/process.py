import openmc
import infermc
import numpy as np


def get_path(coords):

  # Build "path" from LocalCoords for OpenMC
  path = list()

  while coords is not None:

    # If the LocalCoords is at a Universe
    if coords._type == 'universe':
      path.append(coords._universe._id)
      path.append(coords._cell._id)

    # If the LocalCoords is at a Lattice
    else:
      # Add 1 for Fortran indexing
      lat_x = coords._lat_x + 1
      lat_y = coords._lat_y + 1

      # 3D Lattices only
      if coords._lat_z != None:
        lat_z = coords._lat_z + 1

      # 2D Lattices
      else:
        lat_z = 1

      path.append((coords._lattice._id, lat_x, lat_y, lat_z))

    # Traverse LocalCoords linked list to next lowest nested universe
    coords = coords._next

  return tuple(path)



class XSTallyExtractor(object):

  def __init__(self, statepoint=None, summary=None):

    # Initialize TallyExtractor class attributes
    self._statepoint = None
    self._summary = None
    self._openmc_geometry = None
    self._opencg_geometry = None

    # Dictionaries mapping cells/materials to Tallies
    # Keys   - Region ID
    # Values - Tally ID
    self._distribcells_to_tallies = dict()
    self._cells_to_tallies = dict()
    self._universes_to_tallies = dict()
    self._materials_to_tallies = dict()

    # Dictionary mapping Tallies to scores
    # Keys   - Tally ID
    # Values - score list
    self._tallies_to_scores = dict()

    # Dictionary mapping regions to "paths"
    # Keys   - region ID
    # Values - "path" string
    self._all_paths = dict()

    # Nested dictionary of MultiGroupXS objects
    # Keys   - domain type (e.g, 'distribcell', 'material')
    # Values - dictionary for each domain of that type
    self._multigroup_xs = dict()

    if not statepoint is None:
      self.statepoint = statepoint

    if not summary is None:
      self.summary = summary


  @property
  def statepoint(self):
    return self._statepoint


  @statepoint.setter
  def statepoint(self, statepoint):

    self._statepoint = statepoint
    self._statepoint.read_results()
    self._statepoint.compute_ci()

    if self._summary is not None:
      self._statepoint.link_with_summary(self._summary)

      # Build maps to optimize tally lookups
      self._buildTallyMaps()


  @property
  def summary(self):
    return self._summary


  @summary.setter
  def summary(self, summary):

    self._summary = summary
    self._summary.make_opencg_geometry()
    self._summary.opencg_geometry.buildNeighbors()

    # Retrieve the OpenMC Geometry from the statepoint and convert it
    # into an OpenCG geometry object using the compatibility module
    self._openmc_geometry = self._summary.openmc_geometry
    self._opencg_geometry = self._summary.opencg_geometry

    if self._statepoint is not None:
      self._statepoint.link_with_summary(self._summary)

      # Build maps to optimize tally lookups
      self._buildTallyMaps()


  def _buildTallyMaps(self):

    self._buildMaterialTallyMaps()
    self._buildUniverseTallyMaps()
    self._buildCellTallyMaps()
    self._buildDistribcellTallyMaps()


  def _buildMaterialTallyMaps(self):

    # Create a mapping of Tally locations to location IDs to Tally IDs to scores
    for tally_id, tally in self._statepoint._tallies.items():

      # Store a list of the Tally scores
      self._tallies_to_scores[tally._id] = tally._scores

      filters = tally._filters

      for filter in filters:

        if filter._type == 'material':
          material_ids = filter._bins

          # Build maps for all materials to this Tally
          for material_id in material_ids:

            # If this material_id is not already a key, create a list for it
            if not material_id in self._materials_to_tallies.keys():
              self._materials_to_tallies[material_id] = list()

            # Add the Tally's ID to the list for this material
            self._materials_to_tallies[material_id].append(tally._id)


  def _buildUniverseTallyMaps(self):

    # Create a mapping of Tally locations to location IDs to Tally IDs to scores
    for tally_id, tally in self._statepoint._tallies.items():

      # Store a list of the Tally scores
      self._tallies_to_scores[tally._id] = tally._scores

      filters = tally._filters

      for filter in filters:

        if filter._type == 'universe':
          universe_ids = filter._bins

          # Build maps for all universes to this Tally
          for universe_id in universe_ids:

            # If this universe_id is not already a key, create a list for it
            if not universe_id in self._universes_to_tallies.keys():
              self._universes_to_tallies[universe_id] = list()

            # Add the Tally's ID to the list for this universe
            self._universes_to_tallies[universe_id].append(tally._id)


  def _buildCellTallyMaps(self):

    # Create a mapping of Tally locations to location IDs to Tally IDs to scores
    for tally_id, tally in self._statepoint._tallies.items():

      # Store a list of the Tally scores
      self._tallies_to_scores[tally._id] = tally._scores

      filters = tally._filters

      for filter in filters:

        if filter._type == 'cell':
          cell_ids = filter._bins

          # Build maps for all cells to this Tally
          for cell_id in cell_ids:

            # If this distribcell_id is not already a key, create a list for it
            if not cell_id in self._cells_to_tallies.keys():
              self._cells_to_tallies[cell_id] = list()

            # Add the Tally's ID to the list for this cell
            self._cells_to_tallies[cell_id].append(tally._id)


  def _buildDistribcellTallyMaps(self):

    # Create a mapping of Tally locations to location IDs to Tally IDs to scores
    for tally_id, tally in self._statepoint._tallies.items():

      # Store a list of the Tally scores
      self._tallies_to_scores[tally._id] = tally._scores

      filters = tally._filters

      for filter in filters:

        if filter._type == 'distribcell':
          distribcell_ids = filter._bins

          # Build maps for all distribcells to this Tally
          for distribcell_id in distribcell_ids:

            # If this distribcell_id is not already a key, create a list for it
            if not distribcell_id in self._distribcells_to_tallies.keys():
              self._distribcells_to_tallies[distribcell_id] = list()

            # Add the Tally's ID to the list for this distribcell
            self._distribcells_to_tallies[distribcell_id].append(tally._id)


  def buildAllPaths(self):

    # Create a list of "paths" for each unique region in the Geometry
    self._opencg_geometry.initializeCellOffsets()

    num_regions = self._opencg_geometry._num_regions

    for region in range(num_regions):
      coord = self._opencg_geometry.findRegion(region)
      self._all_paths[region] = get_path(coord)


  def getPath(self, region):

    # If this region has not been requested before, memoize its path
    if region not in self._all_paths:
      coord = self._opencg_geometry.findRegion(region)
      self._all_paths[region] = get_path(coord)

    return self._all_paths[region]


  def buildNeighborMaps(self, first_level=0, unique=False):

    distribcell_xs = self._multigroup_xs['distribcell']
    geometry = self._opencg_geometry
    geometry.clearNeighbors()
    geometry.buildNeighbors()
    geometry.countNeighbors()

    for domain_id in distribcell_xs:
      for xs_type in distribcell_xs[domain_id]:

        multigroup_xs = distribcell_xs[domain_id][xs_type]
        multigroup_xs.setUniqueNeighbors(unique)
        subdomains = multigroup_xs.getSubDomains()

        for subdomain in subdomains:

          if unique:
            neighbor = geometry.getUniqueNeighborsHash(subdomain, first_level)
          else:
            neighbor = geometry.getNeighborsHash(subdomain, first_level)

          multigroup_xs.setSubDomainNeighbor(subdomain, neighbor)


  def getTally(self, score, filters, nuclides=[],
               estimator='tracklength', label=''):

    if self._statepoint is None:
      msg = 'Unable to get Tally since statepoint attribute has not been set'
      raise ValueError(msg)

    # Loop over the domain-to-tallies mapping to find the Tally
    tally = None
    tallies = None

    # Determine the Tally domain (e.g, Material, Cell, etc.)
    for filter in filters:

      if filter._type == 'material':
        tallies = self._materials_to_tallies
        domain_id = filter._bins[0]
        break

      elif filter._type == 'cell':
        tallies = self._cells_to_tallies
        domain_id = filter._bins[0]
        break

      elif filter._type == 'distribcell':
        tallies = self._distribcells_to_tallies
        domain_id = filter._bins[0]
        break

      elif filter._type == 'universe':
        tallies = self._universes_to_tallies
        domain_id = filter._bins[0]
        break

    # Iterate over all tallies to find the appropriate one
    for tally_id in tallies[domain_id]:

      # If the Tally score doesn't match
      if not score in self._tallies_to_scores[tally_id]:
        continue

      # Get the Tally from the StatePoint
      test_tally = self._statepoint._tallies[tally_id]

      # If the label parameter was set and doesn't match
      if not label is '' and test_tally._label != label:
        continue

      # If the estimator type ('tracklength' or 'analog') doesn't match
      if test_tally._estimator != estimator:
        continue

      # If the length of the Filters container doesn't match
      if len(filters) != len(test_tally._filters):
        continue

      contains_filters = True

      # Iterate over the Filters requested by the user
      for filter in test_tally._filters:
        if not filter in filters:
          contains_filters = False
          break


      contains_nuclides = True

      '''
      # Iterate over the Nuclides requested by the user
      for nuclide in test_tally._nuclides:
        if not nuclide in nuclides:
          contains_nuclides = False
          break
      '''

      # If the Tally contained all Filters and Nuclides, return the Tally
      if contains_filters and contains_nuclides:
        tally = test_tally
        break

    # If we did not find the Tally, return an error message
    if tally is None:
      msg = 'Unable to get Tally for score {0}'.format(score)
      raise ValueError(msg)

    return tally


  def extractAllMultiGroupXS(self, energy_groups, domain_type='distribcell', corr=False):

    for xs_type in infermc.xs_types:
      self.extractMultiGroupXS(xs_type, energy_groups, domain_type, corr)


  def extractMultiGroupXS(self, xs_type, energy_groups, domain_type='distribcell', corr=False):

    # Add nested dictionary for this domain type if needed
    if not domain_type in self._multigroup_xs.keys():
      self._multigroup_xs[domain_type] = dict()

    # Get a list of the domains for this domain type to iterate over
    if domain_type == 'material':
      domains = self._openmc_geometry.get_all_materials()

    elif domain_type == 'universe':
      domains = self._openmc_geometry.get_all_material_universes()

    elif domain_type == 'cell' or domain_type == 'distribcell':
      domains = self._openmc_geometry.get_all_material_cells()

    # Iterate and create the MultiGroupXS for each domain
    for domain in domains:

      # Add nested dictionary for this domain if needed
      if not domain._id in self._multigroup_xs[domain_type].keys():
        self._multigroup_xs[domain_type][domain._id] = dict()

      # Build the MultiGroupXS for this domain
      xs = self.createMultiGroupXS(xs_type, energy_groups, domain, domain_type, corr)

      # Store a handle to the MultiGroupXS object in the nested dictionary
      self._multigroup_xs[domain_type][domain._id][xs_type] = xs


  def createMultiGroupXS(self, xs_type, energy_groups,
                         domain, domain_type='distribcell', corr=False):

    if self._statepoint is None:
      msg = 'Unable to get cross-sections since the TallyExtractor ' \
            'statepoint attribute has not been set'
      raise ValueError(msg)

    elif self._summary is None:
      msg = 'Unable to get cross-sections since the TallyExtractor ' \
            'summary attribute has not been set'
      raise ValueError(msg)

    # Initialize a list of filters
    filters = list()

    # Create energy and domain filters to search for
    group_edges = energy_groups._group_edges
    filters.append(openmc.Filter(type='energy', bins=group_edges))
    filters.append(openmc.Filter(type=domain_type, bins=domain._id))

    if xs_type == 'total':

      # Get the Tally objects needed to compute the total xs
      flux = self.getTally('flux', filters)
      total = self.getTally('total', filters)

      # Initialize a MultiGroupXS object
      multigroup_xs = infermc.TotalXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['total'] = total

    elif xs_type == 'transport':

      # Get the Tally objects needed to compute the transport xs
      flux = self.getTally('flux', filters, estimator='analog')
      total = self.getTally('total', filters, estimator='analog')
      scatter1 = self.getTally('scatter-1', filters, estimator='analog')

      # Initialize a MultiGroupXS object
      multigroup_xs = infermc.TransportXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['total'] = total
      multigroup_xs._tallies['scatter-1'] = scatter1

    elif xs_type == 'absorption':

      # Get the Tally objects needed to compute the absorption xs
      flux = self.getTally('flux', filters)
      absorption = self.getTally('absorption', filters)

      # Initialize a MultiGroupXS object
      multigroup_xs = infermc.AbsorptionXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['absorption'] = absorption

    elif xs_type == 'capture':

      # Get the Tally objects needed to compute the absorption xs
      flux = self.getTally('flux', filters)
      absorption = self.getTally('absorption', filters)
      fission = self.getTally('fission', filters)

      # Initialize a MultiGroupXS object
      multigroup_xs = infermc.CaptureXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['absorption'] = absorption
      multigroup_xs._tallies['fission'] = fission

    elif xs_type == 'fission':

      # Get the Tally objects needed to compute the fission xs
      flux = self.getTally('flux', filters)
      fission = self.getTally('fission', filters)

      # Initialize a MultiGroupXS object
      multigroup_xs = infermc.FissionXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['fission'] = fission

    elif xs_type == 'nu-fission':

      # Get the Tally objects needed to compute the nu-fission xs
      flux = self.getTally('flux', filters)
      nu_fission = self.getTally('nu-fission', filters)

      # Initialize a MultiGroupXS object
      multigroup_xs = infermc.NuFissionXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['nu-fission'] = nu_fission

    elif xs_type == 'scatter':

      # Get the Tally objects needed to compute the scatter xs
      flux = self.getTally('flux', filters)
      scatter = self.getTally('scatter', filters)

      # Initialize a MultiGroupXS object
      multigroup_xs = infermc.ScatterXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['scatter'] = scatter

    elif xs_type == 'nu-scatter':

      # Get the Tally objects needed to compute the nu-scatter xs
      flux = self.getTally('flux', filters, estimator='analog')
      nu_scatter = self.getTally('nu-scatter', filters, estimator='analog')

      # Initialize a MultiGroupXS object
      multigroup_xs = infermc.NuScatterXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['nu-scatter'] = nu_scatter

    elif xs_type == 'scatter matrix':

      # Get the Tally objects needed to compute the scatter matrix
      flux = self.getTally('flux', filters, estimator='analog')
      scatter1 = self.getTally('scatter-1', filters, estimator='analog')

      filters.append(openmc.Filter(type='energyout', bins=group_edges))
      scatter = self.getTally('scatter', filters, estimator='analog')

      # Initialize a MultiGroupXS object
      multigroup_xs = infermc.ScatterMatrixXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['scatter-1'] = scatter1
      multigroup_xs._tallies['scatter'] = scatter

    elif xs_type == 'nu-scatter matrix':

      # Get the Tally objects needed to compute the scatter matrix
      flux = self.getTally('flux', filters, estimator='analog')
      nu_scatter1 = self.getTally('nu-scatter-1', filters, estimator='analog')

      filters.append(openmc.Filter(type='energyout', bins=group_edges))
      nu_scatter = self.getTally('nu-scatter', filters, estimator='analog')

      # Initialize a MultiGroupXS object
      multigroup_xs = infermc.NuScatterMatrixXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['nu-scatter-1'] = nu_scatter1
      multigroup_xs._tallies['nu-scatter'] = nu_scatter

    elif xs_type == 'chi':

      # Get the Tally objects needed to compute chi
      nu_fission_in = self.getTally('nu-fission', filters, estimator='analog')

      energyout_filter = openmc.Filter(type='energyout', bins=group_edges)
      filters.pop(0)
      filters.append(energyout_filter)
      nu_fission_out = self.getTally('nu-fission', filters, estimator='analog')

      # Initialize a MultiGroupXS object
      multigroup_xs = infermc.Chi(domain, domain_type, energy_groups)
      multigroup_xs._tallies['nu-fission-in'] = nu_fission_in
      multigroup_xs._tallies['nu-fission-out'] = nu_fission_out

    elif xs_type == 'diffusion':
      msg = 'Unable to get diffusion coefficient'
      raise ValueError(msg)

    # Compute the cross-section
    multigroup_xs.computeXS(corr)

    # Build offsets such that a user can query the MultiGroupXS for any region
    if domain_type == 'distribcell':

      multigroup_xs.findDomainOffset()
      domain_offset = multigroup_xs._offset

      # "Cache" a dictionary of region IDs to offsets
      num_regions = self._opencg_geometry._num_regions

      for region in range(num_regions):
        path = self.getPath(region)
        cell_id = path[-1]

        if cell_id == domain._id:
          offset = self._openmc_geometry.get_offset(path, domain_offset)
          multigroup_xs.setSubDomainOffset(region, offset)

    else:
      multigroup_xs.setSubDomainOffset(domain._id, 0)

    return multigroup_xs


  def getMultiGroupXS(self, xs_type, domain, domain_type):

    # Check that MultiGroupXS for the input parameters has been created

    if not domain_type in self._multigroup_xs.keys():
      msg = 'Unable to get cross-section since no cross-sections for ' \
            'domain type {0} have been created'.format(domain_type)
      raise ValueError(msg)

    # Check that the MultiGroupXS corresponding to this domain has been created
    elif not domain in self._multigroup_xs[domain_type].keys():
      msg = 'Unable to get cross-section since no cross-sections for ' \
            'domain {0} {1} have been created'.format(domain_type, domain)
      raise ValueError(msg)

    elif not xs_type in self._multigroup_xs[domain_type][domain].keys():
      msg = 'Unable to get cross-section since no cross-sections of type {0} ' \
            'for domain {1} {2} have been created'.format(xs_type,
                                                          domain_type, domain)
      raise ValueError(msg)

    return self._multigroup_xs[domain_type][domain][xs_type]


  def checkXS(self):

    for domain_type in self._multigroup_xs.keys():
      for domain_id in self._multigroup_xs[domain_type].keys():
        total_xs = self._multigroup_xs[domain_type][domain_id]['total']
        total_xs = total_xs.getXS()

        absorption_xs = self._multigroup_xs[domain_type][domain_id]['absorption']
        scatter_xs = self._multigroup_xs[domain_type][domain_id]['scatter']

        all_xs = absorption_xs.getXS()
        all_xs += scatter_xs.getXS()

        if not np.allclose(total_xs.ravel(), all_xs.ravel()):
          print('The nuclide micro xs {0} in {1} {2} is not equal to the '
                'total macro macro xs'.format(self._xs_type, self._domain_type,
                                              self._domain._id))


  def rebalanceAllScatterMatrices(self, domain_types='all'):

    #FIXME: This does not do error propagation for the uncertainties

    if domain_types == 'all':
      domain_types = self._multigroup_xs.keys()

    for domain_type in self._multigroup_xs.keys():

      if domain_type in domain_types:

        for domain_id in self._multigroup_xs[domain_type].keys():

          # Get dictionary of all MultiGroupXS for this
          all_domain_xs = self._multigroup_xs[domain_type][domain_id]

          if 'scatter matrix' in all_domain_xs.keys():

            # Get MultiGroupXS objects
            transport_xs = all_domain_xs['transport']
            absorption_xs = all_domain_xs['absorption']
            scatter_xs = all_domain_xs['scatter']
            scatter_matrix_xs = all_domain_xs['scatter matrix']
            nu_scatter_matrix_xs = all_domain_xs['nu-scatter matrix']

            # Get cross-section NumPy arrays
            transport = transport_xs.getXS()
            absorption = absorption_xs.getXS()
            scatter_matrix = scatter_matrix_xs.getXS()
            nu_scatter_matrix = nu_scatter_matrix_xs.getXS()

            # Compute scatter cross-section in each subdomain, group, nuclide
            scatter = transport - absorption

            # Compute rebalance factors for each subdomain, group, nuclide
            f = scatter / np.sum(scatter_matrix, axis=2)

            # Update scattering matrices with f factor
            scatter_matrix *= f
            nu_scatter_matrix *= f

            # Convert NaNs to zero
            scatter_matrix = np.nan_to_num(scatter_matrix)
            nu_scatter_matrix = np.nan_to_num(nu_scatter_matrix)

            # Assign rebalanced scattering matrixs to the MultiGroupXS
            scatter_matrix_xs._xs[0,...] = scatter_matrix
            nu_scatter_matrix_xs._xs[0,...] = nu_scatter_matrix


  def getMaxXS(self, xs_type, domain_type, group):

    max_xs = 1e10

    for domain in self._multigroup_xs[domain_type].keys():
      xs = self.getMultiGroupXS(xs_type, domain, domain_type)

      if not xs_type in ['scatter matrix', 'nu-scatter matrix']:
        data = xs.getXS(groups=[group])
      else:
        data = xs.getXS(in_groups=[group], out_groups=[group])

      max_xs = max(max_xs, np.amax(data.ravel()))

    return max_xs


  def getMinXS(self, xs_type, domain_type, group):

    min_xs = -1.

    for domain in self._multigroup_xs[domain_type].keys():
      xs = self.getMultiGroupXS(xs_type, domain, domain_type)

      if not xs_type in ['scatter matrix', 'nu-scatter matrix']:
        data = xs.getXS(groups=[group])
      else:
        data = xs.getXS(in_groups=[group], out_groups=[group])

      min_xs = min(min_xs, np.amin(data.ravel()))

    return min_xs


class MicroXSTallyExtractor(XSTallyExtractor):


  def extractAllMultiGroupXS(self, energy_groups, domain_type='distribcell',
                             nuclides='all', ignore_missing=True, corr=False):

    for xs_type in infermc.xs_types:
      self.extractMultiGroupXS(xs_type, energy_groups, domain_type,
                               nuclides, ignore_missing, corr)


  def extractMultiGroupXS(self, xs_type, energy_groups, domain_type='distribcell',
                          nuclides='all', ignore_missing=True, corr=False):

    # Add nested dictionary for this domain type if needed
    if not domain_type in self._multigroup_xs.keys():
      self._multigroup_xs[domain_type] = dict()

    # Get a list of the domains for this domain type to iterate over
    if domain_type == 'material':
      domains = self._openmc_geometry.get_all_materials()

    elif domain_type == 'universe':
      domains = self._openmc_geometry.get_all_material_universes()

    elif domain_type == 'cell' or domain_type == 'distribcell':
      domains = self._openmc_geometry.get_all_material_cells()

    # Iterate and create the MultiGroupXS for each domain
    for domain in domains:

      try:

        # Build the MultiGroupXS for this domain
        xs = self.createMultiGroupXS(xs_type, energy_groups, domain,
                                     domain_type, nuclides, corr)

        # Add nested dictionary for this domain if needed
        if not domain._id in self._multigroup_xs[domain_type]:
          self._multigroup_xs[domain_type][domain._id] = dict()

        # Store a handle to the MultiGroupXS object in the nested dictionary
        self._multigroup_xs[domain_type][domain._id][xs_type] = xs

      except (KeyError, ValueError):

        if ignore_missing:
          pass

        else:
          msg = 'Unable to build {0} xs for domain {1} {2} since ' \
                'the necessary tallies were not found in the ' \
                'statepoint'.format(xs_type, domain_type, domain._id)
          raise KeyError(msg)


  def createMultiGroupXS(self, xs_type, energy_groups, domain,
                         domain_type='distribcell', nuclides='all', corr=False):

    if self._statepoint is None:
      msg = 'Unable to get cross-sections since the TallyExtractor ' \
            'statepoint attribute has not been set'
      raise ValueError(msg)

    if self._summary is None:
      msg = 'Unable to get cross-sections since the TallyExtractor ' \
            'summary attribute has not been set'
      raise ValueError(msg)

    # Initialize a list of filters
    filters = list()

    # Create energy and domain filters to search for
    group_edges = energy_groups._group_edges
    filters.append(openmc.Filter(type='energy', bins=group_edges))
    filters.append(openmc.Filter(type=domain_type, bins=domain._id))

    if nuclides == 'all':
      all_nuclides = True

    # Extract a list of tuples of Nuclides and number densities (at/b-cm)
    # of all Nuclides in the domain of interest
    nuclides_densities = domain.get_all_nuclides().values()
    densities = list()
    nuclides = list()

    for nuclide_density in nuclides_densities:

      nuclide = nuclide_density[0]
      density = nuclide_density[1]

      if all_nuclides:
        nuclides.append(nuclide)
        densities.append(density)

      elif nuclide in nuclides:
        # Move nuclide to end of list ensure ordering of nuclides is
        # identical to the ordering of the nuclides in the domain
        nuclides.append(nuclides.pop(nuclides.index(nuclide)))
        densities.append(density)

    tot_density = 1.


    if xs_type == 'total':

      # Get the Tally objects needed to compute the total xs
      flux = self.getTally('flux', filters)
      total = self.getTally('total', filters, nuclides)

      # Initialize a MultiGroupXS object
      multigroup_xs = infermc.MicroTotalXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['total'] = total

    elif xs_type == 'transport':

      # Get the Tally objects needed to compute the transport xs
      flux = self.getTally('flux', filters, estimator='analog')
      total = self.getTally('total', filters, nuclides, estimator='analog')
      scatter1 = self.getTally('scatter-1', filters, nuclides, estimator='analog')

      # Initialize a MultiGroupXS object
      multigroup_xs = infermc.MicroTransportXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['total'] = total
      multigroup_xs._tallies['scatter-1'] = scatter1

    elif xs_type == 'absorption':

      # Get the Tally objects needed to compute the absorption xs
      flux = self.getTally('flux', filters)
      absorption = self.getTally('absorption', filters, nuclides)

      # Initialize a MultiGroupXS object
      multigroup_xs = infermc.MicroAbsorptionXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['absorption'] = absorption

    elif xs_type == 'capture':

      # Get the Tally objects needed to compute the absorption xs
      flux = self.getTally('flux', filters)
      absorption = self.getTally('absorption', filters, nuclides)
      fission = self.getTally('fission', filters, nuclides)

      # Initialize a MultiGroupXS object
      multigroup_xs = infermc.MicroCaptureXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['absorption'] = absorption
      multigroup_xs._tallies['fission'] = fission

    elif xs_type == 'fission':

      # Get the Tally objects needed to compute the fission xs
      flux = self.getTally('flux', filters)
      fission = self.getTally('fission', filters, nuclides)

      # Initialize a MultiGroupXS object
      multigroup_xs = infermc.MicroFissionXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['fission'] = fission

    elif xs_type == 'nu-fission':

      # Get the Tally objects needed to compute the nu-fission xs
      flux = self.getTally('flux', filters)
      nu_fission = self.getTally('nu-fission', filters, nuclides)

      # Initialize a MultiGroupXS object
      multigroup_xs = infermc.MicroNuFissionXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['nu-fission'] = nu_fission

    elif xs_type == 'scatter':

      # Get the Tally objects needed to compute the scatter xs
      flux = self.getTally('flux', filters)
      scatter = self.getTally('scatter', filters, nuclides)

      # Initialize a MultiGroupXS object
      multigroup_xs = infermc.MicroScatterXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['scatter'] = scatter

    elif xs_type == 'nu-scatter':

      # Get the Tally objects needed to compute the nu-scatter xs
      flux = self.getTally('flux', filters, estimator='analog')
      nu_scatter = self.getTally('nu-scatter', filters, nuclides, estimator='analog')

      # Initialize a MultiGroupXS object
      multigroup_xs = infermc.MicroNuScatterXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['nu-scatter'] = nu_scatter

    elif xs_type == 'scatter matrix':

      # Get the Tally objects needed to compute the scatter matrix
      flux = self.getTally('flux', filters, estimator='analog')
      scatter1 = self.getTally('scatter-1', filters, estimator='analog')

      filters.append(openmc.Filter(type='energyout', bins=group_edges))
      scatter = self.getTally('scatter', filters, estimator='analog')

      # Initialize a MultiGroupXS object
      multigroup_xs = infermc.MicroScatterMatrixXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['scatter-1'] = scatter1
      multigroup_xs._tallies['scatter'] = scatter

    elif xs_type == 'nu-scatter matrix':

      # Get the Tally objects needed to compute the scatter matrix
      flux = self.getTally('flux', filters, estimator='analog')
      nu_scatter1 = self.getTally('nu-scatter-1', filters, estimator='analog')

      filters.append(openmc.Filter(type='energyout', bins=group_edges))
      nu_scatter = self.getTally('nu-scatter', filters, estimator='analog')

      # Initialize a MultiGroupXS object
      multigroup_xs = infermc.MicroNuScatterMatrixXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['nu-scatter-1'] = nu_scatter1
      multigroup_xs._tallies['nu-scatter'] = nu_scatter

    elif xs_type == 'chi':

      # Get the Tally objects needed to compute chi
      nu_fission_in = self.getTally('nu-fission', filters, nuclides, estimator='analog')

      energyout_filter = openmc.Filter(type='energyout', bins=group_edges)
      filters.pop(0)
      filters.append(energyout_filter)
      nu_fission_out = self.getTally('nu-fission', filters, nuclides, estimator='analog')

      # Initialize a MultiGroupXS object
      multigroup_xs = infermc.MicroChi(domain, domain_type, energy_groups)
      multigroup_xs._tallies['nu-fission-in'] = nu_fission_in
      multigroup_xs._tallies['nu-fission-out'] = nu_fission_out

    elif xs_type == 'diffusion':
      msg = 'Unable to get diffusion coefficient'
      raise ValueError(msg)

    # FIXME!!!! - this does not work for simulations with a subset of nuclides
    # Add Nuclides and densities to the MicroXS
    multigroup_xs.addNuclides(nuclides, densities)
    if all_nuclides:
      multigroup_xs.addNuclide(openmc.Nuclide('total'), tot_density)

    # Compute the cross-section
    multigroup_xs.computeXS(corr)

    # Build offsets such that a user can query the MultiGroupXS for any region
    if domain_type == 'distribcell':

      multigroup_xs.findDomainOffset()
      domain_offset = multigroup_xs._offset

      # "Cache" a dictionary of region IDs to offsets
      num_regions = self._opencg_geometry._num_regions

      for region in range(num_regions):
        path = self.getPath(region)
        cell_id = path[-1]

        if cell_id == domain._id:
          offset = self._openmc_geometry.get_offset(path, domain_offset)
          multigroup_xs.setSubDomainOffset(region, offset)

    else:
      multigroup_xs.setSubDomainOffset(domain._id, 0)

    return multigroup_xs


  def checkXS(self):

    super(MicroXSTallyExtractor, self).checkXS()

    for domain_type in self._multigroup_xs.keys():
      for domain_id in self._multigroup_xs[domain_type].keys():
        for xs_type in self._multigroup_xs[domain_type][domain_id].keys():
          xs = self._multigroup_xs[domain_type][domain_id][xs_type]
          xs.checkXS()


  def getMaxXS(self, xs_type, nuclide, domain_type, group):

    max_xs = -1.

    for domain in self._multigroup_xs[domain_type].keys():
      xs = self.getMultiGroupXS(xs_type, domain, domain_type)

      if not xs.containsNuclide(nuclide):
        continue

      if not xs_type in ['scatter matrix', 'nu-scatter matrix']:
        data = xs.getXS(groups=[group], nuclides=[nuclide])
      else:
        data = xs.getXS(in_groups=[group], out_groups=[group], nuclides=[nuclide])

      max_xs = max(max_xs, np.amax(data.ravel()))

    return max_xs


  def getMinXS(self, xs_type, nuclide, domain_type, group):

    min_xs = 1e10

    for domain in self._multigroup_xs[domain_type].keys():
      xs = self.getMultiGroupXS(xs_type, domain, domain_type)

      if not xs.containsNuclide(nuclide):
        continue

      if not xs_type in ['scatter matrix', 'nu-scatter matrix']:
        data = xs.getXS(groups=[group], nuclides=[nuclide])
      else:
        data = xs.getXS(in_groups=[group], out_groups=[group], nuclides=[nuclide])

      min_xs = min(min_xs, np.amin(data.ravel()))

    return min_xs


  def getMaxRelErr(self, xs_type, nuclide, domain_type, group):

    max_rel_err = -1.

    for domain in self._multigroup_xs[domain_type].keys():
      xs = self.getMultiGroupXS(xs_type, domain, domain_type)

      if not xs.containsNuclide(nuclide):
        continue

      if not xs_type in ['scatter matrix', 'nu-scatter matrix']:
        data = xs.getRelErr(groups=[group], nuclides=[nuclide])
      else:
        data = xs.getRelErr(in_groups=[group], out_groups=[group], nuclides=[nuclide])

      max_rel_err = max(max_rel_err, np.amax(data.ravel()))

    return max_rel_err


  def getMinRelErr(self, xs_type, nuclide, domain_type, group):

    min_rel_err = 1e10

    for domain in self._multigroup_xs[domain_type].keys():
      xs = self.getMultiGroupXS(xs_type, domain, domain_type)

      if not xs.containsNuclide(nuclide):
        continue

      if not xs_type in ['scatter matrix', 'nu-scatter matrix']:
        data = xs.getRelErr(groups=[group], nuclides=[nuclide])
      else:
        data = xs.getRelErr(in_groups=[group], out_groups=[group], nuclides=[nuclide])

      min_rel_err = min(min_rel_err, np.amin(data.ravel()))

    return min_rel_err