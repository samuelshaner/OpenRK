import opencsg
from openmc.opencsg_compatible import get_opencsg_geometry
from infermc.multigroupxs import *


def get_path(coords):

  if not isinstance(coords, opencsg.LocalCoords):
    msg = 'Unable to get the path with {0} which is not an ' \
          'OpenCSG LocalCoords object'.format(coords)
    raise ValueError(msg)

  # Build "path" from LocalCoords
  path = list()

  while coords is not None:

    # If the LocalCoords is at a Universe
    if coords._type == 'universe':
      path.append(coords._universe._id)
      path.append(coords._cell._id)

    # If the LocalCoords is at a Lattice
    else:
      # Add 1 for Fortran indexing
      lat_x = coords._lat_x+1
      lat_y = coords._lat_y+1

      # 3D Lattices only
      if coords._lat_z != None:
        lat_z = coords._lat_z+1

      # 2D Lattices
      else:
        lat_z = 1

      path.append((coords._lattice._id, lat_x, lat_y, lat_z))

    # Traverse LocalCoords linked list to next lowest nested universe
    coords = coords._next

  return path



class XSTallyExtractor(object):

  def __init__(self, statepoint=None):

    # Initialize TallyExtractor class attributes
    self._statepoint = None
    self._openmc_geometry = None
    self._opencsg_geometry = None

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
      self.setStatePoint(statepoint)


  def setStatePoint(self, statepoint):

    if not isinstance(statepoint, openmc.statepoint.StatePoint):
      msg = 'Unable to set the statepoint for a TallyExtractor to {0} ' \
            'since it is not a StatePoint class object'.format(statepoint)
      raise ValueError(msg)

    self._statepoint = statepoint
    self._statepoint.read_results()
    self._statepoint.compute_ci()

    # Retrieve the OpenMC Geometry from the statepoint and convert it
    # into an OpenCSG geometry object using the compatibility module
    self._openmc_geometry = self._statepoint._geometry
    self._opencsg_geometry = get_opencsg_geometry(self._openmc_geometry)

    # Build maps to optimize tally lookups
    self._buildTallyMaps()


  def _buildTallyMaps(self):

    self._buildDistribcellTallyMaps()
    self._buildCellTallyMaps()
    self._buildUniverseTallyMaps()
    self._buildMaterialTallyMaps()
    self._buildAllPaths()


  def _buildMaterialTallyMaps(self):

    # Create a mapping of Tally locations to location IDs to Tally IDs to scores
    for tally_id, tally in self._statepoint._tallies.iteritems():
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
    for tally_id, tally in self._statepoint._tallies.iteritems():
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
    for tally_id, tally in self._statepoint._tallies.iteritems():
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
    for tally_id, tally in self._statepoint._tallies.iteritems():

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


  def _buildAllPaths(self):

    # Create a list of "paths" for each unique region in the Geometry
    self._opencsg_geometry.initializeCellOffsets()

    num_regions = self._opencsg_geometry._num_regions

    for region in range(num_regions):
      coord = self._opencsg_geometry.findRegion(region)
      self._all_paths[region] = get_path(coord)


  def getPath(self, region):

    if not is_integer(region):
      msg = 'Unable to get the path for region {0} which is not an ' \
            'integer'.format(region)
      raise ValueError(msg)

    num_regions = self._opencsg_geometry._num_regions

    if region > num_regions:
      msg = 'Unable to get path for region {0} since it the Geometry only ' \
            'contains {1} regions'.format(region, num_regions)
      raise ValueError(msg)

    elif region < 0:
      msg = 'Unable to get the path for region {0} which is a negative ' \
            'integer'.format(region)
      raise ValueError(msg)

    return self._all_paths[region]


  def getTally(self, score, filters, estimator='tracklength', label=''):

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
      for filter in filters:

        # If the test Tally does not contains this filter, break
        if not filter in test_tally._filters:
          contains_filters = False
          break

      # If the Tally contained all of the filters, then return this Tally
      if contains_filters:
        tally = test_tally
        break

    # If we did not find the Tally, return an error message
    if tally is None:
      msg = 'Unable to get Tally for score {0}'.format(score)
      raise ValueError(msg)

    return tally


  def extractAllMultiGroupXS(self, energy_groups, domain_type='distribcell'):

    for xs_type in xs_types:
      self.extractMultiGroupXS(xs_type, energy_groups, domain_type)


  def extractMultiGroupXS(self, xs_type, energy_groups, domain_type='distribcell'):

    # Add nested dictionary for this domain type if needed
    if not domain_type in self._multigroup_xs.keys():
      self._multigroup_xs[domain_type] = dict()

    # Get a list of the domains for this domain type to iterate over
    if domain_type == 'material':
      domains = self._openmc_geometry.getAllMaterials()

    elif domain_type == 'universe':
      domains = self._openmc_geometry.getAllMaterialUniverses()

    elif domain_type == 'cell' or domain_type == 'distribcell':
      domains = self._openmc_geometry.getAllMaterialCells()

    # Iterate and create the MultiGroupXS for each domain
    for domain in domains:

      # Add nested dictionary for this domain if needed
      if not domain._id in self._multigroup_xs[domain_type].keys():
        self._multigroup_xs[domain_type][domain._id] = dict()

      # Build the MultiGroupXS for this domain
      xs = self.createMultiGroupXS(xs_type, energy_groups, domain, domain_type)

      # Store a handle to the MultiGroupXS object in the nested dictionary
      self._multigroup_xs[domain_type][domain._id][xs_type] = xs



  def createMultiGroupXS(self, xs_type, energy_groups,
                         domain, domain_type='distribcell'):

    if self._statepoint is None:
      msg = 'Unable to get cross-sections since the TallyExtractor ' \
            'statepoint attribute has not been set'
      raise ValueError(msg)

    elif not xs_type in xs_types:
      msg = 'Unable to get {0} cross-sections since it is not a ' \
            'valid cross-section type'.format(xs_type)
      raise ValueError(msg)

    elif not isinstance(energy_groups, EnergyGroups):
      msg = 'Unable to get group cross-sections with energy groups ' \
            '{0} since it is not an EnergyGroups object'.format(energy_groups)
      raise ValueError(msg)

    elif not isinstance(domain, (openmc.Material, openmc.Cell,
                                 openmc.Universe, openmc.Mesh)):
      msg = 'Unable to get group cross-sections for domain ' \
            '{0} since it is not a valid domain'.format(domain)
      raise ValueError(msg)

    elif not domain_type in domain_types:
      msg = 'Unable to get group cross-sections for domain type ' \
            '{0} since it is not a supported domain type'.format(domain_type)
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
      multigroup_xs = TotalXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['total'] = total


    elif xs_type == 'transport':

      # Get the Tally objects needed to compute the transport xs
      flux = self.getTally('flux', filters, estimator='analog')
      total = self.getTally('total', filters, estimator='analog')
      scatter1 = self.getTally('scatter-1', filters, estimator='analog')

      # Initialize a MultiGroupXS object
      multigroup_xs = TransportXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['total'] = total
      multigroup_xs._tallies['scatter-1'] = scatter1


    elif xs_type == 'absorption':

      # Get the Tally objects needed to compute the absorption xs
      flux = self.getTally('flux', filters)
      absorption = self.getTally('absorption', filters)

      # Initialize a MultiGroupXS object
      multigroup_xs = AbsorptionXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['absorption'] = absorption


    elif xs_type == 'fission':

      # Get the Tally objects needed to compute the fission xs
      flux = self.getTally('flux', filters)
      fission = self.getTally('fission', filters)

      # Initialize a MultiGroupXS object
      multigroup_xs = FissionXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['fission'] = fission


    elif xs_type == 'nu-fission':

      # Get the Tally objects needed to compute the nu-fission xs
      flux = self.getTally('flux', filters)
      nu_fission = self.getTally('nu-fission', filters)

      # Initialize a MultiGroupXS object
      multigroup_xs = NuFissionXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['nu-fission'] = nu_fission


    elif xs_type == 'scatter':

      # Get the Tally objects needed to compute the scatter xs
      flux = self.getTally('flux', filters)
      scatter = self.getTally('scatter', filters)

      # Initialize a MultiGroupXS object
      multigroup_xs = ScatterXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['scatter'] = scatter


    elif xs_type == 'nu-scatter':

      # Get the Tally objects needed to compute the nu-scatter xs
      flux = self.getTally('flux', filters, estimator='analog')
      nu_scatter = self.getTally('nu-scatter', filters, estimator='analog')

      # Initialize a MultiGroupXS object
      multigroup_xs = NuScatterXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['nu-scatter'] = nu_scatter


    elif xs_type == 'scatter matrix':

      # Get the Tally objects needed to compute the scatter matrix
      flux = self.getTally('flux', filters, estimator='analog')

      filters.append(openmc.Filter(type='energyout', bins=group_edges))
      nu_scatter = self.getTally('nu-scatter', filters, estimator='analog')

      # Initialize a MultiGroupXS object
      multigroup_xs = ScatterMatrixXS(domain, domain_type, energy_groups)
      multigroup_xs._tallies['flux'] = flux
      multigroup_xs._tallies['nu-scatter'] = nu_scatter


    elif xs_type == 'chi':

      # Get the Tally objects needed to compute chi
      nu_fission_in = self.getTally('nu-fission', filters, estimator='analog')

      energyout_filter = openmc.Filter(type='energyout', bins=group_edges)
      filters.pop(0)
      filters.append(energyout_filter)
      nu_fission_out = self.getTally('nu-fission', filters, estimator='analog')

      # Initialize a MultiGroupXS object
      multigroup_xs = Chi(domain, domain_type, energy_groups)
      multigroup_xs._tallies['nu-fission-in'] = nu_fission_in
      multigroup_xs._tallies['nu-fission-out'] = nu_fission_out


    elif xs_type == 'diffusion':
      msg = 'Unable to get diffusion coefficient'
      raise ValueError(msg)


    # Compute the cross-section
    multigroup_xs.computeXS()



    # Build offsets such that a user can query the MultiGroupXS for any region
    if domain_type == 'distribcell':

      multigroup_xs.findDomainOffset()
      domain_offset = multigroup_xs._offset

      # "Cache" a dictionary of region IDs to offsets
      num_regions = self._opencsg_geometry._num_regions

      for region in range(num_regions):
        path = self.getPath(region)
        cell_id = path[-1]

        if cell_id == domain._id:
          offset = self._openmc_geometry.getOffset(path, domain_offset)
          multigroup_xs.setSubDomainOffset(region, offset)


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
