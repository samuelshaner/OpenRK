import opencsg
from opencsg.checkvalue import *
from statepoint import StatePoint, Tally, score_types           # filter_types


xs_types = ['total',
            'transport',
            'absorption',
            'scatter',
            'nu-scatter',
            'scatter matrix',
            'fission',
            'nu-fission',
            'chi']

#            'diffusion'

domain_types = ['cell',
                'distribcell',
                'universe',
                'material']

#FIXME: mesh domain types


def get_path(coords):

  if not isinstance(coords, opencsg.LocalCoords):
    print('Unable to get the path with %s which is not an '
          'OpenCSG LocalCoords object' % str(coords))

  # Build "path" from LocalCoords
  path = list()

  while coords is not None:

    # If the LocalCoords is at a Universe
    if coords.getType() == 'universe':
      path.append(coords.getUniverse().getId())
      path.append(coords.getCell().getId())

    # If the LocalCoords is at a Lattice
    else:
      # Add 1 for Fortran indexing
      lat_x = coords.getLatticeX()+1
      lat_y = coords.getLatticeY()+1

      # Use z=1 for 3D lattices
      path.append((coords.getLattice().getId(), lat_x, lat_y, 1))

    # Traverse LocalCoords linked list to next lowest nested universe
    coords = coords.getNext()

  return path


class XSTallyExtractor(object):

  def __init__(self, statepoint=None, geometry=None):

    # Initialize TallyExtractor class attributes
    self._statepoint = None
    self._geometry = None

    # Dictionaries mapping cells/materials to Tallies
    # Keys   - Location ID
    # Values - Tally ID
    self._distribcells_to_tallies = dict()
    self._cells_to_tallies = dict()
    self._universes_to_tallies = dict()
    self._materials_to_tallies = dict()

    # Dictionary mapping Tallies to scores
    # Keys   - Tally ID
    # Values - score list
    self._tallies_to_scores = dict()

    # List of all "paths" to each region, indexed by region ID
    self._all_paths = list()

    if not statepoint is None:
      self.setStatePoint(statepoint)

    if not geometry is None:
      self.setGeometry(geometry)


  def getStatePoint(self):
    return self._statepoint


  def getGeometry(self):
    return self._geometry


  def getPath(self, region):

    if not is_integer(region):
      exit('Unable to get the path for region %s which is not an '
           'integer' % str(region))

    if region < 0:
      exit('Unable to get the path for region %d which is a negative '
           'integer' % region)

    if self._geometry is None:
      exit('Unable to get path for region %d since the TallyExtractors '
           'geometry attribute has not been set' % region)

    if region < len(self._geometry.getNumRegions()):
      exit('Unable to get path for region %d since it the Geometry only '
           'contains %d regions' % (region, self._geometry.getNumRegions()))

    return self._all_paths[region]


  def getAllPaths(self):
    return self._all_paths


  def getDistribcellsToTallies(self):
    return self._distribcells_to_tallies


  def getCellsToTallies(self):
    return self._cells_to_tallies


  def getUniversesToTallies(self):
    return self._universes_to_tallies


  def getMaterialsToTallies(self):
    return self._materials_to_tallies


  def setStatePoint(self, statepoint):

    if not isinstance(statepoint, StatePoint):
      exit('Unable to set the statepoint for a TallyExtractor to %s '
           'since it is not a StatePoint class object' % str(statepoint))

    self._statepoint = statepoint
    self._statepoint.read_results()

    # Create a mapping of Tally locations to location IDs to Tally IDs to scores
    for tally in self._statepoint.tallies:

      # Store a list of the Tally scores
      self._tallies_to_scores[tally.id] = tally.scores

      filters = tally.filters

      if 'distribcell' in filters.keys():
        filter = filters['distribcell']
        index = statepoint.geom.cellList.index(filter.bins[0])
        distribcell_id = statepoint.geom.cell[index].userID

        if not distribcell_id in self._distribcells_to_tallies.keys():
          self._distribcells_to_tallies[distribcell_id] = list()

        self._distribcells_to_tallies[distribcell_id].append(tally.id)

      if 'cell' in filters.keys():
        filter = filters['cell']
        index = statepoint.geom.cellList.index(filter.bins[0])
        cell_id = statepoint.geom.cell[index].userID

        if not cell_id in self._cells_to_tallies.keys():
          self._cells_to_tallies[cell_id] = list()

        self._cells_to_tallies[cell_id].append(tally.id)

      if 'universe' in filters.keys():
        filter = filters['universe']
        universe_id = filter.bins[0]

        if not universe_id in self.universes_to_tallies.keys():
          self._universes_to_tallies[universe_id] = list()

        self._universes_to_tallies[universe_id].append(tally.id)

      if 'material' in filters.keys():
        filter = filters['material']
        material_id = filter.bins[0]

        if not material_id in self._materials_to_tallies.keys():
          self._materials_to_tallies[material_id] = list()

        self._materials_to_tallies[material_id].append(tally.id)


  def setGeometry(self, geometry):

    if not isinstance(geometry, opencsg.Geometry):
      exit('Unable to set the geometry for a TallyExtractor to %s '
           'since it is not an OpenCSG Geometry class object' % str(geometry))

    self._geometry = geometry

    num_regions = self._geometry.getNumRegions()

    # Create a list of "paths" for each unique region in the Geometry
    self._geometry.initializeCellOffsets()

    # Compute the volumes for each region in the Geometry
    # First, compute the volume of the bounding box surrounding
    # the Geometry to use as a scaling factor for each subregion
    max_x = self._geometry.getMaxX()
    max_y = self._geometry.getMaxY()
    max_z = self._geometry.getMaxZ()
    min_x = self._geometry.getMinX()
    min_y = self._geometry.getMinY()
    min_z = self._geometry.getMinZ()

    delta_x = max_x - min_x
    delta_y = max_y - min_y
    delta_z = max_z - min_z

    volume = delta_x * delta_y * delta_z

#    self._geometry.setVolume(volume, tolerance=1e-1)

    for region in range(num_regions):
      coord = geometry.findRegion(region)
      self._all_paths.append(get_path(coord))


  def tallyContainsFilter(self, filter, tally):

    # Expect filter as a tuple: (filter type, bins)
    # Expect tally as an statepoint.Tally object

    contains_filter = False

    # Iterate over all of Tally's filters
    for filter_type in tally.filters.keys():

      test_filter = tally.filters[filter_type]

      # Check if the test filter is the same type of filter as
      # the query filter and contains the same bins
      if filter[0] == filter_type and filter[1] == test_filter.bins:
        contains_filter = True
        break

    return contains_filter


  def getTallyScoreIndex(self, score, tally):

    global score_types

    if not score in score_types.values():
      exit('Unable to get the index for score %s since it is an '
           'unsupported score type' % str(score))

    if not isinstance(tally, Tally):
      exit('Unable to get the index for score %s for Tally %s '
           'which is not a Tally object' % (str(score), str(tally)))

    return tally.scores.index(score)


  def getTally(self, score, filters, domain_id, domain='distribcell', estimator=None, label=''):

    global score_types, domain_types

    if self._statepoint is None:
      exit('Unable to get Tally for score %s in %s %s since '
           'the TallyExtractors statepoint attribute has not been set' %
           (str(score), str(domain), str(domain_id)))

    if not score in score_types.values():
      exit('Unable to get Tally for score %s in %s %s since the '
           'score type is not supported' %
           (str(score), str(domain), str(domain_id)))

    if not is_integer(domain_id):
      exit('Unable to get Tally for score %s in %s %s since the '
           'domain ID is not an integer value' %
           (str(score), str(domain), str(domain_id)))

    if domain_id < 0:
      exit('Unable to get Tally for score %s in %s %s since the '
           'domain ID is a negative integer' %
           (str(score), str(domain), str(domain_id)))

    if not domain in domain_types:
      exit('Unable to get Tally for score %s in %s %s since the '
           'domain type is not supported' %
           (str(score), str(domain), str(domain_id)))

    # Loop over the domain-to-tallies mapping to find the Tally
    tally = None
    tallies = None

    if domain == 'distribcell':
      tallies = self._distribcells_to_tallies
    elif domain == 'cell':
      tallies = self._cells_to_tallies
    elif domain == 'universe':
      tallies = self._universes_to_tallies
    elif domain == 'material':
      tallies = self._materials_to_tallies

    # Iterate over all tallies to find the appropriate one
    for tally_id in tallies[domain_id]:

      # Determine if the queried Tally score is associated with this Tally
      if score in self._tallies_to_scores[tally_id]:
        internal_id = self._statepoint.tallyID[tally_id]
        test_tally = self._statepoint.tallies[internal_id]

        #FIXME: Check tally estimator here!!!

        # If the length of the filters container doesn't match the requested
        # filters, continue to next Tally. Subtract one for the Tally domain
        # (e.g., material, cell, distribcell)
        if len(filters) != len(test_tally.filters) - 1:
          continue

        contains_filters = True

        # Iterate over the filters requested by the user
        for filter in filters:

          # If the test Tally does not contains this filter, break
          if not self.tallyContainsFilter(filter, test_tally):
            contains_filters = False
            break

        # If the Tally contained all of the filters, then we can return this Tally
        if contains_filters:
          tally = test_tally
          break

    # If we did not find the Tally, return an error messsage
    if tally is None:
      exit('Unable to get Tally for score %s in %s %d' %
           (score, domain, domain_id))

    return tally


  def getXS(self, xs_type, group_edges, domain_id, domain='distribcell'):

    global xs_types, domain_types

    if self._geometry is None:
      exit('Unable to get cross-sections since the TallyExtractors '
           'geometry attribute has not been set')

    if self._statepoint is None:
      exit('Unable to get cross-sections since the TallyExtractors '
           'statepoint attribute has not been set')

    if not xs_type in xs_types:
      exit('Unable to get %s cross-sections since it is not a '
           'valid cross-section type' % str(xs_type))

    if not isinstance(group_edges, (tuple, list, np.ndarray)):
      exit('Unable to get %s group cross-sections for %s %s since '
           'the group edges is not a Python tuple, list or NumPy array' %
           (str(group_edges), str(domain), str(domain_id)))

    # Determine the number of groups
    num_groups = group_edges.size-1

    if not is_integer(domain_id):
      exit('Unable to get %s group cross-sections for %s %s since the domain '
           'is not an integer' % (num_groups, str(domain), str(domain_id)))

    if domain_id < 0:
      exit('Unable to get %d group cross-sections for %s %d since the domain'
           'is a negative integer' % (num_groups, str(domain), domain_id))

    if not domain in domain_types:
      exit('Unable to get %d group cross-sections for %s %d since the domain '
           'type is not supported' % (num_groups, str(domain), domain_id))


    # Determine the filter for the tally domain
    if domain == 'distribcell':

      # Get the path to the region in the geometry
      coord = self._geometry.findRegion(domain_id)
      path = get_path(coord)

      # Use cell ID for the domain - necessary for the getTally(...) routine
      domain_id = path[-1]

      # Create distribcell filter for the StatePoint object using the path
      domain_filter = (domain, path)

    else:
      domain_filter = (domain, domain_id)


    #FIXME: Deal with diffusion coeff
    #FIXME: Deal with nuclides
    #FIXME: Deal with tracklength vs. analog tallies


    tallies = dict()
    data = dict()
    scores = dict()

    if xs_type == 'transport':

      #FIXME: Use analog tallies
      #FIXME: label='%d groups' % num_groups

      # Get the Tally IDs for the flux and reaction rate needed to compute the xs
      filters = [('energyin', list(group_edges))]
      tallies['flux'] = self.getTally('flux', filters, domain_id, domain)
      tallies['rxn-1'] = self.getTally('total', filters, domain_id, domain)
      tallies['rxn-2'] = self.getTally('scatter-pn', filters, domain_id, domain)

      # Initialize empty arrays for the flux and reaction rate data
      data['flux'] = np.zeros(num_groups)
      data['rxn-1'] = np.zeros(num_groups)
      data['rxn-2'] = np.zeros(num_groups)


    elif xs_type == 'scatter matrix':

      #FIXME: Use analog tallies
      #FIXME: label='%d groups' % num_groups

      # Get the Tally IDs for the flux and reaction rate needed to compute the xs
      filters = [('energyin', list(group_edges))]
      tallies['flux'] = self.getTally('flux', filters, domain_id, domain)

      filters.append(('energyout', list(group_edges)))
      tallies['rxn-1'] = self.getTally('nu-scatter', filters, domain_id, domain)

      # Initialize empty arrays for the flux and reaction rate data
      data['flux'] = np.zeros(num_groups)
      data['rxn-1'] = np.zeros((num_groups, num_groups))


    elif xs_type == 'chi':

      #FIXME: Use analog tallies
      #FIXME: label='%d groups' % num_groups

      # Get the Tally IDs for the flux and reaction rate needed to compute the xs
      filters = [('energyin', list(group_edges))]
      tallies['flux'] = self.getTally('flux', filters, domain_id, domain)
      tallies['rxn-1'] = self.getTally('nu-fission', filters, domain_id, domain)

      filters = [('energyout', list(group_edges))]
      tallies['rxn-2'] = self.getTally('nu-fission', filters, domain_id, domain)

      # Initialize empty arrays for the flux and reaction rate data
      data['flux'] = np.zeros(num_groups)
      data['rxn-1'] = np.zeros(num_groups)
      data['rxn-2'] = np.zeros(num_groups)


    elif xs_type == 'diffusion':
      exit('Unable to get diffusion coefficient')


    else:

      #FIXME: Do not use analog flux
      #FIXME: label='%d groups' % num_groups

      filters = [('energyin', list(group_edges))]

      # Get the Tallies for the flux and reaction rate needed to compute the xs
      tallies['flux'] = self.getTally('flux', filters, domain_id, domain)
      tallies['rxn-1'] = self.getTally(xs_type, filters, domain_id, domain)

      # Initialize empty arrays for the flux and reaction rate data
      data['flux'] = np.zeros(num_groups)
      data['rxn-1'] = np.zeros(num_groups)


    # Get the indices for each of the Tally scores
    for tally_name, tally in tallies.iteritems():
      score = tally.scores[0]                                     #FIXME: Is this safe?? What if a user defines several scores together in the same Tally
      scores[tally_name] = self.getTallyScoreIndex(score, tally)


    # Extract the flux and reaction rate tally averages for each energy group
    for in_group in range(num_groups):

      # Get the flux at this energy group
      flux = tallies['flux']
      filters = [domain_filter, ('energyin', in_group)]
      value = self._statepoint.get_value(flux.id, filters, scores['flux'])
      data['flux'][in_group] = value[0]


      if xs_type == 'transport':

        # Get the total reaction rate at this energy group
        rxn_rate = tallies['rxn-1']
        filters = [domain_filter, ('energyin', in_group)]
        value = self._statepoint.get_value(rxn_rate.id, filters, scores['rxn-1'])
        data['rxn-1'][in_group] = value[0]

        # Get the scatter-P1 reaction rate at this energy group
        rxn_rate = tallies['rxn-2']
        value = self._statepoint.get_value(rxn_rate.id, filters, scores['rxn-2'])
        data['rxn-2'][in_group] = value[0]


      elif xs_type == 'scatter matrix':

        # Need to loop over inner and outer energy groups
        for out_group in range(num_groups):

          filters = [domain_filter, ('energyin', in_group),
                     ('energyout', out_group)]
          rxn_rate = tallies['rxn-1']
          value = self._statepoint.get_value(rxn_rate.id, filters, scores['rxn-1'])
          data['rxn-1'][in_group, out_group] = value[0]


      if xs_type == 'chi':

        # Get the nufission reaction rate at this energy group
        rxn_rate = tallies['rxn-1']
        filters = [domain_filter, ('energyin', in_group)]
        value = self._statepoint.get_value(rxn_rate.id, filters, scores['rxn-1'])
        data['rxn-1'][in_group] = value[0]

        # Get the nufission reaction rate into this energy group
        rxn_rate = tallies['rxn-2']
        filters = [domain_filter, ('energyout', in_group)]
        value = self._statepoint.get_value(rxn_rate.id, filters, scores['rxn-2'])
        data['rxn-2'][in_group] = value[0]


      else:

        # Get the reaction rate at this energy group
        rxn_rate = tallies['rxn-1']
        filters = [domain_filter, ('energyin', in_group)]
        value = self._statepoint.get_value(rxn_rate.id, filters, scores['rxn-1'])
        data['rxn-1'][in_group] = value[0]


    # Compute the cross-section for each energy group and return it
    if xs_type == 'transport':

      xs = (data['rxn-1'] - data['rxn-2']) / data['flux']

      # Replace any negative values with 0.
      #FIXME: This may not be safe!!!
      xs = np.where(xs > 0., xs, 0.)


    elif xs_type == 'chi':

      norm = data['rxn-1'][:].sum()

      if norm == 0.:
        xs = np.zeros(num_groups)

      else:
        xs = data['rxn-2'] / data['rxn-1'][:].sum()

        # Normalize chi to 1.0
        xs /= xs.sum()


    elif xs_type == 'scatter matrix':
      xs = np.transpose(np.transpose(data['rxn-1']) / data['flux'])


    else:
      xs = data['rxn-1'] / data['flux']


    # For any region without flux, convert nan value to zero
    xs = np.nan_to_num(xs)

    return xs