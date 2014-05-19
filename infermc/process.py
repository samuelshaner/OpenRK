import opencsg
from opencsg.checkvalue import *
from statepoint import StatePoint, Tally, score_types           # filter_types

xs_types = ['total',
            'absorption',
            'scatter',
            'scatter matrix',
            'fission',
            'nu-fission',
            'chi']
#            'diffusion',
#            'transport']

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


  def getTallyScoreIndex(self, score, tally):

    global score_types

    if not score in score_types.values():
      exit('Unable to get the index for score %s since it is an '
           'unsupported score type' % str(score))

    if not isinstance(tally, Tally):
      exit('Unable to get the index for score %s for Tally %s '
           'which is not a Tally object' % (str(score), str(tally)))

    return tally.scores.index(score)


  def getTally(self, score, domain_id, domain='distribcell', label=''):

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
    if domain == 'cell':
      tallies = self._cells_to_tallies
    if domain == 'universe':
      tallies = self._universes_to_tallies
    if domain == 'material':
      tallies = self._materials_to_tallies

    # If we reached this point then we did not find the Tally
    if tallies is None:
      exit('Unable to get Tally in %s %d' %
           (score, domain, domain_id))

    for tally_id in tallies[domain_id]:
      if score in self._tallies_to_scores[tally_id]:
        tally = self._statepoint.tallies[self._statepoint.tallyID[tally_id]]

    # If we reached this point then we did not find the Tally
    if tally is None:
      exit('Unable to get Tally for score %s in %s %d' %
           (score, domain, domain_id))

    return tally


  def getXS(self, xs_type, num_groups, domain_id, domain='distribcell'):

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

    if not is_integer(num_groups):
      exit('Unable to get %s group cross-sections for %s %s since '
           'the number of groups is not an integer value' %
           (str(num_groups), str(domain), str(domain_id)))

    if num_groups < 0:
      exit('Unable to get %s group cross-sections for %s %s since '
           'the number of groups is a negative integer' %
           (num_groups, str(domain), str(domain_id)))

    if not is_integer(domain_id):
      exit('Unable to get %s group cross-sections for %s %s since the domain '
           'is not an integer' % (num_groups, str(domain), str(domain_id)))

    if domain_id < 0:
      exit('Unable to get %d group cross-sections for %s %d since the domain'
           'is a negative integer' % (num_groups, str(domain), domain_id))

    if not domain in domain_types:
      exit('Unable to get %d group cross-sections for %s %d since the domain '
           'type is not supported' % (num_groups, str(domain), domain_id))


    #FIXME: Must deal with diffusion coeff, transport xs, and chi here
    #FIXME: Must deal with nuclides here


    # Determine the reaction rate score type for this cross-section type
    rxn_rate_score = xs_type

    if xs_type == 'scatter matrix':
      rxn_rate_score = 'nu-scatter'
    elif xs_type == 'chi':
      rxn_rate_score = 'nu-fission'
    elif xs_type == 'diffusion':
      exit('Unable to get diffusion coefficient')
    elif xs_type == 'transport':
      exit('Unable to get transport xs')


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


    # Get the Tally IDs for the flux and reaction rate needed to compute the xs
    flux_tally = self.getTally('flux', domain_id, domain)                       #FIXME: label='%d groups' % num_groups
    rxn_rate_tally = self.getTally(rxn_rate_score, domain_id, domain)           #FIXME: label='%d groups' % num_groups

    # Initialize empty arrays for the flux and reaction rate data
    flux_data = np.zeros(num_groups)
    rxn_rate_data = np.zeros(num_groups)                                        #FIXME: Doesn't work for chi, scatter matrix

    # Compute the index for flux and reaction rate Tally scores
    flux_index = self.getTallyScoreIndex('flux', flux_tally)
    rxn_rate_index = self.getTallyScoreIndex(rxn_rate_score, rxn_rate_tally)


    # Extract the flux and reaction rate tally averages for each energy group
    for group in range(num_groups):

      filters = [domain_filter, ('energyin', group)]

      # Get the flux at this energy group
      data = self._statepoint.get_value(flux_tally.id, filters, flux_index)
      flux_data[group] = data[0]

      # Get the reaction rate at this energy group
      data = self._statepoint.get_value(rxn_rate_tally.id, filters, rxn_rate_index)
      rxn_rate_data[group] = data[0]


    # Compute the cross-section for this location and return it
    xs = rxn_rate_data / flux_data

    return xs