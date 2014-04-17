import opencsg
from opencsg.checkvalue import *
from statepoint import StatePoint


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


class TallyExtractor(object):

  def __init__(self, statepoint=None, geometry=None):

    # Initialize TallyExtractor class attributes
    self._statepoint = None
    self._geometry = None

    # Dictionary mapping cells to Tallies
    # Keys   - Cell ID
    # Values - Tally ID
    self._cells_to_tallies = dict()

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


  def getCellsToTallies(self):
    return self._cells_to_tallies


  def setStatePoint(self, statepoint):

    if not isinstance(statepoint, StatePoint):
      exit('Unable to set the statepoint for a TallyExtractor to %s '
           'since it is not a StatePoint class object' % str(statepoint))

    self._statepoint = statepoint
    self._statepoint.read_results()

    # Create a mapping of Cell IDs to distribcell Tally IDs
    for tally in self._statepoint.tallies:

      filters = tally.filters

      if 'distribcell' in filters.keys():
        filter = filters['distribcell']
        index = statepoint.geom.cellList.index(filter.bins[0])
        cell_id = statepoint.geom.cell[index].userID
        self._cells_to_tallies[cell_id] = tally.id


  def setGeometry(self, geometry):

    if not isinstance(geometry, opencsg.Geometry):
      exit('Unable to set the geometry for a TallyExtractor to %s '
           'since it is not an OpenCSG Geometry class object' % str(geometry))

    self._geometry = geometry

    # Create a list of "paths" for each unique region in the Geometry
    self._geometry.initializeCellOffsets()

    num_regions = self._geometry.getNumRegions()
    print('The Geometry has %d regions' % num_regions)

    for region in range(num_regions):
      coord = geometry.findRegion(region)
      self._all_paths.append(get_path(coord))


  def getDistribcellTallyData(self, region, score):

    if self._statepoint is None:
      exit('Unable to get tally data for region %s and score %s since '
           'the TallyExtractors statepoint attribute has not been set' %
           (str(region), str(score)))

    if self._geometry is None:
      exit('Unable to get tally data for region %s and score %s since '
           'the TallyExtractors geometry attribute has not been set' %
           (str(region), str(score)))

    coord = self._geometry.findRegion(region)

    path = get_path(coord)
    filters = [('distribcell', path)]

    cell_id = path[-1]
    tally_id = self._cells_to_tallies[cell_id]

    if score is 'flux':
      tally_data = self._statepoint.get_value(tally_id, filters, 0)
    elif score is 'absorption':
      tally_data = self._statepoint.get_value(tally_id, filters, 1)

    return tally_data


