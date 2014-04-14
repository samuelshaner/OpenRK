__author__ = 'Will Boyd'
__email__ = 'wboyd@mit.edu'


from universe import *
from point import Point


class Geometry(object):

  def __init__(self, root_universe=None):

    # Initialize Geometry class attributes
#    self._universes = dict()
#    self._lattices = dict()
    self._root_universe = None
    self._num_regions = 0
    self._volume = np.float64(0.)


  def getRootUniverse(self):
    return self._root_universe


  def getNumRegions(self):
    return self._num_regions


  def getVolume(self):
    return self._volume


  def getMaxX(self):

    if self._root_universe is None:
      exit('Unable to get the maximum x since the Geometry does not '
           'contain the base Universe ID=0')

    return self._root_universe.getMaxX()


  def getMaxY(self):

    if self._root_universe is None:
      exit('Unable to get the maximum y since the Geometry does not '
           'contain the base Universe ID=0')

    return self._root_universe.getMaxY()


  def getMaxZ(self):

    if self._root_universe is None:
      exit('Unable to get the maximum z since the Geometry does not '
           'contain the base Universe ID=0')

    return self._root_universe.getMaxZ()


  def getMinX(self):

    if self._root_universe is None:
      exit('Unable to get the minimum x since the Geometry does not '
           'contain the base Universe ID=0')

    return self._root_universe.getMinX()


  def getMinY(self):

    if self._root_universe is None:
      exit('Unable to get the minimum y since the Geometry does not '
           'contain the base Universe ID=0')

    return self._root_universe.getMinY()


  def getMinZ(self):

    if self._root_universe is None:
      exit('Unable to get the minimum z since the Geometry does not '
           'contain the base Universe ID=0')

    return self._root_universe.getMinZ()


  def setRootUniverse(self, root_universe):

    if not isinstance(root_universe, Universe):
      exit('Unable to set the root Universe for the Geometry since %s is '
           'not a Universe' % str(root_universe))

    if not root_universe.getId() == 0:
      exit('Unable to set the root Universe for the Geometry with a Universe '
           'with ID=%d. The root Universe must have ID=0.' %
           root_universe.getId())

    self._root_universe = root_universe


  def setVolume(self, volume, tolerance=1e-3):

    if not is_float(volume) and not is_integer(volume):
      exit('Unable to set the volume of the Geometry to %s since it is '
           'neither an integer nor a floating point value' % str(volume))

    if not is_float(tolerance):
      exit('Unable to compute the volume of the Geometry to a tolerance '
           'of %s since it is not a floating point value' % str(tolerance))

    if self._root_universe is None:
      exit('Unable to initialize cell offsets since the Geometry does not '
           'contain the base Universe ID=0')

    self._volume = np.float64(volume)

    self._root_universe.computeVolumeFractions(volume=self._volume,
                                               tolerance=tolerance)


  def initializeCellOffsets(self):

    if self._root_universe is None:
      exit('Unable to initialize cell offsets since the Geometry does not '
           'contain the base Universe ID=0')

    self._root_universe.initializeCellOffsets()
    self._num_regions = self._root_universe.getNumRegions()


  def getRegionId(self, x=0., y=0., z=0.):

    coords = self.findCoords(x=x, y=y, z=z)
    region_id = 0

    while not coords is None:

      # The coords is a UnivCoords object
      if isinstance(coords, UnivCoords):
        universe = coords.getUniverse()
        cell = coords.getCell()
        region_id += universe.getCellOffset(cell)

      # The coords is a LatCoords object
      else:
        lattice = coords.getLattice()
        lat_x, lat_y = coords.getLatticeX(), coords.getLatticeY()
        region_id += lattice.getCellOffset(lat_x, lat_y)

      coords = coords.getNext()

    return region_id


  def findCell(self, x=0., y=0., z=0.):

    if self._root_universe is None:
      exit('Unable to find cell since the Geometry does not contain the '
           'base Universe ID=0')

    point = Point(x=x, y=y, z=z)
    localcoords = UnivCoords(point=point)
    localcoords.setUniverse(self._root_universe)

    return self._root_universe.findCell(localcoords=localcoords)




  #def findCell(self, region_id):




  def findRegion(self, region_id=0):

    if self._root_universe is None:
      exit('Unable to find coords since the Geometry does not contain the '
           'base Universe ID=0')

    if not is_integer(region_id):
      exit('Unable to get the path for region_id %s since it is '
           'not an integer value' % str(region_id))

    if region_id < 0:
      exit('Unable to get the path for region_id=%d since it is '
           'a negative integer' % region_id)


    localcoords = UnivCoords(universe=self._root_universe)
    self._root_universe.findRegion(region_id=region_id, univ_coords=localcoords)
    localcoords = localcoords.getHeadNode()
    return localcoords


  def findCoords(self, x=0., y=0., z=0.):

    if self._root_universe is None:
      exit('Unable to find coords since the Geometry does not contain the '
           'base Universe ID=0')

    point = Point(x=x, y=y, z=z)
    localcoords = UnivCoords(point=point)
    localcoords.setUniverse(self._root_universe)
    cell = self._root_universe.findCell(localcoords=localcoords)
    localcoords = localcoords.getHeadNode()

    return localcoords