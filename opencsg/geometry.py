__author__ = 'Will Boyd'
__email__ = 'wboyd@mit.edu'


from universe import *
from point import Point


class Geometry(object):

  def __init__(self):

    # Initialize Geometry class attributes
    self._universes = dict()
    self._lattices = dict()
    self._num_regions = 0
    self._volume = np.float64(0.)


  def getUniverses(self):
    return self._universes


  def getLattices(self):
    return self._lattices


  def getNumRegions(self):
    return self._num_regions


  def getVolume(self):
    return self._volume


  def getMaxX(self):

    if not 0 in self._universes.keys():
      exit('Unable to get the maximum x since the Geometry does not '
           'contain the base Universe ID=0')

    root = self._universes[0]
    return root.getMaxX()


  def getMaxY(self):

    if not 0 in self._universes.keys():
      exit('Unable to get the maximum y since the Geometry does not '
           'contain the base Universe ID=0')

    root = self._universes[0]
    return root.getMaxY()


  def getMaxZ(self):

    if not 0 in self._universes.keys():
      exit('Unable to get the maximum z since the Geometry does not '
           'contain the base Universe ID=0')

    root = self._universes[0]
    return root.getMaxZ()


  def getMinX(self):

    if not 0 in self._universes.keys():
      exit('Unable to get the minimum x since the Geometry does not '
           'contain the base Universe ID=0')

    root = self._universes[0]
    return root.getMinX()


  def getMinY(self):

    if not 0 in self._universes.keys():
      exit('Unable to get the minimum y since the Geometry does not '
           'contain the base Universe ID=0')

    root = self._universes[0]
    return root.getMinY()


  def getMinZ(self):

    if not 0 in self._universes.keys():
      exit('Unable to get the minimum z since the Geometry does not '
           'contain the base Universe ID=0')

    root = self._universes[0]
    return root.getMinZ()


  def addUniverse(self, universe):

    if not isinstance(universe, Universe):
      exit('Unable to add Universe to the Geometry since %s is not '
           'a Universe' % str(universe))

    univ_id = universe.getId()
    self._universes[univ_id] = universe


  def removeUniverse(self, universe):

    if not isinstance(universe, Universe):
      exit('Unable to remove Universe to the Geometry since %s is not '
           'a Universe' % str(universe))

    univ_id = universe.getId()
    if univ_id in self._universes.keys():
      del self._universes[univ_id]


  def addLattice(self, lattice):

    if not isinstance(lattice, Lattice):
      exit('Unable to add Lattice to the Geometry since %s is not '
           'a Lattice' % str(lattice))

    lat_id = lattice.getId()
    self._universes[lat_id] = lattice


  def removeLattice(self, lattice):

    if not isinstance(lattice, Lattice):
      exit('Unable to remove Lattice to the Geometry since %s is not '
           'a Lattice' % str(lattice))

    lat_id = lattice.getId()
    if lat_id in self._lattices.keys():
      del self._lattices[lat_id]


  def setVolume(self, volume, tolerance=1e-3):

    if not is_float(volume) and not is_integer(volume):
      exit('Unable to set the volume of the Geometry to %s since it is '
           'neither an integer nor a floating point value' % str(volume))

    if not is_float(tolerance):
      exit('Unable to compute the volume of the Geometry to a tolerance '
           'of %s since it is not a floating point value' % str(tolerance))

    if not 0 in self._universes.keys():
      exit('Unable to initialize cell offsets since the Geometry does not '
           'contain the base Universe ID=0')

    self._volume = np.float64(volume)

    root = self._universes[0]
    root.computeVolumeFractions(volume=self._volume, tolerance=tolerance)


#  def subdivideCells(self):

#    for universe_id in self._universes:
#      universe = self._universes[universe_id]
#      universe.subdivideCells()


  def initializeCellOffsets(self):

    if not 0 in self._universes.keys():
      exit('Unable to initialize cell offsets since the Geometry does not '
           'contain the base Universe ID=0')

    root = self._universes[0]
    root.initializeCellOffsets()
    self._num_regions = root.getNumRegions()


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

    if not 0 in self._universes.keys():
      exit('Unable to find cell since the Geometry does not contain the '
           'base Universe ID=0')

    root = self._universes[0]
    point = Point(x=x, y=y, z=z)
    localcoords = UnivCoords(point=point)
    localcoords.setUniverse(root)

    return root.findCell(localcoords=localcoords)




  #def findCell(self, region_id):




  def findRegion(self, region_id=0):

    if not 0 in self._universes.keys():
      exit('Unable to find coords since the Geometry does not contain the '
           'base Universe ID=0')

    if not is_integer(region_id):
      exit('Unable to get the path for region_id %s since it is '
           'not an integer value' % str(region_id))

    if region_id < 0:
      exit('Unable to get the path for region_id=%d since it is '
           'a negative integer' % region_id)


    root = self._universes[0]
    localcoords = UnivCoords(universe=root)
    root.findRegion(region_id=region_id, univ_coords=localcoords)
    localcoords = localcoords.getHeadNode()
    return localcoords


  def findCoords(self, x=0., y=0., z=0.):

    if not 0 in self._universes.keys():
      exit('Unable to find coords since the Geometry does not contain the '
           'base Universe ID=0')

    root = self._universes[0]
    point = Point(x=x, y=y, z=z)
    localcoords = UnivCoords(point=point)
    localcoords.setUniverse(root)
    cell = root.findCell(localcoords=localcoords)
    localcoords = localcoords.getHeadNode()

    return localcoords


  def toString(self):

    string = ''

    string += 'Geometry\n'

    for universe in self._universes.values():
      string += universe.toString()

    for lattice in self._lattices:
      string += lattice.toString()

    return string


  def printString(self):
    print(self.toString())