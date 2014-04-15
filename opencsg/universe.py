__author__ = 'Will Boyd'
__email__ = 'wboyd@mit.edu'


from material import Material
from surface import Surface, on_surface_thresh
from point import Point
from checkvalue import *
import numpy as np
import math


# Error threshold for determining how close to the boundary of a Lattice cell
# a Point needs to be to be considered on it
on_lattice_cell_thresh = 1e-12

# Lists of all IDs for all Universes created
universe_ids = list()

# A static variable for auto-generated Universe IDs
auto_universe_id = 10000


class Universe(object):

  def __init__(self, universe_id=None, name=''):

    # Initialize Cell class attributes
    self._id = None
    self._name = None

    # Keys   - Cell IDs
    # Values - Cells
    self._cells = dict()

    # Keys   - Cell IDs
    # Values - Offsets
    self._cell_offsets = dict()
    self._num_regions = 0
    self._volume = np.float("inf")

    self.setId(universe_id)
    self.setName(name)

  def getId(self):
    return self._id


  def getName(self):
    return self._name


  def getCells(self):
    return self._cells


  def getCellOffset(self, cell):

    if not isinstance(cell, Cell):
      exit('Unable to get cell offset from Universe ID=%d for %s which '
           'is not a Cell' % (self._id, str(cell)))

    cell_id = cell.getId()

    if not cell_id in self._cell_offsets.keys():
      exit('Unable to get cell offset from Universe ID=%d for Cell ID=%d '
           'which is not in the Universe' % (self._id, cell_id))

    return self._cell_offsets[cell_id]


  def getCellOffsets(self):
    return self._cell_offsets


  def getNumRegions(self):
    return self._num_regions


  def getVolume(self):
    return self._volume


  def getMaxX(self):

    max_x = -np.float64("inf")

    for cell_id in self._cells:

      cell = self._cells[cell_id]
      if max_x < cell.getMaxX():
        max_x = cell.getMaxX()

    return max_x


  def getMaxY(self):

    max_y = -np.float64("inf")

    for cell_id in self._cells:

      cell = self._cells[cell_id]
      if max_y < cell.getMaxY():
        max_y = cell.getMaxY()

    return max_y


  def getMaxZ(self):

    max_z = -np.float64("inf")

    for cell_id in self._cells:

      cell = self._cells[cell_id]
      if max_z < cell.getMaxZ():
        max_z = cell.getMaxZ()

    return max_z


  def getMinX(self):

    min_x = np.float64("inf")

    for cell_id in self._cells:

      cell = self._cells[cell_id]
      if min_x > cell.getMinX():
        min_x = cell.getMinX()

    return min_x


  def getMinY(self):

    min_y = np.float64("inf")

    for cell_id in self._cells:

      cell = self._cells[cell_id]
      if min_y > cell.getMinY():
        min_y = cell.getMinY()

    return min_y


  def getMinZ(self):

    min_z = np.float64("inf")

    for cell_id in self._cells:

      cell = self._cells[cell_id]
      if min_z > cell.getMinZ():
        min_z = cell.getMinZ()

    return min_z


  def setMaxX(self, max_x):

    if not is_float(max_x) and not is_integer(max_x):
      exit('Unable to set the maximum x for Universe ID=%d to %s since it '
           'is not an integer or floating point value' % (self._id, str(max_x)))

    for cell_id in self._cells.keys():
      cell = self._cells[cell_id]
      cell.setMaxX(max_x=max_x)


  def setMaxY(self, max_y):

    if not is_float(max_y) and not is_integer(max_y):
      exit('Unable to set the maximum y for Universe ID=%d to %s since it '
           'is not an integer or floating point value' % (self._id, str(max_y)))

    for cell_id in self._cells.keys():
      cell = self._cells[cell_id]
      cell.setMaxY(max_y=max_y)


  def setMaxZ(self, max_z):

    if not is_float(max_z) and not is_integer(max_z):
      exit('Unable to set the maximum z for Universe ID=%d to %s since it '
           'is not an integer or floating point value' % (self._id, str(max_z)))

    for cell_id in self._cells.keys():
      cell = self._cells[cell_id]
      cell.setMaxZ(max_z=max_z)


  def setMinX(self, min_x):

    if not is_float(min_x) and not is_integer(min_x):
      exit('Unable to set the minimum x for Universe ID=%d to %s since it '
           'is not an integer or floating point value' % (self._id, str(min_x)))

    for cell_id in self._cells.keys():
      cell = self._cells[cell_id]
      cell.setMinX(min_x=min_x)


  def setMinY(self, min_y):

    if not is_float(min_y) and not is_integer(min_y):
      exit('Unable to set the minimum y for Universe ID=%d to %s since it '
           'is not an integer or floating point value' % (self._id, str(min_y)))

    for cell_id in self._cells.keys():
      cell = self._cells[cell_id]
      cell.setMinY(min_y=min_y)


  def setMinZ(self, min_z):

    if not is_float(min_z) and not is_integer(min_z):
      exit('Unable to set the minimum z for Universe ID=%d to %s since it '
           'is not an integer or floating point value' % (self._id, str(min_z)))

    for cell_id in self._cells.keys():
      cell = self._cells[cell_id]
      cell.setMinZ(min_z=min_z)


  def setId(self, universe_id=None):

    global universe_ids

    if universe_id is None:
      global auto_universe_id
      self._id = auto_universe_id
      universe_ids.append(auto_universe_id)
      auto_universe_id += 1

    # Check that the ID is an integer and wasn't already used
    elif is_integer(universe_id):

      # If the Cell already has an ID, remove it from global list
      if not self._id is None:
        universe_ids.remove(self._id)

      if universe_id in universe_ids:
        exit('Unable to set Universe ID to %s since a Universe with this ID '
             'was already initialized' % str(universe_id))

      if universe_id < 0:
        exit('Unable to set Univeres ID to %d since it must be a '
             'non-negative integer' % universe_id)

      else:
        self._id = universe_id
        universe_ids.append(universe_id)

    else:
      exit('Unable to set Universe ID to a non-integer %s' % str(universe_id))


  def setName(self, name):

    if not is_string(name):
      exit('Unable to set name for Universe ID=%d with a non-string value %s',
           self._id, str(name))

    else:
      self._name = name


  def addCell(self, cell):

    if not isinstance(cell, Cell):
      exit('Unable to add a Cell to Universe ID=%d since %s is not a Cell',
            self._id, str(cell))

    cell_id = cell.getId()

    if not cell_id in self._cells.keys():
      self._cells[cell_id] = cell


  def addCells(self, cells):

    if not isinstance(cells, (list, tuple, np.ndarray)):
      exit('Unable to add Cells to Universe ID=%d since %s is not a Python '
           'tuple/list or NumPy array' % (self._id, str(cells)))

    for i in range(len(cells)):
      self.addCell(cells[i])


  def removeCell(self, cell):

    if not isinstance(cell, Cell):
      exit('Unable to remove a Cell from Universe ID=%d since %s is not a '
           'Cell' % (self._id, str(cell)))

    cell_id = cell.getId()

    if cell_id in self._cells.keys():
      del self._cells[cell_id]


  def clearCells(self):
    self._cells.clear()


  def computeVolumeFractions(self, volume=np.float64(1.), tolerance=1e-3):

    if not is_float(volume):
      exit('Unable to compute volume fractions for Universe ID=%d since '
           'volume=%s is not a floating point value' % (self._id, str(volume)))

    self._volume = volume

    for cell_id in self._cells:

      cell = self._cells[cell_id]
      cell.computeVolumeFraction(volume=self._volume, tolerance=tolerance)

      # Recursively compute volume fractions below this Universe
      fill = cell.getFill()

      if isinstance(fill, Universe):
        volume_fraction = cell.getVolumeFraction()
        fill.computeVolumeFractions(volume=volume*volume_fraction,
                                    tolerance=tolerance)


  def initializeCellOffsets(self):


    # If we have already called this routine, return the total number of regions
    if not self._num_regions == 0:
      return self._num_regions

    # The cell offsets have not yet been initialized
    count = 0

    for cell_id in self._cells.keys():
      self._cell_offsets[cell_id] = count
      cell = self._cells[cell_id]
      count += cell.getNumSubCells()

    self._num_regions = count


  def findCell(self, localcoords):

    if not isinstance(localcoords, LocalCoords):
      exit('Unable to find cell in Universe ID=%d since localcoords %s is not '
           'a LocalCoords' % (self._id, str(localcoords)))

    for cell_id in self._cells:

      cell = self._cells[cell_id]

      if cell.containsPoint(localcoords.getPoint()):
        localcoords.setCell(cell)

        # 'material' type Cell - lowest level, terminate search for Cell
        if cell.getType() == 'material':
            return cell

        # 'fill' type Cell - Cell contains a Universe at a lower level
        # Update coords to next level and continue search
        else:

          fill = cell.getFill()

          if isinstance(fill, Lattice):
            next_coords = LatCoords(localcoords.getPoint())
            next_coords.setLattice(fill)

          elif isinstance(fill, Universe):
            next_coords = UnivCoords(localcoords.getPoint())
            next_coords.setUniverse(fill)

          else:
            exit('Unable to find cell since in Universe ID=%d since Cell ID=%d '
                 'is not filled by a Material, Universe or Lattice' %
                 (self._id, cell.getId()))

          localcoords.setNext(next_coords)
          next_coords.setPrev(localcoords)

          return fill.findCell(next_coords)


#  def findCell(self, region_id):


  def findRegion(self, region_id, univ_coords):

    if not is_integer(region_id):
      exit('Unable to find region_id in Universe ID=%d since %s is '
           'not an integer value' % (self._id, str(region_id)))

    if not isinstance(univ_coords, UnivCoords):
      exit('Unable to find region_id in Universe ID=%d since %s is '
           'not a UnivCoords' % (self._id, str(univ_coords)))

    # Initialize cell and offset
    cell = None
    offset = 0

    # Loop over cells until we reach the one the first one with
    # an offset larger than region_id - return the one prior
    for cell_id in self._cells.keys():

      cell = self._cells[cell_id]

      if self._cell_offsets[cell_id] <= region_id:
        offset = self._cell_offsets[cell_id]

      elif self._cell_offsets[cell_id] > region_id:
        break

    if cell is None:
      exit('Unable to find region_id=%d in Universe '
           'ID=%d' % (self._id, region_id))

    region_id -= offset
    univ_coords.setCell(cell)
    fill = cell.getFill()

    # The next nested level is at a Lattice
    if region_id >= 0 and isinstance(fill, Lattice):
      next_coords = LatCoords(lattice=fill)
      univ_coords.setNext(next_coords)
      next_coords.setPrev(univ_coords)
      fill.findRegion(region_id, next_coords)

    # The next nested level is at a Universe
    elif region_id >= 0 and isinstance(fill, Universe):
      next_coords = UnivCoords(universe=fill)
      univ_coords.setNext(next_coords)
      next_coords.setPrev(univ_coords)
      fill.findRegion(region_id, next_coords)

    # We have found the cell in the base nested universe
    elif region_id == 0 and isinstance(fill, Material):
      return

    # We were unable to find the cell
    else:
      exit('Unable to find cell for region_id=%d in Universe '
           'ID=%d' % (self._id, region_id))


  def toString(self):

    string = ''

    string += 'Universe\n'

    universe_id = '{0: <16}'.format('\tID') + '=\t' + str(self._id)
    string += universe_id + '\n'

    name = '{0: <16}'.format('\tName') + '=\t' + self._name
    string += name + '\n'

    num_regions = '{0: <16}'.format('\t# Regions') + '=\t'
    num_regions += str(self._num_regions)
    string += num_regions + '\n'

    cells = '{0: <16}'.format('\tCells') + '=\t'
    cells += str(self._cells.keys())
    string += cells + '\n'

    return string


  def printString(self):
    print(self.toString())



################################################################################
###################################  Lattice  ##################################
################################################################################


class Lattice(Universe):

  def __init__(self, lattice_id=None, name='', type='rectangular'):

    # Initialize Lattice class attributes
    self._id = None
    self._name = None
    self._type = ''
    self._dimension = None
    self._lower_left = None
    self._width = None
    self._universes = None
    self._cell_offsets = None
    self._num_regions = None

    self.setId(lattice_id)
    self.setName(name)
    self.setType(type)


  def getId(self):
    return self._id


  def getName(self):
    return self._name


  def getType(self):
    return self._type


  def getDimension(self):
    return self._dimension


  def getLowerLeft(self):
    return self._lower_left


  def getWidth(self):
    return self._width


  def getUniverse(self, lat_x, lat_y):

    if not is_integer(lat_x):
      exit('Unable to get Universe from Lattice ID=%d since x=%s is not '
           'a lattice cell' % self._id, str(lat_x))

    if not is_integer(lat_y):
      exit('Unable to get Universe from Lattice ID=%d since y=%s is not '
           'a lattice cell' % self._id, str(lat_y))

    if lat_x < 0 or lat_x > self._dimension[0]:
      exit('Unable to get Universe from Lattice ID=%d since x=%s is '
           'outside the bounds of the lattice cells', self._id, lat_x)

    if lat_y < 0 or lat_y > self._dimension[1]:
      exit('Unable to get Universe from Lattice ID=%d since y=%s is '
           'outside the bounds of the lattice cells', self._id, lat_y)

    return self._universes[lat_x][lat_y]


  def getUniverses(self):
    return self._universes


  def getCellOffset(self, lat_x, lat_y):

    if not is_integer(lat_x):
      exit('Unable to get cell offset from Lattice ID=%d since lat_x=%s is '
           'not a lattice cell' % self._id, str(lat_x))

    if not is_integer(lat_y):
      exit('Unable to get cell offset from Lattice ID=%d since lat_y=%s is '
           'not a lattice cell' % self._id, str(lat_y))

    if lat_x < 0 or lat_x > self._dimension[0]:
      exit('Unable to get Universe from Lattice ID=%d since lat_x=%s is '
           'outside the bounds of the lattice cells', self._id, lat_x)

    if lat_y < 0 or lat_y > self._dimension[1]:
      exit('Unable to get Universe from Lattice ID=%d since lat_y=%s is '
           'outside the bounds of the lattice cells', self._id, lat_y)

    return self._cell_offsets[lat_x][lat_y]


  def getMaxX(self):
    return self._lower_left[0] + self._dimension[0] * self._width[0]


  def getMaxY(self):
    return self._lower_left[1] + self._dimension[1] * self._width[1]


  def getMinX(self):
    return self._lower_left[0]


  def getMinY(self):
    return self._lower_left[1]


  def setId(self, lattice_id=None):

    global universe_ids

    if lattice_id is None:
      global auto_universe_id
      self._id = auto_universe_id
      universe_ids.append(auto_universe_id)
      auto_universe_id += 1

    # Check that the ID is an integer and wasn't already used
    elif is_integer(lattice_id):

      # If the Lattice already has an ID, remove it from global list
      if not self._id is None:
        universe_ids.remove(self._id)

      if lattice_id in universe_ids:
        exit('Unable to set Lattice ID to %s since a Lattice '
             'with this ID was already initialized.' % str(lattice_id))

      if lattice_id < 0:
        exit('Unable to set Lattice ID to %d since it must be a '
             'non-negative integer' % lattice_id)

      else:
        self._id = lattice_id
        universe_ids.append(lattice_id)

    else:
      exit('Unable to set a non-integer Lattice ID %s' % str(lattice_id))


  def setName(self, name):

    if not is_string(name):
      exit('Unable to set name for Lattice ID=%d with a non-string value %s' %
           (self._id, str(name)))

    else:
      self._name = name


  def setType(self, type):

    if not is_string(type):
      exit('Unable to set the type for Lattice ID=%d with a non-string '
           'value %s' % (self._id, str(type)))

    if not type in ['rectangular']:
      exit('Unable to set the type for Lattice ID=%d to %s since it is not '
           'rectangular' % (self._id, type))

    self._type = type


  def setDimension(self, dimension):

    if not isinstance(dimension, (tuple, list, np.ndarray)):
      exit('Unable to set Lattice ID=%d dimension to %s since it is not '
           'a Python tuple/list or NumPy array' % (self._id, str(dimension)))

    if len(dimension) != 2 and len(dimension) != 3:
      exit('Unable to set Lattice ID=%d dimension to %s since it does '
           'not contain 2 or 3 coordinates' % (self._id, str(dimension)))

    for dim in dimension:

      if not isinstance(dim, (int, np.int32, np.int64)):
        exit('Unable to set the dimension for Lattice ID=%d to %s since it '
             'is not an integer' % (self._id, str(dim)))

      if dim < 0:
        exit('Unable to set Lattice ID=%d dimension to %s since it '
             'is a negative value' % (self._id, dim))

    self._dimension = np.zeros(len(dimension), dtype=np.int64)

    for i in range(len(dimension)):
      self._dimension[i] = dimension[i]


  def setLowerLeft(self, lower_left):

    if not isinstance(lower_left, (tuple, list, np.ndarray)):
      exit('Unable to set the lower_left for Lattice ID=%d to %s since it is '
           'not a Python tuple/list or NumPy array' % (self._id, str(lower_left)))

    if len(lower_left) != 2 and len(lower_left) != 3:
      exit('Unable to set the lower_left for Lattice ID=%d to %s since it does '
           'not contain 2 or 3 coordinates' % (self._id, str(lower_left)))

    for dim in lower_left:

      if not isinstance(dim, (int, np.int32, np.int64)) \
        and not isinstance(dim, (float, np.float32, np.float64)):
        exit('Unable to set the lower_left for Lattice ID=%d to %s since it '
             'is not an integer or floating point value' % (self._id, str(dim)))

    self._lower_left = np.zeros(len(lower_left), dtype=np.float64)

    for i in range(len(lower_left)):
      self._lower_left[i] = lower_left[i]


  def setWidth(self, width):

    if not isinstance(width, (tuple, list, np.ndarray)):
      exit('Unable to set the width for Lattice ID=%d to %s since it is not '
           'a Python tuple/list or NumPy array' % (self._id, str(width)))

    if len(width) != 2 and len(width) != 3:
      exit('Unable to set the width for Lattice ID=%d to %s since it does '
           'not contain 2 or 3 coordinates' % (self._id, str(width)))

    for dim in width:

      if not isinstance(dim, (int, np.int32, np.int64)) \
        and not isinstance(dim, (float, np.float32, np.float64)):
        exit('Unable to set the width for Lattice ID=%d to %s since it is '
             'not an integer or floating point value' % (self._id, str(dim)))

      if dim < 0:
        exit('Unable to set the width for Lattice ID=%d to %s since it '
             'is a negative value' % (self._id, dim))

    self._width = width

    self._width = np.zeros(len(width), dtype=np.float64)

    for i in range(len(width)):
      self._width[i] = width[i]


  def setUniverses(self, universes):

    if not isinstance(universes, (tuple, list, np.ndarray)):
      exit('Unable to set the universes for Lattice ID=%d to %s since it is '
           'not a Python tuple/list or NumPy array' % (self._id, str(universes)))

    self._universes = universes

    for i in range(len(self._universes)):
      for j in range(len(self._universes[i])):

        universe = self._universes[i][j]

        universe.setMaxX(self._width[0]/2.)
        universe.setMinX(-self._width[0]/2.)
        universe.setMaxY(self._width[1]/2.)
        universe.setMinY(-self._width[1]/2.)

        if len(self._width) == 3:
          universe.setMaxZ(self._width[2]/2.)
          universe.setMinZ(-self._width[2]/2.)


  def computeVolumeFractions(self, volume=np.float64(1.), tolerance=1e-3):

    if not is_float(volume):
      exit('Unable to compute volume fractions for Lattice ID=%d since '
           'volume=%s is not a floating point value' % (self._id, str(volume)))

    volume_fraction = np.float64(1. / (self._dimension[0] * self._dimension[1]))

    for i in range(len(self._universes)):
      for j in range(len(self._universes[i])):
        universe = self._universes[i][j]
        universe.computeVolumeFractions(volume=(volume * volume_fraction),
                                        tolerance=tolerance)


  def initializeCellOffsets(self):

    # If we have already called this routine, return the total number of regions
    if not self._num_regions is None:
      return self._num_regions

    # First, determine how many dimensions is in the lattice
    num_dims = len(self._dimension)

    # Initialize an array for the cell offsets
    self._cell_offsets = np.zeros(tuple(self._dimension), dtype=np.int64)

    # The cell offsets have not yet been initialized
    count = 0

    if num_dims == 2:

      for i in range(self._dimension[0]):
        for j in range(self._dimension[1]):
          self._cell_offsets[i][j] = count
          self._universes[i][j].initializeCellOffsets()
          count += self._universes[i][j].getNumRegions()

    if num_dims == 3:

      for i in range(self._dimension[0]):
        for j in range(self._dimension[1]):
          for k in range(self._dimension[2]):
            self._cell_offsets[i,j,k] = count
            self._universes[i,j,k].initializeCellOffsets()
            count += self._universes[i,j,k].getNumRegions()


    self._num_regions = count


  def findCell(self, localcoords):

    if not isinstance(localcoords, LocalCoords):
      exit('Unable to find cell in Lattice ID=%d since localcoords %s is not '
           'a LocalCoords' % (self._id, str(localcoords)))

    # Compute the x and y indices for the Lattice cell this coord is in
    point = localcoords.getPoint()
    x = point.getX()
    y = point.getY()

    # Compute the Lattice cell indices
    lat_x = math.floor((x - self._lower_left[0]) / self._width[0])
    lat_y = math.floor((y - self._lower_left[1]) / self._width[1])

    # Check if the LocalCoord is on the Lattice boundaries
    # If so adjust x or y Lattice cell indices

    # Compute the distance to the Lattice cell boundaries
    distance_x = math.fabs(math.fabs(x) - self._dimension[0]*self._width[0]*0.5)
    distance_y = math.fabs(math.fabs(y) - self._dimension[1]*self._width[1]*0.5)

    if distance_x < on_lattice_cell_thresh:
      if x > 0:
        lat_x = self._dimension[0] - 1
      else:
        lat_x = 0

    if distance_y < on_lattice_cell_thresh:
      if y > 0:
        lat_y = self._dimension[1] - 1
      else:
        lat_y = 0

    # If the indices are outside the bound of the Lattice
    if (lat_x < 0 or lat_x >= self._dimension[0]) or \
      (lat_y < 0 or lat_y >= self._dimension[0]):
      exit('Unable to find cell since the lattice indices (%d,%d) are '
           'outside of Lattice ID=%d' % (lat_x, lat_y, self._id))

    # Cast the Lattice indices as integers
    lat_x = int(lat_x)
    lat_y = int(lat_y)

    # Set the Lattice cell indices for the LocalCoords
    localcoords.setLatticeX(lat_x)
    localcoords.setLatticeY(lat_y)

    # Compute local position of Point in the next level Universe
    next_x = x - (self._lower_left[0] + (lat_x + 0.5) * self._width[0])
    next_y = y - (self._lower_left[1] + (lat_y + 0.5) * self._width[1])
    next_point = Point(x=next_x, y=next_y)

    # Get the Universe for this Lattice cell
    universe = self._universes[lat_x][lat_y]

    if isinstance(universe, Universe):
      next_coords = UnivCoords(next_point)
      next_coords.setUniverse(universe)

    elif isinstance(universe, Lattice):
      next_coords = LatCoords(next_point)
      next_coords.setLattice(universe)

    else:
      exit('Unable to find Cell since Lattice ID=%d does not contain a '
           'Universe or Lattice in lattice cell (%d, %d)' %
           (self._id, lat_x, lat_y))

    localcoords.setNext(next_coords)
    next_coords.setPrev(localcoords)

    return universe.findCell(next_coords)


  def findRegion(self, region_id, lat_coords):

    if not is_integer(region_id):
      exit('Unable to find region_id=%d in Lattice ID=%d since %s is '
           'not an integer value' % (region_id, self._id, str(region_id)))

    if not isinstance(lat_coords, LatCoords):
      exit('Unable to find region_id=%d in Lattice ID=%d since %s is '
           'not a LatCoords' % (region_id, self._id, str(lat_coords)))

    # Initialize cell and offset
    universe = None
    offset = 0
    lat_x = 0
    lat_y = 0

    # Loop over cells until we reach the one the first one with
    # an offset larger than region_id - return the one prior
    for i in range(len(self._cell_offsets)):
      for j in range(len(self._cell_offsets[i])):

        universe = self._universes[i][j]

        if self._cell_offsets[i][j] <= region_id:
          offset = self._cell_offsets[i][j]
          lat_x = i
          lat_y = j

        elif self._cell_offsets[i][j] > region_id:
          break


    if universe is None:
      exit('Unable to find region_id=%d for FSR ID=%d in Lattice '
           'ID=%d' % (region_id, self._id))

    region_id -= offset
    lat_coords.setLatticeX(lat_x)
    lat_coords.setLatticeY(lat_y)

    # The next nested level is at a Lattice
    if region_id >= 0 and isinstance(universe, Lattice):
      next_coords = LatCoords(lattice=universe)
      lat_coords.setNext(next_coords)
      next_coords.setPrev(lat_coords)
      universe.findRegion(region_id, next_coords)

    # The next nested level is at a Universe
    elif region_id >= 0 and isinstance(universe, Universe):
      next_coords = UnivCoords(universe=universe)
      lat_coords.setNext(next_coords)
      next_coords.setPrev(lat_coords)
      universe.findRegion(region_id, next_coords)

    # We were unable to find the Cell
    else:
      exit('Unable to find region_id=%d in Lattice '
           'ID=%d' % (region_id, self._id))


  def toString(self):

    string = ''

    string += 'Lattice\n'

    cell_id = '{0: <16}'.format('\tID') + '=\t' + str(self._id)
    string += cell_id + '\n'

    name = '{0: <16}'.format('\tName') + '=\t' + self._name
    string += name + '\n'

    type = '{0: <16}'.format('\tType') + '=\t'
    type += str(self._type)
    string += type + '\n'

    dimension = '{0: <16}'.format('\tDimension') + '=\t'
    dimension += str(self._dimension)
    string += dimension + '\n'

    lower_left = '{0: <16}'.format('\tLower Left') + '=\t'
    lower_left += str(self._lower_left)
    string += lower_left + '\n'

    width = '{0: <16}'.format('\tWidth') + '=\t'
    width += str(self._width)
    string += width + '\n'

    num_regions = '{0: <16}'.format('\t# Regions') + '=\t'
    num_regions += str(self._num_regions)
    string += num_regions + '\n'

    universes = '{0: <16}'.format('\tUniverses') + '\n'

    for i in range(len(self._universes)):
      universes += '\t'

      for j in range(len(self._universes[0])):
        universes += '%s ' % str(self._universes[i][j].getId())

      universes += '\n'

    string += universes.rstrip('\n')

    return string


  def printString(self):
    print(self.toString())



################################################################################
####################################  Cell  ####################################
################################################################################


# Lists of all IDs for all Cells created
cell_ids = list()

# A static variable for auto-generated Cell IDs
auto_cell_id = 10000


class Cell(object):

  def __init__(self, cell_id=None, name='', fill=None):

    # Initialize Cell class attributes
    self._id = None
    self._name = None
    self._fill = None
    self._type = None
    self._num_subcells = None
    self._volume_fraction = np.float64(0.)
    self._volume = np.float64(0.)

    # Keys   - Surface IDs
    # Values - (halfpsace, Surface) tuples
    self._surfaces = dict()

    # Max/min values
    self._max_x = None
    self._max_y = None
    self._max_z = None
    self._min_x = None
    self._min_y = None
    self._min_z = None

    self.setId(cell_id)
    self.setName(name)
    self.setMaxX(np.float64("inf"))
    self.setMaxY(np.float64("inf"))
    self.setMaxZ(np.float64("inf"))
    self.setMinX(-np.float64("inf"))
    self.setMinY(-np.float64("inf"))
    self.setMinZ(-np.float64("inf"))

    if not fill is None:
      self.setFill(fill)


  def getId(self):
    return self._id


  def getName(self):
    return self._name


  def getFill(self):
    return self._fill


  def getType(self):
    return self._type


  def getVolumeFraction(self):
    return self._volume_fraction


  def getVolume(self):
    return self._volume


  def getSurfaces(self):
    return self._surfaces


  def getMaxX(self):
    return self._max_x


  def getMaxY(self):
    return self._max_y


  def getMaxZ(self):
    return self._max_z


  def getMinX(self):
    return self._min_x


  def getMinY(self):
    return self._min_y


  def getMinZ(self):
    return self._min_z


  def setId(self, cell_id=None):

    global cell_ids

    if cell_id is None:
      global auto_cell_id
      self._id = auto_cell_id
      cell_ids.append(auto_cell_id)
      auto_cell_id += 1

    # Check that the ID is an integer and wasn't already used
    elif is_integer(cell_id):

      # If the Cell already has an ID, remove it from global list
      if not self._id is None:
        cell_ids.remove(self._id)

      if cell_id in cell_ids:
        exit('Unable to set Cell ID to %s since a Cell with this ID was '
             'already initialized.' % str(cell_id))

      if cell_id < 0:
        exit('Unable to set Cell ID to %d since it must be a '
             'non-negative integer', cell_id)

      else:
        self._id = cell_id
        cell_ids.append(cell_id)

    else:
      exit('Unable to set Cell ID to a non-integer %s' % str(cell_id))


  def setName(self, name):

    if not is_string(name):
      exit('Unable to set name for Cell ID=%d with a non-string value %s' %
           (self._id, str(name)))

    else:
      self._name = name


  def setFill(self, fill):

    if isinstance(fill, Lattice):
      self.setType('lattice')
    elif isinstance(fill, Universe):
      self.setType('universe')
    elif isinstance(fill, Material):
      self.setType('material')
    else:
      exit('Unable to set fill for Cell ID=%d to %s since it is not a '
           'Universe, Lattice or a Material' % (self._id, str(fill)))

    self._fill = fill


  def setType(self, type):

    if not is_string(type):
      exit('Unable to set the type for Cell ID=%d to %s since it is not '
           'a string' % (self._id, str(type)))

    if not type.lower() in ['universe', 'lattice', 'material']:
      exit('Unable to set the type for Cell ID=%d to %s since it is not '
           'universe, lattice or material' % (self._id, type))

    self._type = type.lower()


  def setMaxX(self, max_x):

    if not is_float(max_x):
      exit('Unable to set the maximum x-coordinate for Cell ID=%d to %s '
           'since it is not a floating point value', self._id, str(max_x))

    self._max_x = max_x


  def setMaxY(self, max_y):

    if not is_float(max_y):
      exit('Unable to set the maximum y-coordinate for Cell ID=%d to %s '
           'since it is not a floating point value', self._id, str(max_y))

    self._max_y = max_y


  def setMaxZ(self, max_z):

    if not is_float(max_z):
      exit('Unable to set the maximum z-coordinate for Cell ID=%d to %s '
           'since it is not a floating point value', self._id, str(max_z))

    self._max_z = max_z


  def setMinX(self, min_x):

    if not is_float(min_x):
      exit('Unable to set the minimum x-coordinate for Cell ID=%d to %s '
           'since it is not a floating point value', self._id, str(min_x))

    self._min_x = min_x


  def setMinY(self, min_y):

    if not is_float(min_y):
      exit('Unable to set the minimum y-coordinate for Cell ID=%d to %s '
           'since it is not a floating point value', self._id, str(min_y))

    self._min_y = min_y


  def setMinZ(self, min_z):

    if not is_float(min_z):
      exit('Unable to set the minimum z-coordinate for Cell ID=%d to %s '
           'since it is not a floating point value', self._id, str(min_z))

    self._min_z = min_z


  def addSurface(self, surface, halfspace):

    if not isinstance(surface, Surface):
      exit('Unable to add a Surface to Cell ID=%d since %s is not a Surface' %
           (self._id, str(surface)))

    if not halfspace in [-1, +1]:
      exit('Unable to add a Surface to Cell ID=%d with halfspace %s since '
           'it is not +/-1' % (self._id, str(halfspace)))

    # Only add the Surface if the Cell does not already contain another
    # Surface of the same type and coefficients (a redundant Surface)
    for surface_id in self._surfaces.keys():

      test_surface = self._surfaces[surface_id][0]
      test_halfspace = self._surfaces[surface_id][1]

      if halfspace == test_halfspace:

        if surface.getType() == test_surface.getType():
          coeffs = surface.getCoeffs()
          test_coeffs = test_surface.getCoeffs()

          for key in coeffs.keys():
            coeff = coeffs[key]
            test_coeff = test_coeffs[key]

          if abs(coeff-test_coeff) < 1e-10:
            return

    surface_id = surface.getId()

    if not surface_id in self._surfaces.keys():
      self._surfaces[surface_id] = (surface, halfspace)

    self.findBoundingBox()


  def addSurfaces(self, surfaces, halfspaces):

    if not isinstance(surfaces, (list, tuple, np.ndarray)):
      exit('Unable to add Surfaces to Cell ID=%d since %s is not a Python '
           'tuple/list or NumPy array' % (self._id, str(surfaces)))

    if not isinstance(halfspaces, (list, tuple, np.ndarray)):
      exit('Unable to add Surfaces to Cell ID=%d since %s is not a Python '
           'tuple/list or NumPy array' % (self._id, str(halfspaces)))

    if len(surfaces) != len(halfspaces):
      exit('Unable to add Surfaces to Cell ID=%d since the number of '
           'Surfaces (%d) and halfspaces (%d) are not equal' %
           (self._id, len(surfaces), len(halfspaces)))

    for i in range(len(surfaces)):
      self.addSurface(surfaces[i], halfspaces[i])


  def removeSurface(self, surface):

    if not isinstance(surface, Surface):
      exit('Unable to remove a surface from Cell ID=%d since %s is not a '
           'Surface' % (self._id, str(surface)))

    surf_id = surface.getId()

    if surf_id in self._surfaces.keys():
      del self._surfaces[surf_id]

    self.findBoundingBox()


  def removeRedundantSurfaces(self):

    for surface_id in self._surfaces.keys():
      surface = self._surfaces[surface_id][0]
      halfspace = self._surfaces[surface_id][1]

      max_x = surface.getMaxX(halfspace=halfspace)
      max_y = surface.getMaxY(halfspace=halfspace)
      max_z = surface.getMaxZ(halfspace=halfspace)

      min_x = surface.getMinX(halfspace=halfspace)
      min_y = surface.getMinY(halfspace=halfspace)
      min_z = surface.getMinZ(halfspace=halfspace)

      delta_box = np.zeros(6)

      if abs(max_x) == np.float("inf") and abs(self._max_x) == np.float("inf"):
        delta_box[0] = np.float("nan")
      else:
        delta_box[0] = max_x - self._max_x

      if abs(min_x) == np.float("inf") and abs(self._min_x) == np.float("inf"):
        delta_box[1] = np.float("nan")
      else:
        delta_box[1] = self._min_x - min_x

      if abs(max_y) == np.float("inf") and abs(self._max_y) == np.float("inf"):
        delta_box[2] = np.float("nan")
      else:
        delta_box[2] = max_y - self._max_y

      if abs(min_y) == np.float("inf") and abs(self._min_y) == np.float("inf"):
        delta_box[3] = np.float("nan")
      else:
        delta_box[3] = self._min_y - min_y

      if abs(max_z) == np.float("inf") and abs(self._max_z) == np.float("inf"):
        delta_box[4] = np.float("nan")
      else:
        delta_box[4] = max_z - self._max_z

      if abs(min_z) == np.float("inf") and abs(self._min_z) == np.float("inf"):
        delta_box[5] = np.float("nan")
      else:
        delta_box[5] = self._min_z - min_z

      to_remove = True

      for delta in delta_box:
        if delta in [np.float("inf"), 0.]:
          to_remove = False


      if to_remove:
        del self._surfaces[surface_id]


  def setVolumeFraction(self, volume_fraction):

    if not is_float(volume_fraction):
      exit('Unable to set the volume fraction for Cell ID=%s to %s since it '
           'is not a floating point value' % str(volume_fraction))

    self._volume_fraction = volume_fraction


  def setVolume(self, volume):

    if not is_float(volume):
      exit('Unable to set the volume for Cell ID=%s to %s since it '
           'is not a floating point value' % str(volume))

    self._volume = volume


  def getNumSubCells(self):

    # If we have already called this routine, return the number of subcells
    if not self._num_subcells is None:
      return self._num_subcells

    # The cell offsets have not yet been initialized - we must compute them
    elif isinstance(self._fill, Material):
      self._num_subcells = 1

    elif isinstance(self._fill, (Universe, Lattice)):
      self._fill.initializeCellOffsets()
      self._num_subcells = self._fill.getNumRegions()

    else:
      exit('Unable to compute the number of subcells for Cell ID=%d since '
           'it is not filled by a Material, Universe or Lattice' % (self._id))

    return self._num_subcells


  def findBoundingBox(self):

    self.setMaxX(np.float64("inf"))
    self.setMaxY(np.float64("inf"))
    self.setMaxZ(np.float64("inf"))
    self.setMinX(-np.float64("inf"))
    self.setMinY(-np.float64("inf"))
    self.setMinZ(-np.float64("inf"))

    for surface_id in self._surfaces:
      surface = self._surfaces[surface_id][0]
      halfspace = self._surfaces[surface_id][1]

      max_x = surface.getMaxX(halfspace=halfspace)
      max_y = surface.getMaxY(halfspace=halfspace)
      max_z = surface.getMaxZ(halfspace=halfspace)

      min_x = surface.getMinX(halfspace=halfspace)
      min_y = surface.getMinY(halfspace=halfspace)
      min_z = surface.getMinZ(halfspace=halfspace)

      if max_x != np.float64("inf") and max_x < self._max_x:
        self.setMaxX(max_x)
      if max_y != np.float64("inf") and max_y < self._max_y:
        self.setMaxY(max_y)
      if max_x != np.float64("inf") and max_z < self._max_z:
        self.setMaxZ(max_z)

      if min_x != -np.float64("inf") and min_x > self._min_x:
        self.setMinX(min_x)
      if min_y != -np.float64("inf") and min_y > self._min_y:
        self.setMinY(min_y)
      if min_z != -np.float64("inf") and min_z > self._min_z:
        self.setMinZ(min_z)

    # If we could not find a bounds for any dimension, readjust
    # it to +/- infinity
    if self._max_x == -np.float64("inf"):
      self.setMaxX(np.float64("inf"))
    if self._max_y == -np.float64("inf"):
      self.setMaxY(np.float64("inf"))
    if self._max_z == -np.float64("inf"):
      self.setMaxZ(np.float64("inf"))

    if self._min_x == np.float64("inf"):
      self.setMinX(-np.float64("inf"))
    if self._min_y == np.float64("inf"):
      self.setMinY(-np.float64("inf"))
    if self._min_z == np.float64("inf"):
      self.setMinZ(-np.float64("inf"))


  def containsPoint(self, point):

    if not isinstance(point, Point):
      exit('Unable to determine if point is in Cell ID=%d since %s is not '
           'a Point' % (self._id, str(point)))

    for surface_id in self._surfaces:

      halfspace = self._surfaces[surface_id][1]
      surface = self._surfaces[surface_id][0]

      # Return false if the Point is not in the correct Surface halfspace
      if (surface.evaluate(point) * halfspace) < -on_surface_thresh:
        return False

    # Return true if the Point is in the correct halfspace for each Surface
    return True


  def computeVolumeFraction(self, volume=np.float(1.),
                            num_samples=10000, tolerance=1e-3):

    # Do not recompute the volume fraction if it was already computed
    if self._volume_fraction != 0.:
      return

    if not is_float(volume):
      exit('Unable to compute the volume fraction for Cell '
           'ID=%d since volume=%s is not a floating point '
           'value' % (self._id, str(volume)))

    if not is_integer(num_samples):
      exit('Unable to compute the volume fraction for Cell '
           'ID=%d since num_samples=%s is not an integer '
           'value' % (self._id, str(num_samples)))

    if not is_float(tolerance):
      exit('Unable to compute the volume fraction for Cell '
           'ID=%d since tolerance=%s is not a floating point '
           'value' % (self._id, str(tolerance)))

    from numpy.random import uniform

    # Initialize the point
    point = Point()

    # Compute the volume/area of the bounding box we sample from
    box_volume = np.float64(1.)

    if self._min_x > -np.float64("inf") and self._max_x < np.float64("inf"):
      box_volume *= (self._max_x - self._min_x)
    if self._min_y > -np.float64("inf") and self._max_y < np.float64("inf"):
      box_volume *= (self._max_y - self._min_y)
    if self._min_z > -np.float64("inf") and self._max_z < np.float64("inf"):
      box_volume *= (self._max_z - self._min_z)

    # Initialize variables
    counter = 0.
    tot_samples = 0.
    uncertainty = np.float64("inf")

    while (uncertainty > tolerance):

      # Increment the total samples
      tot_samples += num_samples

      # Draw random samples for x,y,z coordinates
      x = uniform(low=self._min_x, high=self._max_x, size=num_samples)
      y = uniform(low=self._min_y, high=self._max_y, size=num_samples)
      z = uniform(low=self._min_z, high=self._max_z, size=num_samples)

      for i in range(num_samples):

        point.setX(x[i])
        point.setY(y[i])
        point.setZ(z[i])

        if self.containsPoint(point):
          counter += 1.

      fraction = np.float64(counter / tot_samples)

      # Compute uncertainty
      uncertainty = box_volume * np.sqrt((fraction - fraction**2) / tot_samples)

    # Compute the final volume fraction and volume
    self._volume_fraction = np.float64(counter / tot_samples)
    self._volume = np.float64(self._volume_fraction * volume)


  def clone(self):

    clone = Cell()
    clone.setName(self._name)
    clone.setFill(self._fill)

    for surface_id in self._surfaces.keys():
      surface = self._surfaces[surface_id][0]
      halfspace = self._surfaces[surface_id][1]
      clone.addSurface(surface=surface, halfspace=halfspace)

    clone.setMaxX(self._max_x)
    clone.setMaxY(self._max_y)
    clone.setMaxZ(self._max_z)
    clone.setMinX(self._min_x)
    clone.setMinY(self._min_y)
    clone.setMinZ(self._min_z)
    clone.setVolume(self._volume)
    clone.setVolumeFraction(self._volume_fraction)

    return clone


  def toString(self):

    string = ''

    string += 'Cell\n'

    cell_id = '{0: <16}'.format('\tID') + '=\t' + str(self._id)
    string += cell_id + '\n'

    name = '{0: <16}'.format('\tName') + '=\t' + self._name
    string += name + '\n'

    type = '{0: <16}'.format('\tType') + '=\t'
    type += str(self._type)
    string += type + '\n'

    volume_fraction = '{0: <16}'.format('\tVolume Fraction') + '=\t'
    volume_fraction += str(self._volume_fraction)
    string += volume_fraction + '\n'

    volume = '{0: <16}'.format('\tVolume') + '=\t'
    volume += str(self._volume)
    string += volume + '\n'

    fill = '{0: <16}'.format('\tFill') + '=\t'

    if not self._fill is None:

      if self._type is 'material':
        fill += 'Material ID='
      elif self._type is 'universe':
        fill += 'Universe ID='
      else:
        fill += 'Lattice ID='

      fill += str(self._fill.getId())

    string += fill + '\n'

    num_subcells = '{0: <16}'.format('\t# Regions') + '=\t'
    num_subcells += str(self._num_subcells)
    string += num_subcells + '\n'

    surfaces = '{0: <16}'.format('\tSurfaces') + '=\t'

    for surface_id in self._surfaces.keys():
      halfspace = self._surfaces[surface_id][1]
      surfaces += str(surface_id*halfspace) + ', '

    string += surfaces.rstrip(', ') + '\n'

    return string


  def printString(self):
    print(self.toString())



################################################################################
#################################  LocalCoords  ################################
################################################################################


class LocalCoords(object):

  def __init__(self, point=None, next=None, prev=None):

    self._point = None
    self._type = None
    self._next = None
    self._prev = None

    if not point is None:
      self.setPoint(point)

    if not next is None:
      self.setNext(next)

    if not prev is None:
      self.setPrev(prev)


  def getPoint(self):
    return self._point


  def getType(self):
    return self._type


  def getNext(self):
    return self._next


  def getPrev(self):
    return self._prev


  def setPoint(self, point):

    if not isinstance(point, Point):
      exit('Unable to set the point %s for LocalCoords since it is not '
           'a Point object' % str(point))

    self._point = point


  def setNext(self, next):

    if not isinstance(next, LocalCoords) and not next is None:
      exit('Unable to set the next to %s for LocalCoords since it is not '
           'a LocalCoords object' % str(next))

    self._next = next


  def setPrev(self, prev):

    if not isinstance(prev, LocalCoords) and not prev is None:
      exit('Unable to set the prev to %s for LocalCoords since it is not '
           'a LocalCoords object' % str(prev))

    self._prev = prev


  def getHeadNode(self):

    curr = self
    prev = self.getPrev()

    while not prev is None:
      curr = prev
      prev = curr.getPrev()

    return curr


  def getTailNode(self):

    curr = self
    next = self.getNext()

    while not next is None:
      curr = next
      next = curr.getNext()

    return curr


  def prune(self):

    curr = self.getTailNode()
    next = curr.getPrev()

    # Iterate over LocalCoords beneath this one in the linked list
    while curr != self:
      next = curr.getPrev()
      del curr
      curr = next

    # Set the next LocalCoord in the linked list to null
    self.setNext(None)


  def toString(self):

    string = ''

    string += 'LocalCoords\n'

    type = '{0: <16}'.format('\tType') + '=\t' + str(self._type)
    string += type + '\n'

    point = '{0: <16}'.format('\tPoint') + '=\t'

    if not self._point is None:
        point += str(self._point.getCoords())

    string += point + '\n'

    return string


  def printString(self):
    print(self.toString())



class UnivCoords(LocalCoords):

  def __init__(self, point=None, next=None, prev=None,
               universe=None, cell=None):

    super(UnivCoords, self).__init__(point, next, prev)

    self._type = 'universe'
    self._universe = None
    self._cell = None

    if not universe is None:
      self.setUniverse(universe)

    if not cell is None:
      self.setCell(cell)


  def getUniverse(self):
    return self._universe


  def getCell(self):
    return self._cell


  def setUniverse(self, universe):

    if not isinstance(universe, Universe):
      exit('Unable to set the Universe to %s for LocalCoords since it '
           'is not a Universe' % str(universe))

    self._universe = universe


  def setCell(self, cell):

    if not isinstance(cell, Cell):
      exit('Unable to set the Cell to %s for LocalCoords since it '
           'is not a Cell' % str(cell))

    self._cell = cell


  def toString(self):

    string = super(UnivCoords, self).toString()

    universe_id = '{0: <16}'.format('\tUniverse') + '=\t'
    universe_id += str(self._universe.getId()) + '\n'
    string += universe_id

    cell_id = '{0: <16}'.format('\tCell') + '=\t'\

    if not self._cell is None:
      cell_id += str(self._cell.getId())

    string += cell_id

    return string



class LatCoords(LocalCoords):

  def __init__(self, point=None, next=None, prev=None,
               lattice=None, lattice_x=None, lattice_y=None):

    super(LatCoords, self).__init__(point, next, prev)

    self._type = 'lattice'
    self._lattice = None
    self._lattice_x = None
    self._lattice_y = None

    if not lattice is None:
      self.setLattice(lattice)

    if not lattice_x is None:
      self.setLatticeX(lattice_x)

    if not lattice_y is None:
      self.setLatticeY(lattice_y)


  def getLattice(self):
    return self._lattice


  def getLatticeX(self):
    return self._lattice_x


  def getLatticeY(self):
    return self._lattice_y


  def setLattice(self, lattice):

    if not isinstance(lattice, Lattice):
      exit('Unable to set the Lattice to %s for LocalCoords since it '
           'is not a Lattice' % str(lattice))

    self._lattice = lattice


  def setLatticeX(self, lattice_x):

    if not is_integer(lattice_x):
      exit('Unable to set the Lattice X to %s for LocalCoords since it '
           'is not an integer' % str(lattice_x))

    self._lattice_x = lattice_x


  def setLatticeY(self, lattice_y):

    if not is_integer(lattice_y):
      exit('Unable to set the Lattice Y to %s for LocalCoords since it '
           'is not an integer' % str(lattice_y))

    self._lattice_y = lattice_y


  def toString(self):

    string = super(LatCoords, self).toString()

    lattice_id = '{0: <16}'.format('\tLattice') + '=\t'
    lattice_id += str(self._lattice.getId()) + '\n'
    string += lattice_id

    lattice_x = '{0: <16}'.format('\tLattice X') + '=\t'
    lattice_x += str(self._lattice_x) + '\n'
    string += lattice_x

    lattice_y = '{0: <16}'.format('\tLattice Y') + '=\t'
    lattice_y += str(self._lattice_y)
    string += lattice_y

    return string