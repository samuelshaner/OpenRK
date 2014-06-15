__author__ = 'Will Boyd'
__email__ = 'wboyd@mit.edu'


from material import Material
from surface import Surface, on_surface_thresh
from point import Point
from checkvalue import *
import numpy as np
import math
from collections import OrderedDict


# Error threshold for determining how close to the boundary of a Lattice cell
# a Point needs to be to be considered on it
on_lattice_cell_thresh = 1e-12

# Lists of all IDs for all Universes created
universe_ids = list()

# A static variable for auto-generated Universe IDs
auto_universe_id = 10000

max_float = np.finfo(np.float64).max
min_float = np.finfo(np.float64).min


class Universe(object):

  def __init__(self, universe_id=None, name=''):

    # Initialize Universe class attributes
    self._id = None
    self._name = None

    # Keys   - Cell IDs
    # Values - Cells
    self._cells = dict()

    # Keys   - Cell IDs
    # Values - Offsets
    self._cell_offsets = OrderedDict()
    self._num_regions = 0
    self._volume = np.finfo(np.float64)

    self.setId(universe_id)
    self.setName(name)


  def getAllCells(self):

    cells = dict()

    # Add this Universe's cells to the dictionary
    cells.update(self._cells)

    # Append all Cells in each Cell in the Universe to the dictionary
    for cell_id, cell in self._cells.iteritems():
      cells.update(cell.getAllCells())

    return cells


  def getAllUniverses(self):

    # Get all Cells in this Universe
    cells = self.getAllCells()

    universes = dict()

    # Append all Universes containing each Cell to the dictionary
    for cell_id, cell in cells.iteritems():
      universes.update(cell.getAllUniverses())

    return universes


  def getCellOffset(self, cell):

    if not isinstance(cell, Cell):
      msg = 'Unable to get cell offset from Universe ID={0} for {1} which ' \
            'is not a Cell'.format(self._id, cell)
      raise ValueError(msg)

    cell_id = cell._id

    if not cell_id in self._cell_offsets.keys():
      msg = 'Unable to get cell offset from Universe ID={0} for Cell ID={1} ' \
            'which is not in the Universe'.format(self._id, cell_id)
      raise ValueError(msg)

    return self._cell_offsets[cell_id]


  def getMaxX(self):

    max_x = min_float

    for cell_id in self._cells:

      cell = self._cells[cell_id]
      if max_x < cell.getMaxX():
        max_x = cell.getMaxX()

    return max_x


  def getMaxY(self):

    max_y = min_float

    for cell_id in self._cells:

      cell = self._cells[cell_id]
      if max_y < cell.getMaxY():
        max_y = cell.getMaxY()

    return max_y


  def getMaxZ(self):

    max_z = min_float

    for cell_id in self._cells:

      cell = self._cells[cell_id]
      if max_z < cell.getMaxZ():
        max_z = cell.getMaxZ()

    return max_z


  def getMinX(self):

    min_x = max_float

    for cell_id in self._cells:

      cell = self._cells[cell_id]
      if min_x > cell.getMinX():
        min_x = cell.getMinX()

    return min_x


  def getMinY(self):

    min_y = max_float

    for cell_id in self._cells:

      cell = self._cells[cell_id]
      if min_y > cell.getMinY():
        min_y = cell.getMinY()

    return min_y


  def getMinZ(self):

    min_z = max_float

    for cell_id in self._cells:

      cell = self._cells[cell_id]
      if min_z > cell.getMinZ():
        min_z = cell.getMinZ()

    return min_z


  def setMaxX(self, max_x):

    if not is_float(max_x) and not is_integer(max_x):
      msg = 'Unable to set the maximum x for Universe ID={0} to {1} since it ' \
            'is not an integer or floating point value'.format(self._id, max_x)
      raise ValueError(msg)

    for cell_id in self._cells.keys():
      cell = self._cells[cell_id]
      cell.setMaxX(max_x=max_x)


  def setMaxY(self, max_y):

    if not is_float(max_y) and not is_integer(max_y):
      msg = 'Unable to set the maximum y for Universe ID={1} to {0} since it ' \
           'is not an integer or floating point value'.format(self._id, max_y)
      raise ValueError(msg)

    for cell_id in self._cells.keys():
      cell = self._cells[cell_id]
      cell.setMaxY(max_y=max_y)


  def setMaxZ(self, max_z):

    if not is_float(max_z) and not is_integer(max_z):
      msg = 'Unable to set the maximum z for Universe ID={0} to {1} since it ' \
            'is not an integer or floating point value'.format(self._id, max_z)
      raise ValueError(msg)

    for cell_id in self._cells.keys():
      cell = self._cells[cell_id]
      cell.setMaxZ(max_z=max_z)


  def setMinX(self, min_x):

    if not is_float(min_x) and not is_integer(min_x):
      msg = 'Unable to set the minimum x for Universe ID={0} to {1} since it ' \
            'is not an integer or floating point value'.format(self._id, min_x)
      raise ValueError(msg)

    for cell_id in self._cells.keys():
      cell = self._cells[cell_id]
      cell.setMinX(min_x=min_x)


  def setMinY(self, min_y):

    if not is_float(min_y) and not is_integer(min_y):
      msg = 'Unable to set the minimum y for Universe ID={0} to {1} since it ' \
            'is not an integer or floating point value'.format(self._id, min_y)
      raise ValueError(msg)

    for cell_id in self._cells.keys():
      cell = self._cells[cell_id]
      cell.setMinY(min_y=min_y)


  def setMinZ(self, min_z):

    if not is_float(min_z) and not is_integer(min_z):
      msg = 'Unable to set the minimum z for Universe ID={0} to {1} since it ' \
            'is not an integer or floating point value'.format(self._id, min_z)
      raise ValueError(msg)

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

      # If the Universe already has an ID, remove it from global list
      if not self._id is None:
        universe_ids.remove(self._id)

      if universe_id in universe_ids:
        msg = 'Unable to set Universe ID to {0} since a Universe with this ' \
              'ID was already initialized'.format(universe_id)
        raise ValueError(msg)

      if universe_id < 0:
        msg = 'Unable to set Univeres ID to {0} since it must be a ' \
              'non-negative integer'.format(universe_id)
        raise ValueError(msg)

      else:
        self._id = universe_id
        universe_ids.append(universe_id)

    else:
      msg = 'Unable to set Universe ID to a non-integer {0}'.format(universe_id)
      raise ValueError(msg)


  def setName(self, name):

    if not is_string(name):
      msg = 'Unable to set name for Universe ID={0} with a non-string ' \
            'value {1}'.format(self._id, name)
      raise ValueError(msg)

    else:
      self._name = name


  def addCell(self, cell):

    if not isinstance(cell, Cell):
      msg = 'Unable to add a Cell to Universe ID={0} since {1} is ' \
            'not a Cell'.format(self._id, cell)
      raise ValueError(msg)

    cell_id = cell._id

    if not cell_id in self._cells.keys():
      self._cells[cell_id] = cell


  def addCells(self, cells):

    if not isinstance(cells, (list, tuple, np.ndarray)):
      msg = 'Unable to add Cells to Universe ID={0} since {1} is ' \
            'not a Python tuple/list or NumPy array'.format(self._id, cells)
      raise ValueError(msg)

    for i in range(len(cells)):
      self.addCell(cells[i])


  def removeCell(self, cell):

    if not isinstance(cell, Cell):
      msg = 'Unable to remove a Cell from Universe ID={0} since {1} is not a ' \
            'Cell'.format(self._id, cell)
      raise ValueError(msg)

    cell_id = cell._id

    if cell_id in self._cells.keys():
      del self._cells[cell_id]


  def clearCells(self):
    self._cells.clear()


  def computeVolumeFractions(self, volume=np.float64(1.), tolerance=1e-3):

    if not is_float(volume):
      msg = 'Unable to compute volume fractions for Universe ID={0} since ' \
            'volume={1} is not a floating point value'.format(self._id, volume)
      raise ValueError(msg)

    self._volume = volume

    max_x = self.getMaxX()
    max_y = self.getMaxY()
    max_z = self.getMaxZ()
    min_x = self.getMinX()
    min_y = self.getMinY()
    min_z = self.getMinZ()

    for cell_id in self._cells:

      cell = self._cells[cell_id]

      cell.setMaxX(max_x)
      cell.setMaxY(max_y)
      cell.setMaxZ(max_z)
      cell.setMinX(min_x)
      cell.setMinY(min_y)
      cell.setMinZ(min_z)

      cell.computeVolumeFraction(volume=self._volume, tolerance=tolerance)

      # Recursively compute volume fractions below this Universe
      fill = cell._fill

      if isinstance(fill, Universe):
        volume_fraction = cell._volume_fraction
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
      msg = 'Unable to find cell in Universe ID={0} since localcoords ' \
            '{1} is not a LocalCoords'.format(self._id, localcoords)
      raise ValueError(msg)

    for cell_id in self._cells:

      cell = self._cells[cell_id]

      if cell.containsPoint(localcoords._point):
        localcoords.setCell(cell)

        # 'material' type Cell - lowest level, terminate search for Cell
        if cell._type == 'material':
            return cell

        # 'fill' type Cell - Cell contains a Universe at a lower level
        # Update coords to next level and continue search
        else:

          fill = cell._fill

          if isinstance(fill, Lattice):
            next_coords = LatCoords(localcoords._point)
            next_coords.setLattice(fill)

          elif isinstance(fill, Universe):
            next_coords = UnivCoords(localcoords._point)
            next_coords.setUniverse(fill)

          else:
            msg = 'Unable to find cell since in Universe ID={0} ' \
                  'since Cell ID={1} is not filled by a Material, ' \
                  'Universe or Lattice'.format(self._id, cell._id)
            raise ValueError(msg)

          localcoords.setNext(next_coords)
          next_coords.setPrev(localcoords)

          return fill.findCell(next_coords)


#  def findCell(self, region_id):


  def findRegion(self, region_id, univ_coords):

    if not is_integer(region_id):
      msg = 'Unable to find region_id in Universe ID={0} since {1} is ' \
            'not an integer value'.format(self._id, region_id)
      raise ValueError(msg)

    if not isinstance(univ_coords, UnivCoords):
      msg = 'Unable to find region_id in Universe ID={0} since {1} is ' \
            'not a UnivCoords'.format(self._id, univ_coords)
      raise ValueError(msg)

    # Initialize cell and offset
    cell = None
    offset = 0

    # Loop over cells until we reach the one the first one with
    # an offset larger than region_id - return the one prior
    for cell_id in self._cells.keys():

      if self._cell_offsets[cell_id] <= region_id:
        cell = self._cells[cell_id]
        offset = self._cell_offsets[cell_id]

      else:
        break


    if cell is None:
      msg = 'Unable to find region_id={0} in ' \
            'Universe ID={1}'.format(region_id, self._id)
      raise ValueError(msg)


    region_id -= offset
    univ_coords.setCell(cell)
    fill = cell._fill

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
      msg = 'Unable to find cell for region_id={0} in Universe ' \
            'ID={1}'.format(region_id, self._id)
      raise ValueError(msg)


  def __repr__(self):

    string = 'Universe\n'
    string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self._id)
    string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self._name)
    string += '{0: <16}{1}{2}\n'.format('\t# Regions', '=\t', self._num_regions)
    string += '{0: <16}{1}{2}\n'.format('\tCells', '=\t', self._cells.keys())
    return string



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


  def getUniverse(self, lat_x, lat_y):

    if not is_integer(lat_x):
      msg = 'Unable to get Universe from Lattice ID={0} since x={1} is not ' \
            'a lattice cell'.format(self._id, lat_x)
      raise ValueError(msg)

    if not is_integer(lat_y):
      msg = 'Unable to get Universe from Lattice ID={0} since y={1} is not ' \
           'a lattice cell'.format(self._id, lat_y)
      raise ValueError(msg)

    if lat_x < 0 or lat_x > self._dimension[0]:
      msg = 'Unable to get Universe from Lattice ID={0} since x={1} is ' \
            'outside the bounds of the lattice cells'.format(self._id, lat_x)
      raise ValueError(msg)

    if lat_y < 0 or lat_y > self._dimension[1]:
      msg = 'Unable to get Universe from Lattice ID={0} since y={1} is ' \
            'outside the bounds of the lattice cells'.format(self._id, lat_y)
      raise ValueError(msg)

    return self._universes[lat_x][lat_y]


  def getCellOffset(self, lat_x, lat_y):

    if not is_integer(lat_x):
      msg = 'Unable to get cell offset from Lattice ID={0} since ' \
            'lat_x={1} is not a lattice cell'.format(self._id, lat_x)
      raise ValueError(msg)

    if not is_integer(lat_y):
      msg = 'Unable to get cell offset from Lattice ID={0} since ' \
            'lat_y={1} is not a lattice cell'.format(self._id, lat_y)
      raise ValueError(msg)

    if lat_x < 0 or lat_x > self._dimension[0]:
      msg = 'Unable to get Universe from Lattice ID={0} since lat_x={1}' \
           'is outside the bounds of the lattice cells'.format(self._id, lat_x)
      raise ValueError(msg)

    if lat_y < 0 or lat_y > self._dimension[1]:
      msg = 'Unable to get Universe from Lattice ID={0} since lat_y={1} is ' \
            'outside the bounds of the lattice cells'.format(self._id, lat_y)
      raise ValueError(msg)

    return self._cell_offsets[lat_x][lat_y]


  def getMaxX(self):
    return self._lower_left[0] + self._dimension[0] * self._width[0]


  def getMaxY(self):
    return self._lower_left[1] + self._dimension[1] * self._width[1]


  def getMinX(self):
    return self._lower_left[0]


  def getMinY(self):
    return self._lower_left[1]


  def getUniqueUniverses(self):

    if self._universes is None:
      msg = 'Unable to get unique Universes for Lattice ID={0} since ' \
            'the universes array has not been set'.format(self._id)
      raise ValueError(msg)

    universes = dict()

    for i in range(len(self._universes)):
      for j in range(len(self._universes[i])):

        universe = self._universes[i][j]
        universe_id = universe._id
        universes[universe_id] = universe

    return universes


  def getAllCells(self):

    if self._universes is None:
      msg = 'Unable to get all Cells for Lattice ID={0} since the ' \
            'universes array has not been set'.format(self._id)
      raise ValueError(msg)

    cells = dict()

    for i in range(self._dimension[0]):
      for j in range(self._dimension[1]):
        universe = self._universes[i][j]
        cells.update(universe.getAllCells())

    return cells


  def getAllUniverses(self):

    # Initialize a dictionary of all Universes contained by the Lattice
    # in each nested Universe level
    all_universes = dict()

    # Get all unique Universes contained in each of the lattice cells
    unique_universes = self.getUniqueUniverses()

    # Add the unique Universes filling each Lattice cell
    all_universes.update(unique_universes)

    # Append all Universes containing each cell to the dictionary
    for universe_id, universe in unique_universes.iteritems():
      all_universes.update(universe.getAllUniverses())

    return all_universes


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
        msg = 'Unable to set Lattice ID to {0} since a Lattice ' \
              'with this ID was already initialized'.format(lattice_id)
        raise ValueError(msg)

      if lattice_id < 0:
        msg = 'Unable to set Lattice ID to {0} since it must be a ' \
              'non-negative integer'.format(lattice_id)
        raise ValueError(msg)

      else:
        self._id = lattice_id
        universe_ids.append(lattice_id)

    else:
      msg = 'Unable to set a non-integer Lattice ID {0}'.format(lattice_id)
      raise ValueError(msg)


  def setName(self, name):

    if not is_string(name):
      msg = 'Unable to set name for Lattice ID={0} with a ' \
            'non-string value {1}'.format(self._id, name)
      raise ValueError(msg)

    else:
      self._name = name


  def setType(self, type):

    if not is_string(type):
      msg = 'Unable to set the type for Lattice ID={0} with a non-string ' \
            'value {1}'.format(self._id, type)
      raise ValueError(msg)

    if not type in ['rectangular']:
      msg = 'Unable to set the type for Lattice ID={0} to {1} since it ' \
            'is not rectangular'.format(self._id, type)
      raise ValueError(msg)

    self._type = type


  def setDimension(self, dimension):

    if not isinstance(dimension, (tuple, list, np.ndarray)):
      msg = 'Unable to set Lattice ID={0} dimension to {1} since it is not ' \
            'a Python tuple/list or NumPy array'.format(self._id, dimension)
      raise ValueError(msg)

    if len(dimension) != 2 and len(dimension) != 3:
      msg = 'Unable to set Lattice ID={0} dimension to {1} since it does ' \
            'not contain 2 or 3 coordinates'.format(self._id, dimension)
      raise ValueError(msg)

    for dim in dimension:

      if not isinstance(dim, (int, np.int32, np.int64)):
        msg = 'Unable to set the dimension for Lattice ID={0} to {1} ' \
              'since it is not an integer'.format(self._id, dim)
        raise ValueError(msg)

      if dim < 0:
        msg = 'Unable to set Lattice ID={0} dimension to {1} since it ' \
              'is a negative value'.format(self._id, dim)
        raise ValueError(msg)

    self._dimension = np.zeros(len(dimension), dtype=np.int64)

    for i in range(len(dimension)):
      self._dimension[i] = dimension[i]


  def setLowerLeft(self, lower_left):

    if not isinstance(lower_left, (tuple, list, np.ndarray)):
      msg = 'Unable to set the lower_left for Lattice ID={0} to ' \
            '{1} since it is not a Python tuple/list or ' \
            'NumPy array'.format(self._id, lower_left)
      raise ValueError(msg)

    if len(lower_left) != 2 and len(lower_left) != 3:
      msg = 'Unable to set the lower_left for Lattice ID={0} to ' \
            '{1} since it does not contain 2 or 3 ' \
            'coordinates'.format(self._id, lower_left)
      raise ValueError(msg)

    for dim in lower_left:

      if not is_integer(dim) and not is_float(dim):
        msg = 'Unable to set the lower_left for Lattice ID={0} to ' \
              '{1} since it is not an integer or floating ' \
              'point value'.format(self._id, dim)
        raise ValueError(msg)

    self._lower_left = np.zeros(len(lower_left), dtype=np.float64)

    for i in range(len(lower_left)):
      self._lower_left[i] = lower_left[i]


  def setWidth(self, width):

    if not isinstance(width, (tuple, list, np.ndarray)):
      msg = 'Unable to set the width for Lattice ID={0} to {1} since ' \
            'it is not a Python tuple/list or NumPy ' \
            'array'.format(self._id, width)
      raise ValueError(msg)

    if len(width) != 2 and len(width) != 3:
      msg = 'Unable to set the width for Lattice ID={0} to {1} since it does ' \
            'not contain 2 or 3 coordinates'.format(self._id, width)
      raise ValueError(msg)

    for dim in width:

      if not is_integer(dim) and not is_float(dim):
        msg = 'Unable to set the width for Lattice ID={0} to {1} since it is ' \
              'not an integer or floating point value'.format(self._id, dim)
        raise ValueError(msg)

      if dim < 0:
        msg = 'Unable to set the width for Lattice ID={0} to {1} since it ' \
              'is a negative value'.format(self._id, dim)
        raise ValueError(msg)

    self._width = np.zeros(len(width), dtype=np.float64)

    for i in range(len(width)):
      self._width[i] = width[i]


  def setUniverses(self, universes):

    if not isinstance(universes, (tuple, list, np.ndarray)):
      msg = 'Unable to set the universes for Lattice ID={0} to {1} ' \
            'since it is not a Python tuple/list or ' \
            'NumPy array'.format(self._id, universes)
      raise ValueError(msg)

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
      msg = 'Unable to compute volume fractions for Lattice ID={0} since ' \
            'volume={1} is not a floating point value'.format(self._id, volume)
      raise ValueError(msg)

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
          count += self._universes[i][j]._num_regions

    if num_dims == 3:

      for i in range(self._dimension[0]):
        for j in range(self._dimension[1]):
          for k in range(self._dimension[2]):
            self._cell_offsets[i,j,k] = count
            self._universes[i,j,k].initializeCellOffsets()
            count += self._universes[i,j,k]._num_regions


    self._num_regions = count


  def findCell(self, localcoords):

    if not isinstance(localcoords, LocalCoords):
      msg = 'Unable to find cell in Lattice ID={0} since localcoords ' \
            '{1} is not a LocalCoords'.format(self._id, localcoords)
      raise ValueError(msg)

    # Compute the x and y indices for the Lattice cell this coord is in
    point = localcoords._point
    x = point._coords[0]
    y = point._coords[1]

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
      msg = 'Unable to find cell since the lattice indices ({0},{1}) are ' \
            'outside of Lattice ID={1}'.format(lat_x, lat_y, self._id)
      raise ValueError(msg)

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
      msg = 'Unable to find Cell since Lattice ID={0} does ' \
            'not contain a Universe or Lattice in lattice ' \
            'cell ({1}, {2})'.format(self._id, lat_x, lat_y)
      raise ValueError(msg)

    localcoords.setNext(next_coords)
    next_coords.setPrev(localcoords)

    return universe.findCell(next_coords)


  def findRegion(self, region_id, lat_coords):

    if not is_integer(region_id):
      msg = 'Unable to find region_id={0} in Lattice ID={1} since {2} is ' \
            'not an integer value'.format(region_id, self._id, region_id)
      raise ValueError(msg)

    if not isinstance(lat_coords, LatCoords):
      msg = 'Unable to find region_id={0} in Lattice ID={1} since {2} is ' \
            'not a LatCoords'.format(region_id, self._id, lat_coords)
      raise ValueError(msg)

    # Initialize cell and offset
    universe = None
    offset = 0
    lat_x = 0
    lat_y = 0

    # Loop over cells until we reach the one the first one with
    # an offset larger than region_id - return the one prior
    for i in range(len(self._cell_offsets)):
      for j in range(len(self._cell_offsets[i])):

        if self._cell_offsets[i][j] <= region_id:
          offset = self._cell_offsets[i][j]
          universe = self._universes[i][j]
          lat_x = i
          lat_y = j

        elif self._cell_offsets[i][j] > region_id:
          break


    if universe is None:
      msg = 'Unable to find region_id={0} for FSR ID={1} in Lattice ' \
            'ID={2}'.format(region_id, self._id)
      raise ValueError(msg)

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
      msg = 'Unable to find region_id={0} in Lattice ' \
            'ID={1}'.format(region_id, self._id)
      raise ValueError(msg)


  def __repr__(self):

    string = 'Lattice\n'
    string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self._id)
    string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self._name)
    string += '{0: <16}{1}{2}\n'.format('\tType', '=\t', self._type)
    string += '{0: <16}{1}{2}\n'.format('\tDimension', '=\t', self._dimension)
    string += '{0: <16}{1}{2}\n'.format('\tLower Left', '=\t', self._lower_left)
    string += '{0: <16}{1}{2}\n'.format('\tWidth', '=\t', self._width)
    string += '{0: <16}{1}{2}\n'.format('\t# Regions', '=\t', self._num_regions)

    string += '{0: <16}{1}'.format('\tUniverses', '\n')

    for i in range(len(self._universes)):
      string += '\t'

      for j in range(len(self._universes[0])):
        string += '{0} '.format(self._universes[i][j]._id)

      string += '\n'

    return string



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

    self.setMaxX(max_float)
    self.setMaxY(max_float)
    self.setMaxZ(max_float)
    self.setMinX(min_float)
    self.setMinY(min_float)
    self.setMinZ(min_float)

    if not fill is None:
      self.setFill(fill)


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


  def getAllCells(self):

    if self._fill is None:
      msg = 'Unable to get all Cells from Cell ID={0} since the fill ' \
            'has not been set'.format(self._id)
      raise ValueError(msg)

    cells = dict()

    if self._type == 'universe' or self._type == 'lattice':
      cells.update(self._fill.getAllCells())

    return cells


  def getAllUniverses(self):

    if self._fill is None:
      msg = 'Unable to get all Universes from Cell ID={0} since the fill ' \
            'has not been set'.format(self._id)
      raise ValueError(msg)

    universes = dict()

    if self._type == 'universe' or self._type == 'lattice':
      universes[self._fill._id] = self._fill
      universes.update(self._fill.getAllUniverses())

    return universes


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
        msg = 'Unable to set Cell ID to {0} since a Cell with this ID was ' \
              'already initialized'.format(cell_id)
        raise ValueError(msg)

      if cell_id < 0:
        msg = 'Unable to set Cell ID to {0} since it must be a ' \
              'non-negative integer'.format(cell_id)
        raise ValueError(msg)

      else:
        self._id = cell_id
        cell_ids.append(cell_id)

    else:
      msg = 'Unable to set Cell ID to a non-integer {0}'.format(cell_id)
      raise ValueError(msg)


  def setName(self, name):

    if not is_string(name):
      msg = 'Unable to set name for Cell ID={0} with a ' \
            'non-string value {1}'.format(self._id, name)
      raise ValueError(msg)

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
      msg = 'Unable to set fill for Cell ID={0} to {1} since it is not a ' \
            'Universe, Lattice or a Material'.format(self._id, fill)
      raise ValueError(msg)

    self._fill = fill


  def setType(self, type):

    if not is_string(type):
      msg = 'Unable to set the type for Cell ID={0} to {1} since it is not ' \
            'a string'.format(self._id, type)
      raise ValueError(msg)

    if not type.lower() in ['universe', 'lattice', 'material']:
      msg = 'Unable to set the type for Cell ID={0} to {1} since it is not ' \
            'universe, lattice or material'.format(self._id, type)
      raise ValueError(msg)

    self._type = type.lower()


  def setMaxX(self, max_x):

    if not is_float(max_x):
      msg = 'Unable to set the maximum x-coordinate for Cell ID={0} to {1} ' \
            'since it is not a floating point value'.format(self._id, max_x)
      raise ValueError(msg)

    self._max_x = max_x


  def setMaxY(self, max_y):

    if not is_float(max_y):
      msg = 'Unable to set the maximum y-coordinate for Cell ID={0} to {1} ' \
            'since it is not a floating point value'.format(self._id, max_y)
      raise ValueError(msg)

    self._max_y = max_y


  def setMaxZ(self, max_z):

    if not is_float(max_z):
      msg = 'Unable to set the maximum z-coordinate for Cell ID={0} to {1} ' \
            'since it is not a floating point value'.format(self._id, max_z)
      raise ValueError(msg)

    self._max_z = max_z


  def setMinX(self, min_x):

    if not is_float(min_x):
      msg = 'Unable to set the minimum x-coordinate for Cell ID={0} to {1} ' \
            'since it is not a floating point value'.format(self._id, min_x)
      raise ValueError(msg)

    self._min_x = min_x


  def setMinY(self, min_y):

    if not is_float(min_y):
      msg = 'Unable to set the minimum y-coordinate for Cell ID={0} to {1} ' \
           'since it is not a floating point value'.format(self._id, min_y)
      raise ValueError(msg)

    self._min_y = min_y


  def setMinZ(self, min_z):

    if not is_float(min_z):
      msg = 'Unable to set the minimum z-coordinate for Cell ID={0} to {1} ' \
           'since it is not a floating point value'.format(self._id, min_z)
      raise ValueError(msg)

    self._min_z = min_z


  def addSurface(self, surface, halfspace):

    if not isinstance(surface, Surface):
      msg = 'Unable to add a Surface to Cell ID={0} since {1} is ' \
            'not a Surface'.format(self._id, surface)
      raise ValueError(msg)

    if not halfspace in [-1, +1]:
      msg = 'Unable to add a Surface to Cell ID={0} with halfspace {1} ' \
            'since it is not +/-1'.format(self._id, halfspace)
      raise ValueError(msg)

    # Only add the Surface if the Cell does not already contain another
    # Surface of the same type and coefficients (a redundant Surface)
    for surface_id in self._surfaces.keys():

      test_surface = self._surfaces[surface_id][0]
      test_halfspace = self._surfaces[surface_id][1]

      if halfspace == test_halfspace:

        if surface._type == test_surface._type:
          coeffs = surface._coeffs
          test_coeffs = test_surface._coeffs

          for key in coeffs.keys():
            coeff = coeffs[key]
            test_coeff = test_coeffs[key]

          if abs(coeff-test_coeff) < 1e-10:
            return

    surface_id = surface._id

    if not surface_id in self._surfaces.keys():
      self._surfaces[surface_id] = (surface, halfspace)

    self.findBoundingBox()


  def addSurfaces(self, surfaces, halfspaces):

    if not isinstance(surfaces, (list, tuple, np.ndarray)):
      msg = 'Unable to add Surfaces to Cell ID={0} since {1} is not a Python ' \
            'tuple/list or NumPy array'.format(self._id, surfaces)
      raise ValueError(msg)

    if not isinstance(halfspaces, (list, tuple, np.ndarray)):
      msg = 'Unable to add Surfaces to Cell ID={0} since {1} is not a Python ' \
            'tuple/list or NumPy array'.format(self._id, halfspaces)
      raise ValueError(msg)

    if len(surfaces) != len(halfspaces):
      msg = 'Unable to add Surfaces to Cell ID={0} since the number ' \
            'of Surfaces ({1}) and halfspaces ({2}) are not ' \
            'equal'.format(self._id, len(surfaces), len(halfspaces))
      raise ValueError(msg)

    for i in range(len(surfaces)):
      self.addSurface(surfaces[i], halfspaces[i])


  def removeSurface(self, surface):

    if not isinstance(surface, Surface):
      msg = 'Unable to remove a surface from Cell ID={0} since {1} is not a ' \
            'Surface'.format(self._id, surface)
      raise ValueError(msg)

    surf_id = surface._id

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
      msg = 'Unable to set the volume fraction for Cell ID={0} to {1} since ' \
            'it is not a floating point value'.format(self._id, volume_fraction)
      raise ValueError(msg)

    self._volume_fraction = volume_fraction


  def setVolume(self, volume):

    if not is_float(volume):
      msg = 'Unable to set the volume for Cell ID={0} to {1} since it ' \
            'is not a floating point value'.format(self._id, volume)
      raise ValueError(msg)

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
      self._num_subcells = self._fill._num_regions

    else:
      msg = 'Unable to compute the number of subcells for ' \
            'Cell ID={0} since it is not filled by a Material, ' \
            'Universe or Lattice'.format(self._id)
      raise ValueError(msg)

    return self._num_subcells


  def findBoundingBox(self):

    self.setMaxX(max_float)
    self.setMaxY(max_float)
    self.setMaxZ(max_float)
    self.setMinX(min_float)
    self.setMinY(min_float)
    self.setMinZ(min_float)

    for surface_id in self._surfaces:
      surface = self._surfaces[surface_id][0]
      halfspace = self._surfaces[surface_id][1]

      max_x = surface.getMaxX(halfspace=halfspace)
      max_y = surface.getMaxY(halfspace=halfspace)
      max_z = surface.getMaxZ(halfspace=halfspace)

      min_x = surface.getMinX(halfspace=halfspace)
      min_y = surface.getMinY(halfspace=halfspace)
      min_z = surface.getMinZ(halfspace=halfspace)

      if max_x != max_float and max_x < self._max_x:
        self.setMaxX(max_x)
      if max_y != max_float and max_y < self._max_y:
        self.setMaxY(max_y)
      if max_z != max_float and max_z < self._max_z:
        self.setMaxZ(max_z)

      if min_x != min_float and min_x > self._min_x:
        self.setMinX(min_x)
      if min_y != min_float and min_y > self._min_y:
        self.setMinY(min_y)
      if min_z != min_float and min_z > self._min_z:
        self.setMinZ(min_z)

    # If we could not find a bounds for any dimension, readjust
    # it to +/- infinity
    if self._max_x == min_float:
      self.setMaxX(max_float)
    if self._max_y == min_float:
      self.setMaxY(max_float)
    if self._max_z == min_float:
      self.setMaxZ(max_float)

    if self._min_x == max_float:
      self.setMinX(min_float)
    if self._min_y == max_float:
      self.setMinY(min_float)
    if self._min_z == max_float:
      self.setMinZ(min_float)


  def containsPoint(self, point):

    if not isinstance(point, Point):
      msg = 'Unable to determine if point is in Cell ID={1} since {0} is not ' \
            'a Point'.format(self._id, point)
      raise ValueError(msg)

    for surface_id in self._surfaces:

      halfspace = self._surfaces[surface_id][1]
      surface = self._surfaces[surface_id][0]

      # Return false if the Point is not in the correct Surface halfspace
      if (halfspace * surface.evaluate(point)) < 0.:
        return False

    # Return true if the Point is in the correct halfspace for each Surface
    return True


  def computeVolumeFraction(self, volume=np.float(1.),
                            num_samples=1000, tolerance=1e-3):

    # Do not recompute the volume fraction if it was already computed
    if self._volume_fraction != 0.:
      return

    if not is_float(volume):
      msg = 'Unable to compute the volume fraction for Cell ' \
            'ID={0} since volume={1} is not a floating point ' \
            'value'.format(self._id, volume)
      raise ValueError(msg)

    if not is_integer(num_samples):
      msg = 'Unable to compute the volume fraction for Cell ' \
            'ID={0} since num_samples={1} is not an integer ' \
            'value'.format(self._id, num_samples)
      raise ValueError(msg)

    if not is_float(tolerance):
      msg = 'Unable to compute the volume fraction for Cell ' \
            'ID={0} since tolerance={1} is not a floating point ' \
            'value'.format(self._id, tolerance)
      raise ValueError(msg)

    from numpy.random import uniform

    # Initialize the point
    point = Point()

    # Compute the volume/area of the bounding box we sample from
    box_volume = np.float64(1.)

    if self._min_x > min_float and self._max_x < max_float:
      box_volume *= (self._max_x - self._min_x)
    if self._min_y > min_float and self._max_y < max_float:
      box_volume *= (self._max_y - self._min_y)
    if self._min_z > min_float and self._max_z < max_float:
      box_volume *= (self._max_z - self._min_z)

    # Initialize variables
    counter = 0.
    tot_samples = 0.
    uncertainty = max_float

    while (uncertainty > tolerance):

      # Increment the total samples
      tot_samples += num_samples

      # Draw random samples for x,y,z coordinates
      x = uniform(low=self._min_x, high=self._max_x, size=num_samples)
      y = uniform(low=self._min_y, high=self._max_y, size=num_samples)
      z = uniform(low=self._min_z, high=self._max_z, size=num_samples)

      # Must replace all NaNs with zeros
      x = np.nan_to_num(x)
      y = np.nan_to_num(y)
      z = np.nan_to_num(z)

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


  def __repr__(self):

    string = 'Cell\n'
    string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self._id)
    string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self._name)
    string += '{0: <16}{1}{2}\n'.format('\tType', '=\t', self._type)
    string += '{0: <16}{1}{2}\n'.format('\tVol. Frac.', '=\t', self._volume_fraction)
    string += '{0: <16}{1}{2}\n'.format('\tVolume', '=\t', self._volume)

    string += '{0: <16}{1}'.format('\tFill', '=\t')

    if not self._fill is None:

      if self._type is 'material':
        string += 'Material ID='
      elif self._type is 'universe':
        string += 'Universe ID='
      else:
        string += 'Lattice ID='

      string += str(self._fill._id)

    string += '\n'

    string += '{0: <16}{1}{2}\n'.format('\t# Regions', '=\t', self._num_subcells)

    string += '{0: <16}{1}'.format('\tSurfaces', '=\t')

    for surface_id in self._surfaces.keys():
      halfspace = self._surfaces[surface_id][1]
      string += str(surface_id * halfspace) + ', '

    string = string.rstrip(', ') + '\n'

    return string



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


  def setPoint(self, point):

    if not isinstance(point, Point):
      msg = 'Unable to set the point {0} for LocalCoords since it is not ' \
            'a Point object'.format(point)
      raise ValueError(msg)

    self._point = point


  def setNext(self, next):

    if not isinstance(next, LocalCoords) and not next is None:
      msg = 'Unable to set the next to {0} for LocalCoords since it is not ' \
            'a LocalCoords object'.format(next)
      raise ValueError(msg)

    self._next = next


  def setPrev(self, prev):

    if not isinstance(prev, LocalCoords) and not prev is None:
      msg = 'Unable to set the prev to {0} for LocalCoords since it is not ' \
            'a LocalCoords object'.format(prev)
      raise ValueError(msg)

    self._prev = prev


  def getHeadNode(self):

    curr = self
    prev = self._prev

    while not prev is None:
      curr = prev
      prev = curr._prev

    return curr


  def getTailNode(self):

    curr = self
    next = self._next

    while not next is None:
      curr = next
      next = curr._next

    return curr


  def prune(self):

    curr = self.getTailNode()
    next = curr._prev

    # Iterate over LocalCoords beneath this one in the linked list
    while curr != self:
      next = curr._prev
      del curr
      curr = next

    # Set the next LocalCoord in the linked list to null
    self.setNext(None)


  def __repr__(self):

    string = 'LocalCoords\n'
    string += '{0: <16}{1}{2}\n'.format('\tType', '=\t', self._type)

    string += '{0: <16}{1}'.format('\tPoint', '=\t')

    if not self._point is None:
      string += '{0}'.format(self._point._coords)

    string += '\n'

    return string



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


  def setUniverse(self, universe):

    if not isinstance(universe, Universe):
      msg = 'Unable to set the Universe to {0} for LocalCoords since it ' \
            'is not a Universe'.format(universe)
      raise ValueError(msg)

    self._universe = universe


  def setCell(self, cell):

    if not isinstance(cell, Cell):
      msg = 'Unable to set the Cell to {0} for LocalCoords since it ' \
            'is not a Cell'.format(cell)
      raise ValueError(msg)

    self._cell = cell


  def __repr__(self):

    string = super(UnivCoords, self).__repr__()

    string += '{0: <16}{1}{2}\n'.format('\tUniverse', '=\t', self._universe._id)

    string += '{0: <16}{1}'.format('\tCell', '=\t')

    if not self._cell is None:
      string += '{0}'.format(self._cell._id)

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


  def setLattice(self, lattice):

    if not isinstance(lattice, Lattice):
      msg = 'Unable to set the Lattice to {0} for LocalCoords since it ' \
            'is not a Lattice'.format(lattice)
      raise ValueError(msg)

    self._lattice = lattice


  def setLatticeX(self, lattice_x):

    if not is_integer(lattice_x):
      msg = 'Unable to set the Lattice X to {0} for LocalCoords since it ' \
            'is not an integer'.format(lattice_x)
      raise ValueError(msg)

    self._lattice_x = lattice_x


  def setLatticeY(self, lattice_y):

    if not is_integer(lattice_y):
      msg = 'Unable to set the Lattice Y to {0} for LocalCoords since it ' \
            'is not an integer'.format(lattice_y)
      raise ValueError(msg)

    self._lattice_y = lattice_y


  def __repr__(self):

    string = super(LatCoords, self).__repr__()

    string += '{0: <16}{1}{2}\n'.format('\tLattice', '=\t', self._lattice._id)
    string += '{0: <16}{1}{2}\n'.format('\tLattice X', '=\t', self._lattice_x)
    string += '{0: <16}{1}{2}\n'.format('\tLattice Y', '=\t', self._lattice_y)

    return string