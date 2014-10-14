__author__ = 'Will Boyd'
__email__ = 'wboyd@mit.edu'

import warnings

from opencsg.material import Material
from opencsg.surface import Surface, ON_SURFACE_THRESH
from opencsg.point import Point
from opencsg.checkvalue import *
from collections import OrderedDict
from hashlib import sha1
import numpy as np
from numpy.lib.stride_tricks import as_strided
import copy, math


# Error threshold for determining how close to the boundary of a Lattice cell
# a Point needs to be to be considered on it
ON_LATTICE_CELL_THRESH = 1e-12

# Lists of all IDs for all Universes created
UNIVERSE_IDS = list()

# A static variable for auto-generated Universe IDs
AUTO_UNIVERSE_ID = 10000

def reset_auto_universe_id():
  global AUTO_UNIVERSE_ID, UNIVERSE_IDS
  AUTO_UNIVERSE_ID = 10000
  UNIVERSE_IDS = list()

MAX_FLOAT = np.finfo(np.float64).max
MIN_FLOAT = np.finfo(np.float64).min


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


  def __deepcopy__(self, memo):

    existing = memo.get(id(self))

    # If this is the first time we have tried to copy this object, create a copy
    if existing is None:

      clone = type(self).__new__(type(self))
      clone._id = self._id
      clone._name = self._name

      clone._cells = dict()
      for cell_id in self._cells:
        clone._cells[cell_id] = copy.deepcopy(self._cells[cell_id], memo)

      clone._cell_offsets = copy.deepcopy(self._cell_offsets, memo)
      clone._num_regions = self._num_regions
      clone._volume = self._volume

      memo[id(self)] = clone

      return clone

    # If this object has been copied before, return the first copy made
    else:
      return existing


  def __gt__(self, other):
    return (id(self) > id(other))


  def __ge__(self, other):
    return (id(self) >= id(other))


  def __lt__(self, other):
    return (id(self) < id(other))


  def __le__(self, other):
    return (id(self) <= id(other))


  def getAllCells(self):

    cells = dict()

    # Add this Universe's cells to the dictionary
    cells.update(self._cells)

    # Append all Cells in each Cell in the Universe to the dictionary
    for cell_id, cell in self._cells.items():
      cells.update(cell.getAllCells())

    return cells


  def getAllUniverses(self):

    # Get all Cells in this Universe
    cells = self.getAllCells()

    universes = dict()

    # Append all Universes containing each Cell to the dictionary
    for cell_id, cell in cells.items():
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

    max_x = MIN_FLOAT

    for cell_id in self._cells:

      cell = self._cells[cell_id]
      if max_x < cell.getMaxX():
        max_x = cell.getMaxX()

    return max_x


  def getMaxY(self):

    max_y = MIN_FLOAT

    for cell_id in self._cells:

      cell = self._cells[cell_id]
      if max_y < cell.getMaxY():
        max_y = cell.getMaxY()

    return max_y


  def getMaxZ(self):

    max_z = MIN_FLOAT

    for cell_id in self._cells:

      cell = self._cells[cell_id]
      if max_z < cell.getMaxZ():
        max_z = cell.getMaxZ()

    return max_z


  def getMinX(self):

    min_x = MAX_FLOAT

    for cell_id in self._cells:

      cell = self._cells[cell_id]
      if min_x > cell.getMinX():
        min_x = cell.getMinX()

    return min_x


  def getMinY(self):

    min_y = MAX_FLOAT

    for cell_id in self._cells:

      cell = self._cells[cell_id]
      if min_y > cell.getMinY():
        min_y = cell.getMinY()

    return min_y


  def getMinZ(self):

    min_z = MAX_FLOAT

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

    global UNIVERSE_IDS

    if universe_id is None:
      global AUTO_UNIVERSE_ID
      self._id = AUTO_UNIVERSE_ID
      UNIVERSE_IDS.append(AUTO_UNIVERSE_ID)
      AUTO_UNIVERSE_ID += 1

    # Check that the ID is an integer and wasn't already used
    elif is_integer(universe_id):

      # If the Universe already has an ID, remove it from global list
      if not self._id is None:
        UNIVERSE_IDS.remove(self._id)

      if universe_id in UNIVERSE_IDS:
        msg = 'Unable to set Universe ID to {0} since a Universe with this ' \
              'ID was already initialized'.format(universe_id)
        raise ValueError(msg)

      if universe_id < 0:
        msg = 'Unable to set Univeres ID to {0} since it must be a ' \
              'non-negative integer'.format(universe_id)
        raise ValueError(msg)

      else:
        self._id = universe_id
        UNIVERSE_IDS.append(universe_id)

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


  def containsCell(cell=None, cell_id=None, name=None):

    if not cell is None:
      for cell_id in self._cells.keys():
        if cell == self._cells[cell_id]:
          return True

    if not cell_id is None:
      if cell_id in self._cells.keys():
        return True

    if not name is None:
      for cell_id in self._cells.keys():
        if cell._name == self._cells[cell_id]._name:
          return True

    return False


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


  def buildNeighbors(self):

    # Loop over all of the Universe's Cells
    for cell_id, cell in self._cells.items():

      # Make recursive call for the Cell to build its neighbor Cells
      cell.buildNeighbors()


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

          # Make recursive call to next nested Universe level
          return fill.findCell(next_coords)

    # If we reach this point, we did not find a Cell containing the Point
    return None


  def findRegion(self, region_id, univ_coords):

    if not is_integer(region_id):
      msg = 'Unable to find region_id in Universe ID={0} since {1} is ' \
            'not an integer value'.format(self._id, region_id)
      raise ValueError(msg)

    elif not isinstance(univ_coords, UnivCoords):
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


  def toString(self):
    string = self.__repr__()
    for cell_id, cell in self._cells.items():
      string += cell.toString()
    return string


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
    self._width = None
    self._universes = None
    self._cell_offsets = None
    self._num_regions = None
    self._offset = np.zeros(3, dtype=np.float64)
    self._neighbor_universes = None
    self._neighbor_depth = 1

    self.setId(lattice_id)
    self.setName(name)
    self.setType(type)


  def __deepcopy__(self, memo):

    existing = memo.get(id(self))

    # If this is the first time we have tried to copy this object, create a copy
    if existing is None:

      clone = type(self).__new__(type(self))
      clone._id = self._id
      clone._name = self._name
      clone._type = self._type
      clone._dimension = self._dimension
      clone._width = self._width

      clone._universes = np.empty(self._universes.shape, dtype=Universe)
      for i in range(self._dimension[0]):
        for j in range(self._dimension[1]):
          for k in range(self._dimension[2]):
            universe = self._universes[k][j][i]
            clone._universes[k,j,i] = copy.deepcopy(universe, memo)

      clone._cell_offsets = copy.deepcopy(self._cell_offsets, memo)
      clone._num_regions = self._num_regions
      clone._offset = copy.deepcopy(self._offset, memo)
      clone._neighbor_universes = copy.deepcopy(self._neighbor_universes, memo)
      clone._neighbor_depth = self._neighbor_depth

      memo[id(self)] = clone

      return clone

    # If this object has been copied before, return the first copy made
    else:
      return existing


  def getUniverse(self, lat_x, lat_y, lat_z=None):

    if not is_integer(lat_x):
      msg = 'Unable to get Universe from Lattice ID={0} since x={1} is not ' \
            'a lattice cell'.format(self._id, lat_x)
      raise ValueError(msg)

    elif not is_integer(lat_y):
      msg = 'Unable to get Universe from Lattice ID={0} since y={1} is not ' \
           'a lattice cell'.format(self._id, lat_y)
      raise ValueError(msg)

    elif lat_z != None and not is_integer(lat_y):
      msg = 'Unable to get Universe from Lattice ID={0} since z={1} is not ' \
           'a lattice cell'.format(self._id, lat_z)
      raise ValueError(msg)

    elif lat_x < 0 or lat_x > self._dimension[0]:
      msg = 'Unable to get Universe from Lattice ID={0} since x={1} is ' \
            'outside the bounds of the lattice cells'.format(self._id, lat_x)
      raise ValueError(msg)

    elif lat_y < 0 or lat_y > self._dimension[1]:
      msg = 'Unable to get Universe from Lattice ID={0} since y={1} is ' \
            'outside the bounds of the lattice cells'.format(self._id, lat_y)
      raise ValueError(msg)

    elif lat_z != None and (lat_y < 0 or lat_y > self._dimension[1]):
      msg = 'Unable to get Universe from Lattice ID={0} since z={1} is ' \
            'outside the bounds of the lattice cells'.format(self._id, lat_z)
      raise ValueError(msg)

    if lat_z is None:
      lat_z = 0

    return self._universes[lat_z][lat_y][lat_x]


  def getCellOffset(self, lat_x, lat_y, lat_z=None):

    if not is_integer(lat_x):
      msg = 'Unable to get cell offset from Lattice ID={0} since ' \
            'lat_x={1} is not a lattice cell'.format(self._id, lat_x)
      raise ValueError(msg)

    elif not is_integer(lat_y):
      msg = 'Unable to get cell offset from Lattice ID={0} since ' \
            'lat_y={1} is not a lattice cell'.format(self._id, lat_y)
      raise ValueError(msg)

    elif lat_z != None and not is_integer(lat_y):
      msg = 'Unable to get Universe from Lattice ID={0} since z={1} is not ' \
           'a lattice cell'.format(self._id, lat_z)
      raise ValueError(msg)

    elif lat_x < 0 or lat_x > self._dimension[0]:
      msg = 'Unable to get cell offset from Lattice ID={0} since lat_x={1}' \
           'is outside the bounds of the lattice cells'.format(self._id, lat_x)
      raise ValueError(msg)

    elif lat_y < 0 or lat_y > self._dimension[1]:
      msg = 'Unable to get cell offset from Lattice ID={0} since lat_y={1} is ' \
            'outside the bounds of the lattice cells'.format(self._id, lat_y)
      raise ValueError(msg)

    elif lat_z != None and (lat_y < 0 or lat_y > self._dimension[1]):
      msg = 'Unable to get cell offset from Lattice ID={0} since z={1} is ' \
            'outside the bounds of the lattice cells'.format(self._id, lat_z)
      raise ValueError(msg)

    if lat_z is None:
      lat_z = 0

    return self._cell_offsets[lat_z][lat_y][lat_x]


  def getNeighbors(self, lat_x, lat_y, lat_z=None):
    """Return d-th neighbors of cell (i, j)"""

    if not is_integer(lat_x):
      msg = 'Unable to get neighbor Universes from Lattice ID={0} since ' \
            'x={1} is not a lattice cell'.format(self._id, lat_x)
      raise ValueError(msg)

    elif not is_integer(lat_y):
      msg = 'Unable to get neighbor Universes from Lattice ID={0} since ' \
            'y={1} is not a lattice cell'.format(self._id, lat_y)
      raise ValueError(msg)

    elif lat_z != None and not is_integer(lat_y):
      msg = 'Unable to get neighbor Universes from Lattice ID={0} since ' \
            'z={1} is not a lattice cell'.format(self._id, lat_z)
      raise ValueError(msg)

    elif lat_x < 0 or lat_x > self._dimension[0]:
      msg = 'Unable to get neighbor Universes from Lattice ID={0} since ' \
            'x={1} is outside the bounds of the lattice'.format(self._id, lat_x)
      raise ValueError(msg)

    elif lat_y < 0 or lat_y > self._dimension[1]:
      msg = 'Unable to get neighbor Universes from Lattice ID={0} since ' \
            'y={1} is outside the bounds of the lattice'.format(self._id, lat_y)
      raise ValueError(msg)

    elif lat_z != None and (lat_y < 0 or lat_y > self._dimension[1]):
      msg = 'Unable to get neighbor Universes from Lattice ID={0} since ' \
            'z={1} is outside the bounds of the lattice'.format(self._id, lat_z)
      raise ValueError(msg)

    elif self._universes is None:
      msg = 'Unable to get neighbor Universes for Lattice ID={0} since ' \
            'the universes array has not been set'.format(self._id)
      raise ValueError(msg)

    # Assign the z-coordinate to zero for 2D Lattices
    if lat_z is None:
      lat_z = 0


    depth = self._neighbor_depth

    # Compute indices for Lattice cell entry in the neighbor universes array
    ix = np.clip(lat_z - depth, 0, self._neighbor_universes.shape[0]-1)
    jx = np.clip(lat_y - depth, 0, self._neighbor_universes.shape[1]-1)
    kx = np.clip(lat_x - depth, 0, self._neighbor_universes.shape[2]-1)

    # Compute the neighbor cell starting indices for each dimension
    i0 = max(0, lat_z - depth - ix)
    j0 = max(0, lat_y - depth - jx)
    k0 = max(0, lat_x - depth - kx)

    # Compute the neighbor cell ending indices for each dimension
    i1 = self._neighbor_universes.shape[3] - max(0, depth - lat_z + ix)
    j1 = self._neighbor_universes.shape[4] - max(0, depth - lat_y + jx)
    k1 = self._neighbor_universes.shape[5] - max(0, depth - lat_x + kx)

    # Account for cases where the lattice cell dimensions is 1
    if i0 == 0 and i1 == 0:
      i1 = 1
    if j0 == 0 and j1 == 0:
      j1 = 1
    if k0 == 0 and k1 == 0:
      k1 = 1

    neighbor_universes = self._neighbor_universes[ix,jx,kx][i0:i1,j0:j1,k0:k1]

    # Convert 3D neighbor universes array to a 3D tuple of tuples of tuples

    # Build a 2D tuple of all of the neighbors
    # Iterate over all orderings of the (x,y,z) coordinates to ensure
    # that neighbors tuple is rotationally symmetric across all dimensions

    # x-y-z
    xy_neighbors = neighbor_universes

    # x-z-y
    xz_neighbors = np.swapaxes(neighbor_universes, 0, 1)

    # y-x-z
    yx_neighbors = np.swapaxes(neighbor_universes, 2, 1)

    # y-z-x
    yz_neighbors = np.swapaxes(neighbor_universes, 2, 0)
    yz_neighbors = np.swapaxes(yz_neighbors, 2, 1)

    # z-x-y
    zx_neighbors = np.swapaxes(neighbor_universes, 2, 0)
    zx_neighbors = np.swapaxes(zx_neighbors, 1, 0)

    # z-y-x
    zy_neighbors = np.swapaxes(neighbor_universes, 2, 0)

    neighbors_universes_array = np.array([xy_neighbors.ravel(), xz_neighbors.ravel(),
                                          yx_neighbors.ravel(), yz_neighbors.ravel(),
                                          zx_neighbors.ravel(), zy_neighbors.ravel()])

    return neighbors_universes_array


  def getUniqueNeighbors(self, lat_x, lat_y, lat_z=None):

    # Get the depth x depth x depth array of neighbors and
    # convert to a 1D tuple containing only unique Universes
    neighbor_universes = self.getNeighbors(lat_x, lat_y, lat_z)
    unique_neighbors = np.unique(neighbor_universes)
    return unique_neighbors


  def getNeighborsHash(self, lat_x, lat_y, lat_z=None):

    neighbors_universes = self.getNeighbors(lat_x, lat_y, lat_z)
    neighbors_universes.sort()
    raw_data = neighbors_universes.view()
    neighbors_hash = int(sha1(raw_data).hexdigest(), 16)
    return neighbors_hash


  def getUniqueNeighborsHash(self, lat_x, lat_y, lat_z=None):

    neighbors_universes = self.getUniqueNeighbors(lat_x, lat_y, lat_z)
    neighbors_universes.sort()
    raw_data = neighbors_universes.view()
    neighbors_hash = int(sha1(raw_data).hexdigest(), 16)
    return neighbors_hash


  def getMaxX(self):
    return self._offset[0] + self._dimension[0] / 2.0 * self._width[0]


  def getMaxY(self):
    return self._offset[1] + self._dimension[1] / 2.0 * self._width[1]


  def getMaxZ(self):
    return self._offset[2] + self._dimension[2] / 2.0 * self._width[2]


  def getMinX(self):
    return self._offset[0] - self._dimension[0] / 2.0 * self._width[0]


  def getMinY(self):
    return self._offset[1] - self._dimension[1] / 2.0 * self._width[1]


  def getMinZ(self):
    return self._offset[2] - self._dimension[2] / 2.0 * self._width[2]


  def getUniqueUniverses(self):

    if self._universes is None:
      msg = 'Unable to get unique Universes for Lattice ID={0} since ' \
            'the universes array has not been set'.format(self._id)
      raise ValueError(msg)

    unique_universes = np.unique(self._universes.ravel())
    universes = dict()

    for universe in unique_universes:
      universes[universe._id] = universe

    return universes


  def getAllCells(self):

    cells = dict()
    unique_universes = self.getUniqueUniverses()

    for universe_id, universe in unique_universes.items():
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
    for universe_id, universe in unique_universes.items():
      all_universes.update(universe.getAllUniverses())

    return all_universes


  def setId(self, lattice_id=None):

    global UNIVERSE_IDS

    if lattice_id is None:
      global AUTO_UNIVERSE_ID
      self._id = AUTO_UNIVERSE_ID
      UNIVERSE_IDS.append(AUTO_UNIVERSE_ID)
      AUTO_UNIVERSE_ID += 1

    # Check that the ID is an integer and wasn't already used
    elif is_integer(lattice_id):

      # If the Lattice already has an ID, remove it from global list
      if not self._id is None:
        UNIVERSE_IDS.remove(self._id)

      if lattice_id in UNIVERSE_IDS:
        msg = 'Unable to set Lattice ID to {0} since a Lattice ' \
              'with this ID was already initialized'.format(lattice_id)
        raise ValueError(msg)

      if lattice_id < 0:
        msg = 'Unable to set Lattice ID to {0} since it must be a ' \
              'non-negative integer'.format(lattice_id)
        raise ValueError(msg)

      else:
        self._id = lattice_id
        UNIVERSE_IDS.append(lattice_id)

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

    elif not type in ['rectangular']:
      msg = 'Unable to set the type for Lattice ID={0} to {1} since it ' \
            'is not rectangular'.format(self._id, type)
      raise ValueError(msg)

    self._type = type


  def setOffset(self, offset):
    
    if not isinstance(offset, (tuple, list, np.ndarray)):
      msg = 'Unable to set Lattice ID={0} offset to {1} since it ' \
          'is not a Python tuple/list or NumPy array'.\
          format(self._id, offset)
      raise ValueError(msg)

    elif len(offset) != 3 and len(offset) != 2:
      msg = 'Unable to set Lattice ID={0} offset to {1} since it ' \
            'does not contain 2 or 3 coordinates'.\
            format(self._id, offset)
      raise ValueError(msg)

    for val in offset:

      if not is_float(val):
        msg = 'Unable to set the offset for Lattice ID={0} to {1} ' \
              'since it is not a float'.format(self._id, val)
        raise ValueError(msg)

    for i in range(len(offset)):
      self._offset[i] = offset[i]


  def setDimension(self, dimension):

    if not isinstance(dimension, (tuple, list, np.ndarray)):
      msg = 'Unable to set Lattice ID={0} dimension to {1} since it is not ' \
            'a Python tuple/list or NumPy array'.format(self._id, dimension)
      raise ValueError(msg)

    elif len(dimension) != 2 and len(dimension) != 3:
      msg = 'Unable to set Lattice ID={0} dimension to {1} since it does ' \
            'not contain 2 or 3 coordinates'.format(self._id, dimension)
      raise ValueError(msg)

    for dim in dimension:

      if not isinstance(dim, (int, np.int32, np.int64)):
        msg = 'Unable to set the dimension for Lattice ID={0} to {1} ' \
              'since it is not an integer'.format(self._id, dim)
        raise ValueError(msg)

      elif dim < 0:
        msg = 'Unable to set Lattice ID={0} dimension to {1} since it ' \
              'is a negative value'.format(self._id, dim)
        raise ValueError(msg)

    self._dimension = np.ones(3, dtype=np.int64)

    for i in range(len(dimension)):
      self._dimension[i] = dimension[i]
        

  def setWidth(self, width):

    if not isinstance(width, (tuple, list, np.ndarray)):
      msg = 'Unable to set the width for Lattice ID={0} to {1} since ' \
            'it is not a Python tuple/list or NumPy ' \
            'array'.format(self._id, width)
      raise ValueError(msg)

    elif len(width) != 2 and len(width) != 3:
      msg = 'Unable to set the width for Lattice ID={0} to {1} since it does ' \
            'not contain 2 or 3 coordinates'.format(self._id, width)
      raise ValueError(msg)

    for dim in width:

      if not is_integer(dim) and not is_float(dim):
        msg = 'Unable to set the width for Lattice ID={0} to {1} since it is ' \
              'not an integer or floating point value'.format(self._id, dim)
        raise ValueError(msg)

      elif dim < 0:
        msg = 'Unable to set the width for Lattice ID={0} to {1} since it ' \
              'is a negative value'.format(self._id, dim)
        raise ValueError(msg)

    # Initialize width array to infinity by default
    self._width = np.zeros(3, dtype=np.float64)
    self._width[:] = MAX_FLOAT

    for i in range(len(width)):
      self._width[i] = width[i]


  def setUniverses(self, universes):

    if not isinstance(universes, (tuple, list, np.ndarray)):
      msg = 'Unable to set the universes for Lattice ID={0} to {1} ' \
            'since it is not a Python tuple/list or ' \
            'NumPy array'.format(self._id, universes)
      raise ValueError(msg)

    # if universes was input in 2D -> make 3D
    shape = np.shape(universes)
    if len(shape) == 2:
      universes = [universes]
    
    self._universes = np.asarray(universes, dtype=Universe)

    for i in range(self._dimension[0]):
      for j in range(self._dimension[1]):
        for k in range(self._dimension[2]):

          universe = self._universes[k][j][i]

          if self._width[0] != MAX_FLOAT:
            universe.setMaxX(self._width[0]/2.)
            universe.setMinX(-self._width[0]/2.)

          if self._width[1] != MAX_FLOAT:
            universe.setMaxY(self._width[1]/2.)
            universe.setMinY(-self._width[1]/2.)

          if self._width[2] != MAX_FLOAT:
            universe.setMaxZ(self._width[2]/2.)
            universe.setMinZ(-self._width[2]/2.)


  def computeVolumeFractions(self, volume=np.float64(1.), tolerance=1e-3):

    if not is_float(volume):
      msg = 'Unable to compute volume fractions for Lattice ID={0} since ' \
            'volume={1} is not a floating point value'.format(self._id, volume)
      raise ValueError(msg)

    volume_fraction = np.float64(1. / (self._dimension[0] * self._dimension[1] \
                                         * self._dimension[2]))

    for i in range(self._dimension[0]):
      for j in range(self._dimension[1]):
        for k in range(self._dimension[2]):
          universe = self._universes[k][j][i]
          universe.computeVolumeFractions(volume=(volume * volume_fraction),
                                          tolerance=tolerance)


  def initializeCellOffsets(self):

    # If we have already called this routine, return the total number of regions
    if not self._num_regions is None:
      return self._num_regions

    # Initialize an array for the cell offsets
    self._cell_offsets = np.zeros(tuple(self._dimension[::-1]), dtype=np.int64)

    # The cell offsets have not yet been initialized
    count = 0

    for i in range(self._dimension[0]):
      for j in range(self._dimension[1]):
        for k in range(self._dimension[2]):
          self._cell_offsets[k][j][i] = count
          self._universes[k][j][i].initializeCellOffsets()
          count += self._universes[k][j][i]._num_regions

    self._num_regions = count


  def buildNeighbors(self, depth=1):

    if self._universes is None:
      msg = 'Unable to build neighbor Universes for Lattice ID={0} since ' \
            'its Universes array has not yet been set'.format(self._id)
      raise ValueError(msg)

    elif not is_integer(depth) or depth <= 0:
      msg = 'Unable to build neighbor Universes for Lattice ID={0} to a ' \
            'depth of {1} which is not a positive integer'.format(depth)

    self._neighbor_depth = 1

    # Create 3D depth x depth x depth array of neighbor
    # Universes for each Lattice cell
    self._neighbor_universes = sliding_window(self._universes, 2*depth+1)

    # Iterate over each unique Universe in the Lattice and make recursive
    # call to build their neighbor Cells
    unique_universes = self.getUniqueUniverses()

    for universe_id, universe in unique_universes.items():
      universe.buildNeighbors()


  def findCell(self, localcoords):

    if not isinstance(localcoords, LocalCoords):
      msg = 'Unable to find cell in Lattice ID={0} since localcoords ' \
            '{1} is not a LocalCoords'.format(self._id, localcoords)
      raise ValueError(msg)

    # Compute the x and y indices for the Lattice cell this coord is in
    point = localcoords._point
    x = point._coords[0]
    y = point._coords[1]
    z = point._coords[2]

    # Compute the Lattice cell indices
    lat_x = math.floor((x + self._dimension[0]*self._width[0]*0.5 \
                          - self._offset[0]) / self._width[0])
    lat_y = math.floor((y + self._dimension[1]*self._width[1]*0.5 \
                          + self._offset[1]) / self._width[1])
    if self._dimension[2] == 1:
      lat_z = 0
    else:
      lat_z = math.floor((z + self._dimension[2]*self._width[2]*0.5 \
                          - self._offset[2]) / self._width[2])

    # Check if the LocalCoord is on the Lattice boundaries
    # If so adjust x or y Lattice cell indices

    # Compute the distance to the Lattice cell boundaries
    distance_x = math.fabs(math.fabs(x) - self._dimension[0]*self._width[0]*0.5\
                             - self._offset[0])
    distance_y = math.fabs(math.fabs(y) - self._dimension[1]*self._width[1]*0.5\
                             - self._offset[1])
    distance_z = math.fabs(math.fabs(z) - self._dimension[2]*self._width[2]*0.5\
                             - self._offset[2])

    if distance_x < ON_LATTICE_CELL_THRESH:
      if x > 0:
        lat_x = self._dimension[0] - 1
      else:
        lat_x = 0

    if distance_y < ON_LATTICE_CELL_THRESH:
      if y > 0:
        lat_y = self._dimension[1] - 1
      else:
        lat_y = 0

    if distance_z < ON_LATTICE_CELL_THRESH:
      if z > 0:
        lat_z = self._dimension[2] - 1
      else:
        lat_z = 0

    # Cast the Lattice indices as integers
    lat_x = int(lat_x)
    lat_y = int(lat_y)
    lat_z = int(lat_z)

    if (lat_x < 0 or lat_x >= self._dimension[0]) or \
          (lat_y < 0 or lat_y >= self._dimension[1]) or \
          (lat_z < 0 or lat_z >= self._dimension[2]):
      msg = 'Unable to find cell since the lattice indices ({0},{1},{2}) are ' \
          'outside of Lattice ID={3}'.format(lat_x, lat_y, lat_z, self._id)
      raise ValueError(msg)

    # Set the Lattice cell indices for the LocalCoords
    localcoords.setLatticeX(lat_x)
    localcoords.setLatticeY(lat_y)
    localcoords.setLatticeZ(lat_z)

    # Compute local position of Point in the next level Universe
    next_x = x - (-self._dimension[0]*self._width[0]*0.5 + self._offset[0] \
                    + (lat_x + 0.5) * self._width[0])
    next_y = y - (-self._dimension[1]*self._width[1]*0.5 + self._offset[1] \
                    + (lat_y + 0.5) * self._width[1])
    next_z = z - (-self._dimension[2]*self._width[2]*0.5 + self._offset[2] \
                    + (lat_z + 0.5) * self._width[2])
    next_point = Point(x=next_x, y=next_y, z=next_z)

    universe = self._universes[lat_z][lat_y][lat_x]

    if isinstance(universe, Universe):
      next_coords = UnivCoords(next_point)
      next_coords.setUniverse(universe)

    elif isinstance(universe, Lattice):
      next_coords = LatCoords(next_point)
      next_coords.setLattice(universe)

    else:

      msg = 'Unable to find Cell since Lattice ID={0} does ' \
          'not contain a Universe or Lattice in lattice cell ' \
          '({1}, {2}, {3})'.format(self._id, lat_x, lat_y, lat_z)
      raise ValueError(msg)

    localcoords.setNext(next_coords)
    next_coords.setPrev(localcoords)

    return universe.findCell(next_coords)


  def findRegion(self, region_id, lat_coords):

    if not is_integer(region_id):
      msg = 'Unable to find region_id={0} in Lattice ID={1} since {2} is ' \
            'not an integer value'.format(region_id, self._id, region_id)
      raise ValueError(msg)

    elif not isinstance(lat_coords, LatCoords):
      msg = 'Unable to find region_id={0} in Lattice ID={1} since {2} is ' \
            'not a LatCoords'.format(region_id, self._id, lat_coords)
      raise ValueError(msg)

    # Find Lattice cell indices where region is less than the Cell offset
    indices = np.where(np.swapaxes(self._cell_offsets, 0, 2) <= region_id)

    if indices != None:
      lat_z = indices[2][-1]
      lat_y = indices[1][-1]
      lat_x = indices[0][-1]
      universe = self._universes[lat_z][lat_y][lat_x]
      offset = self._cell_offsets[lat_z][lat_y][lat_x]

    else:
      msg = 'Unable to find region_id={0} for FSR ID={1} in Lattice ' \
            'ID={2}'.format(region_id, self._id)
      raise ValueError(msg)

    region_id -= offset
    lat_coords.setLatticeX(lat_x)
    lat_coords.setLatticeY(lat_y)
    lat_coords.setLatticeZ(lat_z)

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


  def toString(self):
    string = self.__repr__()
    unique_universes = self.getUniqueUniverses()
    for universe_id, universe in unique_universes.items():
      string += universe.toString()
    return string


  def __repr__(self):

    string = 'Lattice\n'
    string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self._id)
    string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self._name)
    string += '{0: <16}{1}{2}\n'.format('\tType', '=\t', self._type)
    string += '{0: <16}{1}{2}\n'.format('\tDimension', '=\t', self._dimension)
    string += '{0: <16}{1}{2}\n'.format('\tOffset', '=\t', self._offset)
    string += '{0: <16}{1}{2}\n'.format('\tWidth', '=\t', self._width)
    string += '{0: <16}{1}{2}\n'.format('\t# Regions', '=\t', self._num_regions)

    string += '{0: <16}{1}'.format('\tUniverses', '\n')

    for i in range(self._dimension[0]):
      string += '\t'

      for j in range(self._dimension[1]):

          for k in range(self._dimension[2]):
            string += '{0} '.format(self._universes[k][j][i]._id)

          string += '\n'

      string += '\n'

    return string



################################################################################
####################################  Cell  ####################################
################################################################################


# Lists of all IDs for all Cells created
CELL_IDS = list()

# A static variable for auto-generated Cell IDs
AUTO_CELL_ID = 10000

def reset_auto_cell_id():
  global AUTO_CELL_ID, CELL_IDS
  AUTO_CELL_ID = 10000
  CELL_IDS = list()


class Cell(object):

  def __init__(self, cell_id=None, name='', fill=None, rot=None):

    # Initialize Cell class attributes
    self._id = None
    self._name = None
    self._fill = None
    self._rot = None
    self._type = None
    self._num_subcells = None
    self._volume_fraction = np.float64(0.)
    self._volume = np.float64(0.)

    self._neighbor_cells = np.empty(shape=(0,), dtype=Cell)

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

    self.setMaxX(MAX_FLOAT)
    self.setMaxY(MAX_FLOAT)
    self.setMaxZ(MAX_FLOAT)
    self.setMinX(MIN_FLOAT)
    self.setMinY(MIN_FLOAT)
    self.setMinZ(MIN_FLOAT)

    if not fill is None:
      self.setFill(fill)

    if not rot is None:
      self.setRotation(rot)


  def __deepcopy__(self, memo):

    existing = memo.get(id(self))

    # If this is the first time we have tried to copy this object, create a copy
    if existing is None:

      clone = type(self).__new__(type(self))
      clone._id = self._id
      clone._name = self._name
      clone._fill = copy.deepcopy(self._fill, memo)
      clone._rot = self._rot
      clone._type = self._type
      clone._num_subcells = self._num_subcells
      clone._volume_fraction = self._volume_fraction
      clone._volume = self._volume
      clone._neighbor_cells = copy.deepcopy(self._neighbor_cells, memo)

      clone._surfaces = dict()
      for surface_id in self._surfaces.keys():
        surface = self._surfaces[surface_id][0]
        halfspace =self._surfaces[surface_id][1]
        clone_surface = copy.deepcopy(surface)
        clone._surfaces[surface_id] = (clone_surface, halfspace)

      clone._max_x = self._max_x
      clone._min_x = self._min_x
      clone._max_y = self._max_y
      clone._min_y = self._min_y
      clone._max_z = self._max_z
      clone._min_z = self._min_z

      memo[id(self)] = clone

      return clone

    # If this object has been copied before, return the first copy made
    else:
      return existing


  def __gt__(self, other):
    return (id(self) > id(other))


  def __ge__(self, other):
    return (id(self) >= id(other))


  def __lt__(self, other):
    return (id(self) < id(other))


  def __le__(self, other):
    return (id(self) <= id(other))


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

    cells = dict()

    if self._type == 'universe' or self._type == 'lattice':
      cells = self._fill.getAllCells()

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


  def getNeighbors(self):
    return self._neighbor_cells


  def getUniqueNeighbors(self):

    # Select only unique Cells and return them as a tuple
    unique_neighbors = np.unique(self._neighbor_cells)
    return unique_neighbors


  def getNeighborsHash(self):
    neighbor_cells = np.copy(self.getNeighbors())
    neighbor_cells.sort()
    raw_data = neighbor_cells.view()
    neighbors_hash = int(sha1(raw_data).hexdigest(), 16)
    return neighbors_hash


  def getUniqueNeighborsHash(self):
    neighbor_cells = np.copy(self.getUniqueNeighbors())
    neighbor_cells.sort()
    raw_data = neighbor_cells.view()
    neighbors_hash = int(sha1(raw_data).hexdigest(), 16)
    return neighbors_hash


  def setId(self, cell_id=None):

    global CELL_IDS

    if cell_id is None:
      global AUTO_CELL_ID
      self._id = AUTO_CELL_ID
      CELL_IDS.append(AUTO_CELL_ID)
      AUTO_CELL_ID += 1

    # Check that the ID is an integer and wasn't already used
    elif is_integer(cell_id):

      # If the Cell already has an ID, remove it from global list
      if not self._id is None:
        CELL_IDS.remove(self._id)

      if cell_id in CELL_IDS:
        msg = 'Unable to set Cell ID to {0} since a Cell with this ID was ' \
              'already initialized'.format(cell_id)
        raise ValueError(msg)

      if cell_id < 0:
        msg = 'Unable to set Cell ID to {0} since it must be a ' \
              'non-negative integer'.format(cell_id)
        raise ValueError(msg)

      else:
        self._id = cell_id
        CELL_IDS.append(cell_id)

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

  def setRotation(self, rot):

    # TODO: error checking

    self._rot = rot

  def setType(self, type):

    if not is_string(type):
      msg = 'Unable to set the type for Cell ID={0} to {1} since it is not ' \
            'a string'.format(self._id, type)
      raise ValueError(msg)

    elif not type.lower() in ['universe', 'lattice', 'material']:
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

    elif not halfspace in [-1, +1]:
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

          match = True
          for key in coeffs.keys():
            coeff = coeffs[key]
            test_coeff = test_coeffs[key]
            if abs(coeff-test_coeff) < 1e-10:
              match = False
              break
          if match:
            warnings.warn('Skipping redundant surface <{0}> ' + \
                          'for cell <{1}>: ' +\
                          'matches <{2}>'.format(surface._name,
                                                 self._name,
                                                 test_surface._name))
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

    elif not isinstance(halfspaces, (list, tuple, np.ndarray)):
      msg = 'Unable to add Surfaces to Cell ID={0} since {1} is not a Python ' \
            'tuple/list or NumPy array'.format(self._id, halfspaces)
      raise ValueError(msg)

    elif len(surfaces) != len(halfspaces):
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


  def addNeighborCell(self, cell):

    if not isinstance(cell, Cell):
      msg = 'Unable to add a neighbor Cell to Cell ID={0} ' \
            'since it {1} is not a Cell object'.format(self._id, cell)
      raise ValueError(msg)

    self._neighbor_cells = np.append(self._neighbor_cells, cell)


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

    self.setMaxX(MAX_FLOAT)
    self.setMaxY(MAX_FLOAT)
    self.setMaxZ(MAX_FLOAT)
    self.setMinX(MIN_FLOAT)
    self.setMinY(MIN_FLOAT)
    self.setMinZ(MIN_FLOAT)

    for surface_id in self._surfaces:
      surface = self._surfaces[surface_id][0]
      halfspace = self._surfaces[surface_id][1]

      max_x = surface.getMaxX(halfspace=halfspace)
      max_y = surface.getMaxY(halfspace=halfspace)
      max_z = surface.getMaxZ(halfspace=halfspace)

      min_x = surface.getMinX(halfspace=halfspace)
      min_y = surface.getMinY(halfspace=halfspace)
      min_z = surface.getMinZ(halfspace=halfspace)

      if max_x != MAX_FLOAT and max_x < self._max_x:
        self.setMaxX(max_x)
      if max_y != MAX_FLOAT and max_y < self._max_y:
        self.setMaxY(max_y)
      if max_z != MAX_FLOAT and max_z < self._max_z:
        self.setMaxZ(max_z)

      if min_x != MIN_FLOAT and min_x > self._min_x:
        self.setMinX(min_x)
      if min_y != MIN_FLOAT and min_y > self._min_y:
        self.setMinY(min_y)
      if min_z != MIN_FLOAT and min_z > self._min_z:
        self.setMinZ(min_z)

    # If we could not find a bounds for any dimension, readjust
    # it to +/- infinity
    if self._max_x == MIN_FLOAT:
      self.setMaxX(MAX_FLOAT)
    if self._max_y == MIN_FLOAT:
      self.setMaxY(MAX_FLOAT)
    if self._max_z == MIN_FLOAT:
      self.setMaxZ(MAX_FLOAT)

    if self._min_x == MAX_FLOAT:
      self.setMinX(MIN_FLOAT)
    if self._min_y == MAX_FLOAT:
      self.setMinY(MIN_FLOAT)
    if self._min_z == MAX_FLOAT:
      self.setMinZ(MIN_FLOAT)


  def buildNeighbors(self):

    # Loop over all of the Surfaces
    for surface_id in self._surfaces.keys():

      # Add this Cell to the Surface's neighbor Cells
      surface = self._surfaces[surface_id][0]
      halfspace = self._surfaces[surface_id][1]
      surface.addNeighborCell(self, halfspace)

    # Make recursive call to the Cell's Fill
    if self._type in ['universe', 'lattice']:
      self._fill.buildNeighbors()


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

    elif not is_integer(num_samples):
      msg = 'Unable to compute the volume fraction for Cell ' \
            'ID={0} since num_samples={1} is not an integer ' \
            'value'.format(self._id, num_samples)
      raise ValueError(msg)

    elif not is_float(tolerance):
      msg = 'Unable to compute the volume fraction for Cell ' \
            'ID={0} since tolerance={1} is not a floating point ' \
            'value'.format(self._id, tolerance)
      raise ValueError(msg)

    from numpy.random import uniform

    # Initialize the point
    point = Point()

    # Compute the volume/area of the bounding box we sample from
    box_volume = np.float64(1.)

    if self._min_x > MIN_FLOAT and self._max_x < MAX_FLOAT:
      box_volume *= (self._max_x - self._min_x)
    if self._min_y > MIN_FLOAT and self._max_y < MAX_FLOAT:
      box_volume *= (self._max_y - self._min_y)
    if self._min_z > MIN_FLOAT and self._max_z < MAX_FLOAT:
      box_volume *= (self._max_z - self._min_z)

    # Initialize variables
    counter = 0.
    tot_samples = 0.
    uncertainty = MAX_FLOAT

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
    clone.setType(self._type)

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
    string = self.__repr__()
    if not isinstance(self._fill, Material):
      string += self._fill.toString()
    return string

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


  def __deepcopy__(self, memo):

    existing = memo.get(id(self))

    # If this is the first time we have tried to copy this object, create a copy
    if existing is None:

      clone = type(self).__new__(type(self))
      clone._point = copy.deepcopy(self._point, memo)
      clone._type = self._type
      clone._next = copy.deepcopy(self._next, memo)
      clone._prev = copy.deepcopy(self._prev, memo)

      memo[id(self)] = clone

      return clone

    # If this object has been copied before, return the first copy made
    else:
      return existing


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


  def getNeighbors(self):

    msg = 'The abstract LocalCoords class does not have the ' \
          'getNeighbors() method implemented. Try again wih UnivCoords ' \
          'and LatCoords subclasses.'
    raise ValueError(msg)


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


  def __deepcopy__(self, memo):

    existing = memo.get(id(self))

    # If this is the first time we have tried to copy this object, create a copy
    if existing is None:

      clone = super(UnivCoords, self).__deepcopy__(self, memo)
      clone._universe = copy.deepcopy(self._universe, memo)
      clone._cell = copy.deepcopy(self._cell, memo)

      memo[id(self)] = clone

      return clone

    # If this object has been copied before, return the first copy made
    else:
      return existing


  def getNeighbors(self, neighbors=None):

    if neighbors is None:
      neighbors = list()

    cells = self._cell.getNeighbors()
    neighbors.append(cells)

    # Make recursive call to next LocalCoords or return
    if self._next is None:
      return neighbors
    else:
      return self._next.getNeighbors(neighbors=neighbors)


  def getUniqueNeighbors(self, neighbors=None):

    if neighbors is None:
      neighbors = list()

    cells = self._cell.getUniqueNeighbors()
    neighbors.append(cells)

    # Make recursive call to next LocalCoords or return
    if self._next is None:
      return tuple(neighbors)
    else:
      return self._next.getUniqueNeighbors(neighbors)


  def getNeighborsHash(self, neighbors=None):

    if neighbors is None:
      neighbors = list()

    cells = self._cell.getNeighborsHash()
    neighbors.append(cells)

    # Make recursive call to next LocalCoords or return
    if self._next is None:
      return hash(tuple(neighbors))
    else:
      return self._next.getNeighborsHash(neighbors=neighbors)


  def getUniqueNeighborsHash(self, neighbors=None):

    if neighbors is None:
      neighbors = list()

    cells = self._cell.getUniqueNeighborsHash()
    neighbors.append(cells)

    # Make recursive call to next LocalCoords or return
    if self._next is None:
      return hash(tuple(neighbors))
    else:
      return self._next.getUniqueNeighborsHash(neighbors=neighbors)


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

  def __init__(self, point=None, next=None, prev=None, lattice=None,
               lat_x=None, lat_y=None, lat_z=None):

    super(LatCoords, self).__init__(point, next, prev)

    self._type = 'lattice'
    self._lattice = None
    self._lat_x = None
    self._lat_y = None
    self._lat_z = None

    if not lattice is None:
      self.setLattice(lattice)

    if not lat_x is None:
      self.setLatticeX(lat_x)

    if not lat_y is None:
      self.setLatticeY(lat_y)

    if not lat_z is None:
      self.setLatticeY(lat_z)


  def __deepcopy__(self, memo):

    existing = memo.get(id(self))

    # If this is the first time we have tried to copy this object, create a copy
    if existing is None:

      clone = super(LatCoords, self).__deepcopy__(self, memo)
      clone._lattice = copy.deepcopy(self._lattice, memo)
      clone._lat_x = self._lat_x
      clone._lat_y = self._lat_y
      clone._lat_z = self._lat_z

      memo[id(self)] = clone

      return clone

    # If this object has been copied before, return the first copy made
    else:
      return existing


  def getNeighbors(self, neighbors=None):

    if neighbors is None:
      neighbors = list()

    # Append the Universes in the neighboring Lattice cells to the
    # neighbors list
    universes = self._lattice.getNeighbors(self._lat_x, self._lat_y, self._lat_z)
    neighbors.append(universes)

    # Make recursive call to next LocalCoords if it exists or return
    if self._next is None:
      return neighbors
    else:
      return self._next.getNeighbors(neighbors=neighbors)


  def getUniqueNeighbors(self, neighbors=None):

    if neighbors is None:
      neighbors = list()

    # Append the Universes in the neighboring Lattice cells to the
    # neighbors list
    universes = self._lattice.getUniqueNeighbors(self._lat_x,
                                                 self._lat_y, self._lat_z)
    neighbors.append(universes)

    # Make recursive call to next LocalCoords if it exists or return
    if self._next is None:
      return neighbors
    else:
      return self._next.getUniqueNeighbors(neighbors)


  def getNeighborsHash(self, neighbors=None):

    if neighbors is None:
      neighbors = list()

    # Append the Universes in the neighboring Lattice cells to the
    # neighbors list
    universes = self._lattice.getNeighborsHash(self._lat_x,
                                               self._lat_y, self._lat_z)
    neighbors.append(universes)

    # Make recursive call to next LocalCoords if it exists or return
    if self._next is None:
      return hash(tuple(neighbors))
    else:
      return self._next.getNeighborsHash(neighbors=neighbors)


  def getUniqueNeighborsHash(self, neighbors=None):

    if neighbors is None:
      neighbors = list()

    # Append the Universes in the neighboring Lattice cells to the
    # neighbors list
    universes = self._lattice.getUniqueNeighborsHash(self._lat_x,
                                                     self._lat_y, self._lat_z)
    neighbors.append(universes)

    # Make recursive call to next LocalCoords if it exists or return
    if self._next is None:
      return hash(tuple(neighbors))
    else:
      return self._next.getUniqueNeighborsHash(neighbors)


  def setLattice(self, lattice):

    if not isinstance(lattice, Lattice):
      msg = 'Unable to set the Lattice to {0} for LocalCoords since it ' \
            'is not a Lattice'.format(lattice)
      raise ValueError(msg)

    self._lattice = lattice


  def setLatticeX(self, lat_x):

    if not is_integer(lat_x):
      msg = 'Unable to set the Lattice X to {0} for LocalCoords since it ' \
            'is not an integer'.format(lat_x)
      raise ValueError(msg)

    self._lat_x = lat_x


  def setLatticeY(self, lat_y):

    if not is_integer(lat_y):
      msg = 'Unable to set the Lattice Y to {0} for LocalCoords since it ' \
            'is not an integer'.format(lat_y)
      raise ValueError(msg)

    self._lat_y = lat_y


  def setLatticeZ(self, lat_z):

    if not is_integer(lat_z):
      msg = 'Unable to set the Lattice Z to {0} for LocalCoords since it ' \
            'is not an integer'.format(lat_z)
      raise ValueError(msg)

    self._lat_z = lat_z


  def __repr__(self):

    string = super(LatCoords, self).__repr__()

    string += '{0: <16}{1}{2}\n'.format('\tLattice', '=\t', self._lattice._id)
    string += '{0: <16}{1}{2}\n'.format('\tLattice X', '=\t', self._lat_x)
    string += '{0: <16}{1}{2}\n'.format('\tLattice Y', '=\t', self._lat_y)
    string += '{0: <16}{1}{2}\n'.format('\tLattice ', '=\t', self._lat_z)

    return string



def sliding_window(arr, window_size):
  """ Construct a sliding window view of the array"""

    # This algorithm is based on a technique for building sliding windows on
    # NumPy arrays using the as_strided(...) routine described here:
    # http://stackoverflow.com/questions/10996769/pixel-neighbors-in-2d-array-image-using-python

  # Correct input parameters if needed
  arr = np.asarray(arr)
  window_size = int(window_size)

  # Determine the shape of the sliding windows array
  shape = (arr.shape[0] - window_size + 1,
           arr.shape[1] - window_size + 1,
           arr.shape[2] - window_size + 1,
           window_size, window_size, window_size)

  # Correct the shape of the sliding windows array to account for dimensions
  # of the original array that are not as large as the window size
  if shape[0] <= 0:
    shape = (1, shape[1], shape[2], arr.shape[0], shape[4], shape[5])
  if shape[1] <= 0:
    shape = (shape[0], 1, shape[2], shape[3], arr.shape[1], shape[5])
  if shape[2] <= 0:
    shape = (shape[0], shape[1], 1, shape[3], shape[4], arr.shape[2])

  # Determine the array of strides for NumPy for the sliding window array
  strides = (arr.shape[1] * arr.shape[2] * arr.itemsize,
             arr.shape[2] * arr.itemsize, arr.itemsize,
             arr.shape[1] * arr.shape[2] * arr.itemsize,
             arr.shape[2] * arr.itemsize, arr.itemsize)

  return as_strided(arr, shape=shape, strides=strides)
