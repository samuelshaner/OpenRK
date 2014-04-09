__author__ = 'Will Boyd'
__email__ = 'wboyd@mit.edu'


from material import Material
from surface import Surface, on_surface_thresh
from localcoords import *
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
    self._num_regions = None

    self.setId(universe_id)
    self.setName(name)

  def getId(self):
    return self._id


  def getName(self):
    return self._name


  def getCells(self):
    return self._cells


  def getCellOffsets(self):
    return self._cell_offsets


  def getNumRegions(self):
    return self._num_regions


  def setId(self, universe_id=None):

    global universe_ids

    if universe_id is None:
      global auto_universe_id
      self._id = auto_universe_id
      universe_ids.append(auto_universe_id)
      auto_universe_id += 1

    # Check that the ID is an integer and wasn't already used
    elif not is_integer(universe_id):

      # If the Cell already has an ID, remove it from global list
      if not self._id is None:
        universe_ids.remove(self._id)

      if universe_id in universe_ids:
        exit('Unable to set Universe ID to %s since a Universe with this ID '
             'was already initialized.', str(universe_id))

      if universe_id < 0:
        exit('Unable to set Univeres ID to %d since it must be a '
             'non-negative integer', universe_id)

      else:
        self._id = universe_id
        universe_ids.append(universe_id)

    else:
      exit('Unable to set Universe ID to a non-integer %s', str(universe_id))


  def setName(self, name):

    if not is_string(name):
      exit('Unable to set name for Universe ID=%d with a non-string value %s',
           self._id, str(name))

    else:
      self._name = name


  def addCell(self, cell):

    if not issubclass(cell, Cell):
      exit('Unable to add a Cell to Universe ID=%d since %s is not a Cell',
            self._id, str(cell))

    cell_id = cell.getId()

    if not cell_id in self._cells.keys():
      self._cells[cell_id] = cell


  def addCells(self, cells):

    if not isinstance(cells, (list, tuple, np.ndarray)):
      exit('Unable to add Cells to Universe ID=%d since %s is not a Python '
           'tuple/list or NumPy array', self._id, str(cells))

    for i in range(len(cells)):
      self.addSurface(cells[i])


  def removeCell(self, cell):

    if not issubclass(cell, Cell):
      exit('Unable to remove a Cell from Universe ID=%d since %s is not a '
           'Cell', self._id, str(cell))

    cell_id = cell.getId()

    if cell_id in self._cells.keys():
      del self._cells[cell_id]


  def initializeCellOffsets(self):


    # If we have already called this routine, return the total number of regions
    if not self._num_regions is None:
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
           'a LocalCoords', self._id, str(localcoords))

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

          if isinstance(fill, Universe):
            next_coords = UnivCoords(localcoords.getPoint())
            next_coords.setUniverse(fill)

          elif isinstance(fill, Lattice):
            next_coords = LatCoords(localcoords.getPoint())
            next_coords.setLattice(fill)

          else:
            exit('Unable to find cell since in Universe ID=%d since Cell ID=%d '
                 'is not filled by a Material, Universe or Lattice',
                 self._id, cell.getId())

          localcoords.setNext(next_coords)
          next_coords.setPrev(localcoords)

          return fill.findCell(next_coords)


#  def findCell(self, region_id):


  def toString(self):

    string = ''

    string += 'Universe\n'

    universe_id = '{0: <16}'.format('\tID') + '=\t' + str(self._id)
    string += universe_id + '\n'

    name = '{0: <16}'.format('\tName') + '=\t' + self._name
    string += name + '\n'

    num_regions = '{0: <16}'.format('\t# Regions') + '=\t' + self._num_regions
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


  def getOutside(self):
    return self._outside


  def getUniverses(self):
    return self._universes


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
             'with this ID was already initialized.', str(lattice_id))

      if lattice_id < 0:
        exit('Unable to set Lattice ID to %d since it must be a '
             'non-negative integer', lattice_id)

      else:
        self._id = lattice_id
        universe_ids.append(lattice_id)

    else:
      exit('Unable to set a non-integer Lattice ID %s', str(lattice_id))


  def setName(self, name):

    if not is_string(name):
      exit('Unable to set name for Lattice ID=%d with a non-string value %s',
           self._id, str(name))

    else:
      self._name = name


  def setType(self, type):

    if not is_string(type):
      exit('Unable to set the type for Lattice ID=%d with a non-string '
           'value %s', self._id, str(type))

    if not type in ['rectangular']:
      exit('Unable to set the type for Lattice ID=%d to %s since it is not '
           'rectangular', self._id, type)

    self._type = type


  def setDimension(self, dimension):

    if not isinstance(dimension, (tuple, list, np.ndarray)):
      exit('Unable to set Lattice ID=%d dimension to %s since it is not '
           'a Python tuple/list or NumPy array', self._id, str(dimension))

    if len(dimension) != 2 and len(dimension) != 3:
      exit('Unable to set Lattice ID=%d dimension to %s since it does '
           'not contain 2 or 3 coordinates', self._id, str(dimension))

    for dim in dimension:

      if not isinstance(dim, (int, np.int32, np.int64)):
        exit('Unable to set the dimension for Lattice ID=%d to %s since it '
             'is not an integer', self._id, str(dim))

      if dim < 0:
        exit('Unable to set Lattice ID=%d dimension to %s since it '
             'is a negative value', self._id, dim)

    self._dimension = np.zeros(len(dimension), dtype=np.int64)

    for i in range(len(dimension)):
      self._dimension[i] = dimension[i]


  def setLowerLeft(self, lower_left):

    if not isinstance(lower_left, (tuple, list, np.ndarray)):
      exit('Unable to set the lower_left for Lattice ID=%d to %s since it is '
           'not a Python tuple/list or NumPy array', self._id, str(lower_left))

    if len(lower_left) != 2 and len(lower_left) != 3:
      exit('Unable to set the lower_left for Lattice ID=%d to %s since it does '
           'not contain 2 or 3 coordinates', self._id, str(lower_left))

    for dim in lower_left:

      if not isinstance(dim, (int, np.int32, np.int64)) \
        and not isinstance(dim, (float, np.float32, np.float64)):
        exit('Unable to set the lower_left for Lattice ID=%d to %s since it '
             'is not an integer or floating point value', self._id, str(dim))

    self._lower_left = np.zeros(len(lower_left), dtype=np.float64)

    for i in range(len(lower_left)):
      self._lower_left[i] = lower_left[i]


  def setWidth(self, width):

    if not isinstance(width, (tuple, list, np.ndarray)):
      exit('Unable to set the width for Lattice ID=%d to %s since it is not '
           'a Python tuple/list or NumPy array', self._id, str(width))

    if len(width) != 2 and len(width) != 3:
      exit('Unable to set the width for Lattice ID=%d to %s since it does '
           'not contain 2 or 3 coordinates', self._id, str(width))

    for dim in width:

      if not isinstance(dim, (int, np.int32, np.int64)) \
        and not isinstance(dim, (float, np.float32, np.float64)):
        exit('Unable to set the width for Lattice ID=%d to %s since it is '
             'not an integer or floating point value', self._id, str(dim))

      if dim < 0:
        exit('Unable to set the width for Lattice ID=%d to %s since it '
             'is a negative value', self._id, dim)

    self._width = width

    self._width = np.zeros(len(width), dtype=np.float64)

    for i in range(len(width)):
      self._with[i] = width[i]


  def setUniverses(self, universes):

    if not isinstance(universes, (tuple, list, np.ndarray)):
      exit('Unable to set the universes for Lattice ID=%d to %s since it is '
           'not a Python tuple/list or NumPy array', self._id, str(universes))

    self._universes = universes


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
          self._cell_offsets[i,j] = count
          self._universes[i,j].initializeCellOffsets()
          count += self._universes[i,j].getNumRegions()

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
           'a LocalCoords', self._id, str(localcoords))

    # Compute the x and y indices for the Lattice cell this coord is in
    point = localcoords.getPoint()
    x = point.getX()
    y = point.getY()

    # Compute the Lattice cell indices
    lat_x = math.floor((x - self._lower_left[0]) / self._width[0])
    lat_y = math.floor((y - self._lower_left[1]) / self._width[1])

    # If the indices are outside the bound of the Lattice
    if (lat_x < 0 or lat_x >= self._dimension[0]) or \
      (lat_y < 0 or lat_y >= self._dimension[0]):
      exit('Unable to find cell since the lattice indices (%d,%d) are '
           'outside of Lattice ID=%d', lat_x, lat_y, self._id)

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
      exit('Unable to find cell since in Lattice ID=%d does not contain a '
           'Universe or Lattice in lattice cell (%d, %d)',
           self._id, lat_x, lat_y)

    localcoords.setNext(next_coords)
    next_coords.setPrev(localcoords)

    return universe.findCell(next_coords)



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

    outside = '{0: <16}'.format('\tOutside') + '=\t'
    outside += str(self._outside)
    string += outside + '\n'

    num_regions = '{0: <16}'.format('\t# Regions') + '=\t' + self._num_regions
    string += num_regions + '\n'

    universes = '{0: <16}'.format('\tUniverses') + '\n'

    for i in range(len(self._universes)):
      universes += '\t'

      for j in range(len(self._universes[0])):
        universes += '%s ' % str(int(self._universes[i][j]))

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

  def __init__(self, cell_id=None, name='', universe=None, fill=None):

    # Initialize Cell class attributes
    self._id = None
    self._name = None
    self._fill = None
    self._type = None
    self._num_subcells = None

    # Keys   - Surface IDs
    # Values - (halfpsace, Surface) tuples
    self._surfaces = dict()

    self.setId(cell_id)
    self.setName(name)

    if not universe is None:
      self.setUniverse(universe)

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


  def getSurfaces(self):
    return self._surfaces


  def setId(self, cell_id=None):

    global cell_ids

    if cell_id is None:
      global auto_cell_id
      self._id = auto_cell_id
      cell_ids.append(auto_cell_id)
      auto_cell_id += 1

    # Check that the ID is an integer and wasn't already used
    elif not is_integer(cell_id):

      # If the Cell already has an ID, remove it from global list
      if not self._id is None:
        cell_ids.remove(self._id)

      if cell_id in cell_ids:
        exit('Unable to set Cell ID to %s since a Cell with this ID was '
             'already initialized.', str(cell_id))

      if cell_id < 0:
        exit('Unable to set Cell ID to %d since it must be a '
             'non-negative integer', cell_id)

      else:
        self._id = cell_id
        cell_ids.append(cell_id)

    else:
      exit('Unable to set Cell ID to a non-integer %s', str(cell_id))


  def setName(self, name):

    if not is_string(name):
      exit('Unable to set name for Cell ID=%d with a non-string value %s',
           self._id, str(name))

    else:
      self._name = name


  def setFill(self, fill):

    if isinstance(fill, Universe):
      self.setType('universe')
    elif isinstance(fill, Lattice):
      self.setType('lattice')
    elif isinstance(fill, Material):
      self.setType('material')
    else:
      exit('Unable to set fill for Cell ID=%d to %s since it is not a '
           'Universe, Lattice or a Material', self._id, str(fill))

    self._fill = fill


  def setType(self, type):

    if not is_string(type):
      exit('Unable to set the type for Cell ID=%d to %s since it is not '
           'a string', self._id, str(type))

    if not type.lower() in ['universe', 'lattice', 'material']:
      exit('Unable to set the type for Cell ID=%d to %s since it is not '
           'universe, lattice or material', self._id, type)

    self._type = type.lower()


  def addSurface(self, surface, halfspace):

    if not issubclass(surface, Surface):
      exit('Unable to add a Surface to Cell ID=%d since %s is not a Surface',
            self._id, str(surface))

    if not halfspace in [-1, +1]:
      exit('Unable to add a Surface to Cell ID=%d with halfspace %s since '
           'it is not +/-1', self._id, str(halfspace))

    surf_id = surface.getId()

    if not surf_id in self._surfaces.keys():
      self._surfaces[surf_id] = (surface, halfspace)


  def addSurfaces(self, surfaces, halfspaces):

    if not isinstance(surfaces, (list, tuple, np.ndarray)):
      exit('Unable to add Surfaces to Cell ID=%d since %s is not a Python '
           'tuple/list or NumPy array', self._id, str(surfaces))

    if not isinstance(halfspaces, (list, tuple, np.ndarray)):
      exit('Unable to add Surfaces to Cell ID=%d since %s is not a Python '
           'tuple/list or NumPy array', self._id, str(halfspaces))

    if len(surfaces) != len(halfspaces):
      exit('Unable to add Surfaces to Cell ID=%d since the number of '
           'Surfaces (%d) and halfspaces (%d) are not equal',
           self._id, len(surfaces), len(halfspaces))

    for i in range(len(surfaces)):
      self.addSurface(surfaces[i], halfspaces[i])


  def removeSurface(self, surface):

    if not issubclass(surface, Surface):
      exit('Unable to remove a surface from Cell ID=%d since %s is not a '
           'Surface', self._id, str(surface))

    surf_id = surface.getId()

    if surf_id in self._surfaces.keys():
      del self._surfaces[surf_id]


  def getNumSubCells(self):

    # If we have already called this routine, return the number of subcells
    if not self._num_subcells is None:
      return self._num_subcells

    # The cell offsets have not yet been initialized - we must compute them
    elif isinstance(self._fill, Material):
      self._num_subcells = 1

    elif isinstance(self._fill, (Universe, Lattice)):
      self._fill.initializeCellOffsets()
      self._num_subcells = self._fill.getNumRgions()

    else:
      exit('Unable to compute the number of subcells for Cell ID=%d since '
           'it is not filled by a Material, Universe or Lattice', self._id)

    return self._num_subcells


  def containsPoint(self, point):

    if not isinstance(point, Point):
      exit('Unable to determine if point is in Cell ID=%d since %s is not '
           'a Point', self._id, str(point))

    for surface_id in self._surfaces:

      halfspace = self._surfaces[surface_id][0]
      surface = self._surfaces[surface_id][1]

      # Return false if the Point is not in the correct Surface halfspace
      if (surface.evaluate(point) * halfspace) < -on_surface_thresh:
        return False

    # Return true if the Point is in the correct halfspace for each Surface
    return True


  def toString(self):

    string = ''

    string += 'Cell\n'

    cell_id = '{0: <16}'.format('\tID') + '=\t' + str(self._id)
    string += cell_id + '\n'

    name = '{0: <16}'.format('\tName') + '=\t' + self._name
    string += name + '\n'

    fill = '{0: <16}'.format('\tFill') + '=\t'
    fill += str(self._fill)
    string += fill + '\n'

    type = '{0: <16}'.format('\tType') + '=\t'
    type += str(self._type)
    string += type + '\n'

    num_subcells = '{0: <16}'.format('\t# Regions') + '=\t' + self._num_subcells
    string += num_subcells + '\n'

    surfaces = '{0: <16}'.format('\tSurfaces') + '=\t'
    surfaces += str(self._surfaces.keys())
    string += surfaces + '\n'

    return string


  def printString(self):
    print(self.toString())