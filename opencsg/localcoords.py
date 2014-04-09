__author__ = 'Will Boyd'
__email__ = 'wboyd@mit.edu'


from universe import *
from checkvalue import *


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
           'a Point object', str(point))

    self._point = point


  def setNext(self, next):

    if not isinstance(next, None) and not issubclass(next, LocalCoords):
      exit('Unable to set the next to %s for LocalCoords since it is not '
           'a LocalCoords object', str(next))

    self._next = next


  def setPrev(self, prev):

    if not isinstance(prev, None) and not issubclass(prev, LocalCoords):
      exit('Unable to set the prev to %s for LocalCoords since it is not '
           'a LocalCoords object', str(prev))

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

    point = '{0: <16}'.format('\tPoint') + '=\t' + str(self._point.getCoords())
    string += point + '\n'

    return string


  def printString(self):
    print(self.toString())



class UnivCoords(LocalCoords):

  def __init__(self, point=None, next=None, prev=None,
               universe=None, cell=None):

    super(self, point, next, prev)

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
           'is not a Universe', str(universe))

    self._universe = universe


  def setCell(self, cell):

    if not isinstance(cell, Cell):
      exit('Unable to set the Cell to %s for LocalCoords since it '
           'is not a Cell', str(cell))

    self._cell = cell


  def toString(self):

    string = super(self)

    universe_id = '{0: <16}'.format('\tUniverse') + '=\t'
    universe_id += str(self._universe.getId()) + '\n'
    string += universe_id

    cell_id = '{0: <16}'.format('\tCell') + '=\t' + str(self._cell.getId())
    string += cell_id

    return string



class LatCoords(LocalCoords):

  def __init__(self, point=None, next=None, prev=None,
               lattice=None, lattice_x=None, lattice_y=None):

    super(self, point, next, prev)

    self._type = 'lattice'
    self._lattice = None
    self._lattice_x = None
    self._lattice_y = None

    if not lattice is None:
      self.setLatticeId(lattice)

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
           'is not a Lattice', str(lattice))

    self._lattice = lattice


  def setLatticeX(self, lattice_x):

    if not is_integer(lattice_x):
      exit('Unable to set the Lattice X to %s for LocalCoords since it '
           'is not an integer', str(lattice_x))

    self._lattice_x = lattice_x


  def setLatticeY(self, lattice_y):

    if not is_integer(lattice_y):
      exit('Unable to set the Lattice Y to %s for LocalCoords since it '
           'is not an integer', str(lattice_y))

    self._lattice_y = lattice_y


  def toString(self):

    string = super(self)

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