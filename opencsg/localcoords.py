__author__ = 'Will Boyd'
__email__ = 'wboyd@mit.edu'


from point import *

coords_types = ['universe', 'lattice']


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

    if not isinstance(next, (None, LocalCoords)):
      exit('Unable to set the next to %s for LocalCoords since it is not '
           'a LocalCoords object', str(next))

    self._next = next


  def setPrev(self, prev):

    if not isinstance(prev, (None, LocalCoords)):
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

    string += 'Surface\n'

    type = '{0: <16}'.format('\tType') + '=\t' + str(self._type)
    string += type + '\n'

    point = '{0: <16}'.format('\tPoint') + '=\t' + str(self._point.getCoords())
    string += point + '\n'

    return string


  def printString(self):
    print(self.toString())



class UnivCoords(LocalCoords):

  def __init__(self, point=None, next=None, prev=None,
               universe_id=None, cell_id=None):

    super(self, point, next, prev)

    self._type = 'universe'
    self._universe = None
    self._cell_id = None

    if not universe_id is None:
      self.setUniverseId(universe_id)

    if not cell_id is None:
      self.setCellId(cell_id)


    def getUniverseId(self):
      return self._universe_id


    def getCellId(self):
      return self._cell_id


    def setUniverseId(self, universe_id):

      if not is_integer(universe_id):
        exit('Unable to set the Universe ID to %s for LocalCoords since it '
             'is not an integer', str(universe_id))

      self._universe_id = universe_id


    def setCellId(self, cell_id):

      if not is_integer(cell_id):
        exit('Unable to set the Cell ID to %s for LocalCoords since it '
             'is not an integer', str(cell_id))

      self._cell_id = cell_id


    def toString(self):

      string = super(self)

      universe_id = '{0: <16}'.format('\tUniverse') + '=\t'
      universe_id += str(self._universe_id) + '\n'
      string += universe_id

      cell_id = '{0: <16}'.format('\tCell') + '=\t' + str(self._cell_id)
      string += cell_id

      return string



class LatCoords(LocalCoords):

  def __init__(self, point=None, next=None, prev=None,
               lattice_id=None, lattice_x=None, lattice_y=None):

    super(self, point, next, prev)

    self._type = 'lattice'
    self._lattice_id = None
    self._lattice_x = None
    self._lattice_y = None

    if not lattice_id is None:
      self.setLatticeId(lattice_id)

    if not lattice_x is None:
      self.setLatticeX(lattice_x)

    if not lattice_y is None:
      self.setLatticeY(lattice_y)


  def getLatticeId(self):
    return self._lattice_id


  def getLatticeX(self):
    return self._lattice_x


  def getLatticeY(self):
    return self._lattice_y


  def setLatticeId(self, lattice_id):

    if not is_integer(lattice_id):
      exit('Unable to set the Lattice ID to %s for LocalCoords since it '
           'is not an integer', str(lattice_id))

    self._lattice_id = lattice_id


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

    lattice_id = '{0: <16}'.format('\tLattuce') + '=\t'
    lattice_id += str(self._lattice_id) + '\n'
    string += lattice_id

    lattice_x = '{0: <16}'.format('\tLattice X') + '=\t'
    lattice_x += str(self._lattice_x) + '\n'
    string += lattice_x

    lattice_y = '{0: <16}'.format('\tLattice Y') + '=\t'
    lattice_y += str(self._lattice_y)
    string += lattice_y

    return string