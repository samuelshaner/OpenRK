import numpy as np
import copy


class Point(object):

  def __init__(self, x=0., y=0., z=0.):

    # Initialize coordinates
    self._coords = np.zeros(3, dtype=np.float64)

    self.setX(x)
    self.setY(y)
    self.setZ(z)


  def __deepcopy__(self, memo):

    existing = memo.get(id(self))

    # If this is the first time we have tried to copy this object, create a copy
    if existing is None:

      clone = type(self).__new__(type(self))
      clone._coords = copy.deepcopy(self._coords, memo)

      memo[id(self)] = clone

      return clone

    # If this object has been copied before, return the first copy made
    else:
      return existing


  def setCoords(self, coords):

    self.setX(coords[0])
    self.setY(coords[1])
    self.setZ(coords[2])

  def setX(self, x):

    self._coords[0] = np.float64(x)


  def setY(self, y):

    self._coords[1] = np.float64(y)


  def setZ(self, z):

    self._coords[2] = np.float64(z)


  def distanceToPoint(self, point):
    delta = self._coords - point._coords
    dist = np.sqrt(np.sum(delta**2))
    return dist


  def __repr__(self):

    string = 'Point\n'
    string += '{0: <16}{1}{2}\n'.format('\tCoords', '=\t', self._coords)
    return string


class Direction(object):

  def __init__(self, u=0., v=0., w=0.):

    # Initialize components
    self._comps = np.zeros(3, dtype=np.float64)

    self.setU(u)
    self.setV(v)
    self.setW(w)

    self._is_normalized = False

  def __deepcopy__(self, memo):

    existing = memo.get(self)

    # If this is the first time we have tried to copy this object, create a copy
    if existing is None:

      clone = type(self).__new__(type(self), memo)
      clone._comps = copy.deepcopy(self._comps, memo)
      clone._is_normalized = copy.deepcopy(self._is_normalized, memo)
      return clone

    # If this object has been copied before, return the first copy made
    else:
      return existing


  def setComps(self, comps):

    self.setU(comps[0])
    self.setV(comps[1])
    self.setW(comps[2])

  def setU(self, u):

    self._comps[0] = np.float64(u)
    self._is_normalized = False

  def setV(self, v):

    self._comps[1] = np.float64(v)
    self._is_normalized = False

  def setW(self, w):

    self._comps[2] = np.float64(w)
    self._is_normalized = False

  def normalize(self):

    # Check if direction is normalized
    if self._is_normalized:
      return self._comps

    # Normalizes direction otherwise
    comps = self._comps
    unit = comps/np.sqrt(np.sum(comps**2))
    self._is_normalized = True
    return unit

  def toPolar(self):

    u, v, w = self._comps
    square_comps = self._comps**2
    r = np.sqrt(np.sum(square_comps))
    phi = np.arctan2(v, u)
    theta = np.arcsin(np.sqrt((square_comps[0] + square_comps[1]))/r)
    
    return np.array([r,phi,theta])

  def __repr__(self):

    string = 'Direction\n'
    string += '{0: <16}{1}{2}\n'.format('\tComponents', '=\t', self._comps)
    return string