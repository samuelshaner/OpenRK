__author__ = 'Will Boyd'
__email__ = 'wboyd@mit.edu'


from opencsg.checkvalue import *
import numpy as np

class Point(object):

  def __init__(self, x=0., y=0., z=0.):

    # Initialize coordinates
    self._coords = np.zeros(3, dtype=np.float64)

    if not x is None:
      self.setX(x)

    if not y is None:
      self.setY(y)

    if not z is None:
      self.setZ(z)


  def setCoords(self, coords):

    if not isinstance(coords, tuple) and len(coords) != 3:
      msg = 'Unable to set coords for point to {0} since it is ' \
            'not a 3D tuple'.format(coords)
      raise ValueError(msg)

    self.setX(coords[0])
    self.setY(coords[1])
    self.setZ(coords[2])

  def setX(self, x):
    if not is_integer(x) and not is_float(x):
      msg = 'Unable to set x coordinate for point to {0} since it is ' \
            'not an integer or floating point value'.format(x)
      raise ValueError(msg)

    self._coords[0] = np.float64(x)


  def setY(self, y):
    if not is_integer(y) and not is_float(y):
      msg = 'Unable to set y coordinate for point to {0} since it is ' \
            'not an integer or floating point value'.format(y)
      raise ValueError(msg)

    self._coords[1] = np.float64(y)


  def setZ(self, z):
    if not is_integer(z) and not is_float(z):
      msg = 'Unable to set z coordinate for point to {0} since it is ' \
            'not an integer or floating point value'.format(z)
      raise ValueError(msg)

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

    if not u is None:
      self.setU(u)

    if not v is None:
      self.setV(v)

    if not w is None:
      self.setW(w)


  def setComps(self, comps):

    if not isinstance(comps, tuple) and len(comps) != 3:
      msg = 'Unable to set comps for direction to {0} since it is ' \
            'not a 3D tuple'.format(comps)
      raise ValueError(msg)

    self.setU(comps[0])
    self.setV(comps[1])
    self.setW(comps[2])

  def setU(self, u):
    if not is_integer(u) and not is_float(u):
      msg = 'Unable to set u component for direction to {0} since it is ' \
            'not an integer or floating point value'.format(u)
      raise ValueError(msg)

    self._comps[0] = np.float64(u)


  def setV(self, v):
    if not is_integer(v) and not is_float(v):
      msg = 'Unable to set v component for point to {0} since it is ' \
            'not an integer or floating point value'.format(v)
      raise ValueError(msg)

    self._comps[1] = np.float64(v)


  def setW(self, w):
    if not is_integer(w) and not is_float(w):
      msg = 'Unable to set w component for point to {0} since it is ' \
            'not an integer or floating point value'.format(w)
      raise ValueError(msg)

    self._comps[2] = np.float64(w)

  def normalize(self):

    comps = self._comps
    unit = comps**2/np.sqrt(np.sum(comps**2))
    unit = unit*(comps/abs(comps))
    return unit

  def toPolar(self):

    u, v, w = self._comps
    r = np.sqrt(np.sum(self._comps**2))
    phi = np.arctan2(v, u)
    theta = np.arcsin(np.sqrt(u**2 + v**2)/r)
    
    return np.array([r,phi,theta])

  def __repr__(self):

    string = 'Direction\n'
    string += '{0: <16}{1}{2}\n'.format('\tComps', '=\t', self._comps)
    return string
