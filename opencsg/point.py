__author__ = 'Will Boyd'
__email__ = 'wboyd@mit.edu'


from checkvalue import *
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
    delta = self._coords - point._coords()
    dist = np.sqrt(np.sum(delta**2))
    return dist


  def __repr__(self):

    string = 'Point\n'
    string += '{0: <16}{1}{2}\n'.format('\tCoords', '=\t', self._coords)
    return string
