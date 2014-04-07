__author__ = 'Will Boyd'
__email__ = 'wboyd@mit.edu'

from checkvalue import *
import numpy as np


def Point(object):

  def __init__(self, x=0., y=0., z=0.):

    # Initialize coordinates
    self._coords = np.zeros(3, dtype=np.float64)

    if not x is None:
      self.setX(x)

    if not y is None:
      self.setY(y)

    if not z is None:
      self.setZ(z)


  def getCoords(self):
    return self._coords


  def getX(self):
    return self._coords[0]


  def getY(self):
    return self._coords[1]


  def getZ(self):
    return self._coords[2]


  def setCoords(self, coords):

    if not isinstance(coords, tuple) and len(coords) != 3:
      exit('Unable to set coords for point to %s since it is '
           'not a 3D tuple', str(coords))

    self.setX(coords[0])
    self.setY(coords[1])
    self.setZ(coords[2])

  def setX(self, x):
    if not is_integer(x) and not is_float(x):
      exit('Unable to set x coordinate for point to %s since it is not an '
           'integer or floating point value', str(x))

    self._coords[0] = np.float64(x)


  def setY(self, y):
    if not is_integer(y) and not is_float(y):
      exit('Unable to set y coordinate for point to %s since it is not an '
           'integer or floating point value', str(y))

    self._coords[1] = np.float64(y)


  def setZ(self, z):
    if not is_integer(z) and not is_float(z):
      exit('Unable to set z coordinate for point to %s since it is not an '
           'integer or floating point value', str(z))

    self._coords[2] = np.float64(z)


  def distanceToPoint(self, point):
    delta = self._coords - point.getCoords()
    dist = np.sqrt(np.sum(delta)**2)
    return dist


  def toString(self):

    string = ''

    string += 'Point\n'

    coords = '{0: <16}'.format('\tCoords') + '=\t' + str(self._coords)
    string += coords

    return string


  def printString(self):
    print(self.toString())
