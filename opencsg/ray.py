__author__ = 'Davis Tran'
__email__ = 'dvtran@mit.edu'


from opencsg import *
from opencsg.point import Segment
import numpy as np
import copy

class Ray(object):

  def __init__(self, point, direction):

    self._point = None
    self._direction = None
    self._segments = np.empty(shape=0, dtype=object)

    if not point is None:
      self.setPoint(point)

    if not direction is None:
      self.setDirection(direction)

  def __deepcopy__(self, memo):

    existing = memo.get(self)

    # If this is the first time we have tried to copy this object, create a copy
    if existing is None:

      clone = type(self).__new__(type(self))
      clone._point = copy.deepcopy(self._point)
      clone._direction = copy.deepcopy(self._direction)
      clone._segments = copy.deepcopy(self._segments)
      return clone

    # If this object has been copied before, return the first copy made
    else:
      return existing


  def setPoint(self, point):

    if not isinstance(point, Point):
      msg = 'Unable to set point for ray to {0} since it is ' \
            'not a point object'.format(point)
      raise ValueError(msg)

    self._point = point


  def setDirection(self, direction):

    if not isinstance(direction, Direction):
      msg = 'Unable to set direction for ray to {0} since it is ' \
            'not a point object'.format(direction)
      raise ValueError(msg)

    self._direction = direction

  def addSegment(self, segment):

    if not isinstance(segment, Segment):
      msg = 'Unable to add segment to ray since it is ' \
            'not a segment object'
      raise ValueError(msg)

    self._segments = np.append(self._segments, segment)


'''
class Segment(object):

  def __init__(self, geometry, start=None, end=None, region_id=None, cell=None):

    self._start = None
    self._end = None
    self._region_id = None
    self._cell = None
    self._length = None

    if not start is None:
      self.setStart(start)

    if not end is None:
      self.setEnd(end)

    if not region_id is None:
      self.setRegion(region_id)
    else:
      x,y,z = self._start._coords
      self.setRegion(geometry.getRegionId(x=x,y=y,z=z))

    if not cell is None:
      self.setCell(cell)
    else:
      x,y,z = self._start._coords
      self.setCell(geometry.findCell(x=x,y=y,z=z))

    if not (start is None and end is None):
      self._length = start.distanceToPoint(end)


  def __deepcopy__(self, memo):

    existing = memo.get(self)

    # If this is the first time we have tried to copy this object, create a copy
    if existing is None:

      clone = type(self).__new__(type(self))
      clone._start = copy.deepcopy(self._start)
      clone._end = copy.deepcopy(self._end)
      clone._region_id = copy.deepcopy(self._region_id)
      clone._cell = copy.deepcopy(self._cell)
      clone._length = copy.deepcopy(self._length)
      return clone

    # If this object has been copied before, return the first copy made
    else:
      return existing


  def setStart(self, start):

    if not isinstance(start, Point):
      msg = 'Unable to set start point for segment to {0} since it is ' \
            'not a point object'.format(start)
      raise ValueError(msg)

    self._start = start

  def setEnd(self, end):

    if not isinstance(end, Point):
      msg = 'Unable to set end point for segment to {0} since it is ' \
            'not a point object'.format(end)
      raise ValueError(msg)

    self._end = end

  def setRegion(self, region_id):
    if not is_integer(region_id):
      msg = 'Unable to set region id for segment to {0} since it is ' \
            'not an integer'.format(region_id)
      raise ValueError(msg)

    self._region_id = region_id

  def setCell(self, cell):
    if not str(type(cell)) == '<class \'opencsg.universe.Cell\'>':
      msg = 'Unable to set cell for segment to {0} since it is ' \
            'not a cell object'.format(cell)
      raise ValueError(msg)

    self._cell = cell

  def getXYCoords(self):
    return [self._start._coords[:2], self._end._coords[:2]]

  def getYZCoords(self):
    return [self._start._coords[1:], self._end._coords[1:]]

  def getXZCoords(self):
    return [self._start._coords[::2], self._end._coords[::2]]

  def getMaterial(self):
    return self._cell._fill

  def __repr__(self):

    string = 'Surface\n'
    string += '{0: <16}{1}{2}\n'.format('\tStart', '=\t', self._start._coords)
    string += '{0: <16}{1}{2}\n'.format('\tEnd', '=\t', self._end._coords)
    string += '{0: <16}{1}{2}\n'.format('\tRegion Id', '=\t', self._region_id)
    string += '{0: <16}{1}{2}\n'.format('\tLength', '=\t', self._length)
    return string
'''