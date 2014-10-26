__author__ = 'Davis Tran'
__email__ = 'dvtran@mit.edu'


from opencsg import *
import numpy as np
import copy
import os
import h5py

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

class Segment(object):

  def __init__(self, geometry=None, start=None, end=None, region_id=None, cell_id=None):

    self._region_id = None
    self._cell_id = None
    self._length = None

    if not region_id is None:
      self.setRegion(region_id)
    elif (not start is None) and (not geometry is None):
      x,y,z = start._coords
      self.setRegion(geometry.getRegionId(x=x,y=y,z=z))

    if not cell_id is None:
      self.setCell(cell_id)
    elif (not start is None) and (not geometry is None):
      x,y,z = start._coords
      self.setCell(geometry.findCell(x=x,y=y,z=z)._id)

    if not (start is None and end is None):
      self._length = start.distanceToPoint(end)


  def __deepcopy__(self, memo):

    existing = memo.get(self)

    # If this is the first time we have tried to copy this object, create a copy
    if existing is None:

      clone = type(self).__new__(type(self))
      clone._region_id = copy.deepcopy(self._region_id)
      clone._cell = copy.deepcopy(self._cell)
      clone._length = copy.deepcopy(self._length)
      return clone

    # If this object has been copied before, return the first copy made
    else:
      return existing

  def setRegion(self, region_id):
    if not is_integer(region_id):
      msg = 'Unable to set region id for segment to {0} since it is ' \
            'not an integer'.format(region_id)
      raise ValueError(msg)

    self._region_id = region_id

  def setCell(self, cell_id):
    if not is_integer(cell_id):
      msg = 'Unable to set cell for segment to {0} since it is ' \
            'not an integer'.format(cell_id)
      raise ValueError(msg)

    self._cell_id = cell_id

  def getMaterial(self):
    return self._cell._fill

  def __repr__(self):

    string = 'Segment\n'
    string += '{0: <16}{1}{2}\n'.format('\tRegion Id', '=\t', self._region_id)
    string += '{0: <16}{1}{2}\n'.format('\tLength', '=\t', self._length)
    return string

def exportRays(rays, directory = 'csg-data/', filename = 'rays-data.h5'):

  # creates folder to contain rays data file if one does not exist
  if not os.path.exists(directory):
    os.makedirs(directory)

  f = h5py.File(directory + filename, 'w')
  f.attrs['Number of Rays'] = len(rays)
  rays_group = f.create_group('Rays')

  # create groups for each ray
  for i in xrange(len(rays)):
    ray_group = rays_group.create_group('Ray (%d)' % (i))
    ray_group.create_dataset('Start Point', data = rays[i]._point._coords)
    ray_group.create_dataset('Direction', data = rays[i]._direction._comps)
    segments = rays[i]._segments
    segments_group = ray_group.create_group('Segments')
    for j in xrange(len(segments)):
      segment = segments_group.create_group('Segment (%d)' % (j))
      segment.create_dataset('Region ID', data = segments[j]._region_id)
      segment.create_dataset('Cell ID', data = segments[j]._cell_id)
      segment.create_dataset('Length', data = segments[j]._length)

  f.close()

def importRays(directory = 'csg-data/', filename = 'rays-data.h5'):

  # checks to see if folder exists and raises error otherwise
  if not os.path.exists(directory + filename):
    msg = 'Could not find file under specified directory and filename.'
    raise ValueError(msg)

  rays = list()

  f = h5py.File(directory + filename, 'r')

  # extract and create rays from file
  for i in xrange(len(f['Rays'])):
    start = Point()
    direction = Direction()
    start.setCoords(f['Rays']['Ray (%d)' % (i)]['Start Point'])
    direction.setComps(f['Rays']['Ray (%d)' % (i)]['Direction'])
    ray = Ray(start, direction)
    for j in xrange(len(f['Rays']['Ray (%d)' % (i)]['Segments'])):
      segment = Segment()
      segment._region_id = f['Rays']['Ray (%d)' % (i)]['Segments']\
                            ['Segment (%d)' % (j)]['Region ID'][...]
      segment._cell_id = f['Rays']['Ray (%d)' % (i)]['Segments']\
                          ['Segment (%d)' % (j)]['Cell ID'][...]
      segment._length = f['Rays']['Ray (%d)' % (i)]['Segments']\
                        ['Segment (%d)' % (j)]['Length'][...]
      ray.addSegment(segment)

    rays.append(ray)

  f.close()

  return rays