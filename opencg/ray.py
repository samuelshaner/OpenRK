__author__ = 'Davis Tran'
__email__ = 'dvtran@mit.edu'


from opencg import *
import numpy as np
import copy
import os
import h5py

class Ray(object):

  def __init__(self, point=None, direction=None):

    # Initializing length 10 array in order avoid many resizings
    self._segments = np.empty(shape=(10), dtype=object)
    self._num_segments = 0

    if point is not None:
      self.setPoint(point)

    if direction is not None:
      self.setDirection(direction)

  def __deepcopy__(self, memo):

    existing = memo.get(self)

    # If this is the first time we have tried to copy this object, create a copy
    if existing is None:

      clone = type(self).__new__(type(self), memo)
      clone._point = copy.deepcopy(self._point, memo)
      clone._direction = copy.deepcopy(self._direction, memo)
      clone._segments = copy.deepcopy(self._segments, memo)
      clone._num_segments = copy.deepcopy(self._num_segments, memo)
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

    self._direction = direction

  def addSegment(self, segment):

    if self._segments.shape[0] <= self._num_segments:
      self._segments = np.append(self._segments, np.empty(shape=(10), dtype = object))

    self._segments[self._num_segments] = segment
    self._num_segments += 1

  def __repr__(self):

    string = 'Ray\n'
    string += '{0: <16}{1}{2}\n'.format('\tStart Point', '=\t', self._point._coords)
    string += '{0: <16}{1}{2}\n'.format('\tDirection', '=\t', self._direction._comps)
    string += '{0: <16}{1}{2}\n'.format('\tNumber of Segments', '=\t', self._num_segments)
    return string

class Segment(object):

  def __init__(self, geometry=None, start=None, end=None, region_id=None, cell_id=None):

    self._region_id = None
    self._cell_id = None
    self._length = None

    if region_id is not None:
      self.setRegion(region_id)
    elif (start is not None) and (geometry is not None):
      x,y,z = start._coords
      self.setRegion(geometry.getRegionId(x=x,y=y,z=z))

    if cell_id is not None:
      self.setCell(cell_id)
    elif (start is not None) and (geometry is not None):
      x,y,z = start._coords
      self.setCell(geometry.findCell(x=x,y=y,z=z)._id)

    if (start is not None) and (end is not None):
      self._length = start.distanceToPoint(end)


  def __deepcopy__(self, memo):

    existing = memo.get(self)

    # If this is the first time we have tried to copy this object, create a copy
    if existing is None:

      clone = type(self).__new__(type(self), memo)
      clone._region_id = copy.deepcopy(self._region_id, memo)
      clone._cell = copy.deepcopy(self._cell_id, memo)
      clone._length = copy.deepcopy(self._length, memo)
      return clone

    # If this object has been copied before, return the first copy made
    else:
      return existing

  def setRegion(self, region_id):

    self._region_id = region_id

  def setCell(self, cell_id):

    self._cell_id = cell_id

  def getMaterial(self, universe):
    return universe._cells[self._cell_id]._fill

  def __repr__(self):

    string = 'Segment\n'
    string += '{0: <16}{1}{2}\n'.format('\tRegion Id', '=\t', self._region_id)
    string += '{0: <16}{1}{2}\n'.format('\tCell Id', '=\t', self._cell_id)
    string += '{0: <16}{1}{2}\n'.format('\tLength', '=\t', self._length)
    return string

def exportRays(rays, directory = 'ray-segments/', filename = 'rays-data.h5'):

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
    for j in xrange(ray._num_segments):
      segment = segments_group.create_group('Segment (%d)' % (j))
      segment.create_dataset('Region ID', data = segments[j]._region_id)
      segment.create_dataset('Cell ID', data = segments[j]._cell_id)
      segment.create_dataset('Length', data = segments[j]._length)

  f.close()

def importRays(directory = 'ray-segments/', filename = 'rays-data.h5'):

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
    ray = Ray(point=start, direction=direction)
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