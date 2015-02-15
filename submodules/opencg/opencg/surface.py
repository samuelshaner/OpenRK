import opencg
from opencg.checkvalue import *
import numpy as np
import copy


# Threshold for determining how close a point must be to a surface to be on it
ON_SURFACE_THRESH = 1e-12

# Threshold for determining if particle is travelling parallel to axis
PARALLEL_TO_AXIS_THRESH = 1e-5

# A static variable for auto-generated Surface UIDs
AUTO_SURFACE_UID = 1

# The Surface boundary conditions
BOUNDARY_TYPES = ['interface', 'vacuum', 'reflective']

# The types of Surfaces
SURF_TYPES = ['plane',
              'x-plane',
              'y-plane',
              'z-plane',
              'x-cylinder',
              'y-cylinder',
              'z-cylinder',
              'sphere',
              'x-squareprism',
              'y-squareprism',
              'z-squareprism']

MAX_FLOAT = np.finfo(np.float64).max
MIN_FLOAT = np.finfo(np.float64).min


class Surface(object):

  def __init__(self, surface_id=None, name='', boundary='interface'):

    # Initialize class attributes
    global AUTO_SURFACE_UID
    self._uid = AUTO_SURFACE_UID
    AUTO_SURFACE_UID += 1

    self._id = None
    self._set_id = False
    self._name = ''
    self._type = ''
    self._boundary_type = ''
    self._neighbor_cells = dict()
    self._neighbor_cells[-1] = set()
    self._neighbor_cells[+1] = set()

    # A dictionary of the quadratic surface coefficients
    # Key - coefficient name
    # Value - coefficient value
    self._coeffs = dict()

    # Max/min values
    self._max_x = MAX_FLOAT
    self._max_y = MAX_FLOAT
    self._max_z = MAX_FLOAT
    self._min_x = MIN_FLOAT
    self._min_y = MIN_FLOAT
    self._min_z = MIN_FLOAT

    if not surface_id is None:
      self.setId(surface_id)
    else:
      self.setId(self._uid)
      self._set_id = False

    self.setName(name)
    self.setBoundaryType(boundary)


  def __deepcopy__(self, memo):

    existing = memo.get(id(self))

    # If this is the first time we have tried to copy this object, create a copy
    if existing is None:

      clone = type(self).__new__(type(self))
      clone._uid = self._uid
      clone._id = self._id
      clone._set_id = self._set_id
      clone._name = self._name
      clone._type = self._type
      clone._boundary_type = self._boundary_type
      clone._coeffs = copy.deepcopy(self._coeffs, memo)

      clone._neighbor_cells = dict()
      clone._neighbor_cells[-1] = set()
      clone._neighbor_cells[+1] = set()

      for halfspace in [+1, -1]:
        for cell in self._neighbor_cells[halfspace]:
          clone._neighbor_cells[halfspace] = copy.deepcopy(cell, memo)

      clone._max_x = self._max_x
      clone._min_x = self._min_x
      clone._max_y = self._max_y
      clone._min_y = self._min_y
      clone._max_z = self._max_z
      clone._min_z = self._min_z

      memo[id(self)] = clone

      return clone

    # If this object has been copied before, return the first copy made
    else:
      return existing


  def __gt__(self, other):
    return (id(self) > id(other))


  def __ge__(self, other):
    return (id(self) >= id(other))


  def __lt__(self, other):
    return (id(self) < id(other))


  def __le__(self, other):
    return (id(self) <= id(other))


  def getMaxX(self, halfspace=None):
    return self._max_x


  def getMaxY(self, halfspace=None):
    return self._max_y


  def getMaxZ(self, halfspace=None):
    return self._max_z


  def getMinX(self, halfspace=None):
    return self._min_x


  def getMinY(self, halfspace=None):
    return self._min_y


  def getMinZ(self, halfspace=None):
    return self._min_z


  def evaluate(self, point):
    if not isinstance(point, opencg.Point):
      msg = 'Unable to evaluate point for Surface ID={0} since the input ' \
            'is not a Point object'.format(self._id)
      raise ValueError(msg)


  def onSurface(self, point):
    if not isinstance(point, opencg.Point):
      msg = 'Unable to determine whether a point is on Surface ID={0} since ' \
            'the input is not a Point object'.format(self._id)
      raise ValueError(msg)

    if np.abs(self.evaluate(point)) < ON_SURFACE_THRESH:
      return True
    else:
      return False


  def setId(self, surface_id):

    # Check that the ID is a non-negative integer
    if is_integer(surface_id):

      if surface_id >= 0:
        self._id = surface_id
        self._set_id = True
      else:
        msg = 'Unable to set Surface ID to {0} since it must be a ' \
              'non-negative integer'.format(surface_id)
        raise ValueError(msg)

    else:
      msg = 'Unable to set a non-integer Surface ID {0}'.format(surface_id)
      raise ValueError(msg)


  def setName(self, name):

    if not is_string(name):
      msg = 'Unable to set name for Surface ID={0} with a non-string ' \
            'value {1}'.format(self._id, name)
      raise ValueError(msg)

    else:
      self._name = name


  def setBoundaryType(self, boundary):

    if not is_string(boundary):
      msg = 'Unable to set boundary type for Surface ID={0} with a ' \
            'non-string value {1}'.format(self._id, boundary)
      raise ValueError(msg)

    elif not boundary in BOUNDARY_TYPES:
      msg = 'Unable to set boundary type for Surface ID={0} to {1} which ' \
            'is not interface, vacuum or reflective'.format(self._id, boundary)
      raise ValueError(msg)

    else:
      self._boundary_type = boundary


  def addNeighborCell(self, cell, halfspace):

    if not 'opencg.universe.Cell' in str(type(cell)):
      msg = 'Unable to add a neighbor Cell to Surface ID={0} ' \
            'since {1} is not a Cell object'.format(self._id, type(cell))
      raise ValueError(msg)

    elif not halfspace in [-1, +1]:
      msg = 'Unable to add a neighbor Cell to Surface ID={0} with ' \
            'halfspace {1} since it is not an +/-1'.format(self._id, halfspace)
      raise ValueError(msg)

    # Add the Cell to the neighbor Cell collection
    self._neighbor_cells[halfspace].add(cell)

    # Update Cells with the neighbor Cells on the opposite Surface halfspace
    for halfspace in [-1, +1]:
      for cell in self._neighbor_cells[halfspace]:
        for neighbor_cell in self._neighbor_cells[-halfspace]:
          cell.addNeighborCell(neighbor_cell)


  def __repr__(self):

    string = 'Surface\n'
    string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self._id)
    string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self._name)
    string += '{0: <16}{1}{2}\n'.format('\tType', '=\t', self._type)
    string += '{0: <16}{1}{2}\n'.format('\tBoundary', '=\t', self._boundary_type)

    string += '{0: <16}\n'.format('\tCoefficients')

    if len(self._coeffs) > 0:

      for coeff in self._coeffs.keys():
        string += '{0: <16}{1}'.format('\t{0}'.format(coeff), '=\t')
        string += '{0: <12}\n'.format(self._coeffs[coeff])

    return string


class Plane(Surface):

  def __init__(self, surface_id=None, name='', boundary='interface',
               A=None, B=None, C=None, D=None):

    # Initialize Plane class attributes
    super(Plane, self).__init__(surface_id, name, boundary)

    self._type = 'plane'
    self._coeffs['A'] = None
    self._coeffs['B'] = None
    self._coeffs['C'] = None
    self._coeffs['D'] = None

    if A is not None:
      self.setA(A)

    if B is not None:
      self.setB(B)

    if C is not None:
      self.setC(C)

    if D is not None:
      self.setD(D)


  def setA(self, A):

    if not is_integer(A) and not is_float(A):
      msg = 'Unable to set A coefficient for Plane ID={0} to a non-integer ' \
            'or floating point value {1}'.format(self._id, A)
      raise ValueError(msg)

    self._coeffs['A'] = np.float64(A)


  def setB(self, B):

    if not is_integer(B) and not is_float(B):
      msg = 'Unable to set B coefficient for Plane ID={0} to a non-integer ' \
            'or floating point value {1}'.format(self._id, B)
      raise ValueError(msg)

    self._coeffs['B'] = np.float64(B)


  def setC(self, C):

    if not is_integer(C) and not is_float(C):
      msg = 'Unable to set C coefficient for Plane ID={0} to a non-integer ' \
            'or floating point value {1}'.format(self._id, C)
      raise ValueError(msg)

    self._coeffs['C'] = np.float64(C)


  def setD(self, D):

    if not is_integer(D) and not is_float(D):
      msg = 'Unable to set D coefficient for Plane ID={0} to a non-integer ' \
            'or floating point value {1}'.format(self._id, D)
      raise ValueError(msg)

    self._coeffs['D'] = np.float64(D)


  def evaluate(self, point):

    super(Plane, self).evaluate(point)

    x, y, z = point._coords

    value = self._coeffs['A'] * x + \
            self._coeffs['B'] * y + \
            self._coeffs['C'] * z + \
            self._coeffs['D']

    return value

  def minSurfaceDist(self, point, direction):

    if self.onSurface(point):
      return 0.

    x, y, z = point._coords
    u, v, w = direction.normalize()
    
    numerator = self._coeffs['D'] - \
                self._coeffs['A'] * x - \
                self._coeffs['B'] * y - \
                self._coeffs['C'] * z
    denominator = self._coeffs['A'] * u + \
                  self._coeffs['B'] * v + \
                  self._coeffs['C'] * w

    if abs(denominator) < ON_SURFACE_THRESH:
      return None

    dist = numerator/denominator

    if dist < 0:
      return None

    intersect = opencg.Point()
    intersect.setCoords((x+dist*u, y+dist*v, z+dist*w))
    dist = point.distanceToPoint(intersect)
    return dist

class XPlane(Plane):

  def __init__(self, surface_id=None, name='',
               boundary='interface', x0=None):

    # Initialize XPlane class attributes
    super(XPlane, self).__init__(surface_id, name, boundary,
                                 A=1., B=0., C=0., D=-x0)

    self._type = 'x-plane'
    self._coeffs['x0'] = None

    if x0 is not None:
      self.setX0(x0)


  def getMaxX(self, halfspace=None):

    if halfspace is None:
      return self._max_x

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the maximum x-coordinate for ' \
              'Surface ID {0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)


      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the maximum x-coordinate for ' \
              'Surface ID={0} since the halfspace is {1} ' \
              'which is not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._max_x

      else:
        return MAX_FLOAT


  def getMinX(self, halfspace=None):

    if halfspace is None:
      return self._min_x

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the minimum x-coordinate for ' \
              'Surface ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the minimum x-coordinate for ' \
              'Surface ID={0} since the halfspace is {1} which ' \
              'is not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return MIN_FLOAT

      else:
        return self._min_x


  def setX0(self, x0):

    if not is_integer(x0) and not is_float(x0):
      msg = 'Unable to set x0 coefficient for XPlane ID={0} to a non-integer ' \
            'or floating point value {1}'.format(self._id, x0)
      raise ValueError(msg)

    self._coeffs['x0'] = np.float64(x0)
    self.setD(-x0)
    self._max_x = np.float64(x0)
    self._min_x = np.float64(x0)


  def minSurfaceDist(self, point, direction):

    super(XPlane, self).minSurfaceDist(point, direction)

    if self.onSurface(point):
      return 0.

    x, y, z = point._coords
    u, v, w = direction.normalize()

    if abs(u) < ON_SURFACE_THRESH:
      return None

    dist = (self._coeffs['x0'] - x)/u

    if dist < 0:
      return None

    intersect = opencg.Point()
    intersect.setCoords((x+dist*u, y+dist*v, z+dist*w))
    dist = point.distanceToPoint(intersect)
    return dist

class YPlane(Plane):

  def __init__(self, surface_id=None, name='',
               boundary='interface', y0=None):

    # Initialize YPlane class attributes
    super(YPlane, self).__init__(surface_id, name, boundary,
                                 A=0., B=1., C=0., D=-y0)

    self._type = 'y-plane'
    self._coeffs['y0'] = None

    if y0 is not None:
      self.setY0(y0)


  def getMaxY(self, halfspace=None):

    if halfspace is None:
      return self._max_y

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the maximum y-coordinate for ' \
              'Surface ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the maximum y-coordinate for ' \
              'Surface ID={0} since the halfspace is {1} ' \
              'which is not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._max_y

      else:
        return MAX_FLOAT


  def getMinY(self, halfspace=None):

    if halfspace is None:
      return self._min_y

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the minimum y-coordinate for ' \
              'XPlane ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the minimum y-coordinate for ' \
              'XPlane ID={0} since the halfspace is {1} which ' \
              'is not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return MIN_FLOAT

      else:
        return self._min_y


  def setY0(self, y0):

    if not is_integer(y0) and not is_float(y0):
      msg = 'Unable to set y0 coefficient for YPlane ID={0} to a ' \
            'non-integer or floating point value {1}'.format(self._id, y0)
      raise ValueError(msg)

    self._coeffs['y0'] = np.float64(y0)
    self.setD(-y0)
    self._max_y = np.float64(y0)
    self._min_y = np.float64(y0)


  def minSurfaceDist(self, point, direction):

    super(YPlane, self).minSurfaceDist(point, direction)

    if self.onSurface(point):
      return 0.

    x, y, z = point._coords
    u, v, w = direction.normalize()

    if abs(v) < ON_SURFACE_THRESH:
      return None

    dist = (self._coeffs['y0'] - y)/v

    if dist < 0:
      return None

    intersect = opencg.Point()
    intersect.setCoords((x+dist*u, y+dist*v, z+dist*w))
    dist = point.distanceToPoint(intersect)
    return dist

class ZPlane(Plane):

  def __init__(self, surface_id=None, name='',
               boundary='interface', z0=None):

    # Initialize ZPlane class attributes
    super(ZPlane, self).__init__(surface_id, name, boundary,
                                 A=0., B=0., C=1., D=-z0)

    self._type = 'z-plane'
    self._coeffs['z0'] = None

    if z0 is not None:
      self.setZ0(z0)


  def getMaxZ(self, halfspace=None):

    if halfspace is None:
      return self._max_z

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the maximum z-coordinate for ' \
              'YPlane ID={0} since the halfspace is a ' \
              'non-integer value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the maximum z-coordinate for ' \
              'YPlane ID={0} since the halfspace is {1} which ' \
              'is not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._max_z

      else:
        return MAX_FLOAT


  def getMinZ(self, halfspace=None):

    if halfspace is None:
      return self._min_z

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the minimum z-coordinate for ' \
              'ZPlane ID={0} since the halfspace is a ' \
              'non-integer value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the minimum z-coordinate for ' \
              'ZPlane ID={0} since the halfspace is {1} which ' \
              'is not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return MIN_FLOAT

      else:
        return self._min_z


  def setZ0(self, z0):

    if not is_integer(z0) and not is_float(z0):
      msg = 'Unable to set z0 coefficient for ZPlane ID={0} to a non-integer ' \
            'or floating point value {1}'.format(self._id, z0)
      raise ValueError(msg)

    self._coeffs['z0'] = np.float64(z0)
    self.setD(-z0)
    self._max_z = z0
    self._min_z = z0


  def minSurfaceDist(self, point, direction):

    super(ZPlane, self).minSurfaceDist(point, direction)

    if self.onSurface(point):
      return 0.

    x, y, z = point._coords
    u, v, w = direction.normalize()

    if abs(w) < ON_SURFACE_THRESH:
      return None

    dist = (self._coeffs['z0'] - z)/w

    if dist < 0:
      return None

    intersect = opencg.Point()
    intersect.setCoords((x+dist*u, y+dist*v, z+dist*w))
    dist = point.distanceToPoint(intersect)
    return dist

class Cylinder(Surface):

  def __init__(self, surface_id=None, name='',
               boundary='interface', R=None):

    # Initialize Cylinder class attributes
    super(Cylinder, self).__init__(surface_id, name, boundary)

    self._coeffs['R'] = None


  def setR(self, R):

    if not is_integer(R) and not is_float(R):
      msg = 'Unable to set R coefficient for Cylinder ID={0} to ' \
            'a non-integer or floating point value {1}'.format(self._id, R)
      raise ValueError(msg)

    self._coeffs['R'] = np.float64(R)



class XCylinder(Cylinder):

  def __init__(self, surface_id=None, name='',
               boundary='interface', y0=None, z0=None, R=None):

    # Initialize XCylinder class attributes
    super(XCylinder, self).__init__(surface_id, name, boundary, R)

    self._type = 'x-cylinder'
    self._coeffs['y0'] = 0.
    self._coeffs['z0'] = 0.
    self._coeffs['R'] = None

    if y0 is not None:
      self.setY0(y0)

    if z0 is not None:
      self.setZ0(z0)

    if R is not None:
      self.setR(R)


  def getMaxY(self, halfspace=None):

    if halfspace is None:
      return self._max_y

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the maximum y-coordinate for ' \
              'XClinder ID={0} since the halfspace is a ' \
              'non-integer value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the maximum y-coordinate for ' \
              'XCylinder ID={0} since the halfspace is {1} ' \
              'which is not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._max_y

      else:
        return MAX_FLOAT


  def getMinY(self, halfspace=None):

    if halfspace is None:
      return self._min_y

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the minimum y-coordinate for ' \
              'XCylinder ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the minimum y-coordinate for ' \
              'XCylinder ID={0} since the halfspace is {1} ' \
              'which is not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._min_y

      else:
        return MIN_FLOAT


  def getMaxZ(self, halfspace=None):

    if halfspace is None:
      return self._max_z

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the maximum z-coordinate for ' \
              'XClinder ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the maximum z-coordinate for ' \
              'XCylinder ID={0} since the halfspace is {1} ' \
              'which is not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._max_z

      else:
        return MAX_FLOAT


  def getMinZ(self, halfspace=None):

    if halfspace is None:
      return self._min_z

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the minimum z-coordinate for ' \
              'XCylinder ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the minimum z-coordinate for ' \
              'XCylinder ID={0} since the halfspace is {1} which ' \
              'is not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._min_z

      else:
        return MIN_FLOAT


  def setY0(self, y0):

    if not is_integer(y0) and not is_float(y0):
      msg = 'Unable to set y0 coefficient for XCylinder ID={0} ' \
            'to a non-integer or floating point value {1}'.format(self._id, y0)
      raise ValueError(msg)

    self._coeffs['y0'] = np.float64(y0)

    if self._coeffs['R'] is not None:
      self._max_y = y0 + self._coeffs['R']
      self._min_y = y0 - self._coeffs['R']


  def setZ0(self, z0):

    if not is_integer(z0) and not is_float(z0):
      msg = 'Unable to set z0 coefficient for XCylinder ID={0} ' \
            'to a non-integer or floating point value {1}'.format(self._id, z0)
      raise ValueError(msg)

    self._coeffs['z0'] = np.float64(z0)

    if self._coeffs['R'] is not None:
      self._max_z = z0 + self._coeffs['R']
      self._min_z = z0 - self._coeffs['R']


  def setR(self, R):

    super(XCylinder, self).setR(R)

    if self._coeffs['y0'] is not None:
      self._max_y = self._coeffs['y0'] + self._coeffs['R']
      self._min_y = self._coeffs['y0'] - self._coeffs['R']

    if self._coeffs['z0'] is not None:
      self._max_z = self._coeffs['z0'] + self._coeffs['R']
      self._min_z = self._coeffs['z0'] - self._coeffs['R']


  def evaluate(self, point):

    super(XCylinder, self).evaluate(point)

    coords = point.getCoords()
    r2 = (coords[1] - self._coeffs['y0'])**2 + \
        (coords[2] - self._coeffs['z0'])**2
    r = np.sqrt(r2)
    return (r - self._coeffs['R'])

  def minSurfaceDist(self, point, direction):

    if self.onSurface(point):
      return 0.

    x, y, z = point._coords
    u, v, w = direction.normalize()

    ybar = y-self._coeffs['y0']
    zbar = z-self._coeffs['z0']
    a = v**2 + w**2
    k = ybar*v + zbar*w
    c = ybar**2 + zbar**2 - self._coeffs['R']**2

    if abs(a) < ON_SURFACE_THRESH or k**2-a*c < 0:
      return None

    if c < 0:
      dist = (-k + np.sqrt(k**2-a*c))/a
      intersect = opencg.Point()
      intersect.setCoords((x+dist*u, y+dist*v, z+dist*w))
      dist = point.distanceToPoint(intersect)
      return dist

    else:
      dist = (-k - np.sqrt(k**2-a*c))/a
      if dist > 0:
        intersect = opencg.Point()
        intersect.setCoords((x+dist*u, y+dist*v, z+dist*w))
        dist = point.distanceToPoint(intersect)
        return dist

      else:
        return None


class YCylinder(Cylinder):

  def __init__(self, surface_id=None, name='',
               boundary='interface', x0=None, z0=None, R=None):

    # Initialize YCylinder class attributes
    super(YCylinder, self).__init__(surface_id, name, boundary, R)

    self._type = 'y-cylinder'
    self._coeffs['x0'] = 0.
    self._coeffs['z0'] = 0.
    self._coeffs['R'] = None

    if x0 is not None:
      self.setX0(x0)

    if z0 is not None:
      self.setZ0(z0)

    if R is not None:
      self.setR(R)


  def getMaxX(self, halfspace=None):

    if halfspace is None:
      return self._max_x

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the maximum x-coordinate for ' \
              'YClinder ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the maximum x-coordinate for ' \
              'YCylinder ID={0} since the halfspace is {1} which ' \
              'is not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._max_x

      else:
        return MAX_FLOAT


  def getMinX(self, halfspace=None):

    if halfspace is None:
      return self._min_x

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the minimum x-coordinate for ' \
              'YCylinder ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the minimum x-coordinate for ' \
              'YCylinder ID={0} since the halfspace is {1} which ' \
              'is not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._min_x

      else:
        return MIN_FLOAT


  def getMaxZ(self, halfspace=None):

    if halfspace is None:
      return self._max_z

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the maximum z-coordinate for ' \
              'YClinder ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the maximum z-coordinate for ' \
              'YCylinder ID={0} since the halfspace is {1} which ' \
              'is not +/-1'.format(self._id, halfspace)

      elif halfspace == -1:
        return self._max_z

      else:
        return MAX_FLOAT


  def getMinZ(self, halfspace=None):

    if halfspace is None:
      return self._min_z

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the minimum z-coordinate for ' \
              'YCylinder ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the minimum z-coordinate for ' \
              'YCylinder ID={0} since the halfspace is {1} which ' \
              'is not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._min_z

      else:
        return MIN_FLOAT


  def setX0(self, x0):

    if not is_integer(x0) and not is_float(x0):
      msg = 'Unable to set x0 coefficient for YCylinder ID={0} to a ' \
            'non-integer or floating point value {1}'.format(self._id, x0)
      raise ValueError(msg)

    self._coeffs['x0'] = np.float64(x0)

    if self._coeffs['R'] is not None:
      self._max_x = x0 + self._coeffs['R']
      self._min_x = x0 - self._coeffs['R']


  def setZ0(self, z0):

    if not is_integer(z0) and not is_float(z0):
      msg = 'Unable to set z0 coefficient for YCylinder ID={0} to a ' \
            'non-integer or floating point value {1}'.format(self._id, z0)
      raise ValueError(msg)

    self._coeffs['z0'] = np.float64(z0)

    if self._coeffs['R'] is not None:
      self._max_z = z0 + self._coeffs['R']
      self._min_z = z0 - self._coeffs['R']


  def setR(self, R):

    super(YCylinder, self).setR(R)

    if self._coeffs['x0'] is not None:
      self._max_x = self._coeffs['x0'] + self._coeffs['R']
      self._min_x = self._coeffs['x0'] - self._coeffs['R']

    if self._coeffs['z0'] is not None:
      self._max_z = self._coeffs['z0'] + self._coeffs['R']
      self._min_z = self._coeffs['z0'] - self._coeffs['R']


  def evaluate(self, point):

    super(YCylinder, self).evaluate(point)

    coords = point._coords
    r2 = (coords[0] - self._coeffs['x0'])**2 + \
        (coords[2] - self._coeffs['z0'])**2
    r = np.sqrt(r2)
    return (r - self._coeffs['R'])

  def minSurfaceDist(self, point, direction):

    if self.onSurface(point):
      return 0.

    x, y, z = point._coords
    u, v, w = direction.normalize()

    xbar = x-self._coeffs['x0']
    zbar = z-self._coeffs['z0']
    a = u**2 + w**2
    k = xbar*u + zbar*w
    c = xbar**2 + zbar**2 - self._coeffs['R']**2

    if abs(a) < ON_SURFACE_THRESH or k**2-a*c < 0:
      return None

    if c < 0:
      dist = (-k + np.sqrt(k**2-a*c))/a
      intersect = opencg.Point()
      intersect.setCoords((x+dist*u, y+dist*v, z+dist*w))
      dist = point.distanceToPoint(intersect)
      return dist

    else:
      dist = (-k - np.sqrt(k**2-a*c))/a
      if dist > 0:
        intersect = opencg.Point()
        intersect.setCoords((x+dist*u, y+dist*v, z+dist*w))
        dist = point.distanceToPoint(intersect)
        return dist

      else:
        return None


class ZCylinder(Cylinder):

  def __init__(self, surface_id=None, name='',
               boundary='interface', x0=None, y0=None, R=None):

    # Initialize ZCylinder class attributes
    super(ZCylinder, self).__init__(surface_id, name, boundary, R)

    self._type = 'z-cylinder'
    self._coeffs['x0'] = 0.
    self._coeffs['y0'] = 0.
    self._coeffs['R'] = None


    if x0 is not None:
      self.setX0(x0)

    if y0 is not None:
      self.setY0(y0)

    if R is not None:
      self.setR(R)


  def getMaxX(self, halfspace=None):

    if halfspace is None:
      return self._max_x

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the maximum x-coordinate for ' \
              'ZCylinder ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the maximum x-coordinate for ' \
              'ZCylinder ID={0} since the halfspace is {1} which ' \
              'is not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._max_x

      else:
        return MAX_FLOAT


  def getMinX(self, halfspace=None):

    if halfspace is None:
      return self._min_x

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the minimum x-coordinate for ' \
              'ZCylinder ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the minimum x-coordinate for ' \
              'ZCylinder ID={0} since the halfspace is {1} which ' \
              'is not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._min_x

      else:
        return MIN_FLOAT


  def getMaxY(self, halfspace=None):

    if halfspace is None:
      return self._max_y

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the maximum y-coordinate for ' \
              'ZClinder ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the maximum y-coordinate for ' \
              'ZCylinder ID={0} since the halfspace is {1} which ' \
              'is not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._max_y

      else:
        return MAX_FLOAT


  def getMinY(self, halfspace=None):

    if halfspace is None:
      return self._min_y

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the minimum y-coordinate for ' \
              'ZCylinder ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the minimum y-coordinate for ' \
              'ZCylinder ID={0} since the halfspace is {1} which ' \
              'is not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._min_y

      else:
        return MIN_FLOAT


  def setX0(self, x0):

    if not is_integer(x0) and not is_float(x0):
      msg = 'Unable to set x0 coefficient for ZCylinder ID={0} to a ' \
            'non-integer value {1}'.format(self._id, x0)
      raise ValueError(msg)

    self._coeffs['x0'] = np.float64(x0)

    if self._coeffs['R'] is not None:
      self._max_x = x0 + self._coeffs['R']
      self._min_x = x0 - self._coeffs['R']


  def setY0(self, y0):

    if not is_integer(y0) and not is_float(y0):
      msg = 'Unable to set y0 coefficient for ZCylinder ID={0} to a ' \
            'non-integer value {1}'.format(self._id, y0)
      raise ValueError(msg)

    self._coeffs['y0'] = np.float64(y0)

    if self._coeffs['R'] is not None:
      self._max_y = y0 + self._coeffs['R']
      self._min_y = y0 - self._coeffs['R']


  def setR(self, R):

    super(ZCylinder, self).setR(R)

    if self._coeffs['x0'] is not None:
      self._max_x = self._coeffs['x0'] + self._coeffs['R']
      self._min_x = self._coeffs['x0'] - self._coeffs['R']

    if self._coeffs['y0'] is not None:
      self._max_y = self._coeffs['y0'] + self._coeffs['R']
      self._min_y = self._coeffs['y0'] - self._coeffs['R']


  def evaluate(self, point):

    super(ZCylinder, self).evaluate(point)

    coords = point._coords
    r2 = (coords[0] - self._coeffs['x0'])**2 + \
        (coords[1] - self._coeffs['y0'])**2
    r = np.sqrt(r2)
    return (r - self._coeffs['R'])

  def minSurfaceDist(self, point, direction):

    x, y, z = point._coords
    u, v, w = direction.normalize()

    xbar = x-self._coeffs['x0']
    ybar = y-self._coeffs['y0']
    a = u**2 + v**2
    k = xbar*u + ybar*v
    c = xbar**2 + ybar**2 - self._coeffs['R']**2

    if abs(a) < ON_SURFACE_THRESH or k**2-a*c < 0:
      return None

    if c < 0:
      dist = (-k + np.sqrt(k**2-a*c))/a
      intersect = opencg.Point()
      intersect.setCoords((x+dist*u, y+dist*v, z+dist*w))
      dist = point.distanceToPoint(intersect)
      return dist

    else:
      dist = (-k - np.sqrt(k**2-a*c))/a
      if dist > 0:
        intersect = opencg.Point()
        intersect.setCoords((x+dist*u, y+dist*v, z+dist*w))
        dist = point.distanceToPoint(intersect)
        return dist

      else:
        return None

class Sphere(Surface):

  def __init__(self, surface_id=None, name='',
               boundary='interface', x0=None, y0=None, z0=None, R=None):

    # Initialize Sphere class attributes
    super(Sphere, self).__init__(surface_id, name, boundary)

    self._type = 'sphere'
    self._coeffs['x0'] = None
    self._coeffs['y0'] = None
    self._coeffs['z0'] = None
    self._coeffs['R'] = None

    if x0 is not None:
      self.setX0(x0)

    if y0 is not None:
      self.setY0(y0)

    if z0 is not None:
      self.setZ0(z0)

    if R is not None:
      self.setZ0(R)


  def getMaxX(self, halfspace=None):

    if halfspace is None:
      return self._max_x

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the maximum x-coordinate for ' \
              'Sphere ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the maximum x-coordinate for ' \
              'Sphere ID={0} since the halfspace is {1} which ' \
              'is not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._max_x

      else:
        return MAX_FLOAT


  def getMinX(self, halfspace=None):

    if halfspace is None:
      return self._min_x

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the minimum x-coordinate for ' \
              'Sphere ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the minimum x-coordinate for ' \
              'Sphere ID={0} since the halfspace is {1} which ' \
              'is not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._min_x

      else:
        return MIN_FLOAT


  def getMaxY(self, halfspace=None):

    if halfspace is None:
      return self._max_y

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the maximum y-coordinate for ' \
              'Sphere ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the maximum y-coordinate for ' \
              'Sphere ID={0} since the halfspace is {1} which is ' \
              'not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._max_y

      else:
        return MAX_FLOAT


  def getMinY(self, halfspace=None):

    if halfspace is None:
      return self._min_y

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the minimum y-coordinate for ' \
              'Sphere ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the minimum y-coordinate for ' \
              'Sphere ID={0} since the halfspace is {1} which is ' \
              'not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._min_y

      else:
        return MIN_FLOAT


  def getMaxZ(self, halfspace=None):

    if halfspace is None:
      return self._max_z

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the maximum z-coordinate for ' \
              'Sphere ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the maximum z-coordinate for ' \
              'Sphere ID={0} since the halfspace is {1} which ' \
              'is not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._max_z

      else:
        return MAX_FLOAT


  def getMinZ(self, halfspace=None):

    if halfspace is None:
      return self._min_z

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the minimum z-coordinate for ' \
              'Sphere ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the minimum z-coordinate for ' \
              'Sphere ID={0} since the halfspace is {1} which ' \
              'is not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._min_z

      else:
        return MIN_FLOAT


  def setX0(self, x0):

    if not is_integer(x0) and not is_float(x0):
      msg = 'Unable to set x0 coefficient for Sphere ID={0} to a non-integer ' \
            'or floating point value {1}'.format(self._id, x0)
      raise ValueError(msg)

    self._coeffs['x0'] = np.float64(x0)

    if self._coeffs['R'] is not None:
      self._max_x = x0 + self._coeffs['R']
      self._min_x = x0 - self._coeffs['R']


  def setY0(self, y0):

    if not is_integer(y0) and not is_float(y0):
      msg = 'Unable to set y0 coefficient for Sphere ID={0} to a non-integer ' \
            'or floating point value {1}'.format(self._id, y0)
      raise ValueError(msg)

    self._coeffs['y0'] = np.float64(y0)

    if self._coeffs['R'] is not None:
      self._max_y = y0 + self._coeffs['R']
      self._min_y = y0 - self._coeffs['R']


  def setZ0(self, z0):

    if not is_integer(z0) and not is_float(z0):
      msg = 'Unable to set z0 coefficient for Sphere ID={0} to a non-integer ' \
            'or floating point value {1}'.format(self._id, z0)
      raise ValueError(msg)

    self._coeffs['z0'] = np.float64(z0)

    if self._coeffs['R'] is not None:
      self._max_z = z0 + self._coeffs['R']
      self._min_z = z0 - self._coeffs['R']


  def setR(self, R):

    if not is_integer(R) and not is_float(R):
      msg = 'Unable to set R coefficient for Sphere ID={0} to a non-integer ' \
            'or floating point value {1}'.format(self._id, R)
      raise ValueError(msg)

    self._coeffs['R'] = np.float64(R)

    if self._coeffs['x0'] is not None:
      self._max_x = self._coeffs['x0'] + R
      self._min_x = self._coeffs['x0'] - R

    if self._coeffs['y0'] is not None:
      self._max_y = self._coeffs['y0'] + R
      self._min_y = self._coeffs['y0'] - R

    if self._coeffs['z0'] is not None:
      self._max_z = self._coeffs['z0'] + R
      self._min_z = self._coeffs['z0'] - R


  def evaluate(self, point):

    super(Sphere, self).evaluate(point)

    coords = point._coords

    R2 = (self._coeffs['x0'] - coords[0])**2 + \
         (self._coeffs['y0'] - coords[1])**2 + \
         (self._coeffs['z0'] - coords[2])**2
    R = np.sqrt(R2)
    return (R - self._coeffs['R'])

  def minSurfaceDist(self, point, direction):

    if self.onSurface(point):
      return 0.

    x, y, z = point._coords
    u, v, w = direction.normalize()

    xbar = x-self._coeffs['x0']
    ybar = y-self._coeffs['y0']
    zbar = z-self._coeffs['z0']
    k = xbar*u + ybar*v + zbar*w
    c = xbar**2 + ybar**2 + zbar**2 - self._coeffs['R']**2

    if k**2-c < 0:
      return None


    if c < 0:
      dist = (-k + np.sqrt(k**2-c))
      intersect = opencg.Point()
      intersect.setCoords((x+dist*u, y+dist*v, z+dist*w))
      dist = point.distanceToPoint(intersect)
      return dist


    else:
      dist1 = (-k + np.sqrt(k**2-c))
      dist2 = (-k - np.sqrt(k**2-c))
      if abs(dist1) < abs(dist2):
        dist = dist1
      else:
        dist = dist2
      intersect = opencg.Point()
      intersect.setCoords((x+dist*u, y+dist*v, z+dist*w))
      dist = point.distanceToPoint(intersect)
      return dist


class SquarePrism(Surface):

  def __init__(self, surface_id=None, name='',
               boundary='interface', R=None):

    # Initialize Cylinder class attributes
    super(SquarePrism, self).__init__(surface_id, name, boundary)

    self._coeffs['R'] = None


  def setR(self, R):

    if not is_integer(R) and not is_float(R):
      msg = 'Unable to set R coefficient for SquarePrism ID={0} to ' \
            'a non-integer or floating point value {1}'.format(self._id, R)
      raise ValueError(msg)

    self._coeffs['R'] = np.float64(R)



class XSquarePrism(SquarePrism):

  def __init__(self, surface_id=None, name='',
               boundary='interface', y0=None, z0=None, R=None):

    # Initialize XSquarePrism class attributes
    super(XSquarePrism, self).__init__(surface_id, name, boundary, R)

    self._type = 'x-squareprism'
    self._coeffs['y0'] = 0.
    self._coeffs['z0'] = 0.
    self._coeffs['R'] = None

    if y0 is not None:
      self.setY0(y0)

    if z0 is not None:
      self.setZ0(z0)

    if R is not None:
      self.setR(R)


  def getMaxY(self, halfspace=None):

    if halfspace is None:
      return self._max_y

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the maximum y-coordinate for ' \
              'XSquarePrism ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the maximum y-coordinate for ' \
              'XSquarePrism ID={0} since the halfspace is {1} which ' \
              'is not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._max_y

      else:
        return MAX_FLOAT


  def getMinY(self, halfspace=None):

    if halfspace is None:
      return self._min_y

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the minimum y-coordinate for ' \
              'XSquarePrism ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the minimum y-coordinate for ' \
              'XSquarePrism ID={0} since the halfspace is {1} which ' \
              'is not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._min_y

      else:
        return MIN_FLOAT


  def getMaxZ(self, halfspace=None):

    if halfspace is None:
      return self._max_z

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the maximum z-coordinate for ' \
              'XSquarePrism ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the maximum z-coordinate for ' \
              'XSquarePrism ID={0} since the halfspace is {1} which is ' \
              'not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._max_z

      else:
        return MAX_FLOAT


  def getMinZ(self, halfspace=None):

    if halfspace is None:
      return self._min_z

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the minimum z-coordinate for ' \
              'XSquarePrism ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the minimum z-coordinate for ' \
              'XSquarePrism ID={0} since the halfspace is {1} which is ' \
              'not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._min_z

      else:
        return MIN_FLOAT


  def setY0(self, y0):

    if not is_integer(y0) and not is_float(y0):
      msg = 'Unable to set x0 coefficient for XSquarePrism ID={0} to ' \
            'a non-integer or floating point value {1}'.format(self._id, y0)
      raise ValueError(msg)

    self._coeffs['x0'] = np.float64(y0)

    if self._coeffs['R'] is not None:
      self._max_x = y0 + self._coeffs['R']
      self._min_x = y0 - self._coeffs['R']


  def setZ0(self, z0):

    if not is_integer(z0) and not is_float(z0):
      msg = 'Unable to set y0 coefficient for XSquarePrism ID={0} to ' \
            'a non-integer or floating point value {1}'.format(self._id, z0)
      raise ValueError(msg)

    self._coeffs['z0'] = np.float64(z0)

    if self._coeffs['R'] is not None:
      self._max_y = z0 + self._coeffs['R']
      self._min_y = z0 - self._coeffs['R']


  def setR(self, R):

    super(XSquarePrism, self).setR(R)

    if self._coeffs['y0'] is not None:
      self._max_y = self._coeffs['y0'] + R
      self._min_y = self._coeffs['y0'] - R

    if self._coeffs['z0'] is not None:
      self._max_z = self._coeffs['z0'] + R
      self._min_z = self._coeffs['z0'] - R


  def evaluate(self, point):

    super(XSquarePrism, self).evaluate(point)

    x, y, z = point._coords

    Ry = abs(self._coeffs['y0'] - y) - self._coeffs['R']
    Rz = abs(self._coeffs['z0'] - z) - self._coeffs['R']

    return max(Ry, Rz)


  def minSurfaceDist(self, point, direction):

    if self.onSurface(point):
      return 0.

    x, y, z = point._coords
    u, v, w = direction.normalize()

    dist_y = (np.sign(v)*self._coeffs['R'] - y)/v
    dist_z = (np.sign(w)*self._coeffs['R'] - z)/w

    if abs(v) < PARALLEL_TO_AXIS_THRESH:
      dist_y = 0.
    if abs(w) < PARALLEL_TO_AXIS_THRESH:
      dist_z = 0.

    for dist in [dist_y, dist_z]:
      new_x = x + dist*u
      new_y = y + dist*v
      new_z = z + dist*w
      if (abs((abs(new_y) - self._coeffs['R'])) < ON_SURFACE_THRESH and
         abs(new_z) < self._coeffs['R']) or \
         (abs((abs(new_z) - self._coeffs['R'])) < ON_SURFACE_THRESH and
          abs(new_y) < self._coeffs['R']):
        intersect = opencg.Point()
        intersect.setCoords((new_x, new_y, new_z))
        dist = point.distanceToPoint(intersect)
        return dist


class YSquarePrism(SquarePrism):

  def __init__(self, surface_id=None, name='',
               boundary='interface', x0=None, z0=None, R=None):

    # Initialize YSquarePrism class attributes
    super(YSquarePrism, self).__init__(surface_id, name, boundary, R)

    self._type = 'y-squareprism'
    self._coeffs['x0'] = 0.
    self._coeffs['z0'] = 0.
    self._coeffs['R'] = None

    if x0 is not None:
      self.setX0(x0)

    if z0 is not None:
      self.setZ0(z0)

    if R is not None:
      self.setR(R)


  def getMaxX(self, halfspace=None):

    if halfspace is None:
      return self._max_x

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the maximum x-coordinate for ' \
              'YSquarePrism ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the maximum x-coordinate for ' \
              'YSquarePrism ID={0} since the halfspace is {1} which ' \
              'is not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._max_x

      else:
        return MAX_FLOAT


  def getMinX(self, halfspace=None):

    if halfspace is None:
      return self._min_x

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the minimum x-coordinate for ' \
              'YSquarePrism ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the minimum x-coordinate for ' \
              'YSquarePrism ID={0} since the halfspace is {1} which ' \
              'is not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._min_x

      else:
        return MIN_FLOAT


  def getMaxZ(self, halfspace=None):

    if halfspace is None:
      return self._max_z

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the maximum z-coordinate for ' \
              'YSquarePrism ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the maximum z-coordinate for ' \
              'YSquarePrism ID={0} since the halfspace is {1} which is ' \
              'not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._max_z

      else:
        return MAX_FLOAT


  def getMinZ(self, halfspace=None):

    if halfspace is None:
      return self._min_z

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the minimum z-coordinate for ' \
              'YSquarePrism ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the minimum z-coordinate for ' \
              'YSquarePrism ID={0} since the halfspace is {1} which is ' \
              'not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._min_z

      else:
        return MIN_FLOAT


  def setX0(self, x0):

    if not is_integer(x0) and not is_float(x0):
      msg = 'Unable to set x0 coefficient for YSquarePrism ID={0} to ' \
            'a non-integer or floating point value {1}'.format(self._id, x0)
      raise ValueError(msg)

    self._coeffs['x0'] = np.float64(x0)

    if self._coeffs['R'] is not None:
      self._max_x = x0 + self._coeffs['R']
      self._min_x = x0 - self._coeffs['R']


  def setZ0(self, z0):

    if not is_integer(z0) and not is_float(z0):
      msg = 'Unable to set y0 coefficient for YSquarePrism ID={0} to ' \
            'a non-integer or floating point value {1}'.format(self._id, z0)
      raise ValueError(msg)

    self._coeffs['z0'] = np.float64(z0)

    if self._coeffs['R'] is not None:
      self._max_y = z0 + self._coeffs['R']
      self._min_y = z0 - self._coeffs['R']


  def setR(self, R):

    super(YSquarePrism, self).setR(R)

    if self._coeffs['x0'] is not None:
      self._max_x = self._coeffs['x0'] + R
      self._min_x = self._coeffs['x0'] - R

    if self._coeffs['z0'] is not None:
      self._max_z = self._coeffs['z0'] + R
      self._min_z = self._coeffs['z0'] - R


  def evaluate(self, point):

    super(YSquarePrism, self).evaluate(point)

    x, y, z = point._coords

    Rx = abs(self._coeffs['x0'] - x) - self._coeffs['R']
    Rz = abs(self._coeffs['z0'] - z) - self._coeffs['R']

    return max(Rx, Rz)

  def minSurfaceDist(self, point, direction):

    if self.onSurface(point):
      return 0.

    x, y, z = point._coords
    u, v, w = direction.normalize()

    dist_x = (np.sign(u)*self._coeffs['R'] - x)/u
    dist_z = (np.sign(w)*self._coeffs['R'] - z)/w

    if abs(u) < PARALLEL_TO_AXIS_THRESH:
      dist_x = 0.
    if abs(w) < PARALLEL_TO_AXIS_THRESH:
      dist_z = 0.

    for dist in [dist_x, dist_z]:
      new_x = x + dist*u
      new_y = y + dist*v
      new_z = z + dist*w
      if (abs((abs(new_x) - self._coeffs['R'])) < ON_SURFACE_THRESH and
         abs(new_z) < self._coeffs['R']) or \
         (abs((abs(new_z) - self._coeffs['R'])) < ON_SURFACE_THRESH and
          abs(new_x) < self._coeffs['R']):
        intersect = opencg.Point()
        intersect.setCoords((new_x, new_y, new_z))
        dist = point.distanceToPoint(intersect)
        return dist


class ZSquarePrism(SquarePrism):

  def __init__(self, surface_id=None, name='',
               boundary='interface', x0=None, y0=None, R=None):

    # Initialize Square class attributes
    super(ZSquarePrism, self).__init__(surface_id, name, boundary, R)

    self._type = 'z-squareprism'
    self._coeffs['x0'] = 0.
    self._coeffs['y0'] = 0.
    self._coeffs['R'] = None

    if x0 is not None:
      self.setX0(x0)

    if y0 is not None:
      self.setY0(y0)

    if R is not None:
      self.setR(R)


  def getMaxX(self, halfspace=None):

    if halfspace is None:
      return self._max_x

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the maximum x-coordinate for ' \
              'ZSquarePrism ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the maximum x-coordinate for ' \
              'ZSquarePrism ID={0} since the halfspace is {1} which ' \
              'is not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._max_x

      else:
        return MAX_FLOAT


  def getMinX(self, halfspace=None):

    if halfspace is None:
      return self._min_x

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the minimum x-coordinate for ' \
              'ZSquarePrism ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the minimum x-coordinate for ' \
              'ZSquarePrism ID={0} since the halfspace is {1} which ' \
              'is not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._min_x

      else:
        return MIN_FLOAT


  def getMaxY(self, halfspace=None):

    if halfspace is None:
      return self._max_y

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the maximum y-coordinate for ' \
              'ZSquarePrism ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the maximum y-coordinate for ' \
              'ZSquarePrism ID={0} since the halfspace is {1} which is ' \
              'not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._max_y

      else:
        return MAX_FLOAT


  def getMinY(self, halfspace=None):

    if halfspace is None:
      return self._min_y

    else:

      if not is_integer(halfspace):
        msg = 'Unable to get the minimum y-coordinate for ' \
              'ZSquarePrism ID={0} since the halfspace is a non-integer ' \
              'value {1}'.format(self._id, halfspace)
        raise ValueError(msg)

      elif not halfspace in [-1, +1]:
        msg = 'Unable to get the minimum y-coordinate for ' \
              'ZSquarePrism ID={0} since the halfspace is {1} which is ' \
              'not +/-1'.format(self._id, halfspace)
        raise ValueError(msg)

      elif halfspace == -1:
        return self._min_y

      else:
        return MIN_FLOAT


  def setX0(self, x0):

    if not is_integer(x0) and not is_float(x0):
      msg = 'Unable to set x0 coefficient for ZSquarePrism ID={0} to ' \
            'a non-integer or floating point value {1}'.format(self._id, x0)
      raise ValueError(msg)

    self._coeffs['x0'] = np.float64(x0)

    if self._coeffs['R'] is not None:
      self._max_x = x0 + self._coeffs['R']
      self._min_x = x0 - self._coeffs['R']


  def setY0(self, y0):

    if not is_integer(y0) and not is_float(y0):
      msg = 'Unable to set y0 coefficient for ZSquarePrism ID={0} to ' \
            'a non-integer or floating point value {1}'.format(self._id, y0)
      raise ValueError(msg)

    self._coeffs['y0'] = np.float64(y0)

    if self._coeffs['R'] is not None:
      self._max_y = y0 + self._coeffs['R']
      self._min_y = y0 - self._coeffs['R']


  def setR(self, R):

    super(ZSquarePrism, self).setR(R)

    if self._coeffs['x0'] is not None:
      self._max_x = self._coeffs['x0'] + R
      self._min_x = self._coeffs['x0'] - R

    if self._coeffs['y0'] is not None:
      self._max_y = self._coeffs['y0'] + R
      self._min_y = self._coeffs['y0'] - R


  def evaluate(self, point):

    super(ZSquarePrism, self).evaluate(point)

    x, y, z = point._coords

    Rx = abs(self._coeffs['x0'] - x) - self._coeffs['R']
    Ry = abs(self._coeffs['y0'] - y) - self._coeffs['R']

    return max(Rx, Ry)


  def minSurfaceDist(self, point, direction):

    if self.onSurface(point):
      return 0.

    x, y, z = point._coords
    u, v, w = direction.normalize()

    dist_x = (np.sign(u)*self._coeffs['R'] - x)/u
    dist_y = (np.sign(v)*self._coeffs['R'] - y)/v

    if abs(u) < PARALLEL_TO_AXIS_THRESH:
      dist_x = 0.
    if abs(v) < PARALLEL_TO_AXIS_THRESH:
      dist_y = 0.

    for dist in [dist_x, dist_y]:
      new_x = x + dist*u
      new_y = y + dist*v
      new_z = z + dist*w
      if (abs((abs(new_x) - self._coeffs['R'])) < ON_SURFACE_THRESH and
         abs(new_y) < self._coeffs['R']) or \
         (abs((abs(new_y) - self._coeffs['R'])) < ON_SURFACE_THRESH and
          abs(new_x) < self._coeffs['R']):
        intersect = opencg.Point()
        intersect.setCoords((new_x, new_y, new_z))
        dist = point.distanceToPoint(intersect)
        return dist