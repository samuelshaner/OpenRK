__author__ = 'Will Boyd'
__email__ = 'wboyd@mit.edu'


from point import Point
from checkvalue import *
import numpy as np


# Threshold for determining how close a point must be to a surface to be on it
on_surface_thresh = 1e-12

# A list of all IDs for all Surfaces created
surface_ids = list()

# A static variable for auto-generated Surface IDs
auto_surface_id = 10000

# The Surface boundary conditions
boundary_types = ['interface', 'vacuum', 'reflective']

# The types of Surfaces
surf_types = ['plane',
              'x-plane',
              'y-plane',
              'z-plane',
              'x-cylinder',
              'y-cylinder',
              'z-cylinder',
              'sphere']

max_float = np.finfo(np.float64).max
min_float = np.finfo(np.float64).min


class Surface(object):

  def __init__(self, surface_id=None, name='', boundary='interface'):

    # Initialize class attributes
    self._id = None
    self._name = ''
    self._type = ''
    self._boundary_type = ''

    # A dictionary of the quadratic surface coefficients
    # Key   - coefficient name
    # Value - coefficient value
    self._coeffs = dict()

    # Max/min values
    self._max_x = np.finfo(np.float64).max
    self._max_y = np.finfo(np.float64).max
    self._max_z = np.finfo(np.float64).max
    self._min_x = np.finfo(np.float64).min
    self._min_y = np.finfo(np.float64).min
    self._min_z = np.finfo(np.float64).min

    self.setId(surface_id)
    self.setName(name)
    self.setBoundaryType(boundary)


  def getId(self):
    return self._id


  def getName(self):
    return self._name


  def getType(self):
    return self._type


  def getBoundaryType(self):
    return self._boundary_type


  def getCoeffs(self):
    return self._coeffs


  def getCoeff(self, coeff):

    if coeff in self._coeffs.keys():
      return self._coeffs[coeff]

    else:
      exit('Unable to return the coeff %s for Surface ID=%d since it '
           'does not contain that coefficient' % (str(coeff), self._id))


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
    if not isinstance(point, Point):
      exit('Unable to evaluate point for Surface ID=%d since the input '
           'is not a Point object' % self._id)


  def onSurface(self, point):
    if not isinstance(point, Point):
      exit('Unable to determine whether a point is on Surface ID=%d since '
           'the input is not a Point object' % self._id)

    if np.abs(self.evaluate(point)) < on_surface_thresh:
      return True
    else:
      return False


  def setId(self, surface_id=None):

    global surface_ids

    if surface_id is None:
      global auto_surface_id
      self._id = auto_surface_id
      surface_ids.append(auto_surface_id)
      auto_surface_id += 1

    # Check that the ID is an integer and wasn't already used
    elif is_integer(surface_id):

      # If the Material already has an ID, remove it from global list
      if not self._id is None:
        surface_ids.remove(self._id)

      if surface_id in surface_ids:
        exit('Unable to set Surface ID to %s since a Surface with '
             'this ID was already initialized' % (str(surface_id)))

      if surface_id < 0:
        exit('Unable to set Surface ID to %d since it must be a '
             'non-negative integer' % surface_id)

      else:
        self._id = surface_id
        surface_ids.append(surface_id)

    else:
      exit('Unable to set a non-integer Surface ID %s' % str(surface_id))


  def setName(self, name):

    if not is_string(name):
      exit('Unable to set name for Surface ID=%d with a non-string '
           'value %s' % (self._id, str(name)))

    else:
      self._name = name


  def setBoundaryType(self, boundary):

    if not is_string(boundary):
      exit('Unable to set boundary type for Surface ID=%d with a '
           'non-string value %s' % (self._id, str(boundary)))

    elif not boundary in boundary_types:
      exit('Unable to set boundary type for Surface ID=%d to %s which '
           'is not interface, vacuum or reflective' % (self._id, boundary))

    else:
      self._boundary_type = boundary


  def toString(self):

    string = ''

    string += 'Surface\n'

    surface_id = '{0: <16}'.format('\tID') + '=\t' + str(self._id)
    string += surface_id + '\n'

    name = '{0: <16}'.format('\tName') + '=\t' + self._name
    string += name + '\n'

    type = '{0: <16}'.format('\tType') + '=\t' + self._type
    string += type + '\n'

    boundary_type = '{0: <16}'.format('\tBoundary')
    boundary_type += '=\t' + self._boundary_type
    string += boundary_type + '\n'

    coeffs = '{0: <16}'.format('\tCoefficients') + '\n'

    if len(self._coeffs) > 0:

      for coeff in self._coeffs.keys():
        coeffs += '{0: <16}'.format('\t%s' % coeff) + '=\t'
        coeffs += '{0: <12}'.format(str(self._coeffs[coeff])) + '\n'

    string += coeffs

    return string


  def printString(self):
    print(self.toString())


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

    if not A is None:
      self.setA(A)

    if not B is None:
      self.setB(B)

    if not C is None:
      self.setC(C)

    if not D is None:
      self.setD(D)


  def getA(self):
    return self._coeffs['A']


  def getB(self):
    return self._coeffs['B']


  def getC(self):
    return self._coeffs['C']


  def getD(self):
    return self._coeffs['D']


  def setA(self, A):

    if not is_integer(A) and not is_float(A):
      exit('Unable to set A coefficient for Plane ID=%d to a non-integer '
           'or floating point value %s' % (self._id, str(A)))

    self._coeffs['A'] = np.float64(A)


  def setB(self, B):

    if not is_integer(B) and not is_float(B):
      exit('Unable to set B coefficient for Plane ID=%d to a non-integer '
           'or floating point value %s' % (self._id, str(B)))

    self._coeffs['B'] = np.float64(B)


  def setC(self, C):

    if not is_integer(C) and not is_float(C):
      exit('Unable to set C coefficient for Plane ID=%d to a non-integer '
           'or floating point value %s' % (self._id, str(C)))

    self._coeffs['C'] = np.float64(C)


  def setD(self, D):

    if not is_integer(D) and not is_float(D):
      exit('Unable to set D coefficient for Plane ID=%d to a non-integer '
           'or floating point value %s' % (self._id, str(D)))

    self._coeffs['D'] = np.float64(D)


  def evaluate(self, point):

    super(Plane, self).evaluate(point)

    x, y, z = point.getX(), point.getY(), point.getZ()

    value = self._coeffs['A'] * x + \
            self._coeffs['B'] * y + \
            self._coeffs['C'] * z + \
            self._coeffs['D']

    return value



class XPlane(Plane):

  def __init__(self, surface_id=None, name='',
               boundary='interface', x0=None):

    # Initialize XPlane class attributes
    super(XPlane, self).__init__(surface_id, name, boundary,
                                 A=1., B=0., C=0., D=-x0)

    self._type = 'x-plane'
    self._coeffs['x0'] = None

    if not x0 is None:
      self.setX0(x0)


  def getX0(self):
    return self._coeffs['x0']


  def getMaxX(self, halfspace=None):

    if halfspace is None:
      return self._max_x

    else:

      if not is_integer(halfspace):
        exit('Unable to get the maximum x-coordinate for Surface ID=%d since '
             'the halfspace is a non-integer value %s',
             self._id, str(halfspace))

      if not halfspace in [-1, +1]:
        exit('Unable to get the maximum x-coordinate for Surface ID=%d since '
             'the halfspace is %s which is not +/-1', self._id, str(halfspace))

      if halfspace == -1:
        return self._max_x

      else:
        return np.finfo(np.float64).max


  def getMinX(self, halfspace=None):

    if halfspace is None:
      return self._min_x

    else:

      if not is_integer(halfspace):
        exit('Unable to get the minimum x-coordinate for Surface ID=%d since '
             'the halfspace is a non-integer value %s',
             self._id, str(halfspace))

      if not halfspace in [-1, +1]:
        exit('Unable to get the minimum x-coordinate for Surface ID=%d since '
             'the halfspace is %s which is not +/-1', self._id, str(halfspace))

      if halfspace == -1:
        return np.finfo(np.float64).min

      else:
        return self._min_x


  def setX0(self, x0):

    if not is_integer(x0) and not is_float(x0):
      exit('Unable to set x0 coefficient for XPlane ID=%d to a non-integer '
           'or floating point value %s' % (self._id, str(x0)))

    self._coeffs['x0'] = np.float64(x0)
    self.setD(-x0)
    self._max_x = np.float64(x0)
    self._min_x = np.float64(x0)



class YPlane(Plane):

  def __init__(self, surface_id=None, name='',
               boundary='interface', y0=None):

    # Initialize YPlane class attributes
    super(YPlane, self).__init__(surface_id, name, boundary,
                                 A=0., B=1., C=0., D=-y0)

    self._type = 'y-plane'
    self._coeffs['y0'] = None

    if not y0 is None:
      self.setY0(y0)


  def getY0(self):
    return self._coeffs['y0']

  def getMaxY(self, halfspace=None):

    if halfspace is None:
      return self._max_y

    else:

      if not is_integer(halfspace):
        exit('Unable to get the maximum y-coordinate for Surface ID=%d since '
             'the halfspace is a non-integer value %s',
             self._id, str(halfspace))

      if not halfspace in [-1, +1]:
        exit('Unable to get the maximum y-coordinate for Surface ID=%d since '
             'the halfspace is %s which is not +/-1', self._id, str(halfspace))

      if halfspace == -1:
        return self._max_y

      else:
        return np.finfo(np.float64).max


  def getMinY(self, halfspace=None):

    if halfspace is None:
      return self._min_y

    else:

      if not is_integer(halfspace):
        exit('Unable to get the minimum y-coordinate for XPlane ID=%d since '
             'the halfspace is a non-integer value %s',
             self._id, str(halfspace))

      if not halfspace in [-1, +1]:
        exit('Unable to get the minimum y-coordinate for XPlane ID=%d since '
             'the halfspace is %s which is not +/-1', self._id, str(halfspace))

      if halfspace == -1:
        return np.finfo(np.float64).min

      else:
        return self._min_y


  def setY0(self, y0):

    if not is_integer(y0) and not is_float(y0):
      exit('Unable to set y0 coefficient for YPlane ID=%d to a non-integer '
           'or floating point value %s' % (self._id, str(y0)))

    self._coeffs['y0'] = np.float64(y0)
    self.setD(-y0)
    self._max_y = np.float64(y0)
    self._min_y = np.float64(y0)


class ZPlane(Plane):

  def __init__(self, surface_id=None, name='',
               boundary='interface', z0=None):

    # Initialize ZPlane class attributes
    super(ZPlane, self).__init__(surface_id, name, boundary,
                                 A=0., B=0., C=1., D=-z0)

    self._type = 'z-plane'
    self._coeffs['z0'] = None

    if not z0 is None:
      self.setZ0(z0)


  def getZ0(self):
    return self._coeffs['z0']


  def getMaxZ(self, halfspace=None):

    if halfspace is None:
      return self._max_z

    else:

      if not is_integer(halfspace):
        exit('Unable to get the maximum z-coordinate for YPlane ID=%d since '
             'the halfspace is a non-integer value %s',
             self._id, str(halfspace))

      if not halfspace in [-1, +1]:
        exit('Unable to get the maximum z-coordinate for YPlane ID=%d since '
             'the halfspace is %s which is not +/-1', self._id, str(halfspace))

      if halfspace == -1:
        return self._max_z

      else:
        return np.finfo(np.float64).max


  def getMinZ(self, halfspace=None):

    if halfspace is None:
      return self._min_z

    else:

      if not is_integer(halfspace):
        exit('Unable to get the minimum z-coordinate for ZPlane ID=%d since '
             'the halfspace is a non-integer value %s',
             self._id, str(halfspace))

      if not halfspace in [-1, +1]:
        exit('Unable to get the minimum z-coordinate for ZPlane ID=%d since '
             'the halfspace is %s which is not +/-1', self._id, str(halfspace))

      if halfspace == -1:
        return np.finfo(np.float64).min

      else:
        return self._min_z


  def setZ0(self, z0):

    if not is_integer(z0) and not is_float(z0):
      exit('Unable to set z0 coefficient for ZPlane ID=%d to a non-integer '
           'or floating point value %s' % (self._id, str(z0)))

    self._coeffs['z0'] = np.float64(z0)
    self.setD(-z0)
    self._max_z = z0
    self._min_z = z0



class Cylinder(Surface):

  def __init__(self, surface_id=None, name='',
               boundary='interface', R=None):

    # Initialize Cylinder class attributes
    super(Cylinder, self).__init__(surface_id, name, boundary)

    self._coeffs['R'] = None


  def getR(self):
    return self._coeffs['R']


  def setR(self, R):

    if not is_integer(R) and not is_float(R):
      exit('Unable to set R coefficient for Cylinder ID=%d to a non-integer '
           'or floating point value %s' % (self._id, str(R)))

    self._coeffs['R'] = np.float64(R)



class XCylinder(Cylinder):

  def __init__(self, surface_id=None, name='',
               boundary='interface', y0=None, z0=None, R=None):

    # Initialize XCylinder class attributes
    super(XCylinder, self).__init__(surface_id, name, boundary, R)

    self._type = 'x-cylinder'
    self._coeffs['y0'] = None
    self._coeffs['z0'] = None

    if not y0 is None:
      self.setY0(y0)

    if not z0 is None:
      self.setZ0(z0)

    if not R is None:
      self.setR(R)


  def getY0(self):
    return self._coeffs['y0']


  def getZ0(self):
    return self._coeffs['z0']


  def getMaxY(self, halfspace=None):

    if halfspace is None:
      return self._max_y

    else:

      if not is_integer(halfspace):
        exit('Unable to get the maximum y-coordinate for XClinder ID=%d since '
             'the halfspace is a non-integer value %s',
             self._id, str(halfspace))

      if not halfspace in [-1, +1]:
        exit('Unable to get the maximum y-coordinate for XCylinder ID=%d since '
             'the halfspace is %s which is not +/-1', self._id, str(halfspace))

      if halfspace == -1:
        return self._max_y

      else:
        return np.finfo(np.float64).max


  def getMinY(self, halfspace=None):

    if halfspace is None:
      return self._min_y

    else:

      if not is_integer(halfspace):
        exit('Unable to get the minimum y-coordinate for XCylinder ID=%d since '
             'the halfspace is a non-integer value %s',
             self._id, str(halfspace))

      if not halfspace in [-1, +1]:
        exit('Unable to get the minimum y-coordinate for XCylinder ID=%d since '
             'the halfspace is %s which is not +/-1', self._id, str(halfspace))

      if halfspace == -1:
        return self._min_y

      else:
        return np.finfo(np.float64).min


  def getMaxZ(self, halfspace=None):

    if halfspace is None:
      return self._max_z

    else:

      if not is_integer(halfspace):
        exit('Unable to get the maximum z-coordinate for XClinder ID=%d since '
             'the halfspace is a non-integer value %s',
             self._id, str(halfspace))

      if not halfspace in [-1, +1]:
        exit('Unable to get the maximum z-coordinate for XCylinder ID=%d since '
             'the halfspace is %s which is not +/-1', self._id, str(halfspace))

      if halfspace == -1:
        return self._max_z

      else:
        return np.finfo(np.float64).max


  def getMinZ(self, halfspace=None):

    if halfspace is None:
      return self._min_z

    else:

      if not is_integer(halfspace):
        exit('Unable to get the minimum z-coordinate for XCylinder ID=%d since '
             'the halfspace is a non-integer value %s',
             self._id, str(halfspace))

      if not halfspace in [-1, +1]:
        exit('Unable to get the minimum z-coordinate for XCylinder ID=%d since '
             'the halfspace is %s which is not +/-1', self._id, str(halfspace))

      if halfspace == -1:
        return self._min_z

      else:
        return np.finfo(np.float64).min


  def setY0(self, y0):

    if not is_integer(y0) and not is_float(y0):
      exit('Unable to set y0 coefficient for XCylinder ID=%d to a non-integer '
           'or floating point value %s' % (self._id, str(y0)))

    self._coeffs['y0'] = np.float64(y0)

    if not self._coeffs['R'] is None:
      self._max_y = y0 + self._coeffs['R']
      self._min_y = y0 - self._coeffs['R']


  def setZ0(self, z0):

    if not is_integer(z0) and not is_float(z0):
      exit('Unable to set z0 coefficient for XCylinder ID=%d to a non-integer '
           'or floating point value %s' % (self._id, str(z0)))

    self._coeffs['z0'] = np.float64(z0)

    if not self._coeffs['R'] is None:
      self._max_z = z0 + self._coeffs['R']
      self._min_z = z0 - self._coeffs['R']


  def setR(self, R):

    super(XCylinder, self).setR(R)

    if not self._coeffs['y0'] is None:
      self._max_y = self._coeffs['y0'] + self._coeffs['R']
      self._min_y = self._coeffs['y0'] - self._coeffs['R']

    if not self._coeffs['z0'] is None:
      self._max_z = self._coeffs['z0'] + self._coeffs['R']
      self._min_z = self._coeffs['z0'] - self._coeffs['R']


  def evaluate(self, point):

    super(XCylinder, self).evaluate(point)

    coords = point.getCoords()
    r2 = (coords[1] - self._coeffs['y0'])**2 + \
        (coords[2] - self._coeffs['z0'])**2
    r = np.sqrt(r2)
    return (r - self._coeffs['R'])


class YCylinder(Cylinder):

  def __init__(self, surface_id=None, name='',
               boundary='interface', x0=None, z0=None, R=None):

    # Initialize YCylinder class attributes
    super(YCylinder, self).__init__(surface_id, name, boundary, R)

    self._type = 'y-cylinder'
    self._coeffs['x0'] = None
    self._coeffs['z0'] = None

    if not x0 is None:
      self.setX0(x0)

    if not z0 is None:
      self.setZ0(z0)

    if not R is None:
      self.setR(R)


  def getX0(self):
    return self._coeffs['x0']


  def getZ0(self):
    return self._coeffs['z0']


  def getMaxX(self, halfspace=None):

    if halfspace is None:
      return self._max_x

    else:

      if not is_integer(halfspace):
        exit('Unable to get the maximum x-coordinate for YClinder ID=%d since '
             'the halfspace is a non-integer value %s',
             self._id, str(halfspace))

      if not halfspace in [-1, +1]:
        exit('Unable to get the maximum x-coordinate for YCylinder ID=%d since '
             'the halfspace is %s which is not +/-1', self._id, str(halfspace))

      if halfspace == -1:
        return self._max_x

      else:
        return np.finfo(np.float64).max


  def getMinX(self, halfspace=None):

    if halfspace is None:
      return self._min_x

    else:

      if not is_integer(halfspace):
        exit('Unable to get the minimum x-coordinate for YCylinder ID=%d since '
             'the halfspace is a non-integer value %s',
             self._id, str(halfspace))

      if not halfspace in [-1, +1]:
        exit('Unable to get the minimum x-coordinate for YCylinder ID=%d since '
             'the halfspace is %s which is not +/-1', self._id, str(halfspace))

      if halfspace == -1:
        return self._min_x

      else:
        return np.finfo(np.float64).min


  def getMaxZ(self, halfspace=None):

    if halfspace is None:
      return self._max_z

    else:

      if not is_integer(halfspace):
        exit('Unable to get the maximum z-coordinate for YClinder ID=%d since '
             'the halfspace is a non-integer value %s',
             self._id, str(halfspace))

      if not halfspace in [-1, +1]:
        exit('Unable to get the maximum z-coordinate for YCylinder ID=%d since '
             'the halfspace is %s which is not +/-1', self._id, str(halfspace))

      if halfspace == -1:
        return self._max_z

      else:
        return np.finfo(np.float64).max


  def getMinZ(self, halfspace=None):

    if halfspace is None:
      return self._min_z

    else:

      if not is_integer(halfspace):
        exit('Unable to get the minimum z-coordinate for YCylinder ID=%d since '
             'the halfspace is a non-integer value %s',
             self._id, str(halfspace))

      if not halfspace in [-1, +1]:
        exit('Unable to get the minimum z-coordinate for YCylinder ID=%d since '
             'the halfspace is %s which is not +/-1', self._id, str(halfspace))

      if halfspace == -1:
        return self._min_z

      else:
        return np.finfo(np.float64).min


  def setX0(self, x0):

    if not is_integer(x0) and not is_float(x0):
      exit('Unable to set x0 coefficient for YCylinder ID=%d to a non-integer '
           'or floating point value %s' % (self._id, str(x0)))

    self._coeffs['x0'] = np.float64(x0)

    if not self._coeffs['R'] is None:
      self._max_x = x0 + self._coeffs['R']
      self._min_x = x0 - self._coeffs['R']


  def setZ0(self, z0):

    if not is_integer(z0) and not is_float(z0):
      exit('Unable to set z0 coefficient for YCylinder ID=%d to a non-integer '
           'or floating point value %s' % (self._id, str(z0)))

    self._coeffs['z0'] = np.float64(z0)

    if not self._coeffs['R'] is None:
      self._max_z = z0 + self._coeffs['R']
      self._min_z = z0 - self._coeffs['R']


  def setR(self, R):

    super(YCylinder, self).setR(R)

    if not self._coeffs['x0'] is None:
      self._max_x = self._coeffs['x0'] + self._coeffs['R']
      self._min_x = self._coeffs['x0'] - self._coeffs['R']

    if not self._coeffs['z0'] is None:
      self._max_z = self._coeffs['z0'] + self._coeffs['R']
      self._min_z = self._coeffs['z0'] - self._coeffs['R']


  def evaluate(self, point):

    super(YCylinder, self).evaluate(point)

    coords = point.getCoords()
    r2 = (coords[0] - self._coeffs['x0'])**2 + \
        (coords[2] - self._coeffs['z0'])**2
    r = np.sqrt(r2)
    return (r - self._coeffs['R'])


class ZCylinder(Cylinder):

  def __init__(self, surface_id=None, name='',
               boundary='interface', x0=None, y0=None, R=None):

    # Initialize ZCylinder class attributes
    super(ZCylinder, self).__init__(surface_id, name, boundary, R)

    self._type = 'z-cylinder'
    self._coeffs['x0'] = None
    self._coeffs['y0'] = None

    if not x0 is None:
      self.setX0(x0)

    if not y0 is None:
      self.setY0(y0)

    if not R is None:
      self.setR(R)


  def getX0(self):
    return self._coeffs['x0']

  def getY0(self):
    return self._coeffs['y0']


  def getMaxX(self, halfspace=None):

    if halfspace is None:
      return self._max_x

    else:

      if not is_integer(halfspace):
        exit('Unable to get the maximum x-coordinate for ZCylinder ID=%d since '
             'the halfspace is a non-integer value %s',
             self._id, str(halfspace))

      if not halfspace in [-1, +1]:
        exit('Unable to get the maximum x-coordinate for ZCylinder ID=%d since '
             'the halfspace is %s which is not +/-1', self._id, str(halfspace))

      if halfspace == -1:
        return self._max_x

      else:
        return np.finfo(np.float64).max


  def getMinX(self, halfspace=None):

    if halfspace is None:
      return self._min_x

    else:

      if not is_integer(halfspace):
        exit('Unable to get the minimum x-coordinate for ZCylinder ID=%d since '
             'the halfspace is a non-integer value %s',
             self._id, str(halfspace))

      if not halfspace in [-1, +1]:
        exit('Unable to get the minimum x-coordinate for ZCylinder ID=%d since '
             'the halfspace is %s which is not +/-1', self._id, str(halfspace))

      if halfspace == -1:
        return self._min_x

      else:
        return np.finfo(np.float64).min


  def getMaxY(self, halfspace=None):

    if halfspace is None:
      return self._max_y

    else:

      if not is_integer(halfspace):
        exit('Unable to get the maximum y-coordinate for ZClinder ID=%d since '
             'the halfspace is a non-integer value %s',
             self._id, str(halfspace))

      if not halfspace in [-1, +1]:
        exit('Unable to get the maximum y-coordinate for ZCylinder ID=%d since '
             'the halfspace is %s which is not +/-1', self._id, str(halfspace))

      if halfspace == -1:
        return self._max_y

      else:
        return np.finfo(np.float64).max


  def getMinY(self, halfspace=None):

    if halfspace is None:
      return self._min_y

    else:

      if not is_integer(halfspace):
        exit('Unable to get the minimum y-coordinate for ZCylinder ID=%d since '
             'the halfspace is a non-integer value %s',
             self._id, str(halfspace))

      if not halfspace in [-1, +1]:
        exit('Unable to get the minimum y-coordinate for ZCylinder ID=%d since '
             'the halfspace is %s which is not +/-1', self._id, str(halfspace))

      if halfspace == -1:
        return self._min_y

      else:
        return np.finfo(np.float64).min


  def setX0(self, x0):

    if not is_integer(x0) and not is_float(x0):
      exit('Unable to set x0 coefficient for ZCylinder ID=%d to a '
         'non-integer value %s' % (self._id, str(x0)))

    self._coeffs['x0'] = np.float64(x0)

    if not self._coeffs['R'] is None:
      self._max_x = x0 + self._coeffs['R']
      self._min_x = x0 - self._coeffs['R']


  def setY0(self, y0):

    if not is_integer(y0) and not is_float(y0):
      exit('Unable to set y0 coefficient for ZCylinder ID=%d to a '
         'non-integer value %s' % (self._id, str(y0)))

    self._coeffs['y0'] = np.float64(y0)

    if not self._coeffs['R'] is None:
      self._max_y = y0 + self._coeffs['R']
      self._min_y = y0 - self._coeffs['R']


  def setR(self, R):

    super(ZCylinder, self).setR(R)

    if not self._coeffs['x0'] is None:
      self._max_x = self._coeffs['x0'] + self._coeffs['R']
      self._min_x = self._coeffs['x0'] - self._coeffs['R']

    if not self._coeffs['y0'] is None:
      self._max_y = self._coeffs['y0'] + self._coeffs['R']
      self._min_y = self._coeffs['y0'] - self._coeffs['R']


  def evaluate(self, point):

    super(ZCylinder, self).evaluate(point)

    coords = point.getCoords()
    r2 = (coords[0] - self._coeffs['x0'])**2 + \
        (coords[1] - self._coeffs['y0'])**2
    r = np.sqrt(r2)
    return (r - self._coeffs['R'])



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

    if not x0 is None:
      self.setX0(x0)

    if not y0 is None:
      self.setY0(y0)

    if not z0 is None:
      self.setZ0(z0)

    if not R is None:
      self.setZ0(R)


  def getX0(self):
    return self._coeffs['x0']


  def getY0(self):
    return self._coeffs['y0']


  def getZ0(self):
    return self._coeffs['z0']


  def getR(self):
    return self._coeffs['R']


  def getMaxX(self, halfspace=None):

    if halfspace is None:
      return self._max_x

    else:

      if not is_integer(halfspace):
        exit('Unable to get the maximum x-coordinate for Sphere ID=%d since '
             'the halfspace is a non-integer value %s',
             self._id, str(halfspace))

      if not halfspace in [-1, +1]:
        exit('Unable to get the maximum x-coordinate for Sphere ID=%d since '
             'the halfspace is %s which is not +/-1', self._id, str(halfspace))

      if halfspace == -1:
        return self._max_x

      else:
        return np.finfo(np.float64).max


  def getMinX(self, halfspace=None):

    if halfspace is None:
      return self._min_x

    else:

      if not is_integer(halfspace):
        exit('Unable to get the minimum x-coordinate for Sphere ID=%d since '
             'the halfspace is a non-integer value %s',
             self._id, str(halfspace))

      if not halfspace in [-1, +1]:
        exit('Unable to get the minimum x-coordinate for Sphere ID=%d since '
             'the halfspace is %s which is not +/-1', self._id, str(halfspace))

      if halfspace == -1:
        return self._min_x

      else:
        return np.finfo(np.float64).min


  def getMaxY(self, halfspace=None):

    if halfspace is None:
      return self._max_y

    else:

      if not is_integer(halfspace):
        exit('Unable to get the maximum y-coordinate for Sphere ID=%d since '
             'the halfspace is a non-integer value %s',
             self._id, str(halfspace))

      if not halfspace in [-1, +1]:
        exit('Unable to get the maximum y-coordinate for Sphere ID=%d since '
             'the halfspace is %s which is not +/-1', self._id, str(halfspace))

      if halfspace == -1:
        return self._max_y

      else:
        return np.finfo(np.float64).max


  def getMinY(self, halfspace=None):

    if halfspace is None:
      return self._min_y

    else:

      if not is_integer(halfspace):
        exit('Unable to get the minimum y-coordinate for Sphere ID=%d since '
             'the halfspace is a non-integer value %s',
             self._id, str(halfspace))

      if not halfspace in [-1, +1]:
        exit('Unable to get the minimum y-coordinate for Sphere ID=%d since '
             'the halfspace is %s which is not +/-1', self._id, str(halfspace))

      if halfspace == -1:
        return self._min_y

      else:
        return np.finfo(np.float64).min


  def getMaxZ(self, halfspace=None):

    if halfspace is None:
      return self._max_z

    else:

      if not is_integer(halfspace):
        exit('Unable to get the maximum z-coordinate for Sphere ID=%d since '
             'the halfspace is a non-integer value %s',
             self._id, str(halfspace))

      if not halfspace in [-1, +1]:
        exit('Unable to get the maximum z-coordinate for Sphere ID=%d since '
             'the halfspace is %s which is not +/-1', self._id, str(halfspace))

      if halfspace == -1:
        return self._max_z

      else:
        return np.finfo(np.float64).max


  def getMinZ(self, halfspace=None):

    if halfspace is None:
      return self._min_z

    else:

      if not is_integer(halfspace):
        exit('Unable to get the minimum z-coordinate for Sphere ID=%d since '
             'the halfspace is a non-integer value %s',
             self._id, str(halfspace))

      if not halfspace in [-1, +1]:
        exit('Unable to get the minimum z-coordinate for Sphere ID=%d since '
             'the halfspace is %s which is not +/-1', self._id, str(halfspace))

      if halfspace == -1:
        return self._min_z

      else:
        return np.finfo(np.float64).min


  def setX0(self, x0):

    if not is_integer(x0) and not is_float(x0):
      exit('Unable to set x0 coefficient for Sphere ID=%d to a non-integer '
           'or floating point value %s' % (self._id, str(x0)))

    self._coeffs['x0'] = np.float64(x0)

    if not self._coeffs['R'] is None:
      self._max_x = x0 + self._coeffs['R']
      self._min_x = x0 - self._coeffs['R']


  def setY0(self, y0):

    if not is_integer(y0) and not is_float(y0):
      exit('Unable to set y0 coefficient for Sphere ID=%d to a non-integer '
           'or floating point value %s' % (self._id, str(y0)))

    self._coeffs['y0'] = np.float64(y0)

    if not self._coeffs['R'] is None:
      self._max_y = y0 + self._coeffs['R']
      self._min_y = y0 - self._coeffs['R']


  def setZ0(self, z0):

    if not is_integer(z0) and not is_float(z0):
      exit('Unable to set z0 coefficient for Sphere ID=%d to a non-integer '
           'or floating point value %s' % (self._id, str(z0)))

    self._coeffs['z0'] = np.float64(z0)

    if not self._coeffs['R'] is None:
      self._max_z = z0 + self._coeffs['R']
      self._min_z = z0 - self._coeffs['R']


  def setR(self, R):

    if not is_integer(R) and not is_float(R):
      exit('Unable to set R coefficient for Sphere ID=%d to a non-integer '
           'or floating point value %s' % (self._id, str(R)))

    self._coeffs['R'] = np.float64(R)

    if not self._coeffs['x0'] is None:
      self._max_x = self._coeffs['x0'] + R
      self._min_x = self._coeffs['x0'] - R

    if not self._coeffs['y0'] is None:
      self._max_y = self._coeffs['y0'] + R
      self._min_y = self._coeffs['y0'] - R

    if not self._coeffs['z0'] is None:
      self._max_z = self._coeffs['z0'] + R
      self._min_z = self._coeffs['z0'] - R


  def evaluate(self, point):

    super(Sphere, self).evaluate(point)

    coords = point.getCoords()

    R2 = (self._coeffs['x0'] - coords[0])**2 + \
         (self._coeffs['y0'] - coords[1])**2 + \
         (self._coeffs['z0'] - coords[2])**2
    R = np.sqrt(R2)
    return (R - self._coeffs['R'])