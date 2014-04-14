__author__ = 'wbinventor'
__email__ = 'wboyd@mit.edu'


from universe import *
from surface import *
from checkvalue import *


# The types of meshes - distance (1D), area (2D) or volume (3D)
spacing_types = ['1D', '2D']    # Add 3D later

# The types of intervals between mesh cells
interval_types = ['linear', 'logarithmic']


class Mesh(object):

  def __init__(self, cell=None):

    # Initialize Mesh class attributes
    self._cell = None
    self._spacing_type = '1D'
    self._interval_type = 'linear'

    if not cell is None:
      self.setCell(cell)


  def getCell(self):
    return self._cell


  def getSpacingType(self):
    return self._spacing_type


  def getIntervalType(self):
    return self._interval_type


  def setCell(self, cell):

    if not isinstance(cell, Cell):
      exit('Unable to set the cell for mesh since %s is not '
           'a Cell' % str(cell))

    self._cell = cell


  def setIntervalType(self, interval_type='linear'):

    if not is_string(interval_type):
      exit('Unable to set the interval type for Mesh to %s since it '
           'is not a string' % str(interval_type))

    if not interval_type in interval_types:
      exit('Unable to set the interval type for Mesh to %s since it '
           'is not a supported type' % interval_type)

    self._interval_type = interval_type


  def setSpacingType(self, spacing_type='1D'):

    if not is_string(spacing_type):
      exit('Unable to set the spacing type for Mesh to %s since it '
           'is not a string' % str(spacing_type))

    if not spacing_type in spacing_types:
      exit('Unable to set the spacing type for Mesh to %s since it '
           'is not a supported type' % spacing_type)

    self._spacing_type = spacing_type


  def subdivideCell(self):

    if self._cell is None:
      exit('Unable to subdivide Mesh cell since it has not been set')



class RadialMesh(Mesh):

  def __init__(self, cell=None):

    super(RadialMesh, self).__init__(cell=cell)

    # Initialize RadialMesh class attributes
    self._num_rings = 0.
    self._max_radius = -np.float("inf")
    self._min_radius = np.float("inf")
    self._with_outer = False
    self._with_inner = False

    if not cell is None:
      self.setCell(cell)


  def getNumRings(self):
    return self._num_rings


  def getMaxRadius(self):
    return self._max_radius


  def getMinRadius(self):
    return self._min_radius


  def getWithOuter(self):
    return self._with_outer


  def getWithInner(self):
    return self._with_inner


  def setNumRings(self, num_rings):

    if not is_integer(num_rings):
      exit('Unable to set the number of rings for RadialMesh to %s '
           'since it is not an integer' % str(num_rings))

    if num_rings < 0:
      exit('Unable to set the number of rings for RadialMesh to %d '
           'since it is a negative integer' % str(num_rings))

    self._num_rings = num_rings


  def setMaxRadius(self, max_radius):

    if not is_float(max_radius) and not is_integer(max_radius):
      exit('Unable to set the max radius for RadialMesh to %s since '
           'it is not an integer or floating point value' % str(max_radius))

    if max_radius < 0.:
      exit('Unable to set the max radius for RadialMesh to %s since '
           'it is a negative value' % str(max_radius))

    self._max_radius = np.float64(max_radius)


  def setMinRadius(self, min_radius):

    if not is_float(min_radius) and not is_integer(min_radius):
      exit('Unable to set the min radius for RadialMesh to %s since '
           'it is not an integer or floating point value' % str(min_radius))

    if min_radius < 0.:
      exit('Unable to set the min radius for RadialMesh to %s since '
           'it is a negative value' % str(min_radius))

    self._min_radius = np.float64(min_radius)


  def setWithOuter(self, with_outer):

    if not isinstance(with_outer, bool):
      exit('Unable to set with outer for RadialMesh to %s since '
           'it is not a boolean value' % str(with_outer))

    self._with_outer = with_outer


  def setWithInner(self, with_inner):

    if not isinstance(with_inner, bool):
      exit('Unable to set with inner for RadialMesh to %s since '
           'it is not a boolean value' % str(with_inner))

    self._with_inner = with_inner


  def subdivideCell(self):

    super(RadialMesh, self).subdivideCell()

    if self._cell is None:
      exit('Unable to subdivide Cell with RadialMesh since no Cell has '
           'been set')

    if self._max_radius == -np.float("inf"):
      exit('Unable to subdivide Cell ID=%s with RadialMesh since the '
           'maximum radius has not been set' % str(self._cell.getId()))

    if self._min_radius == np.float("inf"):
      exit('Unable to subdivide Cell ID=%s with RadialMesh since the '
           'minimum radius has not been set' % str(self._cell.getId()))

    # Create ZCylinders
    cylinders = list()

    if self._interval_type == 'logarithmic':
      exit('Unable to subdivide Cells using a Radial mesh with '
           'logarithmically spaced radii since that feature is not supported')

    if self._spacing_type is '1D':

      # Equally spaced radii
      radii = np.linspace(self._max_radius, self._min_radius, self._num_rings+1)

    elif self._spacing_type is '2D':

      # Equal area radii
      radii = list()
      area = (self._max_radius**2 - self._min_radius**2) / self._num_rings

      # Initialize successively smaller rings
      radii.append(self._max_radius)

      for i in range(self._num_rings):
        delta_area = radii[-1]**2 - area

        if delta_area <= 0.:
          radii.append(0.)

        else:
          radii.append(np.sqrt(delta_area))


    # Create ZCylinders for each radius
    for radius in radii:
      if radius != 0.:
        cylinders.append(ZCylinder(x0=0., y0=0., R=radius))

    # Initialize an empty list of the new subdivided cells
    new_cells = list()

    # Loop over all rings
    for i in range(self._num_rings):

      min_radius = radii[i+1]

      # Create a clone of this cell for this ring
      clone = self._cell.clone()

      # Add outer bounding Surface to the clone
      clone.addSurface(surface=cylinders[i], halfspace=-1)

      # Add non-trivial inner bounding Surface to the clone
      if min_radius != 0.:
        clone.addSurface(surface=cylinders[i+1], halfspace=+1)

      # Add this clone to the new_cells list
      new_cells.append(clone)


    if self._with_outer:
      clone = self._cell.clone()
      clone.addSurface(surface=cylinders[0], halfspace=+1)
      new_cells.append(clone)

    if self._with_inner:
      clone = self._cell.clone()
      clone.addSurface(surface=cylinders[-1], halfspace=-1)
      new_cells.append(clone)


    # MUST REMOVE SURFACES IF AN EQUIVALENT ALREADY EXISTS IN THE CELL!

    return new_cells


class SectorMesh(Mesh):

  def __init__(self, universe=None, cell=None, num_sectors=None):

    super(SectorMesh, self).__init__(universe=universe, cell=cell)

    # Initialize SectorMesh class attributes
    self._num_sectors = 0.

    if not num_sectors is None:
      self.setNumSectors(num_sectors)


  def getNumSectors(self):
    return self._num_sectors


  def setNumSectors(self, num_sectors):

    if not is_integer(num_sectors):
      exit('Unable to set the number of sectors for SectorMesh to %s '
           'since it is not an integer' % str(num_sectors))

    if num_sectors < 0:
      exit('Unable to set the number of rings for SectorMesh to %d '
           'since it is a negative integer' % str(num_sectors))

    self._num_sectors = num_sectors


  def subdivideUniverse(self):

    super(RadialMesh, self).subdivideUniverse()


  def subdivdeCell(self):

    super(RadialMesh, self).subdivideCell()
