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

  def __init__(self, universe=None, cell=None):

    # Initialize Mesh class attributes
    self._universe = None
    self._cell = None
    self._spacing_type = '1D'
    self._interval_type = 'linear'

    if not universe is None:
      self.setUniverse(universe)

    if not cell is None:
      self.setCell(cell)


  def getUniverse(self):
    return self._universe


  def getCell(self):
    return self._cell


  def getSpacingType(self):
    return self._spacing_type


  def getIntervalType(self):
    return self._interval_type


  def setUniverse(self, universe):

    if not isinstance(universe, Universe):
      exit('Unable to set the universe for mesh since %s is not '
           'a Universe' % str(universe))

    self._universe = universe


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


  def setSpacingType(self, spacing_type='linear'):

    if not is_string(spacing_type):
      exit('Unable to set the spacing type for Mesh to %s since it '
           'is not a string' % str(spacing_type))

    if not spacing_type in spacing_types:
      exit('Unable to set the spacing type for Mesh to %s since it '
           'is not a supported type' % spacing_type)

    self._spacing_type = spacing_type


  def subdivideUniverse(self):

    if self._universe is None:
      exit('Unable to subdivide Mesh universe since it has not been set')


  def subdivideCell(self):

    if self._cell is None:
      exit('Unable to subdivide Mesh cell since it has not been set')



class RadialMesh(Mesh):

  def __init__(self, universe=None, cell=None, num_rings=None):

    super(RadialMesh, self).__init__(universe=universe, cell=cell)

    # Initialize RadialMesh class attributes
    self._num_rings = 0.
    self._max_radius = 0.
    self._min_radius = 0.

    if not num_rings is None:
      self.setNumRings(num_rings)


  def getNumRings(self):
    return self._num_rings


  def getMaxRadius(self):
    return self._max_radius


  def getMinRadius(self):
    return self._min_radius


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


  def subdivideUniverse(self):

    super(RadialMesh, self).subdivideUniverse()

    if self._min_radius == 0:
      exit('Unable ot subdivide Universe ID=%d with RadialMesh since the '
           'maximum radius has not been set' % str(self._universe.getId()))

    if self._min_radius == 0:
      exit('Unable ot subdivide Universe ID=%d with RadialMesh since the '
           'minimum radius has not been set' % str(self._universe.getId()))

    # Create ZCylinders
    radii = None
    cylinders = list()

    # Equally spaced radii
    if self._spacing_type is '1D' :
      if self._interval_type is 'linear':
        radii = np.linspace(self._max_radius, self._min_radius, self._num_rings)

      # Logarithmically spaced radii
      else:
        log_max_radius = np.log(self._max_radius)
        log_min_radius = np.log(self._min_radius)
        radii = np.logspace(log_max_radius, log_min_radius, self._num_rings)

    elif self._spacing_type is '2D':

      # Equal area radii
      if self._interval_type is 'linear':
        area = np.pi * (self._max_radius**2 - self._min_radius**2) / self._num_rings

        # Initialize successively smaller rings
        radii.append(self._max_radius)

        for i in range(self._num_rings-1):
          radii.append(np.sqrt(radii[-1]**2 - (area / np.pi)))

      # Equal log(area) spaced radii
      else:

        exit('Unable to create 2D RadialMesh with logarithmically spaced rings')

    # Create ZCylinders for each radius
    for radius in radii:
      if radius != 0.:
        cylinders.append(ZCylinder(x0=0., y0=0., R=radius))

    # Retrieve the Cells from the Universe
    cells = self._universe.getCells()

    # Initialize an empty list of the new subdivided cells
    new_cells = list()

    for cell_id in cells:
      cell = cells[cell_id]

      # Loop over all rings
      for i in range(self._num_rings-1):

        min_radius = radii[i+1]

        # Create a clone of this cell for this ring
        clone = cell.clone()

        # Add outer bounding Surface to the clone
        clone.addSurface(surface=cylinders[i], halfspace=-1)

        # Add non-trivial inner bounding Surface to the clone
        if min_radius != 0.:
          clone.addSurface(surface=cylinders[i+1], halfspace=+1)

        # Add this clone to the new_cells list
        new_cells.append(clone)

        # MUST REMOVE SURFACES IF AN EQUIVALENT ALREADY EXISTS IN THE CELL!




#      for cylinder



  def subdivdeCell(self):

    super(RadialMesh, self).subdivideCell()



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
