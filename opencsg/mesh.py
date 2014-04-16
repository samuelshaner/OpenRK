__author__ = 'wbinventor'
__email__ = 'wboyd@mit.edu'


from universe import *
from surface import *
from checkvalue import *



class RadialMesh(object):

  def __init__(self, num_rings=None):

    super(RadialMesh, self).__init__()

    # Initialize RadialMesh class attributes
    self._num_rings = 0.
    self._max_radius = -np.float("inf")
    self._min_radius = 0.
    self._with_outer = False
    self._with_inner = False

    if not num_rings is None:
      self.setNumRings(num_rings)


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


  def subdivideCell(self, cell, universe=None):

    if not isinstance(cell, Cell):
      exit('Unable to subdivide Cell with RadialMesh since %s is not '
           'a Cell' % str(cell))

    cell.findBoundingBox()

    if self._max_radius == -np.float("inf"):
      self._max_radius = cell.getMaxX()

    # Create container for ZCylinders
    cylinders = list()

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
      cylinders.append(ZCylinder(x0=0., y0=0., R=radius))

    # Initialize an empty list of the new subdivided cells
    new_cells = list()

    # Loop over all rings
    for i in range(self._num_rings):

      min_radius = radii[i+1]

      # Create a clone of this cell for this ring
      clone = cell.clone()

      # Add outer bounding Surface to the clone
      clone.addSurface(surface=cylinders[i], halfspace=-1)

      # Add non-trivial inner bounding Surface to the clone
      if abs(min_radius) > 1e-5:
        clone.addSurface(surface=cylinders[i+1], halfspace=+1)

      # Add this clone to the new_cells list
      new_cells.append(clone)


    if self._with_outer:
      clone = cell.clone()
      clone.addSurface(surface=cylinders[0], halfspace=+1)
      new_cells.append(clone)

    if self._with_inner:
      clone = cell.clone()
      clone.addSurface(surface=cylinders[-1], halfspace=-1)
      new_cells.append(clone)

    for new_cell in new_cells:
      new_cell.findBoundingBox()
      new_cell.removeRedundantSurfaces()

    if isinstance(universe, Universe):
      universe.removeCell(cell)
      universe.addCells(new_cells)

    return new_cells


class SectorMesh(object):

  def __init__(self, num_sectors=None):

    super(SectorMesh, self).__init__()

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


  def subdivideCell(self, cell, universe=None):

    if not isinstance(cell, Cell):
      exit('Unable to subdivide Cell with RadialMesh since %s is not '
           'a Cell' % str(cell))

    cell.findBoundingBox()

    # Initialize an empty list of the new subdivided cells
    new_cells = list()

    if self._num_sectors == 0:
      return new_cells

    # Initialize an empty list for Planes
    planes = list()

    delta_azim = 2. * np.pi / self._num_sectors

    # Create each of the bounding planes for the sector Cells
    for i in range(self._num_sectors):

      # Compute the angle for this plane
      azim_angle = i * delta_azim

      # Instantiate the plane
      A = np.cos(azim_angle)
      B = np.sin(azim_angle)
      planes.append(Plane(A=A, B=B, C=0., D=0.))


    # Create sectors using disjoint halfspaces of pairing Planes
    for i in range(self._num_sectors):

      # Create new Cell clone for this sector Cell
      sector = cell.clone()

      # Add new bounding planar Surfaces to the clone
      sector.addSurface(surface=planes[i], halfspace=+1)

      if self._num_sectors != 2:

        if (i+1) < self._num_sectors:
          sector.addSurface(surface=planes[i+1], halfspace=-1)

        else:
          sector.addSurface(surface=planes[0], halfspace=-1)

      # Store the clone in the container of new sector Cells
      new_cells.append(sector)

    if isinstance(universe, Universe):
      universe.removeCell(cell)
      universe.addCells(new_cells)

    return new_cells


  def subdivideCells(self, cells, universe=None):

    if not isinstance(cells, (tuple, list, np.ndarray)):
      exit('Unable to subdivide cells with a SectorMesh since %s is not '
           'a Python tuple/list or NumPy array' % str(cells))

    # Initialize an empty list of new Cells
    new_cells = list()

    for cell in cells:
      new_cells.extend(self.subdivideCell(cell=cell))

      if isinstance(universe, Universe):
        universe.removeCell(cell)

    universe.addCells(new_cells)

    return new_cells


  def subdivideUniverse(self, universe):

    if not isinstance(universe, Universe):
      exit('Unable to subdivide Universe with a SectorMesh since %s is not '
           'a Universe' % str(universe))

    cells = universe.getCells().values()
    new_cells = self.subdivideCells(cells, universe=universe)
    return new_cells