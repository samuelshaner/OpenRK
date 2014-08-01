__author__ = 'wbinventor'
__email__ = 'wboyd@mit.edu'


from opencsg.universe import *
from opencsg.surface import *
from opencsg.checkvalue import *



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


  def setNumRings(self, num_rings):

    if not is_integer(num_rings):
      msg = 'Unable to set the number of rings for RadialMesh to {0} ' \
            'since it is not an integer'.format(num_rings)
      raise ValueError(msg)

    if num_rings < 0:
      msg = 'Unable to set the number of rings for RadialMesh to {0} ' \
            'since it is a negative integer'.format(num_rings)
      raise ValueError(msg)

    self._num_rings = num_rings


  def setMaxRadius(self, max_radius):

    if not is_float(max_radius) and not is_integer(max_radius):
      msg = 'Unable to set the max radius for RadialMesh to {0} since ' \
            'it is not an integer or floating point value'.format(max_radius)
      raise ValueError(msg)

    if max_radius < 0.:
      msg = 'Unable to set the max radius for RadialMesh to {0} since ' \
            'it is a negative value'.format(max_radius)
      raise ValueError(msg)

    self._max_radius = np.float64(max_radius)


  def setMinRadius(self, min_radius):

    if not is_float(min_radius) and not is_integer(min_radius):
      msg = 'Unable to set the min radius for RadialMesh to {0} since ' \
            'it is not an integer or floating point value'.format(min_radius)
      raise ValueError(msg)

    if min_radius < 0.:
      msg = 'Unable to set the min radius for RadialMesh to {0} since ' \
            'it is a negative value'.format(min_radius)
      raise ValueError(msg)

    self._min_radius = np.float64(min_radius)


  def setWithOuter(self, with_outer):

    if not isinstance(with_outer, bool):
      msg = 'Unable to set with outer for RadialMesh to {0} since ' \
            'it is not a boolean value'.format(with_outer)
      raise ValueError(msg)

    self._with_outer = with_outer


  def setWithInner(self, with_inner):

    if not isinstance(with_inner, bool):
      msg = 'Unable to set with inner for RadialMesh to {0} since ' \
            'it is not a boolean value'.format(with_inner)
      raise ValueError(msg)

    self._with_inner = with_inner


  def subdivideCell(self, cell, universe=None):

    if not isinstance(cell, Cell):
      msg = 'Unable to subdivide Cell with RadialMesh since {0} is not ' \
            'a Cell'.format(cell)
      raise ValueError(msg)

    # Create container for ZCylinders
    cylinders = list()

    if self._num_rings == 0:
      return cylinders

    cell.findBoundingBox()

    if self._max_radius == -np.float("inf"):
      self._max_radius = cell.getMaxX()

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


  def setNumSectors(self, num_sectors):

    if not is_integer(num_sectors):
      msg = 'Unable to set the number of sectors for SectorMesh to {0} ' \
            'since it is not an integer'.format(num_sectors)
      raise ValueError(msg)

    if num_sectors < 0:
      msg = 'Unable to set the number of rings for SectorMesh to {0} ' \
            'since it is a negative integer'.format(num_sectors)
      raise ValueError(msg)

    self._num_sectors = num_sectors


  def subdivideCell(self, cell, universe=None):

    if not isinstance(cell, Cell):
      msg = 'Unable to subdivide Cell with RadialMesh since {0} ' \
            'is not a Cell'.format(cell)
      raise ValueError(msg)

    cell.findBoundingBox()

    # Initialize an empty list of the new subdivided cells
    new_cells = list()

    if self._num_sectors == 0 or self._num_sectors == 1:
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
      msg = 'Unable to subdivide cells with a SectorMesh since {0} ' \
            'is not a Python tuple/list or NumPy array'.format(cells)
      raise ValueError(msg)

    # Initialize an empty list of new Cells
    new_cells = list()

    if self._num_sectors == 0 or self._num_sectors == 1:
      return new_cells

    for cell in cells:
      new_cells.extend(self.subdivideCell(cell=cell))

      if isinstance(universe, Universe):
        universe.removeCell(cell)

    universe.addCells(new_cells)

    return new_cells


  def subdivideUniverse(self, universe):

    if not isinstance(universe, Universe):
      msg = 'Unable to subdivide Universe with a SectorMesh since {0} ' \
            'is not a Universe'.format(universe)
      raise ValueError(msg)

    cells = universe._cells.values()
    new_cells = self.subdivideCells(cells, universe=universe)
    return new_cells
