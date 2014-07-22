__author__ = 'Will Boyd'
__email__ = 'wboyd@mit.edu'


from universe import *
from point import Point


class Geometry(object):

  def __init__(self):

    # Initialize Geometry class attributes
    self._root_universe = None
    self._num_regions = 0
    self._volume = np.float64(0.)

    # A NumPy array of volumes for each region, indexed by Region ID
    self._region_volumes = None


  def buildNeighbors(self):

    # Allocate dictionaries for neighbor Cells for each Surface halfspace
    surface_positive_neighbors = dict()
    surface_negative_neighbors = dict()

    cells = self.getAllCells()

    # Determine the number of Cells sharing each Surface's halfspaces
    for cell_id, cell in cells.iteritems():
      surfaces = cell._surfaces

      for surface_id in surfaces.keys():

        surface = surfaces[surface_id][0]
        halfspace = surfaces[surface_id][1]
        surface.addNeighborCell(cell_id, halfspace)


  def getMaxX(self):

    if self._root_universe is None:
      msg = 'Unable to get the maximum x since the Geometry does not ' \
            'contain the base Universe ID=0'
      raise ValueError(msg)

    return self._root_universe.getMaxX()


  def getMaxY(self):

    if self._root_universe is None:
      msg = 'Unable to get the maximum y since the Geometry does not ' \
            'contain the base Universe ID=0'
      raise ValueError(msg)

    return self._root_universe.getMaxY()


  def getMaxZ(self):

    if self._root_universe is None:
      msg = 'Unable to get the maximum z since the Geometry does not ' \
            'contain the base Universe ID=0'
      raise ValueError(msg)

    return self._root_universe.getMaxZ()


  def getMinX(self):

    if self._root_universe is None:
      msg = 'Unable to get the minimum x since the Geometry does not ' \
            'contain the base Universe ID=0'
      raise ValueError(msg)

    return self._root_universe.getMinX()


  def getMinY(self):

    if self._root_universe is None:
      msg = 'Unable to get the minimum y since the Geometry does not ' \
            'contain the base Universe ID=0'
      raise ValueError(msg)

    return self._root_universe.getMinY()


  def getMinZ(self):

    if self._root_universe is None:
      msg = 'Unable to get the minimum z since the Geometry does not ' \
            'contain the base Universe ID=0'
      raise ValueError(msg)

    return self._root_universe.getMinZ()


  def getBounds(self):

    if self._root_universe is None:
      msg = 'Unable to get the bounds since the Geometry does not ' \
            'contain the base Universe ID=0'
      raise ValueError(msg)

    bounds = [self._root_universe.getMinX(),
              self._root_universe.getMaxX(),
              self._root_universe.getMinY(),
              self._root_universe.getMaxY(),
              self._root_universe.getMinZ(),
              self._root_universe.getMaxZ()]

    return bounds



  def getAllCells(self):

    if self._root_universe is None:
      msg = 'Unable to get all Cells since the Geometry does not ' \
            'contain the base Universe ID=0'
      raise ValueError(msg)

    return self._root_universe.getAllCells()


  def getAllUniverses(self):

    if self._root_universe is None:
      msg = 'Unable to get all Universes since the Geometry does not ' \
            'contain the base Universe ID=0'
      raise ValueError(msg)

    return self._root_universe.getAllUniverses()


  def getAllMaterials(self):

    material_cells = self.getAllMaterialCells()
    materials = set()

    for cell in material_cells:
      materials.add(cell._fill)

    return list(materials)


  def getAllMaterialCells(self):

    all_cells = self.getAllCells()
    material_cells = set()

    for cell_id, cell in all_cells.iteritems():
      if cell._type == 'material':
        material_cells.add(cell)

    return list(material_cells)


  def getAllMaterialUniverses(self):

    all_universes = self.getAllUniverses()
    material_universes = set()

    for universe_id, universe in all_universes.iteritems():

      cells = universe._cells

      for cell_id, cell in cells.iteritems():
        if cell._type == 'material':
          material_universes.add(universe)

    return list(material_universes)


  def setRootUniverse(self, root_universe):

    if not isinstance(root_universe, Universe):
      msg = 'Unable to set the root Universe for the Geometry since {0} is ' \
            'not a Universe'.format(root_universe)
      raise ValueError(msg)

    if not root_universe._id == 0:
      msg = 'Unable to set the root Universe for the Geometry with a ' \
            'Universe with ID={0}. The root Universe must have ' \
            'ID=0.'.format(root_universe._id)
      raise ValueError(msg)

    self._root_universe = root_universe


  def setVolume(self, volume, tolerance=1e-3):

    if not is_float(volume) and not is_integer(volume):
      msg = 'Unable to set the volume of the Geometry to {0} since it is ' \
            'neither an integer nor a floating point value'.format(volume)
      raise ValueError(msg)

    if not is_float(tolerance):
      msg = 'Unable to compute the volume of the Geometry to a tolerance ' \
            'of {0} since it is not a floating point value'.format(tolerance)
      raise ValueError(msg)

    if self._root_universe is None:
      msg = 'Unable to initialize cell offsets since the Geometry does not ' \
            'contain the base Universe ID=0'
      raise ValueError(msg)

    self._volume = np.float64(volume)

    self._root_universe.computeVolumeFractions(volume=self._volume,
                                               tolerance=tolerance)

    # Initialize an array for the region volumes
    self._region_volumes = np.zeros(self._num_regions)

    # Get the region volumes
    for region in range(self._num_regions):
      
      coords = self.findRegion(region)
      coords = coords.getTailNode()

      cell = coords._cell
      volume = cell._volume

      self._region_volumes[region] = volume


  def initializeCellOffsets(self):

    if self._root_universe is None:
      msg = 'Unable to initialize cell offsets since the Geometry ' \
            'does not contain the base Universe ID=0'
      raise ValueError(msg)

    self._root_universe.initializeCellOffsets()
    self._num_regions = self._root_universe._num_regions


  def getRegionId(self, x=0., y=0., z=0.):

    coords = self.findCoords(x=x, y=y, z=z)
    region_id = 0

    while not coords is None:

      # The coords is a UnivCoords object
      if isinstance(coords, UnivCoords):
        universe = coords._universe
        cell = coords._cell
        region_id += universe.getCellOffset(cell)

      # The coords is a LatCoords object
      else:
        lattice = coords._lattice
        lat_x, lat_y = coords._lattice_x, coords._lattice_y
        region_id += lattice.getCellOffset(lat_x, lat_y)

      coords = coords._next

    return region_id


  def findCell(self, x=0., y=0., z=0.):

    if self._root_universe is None:
      msg = 'Unable to find cell since the Geometry does not ' \
            'contain the base Universe ID=0'
      raise ValueError(msg)

    point = Point(x=x, y=y, z=z)
    
    localcoords = UnivCoords(point=point)
    localcoords.setUniverse(self._root_universe)

    return self._root_universe.findCell(localcoords=localcoords)


  def findRegion(self, region_id=0):

    if self._root_universe is None:
      msg = 'Unable to find coords since the Geometry does not ' \
            'contain the base Universe ID=0'
      raise ValueError(msg)

    if not is_integer(region_id):
      msg = 'Unable to get the path for region_id {0} since ' \
            'it is not an integer value'.format(region_id)
      raise ValueError(msg)

    if region_id < 0:
      msg = 'Unable to get the path for region_id={0} since it is ' \
            'a negative integer'.format(region_id)
      raise ValueError(msg)

    localcoords = UnivCoords(universe=self._root_universe)
    self._root_universe.findRegion(region_id=region_id, univ_coords=localcoords)
    localcoords = localcoords.getHeadNode()
    return localcoords


  def findCoords(self, x=0., y=0., z=0.):

    if self._root_universe is None:
      msg = 'Unable to find coords since the Geometry does not ' \
            'contain the base Universe ID=0'
      raise ValueError(msg)

    point = Point(x=x, y=y, z=z)
    localcoords = UnivCoords(point=point)
    localcoords.setUniverse(self._root_universe)
    cell = self._root_universe.findCell(localcoords=localcoords)
    localcoords = localcoords.getHeadNode()

    return localcoords

  def getNearestIntersection(self, point, direction):

    x, y, z = point._coords
    cell = self.findCoords(x=x, y=y, z=z).getTailNode()._cell
    print self.findCoords(x=x,y=y,z=z).getTailNode()
    print cell
    surfaces = cell._surfaces
    points = []

    for surface in surfaces.keys():

      intersect = surfaces[surface][0].getIntersectionPoints(point, direction)

      if intersect is not None:
        points.extend(intersect)

    if points == []:
      return None

    nearestpoint = points[0]
    nearestdist = point.distanceToPoint(points[0])

    for intersect in points:
      if point.distanceToPoint(intersect) < nearestdist:
        nearestpoint = intersect
        nearestdist = point.distanceToPoint(intersect)

    return nearestpoint
