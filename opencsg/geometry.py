__author__ = 'Will Boyd'
__email__ = 'wboyd@mit.edu'


from opencsg.universe import *
from opencsg.point import *

tiny_bit = 1e-5

class Geometry(object):

  def __init__(self):

    # Initialize Geometry class attributes
    self._root_universe = None
    self._num_regions = 0
    self._volume = np.float64(0.)

    # A NumPy array of volumes for each region, indexed by Region ID
    self._region_volumes = None

    # Dictionaries mapping neighbor hashes to consecutive integers
    # Keys    - hashes of the tuples of (unique) neighbors
    # Values  - monotonically consecutive non-negative integers
    self._num_neighbors = 0
    self._neighbor_ids = dict()
    self._num_unique_neighbors = 0
    self._unique_neighbor_ids = dict()

    # Map regions to neighbors
    # Keys    - region IDs
    # Values  - hashes of the tuples of (unique) neighbors
    self._regions_to_neighbors = dict()
    self._regions_to_unique_neighbors = dict()


  def __deepcopy__(self, memo):

    existing = memo.get(self)

    # If this is the first time we have tried to copy this object, create a copy
    if existing is None:

      clone = type(self).__new__(type(self))
      clone._root_universe = copy.deepcopy(self._root_universe)
      clone._num_regions = self._num_regions
      clone._volume = self._volume
      clone._region_volumes = copy.deepcopy(self._region_volumes)
      clone._num_neighbors = self._num_neighbors
      clone._neighbor_ids = copy.deepcopy(self._neighbor_ids)
      clone._num_unique_neighbors = copy.deepcopy(self._num_unique_neighbors)
      clone._unique_neighbor_ids = copy.deepcopy(self._unique_neighbor_ids)
      clone._regions_to_neighbors = copy.deepcopy(self._regions_to_neighbors)
      clone._regions_to_unique_neighbors = copy.deepcopy(self._regions_to_unique_neighbors)

      return clone

    # If this object has been copied before, return the first copy made
    else:
      return existing


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
    materials = dict()

    for cell_id, cell in material_cells.items():
      materials[cell._fill._id] = cell._fill

    return materials


  def getAllMaterialCells(self):

    all_cells = self.getAllCells()
    material_cells = dict()

    for cell_id, cell in all_cells.items():
      if cell._type == 'material':
        material_cells[cell._id] = cell

    return material_cells


  def getAllMaterialUniverses(self):

    all_universes = self.getAllUniverses()
    material_universes = dict()

    for universe_id, universe in all_universes.items():

      # Do not consider Lattices since they are not at the Material level
      if isinstance(universe, Lattice):
        continue

      cells = universe._cells

      for cell_id, cell in cells.items():
        if cell._type == 'material':
          material_universes[universe._id] = universe

    return material_universes


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


  def buildNeighbors(self):

    if self._root_universe is None:
      msg = 'Unable to build neighbor Cells/Universes since the ' \
            'root Universe for the Geometry has not yet been set'
      raise ValueError(msg)

    # If the neighbors have already been built, just return
    if self._num_neighbors > 0:
      return

    self._root_universe.buildNeighbors()

    # Initialize offsets maps
    self.initializeCellOffsets()

    # Build dictionaries mapping neighbor hashes to consecutive integers
    # Keys    - hashes of the tuples of (unique) neighbors
    # Values  - monotonically consecutive non-negative integers
    # Reinitialize each time this routine is called
    self._num_neighbors = 0
    self._neighbor_ids = dict()
    self._num_unique_neighbors = 0
    self._unique_neighbor_ids = dict()

    # Map regions to neighbors
    # Keys    - region IDs
    # Values  - hashes of the tuples of (unique) neighbors
    self._regions_to_neighbors = dict()
    self._regions_to_unique_neighbors = dict()

    for region in range(self._num_regions):

      # Build lists of neighbor Cells/Universes
      neighbors = self.getNeighborsHash(region)
      unique_neighbors = self.getUniqueNeighborsHash(region)

      # Store the hashes to the region-to-neighbor hash maps
      self._regions_to_neighbors[region] = neighbors
      self._regions_to_unique_neighbors[region] = unique_neighbors

      # Add the neighbor hash to the neighbor maps
      if not neighbors in self._neighbor_ids.keys():
        self._neighbor_ids[neighbors] = self._num_neighbors
        self._num_neighbors += 1

      # Add the unique neighbor hash to the unique neighbor maps
      if not unique_neighbors in self._unique_neighbor_ids.keys():
        self._unique_neighbor_ids[unique_neighbors] = self._num_unique_neighbors
        self._num_unique_neighbors += 1


  def getRegionId(self, x=0., y=0., z=0.):

    coords = self.findCoords(x=x, y=y, z=z)
    region_id = 0

    # FIXME!!!
    # If we did not find the coords, return NaN as an error code
    if coords._cell is None:
      return np.nan

    while not coords is None:

      # The coords is a UnivCoords object
      if isinstance(coords, UnivCoords):
        universe = coords._universe
        cell = coords._cell
        region_id += universe.getCellOffset(cell)

      # The coords is a LatCoords object
      else:
        lattice = coords._lattice
        lat_x, lat_y, lat_z = coords._lat_x, coords._lat_y, coords._lat_z
        region_id += lattice.getCellOffset(lat_x, lat_y, lat_z)

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
    if cell is None:
      return None
    surfaces = cell._surfaces
    points = []

    for surface in surfaces.keys():
      intersect = surfaces[surface][0].getIntersectionPoints(point, direction)
      points.extend(intersect)

    while None in points:
      points.remove(None)

    if points == []:
      return None

    nearestpoint = points[0]
    nearestdist = point.distanceToPoint(points[0])

    for intersect in points:
      if point.distanceToPoint(intersect) < nearestdist:
        nearestpoint = intersect
        nearestdist = point.distanceToPoint(intersect)

    return nearestpoint

  def getShortestSegment(self, point, direction):

    intersect = self.getNearestIntersection(point, direction)
    if intersect is None:
      return None

    segment = Segment(self, start=point, end=intersect)
    return segment

  def traceSampleRays(self, num_rays=1000):

    rays = dict()
    bounds = self.getBounds()

    for ray in xrange(num_rays):
      edge = np.random.randint(4)
      if edge == 0:
        x = bounds[edge] + tiny_bit
        y = np.random.uniform(bounds[2], bounds[3])
        z = np.random.uniform(-1e12, 1e12)
      elif edge == 1:
        x = bounds[edge] - tiny_bit
        y = np.random.uniform(bounds[2], bounds[3])
        z = np.random.uniform(-1e12, 1e12)
      elif edge == 2:
        x = np.random.uniform(bounds[0], bounds[1])
        y = bounds[edge] + tiny_bit
        z = np.random.uniform(-1e12, 1e12)
      else:
        x = np.random.uniform(bounds[0], bounds[1])
        y = bounds[edge] - tiny_bit
        z = np.random.uniform(-1e12, 1e12)

      u, v, w = np.random.rand(3)-0.5
      point = Point(x=x, y=y, z=z)
      direction = Direction(u=u, v=v, w=w)
      rays[point] = direction

    segments = []
    while rays != {}:
      points = rays.keys()
      for ray in points:
        intersect = self.getNearestIntersection(ray, rays[ray])
        if intersect is None:
          del rays[ray]
        else:
          segment = Segment(self, start=ray, end=intersect)
          segments.append(segment)
          intersect.setCoords(intersect._coords + tiny_bit*rays[ray]._comps)
          rays[intersect] = rays[ray]
          del rays[ray]
    return segments

  def getNeighbors(self, region_id):
    coords = self.findRegion(region_id)
    return coords.getNeighbors()

  def getUniqueNeighbors(self, region_id):
    coords = self.findRegion(region_id)
    return coords.getUniqueNeighbors()

  def getNeighborsHash(self, region_id):
    coords = self.findRegion(region_id)
    return coords.getNeighborsHash()

  def getUniqueNeighborsHash(self, region_id):
    coords = self.findRegion(region_id)
    return coords.getUniqueNeighborsHash()

