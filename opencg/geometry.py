__author__ = 'Will Boyd'
__email__ = 'wboyd@mit.edu'


import opencg
from opencg.universe import *
from opencg.point import *
from opencg.ray import *
import copy

# Small displacement for moving a point across a surface in ray tracing
TINY_BIT = 1e-10

def reset_auto_ids():
  opencg.reset_auto_material_id()
  opencg.reset_auto_surface_id()
  opencg.reset_auto_cell_id()
  opencg.reset_auto_universe_id()


class Geometry(object):

  def __init__(self):

    # Initialize Geometry class attributes
    self._root_universe = None
    self._num_regions = 0
    self._volume = np.float64(0.)

    # A NumPy array of volumes for each region, indexed by Region ID
    self._region_volumes = None

    self._all_cells = None

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

    # Flag indicating whether or not the neighbors have been built
    self._built_neighbors = False


  def __deepcopy__(self, memo):

    existing = memo.get(id(self))

    # If this is the first time we have tried to copy this object, create a copy
    if existing is None:

      clone = type(self).__new__(type(self))
      clone._root_universe = copy.deepcopy(self._root_universe, memo)
      clone._num_regions = self._num_regions
      clone._volume = self._volume
      clone._region_volumes = copy.deepcopy(self._region_volumes, memo)
      clone._num_neighbors = self._num_neighbors
      clone._neighbor_ids = copy.deepcopy(self._neighbor_ids, memo)
      clone._num_unique_neighbors = copy.deepcopy(self._num_unique_neighbors, memo)
      clone._unique_neighbor_ids = copy.deepcopy(self._unique_neighbor_ids, memo)
      clone._regions_to_neighbors = copy.deepcopy(self._regions_to_neighbors, memo)
      clone._regions_to_unique_neighbors = copy.deepcopy(self._regions_to_unique_neighbors, memo)
      clone._built_neighbors = self._built_neighbors
      clone._all_cells = copy.deepcopy(self._all_cells, memo)

      memo[id(self)] = clone

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
    bounds = [self._root_universe.getMinX(),
              self._root_universe.getMaxX(),
              self._root_universe.getMinY(),
              self._root_universe.getMaxY(),
              self._root_universe.getMinZ(),
              self._root_universe.getMaxZ()]

    return bounds


  def getAllCells(self):

    if self._all_cells is None:
      self._all_cells = self._root_universe.getAllCells()

    return self._all_cells


  def getAllUniverses(self):
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


  def updateBoundingBoxes(self):

    if self._root_universe is None:
      msg = 'Unable to get the bounds since the Geometry does not ' \
            'contain the base Universe ID=0'
      raise ValueError(msg)

    cells = self.getAllCells()

    for cell_id, cell in cells.items():
      cell.findBoundingBox()


  def initializeCellOffsets(self):

    if self._root_universe is None:
      msg = 'Unable to initialize cell offsets since the Geometry ' \
            'does not contain the base Universe ID=0'
      raise ValueError(msg)

    self._root_universe.initializeCellOffsets()
    self._num_regions = self._root_universe._num_regions


  def clearNeighbors(self):
    self._num_neighbors = 0
    self._neighbor_ids = dict()
    self._num_unique_neighbors = 0
    self._unique_neighbor_ids = dict()
    self._regions_to_neighbors = dict()
    self._regions_to_unique_neighbors = dict()
    self._built_neighbors = False


  def buildNeighbors(self):

    if self._root_universe is None:
      msg = 'Unable to build neighbor Cells/Universes since the ' \
            'root Universe for the Geometry has not yet been set'
      raise ValueError(msg)

    # If the neighbors have already been built, just return
    if self._built_neighbors:
      return

    self._root_universe.buildNeighbors()

    # Initialize offsets maps
    self.initializeCellOffsets()

    # Initialize dictionaries mapping neighbor hashes to consecutive integers
    # Keys    - hashes of the tuples of (unique) neighbors
    # Values  - monotonically consecutive non-negative integers
    # Reinitialize each time this routine is called
    self._num_neighbors = 0
    self._neighbor_ids = dict()
    self._num_unique_neighbors = 0
    self._unique_neighbor_ids = dict()

    # Initialize dictionaries mapping regions to neighbors
    # Keys    - region IDs
    # Values  - hashes of the tuples of (unique) neighbors
    self._regions_to_neighbors = dict()
    self._regions_to_unique_neighbors = dict()

    # Set a flag indicating that the neighbors have been built
    self._built_neighbors = True


  def countNeighbors(self):

    for region in range(self._num_regions):
      self.getNeighborsHash(region)
      self.getUniqueNeighborsHash(region)


  def getRegionId(self, x=0., y=0., z=0.):

    coords = self.findCoords(x=x, y=y, z=z)
    region_id = 0

    # FIXME!!!
    # If we did not find the coords, return NaN as an error code
    if coords._cell is None:
      return np.nan

    while coords is not None:

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
    distances = []

    # initialize linked list for point
    next = self.findCoords(x=x, y=y, z=z)

    # Loop through linked list and retrieve intersection distances
    while next is not None:
      if isinstance(next, UnivCoords):
        cell = next._cell
        if cell is not None:
          surfaces = cell._surfaces
          for surface in surfaces.values():
            dist = surface[0].minSurfaceDist(next._point, direction)
            if dist is not None:
              distances.append(dist)

      elif isinstance(next, LatCoords):
        lat = next._lattice
        x_lat, y_lat, z_lat = next._point._coords - lat._offset
        x_lat = x_lat + 0.5*lat._dimension[0]*lat._width[0] - next._lat_x*lat._width[0] - lat._halfwidth[0]
        y_lat = y_lat + 0.5*lat._dimension[1]*lat._width[1] - next._lat_y*lat._width[1] - lat._halfwidth[1]
        z_lat = z_lat + 0.5*lat._dimension[2]*lat._width[2] - next._lat_z*lat._width[2] - lat._halfwidth[2]
        lat_point = Point(x=x_lat, y=y_lat, z=z_lat)
        dist = lat.minSurfaceDist(lat_point, direction)
        if dist is not None:
          distances.append(dist)

      next = next._next

    if distances == []:
      return None


    # find smallest distance
    min_dist = min(distances)


    # sets coordinates for nearest intersection
    nearestpoint = Point()
    poldir = direction.toPolar()
    sines = np.sin(poldir)
    cosines = np.cos(poldir)
    nearestpoint.setX(point._coords[0] + min_dist*sines[2]*cosines[1])
    nearestpoint.setY(point._coords[1] + min_dist*sines[2]*sines[1])
    nearestpoint.setZ(point._coords[2] + min_dist*cosines[2])

    return nearestpoint

  def getShortestSegment(self, point, direction):

    intersect = self.getNearestIntersection(point, direction)
    if intersect is None:
      return None

    segment = Segment(geometry=self, start=point, end=intersect)
    return segment

  def traceRays(self, rays):

    start = Point()
    for ray in rays:

      # sets starting point of ray
      start.setCoords(ray._point._coords)
      direction = ray._direction

      # traces ray until edge of the geometry is reached
      intersect = self.getNearestIntersection(start, direction)

      while intersect is not None:
        segment = Segment(geometry=self, start=start, end=intersect)
        ray.addSegment(segment)

        # adjusts next segment in ray to start at found intersection
        start.setCoords(intersect._coords + TINY_BIT*direction._comps)
        intersect = self.getNearestIntersection(start,direction)

    return rays

  def getNeighbors(self, region_id):
    coords = self.findRegion(region_id)
    return coords.getNeighbors()

  def getUniqueNeighbors(self, region_id):
    coords = self.findRegion(region_id)
    return coords.getUniqueNeighbors()


  def getNeighborsHash(self, region_id, first_level=0):

    # Memoize unique neighbors hash for this region ID
    if region_id not in self._regions_to_neighbors:

      # Compute the neighbor hash for this region ID
      coords = self.findRegion(region_id)
      neighbors = coords.getNeighborsHash(first_level=first_level)

      # Add the neighbor hash to the neighbor maps
      if not neighbors in self._neighbor_ids.keys():
        self._neighbor_ids[neighbors] = self._num_neighbors
        self._num_neighbors += 1

      # Store unique hash to the region-to-neighbor hash maps
      self._regions_to_neighbors[region_id] = self._neighbor_ids[neighbors]

    return self._regions_to_neighbors[region_id]


  def getUniqueNeighborsHash(self, region_id, first_level=0):

    # Memoize unique neighbors hash for this region ID
    if region_id not in self._regions_to_unique_neighbors:

      # Compute the unique neighbor hash for this region ID
      coords = self.findRegion(region_id)
      unique_neighbors = coords.getUniqueNeighborsHash(first_level=first_level)

      # Add the unique neighbor hash to the unique neighbor maps
      if not unique_neighbors in self._unique_neighbor_ids.keys():
        self._unique_neighbor_ids[unique_neighbors] = self._num_unique_neighbors
        self._num_unique_neighbors += 1

      # Store unique neighbor hash to the region-to-unique neighbor hash map
      self._regions_to_unique_neighbors[region_id] = \
        self._unique_neighbor_ids[unique_neighbors]

    return self._regions_to_unique_neighbors[region_id]

  def toString(self):
    string = self._root_universe.toString()
    return string
