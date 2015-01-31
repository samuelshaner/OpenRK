__author__ = 'Samuel Shaner'
__email__ = 'shaner@mit.edu'

# Import modules
from math import *
import numpy as np
from checkvalue import *
from cell import *
from clock import *

# A static variable for auto-generated Mesh UIDs
AUTO_MESH_UID = 1


class Mesh(object):

  def __init__(self, mesh_id=None, name='', width=1.0, height=1.0):
        
    # Initialize class attributes
    global AUTO_MESH_UID
    self._uid = AUTO_MESH_UID
    AUTO_MESH_UID += 1
    self._id = None
    self._set_id = False
    self._name = ''

    # Initialize class attributes
    self._cells = None
    self._x_min = None
    self._x_max = None
    self._y_min = None
    self._y_max = None
    self._offset = np.zeros(2)
    self._num_cells = None
    self._boundaries = np.zeros(4)
    self._num_shape_energy_groups = None
    self._num_amp_energy_groups = None
    self._num_delayed_groups = None
    self._flux = {}
    self._temperature = {}
    self._clock = None
    
    # Set Mesh properties
    self.setXMin(-width/2.0)
    self.setXMax(width/2.0)
    self.setYMin(-height/2.0)
    self.setYMax(height/2.0)
    self.setName(name)

    # Set mesh ID
    if not mesh_id is None:
      self.setId(mesh_id)
    else:
      self.setId(self._uid)
      self._set_id = False


  def setClock(self, clock):

    # Initialize clock
    if not isinstance(clock, Clock):
      msg = 'Unable to initialize Mesh clock since clock input is not of type '\
          'Clock: {0}'.format(clock)
    else:
      self._clock = clock
    

  def setName(self, name):

    if not is_string(name):
      msg = 'Unable to set name for Mesh ID={0} with a non-string ' \
            'value {1}'.format(self._id, (name))
      raise ValueError(msg)

    else:
      self._name = name


  def setId(self, mesh_id):

    # Check that the ID is a non-negative integer
    if is_integer(mesh_id):

      if mesh_id >= 0:
        self._id = mesh_id
        self._set_id = True
      else:
        msg = 'Unable to set Mesh ID to {0} since it must be a ' \
              'non-negative integer'.format(mesh_id)
        raise ValueError(msg)

    else:
       msg = 'Unable to set Mesh ID to non-integer {0}'.format(mesh_id)
       raise ValueError(msg)


  def setXMin(self, x_min):

    if not is_float(x_min):
      msg = 'Unable to set x_min to non-integer {0}'\
          .format(x_min)
      raise ValueError(msg)
        
    else:
      self._x_min = x_min


  def setXMax(self, x_max):

    if not is_float(x_max):
      msg = 'Unable to set x_max to non-integer {0}'\
          .format(x_max)
      raise ValueError(msg)
        
    else:
      self._x_max = x_max


  def setYMin(self, y_min):

    if not is_float(y_min):
      msg = 'Unable to set y_min to non-integer {0}'\
          .format(y_min)
      raise ValueError(msg)
        
    else:
      self._y_min = y_min


  def setYMax(self, y_max):

    if not is_float(y_max):
      msg = 'Unable to set y_max to non-integer {0}'\
          .format(y_max)
      raise ValueError(msg)
        
    else:
      self._y_max = y_max


  def getWidth(self):    
    return (self._x_max - self._x_min)


  def getHeight(self):    
    return (self._y_max - self._y_min)


  def setBoundary(self, side, boundary):

    if not is_integer(side):
      msg = 'Unable to set boundary for non-integer side {0}'\
          .format(side)
      raise ValueError(msg)

    elif not is_integer(boundary):
      msg = 'Unable to set boundary for non-integer boundary {0}'\
          .format(boundary)
      raise ValueError(msg)

    if side < 0 or side > 4:
      msg = 'Unable to set boundary for invalid side {0}'\
          .format(side)
      raise ValueError(msg)

    if boundary < 0 or boundary > 2:
      msg = 'Unable to set boundary for invalid boundary {0}'\
          .format(boundary)
      raise ValueError(msg)

    else:
      self._boundaries[side] = boundary


  def getBoundary(self, side):

    if not is_integer(side):
      msg = 'Unable to get boundary for non-integer side {0}'\
          .format(side)
      raise ValueError(msg)

    elif side < 0 or side > 4:
      msg = 'Unable to set boundary for invalid side {0}'\
          .format(side)
      raise ValueError(msg)

    else:
      return self._boundaries[side]


  def getBounds(self):
    
    return [self._x_min, self._x_max, self._y_min, self._y_max]


  def setNumShapeEnergyGroups(self, num_groups):

    if not is_integer(num_groups):
      msg = 'Unable to set num shape energy groups for non-integer {0}'\
          .format(num_groups)
      raise ValueError(msg)

    elif num_groups < 1:
      msg = 'Unable to set num shape energy groups for non-positive {0}'\
          .format(num_groups)
      raise ValueError(msg)
    
    else:
      self._num_shape_energy_groups = num_groups


  def getNumShapeEnergyGroups(self):

    check_set(self._num_shape_energy_groups, 'Mesh get num shape energy groups', 'self._num_shape_energy_groups')
    return self._num_shape_energy_groups


  def setNumAmpEnergyGroups(self, num_groups):

    if not is_integer(num_groups):
      msg = 'Unable to set num amp energy groups for non-integer {0}'\
          .format(num_groups)
      raise ValueError(msg)

    elif num_groups < 1:
      msg = 'Unable to set num amp energy groups for non-positive {0}'\
          .format(num_groups)
      raise ValueError(msg)
    
    else:
      self._num_amp_energy_groups = num_groups


  def getNumAmpEnergyGroups(self):

    check_set(self._num_amp_energy_groups, 'Mesh get num amp energy groups', 'self._num_amp_energy_groups')
    return self._num_amp_energy_groups


  def setNumDelayedGroups(self, num_groups):

    if not is_integer(num_groups):
      msg = 'Unable to set num delayed groups for non-integer {0}'\
          .format(num_groups)
      raise ValueError(msg)

    elif num_groups < 1:
      msg = 'Unable to set num delayed groups for non-positive {0}'\
          .format(num_groups)
      raise ValueError(msg)
    
    else:
      self._num_delayed_groups = num_groups
 

  def __repr__(self):

    string = 'OpenRK Mesh\n'
    string += ' Name \t\t\t = {0} \n'.format(self._name)
    string += ' ID \t\t\t = {0} \n'.format(self._id)
    string += ' X bounds \t\t = [{0}, {1}] \n'.format(self._x_min, self._x_max)
    string += ' Y bounds \t\t = [{0}, {1}] \n'.format(self._y_min, self._y_max)
    string += ' Width \t\t\t = {0} \n'.format(self._x_max - self._x_min)
    string += ' Height \t\t = {0} \n'.format(self._y_max - self._y_min)
    string += ' Offset \t\t = {0}\n'.format(self._offset)
    string += ' Boundaries \t\t = {0}\n'.format(self._boundaries)
    string += ' Num Shape Groups \t = {0}\n'.format(self._num_shape_energy_groups)
    string += ' Num Amp Groups \t = {0}\n'.format(self._num_amp_energy_groups)
    string += ' Num Delayed Groups \t = {0}\n'.format(self._num_delayed_groups)

    return string


class StructuredMesh(Mesh):

  def __init__(self, mesh_id=None, name='', width=1.0, height=1.0, num_x=1, num_y=1):

    # initialize FunctionalMaterial class attributes
    super(StructuredMesh, self).__init__(mesh_id, name, width, height)
    
    # Initialize class attributes
    self._num_x = None
    self._num_y = None
    self._cell_width  = None
    self._cell_height = None
    
    # Set Mesh properties
    self.setNumX(num_x)
    self.setNumY(num_y)
    self._current = {}


  def getCell(self, x, y):

    # Check input values
    check_is_int(x, 'Structured Mesh {0} get cell'.format(self._name), 'x') 
    check_is_int(y, 'Structured Mesh {0} get cell'.format(self._name), 'y') 

    if self._cells is None:
      msg = 'Cannot get cell from Structured Mesh as cells have not been initialized!'
      raise AttributeError(msg)

    else:
      
      return self._cells[y][x]


  def getNeighborCell(self, x, y, side):

    # Check input values
    check_is_int(x, 'Structured Mesh {0} get neighbor cell'.format(self._name), 'x') 
    check_is_int(y, 'Structured Mesh {0} get neighbor cell'.format(self._name), 'y') 
    check_is_int(side, 'Structured Mesh {0} get neighbor cell'.format(self._name), 'side') 

    if self._cells is None:
      msg = 'Cannot get neighbor cell from Structured Mesh as cells have not been initialized!'
      raise AttributeError(msg)

    elif side < 0 or side > 3:
      msg = 'Cannot get neighbor cell from Structured Mesh with invalid side {0}.'.format(side)
      raise AttributeError(msg)
      
    else:

      neighbor_cell = None

      if side == 0:
        if x != 0:
          neighbor_cell = self._cells[y][x-1]
      elif side == 1:
        if y != 0:
          neighbor_cell = self._cells[y-1][x]
      elif side == 2:
        if x != self._num_x - 1:
          neighbor_cell = self._cells[y][x+1]
      elif side == 3:
        if y != self._num_y - 1:
          neighbor_cell = self._cells[y+1][x]

      return neighbor_cell


  def getNumX(self):

    return self._num_x


  def getNumY(self):

    return self._num_y


  def getCellWidth(self):

    if self._cell_width is None:
      msg = 'Cannot get cell width of Structured Mesh as num cells x has not been set!'
      raise AttributeError(msg)
      
    else:
      
      return self._cell_width


  def getCellHeight(self):

    if self._cell_height is None:
      msg = 'Cannot get cell height of Structured Mesh as num cells y has not been set!'
      raise AttributeError(msg)
      
    else:
      
      return self._cell_height


  def initializeCells(self):
                
    self._cells = np.empty(shape=[self._num_y, self._num_x], dtype=object)
   
    # create cell objects
    for y in range(self._num_y):
      for x in range(self._num_x):
        self._cells[y][x] = CmfdCell()


  def setNumX(self, num_x):

    # Check input values
    check_is_int(num_x, 'Structured mesh num x', 'num x')

    if num_x < 1:
      msg = 'Unable to set num cells x to non-positive {0}'\
          .format(num_x)
      raise ValueError(msg)
        
    else:

      self._num_x = num_x
      self._cell_width = (self._x_max - self._x_min) / num_x    


  def setNumY(self, num_y):

    # Check input values
    check_is_int(num_y, 'Structured mesh num y', 'num y')
    
    if num_y < 1:
      msg = 'Unable to set num cells x to non-positive {0}'\
          .format(num_y)
      raise ValueError(msg)
        
    else:
      self._num_y = num_y
      self._cell_height = (self._y_max - self._y_min) / num_y    


  def findCell(self, x, y):

    # Check input values
    check_is_float_or_int(x, 'Structured mesh find cell', 'x')
    check_is_float_or_int(y, 'Structured mesh find cell', 'y')

    i = floor((x - self._x_min) / self._cell_width)
    j = floor((y - self._y_min) / self._cell_height)

    return j*self._num_x + i


  def __repr__(self):

    string = 'OpenRK StructuredMesh\n'
    string += ' Name \t\t\t = {0} \n'.format(self._name)
    string += ' ID \t\t\t = {0} \n'.format(self._id)
    string += ' X bounds \t\t = [{0}, {1}] \n'.format(self._x_min, self._x_max)
    string += ' Y bounds \t\t = [{0}, {1}] \n'.format(self._y_min, self._y_max)
    string += ' Width \t\t\t = {0} \n'.format(self._x_max - self._x_min)
    string += ' Height \t\t = {0} \n'.format(self._y_max - self._y_min)
    string += ' Cell Width \t\t = {0} \n'.format(self._cell_width)
    string += ' Cell Height \t\t = {0} \n'.format(self._cell_height)
    string += ' Cells X \t\t = {0} \n'.format(self._num_x)
    string += ' Cells Y \t\t = {0} \n'.format(self._num_y)
    string += ' Offset \t\t = {0}\n'.format(self._offset)
    string += ' Boundaries \t\t = {0}\n'.format(self._boundaries)
    string += ' Num Shape Groups \t = {0}\n'.format(self._num_shape_energy_groups)
    string += ' Num Amp Groups \t = {0}\n'.format(self._num_amp_energy_groups)
    string += ' Num Delayed Groups \t = {0}\n'.format(self._num_delayed_groups)

    return string


class AmpMesh(StructuredMesh):

  def __init__(self, mesh_id=None, name='', width=1.0, height=1.0, num_x=1, num_y=1):

    # initialize FunctionalMaterial class attributes
    super(AmpMesh, self).__init__(mesh_id, name, width, height, num_x, num_y)

    # Initialize flux and current arrays
    clock = Clock()    
    for position in clock.getPositions():
      self._flux[position] = np.empty
      self._current[position] = np.empty

  def setNumFSRs(self, num_fsrs):

    self._num_fsrs = num_fsrs


  def initializeCells(self):
                
    self._cells = np.empty(shape=[self._num_y, self._num_x], dtype=object)
   
    # create cell objects
    for y in range(self._num_y):
      for x in range(self._num_x):
        self._cells[y][x] = TcmfdCell()


  def initializeSurfaces(self):

    for y in xrange(self._num_y):
      for x in xrange(self._num_x):
         
        cell = self._cells[y][x]
        
        # Surface 0
        if x == 0:
          cell._surfaces[0] = Surface(cell._material._num_energy_groups)
        else:
          cell._surfaces[0] = self._cells[y][x-1]._surfaces[2]
                    
        # Surface 1
        if y == 0:
          cell._surfaces[1] = Surface(cell._material._num_energy_groups)
        else:
          cell._surfaces[1] = self._cells[y-1][x]._surfaces[3]
      
        # Surface 2 
        cell._surfaces[2] = Surface(cell._material._num_energy_groups)

        # Surface 3
        cell._surfaces[3] = Surface(cell._material._num_energy_groups)


  def getFlux(self, name, cell, group, time='CURRENT'):

    check_clock_position(time, 'AmpMesh flux')
    return self._flux[name][cell*self._num_amp_energy_groups+group]


  def initializeFlux(self):

    clock = Clock()    
    for position in clock.getPositions():
      self._flux[position] = np.zeros(self._num_x*self._num_y*self._num_amp_energy_groups)


  def initializeCurrent(self):

    clock = Clock()    
    for position in clock.getPositions():
      self._current[position] = np.zeros(self._num_x*self._num_y*self._num_amp_energy_groups*4)


  def __repr__(self):

    string = 'OpenRK AmpMesh\n'
    string += ' Name \t\t\t = {0} \n'.format(self._name)
    string += ' ID \t\t\t = {0} \n'.format(self._id)
    string += ' X bounds \t\t = [{0}, {1}] \n'.format(self._x_min, self._x_max)
    string += ' Y bounds \t\t = [{0}, {1}] \n'.format(self._y_min, self._y_max)
    string += ' Width \t\t\t = {0} \n'.format(self._x_max - self._x_min)
    string += ' Height \t\t = {0} \n'.format(self._y_max - self._y_min)
    string += ' Cell Width \t\t = {0} \n'.format(self._cell_width)
    string += ' Cell Height \t\t = {0} \n'.format(self._cell_height)
    string += ' Cells X \t\t = {0} \n'.format(self._num_x)
    string += ' Cells Y \t\t = {0} \n'.format(self._num_y)
    string += ' Offset \t\t = {0}\n'.format(self._offset)
    string += ' Boundaries \t\t = {0}\n'.format(self._boundaries)
    string += ' Num Shape Groups \t = {0}\n'.format(self._num_shape_energy_groups)
    string += ' Num Amp Groups \t = {0}\n'.format(self._num_amp_energy_groups)
    string += ' Num Delayed Groups \t = {0}\n'.format(self._num_delayed_groups)

    return string


class UnstructuredMesh(Mesh):

  def __init__(self, mesh_id=None, name='', width=1.0, height=1.0):

    # initialize FunctionalMaterial class attributes
    super(UnstructuredMesh, self).__init__(mesh_id, name, width, height)
    
    # Initialize class attributes
    self._num_cells = None


  def setNumCells(self, num_cells):

    if not is_integer(num_cells):
      msg = 'Unable to set num cells for non-integer {0}'\
          .format(num_cells)
      raise ValueError(msg)

    elif num_cells < 1:
      msg = 'Unable to set num cells for non-positive {0}'\
          .format(num_cells)
      raise ValueError(msg)
    
    else:
      self._num_cells = num_cells


  def initializeCells(self):
                
    self._cells = np.empty(shape=[self._num_cells], dtype=object)
   
    # create cell objects
    for i in range(self._num_cells):
      self._cells[i] = MOCCell()


  def __repr__(self):

    string = 'OpenRK UnstructuredMesh\n'
    string += ' Name \t\t\t = {0} \n'.format(self._name)
    string += ' ID \t\t\t = {0} \n'.format(self._id)
    string += ' X bounds \t\t = [{0}, {1}] \n'.format(self._x_min, self._x_max)
    string += ' Y bounds \t\t = [{0}, {1}] \n'.format(self._y_min, self._y_max)
    string += ' Width \t\t\t = {0} \n'.format(self._x_max - self._x_min)
    string += ' Height \t\t = {0} \n'.format(self._y_max - self._y_min)
    string += ' Offset \t\t = {0}\n'.format(self._offset)
    string += ' Num cells \t\t = {0}\n'.format(self._num_cells)
    string += ' Boundaries \t\t = {0}\n'.format(self._boundaries)
    string += ' Num Shape Groups \t = {0}\n'.format(self._num_shape_energy_groups)
    string += ' Num Amp Groups \t = {0}\n'.format(self._num_amp_energy_groups)
    string += ' Num Delayed Groups \t = {0}\n'.format(self._num_delayed_groups)

    return string


class UnstructuredShapeMesh(UnstructuredMesh):

  def __init__(self, mesh_id=None, name='', width=1.0, height=1.0):

    # initialize FunctionalMaterial class attributes
    super(MOCMesh, self).__init__(mesh_id, name, width, height)

    # Initialize flux and current arrays
    clock = Clock()    
    for position in clock.getPositions():
      self._flux[position] = np.empty


  def initializeFlux(self):

    # Initialize flux and current arrays
    clock = Clock()    
    for position in clock.getPositions():
      self._flux[position] = np.zeros(self._num_cells*self._num_shape_energy_groups)


  def __repr__(self):

    string = 'OpenRK MOCMesh\n'
    string += ' Name \t\t\t = {0} \n'.format(self._name)
    string += ' ID \t\t\t = {0} \n'.format(self._id)
    string += ' X bounds \t\t = [{0}, {1}] \n'.format(self._x_min, self._x_max)
    string += ' Y bounds \t\t = [{0}, {1}] \n'.format(self._y_min, self._y_max)
    string += ' Width \t\t\t = {0} \n'.format(self._x_max - self._x_min)
    string += ' Height \t\t = {0} \n'.format(self._y_max - self._y_min)
    string += ' Offset \t\t = {0}\n'.format(self._offset)
    string += ' Num cells \t\t = {0}\n'.format(self._num_cells)
    string += ' Boundaries \t\t = {0}\n'.format(self._boundaries)
    string += ' Num Shape Groups \t = {0}\n'.format(self._num_shape_energy_groups)
    string += ' Num Amp Groups \t = {0}\n'.format(self._num_amp_energy_groups)
    string += ' Num Delayed Groups \t = {0}\n'.format(self._num_delayed_groups)

    return string
