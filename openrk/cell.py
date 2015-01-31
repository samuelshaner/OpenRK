__author__ = 'Samuel Shaner'
__email__ = 'shaner@mit.edu'

# Import modules
from math import *
import numpy as np
from checkvalue import *
from clock import *
from material import *
from surface import *

# A static variable for auto-generated Material UIDs
AUTO_CELL_UID = 1

class Cell(object):
    
  def __init__(self, clock=None):

    # Initialize class attributes
    global AUTO_CELL_UID
    self._id = AUTO_CELL_UID
    AUTO_CELL_UID += 1

    # Initialize class attributes
    self._volume = None        
    self._material = None
    self._clock = None
    self._id = None

    # Initialize clock
    if clock is None:
      self._clock = Clock()
    elif not isinstance(clock, Clock):
      msg = 'Unable to initialize Cell since clock input is not of type '\
          'Clock: {0}'.format(clock)
    else:
      self._clock = clock

    # Initialize temperatures
    self._temperature = {}
    for position in self._clock._positions:
      self._temperature[position] =  0.0

    
  def setMaterial(self, material):

    if not isinstance(material, Material):
      msg = 'Unable to set material for Cell ID={0} with a '\
            'non-material value {1}'.format(self._id, material)
      
      raise ValueError(msg)

    else:
      self._material = material


  def getMaterial(self):

    if self._material is None:
      msg = 'Unable to get material for Cell ID={0} since the '\
            'material has not been set'.format(self._id)      
      raise ValueError(msg)

    else:
      return self._material


  def setTemperature(self, temperature, clock_position):

    # A static variable for auto-generated Material UIDs
    global CLOCK_POSITIONS

    # Check input values
    check_is_float_or_int(temperature, 'Cell ID={0} temperature'.format(self._id), 'temperature')

    if clock_position not in CLOCK_POSITIONS:
      msg = 'Clock position {0} not in valid clock positions: {1}'.format(clock_position, CLOCK_POSITIONS)
      raise ValueError(msg)

    else:

      self._temperature[clock_position] = temperature


  def setClock(self, clock):

    self._clock = clock


class MOCCell(Cell):
    
  def __init__(self, clock=None):

    # initialize FunctionalMaterial class attributes
    super(MOCCell, self).__init__(clock)


class CmfdCell(Cell):
    
  def __init__(self, clock=None):

    # initialize FunctionalMaterial class attributes
    super(CmfdCell, self).__init__(clock)

    # Initialize class attributes
    self._surfaces = np.empty(4, dtype=object)
    

  def initializeSurfaces(self):

    # Check input values
    check_set(self._material, 'initialize Cell ID={0} surfaces'.format(self._id), 'self._material')

    for i in range(4):
      self._surfaces[i] = Surface(self._material.getNumEnergyGroups())


  def getSurface(self, side):

    # Check input values
    check_is_int(side, 'Cell ID={0} get surface'.format(self._id), 'side')

    if side < 0 or side > 3:
      msg = 'Cannot get surface for cell ID={0} for invalid side: {1}'\
        .format(self._id, side)
      raise ValueError(msg)

    elif self._surfaces[side] is None:
      msg = 'Cannot get surface for cell ID={0} since the surfaces have not been initialized!'\
        .format(self._id, side)
      raise ValueError(msg)

    else:
        
      return self._surfaces[side]
    

class TcmfdCell(CmfdCell):
    
  def __init__(self, clock=None):

    # initialize FunctionalMaterial class attributes
    super(TcmfdCell, self).__init__(clock)

    # Initialize class attributes
    self._fsrs = None
    self._num_fsrs = None
    self._moc_groups_map = None

  def setNumFSRs(self, num_fsrs):

    # Check input values
    check_is_int(num_fsrs, 'Cell ID={0} number of fsrs'.format(self._id), 'num fsrs')

    # check if group is valid
    if num_fsrs < 1:
      msg = 'Unable to set number of fsrs for Cell ID={0} for non-positive '\
          'number of fsrs {1}.'\
          .format(self._id, num_fsrs)
      raise ValueError(msg)

    else:

      self._num_fsrs = num_fsrs
      self._fsrs = np.zeros(num_fsrs, dtype=int)


