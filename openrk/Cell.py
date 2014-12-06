
from math import *
import numpy as np
from Surface import *
from Material import *

class Cell(object):
    
  def __init__(self):

    # Initialize class attributes
    self._surfaces = np.empty(4, dtype=object)
    self._neighbor_cells = np.empty(4, dtype=object)    
    self._volume = None        
    self._material = None
    self._flux = None
    self._temperature = None
    self._clock = None
    self._id = None

    
  def setMaterial(self, material):

    if not isinstance(material, Material):
      msg = 'Unable to set material for Cell ID={0} with a '\
            'non-material value {1}'.format(self._id, material)

    else:
      self._material = material
      self._flux = np.zeros(material._num_energy_groups)


