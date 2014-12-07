
from math import *
import numpy as np

class Surface(object):
    
  def __init__(self, num_groups):
  
    # Initialize class attributes
    self._current = np.zeros(num_groups)
    self._dif_hat = np.zeros(num_groups)
    self._dif_tilde = np.zeros(num_groups)
  
