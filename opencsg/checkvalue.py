__author__ = 'Will Boyd'
__email__ = 'wboyd@mit.edu'

import numpy as np


def is_integer(val):
  return isinstance(val, (int, np.int32, np.int64))


def is_float(val):
  return isinstance(val, (float, np.float32, np.float64))
  

def is_string(val):
  return isinstance(val, (str, np.str))


def dim_list(val):
  dim = 0
  while isinstance(val, list):
    dim += 1
    val = val[0]

  return dim
