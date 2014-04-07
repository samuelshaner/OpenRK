__author__ = 'Will Boyd'
__email__ = 'wboyd@mit.edu'

import numpy as np


def is_integer(val):
  return isinstance(val, (int, np.int8, np.int16, np.int32, \
                          np.int64, np.int128))


def is_float(val):
  return isinstance(val, (float, np.float16, np.int32, np.int64, \
                          np.float96, np.float128, np.float256))


def is_string(val):
  return isinstance(val, (str, np.str))

