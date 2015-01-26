
from math import *
import numpy as np
from checkvalue import *

# A static variable for auto-generated Surface UIDs
AUTO_SURFACE_UID = 1

class Surface(object):
    
  def __init__(self, num_groups):

    # Initialize class attributes
    global AUTO_SURFACE_UID
    self._id = AUTO_SURFACE_UID
    AUTO_SURFACE_UID += 1
  
    # Check input values
    check_is_int(num_groups, 'Surface num energy groups', 'num energy groups')

    # Initialize class attributes
    self._currents = np.zeros(num_groups)
    self._dif_coef_linear = np.zeros(num_groups)
    self._dif_coef_nonlinear = np.zeros(num_groups)
    self._num_energy_groups = num_groups

  
  def setCurrentByGroup(self, current, group):

    # Check input values
    check_is_float_or_int(current, 'Surface current', 'current')
    check_is_int(group, 'Surface group', 'current')

    # check if group is valid
    if group < 0 or group > self._num_energy_groups-1:
      msg = 'Unable to set current for Surface ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    else:
      
      self._currents[group] = current


  def setDifCoefLinearByGroup(self, dif_coef, group):

    # Check input values
    check_is_float_or_int(dif_coef, 'Surface coef linear', 'dif coef linear')
    check_is_int(group, 'Surface group', 'dif coef linear')

    # check if group is valid
    if group < 0 or group > self._num_energy_groups-1:
      msg = 'Unable to set dif coef linear for Surface ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    else:
      
      self._dif_coef_linear[group] = dif_coef


  def setDifCoefNonlinearByGroup(self, dif_coef, group):

    # Check input values
    check_is_float_or_int(dif_coef, 'Surface coef nonlinear', 'dif coef linear')
    check_is_int(group, 'Surface group', 'dif coef nonlinear')

    # check if group is valid
    if group < 0 or group > self._num_energy_groups-1:
      msg = 'Unable to set dif coef nonlinear for Surface ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    else:
      
      self._dif_coef_nonlinear[group] = dif_coef


  def getDifCoefLinearByGroup(self, group):

    # Check input values
    check_is_int(group, 'Surface group', 'dif coef linear')

    # check if group is valid
    if group < 0 or group > self._num_energy_groups-1:
      msg = 'Unable to get linear dif coef for Surface ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    else:
      
      return self._dif_coef_linear[group]


  def getDifCoefNonlinearByGroup(self, dif_coef):

    # Check input values
    check_is_int(group, 'Surface group', 'dif coef nonlinear')

    # check if group is valid
    if group < 0 or group > self._num_energy_groups-1:
      msg = 'Unable to get nonlinear dif coef for Surface ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    else:
      
      self._dif_coef_nonlinear[group] = dif_coef

