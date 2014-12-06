__author__ = 'Samuel Shaner'
__email__ = 'shaner@mit.edu'


from openrk.checkvalue import *
import numpy as np

# A list of all IDs for all Materials created
MATERIAL_IDS = list()

# A static variable for auto-generated Material IDs
AUTO_MATERIAL_ID = 10000

def reset_auto_material_id():
  global AUTO_MATERIAL_ID, MATERIAL_IDS
  AUTO_MATERIAL_ID = 10000
  MATERIAL_IDS = list()


class Material(object):

  def __init__(self, material_id=None, name='', num_energy_groups=1):

    # Initialize class attributes
    self._id = None
    self._name = ''
    self._sigma_a = None
    self._sigma_f = None
    self._nu_sigma_f = None
    self._sigma_s = None
    self._dif_coef = None

    # Set the Material class attributes
    self.setId(material_id)
    self.setName(name)
    self.setNumEnergyGroups(num_energy_groups)

  def __deepcopy__(self, memo):

    existing = memo.get(id(self))

    # If this is the first time we have tried to copy this object, create a copy
    if existing is None:

      clone = type(self).__new__(type(self))
      clone._id = self._id
      clone._name = self._name

      memo[id(self)] = clone

      return clone

    # If this object has been copied before, return the first copy made
    else:
      return existing


  def __gt__(self, other):
    return (id(self) > id(other))


  def __ge__(self, other):
    return (id(self) >= id(other))


  def __lt__(self, other):
    return (id(self) < id(other))


  def __le__(self, other):
    return (id(self) <= id(other))


  def setId(self, material_id=None):

    global MATERIAL_IDS

    if material_id is None:
      global AUTO_MATERIAL_ID
      self._id = AUTO_MATERIAL_ID
      MATERIAL_IDS.append(AUTO_MATERIAL_ID)
      AUTO_MATERIAL_ID += 1

    # Check that the ID is an integer and wasn't already used
    elif is_integer(material_id):

      # If the Material already has an ID, remove it from global list
      if self._id is not None:
        MATERIAL_IDS.remove(self._id)

      if material_id in MATERIAL_IDS:
        msg = 'Unable to set Material ID to {0} since a Material ' \
              'with this ID was already initialized'.format(material_id)
        raise ValueError(msg)

      if material_id < 0:
        msg = 'Unable to set Material ID to {0} since it must be a ' \
              'non-negative integer'.format(material_id)
        raise ValueError(msg)

      else:
        self._id = material_id
        MATERIAL_IDS.append(material_id)

    else:
      msg = 'Unable to set Material ID to non-integer {0}'.format(material_id)
      raise ValueError(msg)


  def setName(self, name):

    if not is_string(name):
      msg = 'Unable to set name for Material ID={0} with a non-string ' \
            'value {1}'.format(self._id, (name))
      raise ValueError(msg)

    else:
      self._name = name


  def __repr__(self):

    string = 'Material\n'
    string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self._id)
    string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self._name)

    return string


  def setSigmaA(self, sigma_a):

    # check if sigma_a is a list
    if not is_list(sigma_a):
      msg = 'Unable to set absorption xs for Material ID={0} with a '\
            'non-list value {1}'.format(self._id, (sigma_a))
      raise ValueError(msg)

    # check if sigma_a is of length num_groups
    elif len(sigma_a) != self._num_energy_groups:
      msg = 'Unable to set absorption xs for Material ID={0} with {1} groups '\
          'as num energy groups is set to {2}'\
          .format(self._id, len(sigma_a), self._num_energy_groups)
      raise ValueError(msg)

    else:
      if isinstance(sigma_a, list):
        self._sigma_a = np.asarray(sigma_a)
      else:
        self._sigma_a = np.copy(sigma_a)


  def setSigmaF(self, sigma_f):

    # check if sigma_f is a list
    if not is_list(sigma_f):
      msg = 'Unable to set fission xs for Material ID={0} with a '\
            'non-list value {1}'.format(self._id, (sigma_f))
      raise ValueError(msg)

    # check if sigma_f is of length num_groups
    elif len(sigma_f) != self._num_energy_groups:
      msg = 'Unable to set fission xs for Material ID={0} with {1} groups '\
          'as num energy groups is set to {2}'\
          .format(self._id, len(sigma_f), self._num_energy_groups)
      raise ValueError(msg)

    else:
      if isinstance(sigma_f, list):
        self._sigma_f = np.asarray(sigma_f)
      else:
        self._sigma_f = np.copy(sigma_f)


  def setNuSigmaF(self, nu_sigma_f):

    # check if nu_sigma_f is a list
    if not is_list(nu_sigma_f):
      msg = 'Unable to set nu fission xs for Material ID={0} with a '\
            'non-list value {1}'.format(self._id, (nu_sigma_f))
      raise ValueError(msg)

    # check if nu_sigma_f is of length num_groups
    elif len(nu_sigma_f) != self._num_energy_groups:
      msg = 'Unable to set nu fission xs for Material ID={0} with {1} groups '\
          'as num energy groups is set to {2}'\
          .format(self._id, len(nu_sigma_f), self._num_energy_groups)
      raise ValueError(msg)

    else:
      if isinstance(nu_sigma_f, list):
        self._nu_sigma_f = np.asarray(nu_sigma_f)
      else:
        self._nu_sigma_f = np.copy(nu_sigma_f)


  def setSigmaS(self, sigma_s):

    # check if sigma_s is a list
    if not is_list(sigma_s):
      msg = 'Unable to set scattering xs for Material ID={0} with a '\
            'non-list value {1}'.format(self._id, (sigma_s))
      raise ValueError(msg)

    # check if sigma_s is of length num_energy_groups*num_energy_groups
    elif np.shape(sigma_s) != (self._num_energy_groups, self._num_energy_groups):
      msg = 'Unable to set scattering xs for Material ID={0} with {1} groups '\
          'as num energy groups is set to {2}x{3}'\
          .format(self._id, self._num_energy_groups, np.shape(sigma_s)[0], np.shape(sigma_s)[1])
      raise ValueError(msg)

    else:
      if isinstance(sigma_s, list):
        self._sigma_s = np.asarray(sigma_s)
      else:
        self._sigma_s = np.copy(sigma_s)


  def setDifCoef(self, dif_coef):

    # check if dif_coef is a list
    if not is_list(dif_coef):
      msg = 'Unable to set the diffusion coefficients for Material ID={0} '\
          'with a non-list value {1}'.format(self._id, (dif_coef))
      raise ValueError(msg)

    # check if dif_coef is of length num_energy_groups
    elif len(dif_coef) != self._num_energy_groups:
      msg = 'Unable to set the diffusion coefficients for Material ID={0} with'\
          ' {1} groups as num energy groups is set to {2}'\
          .format(self._id, len(dif_coef), self._num_energy_groups)
      raise ValueError(msg)

    else:
      if isinstance(dif_coef, list):
        self._dif_coef = np.asarray(dif_coef)
      else:
        self._dif_coef = np.copy(dif_coef)


  def setSigmaAByGroup(self, sigma_a, group):

    # check if sigma_a is a float
    if not is_float(sigma_a):
      msg = 'Unable to set absorption xs for Material ID={0} with a '\
            'non-float value {1}'.format(self._id, (sigma_a))
      raise ValueError(msg)

    # check if group is valid
    elif group < 1 or group > self._num_energy_groups:
      msg = 'Unable to set absorption xs for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    else:
      self._sigma_a[group] = sigma_a


  def setSigmaFByGroup(self, sigma_f, group):

    # check if sigma_f is a float
    if not is_float(sigma_f):
      msg = 'Unable to set fission xs for Material ID={0} with a '\
            'non-float value {1}'.format(self._id, (sigma_f))
      raise ValueError(msg)

    # check if group is valid
    elif group < 1 or group > self._num_energy_groups:
      msg = 'Unable to set fission xs for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    else:
      self._sigma_f[group] = sigma_f


  def setNuSigmaFByGroup(self, nu_sigma_f, group):

    # check if nu_sigma_f is a list
    if not is_float(nu_sigma_f):
      msg = 'Unable to set nu fission xs for Material ID={0} with a '\
            'non-float value {1}'.format(self._id, (nu_sigma_f))
      raise ValueError(msg)

    # check if group is valid
    elif group < 1 or group > self._num_energy_groups:
      msg = 'Unable to set fission xs for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    else:
      self._nu_sigma_f[group] = nu_sigma_f


  def setSigmaSByGroup(self, sigma_s, group_from, group_to):

    # check if sigma_s is a list
    if not is_float(sigma_s):
      msg = 'Unable to set scattering xs for Material ID={0} with a '\
            'non-float value {1}'.format(self._id, (sigma_s))
      raise ValueError(msg)

    # check if group is valid
    elif group_from < 1 or group_from > self._num_energy_groups or \
    group_to < 1 or group_to > self._num_energy_groups:
      msg = 'Unable to set scattering xs for Material ID={0} for group {1} '\
          'to {2} as num energy groups is set to {3}'\
          .format(self._id, group_to, group_from, self._num_energy_groups)
      raise ValueError(msg)

    else:
      self._sigma_s[group_from][group_to] = sigma_s


  def setDifCoefByGroup(self, dif_coef, group):

    # check if nu_sigma_f is a list
    if not is_float(dif_coef):
      msg = 'Unable to set dif coef for Material ID={0} with a '\
            'non-float value {1}'.format(self._id, (dif_coef))
      raise ValueError(msg)

    # check if group is valid
    elif group < 1 or group > self._num_energy_groups:
      msg = 'Unable to set dif coef for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    else:
      self._dif_coef[group] = dif_coef


  def getSigmaAByGroup(self, group):

    # check if group is valid
    if group < 1 or group > self._num_energy_groups:
      msg = 'Unable to get absorption xs for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    else:
      return self._sigma_a[group]


  def getSigmaFByGroup(self, group):

    # check if group is valid
    if group < 1 or group > self._num_energy_groups:
      msg = 'Unable to get fission xs for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    else:
      return self._sigma_f[group]


  def getNuSigmaFByGroup(self, group):

    # check if group is valid
    if group < 1 or group > self._num_energy_groups:
      msg = 'Unable to get fission xs for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    else:
      return self._nu_sigma_f[group]


  def getSigmaSByGroup(self, group_from, group_to):

    # check if group is valid
    if group_from < 1 or group_from > self._num_energy_groups or \
    group_to < 1 or group_to > self._num_energy_groups:
      msg = 'Unable to get scattering xs for Material ID={0} for group {1} '\
          'to {2} as num energy groups is set to {3}'\
          .format(self._id, group_to, group_from, self._num_energy_groups)
      raise ValueError(msg)

    else:
      return self._sigma_s[group_from][group_to]


  def getDifCoefByGroup(self, group):

    # check if group is valid
    if group < 1 or group > self._num_energy_groups:
      msg = 'Unable to set dif coef for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    else:
      return self._dif_coef[group]


  def setNumEnergyGroups(self, num_energy_groups):

    if is_integer(num_energy_groups):
   
      if num_energy_groups < 1:
        msg = 'Unable to set num energy groups to {0} since it must be an ' \
              ' integer > 0'.format(num_energy_groups)
        raise ValueError(msg)

      else:
        self._num_energy_groups = num_energy_groups
        self._sigma_a = np.zeros(num_energy_groups)
        self._sigma_f = np.zeros(num_energy_groups)
        self._nu_sigma_f = np.zeros(num_energy_groups)
        self._sigma_s = np.zeros((num_energy_groups, num_energy_groups))
        self._dif_coef = np.zeros(num_energy_groups)

    else:
      msg = 'Unable to set num energy groups to non-integer {0}'\
          .format(num_energy_groups)
      raise ValueError(msg)


class FunctionalMaterial(Material):

  def __init__(self, material_id=None, name='', num_energy_groups=1,
               num_time_steps=2):

    # initialize FunctionalMaterial class attributes
    super(FunctionalMaterial, self).__init__(material_id, name, num_energy_groups)

    # Initialize class attributes
    self._num_time_steps = None

    # Set the FunctionalMaterial class attributes
    self.setNumTimeSteps(num_time_steps)
    self.setNumEnergyGroups(num_energy_groups)


  def setSigmaA(self, sigma_a):

    # check if sigma_a is a list
    if not is_list(sigma_a):
      msg = 'Unable to set absorption xs for Material ID={0} with a '\
            'non-list value {1}'.format(self._id, (sigma_a))
      raise ValueError(msg)

    # check if sigma_a is of length num_time_steps x num_energy_groups
    elif np.shape(sigma_a) != (self.num_time_steps, self._num_energy_groups):
      msg = 'Unable to set absorption xs for Material ID={0} with {1} groups '\
          'and {2} time steps. Num groups is {3} and num time steps is {4}.'\
          .format(self._id, np.shape(sigma_a)[0], np.shape(sigma_a)[1], \
                  self._num_energy_groups, self._num_time_steps)
      raise ValueError(msg)

    else:
      if isinstance(sigma_a, list):
        self._sigma_a = np.asarray(sigma_a)
      else:
        self._sigma_a = np.copy(sigma_a)


  def setSigmaF(self, sigma_f):

    # check if sigma_f is a list
    if not is_list(sigma_f):
      msg = 'Unable to set fission xs for Material ID={0} with a '\
            'non-list value {1}'.format(self._id, (sigma_f))
      raise ValueError(msg)

    # check if sigma_f is of length num_time_steps x num_energy_groups
    elif np.shape(sigma_f) != (self.num_time_steps, self._num_energy_groups):
      msg = 'Unable to set fission xs for Material ID={0} with {1} groups '\
          'and {2} time steps. Num groups is {3} and num time steps is {4}.'\
          .format(self._id, np.shape(sigma_f)[0], np.shape(sigma_f)[1], \
                  self._num_energy_groups, self._num_time_steps)
      raise ValueError(msg)

    else:
      if isinstance(sigma_f, list):
        self._sigma_f = np.asarray(sigma_f)
      else:
        self._sigma_f = np.copy(sigma_f)


  def setNuSigmaF(self, nu_sigma_f):

    # check if nu_sigma_f is a list
    if not is_list(sigma_f):
      msg = 'Unable to set nu fission xs for Material ID={0} with a '\
            'non-list value {1}'.format(self._id, (nu_sigma_f))
      raise ValueError(msg)

    # check if nu_sigma_f is of length num_time_steps x num_energy_groups
    elif np.shape(nu_sigma_f) != (self.num_time_steps, self._num_energy_groups):
      msg = 'Unable to set nu fission xs for Material ID={0} with {1} groups '\
          'and {2} time steps. Num groups is {3} and num time steps is {4}.'\
          .format(self._id, np.shape(nu_sigma_f)[0], np.shape(nu_sigma_f)[1], \
                  self._num_energy_groups, self._num_time_steps)
      raise ValueError(msg)

    else:
      if isinstance(nu_sigma_f, list):
        self._nu_sigma_f = np.asarray(nu_sigma_f)
      else:
        self._nu_sigma_f = np.copy(nu_sigma_f)


  def setSigmaS(self, sigma_s):

    # check if sigma_s is a list
    if not is_list(sigma_s):
      msg = 'Unable to set scattering xs for Material ID={0} with a '\
            'non-list value {1}'.format(self._id, (sigma_s))
      raise ValueError(msg)

    # check if sigma_s is of length num_time_steps x num_energy_groups x num_energy_groups
    elif np.shape(nu_sigma_f) != (self.num_time_steps, self._num_energy_groups, self._num_energy_groups):
      msg = 'Unable to set scattering xs for Material ID={0} with {1}/{2} groups '\
          'and {3} time steps. Num groups is {3} and num time steps is {4}.'\
          .format(self._id, np.shape(sigma_s)[0], np.shape(sigma_s)[1], np.shape(sigma_s)[2],
                  self._num_energy_groups, self._num_time_steps)
      raise ValueError(msg)

    else:
      if isinstance(sigma_s, list):
        self._sigma_s = np.asarray(sigma_s)
      else:
        self._sigma_s = np.copy(sigma_s)


  def setDifCoef(self, dif_coef):

    # check if dif_coef is a list
    if not is_list(dif_coef):
      msg = 'Unable to set the diffusion coefficients for Material ID={0} '\
          'with a non-list value {1}'.format(self._id, (dif_coef))
      raise ValueError(msg)

    # check if dif_coef is of length num_time_steps x num_energy_groups
    elif np.shape(dif_coef) != (self.num_time_steps, self._num_energy_groups):
      msg = 'Unable to set dif coef for Material ID={0} with {1} groups '\
          'and {2} time steps. Num groups is {3} and num time steps is {4}.'\
          .format(self._id, np.shape(dif_coef)[0], np.shape(dif_coef)[1], \
                  self._num_energy_groups, self._num_time_steps)
      raise ValueError(msg)

    else:
      if isinstance(dif_coef, list):
        self._dif_coef = np.asarray(dif_coef)
      else:
        self._dif_coef = np.copy(dif_coef)


  def setSigmaAByGroup(self, sigma_a, group, time_step):

    # check if sigma_a is a float
    if not is_float(sigma_a):
      msg = 'Unable to set absorption xs for Material ID={0} with a '\
            'non-float value {1}'.format(self._id, (sigma_a))
      raise ValueError(msg)

    # check if group is valid
    elif group < 1 or group > self._num_energy_groups:
      msg = 'Unable to set absorption xs for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    # check if time_step is valid
    elif time_step < 1 or time_step > self._num_time_steps:
      msg = 'Unable to set absorption xs for Material ID={0} for time step {1} '\
          'as num time steps is set to {2}'\
          .format(self._id, group, self._num_time_steps)
      raise ValueError(msg)

    else:
      self._sigma_a[time_step][group] = sigma_a


  def setSigmaFByGroup(self, sigma_f, group, time_step):

    # check if sigma_f is a float
    if not is_float(sigma_f):
      msg = 'Unable to set fission xs for Material ID={0} with a '\
            'non-float value {1}'.format(self._id, (sigma_f))
      raise ValueError(msg)

    # check if group is valid
    elif group < 1 or group > self._num_energy_groups:
      msg = 'Unable to set fission xs for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    # check if time_step is valid
    elif time_step < 1 or time_step > self._num_time_steps:
      msg = 'Unable to set fission xs for Material ID={0} for time step {1} '\
          'as num time steps is set to {2}'\
          .format(self._id, group, self._num_time_steps)
      raise ValueError(msg)

    else:
      self._sigma_f[time_step][group] = sigma_f


  def setNuSigmaFByGroup(self, nu_sigma_f, group, time_step):

    # check if nu_sigma_f is a list
    if not is_float(nu_sigma_f):
      msg = 'Unable to set nu fission xs for Material ID={0} with a '\
            'non-float value {1}'.format(self._id, (nu_sigma_f))
      raise ValueError(msg)

    # check if group is valid
    elif group < 1 or group > self._num_energy_groups:
      msg = 'Unable to set fission xs for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    # check if time_step is valid
    elif time_step < 1 or time_step > self._num_time_steps:
      msg = 'Unable to set nu fission xs for Material ID={0} for time step {1} '\
          'as num time steps is set to {2}'\
          .format(self._id, group, self._num_time_steps)
      raise ValueError(msg)

    else:
      self._nu_sigma_f[time_step][group] = nu_sigma_f


  def setSigmaSByGroup(self, sigma_s, group_from, group_to, time_step):

    # check if sigma_s is a list
    if not is_float(sigma_s):
      msg = 'Unable to set scattering xs for Material ID={0} with a '\
            'non-float value {1}'.format(self._id, (sigma_s))
      raise ValueError(msg)

    # check if group is valid
    elif group_from < 1 or group_from > self._num_energy_groups or \
    group_to < 1 or group_to > self._num_energy_groups:
      msg = 'Unable to set scattering xs for Material ID={0} for group {1} '\
          'to {2} as num energy groups is set to {3}'\
          .format(self._id, group_to, group_from, self._num_energy_groups)
      raise ValueError(msg)

    # check if time_step is valid
    elif time_step < 1 or time_step > self._num_time_steps:
      msg = 'Unable to set absorption xs for Material ID={0} for time step {1} '\
          'as num time steps is set to {2}'\
          .format(self._id, group, self._num_time_steps)
      raise ValueError(msg)

    else:
      self._sigma_s[time_step][group_from][group_to] = sigma_s


  def setDifCoefByGroup(self, dif_coef, group, time_step):

    # check if nu_sigma_f is a list
    if not is_float(dif_coef):
      msg = 'Unable to set dif coef for Material ID={0} with a '\
            'non-float value {1}'.format(self._id, (dif_coef))
      raise ValueError(msg)

    # check if group is valid
    elif group < 1 or group > self._num_energy_groups:
      msg = 'Unable to set dif coef for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    # check if time_step is valid
    elif time_step < 1 or time_step > self._num_time_steps:
      msg = 'Unable to set dif coef xs for Material ID={0} for time step {1} '\
          'as num time steps is set to {2}'\
          .format(self._id, group, self._num_time_steps)
      raise ValueError(msg)

    else:
      self._dif_coef[time_step][group] = dif_coef


  def getSigmaAByGroup(self, group, time_step):

    # check if group is valid
    if group < 1 or group > self._num_energy_groups:
      msg = 'Unable to get absorption xs for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    # check if time_step is valid
    elif time_step < 1 or time_step > self._num_time_steps:
      msg = 'Unable to get absorption xs for Material ID={0} for time step {1} '\
          'as num time steps is set to {2}'\
          .format(self._id, group, self._num_time_steps)
      raise ValueError(msg)

    else:
      return self._sigma_a[time_step][group]


  def getSigmaFByGroup(self, group, time_step):

    # check if group is valid
    if group < 1 or group > self._num_energy_groups:
      msg = 'Unable to get fission xs for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    # check if time_step is valid
    elif time_step < 1 or time_step > self._num_time_steps:
      msg = 'Unable to set fission xs for Material ID={0} for time step {1} '\
          'as num time steps is set to {2}'\
          .format(self._id, group, self._num_time_steps)
      raise ValueError(msg)

    else:
      return self._sigma_f[time_step][group]


  def getNuSigmaFByGroup(self, group, time_step):

    # check if group is valid
    if group < 1 or group > self._num_energy_groups:
      msg = 'Unable to get fission xs for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    # check if time_step is valid
    elif time_step < 1 or time_step > self._num_time_steps:
      msg = 'Unable to set nu fission xs for Material ID={0} for time step {1} '\
          'as num time steps is set to {2}'\
          .format(self._id, group, self._num_time_steps)
      raise ValueError(msg)

    else:
      return self._nu_sigma_f[time_step][group]


  def getSigmaSByGroup(self, group_from, group_to, time_step):

    # check if group is valid
    if group_from < 1 or group_from > self._num_energy_groups or \
    group_to < 1 or group_to > self._num_energy_groups:
      msg = 'Unable to get scattering xs for Material ID={0} for group {1} '\
          'to {2} as num energy groups is set to {3}'\
          .format(self._id, group_to, group_from, self._num_energy_groups)
      raise ValueError(msg)

    # check if time_step is valid
    elif time_step < 1 or time_step > self._num_time_steps:
      msg = 'Unable to set scattering xs for Material ID={0} for time step {1} '\
          'as num time steps is set to {2}'\
          .format(self._id, group, self._num_time_steps)
      raise ValueError(msg)

    else:
      return self._sigma_s[time_step][group_from][group_to]


  def getDifCoefByGroup(self, group, time_step):

    # check if group is valid
    if group < 1 or group > self._num_energy_groups:
      msg = 'Unable to set dif coef for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    # check if time_step is valid
    elif time_step < 1 or time_step > self._num_time_steps:
      msg = 'Unable to set dif coef xs for Material ID={0} for time step {1} '\
          'as num time steps is set to {2}'\
          .format(self._id, group, self._num_time_steps)
      raise ValueError(msg)

    else:
      return self._dif_coef[time_step][group]


  def setNumTimeSteps(self, num_time_steps):

    if is_integer(num_time_steps):
   
      if num_time_steps < 2:
        msg = 'Unable to set num time steps to {0} since it must be an ' \
              ' integer > 1'.format(num_time_steps)
        raise ValueError(msg)

      else:
        self._num_time_steps = num_time_steps
        self._time_steps = np.zeros(num_time_steps)

    else:
      msg = 'Unable to set num time steps to non-integer {0}'\
          .format(num_time_steps)
      raise ValueError(msg)


  def setNumEnergyGroups(self, num_energy_groups):

    if is_integer(num_energy_groups):
   
      if num_energy_groups < 1:
        msg = 'Unable to set num energy groups to {0} since it must be an ' \
              ' integer > 0'.format(num_energy_groups)
        raise ValueError(msg)

      else:
        self._num_energy_groups = num_energy_groups
        self._sigma_a = np.zeros((num_time_steps, num_energy_groups))
        self._sigma_f = np.zeros((num_time_steps, num_energy_groups))
        self._nu_sigma_f = np.zeros((num_time_steps, num_energy_groups))
        self._sigma_s = np.zeros((num_time_steps, num_energy_groups, \
                                    num_energy_groups))
        self._dif_coef = np.zeros((num_time_steps, num_energy_groups))

    else:
      msg = 'Unable to set num energy groups to non-integer {0}'\
          .format(num_energy_groups)
      raise ValueError(msg)


