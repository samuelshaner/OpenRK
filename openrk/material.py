__author__ = 'Samuel Shaner'
__email__ = 'shaner@mit.edu'

from checkvalue import *
import numpy as np
from clock import *

# A static variable for auto-generated Material UIDs
AUTO_MATERIAL_UID = 1

class Material(object):

  def __init__(self, material_id=None, name=''):

    # Initialize class attributes
    global AUTO_MATERIAL_UID
    self._uid = AUTO_MATERIAL_UID
    AUTO_MATERIAL_UID += 1
    self._id = None
    self._set_id = False
    self._name = ''

    # Initialize cross section attributes    
    self._sigma_t = None
    self._sigma_a = None
    self._sigma_f = None
    self._nu_sigma_f = None
    self._sigma_s = None
    self._dif_coef = None
    self._chi = None
    self._num_energy_groups = None
 
    # Set the Material class attributes
    self.setName(name)

    if not material_id is None:
      self.setId(material_id)
    else:
      self.setId(self._uid)
      self._set_id = False


  def setName(self, name):

    if not is_string(name):
      msg = 'Unable to set name for Material ID={0} with a non-string ' \
            'value {1}'.format(self._id, (name))
      raise ValueError(msg)

    else:
      self._name = name


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


  def setSigmaT(self, sigma_t):

    # check if sigma_t is a list
    if not is_list(sigma_t):
      msg = 'Unable to set total xs for Material ID={0} with a '\
            'non-list value {1}'.format(self._id, (sigma_t))
      raise ValueError(msg)

    # check if sigma_t is of length num_groups
    elif len(sigma_t) != self._num_energy_groups:
      msg = 'Unable to set total xs for Material ID={0} with {1} groups '\
          'as num energy groups is set to {2}'\
          .format(self._id, len(sigma_t), self._num_energy_groups)
      raise ValueError(msg)

    else:
      if isinstance(sigma_t, list):
        self._sigma_t = np.asarray(sigma_t)
      else:
        self._sigma_t = np.copy(sigma_t)


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
    elif len(sigma_s) != self._num_energy_groups*self._num_energy_groups:
      msg = 'Unable to set scattering xs for Material ID={0} with {1} groups '\
          'as num energy groups is set to {2}'\
          .format(self._id, sqrt(len(_sigma_s)), self._num_energy_groups)
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


  def setChi(self, chi):

    # check if chi is a list
    if not is_list(chi):
      msg = 'Unable to set the chi for Material ID={0} '\
          'with a non-list value {1}'.format(self._id, (chi))
      raise ValueError(msg)

    # check if chi is of length num_energy_groups
    elif len(chi) != self._num_energy_groups:
      msg = 'Unable to set the chi for Material ID={0} with'\
          ' {1} groups as num energy groups is set to {2}'\
          .format(self._id, len(chi), self._num_energy_groups)
      raise ValueError(msg)

    else:
      if isinstance(chi, list):
        self._chi = np.asarray(chi)
      else:
        self._chi = np.copy(chi)


  def setSigmaAByGroup(self, sigma_a, group):

    # check if sigma_a is a float
    if not is_float(sigma_a):
      msg = 'Unable to set absorption xs for Material ID={0} with a '\
            'non-float value {1}'.format(self._id, (sigma_a))
      raise ValueError(msg)

    # check if group is valid
    elif group < 0 or group > self._num_energy_groups:
      msg = 'Unable to set absorption xs for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    else:
      self._sigma_a[group] = sigma_a


  def setSigmaTByGroup(self, sigma_t, group):

    # check if sigma_t is a float
    if not is_float(sigma_t):
      msg = 'Unable to set total xs for Material ID={0} with a '\
            'non-float value {1}'.format(self._id, (sigma_t))
      raise ValueError(msg)

    # check if group is valid
    elif group < 0 or group > self._num_energy_groups-1:
      msg = 'Unable to set total xs for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    else:
      self._sigma_t[group] = sigma_t


  def setSigmaFByGroup(self, sigma_f, group):

    # check if sigma_f is a float
    if not is_float(sigma_f):
      msg = 'Unable to set fission xs for Material ID={0} with a '\
            'non-float value {1}'.format(self._id, (sigma_f))
      raise ValueError(msg)

    # check if group is valid
    elif group < 0 or group > self._num_energy_groups-1:
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
    elif group < 0 or group > self._num_energy_groups-1:
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
    elif group_from < 0 or group_from > self._num_energy_groups-1 or \
    group_to < 0 or group_to > self._num_energy_groups-1:
      msg = 'Unable to set scattering xs for Material ID={0} for group {1} '\
          'to {2} as num energy groups is set to {3}'\
          .format(self._id, group_to, group_from, self._num_energy_groups)
      raise ValueError(msg)

    else:
      self._sigma_s[group_from*self._num_energy_groups+group_to] = sigma_s


  def setDifCoefByGroup(self, dif_coef, group):

    # check if dif_coef is a list
    if not is_float(dif_coef):
      msg = 'Unable to set dif coef for Material ID={0} with a '\
            'non-float value {1}'.format(self._id, (dif_coef))
      raise ValueError(msg)

    # check if group is valid
    elif group < 0 or group > self._num_energy_groups-1:
      msg = 'Unable to set dif coef for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    else:
      self._dif_coef[group] = dif_coef


  def setChiByGroup(self, chi, group):

    # check if chi is a list
    if not is_float(dif_coef):
      msg = 'Unable to set chi for Material ID={0} with a '\
            'non-float value {1}'.format(self._id, (chi))
      raise ValueError(msg)

    # check if group is valid
    elif group < 0 or group > self._num_energy_groups-1:
      msg = 'Unable to set chi for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    else:
      self._chi[group] = chi


  def getSigmaAByGroup(self, group):

    # check if group is valid
    if group < 0 or group > self._num_energy_groups-1:
      msg = 'Unable to get absorption xs for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    else:
      return self._sigma_a[group]


  def getSigmaTByGroup(self, group):

    # check if group is valid
    if group < 0 or group > self._num_energy_groups-1:
      msg = 'Unable to get total xs for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    else:
      return self._sigma_t[group]


  def getSigmaFByGroup(self, group):

    # check if group is valid
    if group < 0 or group > self._num_energy_groups-1:
      msg = 'Unable to get fission xs for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    else:
      return self._sigma_f[group]


  def getNuSigmaFByGroup(self, group):

    # check if group is valid
    if group < 0 or group > self._num_energy_groups-1:
      msg = 'Unable to get fission xs for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    else:
      return self._nu_sigma_f[group]


  def getSigmaSByGroup(self, group_from, group_to):

    # check if group is valid
    if group_from < 0 or group_from > self._num_energy_groups-1 or \
    group_to < 0 or group_to > self._num_energy_groups-1:
      msg = 'Unable to get scattering xs for Material ID={0} for group {1} '\
          'to {2} as num energy groups is set to {3}'\
          .format(self._id, group_to, group_from, self._num_energy_groups)
      raise ValueError(msg)

    else:
      return self._sigma_s[group_from*self._num_energy_groups+group_to]


  def getDifCoefByGroup(self, group):

    # check if group is valid
    if group < 0 or group > self._num_energy_groups-1:
      msg = 'Unable to set dif coef for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    else:
      return self._dif_coef[group]


  def getChiByGroup(self, group):

    # check if group is valid
    if group < 0 or group > self._num_energy_groups-1:
      msg = 'Unable to set chi for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    else:
      return self._chi[group]


  def setNumEnergyGroups(self, num_energy_groups):

    if is_integer(num_energy_groups):
   
      if num_energy_groups < 1:
        msg = 'Unable to set num energy groups to {0} since it must be an ' \
              ' integer > 0'.format(num_energy_groups)
        raise ValueError(msg)

      else:
        self._num_energy_groups = num_energy_groups
        self._sigma_t = np.zeros(num_energy_groups)
        self._sigma_a = np.zeros(num_energy_groups)
        self._sigma_f = np.zeros(num_energy_groups)
        self._nu_sigma_f = np.zeros(num_energy_groups)
        self._sigma_s = np.zeros(num_energy_groups*num_energy_groups)
        self._dif_coef = np.zeros(num_energy_groups)
        self._chi = np.zeros(num_energy_groups)

    else:
      msg = 'Unable to set num energy groups to non-integer {0}'\
          .format(num_energy_groups)
      raise ValueError(msg)

  def setId(self, material_id):

    # Check that the ID is a non-negative integer
    if is_integer(material_id):

      if material_id >= 0:
        self._id = material_id
        self._set_id = True
      else:
        msg = 'Unable to set Material ID to {0} since it must be a ' \
              'non-negative integer'.format(material_id)
        raise ValueError(msg)

    else:
       msg = 'Unable to set Material ID to non-integer {0}'.format(material_id)
       raise ValueError(msg)


  def __repr__(self):

    string = 'OpenRK Material\n'
    string += ' Name \t\t= {0} \n'.format(self._name)
    string += ' ID \t\t= {0} \n'.format(self._id)
    string += ' Num Groups \t= {0} \n'.format(self._num_energy_groups)
    string += ' Sigma T \t= {0} \n'.format(self._sigma_t)
    string += ' Sigma A \t= {0} \n'.format(self._sigma_a)
    string += ' Sigma F \t= {0} \n'.format(self._sigma_f)
    string += ' Nu Sigma F \t= {0} \n'.format(self._nu_sigma_f)
    string += ' Dif Coef \t= {0} \n'.format(self._dif_coef)
    string += ' Chi \t\t= {0} \n'.format(self._chi)
    string += ' Sigma s \t= {0} \n'.format(self._sigma_s)

    return string


class FunctionalMaterial(Material):

  def __init__(self, material_id=None, name=''):

    super(FunctionalMaterial, self).__init__(material_id, name)

    # Initialize class attributes
    self._num_time_steps = None
    self._doppler_coefficients = None
    self._decay_constants = None
    self._delayed_fractions = None
    self._precursor_conc = None
    self._velocity = None
    self._energy_per_fission = None

  def setDopplerCoefficients(self, doppler_coefs):

    # Check if necessary variables have been set
    check_set(self._num_energy_groups, 'Doppler Coefs', 'Num Energy Groups')
    check_list_of_floats_or_ints(doppler_coefs, 'Doppler Coefs')

    # check if doppler_coefs is a list
    if not is_list(doppler_coefs):
      msg = 'Unable to set absorption xs for Material ID={0} with a '\
            'non-list value {1}'.format(self._id, (doppler_coefs))
      raise ValueError(msg)

    else:
      if isinstance(doppler_coefs, list):
        self._doppler_coefficients = np.asarray(doppler_coefs)
      else:
        self._doppler_coefficients = np.copy(doppler_coefs)


  def setSigmaA(self, sigma_a):

    # Check if necessary variables have been set
    check_set(self._num_energy_groups, 'Sigma A', 'Num Energy Groups')
    check_set(self._num_time_steps, 'Sigma A', 'Num Time Steps')
    for i in sigma_a:
      check_list_of_floats_or_ints(i, 'Material ID={0} sigma_a'.format(self._id))

    # check if sigma_a is of length num_time_steps x num_energy_groups
    if np.shape(sigma_a) != (self._num_time_steps, self._num_energy_groups):
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


  def setSigmaT(self, sigma_t):

    # Check if necessary variables have been set
    check_set(self._num_energy_groups, 'Sigma T', 'Num Energy Groups')
    check_set(self._num_time_steps, 'Sigma T', 'Num Time Steps')
    for i in sigma_t:
      check_list_of_floats_or_ints(i, 'Material ID={0} sigma_t'.format(self._id))

    # check if sigma_t is of length num_time_steps x num_energy_groups
    if np.shape(sigma_t) != (self._num_time_steps, self._num_energy_groups):
      msg = 'Unable to set total xs for Material ID={0} with {1} groups '\
          'and {2} time steps. Num groups is {3} and num time steps is {4}.'\
          .format(self._id, np.shape(sigma_t)[0], np.shape(sigma_t)[1], \
                  self._num_energy_groups, self._num_time_steps)
      raise ValueError(msg)

    else:
      if isinstance(sigma_t, list):
        self._sigma_t = np.asarray(sigma_t)
      else:
        self._sigma_t = np.copy(sigma_t)


  def setSigmaF(self, sigma_f):

    # Check if necessary variables have been set
    check_set(self._num_energy_groups, 'Sigma F', 'Num Energy Groups')
    check_set(self._num_time_steps, 'Sigma F', 'Num Time Steps')
    for i in sigma_f:
      check_list_of_floats_or_ints(i, 'Material ID={0} sigma_f'.format(self._id))

    # check if sigma_f is of length num_time_steps x num_energy_groups
    if np.shape(sigma_f) != (self._num_time_steps, self._num_energy_groups):
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

    # Check if necessary variables have been set
    check_set(self._num_energy_groups, 'Nu Sigma F', 'Num Energy Groups')
    check_set(self._num_time_steps, 'Nu Sigma F', 'Num Time Steps')
    for i in nu_sigma_f:
      check_list_of_floats_or_ints(i, 'Material ID={0} nu_sigma_f'.format(self._id))

    # check if nu_sigma_f is of length num_time_steps x num_energy_groups
    if np.shape(nu_sigma_f) != (self._num_time_steps, self._num_energy_groups):
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

    # Check if necessary variables have been set
    check_set(self._num_energy_groups, 'Sigma S', 'Num Energy Groups')
    check_set(self._num_time_steps, 'Sigma S', 'Num Time Steps')
    for i in sigma_s:
      check_list_of_floats_or_ints(i, 'Material ID={0} sigma_s'.format(self._id))


    # check if sigma_s is of length num_time_steps x num_energy_groups x num_energy_groups
    if np.shape(sigma_s) != (self._num_time_steps, self._num_energy_groups*self._num_energy_groups):
      msg = 'Unable to set scattering xs for Material ID={0} with {1} groups '\
          'and {3} time steps. Num groups is {3} and num time steps is {4}.'\
          .format(self._id, np.shape(sigma_s)[0], sqrt(np.shape(sigma_s)[1]),
                  self._num_energy_groups, self._num_time_steps)
      raise ValueError(msg)

    else:
      if isinstance(sigma_s, list):
        self._sigma_s = np.asarray(sigma_s)
      else:
        self._sigma_s = np.copy(sigma_s)


  def setDifCoef(self, dif_coef):

    # Check if necessary variables have been set
    check_set(self._num_energy_groups, 'Dif Coefs', 'Num Energy Groups')
    check_set(self._num_time_steps, 'Dif Coefs', 'Num Time Steps')
    for i in dif_coef:
      check_list_of_floats_or_ints(i, 'Material ID={0} dif_coef'.format(self._id))

    # check if dif_coef is of length num_time_steps x num_energy_groups
    if np.shape(dif_coef) != (self._num_time_steps, self._num_energy_groups):
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


  def setChi(self, chi):

    # Check if necessary variables have been set
    check_set(self._num_energy_groups, 'Chi', 'Num Energy Groups')
    check_set(self._num_time_steps, 'Chi', 'Num Time Steps')
    for i in chi:
      check_list_of_floats_or_ints(i, 'Material ID={0} chi'.format(self._id))

    # check if chi is of length num_time_steps x num_energy_groups
    if np.shape(chi) != (self._num_time_steps, self._num_energy_groups):
      msg = 'Unable to set chi for Material ID={0} with {1} groups '\
          'and {2} time steps. Num groups is {3} and num time steps is {4}.'\
          .format(self._id, np.shape(chi)[0], np.shape(chi)[1], \
                  self._num_energy_groups, self._num_time_steps)
      raise ValueError(msg)

    else:
      if isinstance(chi, list):
        self._chi = np.asarray(chi)
      else:
        self._chi = np.copy(chi)


  def setEnergyPerFission(self, energy_per_fission):

    # Check if necessary variables have been set
    check_set(self._num_energy_groups, 'Energy Per Fission', 'Num Energy Groups')
    check_list_of_floats_or_ints(energy_per_fission, 'Material ID={0} energy_per_fission'.format(self._id))

    if isinstance(energy_per_fission, list):
      self._energy_per_fission = np.asarray(energy_per_fission)
    else:
      self._energy_per_fission = np.copy(energy_per_fission)


  def getSigmaAByGroup(self, group, time, temp):

    # Check if necessary variables have been set
    check_set(self._num_energy_groups, 'Material ID={0} sigma_a by group'.format(self._id), 'Num Energy Groups')
    check_set(self._doppler_coefficients, 'Material ID={0} sigma_a by group'.format(self._id), 'Doppler Coefficients')
    check_is_int(group, 'Material ID={0} sigma_a by group'.format(self._id), 'group')
    check_is_float_or_int(time, 'Material ID={0} sigma_a by group'.format(self._id), 'time')
    check_is_float_or_int(group, 'Material ID={0} sigma_a by group'.format(self._id), 'temp')

    # check if group is valid
    if group < 0 or group > self._num_energy_groups-1:
      msg = 'Unable to get absorption xs for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    # check if time_step is valid
    elif time < self._time_steps[0] or time > self._time_steps[-1]:
      msg = 'Unable to get absorption xs for Material ID={0} for group {1} and'\
            ' time {2} outside time window [{3}, {4}]'\
          .format(self._id, group, time, self._time_steps[0], self._time_steps[-1])
      raise ValueError(msg)

    else:

      # find time step window
      (i_left, i_right) = self.getBoundingTimeSteps(time)

      sigma_a_left = self._sigma_a[i_left][group] * self._doppler_coefficients[group] * \
                     (sqrt(temp) - sqrt(300.0))
      sigma_a_right = self._sigma_a[i_right][group] * self._doppler_coefficients[group] * \
                      (sqrt(temp) - sqrt(300.0))

      # compute sigma_a
      sigma_a = sigma_a_left + (time - self._time_steps[i_left]) * \
                (sigma_a_right - sigma_a_left) / \
                (self._time_steps[i_right] - self._time_steps[i_left])

      return sigma_a


  def getSigmaTByGroup(self, group, time, temp):

    # Check if necessary variables have been set
    check_set(self._num_energy_groups, 'Material ID={0} sigma_t by group'.format(self._id), 'Num Energy Groups')
    check_is_int(group, 'Material ID={0} sigma_t by group'.format(self._id), 'group')
    check_is_float_or_int(time, 'Material ID={0} sigma_t by group'.format(self._id), 'time')
    check_is_float_or_int(group, 'Material ID={0} sigma_t by group'.format(self._id), 'temp')

    # check if group is valid
    if group < 0 or group > self._num_energy_groups-1:
      msg = 'Unable to get total xs for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    # check if time_step is valid
    elif time < self._time_steps[0] or time > self._time_steps[-1]:
      msg = 'Unable to get total xs for Material ID={0} for group {1} and '\
            'time {2} outside time window [{3}, {4}]'\
          .format(self._id, group, time, self._time_steps[0], self._time_steps[-1])
      raise ValueError(msg)

    else:

      sigma_t = getSigmaAByGroup(group, time, temp)

      for g in xrange(self._num_energy_groups):
        sigma_t += getSigmaSByGroup(group, g, time, temp)

      return sigma_t


  def getSigmaFByGroup(self, group, time, temp):

    # Check if necessary variables have been set
    check_set(self._num_energy_groups, 'Material ID={0} sigma_f by group'.format(self._id), 'Num Energy Groups')
    check_is_int(group, 'Material ID={0} sigma_f by group'.format(self._id), 'group')
    check_is_float_or_int(time, 'Material ID={0} sigma_f by group'.format(self._id), 'time')
    check_is_float_or_int(group, 'Material ID={0} sigma_f by group'.format(self._id), 'temp')

    # check if group is valid
    if group < 0 or group > self._num_energy_groups-1:
      msg = 'Unable to get fission xs for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    # check if time_step is valid
    elif time < self._time_steps[0] or time > self._time_steps[-1]:
      msg = 'Unable to get fission xs for Material ID={0} for group {1} and '\
            'time {2} outside time window [{3}, {4}]'\
          .format(self._id, group, time, self._time_steps[0], self._time_steps[-1])
      raise ValueError(msg)

    else:

      # find time step window
      (i_left, i_right) = self.getBoundingTimeSteps(time)

      # compute sigma_f
      sigma_f = self._sigma_f[i_left][group] + (time - self._time_steps[i_left]) * \
                (self._sigma_f[i_left][group] - self._sigma_f[i_right][group]) / \
                (self._time_steps[i_right] - self._time_steps[i_left])

      return sigma_f


  def getNuSigmaFByGroup(self, group, time, temp):

    # Check if necessary variables have been set
    check_set(self._num_energy_groups, 'Material ID={0} nu_sigma_f by group'.format(self._id), 'Num Energy Groups')
    check_is_int(group, 'Material ID={0} nu_sigma_f by group'.format(self._id), 'group')
    check_is_float_or_int(time, 'Material ID={0} nu_sigma_f by group'.format(self._id), 'time')
    check_is_float_or_int(group, 'Material ID={0} nu_sigma_f by group'.format(self._id), 'temp')

    # check if group is valid
    if group < 0 or group > self._num_energy_groups-1:
      msg = 'Unable to get nu fission xs for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    # check if time_step is valid
    elif time < self._time_steps[0] or time > self._time_steps[-1]:
      msg = 'Unable to get nu fission xs for Material ID={0} for group {1} and '\
            'time {2} outside time window [{3}, {4}]'\
          .format(self._id, group, time, self._time_steps[0], self._time_steps[-1])
      raise ValueError(msg)

    else:

      # find time step window
      (i_left, i_right) = self.getBoundingTimeSteps(time)

      # compute nu_sigma_f
      nu_sigma_f = self._nu_sigma_f[i_left][group] + (time - self._time_steps[i_left]) * \
                (self._nu_sigma_f[i_left][group] - self._nu_sigma_f[i_right][group]) / \
                (self._time_steps[i_right] - self._time_steps[i_left])


      return nu_sigma_f


  def getSigmaSByGroup(self, group_from, group_to, time, temp):

    # Check if necessary variables have been set
    check_set(self._num_energy_groups, 'Material ID={0} sigma_s by group'.format(self._id), 'Num Energy Groups')
    check_is_int(group_from, 'Material ID={0} sigma_s by group'.format(self._id), 'group_from')
    check_is_int(group_to, 'Material ID={0} sigma_s by group'.format(self._id), 'group_to')
    check_is_float_or_int(time, 'Material ID={0} sigma_s by group'.format(self._id), 'time')
    check_is_float_or_int(temp, 'Material ID={0} sigma_s by group'.format(self._id), 'temp')

    # check if group is valid
    if group_from < 0 or group_from > self._num_energy_groups-1 or \
    group_to < 0 or group_to > self._num_energy_groups-1:
      msg = 'Unable to get scattering xs for Material ID={0} for group {1} '\
          'to {2} as num energy groups is set to {3}'\
          .format(self._id, group_to, group_from, self._num_energy_groups)
      raise ValueError(msg)

    # check if time_step is valid
    elif time < self._time_steps[0] or time > self._time_steps[-1]:
      msg = 'Unable to get scattering xs for Material ID={0} for group {1} to {2} and '\
            'time {3} outside time window [{4}, {5}]'\
          .format(self._id, group_from, group_to, time, self._time_steps[0], self._time_steps[-1])
      raise ValueError(msg)

    else:

      ng = self._num_energy_groups

      # find time step window
      (i_left, i_right) = self.getBoundingTimeSteps(time)

      # compute sigma_s
      sigma_s = self._sigma_s[i_left][group_from*ng+group_to] + (time - self._time_steps[i_left]) * \
                (self._sigma_s[i_left][group_from*ng+group_to] - self._sigma_s[i_right][group_from*ng+group_to]) / \
                (self._time_steps[i_right] - self._time_steps[i_left])

      return sigma_s


  def getDifCoefByGroup(self, group, time, temp):

    # Check if necessary variables have been set
    check_set(self._num_energy_groups, 'Material ID={0} dif_coef by group'.format(self._id), 'Num Energy Groups')
    check_is_int(group, 'Material ID={0} dif_coef by group'.format(self._id), 'group')
    check_is_float_or_int(time, 'Material ID={0} dif_coef by group'.format(self._id), 'time')
    check_is_float_or_int(group, 'Material ID={0} dif_coef by group'.format(self._id), 'temp')

    # check if group is valid
    if group < 0 or group > self._num_energy_groups-1:
      msg = 'Unable to get dif coef for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    # check if time_step is valid
    elif time < self._time_steps[0] or time > self._time_steps[-1]:
      msg = 'Unable to get dif coef for Material ID={0} for group {1} and'\
            ' time {2} outside time window [{3}, {4}]'\
          .format(self._id, group, time, self._time_steps[0], self._time_steps[-1])
      raise ValueError(msg)

    else:

      # find time step window
      (i_left, i_right) = self.getBoundingTimeSteps(time)

      # compute dif_coef
      dif_coef = self._dif_coef[i_left][group] + (time - self._time_steps[i_left]) * \
                (self._dif_coef[i_left][group] - self._dif_coef[i_right][group]) / \
                (self._time_steps[i_right] - self._time_steps[i_left])

      return dif_coef


  def getChiByGroup(self, group, time, temp):

    # Check if necessary variables have been set
    check_set(self._num_energy_groups, 'Material ID={0} chi by group'.format(self._id), 'Num Energy Groups')
    check_is_int(group, 'Material ID={0} chi by group'.format(self._id), 'group')
    check_is_float_or_int(time, 'Material ID={0} chi by group'.format(self._id), 'time')
    check_is_float_or_int(group, 'Material ID={0} chi by group'.format(self._id), 'temp')

    # check if group is valid
    if group < 0 or group > self._num_energy_groups-1:
      msg = 'Unable to get chi for Material ID={0} for group {1} '\
          'as num energy groups is set to {2}'\
          .format(self._id, group, self._num_energy_groups)
      raise ValueError(msg)

    # check if time_step is valid
    elif time < self._time_steps[0] or time > self._time_steps[-1]:
      msg = 'Unable to get chi for Material ID={0} for group {1} and'\
            ' time {2} outside time window [{3}, {4}]'\
          .format(self._id, group, time, self._time_steps[0], self._time_steps[-1])
      raise ValueError(msg)

    else:

      # find time step window
      (i_left, i_right) = self.getBoundingTimeSteps(time)

      # compute chi
      chi = self._chi[i_left][group] + (time - self._time_steps[i_left]) * \
                (self._chi[i_left][group] - self._chi[i_right][group]) / \
                (self._time_steps[i_right] - self._time_steps[i_left])

      return chi


  def setNumTimeSteps(self, num_time_steps):

    # Check if necessary variables have been set
    check_is_int(num_time_steps, 'Material ID={0} num_time_steps'.format(self._id), 'num time steps')    

    if num_time_steps < 2:
      msg = 'Unable to set num time steps to {0} since it must be an ' \
            ' integer > 1'.format(num_time_steps)
      raise ValueError(msg)
      
    else:
      self._num_time_steps = num_time_steps
      self._time_steps = np.zeros(num_time_steps)
      

  def setNumEnergyGroups(self, num_energy_groups):

    # Check if necessary variables have been set
    check_set(self._num_time_steps, 'Material ID={0} num_energy_groups'.format(self._id), 'Num Time Steps')

    if is_integer(num_energy_groups):
   
      if num_energy_groups < 1:
        msg = 'Unable to set num energy groups to {0} since it must be an ' \
              ' integer > 0'.format(num_energy_groups)
        raise ValueError(msg)

      else:
        self._num_energy_groups = num_energy_groups
        self._sigma_t = np.zeros((self._num_time_steps, num_energy_groups))
        self._sigma_a = np.zeros((self._num_time_steps, num_energy_groups))
        self._sigma_f = np.zeros((self._num_time_steps, num_energy_groups))
        self._nu_sigma_f = np.zeros((self._num_time_steps, num_energy_groups))
        self._sigma_s = np.zeros((self._num_time_steps, num_energy_groups*num_energy_groups))
        self._dif_coef = np.zeros((self._num_time_steps, num_energy_groups))
        self._chi = np.zeros((self._num_time_steps, num_energy_groups))
        self._doppler_coefficients = np.zeros(num_energy_groups)
        self._velocity = np.zeros(num_energy_groups)
        self._energy_per_fission = np.zeros(num_energy_groups)

    else:
      msg = 'Unable to set num energy groups to non-integer {0}'\
          .format(num_energy_groups)
      raise ValueError(msg)


  def setNumDelayedGroups(self, num_delayed_groups):

    # Check if necessary variables have been set
    check_is_int(num_delayed_groups, 'Material ID={0} num_delayed_groups'.format(self._id), 'num delayed groups')    

    if num_delayed_groups < 1:
      msg = 'Unable to set num delayed groups to {0} since it must be an ' \
            ' integer > 0'.format(num_delayed_groups)
      raise ValueError(msg)
      
    else:
      self._num_delayed_groups = num_delayed_groups
      self._decay_constants = np.zeros(num_delayed_groups)
      self._delayed_fractions = np.zeros(num_delayed_groups)
      self._precursor_conc = {}
      for i in CLOCK_POSITIONS:
        self._precursor_conc[i] = np.zeros(num_delayed_groups)


  def setTimeSteps(self, time_steps):

    # Check if necessary variables have been set
    check_set(self._num_time_steps, 'Material ID={0} time_steps'.format(self._id), 'Num Time Steps')
    check_list_of_floats_or_ints(time_steps, 'Material ID={0} time_steps'.format(self._id))

    if isinstance(time_steps, list):
      self._time_steps = np.asarray(time_steps)
    else:
      self._time_steps = np.copy(time_steps)


  def setVelocity(self, velocity):

    # Check if necessary variables have been set
    check_set(self._num_energy_groups, 'Material ID={0} velocity'.format(self._id), 'Num Energy Groups')
    check_list_of_floats_or_ints(velocity, 'Material ID={0} velocity'.format(self._id))

    if isinstance(velocity, list):
      self._velocity = np.asarray(velocity)
    else:
      self._velocity = np.copy(velocity)


  def setDecayConstants(self, decay_constants):

    # Check if necessary variables have been set
    check_set(self._num_delayed_groups, 'Material ID={0} decay_constants'.format(self._id), 'Num Delayed Groups')
    check_list_of_floats_or_ints(decay_constants, 'Material ID={0} decay_constants'.format(self._id))

    if isinstance(decay_constants, list):
      self._decay_constants = np.asarray(decay_constants)
    else:
      self._decay_constants = np.copy(decay_constants)
      

  def setDelayedFractions(self, delayed_fractions):

    # Check if necessary variables have been set
    check_set(self._num_delayed_groups, 'Material ID={0} delayed_fractions'.format(self._id), 'Num Delayed Groups')
    check_list_of_floats_or_ints(delayed_fractions, 'Material ID={0} delayed_fractions'.format(self._id))

    if isinstance(delayed_fractions, list):
      self._delayed_fractions = np.asarray(delayed_fractions)
    else:
      self._delayed_fractions = np.copy(delayed_fractions)
      

  def getBoundingTimeSteps(self, time):

    # check if time_step is valid
    if time < self._time_steps[0] or time > self._time_steps[-1]:
      msg = 'Unable to get bounding time steps for Material ID={0} for'\
            ' time {2} outside time window [{3}, {4}]'\
          .format(self._id, group, self._num_time_steps, self._time_steps[0], self._time_steps[-1])
      raise ValueError(msg)

    else:

      i_left = 0
      while time > self._time_steps[i_left+1]:
        i_left += 1

      i_right = i_left + 1
      return (i_left, i_right)


  def __repr__(self):

    global CLOCK_POSITIONS

    string = 'OpenRK FunctionalMaterial\n'
    string += ' Name \t\t\t= {0} \n'.format(self._name)
    string += ' ID \t\t\t= {0} \n'.format(self._id)
    string += ' Num Energy Groups \t= {0} \n'.format(self._num_energy_groups)
    string += ' Num Delayed Groups \t= {0} \n'.format(self._num_delayed_groups)
    string += ' Num Steps \t\t= {0}\n'.format(self._num_time_steps)
    string += ' Time Steps \t\t= {0} \n'.format(self._time_steps)
    string += ' Energy Per Fission \t= {0} \n'.format(self._energy_per_fission)
    string += ' Doppler Coefficients \t= {0} \n'.format(self._doppler_coefficients)
    string += ' Velocity \t\t= {0} \n'.format(self._velocity)
    string += ' Decay Constants \t= {0} \n'.format(self._decay_constants)
    string += ' Delayed Fractions \t= {0} \n'.format(self._delayed_fractions)
    for i in range(self._num_time_steps):
      string += ' Time Step \t\t= {0} \n'.format(i)
      string += ' Sigma T \t\t= {0} \n'.format(self._sigma_t[i])
      string += ' Sigma A \t\t= {0} \n'.format(self._sigma_a[i])
      string += ' Sigma F \t\t= {0} \n'.format(self._sigma_f[i])
      string += ' Nu Sigma F \t\t= {0} \n'.format(self._nu_sigma_f[i])
      string += ' Dif Coef \t\t= {0} \n'.format(self._dif_coef[i])
      string += ' Chi \t\t\t= {0} \n'.format(self._chi[i])
      string += ' Sigma S \t\t= {0} \n'.format(self._sigma_s[i])
    for i in CLOCK_POSITIONS:
      string += ' Precursor Conc ({:^12s}) \t'.format(i)
      string += '= {0} \n'.format(self._precursor_conc[i])

    return string
