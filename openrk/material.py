__author__ = 'Samuel Shaner'
__email__ = 'shaner@mit.edu'

# Import modules
import numpy as np
from clock import Clock
import checkvalue as cv
import math

# A static variable for auto-generated Material UIDs
AUTO_MATERIAL_UID = 1


class Material(object):
    """
    Main Material class which contains information on the cross sections
    and other properties of a material

    Attributes:
      _sigma_a
      _sigma_f
      _sigma_t
      _sigma_s
      _nu_sigma_f
      _chi
      _dif_coef
      _velocity
      _energy_per_fission
      _num_energy_groups
      _name
      _id
      _is_fissionable

    Setter Methods:
      set_name(str name)
      set_energy_per_fission(float energy_per_fission)
      set_sigma_a(list/np.array sigma_a)
      set_sigma_t(list/np.array sigma_t)
      set_sigma_f(list/np.array sigma_f)
      set_sigma_s(list/np.array sigma_s)
      set_nu_sigma_f(list/np.array nu_sigma_f)
      set_chi(list/np.array chi)
      set_dif_coef(list/np.array dif_coef)
      set_velocity(list/np.array velocity)
      set_sigma_a_by_group(float sigma_a, int group)
      set_sigma_t_by_group(float sigma_t, int group)
      set_sigma_f_by_group(float sigma_f, int group)
      set_sigma_s_by_group(float sigma_s, int group_from, int group_to)
      set_nu_sigma_f_by_group(float nu_sigma_f, int group)
      set_chi_by_group(float chi, int group)
      set_dif_coef_by_group(float dif_coef, int group)
      set_velocity_by_group(float velocity, int group)
      set_num_energy_groups(int num_energy_groups)
      set_id(int id)
      set_is_fissionable(bool is_fissionable)

    Getter Methods:
      get_energy_per_fission()
      get_sigma_a_by_group(int group)
      get_sigma_t_by_group(int group)
      get_sigma_f_by_group(int group)
      get_sigma_s_by_group(int group_from, int group_to)
      get_nu_sigma_f_by_group(int group)
      get_chi_by_group(int group)
      get_dif_coef_by_group(int group)
      get_velocity_by_group(int group)
      get_num_energy_groups()
      get_id()
      get_name()
      get_is_fissionable()

    Other Methods:

    """
    def __init__(self, material_id=None, name=''):

        # Initialize class attributes
        global AUTO_MATERIAL_UID
        self._uid = AUTO_MATERIAL_UID
        AUTO_MATERIAL_UID += 1
        self._id = None
        self._set_id = False
        self._name = ''
        self._is_fissionable = False

        # Initialize cross section attributes
        self._sigma_t = None
        self._sigma_a = None
        self._sigma_f = None
        self._nu_sigma_f = None
        self._sigma_s = None
        self._dif_coef = None
        self._velocity = None
        self._chi = None
        self._num_energy_groups = None
        self._energy_per_fission = None

        # Set the Material class attributes
        self.set_name(name)

        if material_id is not None:
            self.set_id(material_id)
        else:
            self.set_id(self._uid)
            self._set_id = False

    def set_is_fissionable(self, is_fissionable):

        cv.check_is_bool(is_fissionable, 'Material is fissionable', 'is_fissionable')
        self._is_fissionable = is_fissionable

    def get_is_fissionable(self):

        return self._is_fissionable

    def set_name(self, name):

        if not cv.is_string(name):
            msg = 'Unable to set name for Material ID={0} with a non-string ' \
                  'value {1}'.format(self._id, name)
            raise ValueError(msg)

        else:
            self._name = name

    def get_name(self):

        return self._name

    def set_energy_per_fission(self, energy_per_fission):

        # Check if necessary variables have been set
        self._energy_per_fission = energy_per_fission

    def get_energy_per_fission(self):

        # Check if necessary variables have been set
        cv.check_set(self._energy_per_fission, 'Energy Per Fission', 'energy per fission')

        return self._energy_per_fission

    def set_sigma_a(self, sigma_a):

        # check if sigma_a is a list
        if not cv.is_list(sigma_a):
            msg = 'Unable to set absorption xs for Material ID={0} with a ' \
                  'non-list value {1}'.format(self._id, sigma_a)
            raise ValueError(msg)

        # check if sigma_a is of length num_groups
        elif len(sigma_a) != self._num_energy_groups:
            msg = 'Unable to set absorption xs for Material ID={0} with {1} groups ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, len(sigma_a), self._num_energy_groups)
            raise ValueError(msg)

        else:
            if isinstance(sigma_a, list):
                sigma_a = np.asarray(sigma_a)

            np.copyto(self._sigma_a, sigma_a)

    def set_sigma_t(self, sigma_t):

        # check if sigma_t is a list
        if not cv.is_list(sigma_t):
            msg = 'Unable to set total xs for Material ID={0} with a ' \
                  'non-list value {1}'.format(self._id, sigma_t)
            raise ValueError(msg)

        # check if sigma_t is of length num_groups
        elif len(sigma_t) != self._num_energy_groups:
            msg = 'Unable to set total xs for Material ID={0} with {1} groups ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, len(sigma_t), self._num_energy_groups)
            raise ValueError(msg)

        else:
            if isinstance(sigma_t, list):
                sigma_t = np.asarray(sigma_t)

            np.copyto(self._sigma_t, sigma_t)

    def set_sigma_f(self, sigma_f):

        # check if sigma_f is a list
        if not cv.is_list(sigma_f):
            msg = 'Unable to set fission xs for Material ID={0} with a ' \
                  'non-list value {1}'.format(self._id, sigma_f)
            raise ValueError(msg)

        # check if sigma_f is of length num_groups
        elif len(sigma_f) != self._num_energy_groups:
            msg = 'Unable to set fission xs for Material ID={0} with {1} groups ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, len(sigma_f), self._num_energy_groups)
            raise ValueError(msg)

        else:
            if isinstance(sigma_f, list):
                sigma_f = np.asarray(sigma_f)

            np.copyto(self._sigma_f, sigma_f)

    def set_nu_sigma_f(self, nu_sigma_f):

        # check if nu_sigma_f is a list
        if not cv.is_list(nu_sigma_f):
            msg = 'Unable to set nu fission xs for Material ID={0} with a ' \
                  'non-list value {1}'.format(self._id, nu_sigma_f)
            raise ValueError(msg)

        # check if nu_sigma_f is of length num_groups
        elif len(nu_sigma_f) != self._num_energy_groups:
            msg = 'Unable to set nu fission xs for Material ID={0} with {1} groups ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, len(nu_sigma_f), self._num_energy_groups)
            raise ValueError(msg)

        else:
            if isinstance(nu_sigma_f, list):
                nu_sigma_f = np.asarray(nu_sigma_f)

            np.copyto(self._nu_sigma_f, nu_sigma_f)

            if np.count_nonzero(self._nu_sigma_f):
                self.set_is_fissionable(True)

    def set_sigma_s(self, sigma_s):

        # check if sigma_s is a list
        if not cv.is_list(sigma_s):
            msg = 'Unable to set scattering xs for Material ID={0} with a ' \
                  'non-list value {1}'.format(self._id, sigma_s)
            raise ValueError(msg)

        # check if sigma_s is of length num_energy_groups*num_energy_groups
        elif len(sigma_s) != self._num_energy_groups * self._num_energy_groups:
            msg = 'Unable to set scattering xs for Material ID={0} with {1} groups ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, math.sqrt(len(sigma_s)), self._num_energy_groups)
            raise ValueError(msg)

        else:
            if isinstance(sigma_s, list):
                sigma_s = np.asarray(sigma_s)

            np.copyto(self._sigma_s, sigma_s)

    def set_dif_coef(self, dif_coef):

        # check if dif_coef is a list
        if not cv.is_list(dif_coef):
            msg = 'Unable to set the diffusion coefficients for Material ID={0} ' \
                  'with a non-list value {1}'.format(self._id, dif_coef)
            raise ValueError(msg)

        # check if dif_coef is of length num_energy_groups
        elif len(dif_coef) != self._num_energy_groups:
            msg = 'Unable to set the diffusion coefficients for Material ID={0} with' \
                  ' {1} groups as num energy groups is set to {2}' \
                .format(self._id, len(dif_coef), self._num_energy_groups)
            raise ValueError(msg)

        else:
            if isinstance(dif_coef, list):
                dif_coef = np.asarray(dif_coef)

            np.copyto(self._dif_coef, dif_coef)

    def set_chi(self, chi):

        # check if chi is a list
        if not cv.is_list(chi):
            msg = 'Unable to set the chi for Material ID={0} ' \
                  'with a non-list value {1}'.format(self._id, chi)
            raise ValueError(msg)

        # check if chi is of length num_energy_groups
        elif len(chi) != self._num_energy_groups:
            msg = 'Unable to set the chi for Material ID={0} with' \
                  ' {1} groups as num energy groups is set to {2}' \
                .format(self._id, len(chi), self._num_energy_groups)
            raise ValueError(msg)

        else:
            if isinstance(chi, list):
                chi = np.asarray(chi)

            np.copyto(self._chi, chi)

    def set_sigma_a_by_group(self, sigma_a, group):

        # check if sigma_a is a float
        if not cv.is_float(sigma_a):
            msg = 'Unable to set absorption xs for Material ID={0} with a ' \
                  'non-float value {1}'.format(self._id, sigma_a)
            raise ValueError(msg)

        # check if group is valid
        elif group < 0 or group > self._num_energy_groups:
            msg = 'Unable to set absorption xs for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            self._sigma_a[group] = sigma_a

    def set_sigma_t_by_group(self, sigma_t, group):

        # check if sigma_t is a float
        if not cv.is_float(sigma_t):
            msg = 'Unable to set total xs for Material ID={0} with a ' \
                  'non-float value {1}'.format(self._id, sigma_t)
            raise ValueError(msg)

        # check if group is valid
        elif group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to set total xs for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            self._sigma_t[group] = sigma_t

    def set_sigma_f_by_group(self, sigma_f, group):

        # check if sigma_f is a float
        if not cv.is_float(sigma_f):
            msg = 'Unable to set fission xs for Material ID={0} with a ' \
                  'non-float value {1}'.format(self._id, sigma_f)
            raise ValueError(msg)

        # check if group is valid
        elif group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to set fission xs for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            self._sigma_f[group] = sigma_f

    def set_nu_sigma_f_by_group(self, nu_sigma_f, group):

        # check if nu_sigma_f is a list
        if not cv.is_float(nu_sigma_f):
            msg = 'Unable to set nu fission xs for Material ID={0} with a ' \
                  'non-float value {1}'.format(self._id, nu_sigma_f)
            raise ValueError(msg)

        # check if group is valid
        elif group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to set fission xs for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            self._nu_sigma_f[group] = nu_sigma_f

            if np.count_nonzero(self._nu_sigma_f):
                self.set_is_fissionable(True)

    def set_sigma_s_by_group(self, sigma_s, group_from, group_to):

        # check if sigma_s is a list
        if not cv.is_float(sigma_s):
            msg = 'Unable to set scattering xs for Material ID={0} with a ' \
                  'non-float value {1}'.format(self._id, sigma_s)
            raise ValueError(msg)

        # check if group is valid
        elif group_from < 0 or group_from > self._num_energy_groups - 1 or \
                group_to < 0 or group_to > self._num_energy_groups - 1:
            msg = 'Unable to set scattering xs for Material ID={0} for group {1} ' \
                'to {2} as num energy groups is set to {3}' \
                .format(self._id, group_to, group_from, self._num_energy_groups)
            raise ValueError(msg)

        else:
            self._sigma_s[group_from * self._num_energy_groups + group_to] = sigma_s

    def set_dif_coef_by_group(self, dif_coef, group):

        # check if dif_coef is a list
        if not cv.is_float(dif_coef):
            msg = 'Unable to set dif coef for Material ID={0} with a ' \
                  'non-float value {1}'.format(self._id, dif_coef)
            raise ValueError(msg)

        # check if group is valid
        elif group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to set dif coef for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            self._dif_coef[group] = dif_coef

    def set_chi_by_group(self, chi, group):

        # check if chi is a list
        if not cv.is_float(chi):
            msg = 'Unable to set chi for Material ID={0} with a ' \
                  'non-float value {1}'.format(self._id, chi)
            raise ValueError(msg)

        # check if group is valid
        elif group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to set chi for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            self._chi[group] = chi

    def get_sigma_a_by_group(self, group, time=None, temp=None):

        # check if group is valid
        if group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to get absorption xs for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            return self._sigma_a[group]

    def get_sigma_t_by_group(self, group, time=None, temp=None):

        # check if group is valid
        if group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to get total xs for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            return self._sigma_t[group]

    def get_sigma_f_by_group(self, group, time=None, temp=None):

        # check if group is valid
        if group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to get fission xs for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            return self._sigma_f[group]

    def get_nu_sigma_f_by_group(self, group, time=None, temp=None):

        # check if group is valid
        if group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to get fission xs for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            return self._nu_sigma_f[group]

    def get_sigma_s_by_group(self, group_from, group_to, time=None, temp=None):

        # check if group is valid
        if group_from < 0 or group_from > self._num_energy_groups - 1 or \
                group_to < 0 or group_to > self._num_energy_groups - 1:
            msg = 'Unable to get scattering xs for Material ID={0} for group {1} ' \
                  'to {2} as num energy groups is set to {3}' \
                .format(self._id, group_to, group_from, self._num_energy_groups)
            raise ValueError(msg)

        else:
            return self._sigma_s[group_from * self._num_energy_groups + group_to]

    def get_dif_coef_by_group(self, group, time=None, temp=None):

        # check if group is valid
        if group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to set dif coef for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            return self._dif_coef[group]

    def get_chi_by_group(self, group, time=None, temp=None):

        # check if group is valid
        if group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to set chi for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            return self._chi[group]

    def get_num_energy_groups(self):

        cv.check_set(self._num_energy_groups, 'Material ID={0} get num energy groups'.
                     format(self._id), 'num energy groups')

        return self._num_energy_groups

    def set_num_energy_groups(self, num_energy_groups):

        if cv.is_integer(num_energy_groups):

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
                self._sigma_s = np.zeros(num_energy_groups * num_energy_groups)
                self._dif_coef = np.zeros(num_energy_groups)
                self._chi = np.zeros(num_energy_groups)
                self._velocity = np.zeros(num_energy_groups)

        else:
            msg = 'Unable to set num energy groups to non-integer {0}' \
                .format(num_energy_groups)
            raise ValueError(msg)

    def set_id(self, material_id):

        # Check that the ID is a non-negative integer
        if cv.is_integer(material_id):

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

    def get_velocity_by_group(self, group):

        # Check if necessary variables have been set
        cv.check_set(self._num_energy_groups, 'Material ID={0} velocity by group'.format(self._id), 'Num Energy Groups')
        cv.check_is_int(group, 'Material ID={0} velocity by group'.format(self._id), 'group')

        # check if group is valid
        if group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to get velocity for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:

            return self._velocity[group]

    def set_velocity(self, velocity):

        # Check if necessary variables have been set
        cv.check_set(self._num_energy_groups, 'Material ID={0} velocity'.format(self._id), 'Num Energy Groups')
        cv.check_list_of_floats_or_ints(velocity, 'Material ID={0} velocity'.format(self._id))

        if isinstance(velocity, list):
            velocity = np.asarray(velocity)

        np.copyto(self._velocity, velocity)

    def set_velocity_by_group(self, velocity, group):

        # check if velocity is a float
        if not cv.is_float(velocity):
            msg = 'Unable to set velocity for Material ID={0} with a ' \
                  'non-float value {1}'.format(self._id, velocity)
            raise ValueError(msg)

        # check if group is valid
        elif group < 0 or group > self._num_energy_groups:
            msg = 'Unable to set velocity for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            self._velocity[group] = velocity

    def get_id(self):

        return self._id

    def clone(self):

        new_material = Material(name=self._name)
        new_material.set_num_energy_groups(self._num_energy_groups)
        new_material.set_sigma_a(self._sigma_a)
        new_material.set_sigma_t(self._sigma_t)
        new_material.set_sigma_f(self._sigma_f)
        new_material.set_nu_sigma_f(self._nu_sigma_f)
        new_material.set_sigma_s(self._sigma_s)
        new_material.set_dif_coef(self._dif_coef)
        new_material.set_chi(self._chi)
        new_material.set_velocity(self._velocity)
        new_material.set_energy_per_fission(self._energy_per_fission)

        return new_material

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
        string += ' Energy Per Fission \t= {0} \n'.format(self._energy_per_fission)

        return string


class TransientMaterial(Material):
    """
    Extension of Material class that adds properties that describe transient materials

    Attributes:
      _sigma_a
      _sigma_f
      _sigma_t
      _sigma_s
      _nu_sigma_f
      _chi
      _dif_coef
      _velocity
      _energy_per_fission
      _num_energy_groups
      _name
      _id
      _decay_constants
      _delayed_fractions
      _precursor_conc
      _clock
      _num_delayed_groups

    Setter Methods:
      set_name(str name)
      set_energy_per_fission(float energy_per_fission)
      set_sigma_a(list/np.array sigma_a, time='CURRENT')
      set_sigma_t(list/np.array sigma_t, time='CURRENT')
      set_sigma_f(list/np.array sigma_f, time='CURRENT')
      set_sigma_s(list/np.array sigma_s, time='CURRENT')
      set_nu_sigma_f(list/np.array nu_sigma_f, time='CURRENT')
      set_chi(list/np.array chi, time='CURRENT')
      set_dif_coef(list/np.array dif_coef, time='CURRENT')
      set_velocity(list/np.array velocity, time='CURRENT')
      set_sigma_a_by_group(float sigma_a, int group, time='CURRENT')
      set_sigma_t_by_group(float sigma_t, int group, time='CURRENT')
      set_sigma_f_by_group(float sigma_f, int group, time='CURRENT')
      set_sigma_s_by_group(float sigma_s, int group_from, int group_to, time='CURRENT')
      set_nu_sigma_f_by_group(float nu_sigma_f, int group, time='CURRENT')
      set_chi_by_group(float chi, int group, time='CURRENT')
      set_dif_coef_by_group(float dif_coef, int group, time='CURRENT')
      set_velocity_by_group(float velocity, int group, time='CURRENT')
      set_num_energy_groups(int num_energy_groups)
      set_id(int id)
      set_decay_contants(list/np.array decay_constants, time='CURRENT')
      set_delayed_fractions(list/np.array delayed_fractions, time='CURRENT')
      set_precursor_conc(list/np.array precursor_conc, time='CURRENT')
      set_decay_contant_by_group(float decay_constant, int group, time='CURRENT')
      set_delayed_fraction_by_group(float delayed_fraction, int group, time='CURRENT')
      set_precursor_conc_by_group(float precursor_conc, int group, time='CURRENT')
      set_clock(Clock clock)
      set_num_delayed_groups(int num_groups)

    Getter Methods:
      get_energy_per_fission()
      get_sigma_a_by_group(int group, time='CURRENT')
      get_sigma_t_by_group(int group, time='CURRENT')
      get_sigma_f_by_group(int group, time='CURRENT')
      get_sigma_s_by_group(int group_from, int group_to, time='CURRENT')
      get_nu_sigma_f_by_group(int group, time='CURRENT')
      get_chi_by_group(int group, time='CURRENT')
      get_dif_coef_by_group(int group, time='CURRENT')
      get_velocity_by_group(int group, time='CURRENT')
      get_num_delayed_groups()
      get_id()
      get_name()
      get_decay_contant_by_group(int group, time='CURRENT')
      get_delayed_fraction_by_group(int group, time='CURRENT')
      get_precursor_conc_by_group(int group, time='CURRENT')

    Other Methods:

    """
    def __init__(self, material_id=None, name='', clock=None):

        super(TransientMaterial, self).__init__(material_id, name)

        # Initialize class attributes
        self._decay_constants = None
        self._delayed_fractions = None
        self._precursor_conc = None
        self._clock = None
        self._num_delayed_groups = None

        # Initialize clock
        if clock is None:
            self.set_clock(Clock())
        else:
            self.set_clock(clock)

    def set_clock(self, clock):

        # Initialize clock
        if not isinstance(clock, Clock):
            msg = 'Unable to initialize FunctionalMaterial clock since clock input is not of type ' \
                  'Clock: {0}'.format(clock)
            raise ValueError(msg)
        else:
            self._clock = clock

    def get_num_delayed_groups(self):

        cv.check_set(self._num_delayed_groups, 'TransientMaterial ID={0} get num delayed groups'.
                     format(self._id), 'num delayed groups')

        return self._num_delayed_groups

    def set_sigma_a(self, sigma_a, time='CURRENT'):

        # check if sigma_a is a list
        if not cv.is_list(sigma_a):
            msg = 'Unable to set absorption xs for Material ID={0} with a ' \
                  'non-list value {1}'.format(self._id, sigma_a)
            raise ValueError(msg)

        # check if sigma_a is of length num_groups
        elif len(sigma_a) != self._num_energy_groups:
            msg = 'Unable to set absorption xs for Material ID={0} with {1} groups ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, len(sigma_a), self._num_energy_groups)
            raise ValueError(msg)

        else:
            if isinstance(sigma_a, list):
                sigma_a = np.asarray(sigma_a)

            np.copyto(self._sigma_a[time], sigma_a)

    def set_sigma_t(self, sigma_t, time='CURRENT'):

        # check if sigma_t is a list
        if not cv.is_list(sigma_t):
            msg = 'Unable to set total xs for Material ID={0} with a ' \
                  'non-list value {1}'.format(self._id, sigma_t)
            raise ValueError(msg)

        # check if sigma_t is of length num_groups
        elif len(sigma_t) != self._num_energy_groups:
            msg = 'Unable to set total xs for Material ID={0} with {1} groups ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, len(sigma_t), self._num_energy_groups)
            raise ValueError(msg)

        else:
            if isinstance(sigma_t, list):
                sigma_t = np.asarray(sigma_t)

            np.copyto(self._sigma_t[time], sigma_t)

    def set_sigma_f(self, sigma_f, time='CURRENT'):

        # check if sigma_f is a list
        if not cv.is_list(sigma_f):
            msg = 'Unable to set fission xs for Material ID={0} with a ' \
                  'non-list value {1}'.format(self._id, sigma_f)
            raise ValueError(msg)

        # check if sigma_f is of length num_groups
        elif len(sigma_f) != self._num_energy_groups:
            msg = 'Unable to set fission xs for Material ID={0} with {1} groups ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, len(sigma_f), self._num_energy_groups)
            raise ValueError(msg)

        else:
            if isinstance(sigma_f, list):
                sigma_f = np.asarray(sigma_f)

            np.copyto(self._sigma_f[time], sigma_f)

    def set_nu_sigma_f(self, nu_sigma_f, time='CURRENT'):

        # check if nu_sigma_f is a list
        if not cv.is_list(nu_sigma_f):
            msg = 'Unable to set nu fission xs for Material ID={0} with a ' \
                  'non-list value {1}'.format(self._id, nu_sigma_f)
            raise ValueError(msg)

        # check if nu_sigma_f is of length num_groups
        elif len(nu_sigma_f) != self._num_energy_groups:
            msg = 'Unable to set nu fission xs for Material ID={0} with {1} groups ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, len(nu_sigma_f), self._num_energy_groups)
            raise ValueError(msg)

        else:
            if isinstance(nu_sigma_f, list):
                nu_sigma_f = np.asarray(nu_sigma_f)

            np.copyto(self._nu_sigma_f[time], nu_sigma_f)

            if np.count_nonzero(self._nu_sigma_f[time]):
                self.set_is_fissionable(True)

    def set_sigma_s(self, sigma_s, time='CURRENT'):

        # check if sigma_s is a list
        if not cv.is_list(sigma_s):
            msg = 'Unable to set scattering xs for Material ID={0} with a ' \
                  'non-list value {1}'.format(self._id, sigma_s)
            raise ValueError(msg)

        # check if sigma_s is of length num_energy_groups*num_energy_groups
        elif len(sigma_s) != self._num_energy_groups * self._num_energy_groups:
            msg = 'Unable to set scattering xs for Material ID={0} with {1} groups ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, math.sqrt(len(sigma_s)), self._num_energy_groups)
            raise ValueError(msg)

        else:
            if isinstance(sigma_s, list):
                sigma_s = np.asarray(sigma_s)

            np.copyto(self._sigma_s[time], sigma_s)

    def set_dif_coef(self, dif_coef, time='CURRENT'):

        # check if dif_coef is a list
        if not cv.is_list(dif_coef):
            msg = 'Unable to set the diffusion coefficients for Material ID={0} ' \
                  'with a non-list value {1}'.format(self._id, dif_coef)
            raise ValueError(msg)

        # check if dif_coef is of length num_energy_groups
        elif len(dif_coef) != self._num_energy_groups:
            msg = 'Unable to set the diffusion coefficients for Material ID={0} with' \
                  ' {1} groups as num energy groups is set to {2}' \
                .format(self._id, len(dif_coef), self._num_energy_groups)
            raise ValueError(msg)

        else:
            if isinstance(dif_coef, list):
                dif_coef = np.asarray(dif_coef)

            np.copyto(self._dif_coef[time], dif_coef)

    def set_chi(self, chi, time='CURRENT'):

        # check if chi is a list
        if not cv.is_list(chi):
            msg = 'Unable to set the chi for Material ID={0} ' \
                  'with a non-list value {1}'.format(self._id, chi)
            raise ValueError(msg)

        # check if chi is of length num_energy_groups
        elif len(chi) != self._num_energy_groups:
            msg = 'Unable to set the chi for Material ID={0} with' \
                  ' {1} groups as num energy groups is set to {2}' \
                .format(self._id, len(chi), self._num_energy_groups)
            raise ValueError(msg)

        else:
            if isinstance(chi, list):
                chi = np.asarray(chi)

            np.copyto(self._chi[time], chi)

    def set_sigma_a_by_group(self, sigma_a, group, time='CURRENT'):

        # check if sigma_a is a float
        if not cv.is_float(sigma_a):
            msg = 'Unable to set absorption xs for Material ID={0} with a ' \
                  'non-float value {1}'.format(self._id, sigma_a)
            raise ValueError(msg)

        # check if group is valid
        elif group < 0 or group > self._num_energy_groups:
            msg = 'Unable to set absorption xs for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            self._sigma_a[time][group] = sigma_a

    def set_sigma_t_by_group(self, sigma_t, group, time='CURRENT'):

        # check if sigma_t is a float
        if not cv.is_float(sigma_t):
            msg = 'Unable to set total xs for Material ID={0} with a ' \
                  'non-float value {1}'.format(self._id, sigma_t)
            raise ValueError(msg)

        # check if group is valid
        elif group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to set total xs for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            self._sigma_t[time][group] = sigma_t

    def set_sigma_f_by_group(self, sigma_f, group, time='CURRENT'):

        # check if sigma_f is a float
        if not cv.is_float(sigma_f):
            msg = 'Unable to set fission xs for Material ID={0} with a ' \
                  'non-float value {1}'.format(self._id, sigma_f)
            raise ValueError(msg)

        # check if group is valid
        elif group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to set fission xs for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            self._sigma_f[time][group] = sigma_f

    def set_nu_sigma_f_by_group(self, nu_sigma_f, group, time='CURRENT'):

        # check if nu_sigma_f is a list
        if not cv.is_float(nu_sigma_f):
            msg = 'Unable to set nu fission xs for Material ID={0} with a ' \
                  'non-float value {1}'.format(self._id, nu_sigma_f)
            raise ValueError(msg)

        # check if group is valid
        elif group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to set fission xs for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            self._nu_sigma_f[time][group] = nu_sigma_f

            if np.count_nonzero(self._nu_sigma_f[time]):
                self.set_is_fissionable(True)

    def set_sigma_s_by_group(self, sigma_s, group_from, group_to, time='CURRENT'):

        # check if sigma_s is a list
        if not cv.is_float(sigma_s):
            msg = 'Unable to set scattering xs for Material ID={0} with a ' \
                  'non-float value {1}'.format(self._id, sigma_s)
            raise ValueError(msg)

        # check if group is valid
        elif group_from < 0 or group_from > self._num_energy_groups - 1 or \
                group_to < 0 or group_to > self._num_energy_groups - 1:
            msg = 'Unable to set scattering xs for Material ID={0} for group {1} ' \
                'to {2} as num energy groups is set to {3}' \
                .format(self._id, group_to, group_from, self._num_energy_groups)
            raise ValueError(msg)

        else:
            self._sigma_s[time][group_from * self._num_energy_groups + group_to] = sigma_s

    def set_dif_coef_by_group(self, dif_coef, group, time='CURRENT'):

        # check if dif_coef is a list
        if not cv.is_float(dif_coef):
            msg = 'Unable to set dif coef for Material ID={0} with a ' \
                  'non-float value {1}'.format(self._id, dif_coef)
            raise ValueError(msg)

        # check if group is valid
        elif group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to set dif coef for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            self._dif_coef[time][group] = dif_coef

    def set_chi_by_group(self, chi, group, time='CURRENT'):

        # check if chi is a list
        if not cv.is_float(chi):
            msg = 'Unable to set chi for Material ID={0} with a ' \
                  'non-float value {1}'.format(self._id, chi)
            raise ValueError(msg)

        # check if group is valid
        elif group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to set chi for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            self._chi[time][group] = chi

    def get_sigma_a_by_group(self, group, time='CURRENT', temp=None):

        # check if group is valid
        if group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to get absorption xs for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            return self._sigma_a[time][group]

    def get_sigma_t_by_group(self, group, time='CURRENT', temp=None):

        # check if group is valid
        if group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to get total xs for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            return self._sigma_t[time][group]

    def get_sigma_f_by_group(self, group, time='CURRENT', temp=None):

        # check if group is valid
        if group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to get fission xs for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            return self._sigma_f[time][group]

    def get_nu_sigma_f_by_group(self, group, time='CURRENT', temp=None):

        # check if group is valid
        if group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to get fission xs for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            return self._nu_sigma_f[time][group]

    def get_sigma_s_by_group(self, group_from, group_to, time='CURRENT', temp=None):

        # check if group is valid
        if group_from < 0 or group_from > self._num_energy_groups - 1 or \
                group_to < 0 or group_to > self._num_energy_groups - 1:
            msg = 'Unable to get scattering xs for Material ID={0} for group {1} ' \
                  'to {2} as num energy groups is set to {3}' \
                .format(self._id, group_to, group_from, self._num_energy_groups)
            raise ValueError(msg)

        else:
            return self._sigma_s[time][group_from * self._num_energy_groups + group_to]

    def get_dif_coef_by_group(self, group, time='CURRENT', temp=None):

        # check if group is valid
        if group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to set dif coef for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            return self._dif_coef[time][group]

    def get_chi_by_group(self, group, time='CURRENT', temp=None):

        # check if group is valid
        if group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to set chi for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            return self._chi[time][group]

    def set_num_energy_groups(self, num_energy_groups):

        if cv.is_integer(num_energy_groups):

            if num_energy_groups < 1:
                msg = 'Unable to set num energy groups to {0} since it must be an ' \
                      ' integer > 0'.format(num_energy_groups)
                raise ValueError(msg)

            else:
                self._num_energy_groups = num_energy_groups
                self._sigma_t = {}
                self._sigma_a = {}
                self._sigma_f = {}
                self._nu_sigma_f = {}
                self._sigma_s = {}
                self._dif_coef = {}
                self._chi = {}
                self._velocity = {}
                for i in self._clock.get_positions():
                    self._sigma_t[i] = np.zeros(num_energy_groups)
                    self._sigma_a[i] = np.zeros(num_energy_groups)
                    self._sigma_f[i] = np.zeros(num_energy_groups)
                    self._nu_sigma_f[i] = np.zeros(num_energy_groups)
                    self._sigma_s[i] = np.zeros(num_energy_groups * num_energy_groups)
                    self._dif_coef[i] = np.zeros(num_energy_groups)
                    self._chi[i] = np.zeros(num_energy_groups)
                    self._velocity[i] = np.zeros(num_energy_groups)

        else:
            msg = 'Unable to set num energy groups to non-integer {0}' \
                .format(num_energy_groups)
            raise ValueError(msg)

    def get_velocity_by_group(self, group, time='CURRENT'):

        # Check if necessary variables have been set
        cv.check_set(self._num_energy_groups, 'Material ID={0} velocity by group'.format(self._id), 'Num Energy Groups')
        cv.check_is_int(group, 'Material ID={0} velocity by group'.format(self._id), 'group')

        # check if group is valid
        if group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to get velocity for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:

            return self._velocity[time][group]

    def set_velocity(self, velocity, time='CURRENT'):

        # Check if necessary variables have been set
        cv.check_set(self._num_energy_groups, 'Material ID={0} velocity'.format(self._id), 'Num Energy Groups')
        cv.check_list_of_floats_or_ints(velocity, 'Material ID={0} velocity'.format(self._id))

        if isinstance(velocity, list):
            velocity = np.asarray(velocity)

        np.copyto(self._velocity[time], velocity)

    def set_velocity_by_group(self, velocity, group, time='CURRENT'):

        # check if velocity is a float
        if not cv.is_float(velocity):
            msg = 'Unable to set velocity for Material ID={0} with a ' \
                  'non-float value {1}'.format(self._id, velocity)
            raise ValueError(msg)

        # check if group is valid
        elif group < 0 or group > self._num_energy_groups:
            msg = 'Unable to set velocity for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            self._velocity[time][group] = velocity

    def get_delayed_fraction_by_group(self, group, time='CURRENT'):

        # Check if necessary variables have been set
        cv.check_set(self._num_delayed_groups, 'Delayed Fraction', 'Num Delayed Groups')
        cv.check_is_int(group, 'Material ID={0} delayed fraction by group'.format(self._id), 'group')
        cv.check_clock_position(time, 'Material ID={0} delayed fraction by group'.format(self._id))

        # check if group is valid
        if group < 0 or group > self._num_delayed_groups - 1:
            msg = 'Unable to get delayed fraction for Material ID={0} for group {1} ' \
                  'as num delayed groups is set to {2}' \
                .format(self._id, group, self._num_delayed_groups)
            raise ValueError(msg)

        else:

            return self._delayed_fractions[time][group]

    def get_decay_constant_by_group(self, group, time='CURRENT'):

        # Check if necessary variables have been set
        cv.check_set(self._num_delayed_groups, 'Decay constant', 'Num Energy Groups')
        cv.check_is_int(group, 'Material ID={0} decay constant by group'.format(self._id), 'group')
        cv.check_clock_position(time, 'Material ID={0} decay constant by group'.format(self._id))

        # check if group is valid
        if group < 0 or group > self._num_delayed_groups - 1:
            msg = 'Unable to get decay constant for Material ID={0} for group {1} ' \
                  'as num delayed groups is set to {2}' \
                .format(self._id, group, self._num_delayed_groups)
            raise ValueError(msg)

        else:

            return self._decay_constants[time][group]

    def get_precursor_conc_by_group(self, group, time='CURRENT'):

        # Check if necessary variables have been set
        cv.check_set(self._num_delayed_groups, 'Precursor Concentration', 'Num Energy Groups')
        cv.check_is_int(group, 'Material ID={0} precursor concentration by group'.format(self._id), 'group')
        cv.check_clock_position(time, 'Material ID={0} precursor concentration by group'.format(self._id))

        # check if group is valid
        if group < 0 or group > self._num_delayed_groups - 1:
            msg = 'Unable to get precursor concentration for Material ID={0} for group {1} ' \
                  'as num delayed groups is set to {2}' \
                .format(self._id, group, self._num_delayed_groups)
            raise ValueError(msg)

        else:

            return self._precursor_conc[time][group]

    def set_num_delayed_groups(self, num_delayed_groups):

        # Check if necessary variables have been set
        cv.check_is_int(num_delayed_groups, 'Material ID={0} num_delayed_groups'.format(self._id), 'num delayed groups')

        if num_delayed_groups < 1:
            msg = 'Unable to set num delayed groups to {0} since it must be an ' \
                  ' integer > 0'.format(num_delayed_groups)
            raise ValueError(msg)

        else:
            self._num_delayed_groups = num_delayed_groups
            self._decay_constants = {}
            self._delayed_fractions = {}
            self._precursor_conc = {}
            for i in self._clock.get_positions():
                self._decay_constants[i] = np.zeros(num_delayed_groups)
                self._delayed_fractions[i] = np.zeros(num_delayed_groups)
                self._precursor_conc[i] = np.zeros(num_delayed_groups)

    def set_decay_constants(self, decay_constants, time='CURRENT'):

        # Check if necessary variables have been set
        cv.check_set(self._num_delayed_groups, 'Material ID={0} decay_constants'.format(self._id), 'Num Delayed Groups')
        cv.check_list_of_floats_or_ints(decay_constants, 'Material ID={0} decay_constants'.format(self._id))
        cv.check_clock_position(time, 'Material ID={0} decay constants'.format(self._id))

        if isinstance(decay_constants, list):
            decay_constants = np.asarray(decay_constants)

        np.copyto(self._decay_constants[time], decay_constants)

    def set_delayed_fractions(self, delayed_fractions, time='CURRENT'):

        # Check if necessary variables have been set
        cv.check_set(self._num_delayed_groups, 'Material ID={0} delayed_fractions'.format(self._id),
                     'Num Delayed Groups')
        cv.check_list_of_floats_or_ints(delayed_fractions, 'Material ID={0} delayed_fractions'.format(self._id))
        cv.check_clock_position(time, 'Material ID={0} delayed fractions'.format(self._id))

        if isinstance(delayed_fractions, list):
            delayed_fractions = np.asarray(delayed_fractions)

        np.copyto(self._delayed_fractions[time], delayed_fractions)

    def set_precursor_conc(self, precursor_conc, time='CURRENT'):

        # Check if necessary variables have been set
        cv.check_set(self._num_delayed_groups, 'Material ID={0} precursor_conc'.format(self._id),
                     'Num Delayed Groups')
        cv.check_list_of_floats_or_ints(precursor_conc, 'Material ID={0} precursor conc'.format(self._id))
        cv.check_clock_position(time, 'Material ID={0} precursor conc'.format(self._id))

        if isinstance(precursor_conc, list):
            precursor_conc = np.asarray(precursor_conc)

        np.copyto(self._precursor_conc[time], precursor_conc)

    def set_decay_constant_by_group(self, decay_constant, group, time='CURRENT'):

        # check if decay_constant is a float
        if not cv.is_float(decay_constant):
            msg = 'Unable to set decay constant for Material ID={0} with a ' \
                  'non-float value {1}'.format(self._id, decay_constant)
            raise ValueError(msg)

        # check if group is valid
        elif group < 0 or group > self._num_energy_groups:
            msg = 'Unable to set decay constant for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            self._decay_constants[time][group] = decay_constant

    def set_delayed_fraction_by_group(self, delayed_fraction, group, time='CURRENT'):

        # check if decay_constant is a float
        if not cv.is_float(delayed_fraction):
            msg = 'Unable to set delayed fraction for Material ID={0} with a ' \
                  'non-float value {1}'.format(self._id, delayed_fraction)
            raise ValueError(msg)

        # check if group is valid
        elif group < 0 or group > self._num_energy_groups:
            msg = 'Unable to set delayed fraction for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            self._delayed_fractions[time][group] = delayed_fraction

    def set_precursor_conc_by_group(self, precursor_conc, group, time='CURRENT'):

        # check if decay_constant is a float
        if not cv.is_float(precursor_conc):
            msg = 'Unable to set precursor concfor Material ID={0} with a ' \
                  'non-float value {1}'.format(self._id, precursor_conc)
            raise ValueError(msg)

        # check if group is valid
        elif group < 0 or group > self._num_energy_groups:
            msg = 'Unable to set precursor conc for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:
            self._precursor_conc[time][group] = precursor_conc

    def clone(self):

        new_material = TransientMaterial(name=self._name, clock=self._clock)
        new_material.set_num_energy_groups(self._num_energy_groups)
        new_material.set_num_delayed_groups(self._num_delayed_groups)
        new_material.set_energy_per_fission(self._energy_per_fission)
        for i in self._clock.get_positions():
            new_material.set_sigma_a(self._sigma_a[i], i)
            new_material.set_sigma_t(self._sigma_t[i], i)
            new_material.set_sigma_f(self._sigma_f[i], i)
            new_material.set_nu_sigma_f(self._nu_sigma_f[i], i)
            new_material.set_sigma_s(self._sigma_s[i], i)
            new_material.set_dif_coef(self._dif_coef[i], i)
            new_material.set_chi(self._chi[i], i)
            new_material.set_velocity(self._velocity[i], i)
            new_material.set_decay_constants(self._decay_constants[i], i)
            new_material.set_delayed_fractions(self._delayed_fractions[i], i)
            new_material.set_precursor_conc(self._precursor_conc[i], i)

        return new_material

    def __repr__(self):

        string = 'OpenRK TransientMaterial\n'
        string += ' Name \t\t\t= {0} \n'.format(self._name)
        string += ' ID \t\t\t= {0} \n'.format(self._id)
        string += ' Num Energy Groups \t= {0} \n'.format(self._num_energy_groups)
        string += ' Num Delayed Groups \t= {0} \n'.format(self._num_delayed_groups)
        string += ' Sigma T \t= {0} \n'.format(self._sigma_t)
        string += ' Sigma A \t= {0} \n'.format(self._sigma_a)
        string += ' Sigma F \t= {0} \n'.format(self._sigma_f)
        string += ' Nu Sigma F \t= {0} \n'.format(self._nu_sigma_f)
        string += ' Dif Coef \t= {0} \n'.format(self._dif_coef)
        string += ' Chi \t\t= {0} \n'.format(self._chi)
        string += ' Sigma s \t= {0} \n'.format(self._sigma_s)
        string += ' Energy Per Fission \t= {0} \n'.format(self._energy_per_fission)
        string += ' Velocity \t\t= {0} \n'.format(self._velocity)
        for i in self._clock.get_positions():
            string += ' Time ({:^12s}) \t'.format(i)
            string += ' Precusur Conc = {0} \n'.format(self._precursor_conc[i])
            string += ' Delayed Fractions = {0} \n'.format(self._delayed_fractions[i])
            string += ' Decay Constants = {0} \n'.format(self._decay_constants[i])

        return string


class FunctionalMaterial(TransientMaterial):
    def __init__(self, material_id=None, name='', clock=None):

        super(FunctionalMaterial, self).__init__(material_id, name)

        # Initialize class attributes
        self._num_time_steps = None
        self._time_steps = None
        self._doppler_coefficients = None

    def set_doppler_coefficients(self, doppler_coefs):

        # Check if necessary variables have been set
        cv.check_set(self._num_energy_groups, 'Doppler Coefs', 'Num Energy Groups')
        cv.check_list_of_floats_or_ints(doppler_coefs, 'Doppler Coefs')

        # check if doppler_coefs is a list
        if not cv.is_list(doppler_coefs):
            msg = 'Unable to set absorption xs for Material ID={0} with a ' \
                  'non-list value {1}'.format(self._id, doppler_coefs)
            raise ValueError(msg)

        else:
            if isinstance(doppler_coefs, list):
                self._doppler_coefficients = np.asarray(doppler_coefs)
            else:
                self._doppler_coefficients = np.copy(doppler_coefs)

    def set_sigma_a(self, sigma_a):

        # Check if necessary variables have been set
        cv.check_set(self._num_energy_groups, 'Sigma A', 'Num Energy Groups')
        cv.check_set(self._num_time_steps, 'Sigma A', 'Num Time Steps')
        for i in sigma_a:
            cv.check_list_of_floats_or_ints(i, 'Material ID={0} sigma_a'.format(self._id))

        # check if sigma_a is of length num_time_steps x num_energy_groups
        if np.shape(sigma_a) != (self._num_time_steps, self._num_energy_groups):
            msg = 'Unable to set absorption xs for Material ID={0} with {1} groups ' \
                  'and {2} time steps. Num groups is {3} and num time steps is {4}.' \
                .format(self._id, np.shape(sigma_a)[0], np.shape(sigma_a)[1],
                        self._num_energy_groups, self._num_time_steps)
            raise ValueError(msg)

        else:
            if isinstance(sigma_a, list):
                self._sigma_a = np.asarray(sigma_a)
            else:
                self._sigma_a = np.copy(sigma_a)

    def set_sigma_t(self, sigma_t):

        # Check if necessary variables have been set
        cv.check_set(self._num_energy_groups, 'Sigma T', 'Num Energy Groups')
        cv.check_set(self._num_time_steps, 'Sigma T', 'Num Time Steps')
        for i in sigma_t:
            cv.check_list_of_floats_or_ints(i, 'Material ID={0} sigma_t'.format(self._id))

        # check if sigma_t is of length num_time_steps x num_energy_groups
        if np.shape(sigma_t) != (self._num_time_steps, self._num_energy_groups):
            msg = 'Unable to set total xs for Material ID={0} with {1} groups ' \
                  'and {2} time steps. Num groups is {3} and num time steps is {4}.' \
                .format(self._id, np.shape(sigma_t)[0], np.shape(sigma_t)[1],
                        self._num_energy_groups, self._num_time_steps)
            raise ValueError(msg)

        else:
            if isinstance(sigma_t, list):
                self._sigma_t = np.asarray(sigma_t)
            else:
                self._sigma_t = np.copy(sigma_t)

    def set_sigma_f(self, sigma_f):

        # Check if necessary variables have been set
        cv.check_set(self._num_energy_groups, 'Sigma F', 'Num Energy Groups')
        cv.check_set(self._num_time_steps, 'Sigma F', 'Num Time Steps')
        for i in sigma_f:
            cv.check_list_of_floats_or_ints(i, 'Material ID={0} sigma_f'.format(self._id))

        # check if sigma_f is of length num_time_steps x num_energy_groups
        if np.shape(sigma_f) != (self._num_time_steps, self._num_energy_groups):
            msg = 'Unable to set fission xs for Material ID={0} with {1} groups ' \
                  'and {2} time steps. Num groups is {3} and num time steps is {4}.' \
                .format(self._id, np.shape(sigma_f)[0], np.shape(sigma_f)[1],
                        self._num_energy_groups, self._num_time_steps)
            raise ValueError(msg)

        else:
            if isinstance(sigma_f, list):
                self._sigma_f = np.asarray(sigma_f)
            else:
                self._sigma_f = np.copy(sigma_f)

    def set_nu_sigma_f(self, nu_sigma_f):

        # Check if necessary variables have been set
        cv.check_set(self._num_energy_groups, 'Nu Sigma F', 'Num Energy Groups')
        cv.check_set(self._num_time_steps, 'Nu Sigma F', 'Num Time Steps')
        for i in nu_sigma_f:
            cv.check_list_of_floats_or_ints(i, 'Material ID={0} nu_sigma_f'.format(self._id))

        # check if nu_sigma_f is of length num_time_steps x num_energy_groups
        if np.shape(nu_sigma_f) != (self._num_time_steps, self._num_energy_groups):
            msg = 'Unable to set nu fission xs for Material ID={0} with {1} groups ' \
                  'and {2} time steps. Num groups is {3} and num time steps is {4}.' \
                .format(self._id, np.shape(nu_sigma_f)[0], np.shape(nu_sigma_f)[1],
                        self._num_energy_groups, self._num_time_steps)
            raise ValueError(msg)

        else:
            if isinstance(nu_sigma_f, list):
                self._nu_sigma_f = np.asarray(nu_sigma_f)
            else:
                self._nu_sigma_f = np.copy(nu_sigma_f)

    def set_sigma_s(self, sigma_s):

        # Check if necessary variables have been set
        cv.check_set(self._num_energy_groups, 'Sigma S', 'Num Energy Groups')
        cv.check_set(self._num_time_steps, 'Sigma S', 'Num Time Steps')
        for i in sigma_s:
            cv.check_list_of_floats_or_ints(i, 'Material ID={0} sigma_s'.format(self._id))

        # check if sigma_s is of length num_time_steps x num_energy_groups x num_energy_groups
        if np.shape(sigma_s) != (self._num_time_steps, self._num_energy_groups * self._num_energy_groups):
            msg = 'Unable to set scattering xs for Material ID={0} with {1} groups ' \
                  'and {3} time steps. Num groups is {3} and num time steps is {4}.' \
                .format(self._id, np.shape(sigma_s)[0], math.sqrt(np.shape(sigma_s)[1]),
                        self._num_energy_groups, self._num_time_steps)
            raise ValueError(msg)

        else:
            if isinstance(sigma_s, list):
                self._sigma_s = np.asarray(sigma_s)
            else:
                self._sigma_s = np.copy(sigma_s)

    def set_dif_coef(self, dif_coef):

        # Check if necessary variables have been set
        cv.check_set(self._num_energy_groups, 'Dif Coefs', 'Num Energy Groups')
        cv.check_set(self._num_time_steps, 'Dif Coefs', 'Num Time Steps')
        for i in dif_coef:
            cv.check_list_of_floats_or_ints(i, 'Material ID={0} dif_coef'.format(self._id))

        # check if dif_coef is of length num_time_steps x num_energy_groups
        if np.shape(dif_coef) != (self._num_time_steps, self._num_energy_groups):
            msg = 'Unable to set dif coef for Material ID={0} with {1} groups ' \
                  'and {2} time steps. Num groups is {3} and num time steps is {4}.' \
                .format(self._id, np.shape(dif_coef)[0], np.shape(dif_coef)[1],
                        self._num_energy_groups, self._num_time_steps)
            raise ValueError(msg)

        else:
            if isinstance(dif_coef, list):
                self._dif_coef = np.asarray(dif_coef)
            else:
                self._dif_coef = np.copy(dif_coef)

    def set_chi(self, chi):

        # Check if necessary variables have been set
        cv.check_set(self._num_energy_groups, 'Chi', 'Num Energy Groups')
        cv.check_set(self._num_time_steps, 'Chi', 'Num Time Steps')
        for i in chi:
            cv.check_list_of_floats_or_ints(i, 'Material ID={0} chi'.format(self._id))

        # check if chi is of length num_time_steps x num_energy_groups
        if np.shape(chi) != (self._num_time_steps, self._num_energy_groups):
            msg = 'Unable to set chi for Material ID={0} with {1} groups ' \
                  'and {2} time steps. Num groups is {3} and num time steps is {4}.' \
                .format(self._id, np.shape(chi)[0], np.shape(chi)[1],
                        self._num_energy_groups, self._num_time_steps)
            raise ValueError(msg)

        else:
            if isinstance(chi, list):
                self._chi = np.asarray(chi)
            else:
                self._chi = np.copy(chi)

    def get_sigma_a_by_group(self, group, time=None, temp=None):

        # Check if necessary variables have been set
        cv.check_set(time, 'FunctionalMaterial ID={0} sigma_a by group'.format(self._id), 'time')
        cv.check_set(temp, 'FunctionalMaterial ID={0} sigma_a by group'.format(self._id), 'temp')
        cv.check_set(self._num_energy_groups, 'Material ID={0} sigma_a by group'.format(self._id), 'Num Energy Groups')
        cv.check_set(self._doppler_coefficients, 'Material ID={0} sigma_a by group'.format(self._id),
                     'Doppler Coefficients')
        cv.check_is_int(group, 'Material ID={0} sigma_a by group'.format(self._id), 'group')
        cv.check_is_float_or_int(time, 'Material ID={0} sigma_a by group'.format(self._id), 'time')
        cv.check_is_float_or_int(group, 'Material ID={0} sigma_a by group'.format(self._id), 'temp')

        # check if group is valid
        if group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to get absorption xs for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        # check if time_step is valid
        elif time < self._time_steps[0] or time > self._time_steps[-1]:
            msg = 'Unable to get absorption xs for Material ID={0} for group {1} and' \
                  ' time {2} outside time window [{3}, {4}]' \
                .format(self._id, group, time, self._time_steps[0], self._time_steps[-1])
            raise ValueError(msg)

        else:

            # find time step window
            (i_left, i_right) = self.get_bounding_time_steps(time)

            sigma_a_left = self._sigma_a[i_left][group] * self._doppler_coefficients[group] * \
                (math.sqrt(temp) - math.sqrt(300.0))
            sigma_a_right = self._sigma_a[i_right][group] * self._doppler_coefficients[group] * \
                (math.sqrt(temp) - math.sqrt(300.0))

            # compute sigma_a
            sigma_a = sigma_a_left + (time - self._time_steps[i_left]) * \
                                     (sigma_a_right - sigma_a_left) / \
                                     (self._time_steps[i_right] - self._time_steps[i_left])

            return sigma_a

    def get_sigma_t_by_group(self, group, time=None, temp=None):

        # Check if necessary variables have been set
        cv.check_set(time, 'FunctionalMaterial ID={0} sigma_t by group'.format(self._id), 'time')
        cv.check_set(temp, 'FunctionalMaterial ID={0} sigma_t by group'.format(self._id), 'temp')
        cv.check_set(self._num_energy_groups, 'Material ID={0} sigma_t by group'.format(self._id), 'Num Energy Groups')
        cv.check_is_int(group, 'Material ID={0} sigma_t by group'.format(self._id), 'group')
        cv.check_is_float_or_int(time, 'Material ID={0} sigma_t by group'.format(self._id), 'time')
        cv.check_is_float_or_int(group, 'Material ID={0} sigma_t by group'.format(self._id), 'temp')

        # check if group is valid
        if group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to get total xs for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        # check if time_step is valid
        elif time < self._time_steps[0] or time > self._time_steps[-1]:
            msg = 'Unable to get total xs for Material ID={0} for group {1} and ' \
                  'time {2} outside time window [{3}, {4}]' \
                .format(self._id, group, time, self._time_steps[0], self._time_steps[-1])
            raise ValueError(msg)

        else:

            sigma_t = self.get_sigma_a_by_group(group, time, temp)

            for g in xrange(self._num_energy_groups):
                sigma_t += self.get_sigma_s_by_group(group, g, time, temp)

            return sigma_t

    def get_sigma_f_by_group(self, group, time=None, temp=None):

        # Check if necessary variables have been set
        cv.check_set(time, 'FunctionalMaterial ID={0} sigma_f by group'.format(self._id), 'time')
        cv.check_set(temp, 'FunctionalMaterial ID={0} sigma_f by group'.format(self._id), 'temp')
        cv.check_set(self._num_energy_groups, 'Material ID={0} sigma_f by group'.format(self._id), 'Num Energy Groups')
        cv.check_is_int(group, 'Material ID={0} sigma_f by group'.format(self._id), 'group')
        cv.check_is_float_or_int(time, 'Material ID={0} sigma_f by group'.format(self._id), 'time')
        cv.check_is_float_or_int(group, 'Material ID={0} sigma_f by group'.format(self._id), 'temp')

        # check if group is valid
        if group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to get fission xs for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        # check if time_step is valid
        elif time < self._time_steps[0] or time > self._time_steps[-1]:
            msg = 'Unable to get fission xs for Material ID={0} for group {1} and ' \
                  'time {2} outside time window [{3}, {4}]' \
                .format(self._id, group, time, self._time_steps[0], self._time_steps[-1])
            raise ValueError(msg)

        else:

            # find time step window
            (i_left, i_right) = self.get_bounding_time_steps(time)

            # compute sigma_f
            sigma_f = self._sigma_f[i_left][group] + (time - self._time_steps[i_left]) * \
                                                     (self._sigma_f[i_left][group] - self._sigma_f[i_right][group]) / \
                                                     (self._time_steps[i_right] - self._time_steps[i_left])

            return sigma_f

    def get_nu_sigma_f_by_group(self, group, time=None, temp=None):

        # Check if necessary variables have been set
        cv.check_set(time, 'FunctionalMaterial ID={0} nu_sigma_f by group'.format(self._id), 'time')
        cv.check_set(temp, 'FunctionalMaterial ID={0} nu_sigma_f by group'.format(self._id), 'temp')
        cv.check_set(self._num_energy_groups, 'Material ID={0} nu_sigma_f by group'.format(self._id),
                     'Num Energy Groups')
        cv.check_is_int(group, 'Material ID={0} nu_sigma_f by group'.format(self._id), 'group')
        cv.check_is_float_or_int(time, 'Material ID={0} nu_sigma_f by group'.format(self._id), 'time')
        cv.check_is_float_or_int(group, 'Material ID={0} nu_sigma_f by group'.format(self._id), 'temp')

        # check if group is valid
        if group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to get nu fission xs for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        # check if time_step is valid
        elif time < self._time_steps[0] or time > self._time_steps[-1]:
            msg = 'Unable to get nu fission xs for Material ID={0} for group {1} and ' \
                  'time {2} outside time window [{3}, {4}]' \
                .format(self._id, group, time, self._time_steps[0], self._time_steps[-1])
            raise ValueError(msg)

        else:

            # find time step window
            (i_left, i_right) = self.get_bounding_time_steps(time)

            # compute nu_sigma_f
            nu_sigma_f = self._nu_sigma_f[i_left][group] + (time - self._time_steps[i_left]) * \
                                                           (self._nu_sigma_f[i_left][group] - self._nu_sigma_f[i_right][
                                                               group]) / \
                                                           (self._time_steps[i_right] - self._time_steps[i_left])

            return nu_sigma_f

    def get_sigma_s_by_group(self, group_from, group_to, time=None, temp=None):

        # Check if necessary variables have been set
        cv.check_set(time, 'FunctionalMaterial ID={0} sigma_s by group'.format(self._id), 'time')
        cv.check_set(temp, 'FunctionalMaterial ID={0} sigma_s by group'.format(self._id), 'temp')
        cv.check_set(self._num_energy_groups, 'Material ID={0} sigma_s by group'.format(self._id), 'Num Energy Groups')
        cv.check_is_int(group_from, 'Material ID={0} sigma_s by group'.format(self._id), 'group_from')
        cv.check_is_int(group_to, 'Material ID={0} sigma_s by group'.format(self._id), 'group_to')
        cv.check_is_float_or_int(time, 'Material ID={0} sigma_s by group'.format(self._id), 'time')
        cv.check_is_float_or_int(temp, 'Material ID={0} sigma_s by group'.format(self._id), 'temp')

        # check if group is valid
        if group_from < 0 or group_from > self._num_energy_groups - 1 or \
                group_to < 0 or group_to > self._num_energy_groups - 1:
            msg = 'Unable to get scattering xs for Material ID={0} for group {1} ' \
                  'to {2} as num energy groups is set to {3}' \
                .format(self._id, group_to, group_from, self._num_energy_groups)
            raise ValueError(msg)

        # check if time_step is valid
        elif time < self._time_steps[0] or time > self._time_steps[-1]:
            msg = 'Unable to get scattering xs for Material ID={0} for group {1} to {2} and ' \
                  'time {3} outside time window [{4}, {5}]' \
                .format(self._id, group_from, group_to, time, self._time_steps[0], self._time_steps[-1])
            raise ValueError(msg)

        else:

            ng = self._num_energy_groups

            # find time step window
            (i_left, i_right) = self.get_bounding_time_steps(time)

            # compute sigma_s
            sigma_s = self._sigma_s[i_left][group_from * ng + group_to] + (time - self._time_steps[i_left]) * \
                                                                          (self._sigma_s[i_left][
                                                                           group_from * ng + group_to] -
                                                                           self._sigma_s[i_right][
                                                                           group_from * ng + group_to]) / \
                                                                          (self._time_steps[i_right] - self._time_steps[
                                                                           i_left])

            return sigma_s

    def get_dif_coef_by_group(self, group, time=None, temp=None):

        # Check if necessary variables have been set
        cv.check_set(time, 'FunctionalMaterial ID={0} dif_coef by group'.format(self._id), 'time')
        cv.check_set(temp, 'FunctionalMaterial ID={0} dif_coef by group'.format(self._id), 'temp')
        cv.check_set(self._num_energy_groups, 'Material ID={0} dif_coef by group'.format(self._id), 'Num Energy Groups')
        cv.check_is_int(group, 'Material ID={0} dif_coef by group'.format(self._id), 'group')
        cv.check_is_float_or_int(time, 'Material ID={0} dif_coef by group'.format(self._id), 'time')
        cv.check_is_float_or_int(group, 'Material ID={0} dif_coef by group'.format(self._id), 'temp')

        # check if group is valid
        if group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to get dif coef for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        # check if time_step is valid
        elif time < self._time_steps[0] or time > self._time_steps[-1]:
            msg = 'Unable to get dif coef for Material ID={0} for group {1} and' \
                  ' time {2} outside time window [{3}, {4}]' \
                .format(self._id, group, time, self._time_steps[0], self._time_steps[-1])
            raise ValueError(msg)

        else:

            # find time step window
            (i_left, i_right) = self.get_bounding_time_steps(time)

            # compute dif_coef
            dif_coef = self._dif_coef[i_left][group] + (time - self._time_steps[i_left]) * \
                                                       (self._dif_coef[i_left][group] - self._dif_coef[i_right][group]) / \
                                                       (self._time_steps[i_right] - self._time_steps[i_left])

            return dif_coef

    def get_chi_by_group(self, group, time=None, temp=None):

        # Check if necessary variables have been set
        cv.check_set(time, 'FunctionalMaterial ID={0} chi by group'.format(self._id), 'time')
        cv.check_set(temp, 'FunctionalMaterial ID={0} chi by group'.format(self._id), 'temp')
        cv.check_set(self._num_energy_groups, 'Material ID={0} chi by group'.format(self._id), 'Num Energy Groups')
        cv.check_is_int(group, 'Material ID={0} chi by group'.format(self._id), 'group')
        cv.check_is_float_or_int(time, 'Material ID={0} chi by group'.format(self._id), 'time')
        cv.check_is_float_or_int(group, 'Material ID={0} chi by group'.format(self._id), 'temp')

        # check if group is valid
        if group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to get chi for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        # check if time_step is valid
        elif time < self._time_steps[0] or time > self._time_steps[-1]:
            msg = 'Unable to get chi for Material ID={0} for group {1} and' \
                  ' time {2} outside time window [{3}, {4}]' \
                .format(self._id, group, time, self._time_steps[0], self._time_steps[-1])
            raise ValueError(msg)

        else:

            # find time step window
            (i_left, i_right) = self.get_bounding_time_steps(time)

            # compute chi
            chi = self._chi[i_left][group] + (time - self._time_steps[i_left]) * \
                                             (self._chi[i_left][group] - self._chi[i_right][group]) / \
                                             (self._time_steps[i_right] - self._time_steps[i_left])

            return chi

    def get_doppler_coefficient_by_group(self, group):

        # Check if necessary variables have been set
        cv.check_set(self._num_energy_groups, 'Doppler Coef', 'Num Energy Groups')
        cv.check_is_int(group, 'Material ID={0} doppler coefficient by group'.format(self._id), 'group')

        # check if group is valid
        if group < 0 or group > self._num_energy_groups - 1:
            msg = 'Unable to get doppler coefficient for Material ID={0} for group {1} ' \
                  'as num energy groups is set to {2}' \
                .format(self._id, group, self._num_energy_groups)
            raise ValueError(msg)

        else:

            return self._doppler_coefficients[group]

    def set_num_time_steps(self, num_time_steps):

        # Check if necessary variables have been set
        cv.check_is_int(num_time_steps, 'Material ID={0} num_time_steps'.format(self._id), 'num time steps')

        if num_time_steps < 2:
            msg = 'Unable to set num time steps to {0} since it must be an ' \
                  ' integer > 1'.format(num_time_steps)
            raise ValueError(msg)

        else:
            self._num_time_steps = num_time_steps
            self._time_steps = np.zeros(num_time_steps)

    def set_num_energy_groups(self, num_energy_groups):

        # Check if necessary variables have been set
        cv.check_set(self._num_time_steps, 'Material ID={0} num_energy_groups'.format(self._id), 'Num Time Steps')

        if cv.is_integer(num_energy_groups):

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
                self._sigma_s = np.zeros((self._num_time_steps, num_energy_groups * num_energy_groups))
                self._dif_coef = np.zeros((self._num_time_steps, num_energy_groups))
                self._chi = np.zeros((self._num_time_steps, num_energy_groups))
                self._doppler_coefficients = np.zeros(num_energy_groups)
                self._velocity = np.zeros(num_energy_groups)

        else:
            msg = 'Unable to set num energy groups to non-integer {0}' \
                .format(num_energy_groups)
            raise ValueError(msg)

    def set_time_steps(self, time_steps):

        # Check if necessary variables have been set
        cv.check_set(self._num_time_steps, 'Material ID={0} time_steps'.format(self._id), 'Num Time Steps')
        cv.check_list_of_floats_or_ints(time_steps, 'Material ID={0} time_steps'.format(self._id))

        if isinstance(time_steps, list):
            self._time_steps = np.asarray(time_steps)
        else:
            self._time_steps = np.copy(time_steps)

    def get_bounding_time_steps(self, time):

        # check if time_step is valid
        if time < self._time_steps[0] or time > self._time_steps[-1]:
            msg = 'Unable to get bounding time steps for Material ID={0} for' \
                  ' time {1} outside time window [{2}, {3}]' \
                .format(self._id, self._num_time_steps, self._time_steps[0], self._time_steps[-1])
            raise ValueError(msg)

        else:

            i_left = 0
            while time > self._time_steps[i_left + 1]:
                i_left += 1

            i_right = i_left + 1
            return i_left, i_right

    def __repr__(self):

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
        for i in self._clock.get_positions():
            string += ' Precursor Conc ({:^12s}) \t'.format(i)
            string += '= {0} \n'.format(self._precursor_conc[i])

        return string