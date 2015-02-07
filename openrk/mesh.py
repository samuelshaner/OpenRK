__author__ = 'Samuel Shaner'
__email__ = 'shaner@mit.edu'

# Import modules
import numpy as np
import checkvalue as cv
from clock import Clock
from math import floor, exp, cos, pi, asin
from material import TransientMaterial, FunctionalMaterial

# A static variable for auto-generated Mesh UIDs
AUTO_MESH_UID = 1


class Mesh(object):
    """
    Main Mesh class which contains information about the geometry and
    contains field variables for various material properties

    Attributes:
      _materials               = array of Material objects
      _dif_linear              = map of linear surface diffusion coefficients
      _dif_nonlinear           = map of nonlinear surface diffusion coefficients
      _x_min                   = minimum x-value
      _y_min                   = minimum y-value
      _x_max                   = maximum x-value
      _y_max                   = maximum y-value
      _offset                  = offset of origin from center of geometry
      _num_cells               = number of cells in mesh
      _boundaries              = boundary condition for the four sides (1 - vacuum; 0 - reflective)
      _num_shape_energy_groups = number of energy groups in shape solve
      _num_amp_energy_groups   = number of energy groups in amp solve
      _num_delayed_groups      = number of delayed groups
      _clock                   = clock indicating times at various states
      _flux                    = map of flux at various states
      _temperature             = map of temperature at various states
      _power                   = map of power at various states
      _current                 = map of current at various states
      _k_eff_0                 = the initial eigenvalue
      _buckling                = the axial buckling
      _decay_constants         = the delayed neutron precursor decay constants
      _delayed_fractions       = the delayed neutron precursor fractional neutron yields

    Setter Methods:
      set_clock(Clock clock)
      set_name(str name)
      set_id(int mesh_id)
      set_x_min(float x_min)
      set_y_min(float y_min)
      set_x_max(float x_max)
      set_y_max(float y_max)
      set_boundary(int side, int boundary)
      set_num_shape_energy_groups(int num_groups)
      set_num_amp_energy_groups(int num_groups)
      set_num_delayed_groups(int num_groups)

    Getter Methods:
      get_width()
      get_height()
      get_boundary(int side)
      get_bounds()
      get_material(int cell)
      get_num_shape_energy_groups()
      get_num_amp_energy_groups()
      get_num_delayed_groups()
      get_flux(str time='CURRENT')
      get_flux_by_value(int cell, int group, str time='CURRENT')
      get_temperature(str time='CURRENT')
      get_temperature_by_value(int cell, str time='CURRENT')
      get_current(str time='CURRENT')
      get_current_by_value(int cell, int group, int side, str time='CURRENT')
      get_power(str time='CURRENT')
      get_power_by_value(int cell, str time='CURRENT')

    Other Methods:
      intialize_field_variables()
      compute_power(str time='CURRENT')
      copy_flux(str time_from, str time_to)
      copy_temperature(str time_from, str time_to)
      copy_current(str time_from, str time_to)
      copy_power(str time_from, str time_to)
      uniform_refine(int num_refines)
      uniquify_materials()

    """
    def __init__(self, mesh_id=None, name='', width=1.0, height=1.0):

        # Initialize class attributes
        global AUTO_MESH_UID
        self._uid = AUTO_MESH_UID
        AUTO_MESH_UID += 1
        self._id = None
        self._set_id = False
        self._name = ''
        self._k_eff_0 = None

        # Initialize class attributes
        self._materials = None
        self._x_min = None
        self._x_max = None
        self._y_min = None
        self._y_max = None
        self._offset = np.zeros(2)
        self._num_cells = None
        self._boundaries = np.zeros(4, dtype=np.int)
        self._num_shape_energy_groups = None
        self._num_amp_energy_groups = None
        self._num_delayed_groups = None
        self._clock = None
        self._flux = {}
        self._temperature = {}
        self._power = {}
        self._fuel_volume = None
        self._buckling = None
        self._decay_constants = None
        self._delayed_fractions = None

        # Set Mesh properties
        self.set_x_min(-width / 2.0)
        self.set_x_max(width / 2.0)
        self.set_y_min(-height / 2.0)
        self.set_y_max(height / 2.0)
        self.set_name(name)

        # Initialize flux and current arrays
        clock = Clock()
        for position in clock.get_positions():
            self._flux[position] = np.empty
            self._temperature[position] = np.empty
            self._power[position] = np.empty

        # Set mesh ID
        if mesh_id is not None:
            self.set_id(mesh_id)
        else:
            self.set_id(self._uid)
            self._set_id = False

    def set_k_eff_0(self, k_eff_0):

        cv.check_is_float(k_eff_0, 'Material initial k_eff', 'k_eff_0')
        self._k_eff_0 = k_eff_0

    def get_k_eff_0(self):

        cv.check_set(self._k_eff_0, 'initial k_eff', 'k_eff_0')
        return self._k_eff_0

    def set_clock(self, clock):

        # Initialize clock
        if not isinstance(clock, Clock):
            msg = 'Unable to initialize Mesh clock since clock input is not of type ' \
                  'Clock: {0}'.format(clock)
            raise ValueError(msg)
        else:
            self._clock = clock

            if self._materials is not None:
                for mat in self._materials:
                    if mat is not None and isinstance(mat, TransientMaterial):
                        mat.set_clock(clock)

    def get_clock(self):

        return self._clock

    def set_name(self, name):

        if not cv.is_string(name):
            msg = 'Unable to set name for Mesh ID={0} with a non-string ' \
                  'value {1}'.format(self._id, name)
            raise ValueError(msg)

        else:
            self._name = name

    def set_id(self, mesh_id):

        # Check that the ID is a non-negative integer
        if cv.is_integer(mesh_id):

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

    def set_x_min(self, x_min):

        if not cv.is_float(x_min):
            msg = 'Unable to set x_min to non-integer {0}' \
                .format(x_min)
            raise ValueError(msg)

        else:
            self._x_min = x_min

    def set_x_max(self, x_max):

        if not cv.is_float(x_max):
            msg = 'Unable to set x_max to non-integer {0}' \
                .format(x_max)
            raise ValueError(msg)

        else:
            self._x_max = x_max

    def set_y_min(self, y_min):

        if not cv.is_float(y_min):
            msg = 'Unable to set y_min to non-integer {0}' \
                .format(y_min)
            raise ValueError(msg)

        else:
            self._y_min = y_min

    def set_y_max(self, y_max):

        if not cv.is_float(y_max):
            msg = 'Unable to set y_max to non-integer {0}' \
                .format(y_max)
            raise ValueError(msg)

        else:
            self._y_max = y_max

    def get_width(self):
        return self._x_max - self._x_min

    def get_height(self):
        return self._y_max - self._y_min

    def set_boundary(self, side, boundary):

        if not cv.is_integer(side):
            msg = 'Unable to set boundary for non-integer side {0}' \
                .format(side)
            raise ValueError(msg)

        elif not cv.is_integer(boundary):
            msg = 'Unable to set boundary for non-integer boundary {0}' \
                .format(boundary)
            raise ValueError(msg)

        elif side < 0 or side > 4:
            msg = 'Unable to set boundary for invalid side {0}' \
                .format(side)
            raise ValueError(msg)

        elif boundary < 0 or boundary > 1:
            msg = 'Unable to set boundary for invalid boundary {0}' \
                .format(boundary)
            raise ValueError(msg)

        else:
            self._boundaries[side] = boundary

    def get_boundary(self, side):

        if not cv.is_integer(side):
            msg = 'Unable to get boundary for non-integer side {0}' \
                .format(side)
            raise ValueError(msg)

        elif side < 0 or side > 4:
            msg = 'Unable to set boundary for invalid side {0}' \
                .format(side)
            raise ValueError(msg)

        else:
            return self._boundaries[side]

    def get_bounds(self):

        return [self._x_min, self._x_max, self._y_min, self._y_max]

    def set_num_shape_energy_groups(self, num_groups):

        if not cv.is_integer(num_groups):
            msg = 'Unable to set num shape energy groups for non-integer {0}' \
                .format(num_groups)
            raise ValueError(msg)

        elif num_groups < 1:
            msg = 'Unable to set num shape energy groups for non-positive {0}' \
                .format(num_groups)
            raise ValueError(msg)

        else:
            self._num_shape_energy_groups = num_groups

    def get_num_shape_energy_groups(self):

        cv.check_set(self._num_shape_energy_groups, 'Mesh get num shape energy groups', 'self._num_shape_energy_groups')
        return self._num_shape_energy_groups

    def set_num_amp_energy_groups(self, num_groups):

        if not cv.is_integer(num_groups):
            msg = 'Unable to set num amp energy groups for non-integer {0}' \
                .format(num_groups)
            raise ValueError(msg)

        elif num_groups < 1:
            msg = 'Unable to set num amp energy groups for non-positive {0}' \
                .format(num_groups)
            raise ValueError(msg)

        else:
            self._num_amp_energy_groups = num_groups

    def get_num_amp_energy_groups(self):

        cv.check_set(self._num_amp_energy_groups, 'Mesh get num amp energy groups', 'self._num_amp_energy_groups')
        return self._num_amp_energy_groups

    def set_num_delayed_groups(self, num_groups):

        if not cv.is_integer(num_groups):
            msg = 'Unable to set num delayed groups for non-integer {0}' \
                .format(num_groups)
            raise ValueError(msg)

        elif num_groups < 1:
            msg = 'Unable to set num delayed groups for non-positive {0}' \
                .format(num_groups)
            raise ValueError(msg)

        else:
            self._num_delayed_groups = num_groups

    def get_num_delayed_groups(self):

        cv.check_set(self._num_delayed_groups, 'Mesh get num delayed groups', 'self._num_delayed_groups')
        return self._num_delayed_groups

    def get_flux(self, time='CURRENT'):

        cv.check_clock_position(time, 'Mesh get flux')
        return self._flux[time]

    def get_temperature(self, time='CURRENT'):

        cv.check_clock_position(time, 'Mesh get temperature')
        return self._temperature[time]

    def get_power(self, time='CURRENT'):

        cv.check_clock_position(time, 'Mesh get power')
        return self._power[time]

    def get_temperature_by_value(self, cell, time='CURRENT'):

        cv.check_clock_position(time, 'AmpMesh temperature by value')
        return self._temperature[time][cell]

    def get_power_by_value(self, cell, time='CURRENT'):

        cv.check_clock_position(time, 'AmpMesh power by value')
        return self._power[time][cell]

    def copy_flux(self, time_from, time_to):

        np.copyto(self._flux[time_to], self._flux[time_from])

    def copy_temperature(self, time_from, time_to):

        np.copyto(self._temperature[time_to], self._temperature[time_from])

    def copy_power(self, time_from, time_to):

        np.copyto(self._power[time_to], self._power[time_from])

    def compute_flux_l2_norm(self, time_1='CURRENT', time_2='FORWARD_IN_OLD'):

        return np.linalg.norm((self._flux[time_1] - self._flux[time_2]) / self._flux[time_1])

    def set_temperature(self, temperature, time='CURRENT'):

        if cv.is_list(temperature):
            if len(temperature) != self._num_x * self._num_y:
                msg = 'unable to set Mesh temperature field since input temperature field '\
                    'is of length {0} and there are {1} cells'.format(len(temperature),
                                                                      self._num_x * self._num_y)
                raise ValueError(msg)

            else:
                if isinstance(temperature, list):
                    temperature = np.asarray(temperature)

                np.copyto(self._temperature[time], temperature)

        else:
            cv.check_is_float_or_int(temperature, 'set Mesh temperature', 'temperature')
            self._temperature[time].fill(temperature)

    def set_temperature_by_value(self, temp, cell, time='CURRENT'):

        cv.check_clock_position(time, 'AmpMesh temp by value')
        self._temperature[time][cell] = temp

    def interpolate_flux(self, time_begin='PREVIOUS_OUT', time_end='FORWARD_OUT', time='CURRENT'):

        dt = self._clock.get_time(time_end) - self._clock.get_time(time_begin)
        wt_begin = (self._clock.get_time(time_end) - self._clock.get_time(time)) / dt
        wt_end = (self._clock.get_time(time) - self._clock.get_time(time_begin)) / dt
        self.copy_flux(time_begin, time)
        self._flux[time] *= wt_begin
        self._flux[time] += wt_end * self._flux[time_end]

    def set_material(self, material, cell_id):

        self._materials[cell_id] = material

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
        self._cell_width = None
        self._cell_height = None
        self._current = {}
        self._dif_linear = {}
        self._dif_nonlinear = {}

        # Set Mesh properties
        self.set_num_x(num_x)
        self.set_num_y(num_y)

    def set_decay_constants(self, decay_constants):

        # check if decay constants is a list
        if not cv.is_list(decay_constants):
            msg = 'Unable to set decay constants for StructuredMesh with a ' \
                  'non-list value {0}'.format(decay_constants)
            raise ValueError(msg)

        # check if decay constants is of length num_groups
        elif len(decay_constants) != self._num_delayed_groups:
            msg = 'Unable to set decay constants for StructuredMesh with {0} groups ' \
                  'as num delayed groups is set to {1}' \
                .format(len(decay_constants), self._num_delayed_groups)
            raise ValueError(msg)

        else:
            if isinstance(decay_constants, list):
                decay_constants = np.asarray(decay_constants)

            if self._decay_constants is None:
                self._decay_constants = np.zeros(self._num_delayed_groups)

            np.copyto(self._decay_constants, decay_constants)

    def set_delayed_fractions(self, delayed_fractions):

        # check if decay constants is a list
        if not cv.is_list(delayed_fractions):
            msg = 'Unable to set delayed fractions for StructuredMesh with a ' \
                  'non-list value {0}'.format(delayed_fractions)
            raise ValueError(msg)

        # check if delayed fractions is of length num_groups
        elif len(delayed_fractions) != self._num_delayed_groups:
            msg = 'Unable to set delayed fractions for StructuredMesh with {0} groups ' \
                  'as num delayed groups is set to {1}' \
                .format(len(delayed_fractions), self._num_delayed_groups)
            raise ValueError(msg)

        else:
            if isinstance(delayed_fractions, list):
                delayed_fractions = np.asarray(delayed_fractions)

            if self._delayed_fractions is None:
                self._delayed_fractions = np.zeros(self._num_delayed_groups)

            np.copyto(self._delayed_fractions, delayed_fractions)

    def get_buckling(self):

        return self._buckling

    def get_decay_constants(self):

        return self._decay_constants

    def get_delayed_fractions(self):

        return self._delayed_fractions

    def get_decay_constant_by_group(self, group):

        # check if group is valid
        if group < 0 or group > self._num_delayed_groups - 1:
            msg = 'Unable to get decay constant for StructuredMesh for group {0} ' \
                  'as num delayed groups is set to {1}' \
                .format(group, self._num_delayed_groups)
            raise ValueError(msg)

        else:
            return self._decay_constants[group]

    def get_delayed_fraction_by_group(self, group):

        # check if group is valid
        if group < 0 or group > self._num_delayed_groups - 1:
            msg = 'Unable to get delayed fraction for StructuredMesh for group {0} ' \
                  'as num delayed groups is set to {1}' \
                .format(group, self._num_delayed_groups)
            raise ValueError(msg)

        else:
            return self._delayed_fractions[group]

    def compute_fuel_volume(self):

        self._fuel_volume = 0.0
        for material in self._materials:
            if material.get_is_fissionable():
                self._fuel_volume += 1

        self._fuel_volume *= self.get_cell_volume()

    def get_current(self, time='CURRENT'):

        cv.check_clock_position(time, 'StructuredMesh get current')
        return self._current[time]

    def get_dif_linear(self, time='CURRENT'):

        cv.check_clock_position(time, 'StructuredMesh get dif_linear')
        return self._dif_linear[time]

    def get_dif_nonlinear(self, time='CURRENT'):

        cv.check_clock_position(time, 'StructuredMesh get dif_nonlinear')
        return self._dif_nonlinear[time]

    def copy_current(self, time_from, time_to):

        np.copyto(self._current[time_to], self._current[time_from])

    def copy_dif_linear(self, time_from, time_to):

        np.copyto(self._dif_linear[time_to], self._dif_linear[time_from])

    def copy_dif_nonlinear(self, time_from, time_to):

        np.copyto(self._dif_nonlinear[time_to], self._dif_nonlinear[time_from])

    def interpolate_dif_nonlinear(self, time_begin='PREVIOUS_OUT', time_end='FORWARD_OUT', time='CURRENT'):

        dt = self._clock.get_time(time_end) - self._clock.get_time(time_begin)
        wt_begin = (self._clock.get_time(time_end) - self._clock.get_time(time)) / dt
        wt_end = (self._clock.get_time(time) - self._clock.get_time(time_begin)) / dt
        self.copy_dif_nonlinear(time_begin, time)
        self._dif_nonlinear[time] *= wt_begin
        self._dif_nonlinear[time] += wt_end * self._dif_nonlinear[time_end]

    def scale_flux(self, scale_val, time='CURRENT'):

        self._flux[time] = self._flux[time] * scale_val

    def scale_current(self, scale_val, time='CURRENT'):

        self._current[time] = self._current[time] * scale_val

    def get_average_power(self, time='CURRENT'):

        if self._fuel_volume is None:
            self.compute_fuel_volume()

        self.compute_power(time)

        return np.sum(self._power[time]) * self.get_cell_volume() / self._fuel_volume

    def get_neighbor_cell(self, x, y, side):

        # Check input values
        cv.check_is_int(x, 'Structured Mesh {0} get neighbor cell'.format(self._name), 'x')
        cv.check_is_int(y, 'Structured Mesh {0} get neighbor cell'.format(self._name), 'y')
        cv.check_is_int(side, 'Structured Mesh {0} get neighbor cell'.format(self._name), 'side')

        if side < 0 or side > 3:
            msg = 'Cannot get neighbor cell from Structured Mesh with invalid side {0}.'.format(side)
            raise AttributeError(msg)

        else:

            neighbor_cell = None

            if side == 0:
                if x != 0:
                    neighbor_cell = y*self._num_x + x - 1
            elif side == 1:
                if y != 0:
                    neighbor_cell = (y - 1)*self._num_x + x
            elif side == 2:
                if x != self._num_x - 1:
                    neighbor_cell = y*self._num_x + x + 1
            elif side == 3:
                if y != self._num_y - 1:
                    neighbor_cell = (y + 1) * self._num_x + x

            return neighbor_cell

    def get_neighbor_material(self, x, y, side):

        neighbor_cell = self.get_neighbor_cell(x, y, side)

        if neighbor_cell is None:
            return None
        else:
            return self._materials[neighbor_cell]

    def get_num_x(self):

        return self._num_x

    def get_num_y(self):

        return self._num_y

    def get_cell_width(self):

        if self._cell_width is None:
            msg = 'Cannot get cell width of Structured Mesh as num cells x has not been set!'
            raise AttributeError(msg)

        else:

            return self._cell_width

    def get_cell_height(self):

        if self._cell_height is None:
            msg = 'Cannot get cell height of Structured Mesh as num cells y has not been set!'
            raise AttributeError(msg)

        else:

            return self._cell_height

    def get_cell_volume(self):

        return self.get_cell_width() * self.get_cell_height()

    def set_num_x(self, num_x):

        # Check input values
        cv.check_is_int(num_x, 'Structured mesh num x', 'num x')

        if num_x < 1:
            msg = 'Unable to set num cells x to non-positive {0}' \
                .format(num_x)
            raise ValueError(msg)

        else:

            self._num_x = num_x
            self._cell_width = (self._x_max - self._x_min) / num_x

    def set_num_y(self, num_y):

        # Check input values
        cv.check_is_int(num_y, 'Structured mesh num y', 'num y')

        if num_y < 1:
            msg = 'Unable to set num cells x to non-positive {0}' \
                .format(num_y)
            raise ValueError(msg)

        else:
            self._num_y = num_y
            self._cell_height = (self._y_max - self._y_min) / num_y

    def find_cell(self, x, y):

        # Check input values
        cv.check_is_float_or_int(x, 'Structured mesh find cell', 'x')
        cv.check_is_float_or_int(y, 'Structured mesh find cell', 'y')

        i = floor((x - self._x_min) / self._cell_width)
        j = floor((y - self._y_min) / self._cell_height)

        return int(j * self._num_x + i)

    def get_material(self, cell_id):

        return self._materials[cell_id]

    def uniform_refine(self, num_refines):

        cv.check_is_int(num_refines, 'Structured mesh number of refinements', 'num refines')

        # Clone the current mesh
        mesh = self.clone()

        # reset the num_x and num_y
        mesh.set_num_x(self._num_x * num_refines)
        mesh.set_num_y(self._num_y * num_refines)
        mesh.initialize()

        for j in xrange(self._num_y * num_refines):
            jj = j / num_refines
            for i in xrange(self._num_x * num_refines):
                ii = i / num_refines
                mesh.set_material(self.get_material(jj*self._num_x + ii), j*self._num_x*num_refines + i)
                mesh.set_temperature_by_value(self._temperature['CURRENT'][jj*self._num_x + ii],
                                              j*self._num_x*num_refines + i)

        return mesh

    def uniquify_materials(self):

        for i in xrange(self._num_x * self._num_y):
            new_material = self.get_material(i).clone()
            self._materials[i] = new_material

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

        self._shape_map = None
        self._shape_mesh = None
        self._shape_groups = None
        self._optically_thick = False

        # Initialize flux and current arrays
        clock = Clock()
        for position in clock.get_positions():
            self._current[position] = np.empty

    def set_optically_thick(self, optically_thick):

        self._optically_thick = optically_thick

    def set_buckling(self, buckling):

        # check if buckling is a list
        if not cv.is_list(buckling):
            msg = 'Unable to set buckling for AmpMesh with a ' \
                  'non-list value {0}'.format(buckling)
            raise ValueError(msg)

        # check if buckling is of length num_groups
        elif len(buckling) != self._num_amp_energy_groups:
            msg = 'Unable to set buckling for AmpMesh with {0} groups ' \
                  'as num amp energy groups is set to {1}' \
                .format(len(buckling), self._num_amp_energy_groups)
            raise ValueError(msg)

        else:
            if isinstance(buckling, list):
                buckling = np.asarray(buckling)

            if self._buckling is None:
                self._buckling = np.zeros(self._num_amp_energy_groups)

            np.copyto(self._buckling, buckling)

    def get_buckling_by_group(self, group):

        # check if group is valid
        if group < 0 or group > self._num_amp_energy_groups - 1:
            msg = 'Unable to get buckling for AmpMesh for group {0} ' \
                  'as num amp energy groups is set to {1}' \
                .format(group, self._num_amp_energy_groups)
            raise ValueError(msg)

        else:
            return self._buckling[group]

    def initialize(self):

        self._materials = np.empty(shape=[self._num_y*self._num_x], dtype=object)

        # create material objects
        for y in range(self._num_y):
            for x in range(self._num_x):
                material = TransientMaterial()
                material.set_num_energy_groups(self._num_amp_energy_groups)
                material.set_num_delayed_groups(self._num_delayed_groups)
                self._materials[y*self._num_x + x] = material

        # Create map of field variables and diffusion coefficients
        clock = Clock()
        for position in clock.get_positions():
            self._dif_linear[position] = np.zeros(self._num_x * self._num_y * self._num_amp_energy_groups * 4)
            self._dif_nonlinear[position] = np.zeros(self._num_x * self._num_y * self._num_amp_energy_groups * 4)
            self._current[position] = np.zeros(self._num_x * self._num_y * self._num_amp_energy_groups * 4)
            self._flux[position] = np.ones(self._num_x * self._num_y * self._num_amp_energy_groups)
            self._temperature[position] = np.ones(self._num_x * self._num_y) * 400
            self._power[position] = np.zeros(self._num_x * self._num_y)

        if self._buckling is None:
            self._buckling = np.zeros(self._num_amp_energy_groups)

        if self._delayed_fractions is None or self._decay_constants is None:
            msg = 'Cannot initialize AmpMesh without setting the delayed fractions and decay constants'
            raise ValueError(msg)

    def set_shape_groups(self, shape_groups):

        self._shape_groups = shape_groups

    def set_shape_mesh(self, shape_mesh):

        if not isinstance(shape_mesh, StructuredShapeMesh):
            msg = 'Cannot generate AmpMesh to Shape Mesh map with a {0}'\
                ' type Mesh'.format(type(shape_mesh))
            raise ValueError(msg)

        num_refines = shape_mesh.get_num_x() / self._num_x
        self._shape_mesh = shape_mesh

        # Create list of lists to store mesh map
        self._shape_map = []
        for j in xrange(self._num_y):
            for i in xrange(self._num_x):
                self._shape_map.append([])

        # Populate shape mesh map
        for j in xrange(self._num_y * num_refines):
            jj = j / num_refines
            for i in xrange(self._num_x * num_refines):
                ii = i / num_refines
                amp_cell = jj*self._num_x+ii
                shape_cell = j*self._num_y*num_refines + i
                self._shape_map[amp_cell].append(shape_cell)

    def condense_materials(self, time='CURRENT', save_flux=False):

        shape_cell_volume = self._shape_mesh.get_cell_width() * self._shape_mesh.get_cell_height()
        amp_cell_volume = self.get_cell_width() * self.get_cell_height()
        temps = self._shape_mesh.get_temperature(time)

        for i in xrange(self._num_y*self._num_x):

            # amp material
            amp_mat = self.get_material(i)

            for g in xrange(self._num_amp_energy_groups):

                # initialize variables
                sigma_a = 0.0
                sigma_t = 0.0
                sigma_f = 0.0
                nu_sigma_f = 0.0
                sigma_s = np.zeros(self._num_amp_energy_groups)
                chi = np.zeros(self._num_amp_energy_groups)
                dif_coef = 0.0
                rxn = 0.0
                production = 0.0
                velocity = 0.0

                # Condense chi
                for ii in self._shape_map[i]:

                    shape_mat = self._shape_mesh.get_material(ii)
                    temp = temps[ii]

                    for e in xrange(self._num_amp_energy_groups):
                        chi_tally = 0.0

                        for ee in self._shape_groups[e]:
                            chi_tally += shape_mat.get_chi_by_group(ee, time, temp)

                        for h in xrange(self._num_shape_energy_groups):
                            chi[e] += chi_tally * shape_mat.get_nu_sigma_f_by_group(h, time, temp) \
                                * self._shape_mesh.get_flux_by_value(ii, h, time) * shape_cell_volume
                            production += chi_tally * shape_mat.get_nu_sigma_f_by_group(h, time, temp) \
                                * self._shape_mesh.get_flux_by_value(ii, h, time) * shape_cell_volume

                for gg in self._shape_groups[g]:

                    rxn_group = 0.0
                    volume = 0.0
                    sigma_t_group = 0.0

                    for ii in self._shape_map[i]:

                        shape_mat = self._shape_mesh.get_material(ii)
                        temp = temps[ii]

                        # accumulate cross section tallies
                        shape = self._shape_mesh.get_flux_by_value(ii, gg, time)
                        sigma_a += shape_mat.get_sigma_a_by_group(gg, time, temp) * shape * shape_cell_volume
                        sigma_t += shape_mat.get_sigma_t_by_group(gg, time, temp) * shape * shape_cell_volume
                        sigma_f += shape_mat.get_sigma_f_by_group(gg, time, temp) * shape * shape_cell_volume
                        nu_sigma_f += shape_mat.get_nu_sigma_f_by_group(gg, time, temp) * shape * shape_cell_volume
                        rxn += shape * shape_cell_volume
                        rxn_group += shape * shape_cell_volume
                        volume += shape_cell_volume
                        sigma_t_group += shape_mat.get_sigma_t_by_group(gg, time, temp) * shape * shape_cell_volume
                        velocity += 1.0 / shape_mat.get_velocity_by_group(gg, time, temp) * shape * shape_cell_volume

                        for h in xrange(self._num_shape_energy_groups):
                            sigma_s[self._shape_mesh.get_amp_group(h)] += shape_mat.get_sigma_s_by_group(gg, h, time, temp) * \
                                shape * shape_cell_volume

                    dif_coef += rxn_group / (3.0 * sigma_t_group / rxn_group)

                # set amp mesh xs and flux
                amp_mat.set_sigma_a_by_group(sigma_a / rxn, g, time)
                amp_mat.set_sigma_t_by_group(sigma_t / rxn, g, time)
                amp_mat.set_sigma_f_by_group(sigma_f / rxn, g, time)
                amp_mat.set_nu_sigma_f_by_group(nu_sigma_f / rxn, g, time)
                amp_mat.set_dif_coef_by_group(dif_coef / rxn, g, time)
                amp_mat.set_velocity_by_group(1.0 / (velocity / rxn), g, time)

                if save_flux:
                    self._flux[time][i*self._num_amp_energy_groups+g] = rxn / volume

                # set chi
                if production != 0.0:
                    amp_mat.set_chi_by_group(chi[g] / production, g, time)
                else:
                    amp_mat.set_chi_by_group(0.0, g, time)

                # set scattering xs
                for e in xrange(self._num_amp_energy_groups):
                    amp_mat.set_sigma_s_by_group(sigma_s[e] / rxn, g, e, time)

            # Condense delayed neutron precursors
            for d in xrange(self._num_delayed_groups):
                precursor_conc = 0.0

                for ii in self._shape_map[i]:
                    if isinstance(self._shape_mesh.get_material(ii), TransientMaterial):
                        precursor_conc += self._shape_mesh.get_material(ii).get_precursor_conc_by_group(d, time) \
                            * shape_cell_volume

                amp_mat.set_precursor_conc_by_group(precursor_conc / amp_cell_volume, d, time)

    def get_flux_by_value(self, cell, group, time='CURRENT'):

        cv.check_clock_position(time, 'AmpMesh flux by value')
        return self._flux[time][cell * self._num_amp_energy_groups + group]

    def get_current_by_value(self, cell, group, surface, time='CURRENT'):

        cv.check_clock_position(time, 'AmpMesh current by value')
        return self._current[time][(cell*4 + surface) * self._num_amp_energy_groups + group]

    def get_dif_linear_by_value(self, cell, group, surface, time='CURRENT'):

        cv.check_clock_position(time, 'AmpMesh dif linear by value')
        return self._dif_linear[time][(cell*4 + surface) * self._num_amp_energy_groups + group]

    def get_dif_nonlinear_by_value(self, cell, group, surface, time='CURRENT'):

        cv.check_clock_position(time, 'AmpMesh dif nonlinear by value')
        return self._dif_nonlinear[time][(cell*4 + surface) * self._num_amp_energy_groups + group]

    def set_flux_by_value(self, flux, cell, group, time='CURRENT'):

        cv.check_clock_position(time, 'AmpMesh flux by value')
        self._flux[time][cell * self._num_amp_energy_groups + group] = flux

    def set_current_by_value(self, current, cell, group, surface, time='CURRENT'):

        cv.check_clock_position(time, 'AmpMesh current by value')
        self._current[time][(cell*4 + surface) * self._num_amp_energy_groups + group] = current

    def set_dif_linear_by_value(self, dif_linear, cell, group, surface, time='CURRENT'):

        cv.check_clock_position(time, 'AmpMesh dif linear by value')
        self._dif_linear[time][(cell*4 + surface) * self._num_amp_energy_groups + group] = dif_linear

    def set_dif_nonlinear_by_value(self, dif_nonlinear, cell, group, surface, time='CURRENT'):

        cv.check_clock_position(time, 'AmpMesh dif nonlinear by value')
        self._dif_nonlinear[time][(cell*4 + surface) * self._num_amp_energy_groups + group] = dif_nonlinear

    def compute_power(self, time='CURRENT'):

        cv.check_clock_position(time, 'AmpMesh compute power')

        for y in xrange(self._num_y):
            for x in xrange(self._num_x):

                # Get cell properties
                cell_id = y*self._num_x + x
                material = self._materials[cell_id]
                fission_rate = 0.0

                # Compute the fission rate
                for g in xrange(self._num_amp_energy_groups):
                    fission_rate += material.get_sigma_f_by_group(g) * self.get_flux_by_value(cell_id, g, time)

                self._power[time][cell_id] = fission_rate * material.get_energy_per_fission()

    def clone(self):

        mesh = AmpMesh(name=self._name, width=self.get_width(), height=self.get_height(),
                       num_x=self._num_x, num_y=self._num_y)

        mesh.set_num_amp_energy_groups(self._num_amp_energy_groups)
        mesh.set_num_shape_energy_groups(self._num_shape_energy_groups)
        mesh.set_num_delayed_groups(self._num_delayed_groups)

        if self._buckling is not None:
            mesh.set_buckling(self._buckling)

        if self._delayed_fractions is not None:
            mesh.set_delayed_fractions(self._delayed_fractions)

        if self._decay_constants is not None:
            mesh.set_decay_constants(self._decay_constants)

        for s in xrange(4):
            mesh.set_boundary(s, self.get_boundary(s))

        return mesh

    def compute_current(self, time='CURRENT'):

        num_refines = self._shape_mesh.get_num_x() / self._num_x
        sm = self._shape_mesh
        sm_nx = sm.get_num_x()
        sm_cw = sm.get_cell_width()
        sm_ch = sm.get_cell_height()
        temps = sm.get_temperature(time)

        for i in xrange(self._num_y*self._num_x):

            for g in xrange(self._num_amp_energy_groups):

                for s in xrange(4):

                    current = 0.0

                    for ii in self._shape_map[i]:

                        mat = sm.get_material(ii)
                        temp = temps[ii]
                        cell_next = sm.get_neighbor_cell(ii % sm_nx, ii / sm_nx, s)
                        temp_next = temps[cell_next]
                        mat_next = sm.get_neighbor_material(ii % sm_nx, ii / sm_nx, s)

                        if s == 0 and (ii % sm_nx) % num_refines == 0:

                            for gg in self._shape_groups[g]:
                                flux = sm.get_flux_by_value(ii, gg, time)
                                d = mat.get_dif_coef_by_group(gg, time, temp)
                                length_perpen = sm_cw
                                length = sm_ch

                                if mat_next is None:
                                    dif_linear = 2 * d / length_perpen / (1 + 4 * d / length_perpen)
                                    dif_linear *= sm.get_boundary(s)
                                    current += - dif_linear * flux * length

                                else:
                                    d_next = mat_next.get_dif_coef_by_group(gg, time, temp_next)
                                    dif_linear = 2 * d * d_next / (length_perpen * d + length_perpen * d_next)
                                    flux_next = sm.get_flux_by_value(cell_next, gg, time)
                                    current += - dif_linear * (flux - flux_next) * length

                        elif s == 1 and (ii / sm_nx) % num_refines == 0:

                            for gg in self._shape_groups[g]:
                                flux = sm.get_flux_by_value(ii, gg, time)
                                d = mat.get_dif_coef_by_group(gg, time, temp)
                                length_perpen = sm_ch
                                length = sm_cw

                                if mat_next is None:
                                    dif_linear = 2 * d / length_perpen / (1 + 4 * d / length_perpen)
                                    dif_linear *= sm.get_boundary(s)
                                    current += - dif_linear * flux * length

                                else:
                                    d_next = mat_next.get_dif_coef_by_group(gg, time, temp_next)
                                    dif_linear = 2 * d * d_next / (length_perpen * d + length_perpen * d_next)
                                    flux_next = sm.get_flux_by_value(cell_next, gg, time)
                                    current += - dif_linear * (flux - flux_next) * length

                        elif s == 2 and (ii % sm_nx) % num_refines == num_refines - 1:

                            for gg in self._shape_groups[g]:
                                flux = sm.get_flux_by_value(ii, gg, time)
                                d = mat.get_dif_coef_by_group(gg, time, temp)
                                length_perpen = sm_cw
                                length = sm_ch

                                if mat_next is None:
                                    dif_linear = 2 * d / length_perpen / (1 + 4 * d / length_perpen)
                                    dif_linear *= sm.get_boundary(s)
                                    current += dif_linear * flux * length

                                else:
                                    d_next = mat_next.get_dif_coef_by_group(gg, time, temp_next)
                                    dif_linear = 2 * d * d_next / (length_perpen * d + length_perpen * d_next)
                                    flux_next = sm.get_flux_by_value(cell_next, gg, time)
                                    current += - dif_linear * (flux_next - flux) * length

                        elif s == 3 and (ii / sm_nx) % num_refines == num_refines - 1:

                            for gg in self._shape_groups[g]:
                                flux = sm.get_flux_by_value(ii, gg, time)
                                d = mat.get_dif_coef_by_group(gg, time, temp)
                                length_perpen = sm_ch
                                length = sm_cw

                                if mat_next is None:
                                    dif_linear = 2 * d / length_perpen / (1 + 4 * d / length_perpen)
                                    dif_linear *= sm.get_boundary(s)
                                    current += dif_linear * flux * length

                                else:
                                    d_next = mat_next.get_dif_coef_by_group(gg, time, temp_next)
                                    dif_linear = 2 * d * d_next / (length_perpen * d + length_perpen * d_next)
                                    flux_next = sm.get_flux_by_value(cell_next, gg, time)
                                    current += - dif_linear * (flux_next - flux) * length

                    # set the current for this surface
                    self._current[time][(i*4 + s) * self._num_amp_energy_groups + g] = current

    def compute_dif_correct(self, dif_coef, length):

        if self._optically_thick:

            mu = cos(asin(0.798184))
            expon = exp(- length / (3 * dif_coef * mu))
            alpha = (1 + expon) / (1 - expon) - 2 * 3 * dif_coef * mu / length
            rho = mu * alpha
            f = 1 + length * rho / (2 * dif_coef)

        else:
            f = 1.0

        return f

    def compute_dif_coefs(self, time='CURRENT'):

        nx = self.get_num_x()
        ny = self.get_num_y()
        ng = self.get_num_amp_energy_groups()
        width = self.get_cell_width()
        height = self.get_cell_height()
        temps = self.get_temperature(time)

        for y in xrange(ny):
            for x in xrange(nx):

                cell = y*nx+x
                temp = temps[cell]

                for s in xrange(4):

                    cell_next = self.get_neighbor_cell(x, y, s)

                    if s == 0 or s == 1:
                        sense = -1.0
                    else:
                        sense = 1.0

                    if s == 0 or s == 2:
                        length = height
                        length_perpen = width
                    else:
                        length = width
                        length_perpen = height

                    for g in xrange(ng):

                        dif_coef = self.get_material(cell).get_dif_coef_by_group(g, time, temp)
                        current = self.get_current_by_value(cell, g, s, time)
                        flux = self.get_flux_by_value(cell, g, time)
                        f = self.compute_dif_correct(dif_coef, length_perpen)

                        if cell_next is None:

                            dif_linear = 2 * dif_coef * f / length_perpen / (1.0 + 4.0 * dif_coef * f / length_perpen)
                            dif_nonlinear = (sense * dif_linear * flux - current / length) / flux
                            dif_linear *= self.get_boundary(s)
                            dif_nonlinear *= self.get_boundary(s)
                        else:

                            flux_next = self.get_flux_by_value(cell_next, g, time)
                            dif_coef_next = self.get_material(cell_next).\
                                get_dif_coef_by_group(g, time, temps[cell_next])
                            f_next = self.compute_dif_correct(dif_coef_next, length_perpen)
                            dif_linear = 2.0 * dif_coef * f * dif_coef_next * f_next / \
                                (length_perpen * dif_coef * f + length_perpen * dif_coef_next * f_next)
                            dif_nonlinear = - (sense * dif_linear * (flux_next - flux) +
                                               current / length) / (flux_next + flux)

                        if dif_nonlinear > dif_linear:
                            if sense == -1.0:
                                if dif_nonlinear > 0.0:
                                    dif_linear = - current / (2 * flux)
                                    dif_nonlinear = - current / (2 * flux)
                                else:
                                    dif_linear = current / (2 * flux_next)
                                    dif_nonlinear = - current / (2 * flux_next)
                            else:
                                if dif_nonlinear > 0.0:
                                    dif_linear = - current / (2 * flux_next)
                                    dif_nonlinear = - current / (2 * flux_next)
                                else:
                                    dif_linear = current / (2 * flux)
                                    dif_nonlinear = - current / (2 * flux)

                        self.set_dif_linear_by_value(dif_linear, cell, g, s, time)
                        self.set_dif_nonlinear_by_value(dif_nonlinear, cell, g, s, time)

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


class StructuredShapeMesh(StructuredMesh):
    def __init__(self, mesh_id=None, name='', width=1.0, height=1.0, num_x=1, num_y=1):

        # initialize FunctionalMaterial class attributes
        super(StructuredShapeMesh, self).__init__(mesh_id, name, width, height, num_x, num_y)

        self._amp_map = None
        self._amp_mesh = None
        self._amp_groups = None

    def set_buckling(self, buckling):

        # check if buckling is a list
        if not cv.is_list(buckling):
            msg = 'Unable to set buckling for StructuredShapeMesh with a ' \
                  'non-list value {0}'.format(buckling)
            raise ValueError(msg)

        # check if buckling is of length num_groups
        elif len(buckling) != self._num_shape_energy_groups:
            msg = 'Unable to set buckling for StructuredShapeMesh with {0} groups ' \
                  'as num shape energy groups is set to {1}' \
                .format(len(buckling), self._num_shape_energy_groups)
            raise ValueError(msg)

        else:
            if isinstance(buckling, list):
                buckling = np.asarray(buckling)

            if self._buckling is None:
                self._buckling = np.zeros(self._num_shape_energy_groups)

            np.copyto(self._buckling, buckling)

    def get_buckling_by_group(self, group):

        # check if group is valid
        if group < 0 or group > self._num_shape_energy_groups - 1:
            msg = 'Unable to get buckling for AmpMesh for group {0} ' \
                  'as num shape energy groups is set to {1}' \
                .format(group, self._num_shape_energy_groups)
            raise ValueError(msg)

        else:
            return self._buckling[group]

    def initialize(self):

        self._materials = np.empty(shape=[self._num_y*self._num_x], dtype=object)

        # Create map of field variables and diffusion coefficients
        clock = Clock()
        for position in clock.get_positions():
            self._dif_linear[position] = np.zeros(self._num_x * self._num_y * self._num_shape_energy_groups * 4)
            self._dif_nonlinear[position] = np.zeros(self._num_x * self._num_y * self._num_shape_energy_groups * 4)
            self._current[position] = np.zeros(self._num_x * self._num_y * self._num_shape_energy_groups * 4)
            self._flux[position] = np.ones(self._num_x * self._num_y * self._num_shape_energy_groups)
            self._temperature[position] = np.ones(self._num_x * self._num_y) * 400
            self._power[position] = np.zeros(self._num_x * self._num_y)

        if self._buckling is None:
            self._buckling = np.zeros(self._num_shape_energy_groups)

        if self._delayed_fractions is None or self._decay_constants is None:
            msg = 'Cannot initialize StructuredShapeMesh without setting the delayed fractions and decay constants'
            raise ValueError(msg)

    def synthesize_flux(self, time='CURRENT'):

        ngs = self._num_shape_energy_groups
        nga = self._num_amp_energy_groups
        shape_previous = np.zeros(self._num_x * self._num_y * ngs)
        shape_forward = np.zeros(self._num_x * self._num_y * ngs)
        shape_current = np.zeros(self._num_x * self._num_y * ngs)

        # Compute PREVIOUS_OUT and FORWARD_OUT shapes
        for i in xrange(self._num_x*self._num_y):
            for g in xrange(ngs):
                shape_previous[i*ngs+g] = self._flux['PREVIOUS_OUT'][i*ngs+g] / \
                    self._amp_mesh.get_flux('PREVIOUS_OUT')[self._amp_map[i]*nga+self._amp_groups[g]]
                shape_forward[i*ngs+g] = self._flux['FORWARD_OUT'][i*ngs+g] / \
                    self._amp_mesh.get_flux('FORWARD_OUT')[self._amp_map[i]*nga+self._amp_groups[g]]

        # Interpolate to get shape at time
        dt = self._clock.get_time('FORWARD_OUT') - self._clock.get_time('PREVIOUS_OUT')
        wt_begin = (self._clock.get_time('FORWARD_OUT') - self._clock.get_time(time)) / dt
        wt_end = (self._clock.get_time(time) - self._clock.get_time('PREVIOUS_OUT')) / dt
        np.copyto(shape_current, shape_previous)
        shape_current *= wt_begin
        shape_current += wt_end * shape_forward

        # Reconstruct flux at time
        for i in xrange(self._num_x*self._num_y):
            for g in xrange(ngs):
                self._flux[time][i*ngs+g] = shape_current[i*ngs+g] * \
                    self._amp_mesh.get_flux(time)[self._amp_map[i]*nga+self._amp_groups[g]]

    def reconstruct_flux(self, time='FORWARD_OUT_OLD', time_shape='PREVIOUS_OUT', time_amp='CURRENT'):

        ngs = self._num_shape_energy_groups
        nga = self._num_amp_energy_groups
        shape = np.zeros(self._num_x * self._num_y * ngs)

        # Compute the shape at time_shape
        for i in xrange(self._num_x*self._num_y):
            for g in xrange(ngs):
                shape[i*ngs+g] = self._flux[time_shape][i*ngs+g] / \
                    self._amp_mesh.get_flux(time_shape)[self._amp_map[i]*nga+self._amp_groups[g]]

        # Reconstruct flux at time using amp at time_amp and shape at time_shape
        for i in xrange(self._num_x*self._num_y):
            for g in xrange(ngs):
                self._flux[time][i*ngs+g] = shape[i*ngs+g] * \
                    self._amp_mesh.get_flux(time_amp)[self._amp_map[i]*nga+self._amp_groups[g]]

    def set_amp_groups(self, amp_groups):

        if isinstance(amp_groups, list):
            amp_groups = np.asarray(amp_groups)

        self._amp_groups = np.copy(amp_groups)

    def get_amp_group(self, group):

        return int(self._amp_groups[group])

    def set_amp_mesh(self, amp_mesh):

        if not isinstance(amp_mesh, AmpMesh):
            msg = 'Cannot set the StructuredShapeMesh amp mesh with non-AmpMesh '\
                'input: {0}'.format(amp_mesh)
            raise ValueError(msg)

        num_refines = self._num_x / amp_mesh.get_num_x()
        self._amp_mesh = amp_mesh
        self._amp_map = np.zeros(self._num_x*self._num_y, dtype=np.int)

        # Populate amp mesh map
        for j in xrange(self._num_y):
            jj = j / num_refines
            for i in xrange(self._num_x):
                ii = i / num_refines
                amp_cell = jj*self._num_x/num_refines+ii
                shape_cell = j*self._num_y + i
                self._amp_map[shape_cell] = amp_cell

    def get_flux_by_value(self, cell, group, time='CURRENT'):

        cv.check_clock_position(time, 'StructuredShapeMesh flux by value')
        return self._flux[time][cell * self._num_shape_energy_groups + group]

    def get_current_by_value(self, cell, group, surface, time='CURRENT'):

        cv.check_clock_position(time, 'StructuredShapeMesh current by value')
        return self._current[time][(cell*4 + surface) * self._num_shape_energy_groups + group]

    def get_dif_linear_by_value(self, cell, group, surface, time='CURRENT'):

        cv.check_clock_position(time, 'StructuredShapeMesh dif linear by value')
        return self._dif_linear[time][(cell*4 + surface) * self._num_shape_energy_groups + group]

    def get_dif_nonlinear_by_value(self, cell, group, surface, time='CURRENT'):

        cv.check_clock_position(time, 'StructuredShapeMesh dif nonlinear by value')
        return self._dif_nonlinear[time][(cell*4 + surface) * self._num_shape_energy_groups + group]

    def set_flux_by_value(self, flux, cell, group, time='CURRENT'):

        cv.check_clock_position(time, 'StructuredShapeMesh flux by value')
        self._flux[time][cell * self._num_shape_energy_groups + group] = flux

    def set_current_by_value(self, current, cell, group, surface, time='CURRENT'):

        cv.check_clock_position(time, 'StructuredShapeMesh current by value')
        self._current[time][(cell*4 + surface) * self._num_shape_energy_groups + group] = current

    def set_dif_linear_by_value(self, dif_linear, cell, group, surface, time='CURRENT'):

        cv.check_clock_position(time, 'StructuredShapeMesh dif linear by value')
        self._dif_linear[time][(cell*4 + surface) * self._num_shape_energy_groups + group] = dif_linear

    def set_dif_nonlinear_by_value(self, dif_nonlinear, cell, group, surface, time='CURRENT'):

        cv.check_clock_position(time, 'StructuredShapeMesh dif nonlinear by value')
        self._dif_nonlinear[time][(cell*4 + surface) * self._num_shape_energy_groups + group] = dif_nonlinear

    def clone(self):

        mesh = StructuredShapeMesh(name=self._name, width=self.get_width(), height=self.get_height(),
                                   num_x=self._num_x, num_y=self._num_y)

        mesh.set_num_amp_energy_groups(self._num_amp_energy_groups)
        mesh.set_num_shape_energy_groups(self._num_shape_energy_groups)
        mesh.set_num_delayed_groups(self._num_delayed_groups)

        if self._buckling is not None:
            mesh.set_buckling(self._buckling)

        if self._delayed_fractions is not None:
            mesh.set_delayed_fractions(self._delayed_fractions)

        if self._decay_constants is not None:
            mesh.set_decay_constants(self._decay_constants)

        for s in xrange(4):
            mesh.set_boundary(s, self.get_boundary(s))

        return mesh

    def compute_power(self, time='CURRENT'):

        cv.check_clock_position(time, 'StructuredShapeMesh compute power')
        temps = self.get_temperature(time)

        for y in xrange(self._num_y):
            for x in xrange(self._num_x):

                # Get cell properties
                cell_id = y*self._num_x + x
                material = self._materials[y*self._num_x + x]
                fission_rate = 0.0
                temp = temps[cell_id]

                if isinstance(material, TransientMaterial):

                    # Compute the fission rate
                    for g in xrange(self._num_shape_energy_groups):
                        fission_rate += material.get_sigma_f_by_group(g, time, temp) * self.get_flux_by_value(cell_id, g, time)

                    self._power[time][cell_id] = fission_rate * material.get_energy_per_fission()

    def compute_initial_precursor_conc(self, time='CURRENT'):

        cv.check_clock_position(time, 'StructuredShapeMesh compute initial precursor conc')
        temps = self._temperature[time]

        for y in xrange(self._num_y):
            for x in xrange(self._num_x):

                # Get cell properties
                cell_id = y*self._num_x + x
                material = self._materials[y*self._num_x + x]
                fission_rate = 0.0
                temp = temps[cell_id]

                if isinstance(material, TransientMaterial):

                    # Compute the fission rate
                    for g in xrange(self._num_shape_energy_groups):
                        fission_rate += material.get_nu_sigma_f_by_group(g, time, temp) * self.get_flux_by_value(cell_id, g, time)

                    for d in xrange(material.get_num_delayed_groups()):
                        conc = fission_rate * self.get_delayed_fraction_by_group(d) \
                            / self.get_k_eff_0() / self.get_decay_constant_by_group(d)
                        material.set_precursor_conc_by_group(conc, d, time)

    def integrate_precursor_conc(self, time_from='CURRENT', time_to='FORWARD_OUT'):

        cv.check_clock_position(time_from, 'StructuredShapeMesh integrate precursor conc')
        cv.check_clock_position(time_to, 'StructuredShapeMesh integrate precursor conc')
        dt = self._clock.get_time(time_to) - self._clock.get_time(time_from)
        temps_from = self._temperature[time_from]
        temps_to = self._temperature[time_to]

        for y in xrange(self._num_y):
            for x in xrange(self._num_x):

                # Get cell properties
                cell_id = y*self._num_x + x
                material = self._materials[y*self._num_x + x]
                fission_rate_from = 0.0
                fission_rate_to = 0.0
                temp_from = temps_from[cell_id]
                temp_to = temps_to[cell_id]

                if isinstance(material, TransientMaterial):

                    # Compute the fission rate
                    for g in xrange(self._num_shape_energy_groups):
                        fission_rate_from += material.get_nu_sigma_f_by_group(g, time_from, temp_from) * \
                            self.get_flux_by_value(cell_id, g, time_from)
                        fission_rate_to += material.get_nu_sigma_f_by_group(g, time_to, temp_to) * \
                            self.get_flux_by_value(cell_id, g, time_to)

                    for d in xrange(material.get_num_delayed_groups()):
                        delayed_fraction = self.get_delayed_fraction_by_group(d)
                        decay_constant = self.get_decay_constant_by_group(d)
                        k1 = exp(- decay_constant * dt)
                        k2 = 1.0 - (1.0 - k1) / (decay_constant * dt)
                        k3 = k1 + (k2 - 1.0)

                        conc = material.get_precursor_conc_by_group(d, time_from)
                        new_conc = k1 * conc + k2 * delayed_fraction / decay_constant / self.get_k_eff_0() \
                            * fission_rate_to - k3 * delayed_fraction / decay_constant / self.get_k_eff_0() \
                            * fission_rate_from

                        material.set_precursor_conc_by_group(new_conc, d, time_to)

    def integrate_temperature(self, time_from='CURRENT', time_to='FORWARD_OUT'):

        cv.check_clock_position(time_from, 'StructuredShapeMesh integrate temperature')
        cv.check_clock_position(time_to, 'StructuredShapeMesh integrate temperature')
        dt = self._clock.get_time(time_to) - self._clock.get_time(time_from)
        temps_from = self._temperature[time_from]
        temps_to = self._temperature[time_to]

        for y in xrange(self._num_y):
            for x in xrange(self._num_x):

                # Get cell properties
                cell_id = y*self._num_x + x
                material = self._materials[y*self._num_x + x]
                fission_rate_from = 0.0
                fission_rate_to = 0.0
                temp_from = temps_from[cell_id]
                temp_to = temps_to[cell_id]

                if isinstance(material, FunctionalMaterial):

                    # Compute the fission rate
                    for g in xrange(self._num_shape_energy_groups):
                        fission_rate_from += material.get_nu_sigma_f_by_group(g, time_from, temp_from) * \
                            self.get_flux_by_value(cell_id, g, time_from)
                        fission_rate_to += material.get_nu_sigma_f_by_group(g, time_to, temp_to) * \
                            self.get_flux_by_value(cell_id, g, time_to)

                    self._temperature[time_to][cell_id] = self._temperature[time_from][cell_id] + dt * 0.5 * \
                        (fission_rate_from + fission_rate_to) * material.get_temperature_conversion_factor()

    def compute_dif_coefs(self, time='CURRENT'):

        nx = self.get_num_x()
        ny = self.get_num_y()
        ng = self.get_num_shape_energy_groups()
        width = self.get_cell_width()
        height = self.get_cell_height()
        temps = self.get_temperature(time)

        for y in xrange(ny):
            for x in xrange(nx):

                mat = self.get_material(y*nx+x)

                for side in xrange(4):

                    cell_next = self.get_neighbor_cell(x, y, side)
                    mat_next = self.get_neighbor_material(x, y, side)

                    for e in xrange(ng):

                        d = mat.get_dif_coef_by_group(e, time, temps[y*nx+x])

                        # set the length of the surface parallel to and perpendicular from surface
                        if side == 0 or side == 2:
                            length_perpen = width
                        elif side == 1 or side == 3:
                            length_perpen = height

                        if mat_next is None:
                            dif_linear = 2 * d / length_perpen / (1 + 4 * d / length_perpen)
                            dif_linear *= self.get_boundary(side)

                        else:

                            if side == 0 or side == 2:
                                next_length_perpen = width
                            else:
                                next_length_perpen = height

                            d_next = mat_next.get_dif_coef_by_group(e, time, temps[cell_next])
                            dif_linear = 2 * d * d_next / (length_perpen * d + next_length_perpen * d_next)

                        self.set_dif_linear_by_value(dif_linear, y*nx+x, e, side, time)

    def __repr__(self):

        string = 'OpenRK StructuredShapeMesh\n'
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

    def set_num_cells(self, num_cells):

        if not cv.is_integer(num_cells):
            msg = 'Unable to set num cells for non-integer {0}' \
                .format(num_cells)
            raise ValueError(msg)

        elif num_cells < 1:
            msg = 'Unable to set num cells for non-positive {0}' \
                .format(num_cells)
            raise ValueError(msg)

        else:
            self._num_cells = num_cells

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
        super(UnstructuredShapeMesh, self).__init__(mesh_id, name, width, height)

        self._volumes = []

    def get_flux_by_value(self, cell, group, time='CURRENT'):

        cv.check_clock_position(time, 'UnstructuredShapeMesh flux by value')
        return self._flux[time][cell * self._num_shape_energy_groups + group]

    def get_temperature_by_value(self, cell, time='CURRENT'):

        cv.check_clock_position(time, 'UnstructuredShapeMesh temperature by value')
        return self._temperature[time][cell]

    def get_power_by_value(self, cell, time='CURRENT'):

        cv.check_clock_position(time, 'UnstructuredShapeMesh power by value')
        return self._power[time][cell]

    def initialize_field_variables(self):

        self._volumes = np.zeros(self._num_x * self._num_y)

        clock = Clock()
        for position in clock.get_positions():
            self._flux[position] = np.ones(self._num_x * self._num_y * self._num_shape_energy_groups)
            self._temperature[position] = np.ones(self._num_x * self._num_y) * 400
            self._power[position] = np.zeros(self._num_x * self._num_y)

    def compute_power(self, time='CURRENT'):

        cv.check_clock_position(time, 'UnstructuredShapeMesh compute power')

        for i in xrange(self._num_cells):

            # Get cell properties
            material = self._materials[i]
            fission_rate = 0.0

            # Compute the fission rate
            for g in xrange(self._num_shape_energy_groups):
                fission_rate += material.get_sigma_f_by_group(g) * self.get_flux_by_value(i, g, time)

            self._power[time][i] = fission_rate * material.get_energy_per_fission()

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