__author__ = 'Samuel Shaner'
__email__ = 'shaner@mit.edu'

# Import modules
import numpy as np
import checkvalue as cv
from cell import Cell, CmfdCell
from clock import Clock
from math import floor
from surface import Surface
from material import TransientMaterial
import openrk

# A static variable for auto-generated Mesh UIDs
AUTO_MESH_UID = 1


class Mesh(object):
    """
    Main Mesh class which contains information about the geometry and
    contains field variables for various material properties

    Attributes:
      _cells                   = array of Cell objects
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
      uniqueify_materials()

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
        self._cells = None
        self._x_min = None
        self._x_max = None
        self._y_min = None
        self._y_max = None
        self._offset = np.zeros(2)
        self._num_cells = None
        self._boundaries = np.zeros(4)
        self._num_shape_energy_groups = None
        self._num_amp_energy_groups = None
        self._num_delayed_groups = None
        self._clock = None
        self._flux = {}
        self._temperature = {}
        self._power = {}

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

        if side < 0 or side > 4:
            msg = 'Unable to set boundary for invalid side {0}' \
                .format(side)
            raise ValueError(msg)

        if boundary < 0 or boundary > 2:
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

        # Set Mesh properties
        self.set_num_x(num_x)
        self.set_num_y(num_y)

    def get_cell(self, x, y):

        # Check input values
        cv.check_is_int(x, 'Structured Mesh {0} get cell'.format(self._name), 'x')
        cv.check_is_int(y, 'Structured Mesh {0} get cell'.format(self._name), 'y')

        if self._cells is None:
            msg = 'Cannot get cell from Structured Mesh as cells have not been initialized!'
            raise AttributeError(msg)

        else:

            return self._cells[y*self._num_x+x]

    def get_neighbor_cell(self, x, y, side):

        # Check input values
        cv.check_is_int(x, 'Structured Mesh {0} get neighbor cell'.format(self._name), 'x')
        cv.check_is_int(y, 'Structured Mesh {0} get neighbor cell'.format(self._name), 'y')
        cv.check_is_int(side, 'Structured Mesh {0} get neighbor cell'.format(self._name), 'side')

        if self._cells is None:
            msg = 'Cannot get neighbor cell from Structured Mesh as materials have not been initialized!'
            raise AttributeError(msg)

        elif side < 0 or side > 3:
            msg = 'Cannot get neighbor cell from Structured Mesh with invalid side {0}.'.format(side)
            raise AttributeError(msg)

        else:

            neighbor_cell = None

            if side == 0:
                if x != 0:
                    neighbor_cell = self._cells[y*self._num_x + x - 1]
            elif side == 1:
                if y != 0:
                    neighbor_cell = self._cells[(y - 1)*self._num_x + x]
            elif side == 2:
                if x != self._num_x - 1:
                    neighbor_cell = self._cells[y*self._num_x + x + 1]
            elif side == 3:
                if y != self._num_y - 1:
                    neighbor_cell = self._cells[(y + 1) * self._num_x + x]

            return neighbor_cell

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

        return self._cells[cell_id].get_material()

    def uniform_refine(self, num_refines):

        cv.check_is_int(num_refines, 'Structured mesh number of refinements', 'num refines')

        # Clone the current mesh
        mesh = self.clone()

        # reset the num_x and num_y
        mesh.set_num_x(self._num_x * num_refines)
        mesh.set_num_y(self._num_y * num_refines)
        mesh.initialize_cells()
        mesh.set_num_amp_energy_groups(self._num_amp_energy_groups)
        mesh.set_num_shape_energy_groups(self._num_shape_energy_groups)

        for j in xrange(self._num_y * num_refines):
            jj = j / num_refines
            for i in xrange(self._num_x * num_refines):
                ii = i / num_refines
                mesh._cells[j*self._num_x*num_refines + i].\
                    set_material(self._cells[jj*self._num_x + ii].get_material())

        return mesh

    def uniquify_materials(self):

        for i in xrange(self._num_x * self._num_y):
            new_material = self.get_material(i).clone()
            self._cells[i].set_material(new_material)

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

        self._num_fsrs = None
        self._current = {}
        self._shape_map = None
        self._shape_mesh = None
        self._shape_groups = None

        # Initialize flux and current arrays
        clock = Clock()
        for position in clock.get_positions():
            self._current[position] = np.empty

    def set_num_fsrs(self, num_fsrs):

        self._num_fsrs = num_fsrs

    def initialize_cells(self):

        self._cells = np.empty(shape=[self._num_y*self._num_x], dtype=object)

        # create cell objects
        for y in range(self._num_y):
            for x in range(self._num_x):
                self._cells[y*self._num_x + x] = CmfdCell()
                material = TransientMaterial()
                material.set_num_energy_groups(self._num_amp_energy_groups)
                material.set_num_delayed_groups(self._num_delayed_groups)
                self._cells[y*self._num_x + x].set_material(material)

        self.initialize_surfaces()
        self.initialize_field_variables()

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

    def condense_materials(self, time='CURRENT'):

        shape_cell_volume = self._shape_mesh.get_cell_width() * self._shape_mesh.get_cell_height()
        amp_cell_volume = self.get_cell_width() * self.get_cell_height()

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

                # Condense chi
                for ii in self._shape_map[i]:

                    shape_mat = self._shape_mesh.get_material(ii)

                    for e in xrange(self._num_amp_energy_groups):
                        chi_tally = 0.0

                        for ee in self._shape_groups[e]:
                            chi_tally += shape_mat.get_chi_by_group(ee, time)

                        for h in xrange(self._num_shape_energy_groups):
                            chi[e] += chi_tally * shape_mat.get_nu_sigma_f_by_group(h, time) \
                                * self._shape_mesh.get_flux_by_value(ii, h, time) * shape_cell_volume
                            production += chi_tally * shape_mat.get_nu_sigma_f_by_group(h, time) \
                                * self._shape_mesh.get_flux_by_value(ii, h, time) * shape_cell_volume

                for gg in self._shape_groups[g]:

                    rxn_group = 0.0
                    volume = 0.0
                    sigma_t_group = 0.0

                    for ii in self._shape_map[i]:

                        shape_mat = self._shape_mesh.get_material(ii)

                        # accumulate cross section tallies
                        shape = self._shape_mesh.get_flux_by_value(ii, gg, time)
                        sigma_a += shape_mat.get_sigma_a_by_group(gg, time) * shape * shape_cell_volume
                        sigma_t += shape_mat.get_sigma_t_by_group(gg, time) * shape * shape_cell_volume
                        sigma_f += shape_mat.get_sigma_f_by_group(gg, time) * shape * shape_cell_volume
                        nu_sigma_f += shape_mat.get_nu_sigma_f_by_group(gg, time) * shape * shape_cell_volume
                        rxn += shape * shape_cell_volume
                        rxn_group += shape * shape_cell_volume
                        volume += shape_cell_volume
                        sigma_t_group += shape_mat.get_sigma_t_by_group(gg, time) * shape * shape_cell_volume

                        for h in xrange(self._num_shape_energy_groups):
                            sigma_s[self._shape_mesh.get_amp_group(h)] += shape_mat.get_sigma_s_by_group(h, gg, time) * \
                                shape * shape_cell_volume

                    dif_coef += rxn_group / (3.0 * sigma_t_group / rxn_group)

                # set amp mesh xs and flux
                amp_mat.set_sigma_a_by_group(sigma_a / rxn, g, time)
                amp_mat.set_sigma_t_by_group(sigma_t / rxn, g, time)
                amp_mat.set_sigma_f_by_group(sigma_f / rxn, g, time)
                amp_mat.set_nu_sigma_f_by_group(nu_sigma_f / rxn, g, time)
                amp_mat.set_dif_coef_by_group(dif_coef / rxn, g, time)
                self._flux[time][i*self._num_amp_energy_groups+g] = rxn / volume

                # set chi
                if production != 0.0:
                    amp_mat.set_chi_by_group(chi[g] / production, g, time)
                else:
                    amp_mat.set_chi_by_group(0.0, g, time)

                # set scattering xs
                for e in xrange(self._num_amp_energy_groups):
                    amp_mat.set_sigma_s_by_group(sigma_s[e], g, e, time)

            # Condense delayed neutron precursors
            for d in xrange(self._num_delayed_groups):
                precursor_conc = 0.0

                for ii in self._shape_map[i]:
                    if isinstance(self._shape_mesh.get_material(ii), TransientMaterial):
                        precursor_conc += self._shape_mesh.get_material(ii).get_precursor_conc_by_group(d, time) \
                            * shape_cell_volume

                self.get_material(i).set_precursor_conc_by_group(precursor_conc / amp_cell_volume, d, time)

    def initialize_surfaces(self):

        for y in xrange(self._num_y):
            for x in xrange(self._num_x):

                cell = self._cells[y*self._num_x + x]

                # Surface 0
                if x == 0:
                    cell._surfaces[0] = Surface(cell.get_material().get_num_energy_groups())
                else:
                    cell._surfaces[0] = self._cells[y*self._num_x + x - 1]._surfaces[2]

                # Surface 1
                if y == 0:
                    cell._surfaces[1] = Surface(cell.get_material().get_num_energy_groups())
                else:
                    cell._surfaces[1] = self._cells[(y - 1)*self._num_x + x]._surfaces[3]

                # Surface 2
                cell._surfaces[2] = Surface(cell.get_material().get_num_energy_groups())

                # Surface 3
                cell._surfaces[3] = Surface(cell.get_material().get_num_energy_groups())

    def get_flux_by_value(self, cell, group, time='CURRENT'):

        cv.check_clock_position(time, 'AmpMesh flux by value')
        return self._flux[time][cell * self._num_amp_energy_groups + group]

    def get_temperature_by_value(self, cell, time='CURRENT'):

        cv.check_clock_position(time, 'AmpMesh temperature by value')
        return self._temperature[time][cell]

    def get_current_by_value(self, cell, group, surface, time='CURRENT'):

        cv.check_clock_position(time, 'AmpMesh current by value')
        return self._current[time][(cell*4 + surface) * self._num_amp_energy_groups + group]

    def get_power_by_value(self, cell, time='CURRENT'):

        cv.check_clock_position(time, 'AmpMesh power by value')
        return self._power[time][cell]

    def get_current(self, time='CURRENT'):

        cv.check_clock_position(time, 'AmpMesh get current')
        return self._current[time]

    def copy_flux(self, time_from, time_to):

        openrk.vector_copy(self.get_flux[time_from], self.get_flux[time_to])

    def copy_temperature(self, time_from, time_to):

        openrk.vector_copy(self.get_temperature[time_from], self.get_temperature[time_to])

    def copy_power(self, time_from, time_to):

        openrk.vector_copy(self.get_power[time_from], self.get_power[time_to])

    def copy_current(self, time_from, time_to):

        openrk.vector_copy(self.get_current[time_from], self.get_current[time_to])

    def compute_power(self, time='CURRENT'):

        cv.check_clock_position(time, 'AmpMesh compute power')
        cell_volume = self.get_cell_width() * self.get_cell_height()

        for y in xrange(self._num_y):
            for x in xrange(self._num_x):

                # Get cell properties
                cell_id = y*self._num_x + x
                material = self._cells[cell_id].get_material()
                fission_rate = 0.0

                # Compute the fission rate
                for g in xrange(self._num_amp_energy_groups):
                    fission_rate += material.get_sigma_f_by_group(g) * self.get_flux_by_value(cell_id, g, time)

                self._power[time][cell_id] = fission_rate * material.get_energy_per_fission() * cell_volume

    def initialize_field_variables(self):

        clock = Clock()
        for position in clock.get_positions():
            self._flux[position] = np.ones(self._num_x * self._num_y * self._num_amp_energy_groups)
            self._current[position] = np.zeros(self._num_x * self._num_y * self._num_amp_energy_groups * 4)
            self._temperature[position] = np.zeros(self._num_x * self._num_y)
            self._power[position] = np.zeros(self._num_x * self._num_y)

    def clone(self):

        mesh = AmpMesh(name=self._name, width=self.get_width(), height=self.get_height(),
                       num_x=self._num_x, num_y=self._num_y)
        return mesh

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

    def initialize_cells(self):

        self._cells = np.empty(shape=[self._num_y*self._num_x], dtype=object)

        # create cell objects
        for y in range(self._num_y):
            for x in range(self._num_x):
                self._cells[y*self._num_x + x] = CmfdCell()

    def initialize_surfaces(self):

        for y in xrange(self._num_y):
            for x in xrange(self._num_x):

                cell = self._cells[y*self._num_x + x]

                # Surface 0
                if x == 0:
                    cell._surfaces[0] = Surface(cell.get_material().get_num_energy_groups())
                else:
                    cell._surfaces[0] = self._cells[y*self._num_x + x - 1]._surfaces[2]

                # Surface 1
                if y == 0:
                    cell._surfaces[1] = Surface(cell.get_material().get_num_energy_groups())
                else:
                    cell._surfaces[1] = self._cells[(y - 1)*self._num_x + x]._surfaces[3]

                # Surface 2
                cell._surfaces[2] = Surface(cell.get_material().get_num_energy_groups())

                # Surface 3
                cell._surfaces[3] = Surface(cell.get_material().get_num_energy_groups())

    def get_flux_by_value(self, cell, group, time='CURRENT'):

        cv.check_clock_position(time, 'StructuredShapeMesh flux by value')
        return self._flux[time][cell * self._num_shape_energy_groups + group]

    def get_temperature_by_value(self, cell, time='CURRENT'):

        cv.check_clock_position(time, 'StructuredShapeMesh temperature by value')
        return self._temperature[time][cell]

    def get_power_by_value(self, cell, time='CURRENT'):

        cv.check_clock_position(time, 'StructuredShapeMesh power by value')
        return self._power[time][cell]

    def initialize_field_variables(self):

        clock = Clock()
        for position in clock.get_positions():
            self._flux[position] = np.ones(self._num_x * self._num_y * self._num_shape_energy_groups)
            self._temperature[position] = np.zeros(self._num_x * self._num_y)
            self._power[position] = np.zeros(self._num_x * self._num_y)

    def clone(self):

        mesh = StructuredShapeMesh(name=self._name, width=self.get_width(), height=self.get_height(),
                                   num_x=self._num_x, num_y=self._num_y)
        return mesh

    def compute_power(self, time='CURRENT'):

        cv.check_clock_position(time, 'StructuredShapeMesh compute power')
        cell_volume = self.get_cell_width() * self.get_cell_height()

        for y in xrange(self._num_y):
            for x in xrange(self._num_x):

                # Get cell properties
                cell_id = y*self._num_x + x
                material = self._cells[y*self._num_x + x].get_material()
                fission_rate = 0.0

                # Compute the fission rate
                for g in xrange(self._num_shape_energy_groups):
                    fission_rate += material.get_sigma_f_by_group(g) * self.get_flux_by_value(cell_id, g, time)

                self._power[time][cell_id] = fission_rate * material.get_energy_per_fission() * cell_volume

    def compute_initial_precursor_conc(self, time='CURRENT'):

        cv.check_clock_position(time, 'StructuredShapeMesh compute initial precursor conc')

        for y in xrange(self._num_y):
            for x in xrange(self._num_x):

                # Get cell properties
                cell_id = y*self._num_x + x
                material = self._cells[y*self._num_x + x].get_material()
                fission_rate = 0.0

                if isinstance(material, TransientMaterial):

                    # Compute the fission rate
                    for g in xrange(self._num_shape_energy_groups):
                        fission_rate += material.get_nu_sigma_f_by_group(g) * self.get_flux_by_value(cell_id, g, time)

                    for d in xrange(material.get_num_delayed_groups()):
                        conc = fission_rate * material.get_delayed_fraction_by_group(d, time) \
                            / self.get_k_eff_0() / material.get_decay_constant_by_group(d, time)
                        material.set_precursor_conc_by_group(conc, d, time)

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

    def initialize_cells(self):

        self._cells = np.empty(shape=[self._num_cells], dtype=object)

        # create cell objects
        for i in range(self._num_cells):
            self._cells[i] = Cell()

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
            self._temperature[position] = np.zeros(self._num_x * self._num_y)
            self._power[position] = np.zeros(self._num_x * self._num_y)

    def compute_power(self, time='CURRENT'):

        cv.check_clock_position(time, 'UnstructuredShapeMesh compute power')

        for i in xrange(self._num_cells):

            # Get cell properties
            material = self._cells[i].get_material()
            fission_rate = 0.0

            # Compute the fission rate
            for g in xrange(self._num_shape_energy_groups):
                fission_rate += material.get_sigma_f_by_group(g) * self.get_flux_by_value(i, g, time)

            self._power[time][i] = fission_rate * material.get_energy_per_fission() * self._volumes[i]

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