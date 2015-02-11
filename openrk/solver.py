__author__ = 'Samuel Shaner'
__email__ = 'shaner@mit.edu'

# Import modules
import numpy as np
import checkvalue as cv
import openrk
from material import TransientMaterial
from numba import jit


class Solver(object):
    def __init__(self):
        self._A_shape = None
        self._AM_shape = None
        self._M_shape = None
        self._b_shape = None
        self._A_amp = None
        self._AM_amp = None
        self._M_amp = None
        self._b_amp = None
        self._num_threads = None
        self._k_eff = None

        openrk.setNumThreads(1)

    def set_num_threads(self, num_threads):
        # Check input values
        cv.check_is_int(num_threads, 'Solver number of threads', 'num threads')

        self._num_threads = num_threads


class CmfdSolver(Solver):
    def __init__(self, shape_mesh, amp_mesh=None):

        super(CmfdSolver, self).__init__()

        self._shape_mesh = shape_mesh
        nc = self._shape_mesh.get_num_x() * self._shape_mesh.get_num_y()
        ng = self._shape_mesh.get_num_shape_energy_groups()
        self._A_shape = np.zeros((nc, ng * (ng + 4)))
        self._AM_shape = np.zeros((nc, ng * (ng + 4)))
        self._M_shape = np.zeros((nc, ng * ng))
        self._b_shape = np.zeros(nc * ng)

        if amp_mesh is not None:
            self._amp_mesh = amp_mesh
            nc = self._amp_mesh.get_num_x() * self._amp_mesh.get_num_y()
            ng = self._shape_mesh.get_num_amp_energy_groups()
            self._A_amp = np.zeros((nc, ng * (ng + 4)))
            self._AM_amp = np.zeros((nc, ng * (ng + 4)))
            self._M_amp = np.zeros((nc, ng * ng))
            self._b_amp = np.zeros(nc * ng)
        else:
            self._amp_mesh = None
            self._A_amp = None
            self._AM_amp = None
            self._M_amp = None
            self._b_amp = None

    def get_am_shape(self):

        return self._AM_shape

    def get_am_amp(self):

        return self._AM_amp

    def get_a_shape(self):

        return self._A_shape

    def get_a_amp(self):

        return self._A_amp

    def get_m_shape(self):

        return self._M_shape

    def get_m_amp(self):

        return self._M_amp

    def get_b_shape(self):

        return self._b_shape

    def get_b_amp(self):

        return self._b_amp

    def make_am_shape_initial(self, time='CURRENT'):

        nx = self._shape_mesh.get_num_x()
        ny = self._shape_mesh.get_num_y()
        ng = self._shape_mesh.get_num_shape_energy_groups()
        height = self._shape_mesh.get_cell_height()
        width = self._shape_mesh.get_cell_width()
        volume = width * height
        mesh = self._shape_mesh
        temps = mesh.get_temperature(time)

        # reinitialize matrices to zero
        self._A_shape.fill(0.0)
        self._M_shape.fill(0.0)

        for y in xrange(ny):
            for x in xrange(nx):

                cell = y*nx+x
                temp = temps[cell]

                for e in range(ng):

                    material = self._shape_mesh.get_material(y*nx+x)

                    # absorption term on diagonal
                    self._A_shape[y * nx + x][e * (ng + 4) + e + 2] += material.get_sigma_a_by_group(e, time, temp) \
                        * volume

                    # buckling term on diagonal
                    self._A_shape[y * nx + x][e * (ng + 4) + e + 2] += material.get_dif_coef_by_group(e, time, temp) * \
                        self._shape_mesh.get_buckling_by_group(e) * volume

                    # out scattering term on diagonal
                    for g in xrange(ng):
                        if e != g:
                            self._A_shape[y * nx + x][e * (ng + 4) + e + 2] += \
                                material.get_sigma_s_by_group(e, g, time, temp) * volume

                    # fission terms on diagonal
                    for g in xrange(ng):
                        self._M_shape[y * nx + x][e * ng + g] = material.get_chi_by_group(e, time, temp) * \
                            material.get_nu_sigma_f_by_group(g) * volume

                    # in scattering terms on off diagonals
                    for g in xrange(ng):
                        if e != g:
                            self._A_shape[y * nx + x][e * (ng + 4) + g + 2] -= \
                                material.get_sigma_s_by_group(g, e, time, temp) * volume

                    # RIGHT SURFACE

                    # transport term on diagonal
                    self._A_shape[y * nx + x][e * (ng + 4) + e + 2] += mesh.get_dif_linear_by_value(cell, e, 2) * height

                    # transport terms on off diagonals
                    if x != nx - 1:
                        self._A_shape[y * nx + x][e * (ng + 4) + ng + 2] -= mesh.get_dif_linear_by_value(cell, e, 2) \
                            * height

                    # LEFT SURFACE

                    # transport term on diagonal
                    self._A_shape[y * nx + x][e * (ng + 4) + e + 2] += mesh.get_dif_linear_by_value(cell, e, 0) * height

                    # transport terms on off diagonals
                    if x != 0:
                        self._A_shape[y * nx + x][e * (ng + 4)] -= mesh.get_dif_linear_by_value(cell, e, 0) * height

                    # BOTTOM SURFACE

                    # transport term on diagonal
                    self._A_shape[y * nx + x][e * (ng + 4) + e + 2] += mesh.get_dif_linear_by_value(cell, e, 1) * width

                    # transport terms on off diagonals
                    if y != 0:
                        self._A_shape[y * nx + x][e * (ng + 4) + 1] -= mesh.get_dif_linear_by_value(cell, e, 1) * width

                    # TOP SURFACE

                    # transport term on diagonal
                    self._A_shape[y * nx + x][e * (ng + 4) + e + 2] += mesh.get_dif_linear_by_value(cell, e, 3) * width

                    # transport terms on off diagonals
                    if y != ny - 1:
                        self._A_shape[y * nx + x][e * (ng + 4) + ng + 3] -= mesh.get_dif_linear_by_value(cell, e, 3) \
                            * width

    @jit
    def compute_initial_shape(self, tol):

        # Compute the surface diffusion coefficients
        self._shape_mesh.compute_dif_coefs('CURRENT')

        # Create the production and loss matrices
        self.make_am_shape_initial('CURRENT')

        # Initialize the old source and a temporary flux variable
        nx = self._shape_mesh.get_num_x()
        ny = self._shape_mesh.get_num_y()
        ng = self._shape_mesh.get_num_shape_energy_groups()
        old_source = np.zeros(nx * ny * ng)
        flux_temp = np.zeros(nx * ny * ng)

        # solve the eigenvalue problem
        self._k_eff = openrk.eigenvalueSolve(self._A_shape, self._M_shape, self._shape_mesh.get_flux('CURRENT'),
                                             self._b_shape, old_source, flux_temp, ng, nx, ny, tol)

    @jit
    def make_am_amp(self, wt=0.5):

        mesh = self._amp_mesh
        am = self._AM_amp
        b = self._b_amp

        nx = mesh.get_num_x()
        ny = mesh.get_num_y()
        ng = mesh.get_num_amp_energy_groups()
        height = mesh.get_cell_height()
        width = mesh.get_cell_width()
        volume = width * height
        flux = mesh.get_flux('PREVIOUS_IN')
        dt = mesh.get_clock().get_dt_inner()
        temps = mesh.get_temperature('CURRENT')
        temps_prev = mesh.get_temperature('PREVIOUS_IN')

        # reinitialize matrices to zero
        am.fill(0.0)
        b.fill(0.0)

        for y in xrange(ny):
            for x in xrange(nx):

                cell = y*nx+x
                temp = temps[cell]
                temp_prev = temps_prev[cell]
                material = mesh.get_material(cell)

                beta = 0.0
                beta_prev = 0.0
                for e in xrange(mesh.get_num_delayed_groups()):
                    beta += mesh.get_delayed_fraction_by_group(e)
                    beta_prev += mesh.get_delayed_fraction_by_group(e)

                for e in range(ng):

                    # time absorption term on diagonal
                    am[cell][e * (ng + 4) + e + 2] += 1.0 / dt / \
                        material.get_velocity_by_group(e, 'CURRENT', temp) * volume
                    b[cell*ng+e] += flux[cell*ng+e] / dt / \
                        material.get_velocity_by_group(e, 'PREVIOUS_IN', temp_prev) * volume

                    # delayed neutron precursors
                    for d in xrange(mesh.get_num_delayed_groups()):
                        b[cell*ng+e] += wt * material.get_chi_by_group(e, 'CURRENT', temp) * \
                            mesh.get_decay_constant_by_group(d) * \
                            material.get_precursor_conc_by_group(d, 'CURRENT') * volume
                        b[cell*ng+e] += (1-wt) * material.get_chi_by_group(e, 'PREVIOUS_IN', temp_prev) * \
                            mesh.get_decay_constant_by_group(d) * \
                            material.get_precursor_conc_by_group(d, 'PREVIOUS_IN') * volume

                    # absorption term on diagonal
                    am[cell][e * (ng + 4) + e + 2] += volume * \
                        wt * material.get_sigma_a_by_group(e, 'CURRENT', temp)
                    b[cell*ng+e] -= volume * \
                        (1-wt) * material.get_sigma_a_by_group(e, 'PREVIOUS_IN', temp_prev) * flux[cell*ng+e]

                    # buckling term on diagonal
                    am[cell][e * (ng + 4) + e + 2] += wt * material.get_dif_coef_by_group(e, 'CURRENT', temp) * \
                        self._shape_mesh.get_buckling_by_group(e) * volume
                    b[cell*ng+e] -= (1-wt) * material.get_dif_coef_by_group(e, 'PREVIOUS_IN', temp_prev) * \
                        self._shape_mesh.get_buckling_by_group(e) * volume * flux[cell*ng+e]

                    # out scattering term on diagonal
                    for g in xrange(ng):
                        if e != g:
                            am[cell][e * (ng + 4) + e + 2] += volume * \
                                wt * material.get_sigma_s_by_group(e, g, 'CURRENT', temp)
                            b[cell*ng+e] -= volume * \
                                (1-wt) * material.get_sigma_s_by_group(e, g, 'PREVIOUS_IN', temp_prev) * flux[cell*ng+e]

                    # fission terms on diagonal
                    for g in xrange(ng):
                        am[cell][e * (ng+4) + g + 2] -= wt * (1 - beta) * material.get_chi_by_group(e, 'CURRENT', temp) * \
                            material.get_nu_sigma_f_by_group(g, 'CURRENT', temp) / mesh.get_k_eff_0() * volume
                        b[cell*ng+e] += (1-wt) * (1 - beta_prev) * material.get_chi_by_group(e, 'PREVIOUS_IN', temp_prev) * \
                            material.get_nu_sigma_f_by_group(g, 'PREVIOUS_IN', temp_prev) / mesh.get_k_eff_0() \
                            * flux[cell*ng+g] * volume

                    # in scattering terms on off diagonals
                    for g in xrange(ng):
                        if e != g:
                            am[cell][e * (ng + 4) + g + 2] -= volume * \
                                wt * material.get_sigma_s_by_group(g, e, 'CURRENT', temp)
                            b[cell*ng+e] += volume * \
                                (1-wt) * material.get_sigma_s_by_group(g, e, 'PREVIOUS_IN', temp_prev) * flux[cell*ng+g]

                    # TRANSPORT TO ADJACENT CELLS

                    # RIGHT SURFACE

                    # transport term on diagonal
                    am[cell][e * (ng + 4) + e + 2] += wt * \
                        (mesh.get_dif_linear_by_value(cell, e, 2, 'CURRENT') -
                         mesh.get_dif_nonlinear_by_value(cell, e, 2, 'CURRENT')) * height
                    b[cell*ng+e] -= (1-wt) * \
                        (mesh.get_dif_linear_by_value(cell, e, 2, 'PREVIOUS_IN') -
                         mesh.get_dif_nonlinear_by_value(cell, e, 2, 'PREVIOUS_IN')) * height * flux[cell*ng+e]

                    # transport terms on off diagonals
                    if x != nx - 1:
                        am[cell][e * (ng + 4) + ng + 2] -= wt * \
                            (mesh.get_dif_linear_by_value(cell, e, 2, 'CURRENT') +
                             mesh.get_dif_nonlinear_by_value(cell, e, 2, 'CURRENT')) * height
                        b[cell*ng+e] += (1-wt) * \
                            (mesh.get_dif_linear_by_value(cell, e, 2, 'PREVIOUS_IN') +
                             mesh.get_dif_nonlinear_by_value(cell, e, 2, 'PREVIOUS_IN')) * height * flux[(cell+1)*ng+e]

                    # LEFT SURFACE

                    # transport term on diagonal
                    am[cell][e * (ng + 4) + e + 2] += wt * \
                        (mesh.get_dif_linear_by_value(cell, e, 0, 'CURRENT') +
                         mesh.get_dif_nonlinear_by_value(cell, e, 0, 'CURRENT')) * height
                    b[cell*ng+e] -= (1-wt) * \
                        (mesh.get_dif_linear_by_value(cell, e, 0, 'PREVIOUS_IN') +
                         mesh.get_dif_nonlinear_by_value(cell, e, 0, 'PREVIOUS_IN')) * height * flux[cell*ng+e]

                    # transport terms on off diagonals
                    if x != 0:
                        am[cell][e * (ng + 4)] -= wt * \
                            (mesh.get_dif_linear_by_value(cell, e, 0, 'CURRENT') -
                             mesh.get_dif_nonlinear_by_value(cell, e, 0, 'CURRENT')) * height
                        b[cell*ng+e] += (1-wt) * \
                            (mesh.get_dif_linear_by_value(cell, e, 0, 'PREVIOUS_IN') -
                             mesh.get_dif_nonlinear_by_value(cell, e, 0, 'PREVIOUS_IN')) * height * flux[(cell-1)*ng+e]

                    # BOTTOM SURFACE

                    # transport term on diagonal
                    am[cell][e * (ng + 4) + e + 2] += wt * \
                        (mesh.get_dif_linear_by_value(cell, e, 1, 'CURRENT') +
                         mesh.get_dif_nonlinear_by_value(cell, e, 1, 'CURRENT')) * width
                    b[cell*ng+e] -= (1-wt) * \
                        (mesh.get_dif_linear_by_value(cell, e, 1, 'PREVIOUS_IN') +
                         mesh.get_dif_nonlinear_by_value(cell, e, 1, 'PREVIOUS_IN')) * width * flux[cell*ng+e]

                    # transport terms on off diagonals
                    if y != 0:
                        am[cell][e * (ng + 4) + 1] -= wt * \
                            (mesh.get_dif_linear_by_value(cell, e, 1, 'CURRENT') -
                             mesh.get_dif_nonlinear_by_value(cell, e, 1, 'CURRENT')) * width
                        b[cell*ng+e] += (1-wt) * \
                            (mesh.get_dif_linear_by_value(cell, e, 1, 'PREVIOUS_IN') -
                             mesh.get_dif_nonlinear_by_value(cell, e, 1, 'PREVIOUS_IN')) * width * flux[(cell-nx)*ng+e]

                    # TOP SURFACE

                    # transport term on diagonal
                    am[cell][e * (ng + 4) + e + 2] += wt * \
                        (mesh.get_dif_linear_by_value(cell, e, 3, 'CURRENT') -
                         mesh.get_dif_nonlinear_by_value(cell, e, 3, 'CURRENT')) * width
                    b[cell*ng+e] -= (1-wt) * \
                        (mesh.get_dif_linear_by_value(cell, e, 3, 'PREVIOUS_IN') -
                         mesh.get_dif_nonlinear_by_value(cell, e, 3, 'PREVIOUS_IN')) * width * flux[cell*ng+e]

                    # transport terms on off diagonals
                    if y != ny - 1:
                        am[cell][e * (ng + 4) + ng + 3] -= wt * \
                            (mesh.get_dif_linear_by_value(cell, e, 3, 'CURRENT') +
                             mesh.get_dif_nonlinear_by_value(cell, e, 3, 'CURRENT')) * width
                        b[cell*ng+e] += (1-wt) * \
                            (mesh.get_dif_linear_by_value(cell, e, 3, 'PREVIOUS_IN') +
                             mesh.get_dif_nonlinear_by_value(cell, e, 3, 'PREVIOUS_IN')) * width * flux[(cell+nx)*ng+e]

    @jit
    def make_am_shape(self, wt=0.5):

        mesh = self._shape_mesh
        am = self._AM_shape
        b = self._b_shape

        nx = mesh.get_num_x()
        ny = mesh.get_num_y()
        ng = mesh.get_num_shape_energy_groups()
        height = mesh.get_cell_height()
        width = mesh.get_cell_width()
        volume = width * height
        flux = mesh.get_flux('PREVIOUS_OUT')
        dt = mesh.get_clock().get_dt_outer()
        temps_prev = mesh.get_temperature('PREVIOUS_OUT')
        temps = mesh.get_temperature('FORWARD_OUT')

        # reinitialize matrices to zero
        am.fill(0.0)
        b.fill(0.0)

        for y in xrange(ny):
            for x in xrange(nx):

                cell = y*nx+x
                temp = temps[cell]
                temp_prev = temps_prev[cell]
                material = mesh.get_material(y*nx+x)

                beta = 0.0
                beta_prev = 0.0
                for e in xrange(mesh.get_num_delayed_groups()):
                    beta += mesh.get_delayed_fraction_by_group(e)
                    beta_prev += mesh.get_delayed_fraction_by_group(e)

                for e in range(ng):

                    # time absorption term on diagonal
                    am[cell][e * (ng + 4) + e + 2] += 1.0 / dt / \
                        material.get_velocity_by_group(e, 'FORWARD_OUT', temp) * volume
                    b[cell*ng+e] += flux[cell*ng+e] / dt / \
                        material.get_velocity_by_group(e, 'PREVIOUS_OUT', temp_prev) * volume

                    # delayed neutron precursors
                    if isinstance(material, TransientMaterial):
                        for d in xrange(mesh.get_num_delayed_groups()):
                            b[cell*ng+e] += wt * material.get_chi_by_group(e, 'FORWARD_OUT', temp) * \
                                mesh.get_decay_constant_by_group(d) * \
                                material.get_precursor_conc_by_group(d, 'FORWARD_OUT') * volume
                            b[cell*ng+e] += (1-wt) * material.get_chi_by_group(e, 'PREVIOUS_OUT', temp_prev) * \
                                mesh.get_decay_constant_by_group(d) * \
                                material.get_precursor_conc_by_group(d, 'PREVIOUS_OUT') * volume

                    # absorption term on diagonal
                    am[cell][e * (ng + 4) + e + 2] += volume * \
                        wt * material.get_sigma_a_by_group(e, 'FORWARD_OUT', temp)
                    b[cell*ng+e] -= volume * \
                        (1-wt) * material.get_sigma_a_by_group(e, 'PREVIOUS_OUT', temp_prev) * flux[cell*ng+e]

                    # buckling term on diagonal
                    am[cell][e * (ng + 4) + e + 2] += wt * material.get_dif_coef_by_group(e, 'FORWARD_OUT', temp) * \
                        volume * self._amp_mesh.get_buckling_by_group(e)
                    b[cell*ng+e] -= (1-wt) * material.get_dif_coef_by_group(e, 'PREVIOUS_OUT', temp_prev) * \
                        flux[cell*ng+e] * volume * self._amp_mesh.get_buckling_by_group(e)

                    # out scattering term on diagonl
                    for g in xrange(ng):
                        if e != g:
                            am[cell][e * (ng + 4) + e + 2] += volume * \
                                wt * material.get_sigma_s_by_group(e, g, 'FORWARD_OUT', temp)
                            b[cell*ng+e] -= volume * \
                                (1-wt) * material.get_sigma_s_by_group(e, g, 'PREVIOUS_OUT', temp_prev) * flux[cell*ng+e]

                    # fission terms on diagonal
                    for g in xrange(ng):
                        am[cell][e * (ng+4) + g + 2] -= wt * (1 - beta) * material.get_chi_by_group(e, 'FORWARD_OUT', temp) * \
                            material.get_nu_sigma_f_by_group(g, 'FORWARD_OUT', temp) / mesh.get_k_eff_0() * volume
                        b[cell*ng+e] += (1-wt) * (1 - beta_prev) * material.get_chi_by_group(e, 'PREVIOUS_OUT', temp_prev) * \
                            material.get_nu_sigma_f_by_group(g, 'PREVIOUS_OUT', temp_prev) / mesh.get_k_eff_0() * \
                            flux[cell*ng+g] * volume

                    # in scattering terms on off diagonals
                    for g in xrange(ng):
                        if e != g:
                            am[cell][e * (ng + 4) + g + 2] -= \
                                wt * material.get_sigma_s_by_group(g, e, 'FORWARD_OUT', temp) * volume
                            b[cell*ng+e] += volume * \
                                (1-wt) * material.get_sigma_s_by_group(g, e, 'PREVIOUS_OUT', temp_prev) * flux[cell*ng+g]

                    # TRANSPORT TO ADJACENT CELLS

                    # RIGHT SURFACE

                    # transport term on diagonal
                    am[cell][e * (ng + 4) + e + 2] += wt * \
                        mesh.get_dif_linear_by_value(cell, e, 2, 'FORWARD_OUT') * height
                    b[cell*ng+e] -= (1-wt) * \
                        mesh.get_dif_linear_by_value(cell, e, 2, 'PREVIOUS_OUT') * height * flux[cell*ng+e]

                    # transport terms on off diagonals
                    if x != nx - 1:
                        am[cell][e * (ng + 4) + ng + 2] -= wt * \
                            mesh.get_dif_linear_by_value(cell, e, 2, 'FORWARD_OUT') * height
                        b[cell*ng+e] += (1-wt) * \
                            mesh.get_dif_linear_by_value(cell, e, 2, 'PREVIOUS_OUT') * height * flux[(cell+1)*ng+e]

                    # LEFT SURFACE

                    # transport term on diagonal
                    am[cell][e * (ng + 4) + e + 2] += wt * \
                        mesh.get_dif_linear_by_value(cell, e, 0, 'FORWARD_OUT') * height
                    b[cell*ng+e] -= (1-wt) * \
                        mesh.get_dif_linear_by_value(cell, e, 0, 'PREVIOUS_OUT') * height * flux[cell*ng+e]

                    # transport terms on off diagonals
                    if x != 0:
                        am[cell][e * (ng + 4)] -= wt * \
                            mesh.get_dif_linear_by_value(cell, e, 0, 'FORWARD_OUT') * height
                        b[cell*ng+e] += (1-wt) * \
                            mesh.get_dif_linear_by_value(cell, e, 0, 'PREVIOUS_OUT') * height * flux[(cell-1)*ng+e]

                    # BOTTOM SURFACE

                    # transport term on diagonal
                    am[cell][e * (ng + 4) + e + 2] += wt * \
                        mesh.get_dif_linear_by_value(cell, e, 1, 'FORWARD_OUT') * width
                    b[cell*ng+e] -= (1-wt) * \
                        mesh.get_dif_linear_by_value(cell, e, 1, 'PREVIOUS_OUT') * width * flux[cell*ng+e]

                    # transport terms on off diagonals
                    if y != 0:
                        am[cell][e * (ng + 4) + 1] -= wt * \
                            mesh.get_dif_linear_by_value(cell, e, 1, 'FORWARD_OUT') * width
                        b[cell*ng+e] += (1-wt) * \
                            mesh.get_dif_linear_by_value(cell, e, 1, 'PREVIOUS_OUT') * width * flux[(cell-nx)*ng+e]

                    # TOP SURFACE

                    # transport term on diagonal
                    am[cell][e * (ng + 4) + e + 2] += wt * \
                        mesh.get_dif_linear_by_value(cell, e, 3, 'FORWARD_OUT') * width
                    b[cell*ng+e] -= (1-wt) * \
                        mesh.get_dif_linear_by_value(cell, e, 3, 'PREVIOUS_OUT') * width * flux[cell*ng+e]

                    # transport terms on off diagonals
                    if y != ny - 1:
                        am[cell][e * (ng + 4) + ng + 3] -= wt * \
                            mesh.get_dif_linear_by_value(cell, e, 3, 'FORWARD_OUT') * width
                        b[cell*ng+e] += (1-wt) * \
                            mesh.get_dif_linear_by_value(cell, e, 3, 'PREVIOUS_OUT') * width * flux[(cell+nx)*ng+e]