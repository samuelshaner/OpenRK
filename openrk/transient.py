__author__ = 'Samuel Shaner'
__email__ = 'shaner@mit.edu'

# Import modules
import numpy as np
import checkvalue as cv
from solver import Solver, CmfdSolver
from clock import Clock
from mesh import Mesh
import openrk
from material import TransientMaterial


class Transient(object):
    """
    This class contains information on the transient problem setup
    and methods for stepping forward in time.

    Attributes:
      _initial_power
      _k_eff_0
      _name
      _clock
      _id
      _solver
      _shape_mesh
      _amp_mesh

    Setter Methods:
      set_name(str name)
      set_id(int mesh_id)
      set_initial_power(float initial_power)
      set_solver(Solver solver)
      set_shape_mesh(Mesh mesh)
      set_amp_mesh(AmpMesh mesh)
      set_clock()

    Getter Methods:
      get_name()
      get_id()
      get_initial_power()
      get_solver()
      get_shape_mesh()
      get_amp_mesh()

    Other Methods:
      compute_initial_shape()
      take_outer_step()
      take_inner_step()
    """
    def __init__(self):

        # Initialize class attributes
        self._initial_power = None
        self._clock = None
        self._k_eff_0 = None
        self._name = None
        self._solver = None
        self._shape_mesh = None
        self._amp_mesh = None

    def set_initial_power(self, initial_power):

        self._initial_power = initial_power

    def set_clock(self, clock):

        if not isinstance(clock, Clock):
            msg = 'Cannot set clock for Transient since input is not of type Clock: {0}'.format(clock)
            raise ValueError(msg)
        else:
            self._clock = clock

    def set_name(self, name):

        if not cv.is_string(name):
            msg = 'Unable to set name for Transient with a non-string ' \
                  'value {1}'.format(name)
            raise ValueError(msg)

        else:
            self._name = name

    def set_solver(self, solver):

        if not isinstance(solver, Solver):
            msg = 'Cannot set solver for Transient since input is not of type Solver: {0}'.format(solver)
            raise ValueError(msg)
        else:
            self._solver = solver

    def set_shape_mesh(self, mesh):

        if not isinstance(mesh, Mesh):
            msg = 'Cannot set shape mesh for Transient since input is not of type Mesh: {0}'.format(mesh)
            raise ValueError(msg)
        else:
            self._shape_mesh = mesh

    def set_amp_mesh(self, mesh):

        if not isinstance(mesh, Mesh):
            msg = 'Cannot set amp mesh for Transient since input is not of type Mesh: {0}'.format(mesh)
            raise ValueError(msg)
        else:
            self._amp_mesh = mesh

    def compute_initial_shape(self):

        # Give the amp and shape mesh a common clock
        self._amp_mesh.set_clock(self._clock)
        self._shape_mesh.set_clock(self._clock)

        # Compute the initial shape
        self._solver.compute_initial_shape(1.e-12)

        # Normalize the flux by the user-set initial power
        self._shape_mesh.set_k_eff_0(self._solver._k_eff)
        self._amp_mesh.set_k_eff_0(self._solver._k_eff)
        initial_power = self._shape_mesh.get_average_power()
        self._shape_mesh.scale_flux(self._initial_power / initial_power)

        # Compute hte initial precursor concentrations on the shape mesh
        self._shape_mesh.compute_initial_precursor_conc()

        # Condense the materials and compute the dif coefs for the amp mesh
        self._amp_mesh.condense_materials(time='CURRENT', save_flux=True)
        self._amp_mesh.compute_current()
        self._solver.compute_dif_coefs_amp()

        # broadcast current properties to all other times
        self.broadcast_to_all()

    def take_inner_step(self):

        # Increment CURRENT time by dt_inner
        self._clock.take_inner_step()

        # Interpolate the nonlinear diffusion coefficients
        self._amp_mesh.interpolate_dif_nonlinear(time_begin='PREVIOUS_OUT', time_end='FORWARD_OUT', time='CURRENT')

        #  Synthesize fine mesh flux at CURRENT time
        self._shape_mesh.synthesize_flux(time='CURRENT')

        # Compute the temperatures at the CURRENT time
        self._shape_mesh.integrate_temperature(time_from='PREVIOUS_IN', time_to='CURRENT')

        # Compute the precursor concentrations at the CURRENT time
        self._shape_mesh.integrate_precursor_conc(time_from='PREVIOUS_IN', time_to='CURRENT')

        # Condense precursor conc and xs to amp mesh
        self._amp_mesh.condense_materials(time='CURRENT')

        converged = False
        ng = self._amp_mesh.get_num_amp_energy_groups()
        nx = self._amp_mesh.get_num_x()
        ny = self._amp_mesh.get_num_y()
        flux_temp = np.zeros(ng*nx*ny)
        tol = 1.e-5

        while not converged:

            # Copy previous flux guess to FORWARD_IN_OLD
            self._amp_mesh.copy_flux(time_from='CURRENT', time_to='FORWARD_IN_OLD')

            # Update the amplitude at CURRENT time
            self._solver.make_am_amp(wt=0.5)
            openrk.linearSolve(self._solver.get_am_amp(), self._amp_mesh.get_flux('CURRENT'),
                               self._solver.get_b_amp(), flux_temp, nx, ny, ng, 1.e-8)

            #  Synthesize fine mesh flux at CURRENT time
            self._shape_mesh.synthesize_flux(time='CURRENT')

            # Compute the temperatures at the CURRENT time
            self._shape_mesh.integrate_temperature(time_from='PREVIOUS_IN', time_to='CURRENT')

            # Compute the precursor concentrations at the CURRENT time
            self._shape_mesh.integrate_precursor_conc(time_from='PREVIOUS_IN', time_to='CURRENT')

            # Condense precursor conc and xs to amp mesh
            self._amp_mesh.condense_materials(time='CURRENT')

            residual = self._amp_mesh.compute_flux_l2_norm(time_1='CURRENT', time_2='FORWARD_IN_OLD')

            power = self._shape_mesh.get_average_power()
            print 'TIME = {:1.4f}, POWER = {:.6e}, RESIDUAL = {:.6e}'.format(self._clock.get_time('CURRENT'), power, residual)

            if residual < tol:
                converged = True

        self.broadcast_to_one(time_from='CURRENT', time_to='PREVIOUS_IN')

    def take_outer_step(self):

        # Increment FORWARD_OUT time by dt_outer
        self._clock.take_outer_step()

        # broadcast the current state to all active times
        self.broadcast_to_active('FORWARD_OUT')
        ng = self._shape_mesh.get_num_shape_energy_groups()
        nx = self._shape_mesh.get_num_x()
        ny = self._shape_mesh.get_num_y()
        flux_temp = np.zeros(ng*nx*ny)
        tol = 1.e-6

        # Prolongate the amplitude
        while self._clock.get_time('CURRENT') - self._clock.get_time('FORWARD_OUT') < -1.e-6:
            self.take_inner_step()

        # Reconstruct the prolongated fine mesh flux and save to FORWARD_OUT
        self._shape_mesh.reconstruct_flux(time='CURRENT', time_shape='FORWARD_OUT', time_amp='CURRENT')
        self.broadcast_to_one(time_from='CURRENT', time_to='FORWARD_OUT')
        self.broadcast_to_one(time_from='PREVIOUS_OUT', time_to='CURRENT')
        self.broadcast_to_one(time_from='PREVIOUS_OUT', time_to='PREVIOUS_IN')

        # Iteratively take outer step until converged
        while True:

            # Copy current approx of fine mesh flux at FORWARD_OUT to FORWARD_OUT_OLD for convergence check
            self._shape_mesh.copy_flux(time_from='FORWARD_OUT', time_to='FORWARD_OUT_OLD')

            # Update the dif coefs at the FORWARD_OUT time
            self._solver.compute_surface_dif_coefs_shape(time='FORWARD_OUT')

            # Recreate the shape matrices
            self._solver.make_am_shape(0.5)

            # Solve for the updated flux at FORWARD_OUT
            openrk.linearSolve(self._solver.get_am_shape(), self._shape_mesh.get_flux('FORWARD_OUT'),
                               self._solver.get_b_shape(), flux_temp, nx, ny, ng, 1.e-8)

            # Update the amp mesh dif coefs at FORWARD_OUT
            self._amp_mesh.compute_current(time='FORWARD_OUT')
            self._amp_mesh.condense_materials(time='FORWARD_OUT', save_flux=True)
            self._solver.compute_dif_coefs_amp(time='FORWARD_OUT')

            # Reset CURRENT time
            self._clock.reset_to_previous_outer_step()

            # March the amplitude forward until end of outer step reached
            while self._clock.get_time('CURRENT') - self._clock.get_time('FORWARD_OUT') < -1.e-6:
                self.take_inner_step()

            # Copy CURRENT amp to FORWARD_OUT
            self._shape_mesh.reconstruct_flux(time='CURRENT', time_shape='FORWARD_OUT', time_amp='CURRENT')
            self.broadcast_to_one(time_from='CURRENT', time_to='FORWARD_OUT')

            # Check for convergence
            residual = self._shape_mesh.compute_flux_l2_norm(time_1='FORWARD_OUT', time_2='FORWARD_OUT_OLD')
            print 'outer residual ' + str(residual)

            if residual < tol:
                break
            else:
                self.broadcast_to_one(time_from='PREVIOUS_OUT', time_to='CURRENT')
                self.broadcast_to_one(time_from='PREVIOUS_OUT', time_to='PREVIOUS_IN')

    def take_outer_step_only(self):

        # Increment FORWARD_OUT time by dt_outer
        self._clock.take_outer_step()

        # broadcast the current state to all active times
        self.broadcast_to_active('FORWARD_OUT')
        converged = False
        ng = self._shape_mesh.get_num_shape_energy_groups()
        nx = self._shape_mesh.get_num_x()
        ny = self._shape_mesh.get_num_y()
        flux_temp = np.zeros(ng*nx*ny)
        tol = 1.e-8

        # Iteratively take outer step until converged
        while not converged:

            # Compute the temperatures at the CURRENT time
            self._shape_mesh.integrate_temperature(time_from='PREVIOUS_OUT', time_to='FORWARD_OUT')

            # Compute the precursor concentrations at the CURRENT time
            self._shape_mesh.integrate_precursor_conc(time_from='PREVIOUS_OUT', time_to='FORWARD_OUT')

            # Copy current approx of fine mesh flux at FORWARD_OUT to FORWARD_OUT_OLD for convergence check
            self._shape_mesh.copy_flux(time_from='FORWARD_OUT', time_to='FORWARD_OUT_OLD')

            # Update the dif coefs at the FORWARD_OUT time
            self._solver.compute_surface_dif_coefs_shape(time='FORWARD_OUT')

            # Recreate the shape matrices
            self._solver.make_am_shape(1.0)

            # Solve for the updated flux at FORWARD_OUT
            openrk.linearSolve(self._solver.get_am_shape(), self._shape_mesh.get_flux('FORWARD_OUT'),
                               self._solver.get_b_shape(), flux_temp, nx, ny, ng, 1.e-8)

            # Check for convergence
            residual = self._shape_mesh.compute_flux_l2_norm(time_1='FORWARD_OUT', time_2='FORWARD_OUT_OLD')
            power = self._shape_mesh.get_average_power('FORWARD_OUT')
            print 'TIME = ' + str(self._clock.get_time('FORWARD_OUT')) + ' POWER = ' + str(power) + ' RESIDUAL = ' + str(residual)

            if residual < tol:
                converged = True
                self._clock.set_time(position_from='FORWARD_OUT', position_to='CURRENT')

    def broadcast_to_active(self, time='CURRENT'):

        for position in ['PREVIOUS_OUT', 'PREVIOUS_IN', 'CURRENT', 'FORWARD_IN_OLD', 'FORWARD_OUT', 'FORWARD_OUT_OLD']:

            self._amp_mesh.copy_flux(time, position)
            self._amp_mesh.copy_current(time, position)
            self._amp_mesh.copy_temperature(time, position)
            self._amp_mesh.copy_dif_linear(time, position)
            self._amp_mesh.copy_dif_nonlinear(time, position)

            self._shape_mesh.copy_flux(time, position)
            self._shape_mesh.copy_current(time, position)
            self._shape_mesh.copy_temperature(time, position)
            self._shape_mesh.copy_dif_linear(time, position)

            for i in xrange(self._amp_mesh.get_num_x() * self._amp_mesh.get_num_y()):
                mat = self._amp_mesh.get_material(i)
                if isinstance(mat, TransientMaterial):
                    mat.broadcast(time, position)

            for i in xrange(self._shape_mesh.get_num_x() * self._shape_mesh.get_num_y()):
                mat = self._shape_mesh.get_material(i)
                if isinstance(mat, TransientMaterial):
                    mat.broadcast(time, position)

    def broadcast_to_all(self, time='CURRENT'):

        for position in self._clock.get_positions():

            self._amp_mesh.copy_flux(time, position)
            self._amp_mesh.copy_current(time, position)
            self._amp_mesh.copy_temperature(time, position)
            self._amp_mesh.copy_dif_linear(time, position)
            self._amp_mesh.copy_dif_nonlinear(time, position)

            self._shape_mesh.copy_flux(time, position)
            self._shape_mesh.copy_current(time, position)
            self._shape_mesh.copy_temperature(time, position)
            self._shape_mesh.copy_dif_linear(time, position)

            for i in xrange(self._amp_mesh.get_num_x() * self._amp_mesh.get_num_y()):
                mat = self._amp_mesh.get_material(i)
                if isinstance(mat, TransientMaterial):
                    mat.broadcast(time, position)

            for i in xrange(self._shape_mesh.get_num_x() * self._shape_mesh.get_num_y()):
                mat = self._shape_mesh.get_material(i)
                if isinstance(mat, TransientMaterial):
                    mat.broadcast(time, position)

    def broadcast_to_one(self, time_from, time_to):

        self._amp_mesh.copy_flux(time_from, time_to)
        self._amp_mesh.copy_current(time_from, time_to)
        self._amp_mesh.copy_temperature(time_from, time_to)
        self._amp_mesh.copy_dif_linear(time_from, time_to)
        self._amp_mesh.copy_dif_nonlinear(time_from, time_to)

        self._shape_mesh.copy_flux(time_from, time_to)
        self._shape_mesh.copy_current(time_from, time_to)
        self._shape_mesh.copy_temperature(time_from, time_to)
        self._shape_mesh.copy_dif_linear(time_from, time_to)

        for i in xrange(self._amp_mesh.get_num_x() * self._amp_mesh.get_num_y()):
            mat = self._amp_mesh.get_material(i)
            if isinstance(mat, TransientMaterial):
                mat.broadcast(time_from, time_to)

        for i in xrange(self._shape_mesh.get_num_x() * self._shape_mesh.get_num_y()):
            mat = self._shape_mesh.get_material(i)
            if isinstance(mat, TransientMaterial):
                mat.broadcast(time_from, time_to)