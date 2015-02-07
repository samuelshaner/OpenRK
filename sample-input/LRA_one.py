import openrk as rk
import random

mesh = rk.mesh.StructuredShapeMesh(name='shape mesh', width=1.0, height=1.0, num_x=1, num_y=1)
mesh.set_num_amp_energy_groups(2)
mesh.set_num_shape_energy_groups(2)
mesh.set_num_delayed_groups(2)
mesh.set_buckling([1.e-4, 1.e-4])
mesh.set_delayed_fractions([0.0054, 0.001087])
mesh.set_decay_constants([0.00654, 1.35])
mesh.initialize()
mesh.set_temperature(300)

# create fuel 1 blade in
fuel1bin = rk.material.FunctionalMaterial(name='fuel 1 blade in')
fuel1bin.set_num_time_steps(3)
fuel1bin.set_num_energy_groups(2)
fuel1bin.set_num_delayed_groups(2)
fuel1bin.set_time_steps([0.0, 2.0, 3.0])
fuel1bin.set_sigma_a([[0.008252, 0.1003], [0.008252, 0.1003], [0.008252, 0.1003]])
fuel1bin.set_dif_coef([[1.255, 0.211], [1.255, 0.211], [1.255, 0.211]])
fuel1bin.set_sigma_t([[1.0/(3.0*1.255), 1.0/(3.0*0.211)], [1.0/(3.0*1.255), 1.0/(3.0*0.211)], [1.0/(3.0*1.255), 1.0/(3.0*0.211)]])
fuel1bin.set_nu_sigma_f([[0.004602, 0.1091], [0.004602, 0.1091], [0.004602, 0.1091]])
fuel1bin.set_sigma_f([[0.004602/2.43, 0.1091/2.43], [0.004602/2.43, 0.1091/2.43], [0.004602/2.43, 0.1091/2.43]])
fuel1bin.set_chi([[1.0, 0.0], [1.0, 0.0], [1.0, 0.0]])
fuel1bin.set_sigma_s([[0.0, 0.02533, 0.0, 0.0], [0.0, 0.02533, 0.0, 0.0], [0.0, 0.02533, 0.0, 0.0]])
fuel1bin.set_doppler_coefficients([0.0, 0.0])
fuel1bin.set_energy_per_fission(3.204e-11)
fuel1bin.set_velocity([[3.e7, 3.e5], [3.e7, 3.e5], [3.e7, 3.e5]])
fuel1bin.set_temperature_conversion_factor(3.83e-11)



mesh.set_material(fuel1bin, 0)

# Create and initialize the amplitude mesh
amp_mesh = rk.mesh.AmpMesh(name='amp mesh', width=1.0, height=1.0, num_x=1, num_y=1)
amp_mesh.set_num_amp_energy_groups(2)
amp_mesh.set_num_shape_energy_groups(2)
amp_mesh.set_num_delayed_groups(2)
amp_mesh.set_buckling([1.e-4, 1.e-4])
amp_mesh.set_delayed_fractions([0.0054, 0.001087])
amp_mesh.set_decay_constants([0.00654, 1.35])
amp_mesh.initialize()
amp_mesh.set_shape_mesh(mesh)
amp_mesh.set_shape_groups([[0], [1]])
mesh.set_amp_mesh(amp_mesh)
mesh.set_amp_groups([0, 1])

# Solve diffusion problem
solver = rk.solver.CmfdSolver(mesh, amp_mesh)
solver.set_num_threads(1)

transient = rk.transient.Transient()
transient.set_clock(rk.clock.Clock())
transient.set_shape_mesh(mesh)
transient.set_amp_mesh(amp_mesh)
transient.set_solver(solver)
transient.set_initial_power(1.e-6)
transient.compute_initial_shape()
for i in xrange(100):
    transient.take_outer_step()

#rk.plotter.plot_precursor_conc(amp_mesh, name='amp-precursor-conc')
#rk.plotter.plot_flux(amp_mesh, name='amp-flux')
#rk.plotter.plot_flux(mesh)
#rk.plotter.plot_power(mesh)
#rk.plotter.plot_temperature(mesh)
#rk.plotter.plot_materials(mesh)
#rk.plotter.plot_precursor_conc(mesh)
#k.plotter.plot_sigma_a(mesh)

