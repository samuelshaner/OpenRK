import openrk as rk
import random

mesh = rk.mesh.StructuredShapeMesh(name='shape mesh', width=165.0, height=165.0, num_x=11, num_y=11)
mesh.set_boundary(2, 1)
mesh.set_boundary(3, 1)
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
fuel1bin.set_doppler_coefficients([3.034e-3, 0.0])
fuel1bin.set_energy_per_fission(3.204e-11)
fuel1bin.set_velocity([[3.e7, 3.e5], [3.e7, 3.e5], [3.e7, 3.e5]])
fuel1bin.set_temperature_conversion_factor(3.83e-11)

# create fuel 1 blade out
fuel1bo = rk.material.FunctionalMaterial(name='fuel 1 blade out')
fuel1bo.set_num_time_steps(3)
fuel1bo.set_num_energy_groups(2)
fuel1bo.set_num_delayed_groups(2)
fuel1bo.set_time_steps([0.0, 2.0, 3.0])
fuel1bo.set_sigma_a([[0.007181, 0.07047], [0.007181, 0.07047], [0.007181, 0.07047]])
fuel1bo.set_dif_coef([[1.268, 0.1902], [1.268, 0.1902], [1.268, 0.1902]])
fuel1bo.set_sigma_t([[1.0/(3.0*1.268), 1.0/(3.0*0.1902)], [1.0/(3.0*1.268), 1.0/(3.0*0.1902)], [1.0/(3.0*1.268), 1.0/(3.0*0.1902)]])
fuel1bo.set_nu_sigma_f([[0.004609, 0.08675], [0.004609, 0.08675], [0.004609, 0.08675]])
fuel1bo.set_sigma_f([[0.004609/2.43, 0.08675/2.43], [0.004609/2.43, 0.08675/2.43], [0.004609/2.43, 0.08675/2.43]])
fuel1bo.set_chi([[1.0, 0.0], [1.0, 0.0], [1.0, 0.0]])
fuel1bo.set_sigma_s([[0.0, 0.02767, 0.0, 0.0], [0.0, 0.02767, 0.0, 0.0], [0.0, 0.02767, 0.0, 0.0]])
fuel1bo.set_doppler_coefficients([3.034e-3, 0.0])
fuel1bo.set_energy_per_fission(3.204e-11)
fuel1bo.set_velocity([[3.e7, 3.e5], [3.e7, 3.e5], [3.e7, 3.e5]])
fuel1bo.set_temperature_conversion_factor(3.83e-11)

# create fuel 2 blade in
fuel2bin = rk.material.FunctionalMaterial(name='fuel 2 blade in')
fuel2bin.set_num_time_steps(3)
fuel2bin.set_num_energy_groups(2)
fuel2bin.set_num_delayed_groups(2)
fuel2bin.set_time_steps([0.0, 2.0, 3.0])
fuel2bin.set_sigma_a([[0.008002, 0.08344], [0.008002, 0.08344], [0.008002, 0.08344]])
fuel2bin.set_dif_coef([[1.259, 0.2091], [1.259, 0.2091], [1.259, 0.2091]])
fuel2bin.set_sigma_t([[1.0/(3.0*1.259), 1.0/(3.0*0.2091)], [1.0/(3.0*1.259), 1.0/(3.0*0.2091)], [1.0/(3.0*1.259), 1.0/(3.0*0.2091)]])
fuel2bin.set_nu_sigma_f([[0.004663, 0.1021], [0.004663, 0.1021], [0.004663, 0.1021]])
fuel2bin.set_sigma_f([[0.004663/2.43, 0.1021/2.43], [0.004663/2.43, 0.1021/2.43], [0.004663/2.43, 0.1021/2.43]])
fuel2bin.set_chi([[1.0, 0.0], [1.0, 0.0], [1.0, 0.0]])
fuel2bin.set_sigma_s([[0.0, 0.02617, 0.0, 0.0], [0.0, 0.02617, 0.0, 0.0], [0.0, 0.02617, 0.0, 0.0]])
fuel2bin.set_doppler_coefficients([3.034e-3, 0.0])
fuel2bin.set_energy_per_fission(3.204e-11)
fuel2bin.set_velocity([[3.e7, 3.e5], [3.e7, 3.e5], [3.e7, 3.e5]])
fuel2bin.set_temperature_conversion_factor(3.83e-11)

# create fuel 2 blade out
fuel2bo = rk.material.FunctionalMaterial(name='fuel 2 blade out')
fuel2bo.set_num_time_steps(3)
fuel2bo.set_num_energy_groups(2)
fuel2bo.set_num_delayed_groups(2)
fuel2bo.set_time_steps([0.0, 2.0, 3.0])
fuel2bo.set_sigma_a([[0.008002, 0.073324], [0.008002, 0.073324], [0.008002, 0.073324]])
fuel2bo.set_dif_coef([[1.259, 0.2091], [1.259, 0.2091], [1.259, 0.2091]])
fuel2bo.set_sigma_t([[1.0/(3.0*1.259), 1.0/(3.0*0.2091)], [1.0/(3.0*1.259), 1.0/(3.0*0.2091)], [1.0/(3.0*1.259), 1.0/(3.0*0.2091)]])
fuel2bo.set_nu_sigma_f([[0.004663, 0.1021], [0.004663, 0.1021], [0.004663, 0.1021]])
fuel2bo.set_sigma_f([[0.004663/2.43, 0.1021/2.43], [0.004663/2.43, 0.1021/2.43], [0.004663/2.43, 0.1021/2.43]])
fuel2bo.set_chi([[1.0, 0.0], [1.0, 0.0], [1.0, 0.0]])
fuel2bo.set_sigma_s([[0.0, 0.02617, 0.0, 0.0], [0.0, 0.02617, 0.0, 0.0], [0.0, 0.02617, 0.0, 0.0]])
fuel2bo.set_doppler_coefficients([3.034e-3, 0.0])
fuel2bo.set_energy_per_fission(3.204e-11)
fuel2bo.set_velocity([[3.e7, 3.e5], [3.e7, 3.e5], [3.e7, 3.e5]])
fuel2bo.set_temperature_conversion_factor(3.83e-11)

# create fuel 2 blade in then out
fuel2bino = rk.material.FunctionalMaterial(name='fuel 2 blade in then out')
fuel2bino.set_num_time_steps(3)
fuel2bino.set_num_energy_groups(2)
fuel2bino.set_num_delayed_groups(2)
fuel2bino.set_time_steps([0.0, 2.0, 3.0])
fuel2bino.set_sigma_a([[0.008002, 0.08344], [0.008002, 0.073324], [0.008002, 0.073324]])
fuel2bino.set_dif_coef([[1.259, 0.2091], [1.259, 0.2091], [1.259, 0.2091]])
fuel2bino.set_sigma_t([[1.0/(3.0*1.259), 1.0/(3.0*0.2091)], [1.0/(3.0*1.259), 1.0/(3.0*0.2091)], [1.0/(3.0*1.259), 1.0/(3.0*0.2091)]])
fuel2bino.set_nu_sigma_f([[0.004663, 0.1021], [0.004663, 0.1021], [0.004663, 0.1021]])
fuel2bino.set_sigma_f([[0.004663/2.43, 0.1021/2.43], [0.004663/2.43, 0.1021/2.43], [0.004663/2.43, 0.1021/2.43]])
fuel2bino.set_chi([[1.0, 0.0], [1.0, 0.0], [1.0, 0.0]])
fuel2bino.set_sigma_s([[0.0, 0.02617, 0.0, 0.0], [0.0, 0.02617, 0.0, 0.0], [0.0, 0.02617, 0.0, 0.0]])
fuel2bino.set_doppler_coefficients([3.034e-3, 0.0])
fuel2bino.set_energy_per_fission(3.204e-11)
fuel2bino.set_velocity([[3.e7, 3.e5], [3.e7, 3.e5], [3.e7, 3.e5]])
fuel2bino.set_temperature_conversion_factor(3.83e-11)

# create reflector
reflector = rk.material.Material(name='reflector')
reflector.set_num_energy_groups(2)
reflector.set_sigma_a([0.000603, 0.01911])
reflector.set_dif_coef([1.257, 0.1592])
reflector.set_sigma_t([1.0/(3.0*1.257), 1.0/(3.0*0.1592)])
reflector.set_nu_sigma_f([0.0, 0.0])
reflector.set_sigma_f([0.0, 0.0])
reflector.set_chi([1.0, 0.0])
reflector.set_sigma_s([0.0, 0.04754, 0.0, 0.0])
reflector.set_velocity([3.e7, 3.e5])

nx = mesh.get_num_x()

for i in xrange(9, 11):
    for cell_id in xrange(i*nx, (i+1)*nx):
        mesh.set_material(reflector, cell_id)

for i in xrange(9):
    for cell_id in xrange(i*nx+7, (i+1)*nx):
        mesh.set_material(reflector, cell_id)

for i in xrange(7, 9):
    for cell_id in xrange(i*nx, i*nx+7):
        mesh.set_material(fuel2bin, cell_id)

for i in xrange(5, 7):
    for cell_id in xrange(i*nx+7, i*nx+9):
        mesh.set_material(fuel2bino, cell_id)

for i in xrange(5):
    for cell_id in xrange(i*nx+7, i*nx+9):
        mesh.set_material(fuel2bin, cell_id)

mesh.set_material(fuel1bo, 6*nx)
mesh.set_material(fuel1bo, 5*nx)

for i in xrange(5, 7):
    for cell_id in xrange(i*nx+5, i*nx+7):
        mesh.set_material(fuel1bo, cell_id)

mesh.set_material(fuel1bo, 0)
for cell_id in xrange(5, 7):
    mesh.set_material(fuel1bo, cell_id)

mesh.set_material(fuel2bo, 7*nx+7)

for i in xrange(5, 7):
    for cell_id in xrange(i*nx+1, i*nx+5):
        mesh.set_material(fuel1bin, cell_id)

for i in xrange(1, 5):
    for cell_id in xrange(i*nx, i*nx+7):
        mesh.set_material(fuel1bin, cell_id)

for cell_id in xrange(1, 5):
    mesh.set_material(fuel1bin, cell_id)

# refine mesh and uniquify materials
mesh = mesh.uniform_refine(3)
mesh.uniquify_materials()

# Create and initialize the amplitude mesh
amp_mesh = rk.mesh.AmpMesh(name='amp mesh', width=165.0, height=165.0, num_x=11, num_y=11)
amp_mesh.set_num_amp_energy_groups(2)
amp_mesh.set_num_shape_energy_groups(2)
amp_mesh.set_num_delayed_groups(2)
amp_mesh.set_optically_thick(True)
amp_mesh.set_boundary(2, 1)
amp_mesh.set_boundary(3, 1)
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
transient.set_clock(rk.clock.Clock(dt_inner=2.5e-2, dt_outer=1.e-1))
transient.set_shape_mesh(mesh)
transient.set_amp_mesh(amp_mesh)
transient.set_solver(solver)
transient.set_initial_power(1.e-6)
transient.compute_initial_shape()

for i in xrange(30):
    transient.take_outer_step()
    rk.plotter.plot_power(mesh, name='mesh-power-{:.4f}s'.format(mesh.get_clock().get_time('CURRENT')))
    rk.plotter.plot_temperature(mesh, name='mesh-temp-{:.4f}s'.format(mesh.get_clock().get_time('CURRENT')))

#rk.plotter.plot_precursor_conc(amp_mesh, name='amp-precursor-conc')
#rk.plotter.plot_flux(amp_mesh, name='amp-flux')
#rk.plotter.plot_flux(mesh)
#rk.plotter.plot_power(mesh)
#rk.plotter.plot_temperature(mesh)
#rk.plotter.plot_materials(mesh)
#rk.plotter.plot_precursor_conc(mesh)
#k.plotter.plot_sigma_a(mesh)

