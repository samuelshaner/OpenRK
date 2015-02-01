import openrk as rk

mesh = rk.mesh.StructuredShapeMesh(name='shape mesh', width=165.0, height=165.0, num_x=11, num_y=11)
mesh.set_boundary(2, 1)
mesh.set_boundary(3, 1)
mesh.initialize_cells()
mesh.set_num_amp_energy_groups(1)
mesh.set_num_shape_energy_groups(2)
mesh.set_num_delayed_groups(2)

# create fuel 1 blade in
fuel1bin = rk.material.TransientMaterial(name='fuel 1 blade in')
fuel1bin.set_num_energy_groups(2)
fuel1bin.set_num_delayed_groups(2)
fuel1bin.set_sigma_a([0.0083775, 0.1003211])
fuel1bin.set_dif_coef([1.255, 0.211])
fuel1bin.set_sigma_t([1.0/(3.0*1.255), 1.0/(3.0*0.211)])
fuel1bin.set_nu_sigma_f([0.004602, 0.1091])
fuel1bin.set_sigma_f([0.004602, 0.1091])
fuel1bin.set_chi([1.0, 0.0])
fuel1bin.set_sigma_s([0.0, 0.02533, 0.0, 0.0])
fuel1bin.set_energy_per_fission(200.0e6)
fuel1bin.set_velocity([3.e7, 3.e5])
fuel1bin.set_delayed_fractions([0.0054, 0.001087])
fuel1bin.set_decay_constants([0.00654, 1.35])


# create fuel 1 blade out
fuel1bo = rk.material.TransientMaterial(name='fuel 1 blade out')
fuel1bo.set_num_energy_groups(2)
fuel1bo.set_num_delayed_groups(2)
fuel1bo.set_sigma_a([0.0073078, 0.07048902])
fuel1bo.set_dif_coef([1.268, 0.1902])
fuel1bo.set_sigma_t([1.0/(3.0*1.268), 1.0/(3.0*0.1902)])
fuel1bo.set_nu_sigma_f([0.004609, 0.08675])
fuel1bo.set_sigma_f([0.004609, 0.08675])
fuel1bo.set_chi([1.0, 0.0])
fuel1bo.set_sigma_s([0.0, 0.02767, 0.0, 0.0])
fuel1bo.set_energy_per_fission(200.0e6)
fuel1bo.set_velocity([3.e7, 3.e5])
fuel1bo.set_delayed_fractions([0.0054, 0.001087])
fuel1bo.set_decay_constants([0.00654, 1.35])

# create fuel 2 blade in
fuel2bin = rk.material.TransientMaterial(name='fuel 2 blade in')
fuel2bin.set_num_energy_groups(2)
fuel2bin.set_num_delayed_groups(2)
fuel2bin.set_sigma_a([0.0081279, 0.08346091])
fuel2bin.set_dif_coef([1.259, 0.2091])
fuel2bin.set_sigma_t([1.0/(3.0*1.259), 1.0/(3.0*0.2091)])
fuel2bin.set_nu_sigma_f([0.004663, 0.1021])
fuel2bin.set_sigma_f([0.004663, 0.1021])
fuel2bin.set_chi([1.0, 0.0])
fuel2bin.set_sigma_s([0.0, 0.02617, 0.0, 0.0])
fuel2bin.set_energy_per_fission(200.0e6)
fuel2bin.set_velocity([3.e7, 3.e5])
fuel2bin.set_delayed_fractions([0.0054, 0.001087])
fuel2bin.set_decay_constants([0.00654, 1.35])

# create fuel 2 blade out
fuel2bo = rk.material.TransientMaterial(name='fuel 2 blade out')
fuel2bo.set_num_energy_groups(2)
fuel2bo.set_num_delayed_groups(2)
fuel2bo.set_sigma_a([0.0081279, 0.07334491])
fuel2bo.set_dif_coef([1.259, 0.2091])
fuel2bo.set_sigma_t([1.0/(3.0*1.259), 1.0/(3.0*0.2091)])
fuel2bo.set_nu_sigma_f([0.004663, 0.1021])
fuel2bo.set_sigma_f([0.004663, 0.1021])
fuel2bo.set_chi([1.0, 0.0])
fuel2bo.set_sigma_s([0.0, 0.02617, 0.0, 0.0])
fuel2bo.set_energy_per_fission(200.0e6)
fuel2bo.set_velocity([3.e7, 3.e5])
fuel2bo.set_delayed_fractions([0.0054, 0.001087])
fuel2bo.set_decay_constants([0.00654, 1.35])

# create reflector
reflector = rk.material.Material(name='reflector')
reflector.set_num_energy_groups(2)
reflector.set_sigma_a([0.0007291, 0.01912592])
reflector.set_dif_coef([1.257, 0.1592])
reflector.set_sigma_t([1.0/(3.0*1.257), 1.0/(3.0*0.1592)])
reflector.set_nu_sigma_f([0.0, 0.0])
reflector.set_sigma_f([0.0, 0.0])
reflector.set_chi([1.0, 0.0])
reflector.set_sigma_s([0.0, 0.04754, 0.0, 0.0])
reflector.set_energy_per_fission(200.0e6)
reflector.set_velocity([3.e7, 3.e5])

for i in xrange(9, 11):
    for cell in mesh._cells[i*mesh._num_x:(i+1)*mesh._num_x]:
        cell.set_material(reflector)

for i in xrange(9):
    for cell in mesh._cells[i*mesh._num_x+7:(i+1)*mesh._num_x]:
        cell.set_material(reflector)

for i in xrange(7, 9):
    for cell in mesh._cells[i*mesh._num_x:i*mesh._num_x+7]:
        cell.set_material(fuel2bin)

for i in xrange(7):
    for cell in mesh._cells[i*mesh._num_x+7:i*mesh._num_x+9]:
        cell.set_material(fuel2bin)

mesh._cells[6*mesh._num_x].set_material(fuel1bo)
mesh._cells[5*mesh._num_x].set_material(fuel1bo)

for i in xrange(5, 7):
    for cell in mesh._cells[i*mesh._num_x+5:i*mesh._num_x+7]:
        cell.set_material(fuel1bo)

mesh._cells[0].set_material(fuel1bo)
for cell in mesh._cells[5:7]:
    cell.set_material(fuel1bo)

mesh._cells[7*mesh._num_x+7].set_material(fuel2bo)

for i in xrange(5, 7):
    for cell in mesh._cells[i*mesh._num_x+1:i*mesh._num_x+5]:
        cell.set_material(fuel1bin)

for i in xrange(1, 5):
    for cell in mesh._cells[i*mesh._num_x:i*mesh._num_x+7]:
        cell.set_material(fuel1bin)

for cell in mesh._cells[1:5]:
    cell.set_material(fuel1bin)

mesh = mesh.uniform_refine(8)
mesh.uniquify_materials()
mesh.initialize_surfaces()
mesh.initialize_field_variables()
solver = rk.solver.CmfdSolver(mesh)
solver.set_num_threads(4)
solver.solve(1.e-6)

mesh.set_k_eff_0(solver._k_eff)
mesh.compute_initial_precursor_conc()

amp_mesh = rk.mesh.AmpMesh(name='amp mesh', width=165.0, height=165.0, num_x=11, num_y=11)
amp_mesh.set_num_amp_energy_groups(1)
amp_mesh.set_num_shape_energy_groups(2)
amp_mesh.set_num_delayed_groups(2)
amp_mesh.set_boundary(2, 1)
amp_mesh.set_boundary(3, 1)
amp_mesh.initialize_cells()
amp_mesh.set_shape_mesh(mesh)
amp_mesh.set_shape_groups([[0, 1]])

mesh.set_amp_mesh(amp_mesh)
mesh.set_amp_groups([0, 0])

amp_mesh.condense_materials()
rk.plotter.plot_precursor_conc(amp_mesh, name='amp-precursor-conc')
rk.plotter.plot_flux(amp_mesh, name='amp-flux')

mesh.compute_power()
rk.plotter.plot_flux(mesh)
rk.plotter.plot_power(mesh)
rk.plotter.plot_materials(mesh)
rk.plotter.plot_precursor_conc(mesh)