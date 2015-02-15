import opencg, openmc, openmoc
from openmoc.process import store_simulation_state
from openmc.statepoint import StatePoint
from openmc.summary import Summary
from datasets.energy_groups import group_structures
from geometry import geometry
from infermc.build import MicroXSTallyFactory, xs_types
from infermc.process import MicroXSTallyExtractor
from openmc.opencg_compatible import get_openmc_geometry
from openmoc.compatible.opencg_compatible import get_openmoc_geometry
import matplotlib.pyplot as plt
import numpy, h5py


################################################################################
####################   Simulation Input File Parameters   ######################
################################################################################

# OpenMC simulation parameters
batches = 20 #100
inactive = 5
particles = 50000
structures = [2,4,8,12,16]#,25,40,70]

# Initialize array to contain all data
kinf = numpy.zeros((len(structures), batches-inactive-4), dtype=numpy.float64)

################################################################################
##################   Exporting to OpenMC settings.xml File  ####################
################################################################################

settings_file = openmc.SettingsFile()
settings_file.set_batches(batches)
settings_file.set_inactive(inactive)
settings_file.set_particles(particles)
settings_file.set_statepoint_interval(1)
settings_file.set_ptables(True)
settings_file.set_output({'tallies': False, 'summary': True})
source_bounds = [geometry.getMinX(), geometry.getMinY(), geometry.getMinZ(), \
                 geometry.getMaxX(), geometry.getMaxY(), geometry.getMaxZ()]
settings_file.set_source_space('box', source_bounds)
settings_file.export_to_xml()


################################################################################
########################   Creating OpenMOC Geometry  ##########################
################################################################################

water_mesh = opencg.RadialMesh()
water_mesh.setNumRings(3)
water_mesh.setMaxRadius(0.7)
water_mesh.setMinRadius(0.4572)
water_mesh.setWithOuter(True)

fuel_mesh = opencg.RadialMesh()
fuel_mesh.setNumRings(5)
fuel_mesh.setMaxRadius(0.39218)
fuel_mesh.setMinRadius(0.)
fuel_mesh.setWithOuter(False)

mesh = opencg.SectorMesh(num_sectors=8)

universes = geometry.getAllMaterialUniverses()

cells = geometry._root_universe.getAllCells()
for cell_id, cell in cells.items():
  if cell._type == 'material':
    if cell._fill._id == 10003 and cell._name == 'Moderator':
      new_cells = water_mesh.subdivideCell(cell=cell, universe=universes[10000])
    if cell._fill._id == 10003 and cell._name == 'Outer Water':
      new_cells = water_mesh.subdivideCell(cell=cell, universe=universes[10001])
    if cell._fill._id == 10000:
      new_cells = fuel_mesh.subdivideCell(cell=cell, universe=universes[10000])

mesh.subdivideUniverse(universe=universes[10000])
mesh.subdivideUniverse(universe=universes[10001])
mesh.subdivideUniverse(universe=universes[10002])

openmoc.reset_auto_ids()
openmoc_geometry = get_openmoc_geometry(geometry)
cells = openmoc_geometry.getRootUniverse().getAllCells()

#####################   Parametric Sweep Over Energy Groups ####################

for i, num_groups in enumerate(structures):

  print('testing {0} groups'.format(num_groups))

  groups = group_structures['CASMO']['{0}-group'.format(num_groups)]

  '''

  ##################   Exporting to OpenMC tallies.xml File  ###################

  openmc_geometry = get_openmc_geometry(geometry)
  tally_factory = MicroXSTallyFactory(openmc_geometry)

  tally_factory.createAllMultiGroupXS(groups, domain_type='material')
  tally_factory.createTalliesFile()


  ###############################   Running OpenMC  ############################

  print('running openmc...')

  executor = openmc.Executor()
  executor.run_simulation(output=True, mpi_procs=12)
  '''

  ########################   Extracting Cross-Sections  ########################

  # Initialize handle on the OpenMC summary file
  openmc.reset_auto_ids()
  summary = Summary('summary.h5')

  for j, batch in enumerate(range(inactive+5, batches+1, 1)):

    print('batch {0}...'.format(batch))

    '''
    openmc.reset_auto_ids()

    # Initialize handle on the OpenMC statepoint file
    filename = 'statepoint.{0:02}.h5'.format(batch)
    statepoint = StatePoint('statepoint.{0:02}.h5'.format(batch))

    micro_extractor = MicroXSTallyExtractor(statepoint, summary)
    micro_extractor.extractAllMultiGroupXS(groups, 'material')

    materials = summary.openmc_geometry.get_all_materials()
    '''

    # DUMP-TO-FILE and PRINT XS
    filename = 'mgxs-batch-{0}-groups-{1}'.format(batch, num_groups)

    '''
    for material in materials:
      for xs_type in xs_types:
        xs = micro_extractor._multigroup_xs['material'][material._id][xs_type]
        xs.exportResults(filename=filename)

    statepoint.close()
    del statepoint
    '''

    ###################   Injecting Cross-Sections into OpenMOC  #################

    print('opening {0} ...'.format(filename))

    f = h5py.File('multigroupxs/{0}.h5'.format(filename), 'r')

    mats = f['material']

    for mat_key in mats.keys():

      print('initializing material {0} ...'.format(mat_key))

      material_id = int(mat_key.split(' ')[1])
      material_group = mats[mat_key]

      openmoc_material = openmoc.Material(material_id)
      openmoc_material.setNumEnergyGroups(num_groups)
      openmoc_material.thisown = 0   # FIXME: Can SWIG do this on its own??

      macro_xs = dict()
      macro_xs['total'] = numpy.zeros(num_groups)
      macro_xs['transport'] = numpy.zeros(num_groups)
      macro_xs['scatter matrix'] = numpy.zeros((num_groups, num_groups))
      macro_xs['absorption'] = numpy.zeros(num_groups)
      macro_xs['fission'] = numpy.zeros(num_groups)
      macro_xs['nu-fission'] = numpy.zeros(num_groups)
      macro_xs['chi'] = numpy.zeros(num_groups)

      #####################  Calculate Macro Cross-Sections  #####################

      for nuclide in material_group.keys():
        nuclide_group = material_group[nuclide]
        density = nuclide_group['density'][0]

        if nuclide == 'total':
          continue

        for rxn_type in nuclide_group.keys():
          if rxn_type in macro_xs.keys():
            micro_xs = nuclide_group[rxn_type]['average'][...]
            macro_xs[rxn_type] += micro_xs * density

      # Manually set transport to absorption + scattering
      macro_xs['transport'] = macro_xs['absorption'] + \
                              numpy.sum(macro_xs['scatter matrix'], axis=1)

      #########################  Create OpenMOC Materials  #######################

      if material_id != 10004 and material_id != 10006:
        openmoc_material.setSigmaT(macro_xs['transport'])
      else:
        openmoc_material.setSigmaT(macro_xs['total'])

      openmoc_material.setSigmaA(macro_xs['absorption'])
      openmoc_material.setSigmaF(macro_xs['fission'])
      openmoc_material.setNuSigmaF(macro_xs['nu-fission'])
      openmoc_material.setSigmaS(macro_xs['scatter matrix'].ravel())
      openmoc_material.setChi(macro_xs['chi'])

      ######################  Set Materials for OpenMOC Cells  ###################

      old_material = None

      for cell_id, cell in cells.items():
        if cell.getType() == openmoc.MATERIAL:
          cell = openmoc.castCellToCellBasic(cell)
          if cell.getMaterial().getId() == material_id:
            old_material = cell.getMaterial()
            cell.setMaterial(openmoc_material)

      # Delete the old Materials to save on memory
      del old_material

    # Close the HDF5 multi-group cross-section file
    f.close()

    #############################   Running OpenMOC  ###########################

    print('running openmoc...')

    cmfd = openmoc.Cmfd()
    cmfd.setLatticeStructure(17,17)
    openmoc_geometry.setCmfd(cmfd)

    openmoc_geometry.initializeFlatSourceRegions()
    track_generator = openmoc.TrackGenerator(openmoc_geometry, 128, 0.05)
    track_generator.generateTracks()

    import openmoc.cuda
    openmoc.cuda.attach_gpu(1)
    solver = openmoc.cuda.GPUSolver(openmoc_geometry, track_generator)
    solver.setSourceConvergenceThreshold(1E-5)
    solver.convergeSource(1000)
    solver.printTimerReport()

#    import openmoc.cuda
#    openmoc.cuda.attach_gpu(1)
#    solver = openmoc.cuda.GPUSolver(openmoc_geometry, track_generator)
#    solver.setSourceConvergenceThreshold(1E-6)
#    solver.convergeSource(1000)
#    solver.printTimerReport()

    store_simulation_state(solver, use_hdf5=True,
                           fission_rates=True,
                           filename='sim-state-{0}-group-2'.format(num_groups),
                           note='batch-{0}'.format(batch))

    print('stored simulation state...')

    kinf[i,j] = solver.getKeff()

    del track_generator, solver

###############################   Plot k-inf Error  ############################

for i, num_groups in enumerate(structures):

  filename = 'simulation-states/sim-state-{0}-group.h5'.format(num_groups)
  f = h5py.File(filename, 'r')
  date_key = f.keys()[0]
  date_group = f[date_key]

  for time in date_group.keys():
    time_group = date_group[time]
    note = time_group.attrs['note']
    batch = int(note.split('-')[1])
    keff = time_group['keff'][...]
    kinf[i,batch-10] = keff

  f.close()


kinf_ref = numpy.zeros(batches-inactive-4)
kinf_std_dev = numpy.zeros(batches-inactive-4)

for i, batch in enumerate(range(inactive+5, batches+1, 1)):
  statepoint = StatePoint('statepoint.{0:03}.h5'.format(batch))
  statepoint.read_results()
  kinf_ref[i] = statepoint.k_combined[0]
  kinf_std_dev[:] = statepoint.k_combined[1]
  statepoint.close()
  del statepoint

print('reference openmc kinf: {0}'.format(kinf_ref[-1]))

kinf_std_dev *= 1E5

fix, ax = plt.subplots(1)
legend = list()
batches = numpy.arange(inactive+5, batches+1)

plt.plot(batches, (kinf_ref-kinf_ref[-1])*1E5, linewidth=2)
legend.append('openmc')

for i, num_groups in enumerate(structures):
  print('plotting {0} groups'.format(num_groups))
  kinf_err = (kinf[i,:] - kinf_ref[-1]) * 1E5
  plt.plot(batches, kinf_err, linewidth=2)
  legend.append('{0}-group'.format(num_groups))

ax.fill_between(batches, +kinf_std_dev, -kinf_std_dev,
                facecolor='blue', alpha=0.3)

plt.xlabel('Batch #')
plt.ylabel('Error [pcm]')
plt.title('1.6% Enr. k-inf Error')
plt.legend(legend)
plt.grid()
plt.savefig('k-inf-err.png')