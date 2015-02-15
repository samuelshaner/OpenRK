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
import numpy


################################################################################
####################   Simulation Input File Parameters   ######################
################################################################################

# OpenMC simulation parameters
batches = 100
inactive = 5
particles = 250000
structures = [1,2,4,8,12,16,25,40,70]

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


###############################################################################
##################   Exporting to OpenMC geometry.xml File  ###################
###############################################################################

openmc_geometry = get_openmc_geometry(geometry)
geometry_file = openmc.GeometryFile()
geometry_file.set_geometry(openmc_geometry)
geometry_file.export_to_xml()


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
universe = universes.values()[0]

cells = geometry._root_universe.getAllCells()
for cell_id, cell in cells.items():
  if cell._type == 'material':
    if cell._fill._name == 'Borated Water':
      new_cells = water_mesh.subdivideCell(cell=cell, universe=universe)
    if cell._fill._name == '1.6% Fuel':
      new_cells = fuel_mesh.subdivideCell(cell=cell, universe=universe)

mesh.subdivideUniverse(universe=universe)
geometry.assignAutoIds()

# Get the Gap's ID
materials = geometry.getAllMaterials()
for material_id, material in materials.items():
  if material._name == 'Gap':
    gap_id = material_id

tally_factory = MicroXSTallyFactory(openmc_geometry)

for i, num_groups in enumerate(structures):
  groups = group_structures['CASMO']['{0}-group'.format(num_groups)]
  tally_factory.createAllMultiGroupXS(groups, domain_type='material')

tally_factory.createTalliesFile()

print('running openmc...')

executor = openmc.Executor()
executor.run_simulation(output=True, mpi_procs=8)


#####################   Parametric Sweep Over Energy Groups ####################

for i, num_groups in enumerate(structures):

  print('testing {0} groups'.format(num_groups))

  groups = group_structures['CASMO']['{0}-group'.format(num_groups)]

  ########################   Extracting Cross-Sections  ########################

  # Initialize handle on the OpenMC summary file
  openmc.reset_auto_ids()
  summary = Summary('summary.h5')

  for j, batch in enumerate(range(inactive+5, batches+1, 1)):

    print('batch {0}...'.format(batch))

    openmc.reset_auto_ids()

    openmoc_geometry = get_openmoc_geometry(geometry)
    cells = openmoc_geometry.getRootUniverse().getAllCells()

    # Initialize handle on the OpenMC statepoint file
    filename = 'statepoint.{0:03}.h5'.format(batch)
    statepoint = StatePoint(filename)

    micro_extractor = MicroXSTallyExtractor(statepoint, summary)
    micro_extractor.extractAllMultiGroupXS(groups, 'material')

    micro_extractor.rebalanceAllScatterMatrices()

    all_materials_xs = micro_extractor._multigroup_xs['material']

    for material_id in all_materials_xs.keys():

      print('initializing material {0} ...'.format(material_id))

      material_xs = all_materials_xs[material_id]
      macro_xs = dict()

      for xs_type in xs_types:
        multigroup_xs = material_xs[xs_type]
        macro_xs[xs_type] = multigroup_xs.getMacroXS()

      #########################  Create OpenMOC Materials  #######################
      openmoc_material = openmoc.Material(material_id)
      openmoc_material.setNumEnergyGroups(num_groups)
      openmoc_material.thisown = 0   # FIXME: Can SWIG do this on its own??

      if material_id != gap_id:
        openmoc_material.setSigmaT(macro_xs['transport'])
      else:
        openmoc_material.setSigmaT(macro_xs['total'])
      openmoc_material.setSigmaA(macro_xs['absorption'])
      openmoc_material.setSigmaF(macro_xs['fission'])
      openmoc_material.setNuSigmaF(macro_xs['nu-fission'])
      openmoc_material.setSigmaS(macro_xs['nu-scatter matrix'].ravel())
      openmoc_material.setChi(macro_xs['chi'])

      ######################  Set Materials for OpenMOC Cells  ###################

      for cell_id, cell in cells.items():
        if cell.getType() == openmoc.MATERIAL:
          cell = openmoc.castCellToCellBasic(cell)
          if cell.getMaterial().getId() == material_id:
            cell.setMaterial(openmoc_material)

    #############################   Running OpenMOC  ###########################

    print('running openmoc...')

    openmoc_geometry.initializeFlatSourceRegions()
    track_generator = openmoc.TrackGenerator(openmoc_geometry, 128, 0.05)
    track_generator.generateTracks()

    solver = openmoc.CPUSolver(openmoc_geometry, track_generator)
    solver.setSourceConvergenceThreshold(1E-6)
    solver.convergeSource(1000)
    solver.printTimerReport()

    store_simulation_state(solver, use_hdf5=True,
                           filename='sim-state-{0}-group'.format(num_groups),
                           note='batch-{0}'.format(batch))

    print('stored simulation state...')

    kinf[i,j] = solver.getKeff()


###############################   Plot k-inf Error  ############################

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
plt.xlim((20,100))
plt.legend(legend)
plt.grid()
plt.xlim(-400,400)
plt.savefig('k-inf-err.png')