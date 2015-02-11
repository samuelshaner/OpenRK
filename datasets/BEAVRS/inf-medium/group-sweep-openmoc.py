import openmc, openmoc
from openmoc.process import store_simulation_state
from openmc.statepoint import StatePoint
from openmc.summary import Summary
from datasets.energy_groups import group_structures
from infermc.build import MicroXSTallyFactory, xs_types
from infermc.process import MicroXSTallyExtractor
from openmc.opencg_compatible import get_opencg_geometry
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

###############################################################################
#                 Exporting to OpenMC materials.xml File
###############################################################################

# Instantiate some Nuclides
h1 = openmc.Nuclide('H-1')
o16 = openmc.Nuclide('O-16')
u235 = openmc.Nuclide('U-235')
u238 = openmc.Nuclide('U-238')
zr90 = openmc.Nuclide('Zr-90')

# Instantiate a Material and register the Nuclides
inf_medium = openmc.Material(name='moderator')
inf_medium.set_density('g/cc', 5.)
inf_medium.add_nuclide(h1,  0.028999667)
inf_medium.add_nuclide(o16, 0.01450188)
inf_medium.add_nuclide(u235, 0.000114142)
inf_medium.add_nuclide(u238, 0.006886019)
inf_medium.add_nuclide(zr90, 0.002116053)

# Instantiate a MaterialsFile, register all Materials, and export to XML
materials_file = openmc.MaterialsFile()
materials_file.set_default_xs('71c')
materials_file.add_material(inf_medium)
materials_file.export_to_xml()


###############################################################################
#                 Exporting to OpenMC geometry.xml File
###############################################################################

# Instantiate boundary Planes
left = openmc.XPlane(bc_type='reflective', x0=-0.63)
right = openmc.XPlane(bc_type='reflective', x0=0.63)
bottom = openmc.YPlane(bc_type='reflective', y0=-0.63)
top = openmc.YPlane(bc_type='reflective', y0=0.63)

# Instantiate Cells
cell = openmc.Cell(cell_id=1, name='cell')

# Register Surfaces with Cells
cell.add_surface(surface=left, halfspace=+1)
cell.add_surface(surface=right, halfspace=-1)
cell.add_surface(surface=bottom, halfspace=+1)
cell.add_surface(surface=top, halfspace=-1)

# Register Material with Cells
cell.set_fill(inf_medium)

# Instantiate Universe
root = openmc.Universe(universe_id=0, name='root universe')

# Register Cell with Universe
root.add_cell(cell)

# Instantiate a Geometry and register the root Universe
openmc_geometry = openmc.Geometry()
openmc_geometry.set_root_universe(root)

opencg_geometry = get_opencg_geometry(openmc_geometry)
geometry_file = openmc.GeometryFile()
geometry_file.set_geometry(openmc_geometry)
geometry_file.export_to_xml()

#opencg_geometry.assignAutoIds()


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
source_bounds = [-0.63, -0.63, -0.63, +0.63, +0.63, +0.63]
settings_file.set_source_space('box', source_bounds)
settings_file.export_to_xml()

tally_factory = MicroXSTallyFactory(openmc_geometry)

for i, num_groups in enumerate(structures):
  groups = group_structures['CASMO']['{0}-group'.format(num_groups)]
  tally_factory.createAllMultiGroupXS(groups, domain_type='material')

tally_factory.createTalliesFile()

###############################   Running OpenMC  ############################

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

    openmoc_geometry = get_openmoc_geometry(opencg_geometry)
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

      openmoc_material.setSigmaT(macro_xs['transport'])
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
    solver.setSourceConvergenceThreshold(1E-7)
    solver.setNumThreads(4)
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
plt.ylim(-400, 400)
plt.savefig('k-inf-err-nu-scatt.png')
