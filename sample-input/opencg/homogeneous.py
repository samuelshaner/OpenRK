import opencg, openmc, openrk, infermc
from openmc.statepoint import StatePoint
from openmc.summary import Summary
from geometry import geometry, mat
from infermc.build import MicroXSTallyFactory, xs_types
from infermc.process import MicroXSTallyExtractor
from openmc.opencg_compatible import get_openmc_geometry, get_openmc_mesh
#from openrk.compatible import get_openrk_currents, get_openrk_amp_mesh, get_openrk_shape_mesh
import matplotlib.pyplot as plt
import numpy


################################################################################
####################   Simulation Input File Parameters   ######################
################################################################################

# OpenMC simulation parameters
batches = 30
inactive = 10
particles = 1000

# Transient simulation parameters
num_x = 1
num_y = 1
num_z = 1
num_groups = 2
mesh_level = 0

###############################################################################
##################   Exporting to OpenMC materials.xml File  ##################
###############################################################################

# OpenMC/OpenCG - set openmc materials
materials_file = openmc.MaterialsFile()
materials_file.set_default_xs('71c')
materials_file.add_material(mat)
materials_file.export_to_xml()

################################################################################
##################   Exporting to OpenMC settings.xml File  ####################
################################################################################

# OpenMC/OpenCG - set openmc simulation settings
settings_file = openmc.SettingsFile()
settings_file.set_batches(batches)
settings_file.set_inactive(inactive)
settings_file.set_particles(particles)
settings_file.set_statepoint_interval(10)
settings_file.set_ptables(True)
settings_file.set_output({'tallies': True, 'summary': True})
source_bounds = [geometry.getMinX(), geometry.getMinY(), geometry.getMinZ(), \
                 geometry.getMaxX(), geometry.getMaxY(), geometry.getMaxZ()]
settings_file.set_source_space('box', source_bounds)
settings_file.export_to_xml()


###############################################################################
##################   Exporting to OpenMC geometry.xml File  ###################
###############################################################################

# OpenMC/OpenCG - set openmc geometry
openmc_geometry = get_openmc_geometry(geometry)
geometry_file = openmc.GeometryFile()
geometry_file.set_geometry(openmc_geometry)
geometry_file.export_to_xml()

################################################################################
#############################   Create Tallies   ###############################
################################################################################


# InferMC - Create XS Tally Factory and get tally file
xs_tally_factory = MicroXSTallyFactory(openmc_geometry)
tally_file = xs_tally_factory.getTalliesFile()

# OpenCG/OpenMC - Extract mesh from OpenCG geometry
opencg_mesh = geometry.extractMesh(mesh_level)
openmc_mesh = get_openmc_mesh(opencg_mesh)

# OpenMC - Create currents tally
current_tally = openmc.Tally()
mesh_filter = openmc.Filter(type='mesh', bins=openmc_mesh._id)
energy_filter = openmc.Filter(type='energy', bins=[0., 0.1e-1, 20.])
current_tally.add_filter(mesh_filter)
current_tally.add_filter(energy_filter)
current_tally.add_score('current')
tally_file.add_mesh(openmc_mesh)
tally_file.add_tally(current_tally)

# InferMC - Create xs tallies
groups = infermc.EnergyGroups()
groups.group_edges = numpy.array([0., 0.1e-1, 20.])
xs_tally_factory.createAllMultiGroupXS(groups, domain_type='material')
xs_tally_factory.createTalliesFile()

# InferMC - Export tallies file to xml
xs_tally_factory.exportTalliesFile()

################################################################################
##############################   Run OpenMC  ###################################
################################################################################

# OpenMC - Execute OpenMC
print('running openmc...')
executor = openmc.Executor()
executor.run_simulation(output=True, mpi_procs=1)

################################################################################
#########################   Extract XS, Flux, and Current  #####################
################################################################################

# OpenMC - Initialize handle on the OpenMC summary file
openmc.reset_auto_ids()
summary = Summary('summary.h5')
  
# OpenMC - Initialize handle on the OpenMC statepoint file
filename = 'statepoint.30.h5'
statepoint = StatePoint(filename)
statepoint.read_results()

# OpenMC/OpenRK - Extract the surface currents
current_tally = statepoint._tallies[10000]
#openrk_currents = get_openrk_currents(current_tally, [1,1,1], 2)

# OpenMC/OpenRK - Extract the surface currents
current_tally = statepoint._tallies[10000]
#openrk_currents = get_openrk_currents(current_tally, [1,1,1], 2)

# InferMC - Extact the macroscopic multigroup xs
micro_extractor = MicroXSTallyExtractor(statepoint, summary)
micro_extractor.extractAllMultiGroupXS(groups, 'material')
micro_extractor.rebalanceAllScatterMatrices()
all_materials_xs = micro_extractor._multigroup_xs['material']
material_id = all_materials_xs.keys()[0]
material_xs = all_materials_xs[material_id]
macro_xs = dict()

for xs_type in xs_types:
  multigroup_xs = material_xs[xs_type]
  macro_xs[xs_type] = multigroup_xs.getMacroXS()

################################################################################
####################   Set up the OpenRK Transient Materials  ##################
################################################################################

# OpenRK - Create openrk material
openrk_material = openrk.Material(material_id)
num_groups = 2
openrk_material.setNumEnergyGroups(num_groups)
openrk_material.setNumDelayedGroups(num_groups)
openrk_material.setSigmaT(macro_xs['transport'])
openrk_material.setSigmaA(macro_xs['absorption'])
openrk_material.setDifCoef([1.0/(3.0*macro_xs['transport'][0]), 1.0/(3.0*macro_xs['transport'][1])]) 
openrk_material.setSigmaF(macro_xs['fission'])
openrk_material.setNuSigmaF(macro_xs['nu-fission'])
openrk_material.setSigmaS(macro_xs['nu-scatter matrix'].ravel())
openrk_material.setChi(macro_xs['chi'])
openrk_material.setVelocity([3.e7, 3.e5])
openrk_material.setTemperatureConversionFactor(3.83e-11)
openrk_material.setEnergyPerFission(3.204e-11)

################################################################################
####################   Set up the OpenRK Transient Meshes  #####################
################################################################################



# OpenRK - Create unstructured OpenRK shape mesh
shape_mesh = get_openrk_shape_mesh(opencg_mesh)
shape_mesh.setNumAmpEnergyGroups(2)
shape_mesh.setNumShapeEnergyGroups(2)
shape_mesh.setNumDelayedGroups(2)
shape_mesh.setDelayedFractions([0.0054, 0.001087])
shape_mesh.setDecayConstants([0.00654, 1.35])
shape_mesh.initialize()
shape_mesh.setTemperature(300)
shape_mesh.setMaterial(openrk_material, 0)
shape_mesh.setFlux()

# OpenRK - Create OpenRK amplitude mesh
amp_mesh = get_openrk_amp_mesh(opencg_mesh)
amp_mesh.setNumAmpEnergyGroups(2)
amp_mesh.setNumShapeEnergyGroups(2)
amp_mesh.setNumDelayedGroups(2)
amp_mesh.setDelayedFractions([0.0054, 0.001087])
amp_mesh.setDecayConstants([0.00654, 1.35])
amp_mesh.initialize()
amp_mesh.setShapeMesh(shape_mesh)
amp_mesh.setGroupStructure([0, 1, 2])
amp_mesh.setCurrents(openrk_currents)
shape_mesh.setAmpMesh(amp_mesh)
shape_mesh.setGroupStructure([0, 1])

# OpenRK - Create diffusion solver
solver = openrk.Solver(shape_mesh, amp_mesh)

# OpenRK - Create transient object
transient = openrk.Transient()
clock = openrk.Clock(dt_inner=1.e-3, dt_outer=1.e-1)
transient.setClock(clock)
transient.setShapeMesh(shape_mesh)
transient.setAmpMesh(amp_mesh)
transient.setSolver(solver)
transient.setInitialPower(1.e-6)

# OpenRK - Use the xs from openmc to solve a steady-state diffusion problem
transient.computeInitialShape()

# Solve the transient problem with constant cross sections to check
# that it holds steady state
# while not clock.isTransientComplete():
#     propagate amplitude forward with openmc
#
#     recompute openmc shape
#     transfer openmc parameters to openrk
#     

