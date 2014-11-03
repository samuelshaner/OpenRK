from datasets.energy_groups import group_structures
from openmc import *
from datasets.BEAVRS.materials import openmc_materials, nuclides
from geometry import geometry
from infermc.build import *


###############################################################################
###################   Simulation Input File Parameters   ######################
###############################################################################

# OpenMC simulation parameters
batches = 250
inactive = 50
particles = 100000


###############################################################################
##################   Exporting to OpenMC geometry.xml File  ###################
###############################################################################

# Create geometry.xml
openmc_geometry = get_openmc_geometry(geometry)
geometry_file = openmc.GeometryFile()
geometry_file.set_geometry(openmc_geometry)
geometry_file.export_to_xml()


###############################################################################
##################   Exporting to OpenMC materials.xml File  ##################
###############################################################################

materials_file = MaterialsFile()
materials_file.add_materials(openmc_materials.values())
materials_file.export_to_xml()


###############################################################################
##################   Exporting to OpenMC settings.xml File  ###################
###############################################################################

settings_file = SettingsFile()
settings_file.set_batches(batches)
settings_file.set_inactive(inactive)
settings_file.set_particles(particles)
settings_file.set_statepoint_interval(5)
settings_file.set_output({'tallies': False})
settings_file.set_ptables(True)
source_bounds = [geometry.getMinX(), geometry.getMinY(), geometry.getMinZ(), \
                 geometry.getMaxX(), geometry.getMaxY(), geometry.getMaxZ()]
settings_file.set_source_space('box', source_bounds)
settings_file.set_threads(8)
settings_file.export_to_xml()


###############################################################################
##################   Exporting to OpenMC plots.xml File  ######################
###############################################################################

plot = Plot(plot_id=1)
plot.set_width(width=[geometry.getMaxX()-geometry.getMinX(),
                     geometry.getMaxY()-geometry.getMinY()])
plot_file = PlotsFile()
plot_file.add_plot(plot)
plot_file.export_to_xml()


###############################################################################
##################   Exporting to OpenMC tallies.xml File  ####################
###############################################################################

tally_factory = XSTallyFactory(openmc_geometry)
micro_tally_factory = MicroXSTallyFactory(openmc_geometry)

groups = group_structures['CASMO']['2-group']

nuclides = openmc_geometry.get_all_nuclides()
tally_nuclides = [nuclides['H-1'], nuclides['B-10'], nuclides['O-16'],
                  nuclides['U-235'], nuclides['U-238'], nuclides['Zr-90'],
                  nuclides['Fe-56']]

micro_tally_factory.createAllMultiGroupXS(groups, domain_type='distribcell',
                                          nuclides=tally_nuclides)
micro_tally_factory.createAllMultiGroupXS(groups, domain_type='cell',
                                          nuclides=tally_nuclides)

micro_tally_factory.createTalliesFile()
