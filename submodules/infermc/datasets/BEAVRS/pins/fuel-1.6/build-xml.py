from openmc import *
from datasets.BEAVRS.materials import openmc_materials
from datasets.energy_groups import group_structures
from geometry import geometry
from infermc.build import *
import opencg.plotter as plotter


###############################################################################
###################   Simulation Input File Parameters   ######################
###############################################################################

# OpenMC simulation parameters
batches = 30
inactive = 10
particles = 1000


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
materials_file.set_default_xs('71c')
materials_file.export_to_xml()


###############################################################################
##################   Exporting to OpenMC settings.xml File  ###################
###############################################################################

settings_file = SettingsFile()
settings_file.set_batches(batches)
settings_file.set_inactive(inactive)
settings_file.set_particles(particles)
settings_file.set_statepoint_interval(5)
settings_file.set_ptables(True)
settings_file.set_output({'tallies': False, 'summary': True})
source_bounds = [geometry.getMinX(), geometry.getMinY(), geometry.getMinZ(), \
                 geometry.getMaxX(), geometry.getMaxY(), geometry.getMaxZ()]
settings_file.set_source_space('box', source_bounds)
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

#micro_tally_factory.createAllMultiGroupXS(groups, domain_type='distribcell')
micro_tally_factory.createAllMultiGroupXS(groups, domain_type='material')
micro_tally_factory.createTalliesFile()


###############################################################################
#########################   Plotting the Geometry  ############################
###############################################################################

#plotter.plot_neighbor_cells(geometry)
#plotter.plot_materials(geometry)