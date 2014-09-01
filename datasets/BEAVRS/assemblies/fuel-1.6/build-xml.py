from datasets.energy_groups import group_structures
from openmc import *
from datasets.BEAVRS.materials import openmc_materials
from geometry import geometry
from infermc.build import *
import opencsg.plotter as plotter


###############################################################################
###################   Simulation Input File Parameters   ######################
###############################################################################

# OpenMC simulation parameters
batches = 30
inactive = 5
particles = 100


###############################################################################
##################   Exporting to OpenMC geometry.xml File  ###################
###############################################################################

# Create geometry.xml
openmc_geometry = get_openmc_geometry(geometry)
geometry_file = openmc.GeometryFile()
geometry_file.setGeometry(openmc_geometry)
geometry_file.exportToXML()


###############################################################################
##################   Exporting to OpenMC materials.xml File  ##################
###############################################################################

materials_file = MaterialsFile()
materials_file.addMaterials(openmc_materials.values())
materials_file.exportToXML()


###############################################################################
##################   Exporting to OpenMC settings.xml File  ###################
###############################################################################

settings_file = SettingsFile()
settings_file.setBatches(batches)
settings_file.setInactive(inactive)
settings_file.setParticles(particles)
settings_file.setStatepointInterval(5)
settings_file.setOutput({'tallies': False})
source_bounds = [geometry.getMinX(), geometry.getMinY(), geometry.getMinZ(), \
                 geometry.getMaxX(), geometry.getMaxY(), geometry.getMaxZ()]
settings_file.setSourceSpace('box', source_bounds)
settings_file.exportToXML()


###############################################################################
##################   Exporting to OpenMC plots.xml File  ######################
###############################################################################

plot = Plot(plot_id=1)
plot.setWidth(width=[geometry.getMaxX()-geometry.getMinX(),
                     geometry.getMaxY()-geometry.getMinY()])
plot_file = PlotsFile()
plot_file.addPlot(plot)
plot_file.exportToXML()


###############################################################################
##################   Exporting to OpenMC tallies.xml File  ####################
###############################################################################

tally_factory = XSTallyFactory(openmc_geometry)
micro_tally_factory = MicroXSTallyFactory(openmc_geometry)

groups = group_structures['CASMO']['2-group']

tally_factory.createAllMultiGroupXS(groups, domain_type='distribcell')
tally_factory.createAllMultiGroupXS(groups, domain_type='material')
tally_factory.createAllMultiGroupXS(groups, domain_type='cell')
tally_factory.createAllMultiGroupXS(groups, domain_type='universe')
tally_factory.createTalliesFile()

'''
micro_tally_factory.createAllMultiGroupXS(groups, domain_type='distribcell')
micro_tally_factory.createAllMultiGroupXS(groups, domain_type='material')
micro_tally_factory.createAllMultiGroupXS(groups, domain_type='cell')
micro_tally_factory.createAllMultiGroupXS(groups, domain_type='universe')
micro_tally_factory.createTalliesFile()
'''


###############################################################################
#########################   Plotting the Geometry  ############################
###############################################################################

#plotter.plot_neighbor_cells(geometry)
#plotter.plot_neighbor_cells(geometry, unique=True)
#plotter.plot_regions(geometry)
#plotter.plot_materials(geometry)
#plotter.plot_cells(geometry)
