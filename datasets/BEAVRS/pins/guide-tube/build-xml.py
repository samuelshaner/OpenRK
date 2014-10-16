from datasets.energy_groups import group_structures
from openmc import *
from datasets.BEAVRS.materials import openmc_materials
from geometry import geometry
from infermc.build import *


###############################################################################
###################   Simulation Input File Parameters   ######################
###############################################################################

# OpenMC simulation parameters
batches = 15
inactive = 5
particles = 1000


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
materials_file.setDefaultXS('70c')
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
settings_file.setPTables(True)
settings_file.setOutput({'tallies': False})
settings_file.setSourceSpace('box', geometry.getBounds())
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

groups = group_structures['CASMO']['8-group']

tally_factory.createAllXS(groups, domain_type='distribcell')
tally_factory.createAllXS(groups, domain_type='material')
tally_factory.createAllXS(groups, domain_type='cell')
tally_factory.createAllXS(groups, domain_type='universe')

tally_factory.createTalliesFile()
