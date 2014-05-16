from openmc.input.settings import SettingsFile
from openmc.input.tallies import TalliesFile, Tally
from openmc.input.plots import PlotsFile, Plot
from openmc.input.material import MaterialsFile
from openmc.input.opencsg_compatible import create_geometry_xml
from datasets.BEAVRS.materials import openmc_materials
from geometry import geometry, slice_height, pin_pitch
import numpy as np
from infermc.build import *
from sets import Set


###############################################################################
###################   Simulation Input File Parameters   ######################
###############################################################################

# OpenMC simulation parameters
batches = 25
inactive = 10
particles = 10000

# Plotting parameters
pixels = 1000


###############################################################################
##################   Exporting to OpenMC geometry.xml File  ###################
###############################################################################

create_geometry_xml(geometry)


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

source = [-pin_pitch/2., -pin_pitch/2., -slice_height/2.,
          pin_pitch/2., pin_pitch/2., slice_height/2.]
settings_file.setSourceSpace(type='box', params=source)

settings_file.setEntropyDimension([17,17,1])
lower_left = [geometry.getMinX(), geometry.getMinY(), geometry.getMinZ()]
upper_right = [geometry.getMaxX(), geometry.getMaxY(), geometry.getMaxZ()]
settings_file.setEntropyLowerLeft(lower_left)
settings_file.setEntropyUpperRight(upper_right)

settings_file.exportToXML()


###############################################################################
##################   Exporting to OpenMC plots.xml File  ######################
###############################################################################

plot = Plot(plot_id=1)
plot.setWidth(width=[geometry.getMaxX()-geometry.getMinX(),
                     geometry.getMaxY()-geometry.getMinY()])
plot.setOrigin([0., 0., 0.])
plot.setPixels([pixels, pixels])

plot_file = PlotsFile()
plot_file.addPlot(plot)
plot_file.exportToXML()


###############################################################################
##################   Exporting to OpenMC tallies.xml File  ####################
###############################################################################

cells = geometry.getAllCells()
tallies_file = TalliesFile()

energy_groups = EnergyGroups()
energy_groups.setGroupEdges(np.array([0.0, 0.625, 10000000.]))

for cell_id in cells.keys():

  cell = cells[cell_id]

  if cell.getType() == 'material':

    total_xs = TotalXS()
    abs_xs = AbsorptionXS()
    nufiss_xs = NuFissionXS()
    scat_xs = ScatterMatrixXS()
    chi = Chi()

    total_xs.setEnergyGroups(energy_groups)
    abs_xs.setEnergyGroups(energy_groups)
    nufiss_xs.setEnergyGroups(energy_groups)
    scat_xs.setEnergyGroups(energy_groups)
    chi.setEnergyGroups(energy_groups)

    total_xs.createTallies()
    abs_xs.createTallies()
    nufiss_xs.createTallies()
    scat_xs.createTallies()
    chi.createTallies()

    tallies = Set()
    tallies.union_update(Set(total_xs.getTallies()))
    tallies.union_update(Set(abs_xs.getTallies()))
    tallies.union_update(Set(nufiss_xs.getTallies()))
    tallies.union_update(Set(scat_xs.getTallies()))
    tallies.union_update(Set(chi.getTallies()))

    for tally in tallies:
      tally.addFilter(type='distribcell', bins=cell_id)
      tallies_file.addTally(tally)

tallies_file.exportToXML()