from openmc.input.settings import SettingsFile
from openmc.input.tallies import TalliesFile, Tally
from openmc.input.plots import PlotsFile, Plot
from openmc.input.material import MaterialsFile
from openmc.input.opencsg_compatible import create_geometry_xml
import opencsg
from datasets.BEAVRS.materials import openmc_materials
from datasets.BEAVRS.lattices import *
import numpy as np


###############################################################################
###################   Simulation Input File Parameters   ######################
###############################################################################

# Get the appropriate lattice from the lattices module
lattice1 = lattices['2.4% Fuel - 0BA']
lattice2 = lattices['3.1% Fuel - 0BA']

# Discretization of pin cells
fuel_rings = 3
mod_rings = 3
sectors = 4

# Height of the axial slice
slice_height = 10.

# OpenMC simulation parameters
batches = 20
inactive = 10
particles = 1000

# Plotting parameters
pixels = 1000


###############################################################################
######################   Creating Bounding Surfaces   #########################
###############################################################################

print('Creating the bounding Surfaces...')

boundaries = dict()

width = lattice_width * 3.

boundaries['X-Min'] = opencsg.XPlane(x0=-width / 2.)
boundaries['X-Max'] = opencsg.XPlane(x0=width / 2.)
boundaries['Y-Min'] = opencsg.YPlane(y0=-width / 2.)
boundaries['Y-Max'] = opencsg.YPlane(y0=width / 2.)
boundaries['Z-Min'] = opencsg.ZPlane(z0=-slice_height / 2.)
boundaries['Z-Max'] = opencsg.ZPlane(z0=slice_height / 2.)

for index in boundaries.keys():
  boundaries[index].setBoundaryType('reflective')


###############################################################################
#########################   Discretize Pin Cells  #############################
###############################################################################

print('Discretizing the Cells...')

sectormesh = opencsg.SectorMesh(num_sectors=sectors)
radialmesh = opencsg.RadialMesh(num_rings=fuel_rings)

# Loop over all pin types
for pin in pin_cells.keys():

  universe = pin_cells[pin]
  sectormesh.subdivideUniverse(universe=universe)
  cells = universe.getCells()

  for cell_id in cells.keys():
    cell = cells[cell_id]

    if 'Fuel' in cell.getName():
      radialmesh.subdivideCell(cell=cell, universe=universe)


###############################################################################
#####################   Creating Colorset Lattice   ###########################
###############################################################################

cell1 = opencsg.Cell(name=lattice1.getName(), fill=lattice1)
cell2 = opencsg.Cell(name=lattice2.getName(), fill=lattice2)

universe1 = opencsg.Universe(name=cell1.getName())
universe1.addCell(cell1)

universe2 = opencsg.Universe(name=cell2.getName())
universe2.addCell(cell2)

lattice = opencsg.Lattice(name='3x3 Lattice')
lattice.setDimension((3, 3))
lattice.setWidth((lattice_width, lattice_width))
lattice.setLowerLeft((-width/2., -width/2.))
lattice.setUniverses([[universe2, universe1, universe2],
                      [universe1, universe2, universe1],
                      [universe2, universe1, universe2]])


###############################################################################
######################   Creating Root Universe   #############################
###############################################################################

print('Creating the root Universe...')

# Root Cell fills the root Universe which encapsulates the complete Geometry
root_universe = opencsg.Universe(universe_id=0, name='Root Universe')
root_cell = opencsg.Cell(name='Root Cell', fill=lattice)
root_universe.addCell(root_cell)

# Append the bounding surfaces for the Geometry to the root Cell
root_cell.addSurface(surface=boundaries['X-Min'], halfspace=+1)
root_cell.addSurface(surface=boundaries['X-Max'], halfspace=-1)
root_cell.addSurface(surface=boundaries['Y-Min'], halfspace=+1)
root_cell.addSurface(surface=boundaries['Y-Max'], halfspace=-1)
root_cell.addSurface(surface=boundaries['Z-Min'], halfspace=+1)
root_cell.addSurface(surface=boundaries['Z-Max'], halfspace=-1)


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

print('Creating the Geometry...')

geometry = opencsg.Geometry()
geometry.setRootUniverse(root_universe)


###############################################################################
###################   Exporting to OpenMC XML Input Files  ####################
###############################################################################

print('Exporting to OpenMC XML Files...')

# geometry.xml
create_geometry_xml(geometry)

# material.xml
materials_file = MaterialsFile()
materials_file.addMaterials(openmc_materials.values())
materials_file.exportToXML()

# settings.xml
settings_file = SettingsFile()
settings_file.setBatches(batches)
settings_file.setInactive(inactive)
settings_file.setParticles(particles)
settings_file.setStatepointInterval(5)

source = [-width/2., -width/2., -slice_height/2.,
          width/2., width/2., slice_height/2.]
settings_file.setSourceSpace(type='box', params=source)

settings_file.exportToXML()

# plots.xml
plot = Plot(plot_id=1)
plot.setWidth(width=[geometry.getMaxX()-geometry.getMinX(),
                     geometry.getMaxY()-geometry.getMinY()])
plot.setOrigin([0., 0., 0.])
plot.setPixels([pixels, pixels])

plot_file = PlotsFile()
plot_file.addPlot(plot)
plot_file.exportToXML()

# tallies.xml
cells = geometry.getAllCells()
tallies_file = TalliesFile()
scores = ['flux']
bins = np.array([0.0, 0.625, 10000000.])

for cell_id in cells.keys():
  cell = cells[cell_id]

  if cell.getType() == 'material':
    tally = Tally()
    tally.addFilter(type='distribcell', bins=cell_id)
    tally.addFilter(type='energy', bins=bins)

    for score in scores:
      tally.addScore(score=score)

    tallies_file.addTally(tally)

tallies_file.exportToXML()