from openmc.input.settings import SettingsFile
from openmc.input.tallies import TalliesFile, Tally
from openmc.input.plots import PlotsFile, Plot
from openmc.input.material import MaterialsFile
from openmc.input.opencsg_compatible import create_geometry_xml
import opencsg
from datasets.BEAVRS.materials import openmc_materials, opencsg_materials
from datasets.BEAVRS.lattices import *
import numpy as np


# NOTE - This is a quarter core BEAVRS geometry.


###############################################################################
###################   Simulation Input File Parameters   ######################
###############################################################################

# Get the Lattice IDs from the lattices module
lat1 = lattices['1.6% Fuel - 0BA']
lat2 = lattices['2.4% Fuel - 0BA']
lat3 = lattices['2.4% Fuel - 12BA']
lat4 = lattices['2.4% Fuel - 16BA']
lat5 = lattices['3.1% Fuel - 0BA']
lat6 = lattices['3.1% Fuel - 6BA']
lat7 = lattices['3.1% Fuel - 15BA']
lat8 = lattices['3.1% Fuel - 16BA']
lat9 = lattices['3.1% Fuel - 20BA']

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
pixels = 2000


###############################################################################
######################   Creating Bounding Surfaces   #########################
###############################################################################

print('Creating the bounding Surfaces...')

boundaries = dict()

width = lattice_width * 10.

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

cell1 = opencsg.Cell(name=lat1.getName(), fill=lat1)
cell2 = opencsg.Cell(name=lat2.getName(), fill=lat2)
cell3 = opencsg.Cell(name=lat3.getName(), fill=lat3)
cell4 = opencsg.Cell(name=lat4.getName(), fill=lat4)
cell5 = opencsg.Cell(name=lat5.getName(), fill=lat5)
cell6 = opencsg.Cell(name=lat6.getName(), fill=lat6)
cell7 = opencsg.Cell(name=lat7.getName(), fill=lat7)
cell8 = opencsg.Cell(name=lat8.getName(), fill=lat8)
cell9 = opencsg.Cell(name=lat9.getName(), fill=lat9)

# Create a pure water Cell
water = opencsg_materials['Borated Water']
cell10 = opencsg.Cell(name='Water cell', fill=water)

u1 = opencsg.Universe(name=cell1.getName())
u1.addCell(cell1)

u2 = opencsg.Universe(name=cell2.getName())
u2.addCell(cell2)

u3 = opencsg.Universe(name=cell3.getName())
u3.addCell(cell3)

u4 = opencsg.Universe(name=cell4.getName())
u4.addCell(cell4)

u5 = opencsg.Universe(name=cell5.getName())
u5.addCell(cell5)

u6 = opencsg.Universe(name=cell6.getName())
u6.addCell(cell6)

u7 = opencsg.Universe(name=cell7.getName())
u7.addCell(cell7)

u8 = opencsg.Universe(name=cell8.getName())
u8.addCell(cell8)

u9 = opencsg.Universe(name=cell9.getName())
u9.addCell(cell9)

u10 = opencsg.Universe(name=cell10.getName())
u10.addCell(cell10)

lattice = Lattice(name='Quarter Core')
lattice.setDimension((10, 10))
lattice.setWidth((lattice_width, lattice_width))
lattice.setLowerLeft((-width/2., -width/2.))
lattice.setUniverses([[u1, u4, u1, u3, u1, u4, u1, u6, u10, u10],
                      [u4, u1, u3, u1, u3, u1, u9, u5, u10, u10],
                      [u1, u3, u1, u3, u1, u4, u1, u6, u10, u10],
                      [u3, u1, u3, u1, u4, u1, u8, u5, u10, u10],
                      [u1, u3, u1, u4, u2, u4, u5, u10, u10, u10],
                      [u4, u1, u4, u1, u4, u7, u5, u10, u10, u10],
                      [u1, u9, u1, u8, u5, u5, u10, u10, u10, u10],
                      [u6, u5, u6, u5, u10, u10, u10, u10, u10, u10],
                      [u10, u10, u10, u10, u10, u10, u10, u10, u10, u10],
                      [u10, u10, u10, u10, u10, u10, u10, u10, u10, u10]])


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

#from statepoint import StatePoint
#import openmoc.compatible.plotter.plotter as plot

#sp = StatePoint('statepoint.20.h5')
#sp.read_results()

#plot.plot_fluxes(geometry, sp, energies=[0,1], gridsize=250)