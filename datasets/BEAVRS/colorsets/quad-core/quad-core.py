from openmoc.compatible.openmc import *
import openmoc.plotter as plotter
from lattices import *
import materials
import surfaces
import pincells
from input.settings import SettingsFile
from input.tallies import TalliesFile, Tally
from input.plots import PlotsFile, Plot


# NOTE - This is a quarter core BEAVRS geometry.


###############################################################################
###################   Simulation Input File Parameters   ######################
###############################################################################

# Get the Lattice IDs from the lattices module
id1 = lattices['1.6% Fuel - 0BA'].getId()
id2 = lattices['2.4% Fuel - 0BA'].getId()
id3 = lattices['2.4% Fuel - 12BA'].getId()
id4 = lattices['2.4% Fuel - 16BA'].getId()
id5 = lattices['3.1% Fuel - 0BA'].getId()
id6 = lattices['3.1% Fuel - 6BA'].getId()
id7 = lattices['3.1% Fuel - 15BA'].getId()
id8 = lattices['3.1% Fuel - 16BA'].getId()
id9 = lattices['3.1% Fuel - 20BA'].getId()

# Discretization of pin cells
rings = 3
sectors = 4

# Height of the axial slice
slice_height = 10.

# OpenMC simulation parameters
batches = 11
inactive = 10
particles = 10000

# Plotting parameters
pixels = 2000


###############################################################################
######################   Creating Bounding Surfaces   #########################
###############################################################################

log.py_printf('NORMAL', 'Creating the bounding Surfaces...')

boundaries = dict()

width = lattice_width * 10.

boundaries['X-Min'] = XPlane(x=-width / 2.)
boundaries['X-Max'] = XPlane(x=width / 2.)
boundaries['Y-Min'] = YPlane(y=-width / 2.)
boundaries['Y-Max'] = YPlane(y=width / 2.)
boundaries['Z-Min'] = ZPlane(z=-slice_height / 2.)
boundaries['Z-Max'] = ZPlane(z=slice_height / 2.)

for index in boundaries.keys():
  boundaries[index].setBoundaryType(REFLECTIVE)
  surfaces.surfaces[index] = boundaries[index]


###############################################################################
#########################   Discretize Pin Cells  #############################
###############################################################################

# Loop over all pin types
for pin in pincells.pincells.keys():

  # Loop over all cells for this pin type
  for cell in pincells.pincells[pin].keys():

    # Ignore the Universe entry
    if cell == 'Universe':
      continue

    # Set the number of sectors
    pincells.pincells[pin][cell].setNumSectors(sectors)

    # Set the number of rings in the fuel
    if 'Fuel' in cell:
      pincells.pincells[pin][cell].setNumRings(rings)


###############################################################################
#####################   Creating Colorset Lattice   ###########################
###############################################################################

cell1 = CellFill(universe=universe_id(), universe_fill=id1)
cell2 = CellFill(universe=universe_id(), universe_fill=id2)
cell3 = CellFill(universe=universe_id(), universe_fill=id3)
cell4 = CellFill(universe=universe_id(), universe_fill=id4)
cell5 = CellFill(universe=universe_id(), universe_fill=id5)
cell6 = CellFill(universe=universe_id(), universe_fill=id6)
cell7 = CellFill(universe=universe_id(), universe_fill=id7)
cell8 = CellFill(universe=universe_id(), universe_fill=id8)
cell9 = CellFill(universe=universe_id(), universe_fill=id9)

# Create a pure water Cell
water_id = materials.materials['Borated Water'].getId()
cell10 = CellBasic(universe=universe_id(), material=water_id)

id1 = cell1.getUniverseId()
id2 = cell2.getUniverseId()
id3 = cell3.getUniverseId()
id4 = cell4.getUniverseId()
id5 = cell5.getUniverseId()
id6 = cell6.getUniverseId()
id7 = cell7.getUniverseId()
id8 = cell8.getUniverseId()
id9 = cell9.getUniverseId()
id10 = cell10.getUniverseId()

lattice = Lattice(id=universe_id(), width_x=lattice_width, width_y=lattice_width)
lattice.setLatticeCells([[id1, id4, id1, id3, id1, id4, id1, id6, id10, id10],
                         [id4, id1, id3, id1, id3, id1, id9, id5, id10, id10],
                         [id1, id3, id1, id3, id1, id4, id1, id6, id10, id10],
                         [id3, id1, id3, id1, id4, id1, id8, id5, id10, id10],
                         [id1, id3, id1, id4, id2, id4, id5, id10, id10, id10],
                         [id4, id1, id4, id1, id4, id7, id5, id10, id10, id10],
                         [id1, id9, id1, id8, id5, id5, id10, id10, id10, id10],
                         [id6, id5, id6, id5, id10, id10, id10, id10, id10, id10],
                         [id10, id10, id10, id10, id10, id10, id10, id10, id10, id10],
                         [id10, id10, id10, id10, id10, id10, id10, id10, id10, id10]])

lattices['Quarter Core'] = lattice


###############################################################################
######################   Creating Root Universe   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating the root Universe...')

# Root cell encapsulates the full geometry
root = CellFill(universe=0, universe_fill=lattice.getId())
root.addSurface(halfspace=+1, surface=boundaries['X-Min'])
root.addSurface(halfspace=-1, surface=boundaries['X-Max'])
root.addSurface(halfspace=+1, surface=boundaries['Y-Min'])
root.addSurface(halfspace=-1, surface=boundaries['Y-Max'])
root.addSurface(halfspace=+1, surface=boundaries['Z-Min'])
root.addSurface(halfspace=-1, surface=boundaries['Z-Max'])

# Append the root to the list of all Cells
pincells.cells.append(root)


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating the Geometry...')

geometry = Geometry()

# Add all Materials to the Geometry
for material in materials.materials.values():
  geometry.addMaterial(material)

# Add all Surfaces to the Geometry
for surface in surfaces.surfaces.values():
  geometry.addSurface(surface)

for cell in pincells.pincells['1.6% Fuel'].values():
  if not isinstance(cell, int):
    geometry.addCell(cell)

for cell in pincells.pincells['2.4% Fuel'].values():
  if not isinstance(cell, int):
    geometry.addCell(cell)

for cell in pincells.pincells['3.1% Fuel'].values():
  if not isinstance(cell, int):
    geometry.addCell(cell)

for cell in pincells.pincells['Guide Tube'].values():
  if not isinstance(cell, int):
    geometry.addCell(cell)

for cell in pincells.pincells['Instrument Tube'].values():
  if not isinstance(cell, int):
    geometry.addCell(cell)

for cell in pincells.pincells['Burnable Absorber'].values():
  if not isinstance(cell, int):
    geometry.addCell(cell)

geometry.addCell(root)
geometry.addCell(cell1)
geometry.addCell(cell2)
geometry.addCell(cell3)
geometry.addCell(cell4)
geometry.addCell(cell5)
geometry.addCell(cell6)
geometry.addCell(cell7)
geometry.addCell(cell8)
geometry.addCell(cell9)
geometry.addCell(cell10)

# Add all Lattices to the Geometry
for lattice in lattices.values():
  geometry.addLattice(lattice)


###############################################################################
###################   Exporting to OpenMC XML Input Files  ####################
###############################################################################

log.py_printf('NORMAL', 'Exporting to OpenMC XML Files...')

# settings.xml
settings_file = SettingsFile()
settings_file.setBatches(batches)
settings_file.setInactive(inactive)
settings_file.setParticles(particles)
settings_file.setStatepointInterval(5)

source = [-pin_pitch/2., -pin_pitch/2., -slice_height/2.,
          pin_pitch/2., pin_pitch/2., slice_height/2.]
settings_file.setSourceSpace(type='box', params=source)

settings_file.exportToXML()

# plots.xml
plot = Plot(plot_id=1)
plot.setWidth(width=[geometry.getXMax()-geometry.getXMin(),
                     geometry.getXMax()-geometry.getXMin()])
plot.setOrigin([0., 0., 0.])
plot.setPixels([pixels, pixels])

plot_file = PlotsFile()
plot_file.addPlot(plot)
plot_file.exportToXML()

# materials.xml and geometry.xml
create_input_files(geometry)

# Create a plot using OpenMOC's plotting module
plotter.plot_cells(geometry, gridsize=1000)