from openmoc import *
import openmoc.log as log
from openmoc.compatible.openmc import create_input_files
from datasets.BEAVRS.materials import materials
from datasets.BEAVRS.templates import lattice_width
from datasets.BEAVRS.geometry.openmoc import *
from openmc.input.settings import SettingsFile
from openmc.input.plots import PlotsFile, Plot


# NOTE - This is a 3x3 geometry of 1.6% and 3.1% (16BA) enriched fuel assemblies


###############################################################################
###################   Simulation Input File Parameters   ######################
###############################################################################

# Get the appropriate lattice from the lattices module
lattice1 = lattices['1.6% Fuel - 0BA']
lattice2 = lattices['3.1% Fuel - 16BA']
lattice1_id = lattice1.getId()
lattice2_id = lattice2.getId()

# Discretization of pin cells
rings = 3
sectors = 4

# Height of the axial slice
slice_height = 10.

# OpenMC simulation parameters
batches = 25
inactive = 10
particles = 10000

# Plotting parameters
pixels = 1000


###############################################################################
######################   Creating Bounding Surfaces   #########################
###############################################################################

log.py_printf('NORMAL', 'Creating the bounding Surfaces...')

boundaries = dict()

width = lattice_width * 3.

boundaries['X-Min'] = XPlane(x=-width / 2.)
boundaries['X-Max'] = XPlane(x=width / 2.)
boundaries['Y-Min'] = YPlane(y=-width / 2.)
boundaries['Y-Max'] = YPlane(y=width / 2.)
boundaries['Z-Min'] = ZPlane(z=-slice_height / 2.)
boundaries['Z-Max'] = ZPlane(z=slice_height / 2.)

for index in boundaries.keys():
  boundaries[index].setBoundaryType(REFLECTIVE)
  surfaces[index] = boundaries[index]


###############################################################################
#########################   Discretize Pin Cells  #############################
###############################################################################

# Loop over all pin types
for pin in pincells.keys():

  # Loop over all cells for this pin type
  for cell in pincells[pin].keys():

    # Ignore the Universe entry
    if cell == 'Universe':
      continue

    # Set the number of sectors
    pincells[pin][cell].setNumSectors(sectors)

    # Set the number of rings in the fuel
    if 'Fuel' in cell:
      pincells[pin][cell].setNumRings(rings)


###############################################################################
#####################   Creating Colorset Lattice   ###########################
###############################################################################

cell1 = CellFill(universe=universe_id(), universe_fill=lattice1_id)
cell2 = CellFill(universe=universe_id(), universe_fill=lattice2_id)

cell1_id = cell1.getUniverseId()
cell2_id = cell2.getUniverseId()

cells.append(cell1)
cells.append(cell2)

lattice = Lattice(id=universe_id(), width_x=lattice_width, width_y=lattice_width)
lattice.setLatticeCells([[cell2_id, cell1_id, cell2_id],
                         [cell1_id, cell2_id, cell1_id],
                         [cell2_id, cell1_id, cell2_id]])


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
cells.append(root)


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating the Geometry...')

geometry = Geometry()

# Add all Materials to the Geometry
for material in materials.values():
  geometry.addMaterial(material)

# Add all Surfaces to the Geometry
for surface in surfaces.values():
  geometry.addSurface(surface)

for cell in pincells['1.6% Fuel'].values():
  if not isinstance(cell, int):
    geometry.addCell(cell)

for cell in pincells['3.1% Fuel'].values():
  if not isinstance(cell, int):
    geometry.addCell(cell)

for cell in pincells['Guide Tube'].values():
  if not isinstance(cell, int):
    geometry.addCell(cell)

for cell in pincells['Instrument Tube'].values():
  if not isinstance(cell, int):
    geometry.addCell(cell)

for cell in pincells['Burnable Absorber'].values():
  if not isinstance(cell, int):
    geometry.addCell(cell)

geometry.addCell(root)
geometry.addCell(cell1)
geometry.addCell(cell2)

# Add all Lattices to the Geometry
geometry.addLattice(lattice1)
geometry.addLattice(lattice2)
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

source = [-width/2., -width/2., -slice_height/2.,
          width/2., width/2., slice_height/2.]
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