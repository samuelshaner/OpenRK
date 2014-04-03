from openmoc.compatible.openmc import *
import openmoc.plotter as plotter
from lattices import *
import materials
import surfaces
import pincells
from input.settings import SettingsFile
from input.tallies import TalliesFile, Tally
from input.plots import PlotsFile, Plot


# NOTE - This is the 3x3 geometry of 1.6% and 2.4% (16BA) assemblies found at
#        the center of the BEAVRS core.


###############################################################################
###################   Simulation Input File Parameters   ######################
###############################################################################

# Get the appropriate lattice from the lattices module
lattice1 = lattices['1.6% Fuel - 0BA']
lattice2 = lattices['2.4% Fuel - 16BA']
lattice1_id = lattice1.getId()
lattice2_id = lattice2.getId()

# Discretization of pin cells
rings = 3
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

cell1 = CellFill(universe=universe_id(), universe_fill=lattice1_id)
cell2 = CellFill(universe=universe_id(), universe_fill=lattice2_id)

cell1_id = cell1.getUniverseId()
cell2_id = cell2.getUniverseId()

pincells.cells.append(cell1)
pincells.cells.append(cell2)

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

# Add all Lattices to the Geometry
geometry.addLattice(lattice1)
geometry.addLattice(lattice2)
geometry.addLattice(lattice)


###############################################################################
###################   Exporting to OpenMC XML Input Files  ####################
###############################################################################

log.py_printf('NORMAL', 'Exporting to OpenMC XML Files...')

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


# tallies.xml
#num_cells = geometry.getNumCells()
#cell_ids = geometry.getCellIds(num_cells)

#tallies_file = openmc.TalliesFile()
#scores = ['flux', 'total', 'nu-fission']

#print cell_ids

#for index, cell_id in enumerate(cell_ids):

#  if cell_id == 10027:

#    print index, cell_id

#    tally = openmc.Tally(label='test')
#    tally.addFilter(type='distribcell', bins=cell_id)
#    tally.addFilter(type='energy', bins=[0.0, 0.625, 10000000.])

#    for score in scores:
#      tally.addScore(score=score)

#    tallies_file.addTallySubelement(tally)

#tallies_file.exportToXML()


# This must be called after the
# hmm = sp.get_value(0, [('distribcell',2)], 0)[0]


# tallies.xml
#tally = openmc.Tally(label='test')
#tally.addFilter(type='distribcell', bins=10034)
#tally.addScore(score='flux')
#tallies_file = openmc.TalliesFile()
#tallies_file.addTallySubelement(tally)
#tallies_file.exportToXML()


#from statepoint import StatePoint
#import openmoc.compatible.openmc.plotter.plotter as plot

#sp = StatePoint('statepoint.15.h5')
#sp.read_results()

#plot.plot_fluxes(geometry, sp, energies=[0], gridsize=100)