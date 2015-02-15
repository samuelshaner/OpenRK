import opencg
from datasets.BEAVRS.lattices import *


###############################################################################
##################   Geometry Discretization Parameters   #####################
###############################################################################

# Discretization of pin cells
fuel_rings = 0
mod_rings = 0
sectors = 0

# Height of the axial slice
slice_height = 10.


###############################################################################
######################   Creating Bounding Surfaces   #########################
###############################################################################

print('Creating the bounding Surfaces...')

boundaries = dict()

width = lattice_width * 3.

boundaries['Box'] = opencg.ZSquarePrism(boundary='reflective', R=width/2.)
boundaries['Z-Min'] = opencg.ZPlane(z0=-slice_height / 2.)
boundaries['Z-Max'] = opencg.ZPlane(z0=slice_height / 2.)

for index in boundaries.keys():
  boundaries[index].setBoundaryType('reflective')


###############################################################################
#########################   Discretize Pin Cells  #############################
###############################################################################

print('Discretizing the Cells...')

sectormesh = opencg.SectorMesh(num_sectors=sectors)
radialmesh = opencg.RadialMesh(num_rings=fuel_rings)

# Loop over all pin types
for pin in pin_cells.keys():

  universe = pin_cells[pin]
  sectormesh.subdivideUniverse(universe=universe)
  cells = universe._cells

  for cell_id in cells.keys():
    cell = cells[cell_id]

    if 'Fuel' in cell._name:
      radialmesh.subdivideCell(cell=cell, universe=universe)


###############################################################################
#####################   Creating Colorset Lattice   ###########################
###############################################################################

# Get the appropriate lattice from the lattices module
lattice1 = lattices['1.6% Fuel - 0BA']
lattice2 = lattices['2.4% Fuel - 16BA']

cell1 = opencg.Cell(name=lattice1._name, fill=lattice1)
cell2 = opencg.Cell(name=lattice2._name, fill=lattice2)

universe1 = opencg.Universe(name=cell1._name)
universe1.addCell(cell1)

universe2 = opencg.Universe(name=cell2._name)
universe2.addCell(cell2)

lattice = opencg.Lattice(name='3x3 Lattice')
lattice.setDimension((3, 3))
lattice.setWidth((lattice_width, lattice_width))
lattice.setUniverses([[universe2, universe1, universe2],
                      [universe1, universe2, universe1],
                      [universe2, universe1, universe2]])


###############################################################################
######################   Creating Root Universe   #############################
###############################################################################

print('Creating the root Universe...')

# Root Cell fills the root Universe which encapsulates the complete Geometry
root_universe = opencg.Universe(universe_id=0, name='Root Universe')
root_cell = opencg.Cell(name='Root Cell', fill=lattice)
root_universe.addCell(root_cell)

# Append the bounding surfaces for the Geometry to the root Cell
root_cell.addSurface(surface=boundaries['Box'], halfspace=-1)
root_cell.addSurface(surface=boundaries['Z-Min'], halfspace=+1)
root_cell.addSurface(surface=boundaries['Z-Max'], halfspace=-1)


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

print('Creating the Geometry...')

geometry = opencg.Geometry()
geometry.setRootUniverse(root_universe)

geometry.initializeCellOffsets()
num_regions = geometry._num_regions
print('# regions = %d' % num_regions)
