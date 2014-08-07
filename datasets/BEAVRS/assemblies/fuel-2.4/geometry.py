import opencsg
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

boundaries['Box'] = opencsg.ZSquarePrism(boundary='reflective', R=lattice_width/2.)
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
  cells = universe._cells

  for cell_id in cells.keys():
    cell = cells[cell_id]

    if 'Fuel' in cell._name:
      radialmesh.subdivideCell(cell=cell, universe=universe)


###############################################################################
######################   Creating Root Universe   #############################
###############################################################################

print('Creating the root Universe...')

# Get the appropriate lattice from the lattices module
lattice = lattices['2.4% Fuel - 0BA']

# Root Cell fills the root Universe which encapsulates the complete Geometry
root_universe = opencsg.Universe(universe_id=0, name='Root Universe')
root_cell = opencsg.Cell(name='Root Cell', fill=lattice)
root_universe.addCell(root_cell)

# Append the bounding surfaces for the Geometry to the root Cell
root_cell.addSurface(surface=boundaries['Box'], halfspace=-1)
root_cell.addSurface(surface=boundaries['Z-Min'], halfspace=+1)
root_cell.addSurface(surface=boundaries['Z-Max'], halfspace=-1)

###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

print('Creating the Geometry...')

geometry = opencsg.Geometry()
geometry.setRootUniverse(root_universe)

geometry.initializeCellOffsets()
num_regions = geometry._num_regions
print('# regions = %d' % num_regions)