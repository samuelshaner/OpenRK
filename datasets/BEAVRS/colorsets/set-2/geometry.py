import opencsg
from datasets.BEAVRS.lattices import *


###############################################################################
##################   Geometry Discretization Parameters   #####################
###############################################################################

# Discretization of pin cells
fuel_rings = 3
mod_rings = 3
sectors = 4

# Height of the axial slice
slice_height = 10.


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

# Get the appropriate lattice from the lattices module
lattice1 = lattices['1.6% Fuel - 0BA']
lattice2 = lattices['3.1% Fuel - 16BA']

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