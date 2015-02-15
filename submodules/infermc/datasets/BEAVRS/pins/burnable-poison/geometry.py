import opencg
from datasets.BEAVRS.pin_cells import *
from datasets.BEAVRS.lattices import pin_pitch


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

boundaries['Box'] = opencg.ZSquarePrism(boundary='reflective', R=pin_pitch/2.)
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

universe = pin_cells['Burnable Absorber']
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

# Root Cell fills the root Universe which encapsulates the complete Geometry
root_universe = opencg.Universe(universe_id=0, name='Root Universe')
root_cell = opencg.Cell(name='Root Cell', fill=universe)
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