import opencg
from datasets.BEAVRS.materials import openmc_materials, opencsg_materials
from datasets.BEAVRS.lattices import *


# NOTE - This is a quarter core BEAVRS geometry.


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

width = lattice_width * 10.

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

cell1 = opencg.Cell(name=lat1._name, fill=lat1)
cell2 = opencg.Cell(name=lat2._name, fill=lat2)
cell3 = opencg.Cell(name=lat3._name, fill=lat3)
cell4 = opencg.Cell(name=lat4._name, fill=lat4)
cell5 = opencg.Cell(name=lat5._name, fill=lat5)
cell6 = opencg.Cell(name=lat6._name, fill=lat6)
cell7 = opencg.Cell(name=lat7._name, fill=lat7)
cell8 = opencg.Cell(name=lat8._name, fill=lat8)
cell9 = opencg.Cell(name=lat9._name, fill=lat9)

# Create a pure water Cell
water = opencsg_materials['Borated Water']
cell10 = opencg.Cell(name='Water cell', fill=water)

u1 = opencg.Universe(name=cell1._name)
u1.addCell(cell1)

u2 = opencg.Universe(name=cell2._name)
u2.addCell(cell2)

u3 = opencg.Universe(name=cell3._name)
u3.addCell(cell3)

u4 = opencg.Universe(name=cell4._name)
u4.addCell(cell4)

u5 = opencg.Universe(name=cell5._name)
u5.addCell(cell5)

u6 = opencg.Universe(name=cell6._name)
u6.addCell(cell6)

u7 = opencg.Universe(name=cell7._name)
u7.addCell(cell7)

u8 = opencg.Universe(name=cell8._name)
u8.addCell(cell8)

u9 = opencg.Universe(name=cell9._name)
u9.addCell(cell9)

u10 = opencg.Universe(name=cell10._name)
u10.addCell(cell10)

lattice = Lattice(name='Quarter Core')
lattice.setDimension((10, 10))
lattice.setWidth((lattice_width, lattice_width))
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