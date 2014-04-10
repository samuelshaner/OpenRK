from opencsg import *
import opencsg.plotter as plotter


###############################################################################
###########################   Creating Materials   ############################
###############################################################################

print('Creating Materials...')

uo2 = Material(name='UO2')
water = Material(name='Water')


###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

print('Creating Surfaces...')

cylinder = ZCylinder(x0=0.0, y0=0.0, R=1.0)
left = XPlane(x0=-2.0)
right = XPlane(x0=2.0)
top = YPlane(y0=2.0)
bottom = YPlane(y0=-2.0)

left.setBoundaryType('reflective')
right.setBoundaryType('reflective')
top.setBoundaryType('reflective')
bottom.setBoundaryType('reflective')


###############################################################################
###########################   Creating Universes  #############################
###############################################################################

print('Creating Universes...')

universes = list()
universes.append(Universe(name='pin'))
universes.append(Universe(universe_id=0, name='root'))


###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

print('Creating Cells...')

cells = list()
cells.append(Cell(name='fuel', fill=uo2))
cells.append(Cell(name='water', fill=water))
cells.append(Cell(name='root'))

cells[0].addSurface(halfspace=-1, surface=cylinder)
cells[1].addSurface(halfspace=+1, surface=cylinder)
cells[2].addSurface(halfspace=+1, surface=left)
cells[2].addSurface(halfspace=-1, surface=right)
cells[2].addSurface(halfspace=+1, surface=bottom)
cells[2].addSurface(halfspace=-1, surface=top)

universes[0].addCells(cells[0:2])
universes[1].addCell(cells[2])


###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

print('Creating Lattices...')

lattice = Lattice(name='1x1')
lattice.setWidth((4.0, 4.0))
lattice.setLowerLeft((-2.0, -2.0))
lattice.setDimension((1, 1))
lattice.setUniverses([[universes[0]]])

cells[2].setFill(lattice)


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

print('Creating Geometry...')

geometry = Geometry()

for universe in universes: geometry.addUniverse(universe)
geometry.addLattice(lattice)

geometry.initializeCellOffsets()
geometry.setVolume(16.)

plotter.plot_cells(geometry)
plotter.plot_materials(geometry)
plotter.plot_regions(geometry)