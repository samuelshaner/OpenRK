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


mesh = RadialMesh(cell=cells[0])
mesh.setSpacingType('2D')
mesh.setNumRings(5)
mesh.setMaxRadius(1.0)
mesh.setMinRadius(0.0)
new_cells = mesh.subdivideCell()

universes[0].removeCell(cells[0])
universes[0].addCells(new_cells)

mesh = RadialMesh(cell=cells[1])
mesh.setSpacingType('2D')
mesh.setNumRings(5)
mesh.setMaxRadius(2.3)
mesh.setMinRadius(1.0)
mesh.setWithOuter(True)
new_cells = mesh.subdivideCell()

universes[0].removeCell(cells[1])
universes[0].addCells(new_cells)


cells = universes[0].getCells()
for cell_id in cells:
  cell = cells[cell_id]
  cell.printString()




###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

print('Creating Geometry...')

geometry = Geometry()
geometry.setRootUniverse(universes[1])

universes[0].printString()
universes[1].printString()

geometry.initializeCellOffsets()
#geometry.setVolume(volume=16., tolerance=1e-2)


###############################################################################
##########################   Plotting the Geometry   ##########################
###############################################################################

print('Plotting Geometry...')

plotter.plot_cells(geometry)
# #plotter.plot_materials(geometry)
#plotter.plot_regions(geometry)