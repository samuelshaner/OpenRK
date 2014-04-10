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

planes = list()
cylinders = list()
planes.append(XPlane(x0=-2.0))
planes.append(XPlane(x0=2.0))
planes.append(YPlane(y0=-2.0))
planes.append(YPlane(y0=2.0))

cylinders.append(ZCylinder(x0=0.0, y0=0.0, R=0.4))
cylinders.append(ZCylinder(x0=0.0, y0=0.0, R=0.3))
cylinders.append(ZCylinder(x0=0.0, y0=0.0, R=0.2))

for plane in planes: plane.setBoundaryType('reflective')


###############################################################################
###########################   Creating Universes  #############################
###############################################################################

print('Creating Universes...')

universes = list()
universes.append(Universe(name='pin 1'))
universes.append(Universe(name='pin 2'))
universes.append(Universe(name='pin 3'))
universes.append(Universe(name='small lattice'))
universes.append(Universe(universe_id=0, name='root'))


###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

print('Creating Cells...')

cells = list()
cells.append(Cell(name='pin 1 fuel', fill=uo2))
cells.append(Cell(name='pin 1 water', fill=water))
cells.append(Cell(name='pin 2 fuel', fill=uo2))
cells.append(Cell(name='pin 2 water', fill=water))
cells.append(Cell(name='pin 3 fuel', fill=uo2))
cells.append(Cell(name='pin 3 water', fill=water))
cells.append(Cell(name='small lattice'))
cells.append(Cell(name='root'))

cells[0].addSurface(halfspace=-1, surface=cylinders[0])
cells[1].addSurface(halfspace=+1, surface=cylinders[0])
cells[2].addSurface(halfspace=-1, surface=cylinders[1])
cells[3].addSurface(halfspace=+1, surface=cylinders[1])
cells[4].addSurface(halfspace=-1, surface=cylinders[2])
cells[5].addSurface(halfspace=+1, surface=cylinders[2])

cells[7].addSurface(halfspace=+1, surface=planes[0])
cells[7].addSurface(halfspace=-1, surface=planes[1])
cells[7].addSurface(halfspace=+1, surface=planes[2])
cells[7].addSurface(halfspace=-1, surface=planes[3])

universes[0].addCells(cells[0:2])
universes[1].addCells(cells[2:4])
universes[2].addCells(cells[4:6])
universes[3].addCell(cells[6])
universes[4].addCell(cells[7])


###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

print('Creating Lattices...')

# Initialize lattices
lattices = list()
lattices.append(Lattice(name='2x2 assembly'))
lattices.append(Lattice(name='2x2 core'))

# 2x2 pin fuel assembly
lattices[0].setWidth((1.0, 1.0))
lattices[0].setLowerLeft((-1.0, -1.0))
lattices[0].setDimension((2, 2))

pin_template = [[1, 2], [1, 3]]

for i in range(len(pin_template)):
  for j in range(len(pin_template[i])):
    pin_template[i][j] = universes[pin_template[i][j]-1]

lattices[0].setUniverses(pin_template)
cells[6].setFill(lattices[0])


# 2x2 fuel assembly core
lattices[1].setWidth((2.0, 2.0))
lattices[1].setLowerLeft((-2.0, -2.0))
lattices[1].setDimension((2, 2))

assembly_template = [[universes[3], universes[3]],
                     [universes[3], universes[3]]]

lattices[1].setUniverses(assembly_template)
cells[7].setFill(lattices[1])


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

print('Creating Geometry...')

geometry = Geometry()

for universe in universes: geometry.addUniverse(universe)
for lattice in lattices: geometry.addLattice(lattice)

geometry.initializeCellOffsets()
geometry.setVolume(16., tolerance=1e-2)

plotter.plot_cells(geometry)
plotter.plot_materials(geometry)
plotter.plot_regions(geometry)
