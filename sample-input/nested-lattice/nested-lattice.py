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

cylinders = list()
cylinders.append(ZCylinder(x0=0.0, y0=0.0, R=0.4))
cylinders.append(ZCylinder(x0=0.0, y0=0.0, R=0.3))
cylinders.append(ZCylinder(x0=0.0, y0=0.0, R=0.2))

square = Square(x0=0., y0=0., R=2.)
square.setBoundaryType('reflective')


###############################################################################
###########################   Creating Universes  #############################
###############################################################################

print('Creating Universes...')

universes = list()
universes.append(Universe(name='Pin 1'))
universes.append(Universe(name='Pin 2'))
universes.append(Universe(name='Pin 3'))
universes.append(Universe(name='Small Lattice'))
universes.append(Universe(universe_id=0, name='Root Universe'))


###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

print('Creating Cells...')

cells = list()
cells.append(Cell(name='Pin 1 Fuel', fill=uo2))
cells.append(Cell(name='Pin 1 Water', fill=water))
cells.append(Cell(name='Pin 2 Fuel', fill=uo2))
cells.append(Cell(name='Pin 2 Water', fill=water))
cells.append(Cell(name='Pin 3 Fuel', fill=uo2))
cells.append(Cell(name='Pin 3 Water', fill=water))
cells.append(Cell(name='Small Lattice'))
cells.append(Cell(name='Root Cell'))

cells[0].addSurface(halfspace=-1, surface=cylinders[0])
cells[1].addSurface(halfspace=+1, surface=cylinders[0])
cells[2].addSurface(halfspace=-1, surface=cylinders[1])
cells[3].addSurface(halfspace=+1, surface=cylinders[1])
cells[4].addSurface(halfspace=-1, surface=cylinders[2])
cells[5].addSurface(halfspace=+1, surface=cylinders[2])

cells[7].addSurface(halfspace=-1, surface=square)

universes[0].addCells(cells[0:2])
universes[1].addCells(cells[2:4])
universes[2].addCells(cells[4:6])
universes[3].addCell(cells[6])
universes[4].addCell(cells[7])


###############################################################################
###########################   Meshing the Cells   #############################
###############################################################################

print('Meshing the Cells...')

mesh = RadialMesh()
mesh.setNumRings(3)
mesh.setMaxRadius(cells[0].getMaxX())
mesh.setMinRadius(0.0)
new_cells = mesh.subdivideCell(cell=cells[0], universe=universes[0])

mesh = RadialMesh()
mesh.setNumRings(3)
mesh.setMaxRadius(cells[2].getMaxX())
mesh.setMinRadius(0.0)
new_cells = mesh.subdivideCell(cell=cells[2], universe=universes[1])

mesh = RadialMesh()
mesh.setNumRings(3)
mesh.setMaxRadius(cells[4].getMaxX())
mesh.setMinRadius(0.0)
new_cells = mesh.subdivideCell(cell=cells[4], universe=universes[2])

mesh = SectorMesh(num_sectors=8)
mesh.subdivideUniverse(universe=universes[0])

mesh = SectorMesh(num_sectors=8)
mesh.subdivideUniverse(universe=universes[1])

mesh = SectorMesh(num_sectors=8)
mesh.subdivideUniverse(universe=universes[2])


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
lattices[0].setDimension((2, 2))

pin_template = [[1, 2], [1, 3]]

for i in range(len(pin_template)):
  for j in range(len(pin_template[i])):
    pin_template[i][j] = universes[pin_template[i][j]-1]

lattices[0].setUniverses(pin_template)
cells[6].setFill(lattices[0])


# 2x2 fuel assembly core
lattices[1].setWidth((2.0, 2.0))
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
geometry.setRootUniverse(universes[4])

geometry.initializeCellOffsets()
geometry.setVolume(volume=16., tolerance=1e-1)


###############################################################################
##########################   Plotting the Geometry   ##########################
###############################################################################

print('Plotting Geometry...')

plotter.plot_cells(geometry)
plotter.plot_materials(geometry)
plotter.plot_regions(geometry)
#plotter.plot_neighbor_cells(geometry)
