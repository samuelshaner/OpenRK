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

cylinder = ZCylinder(x0=0.0, y0=0.0, R=0.7)
left = XPlane(boundary='reflective', x0=-4.0)
right = XPlane(boundary='reflective', x0=4.0)
top = YPlane(boundary='reflective', y0=4.0)
bottom = YPlane(boundary='reflective', y0=-4.0)

###############################################################################
###########################   Creating Universes  #############################
###############################################################################

print('Creating Universes...')

universes = list()
universes.append(Universe(name='Pin'))
universes.append(Universe(name='Small Lattice'))
universes.append(Universe(universe_id=0, name='Root Universe'))


###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

print('Creating Cells...')

cells = list()
cells.append(Cell(name='Fuel', fill=uo2))
cells.append(Cell(name='Water', fill=water))
cells.append(Cell(name='Small Lattice'))
cells.append(Cell(name='Root Cell'))

cells[0].addSurface(halfspace=-1, surface=cylinder)
cells[1].addSurface(halfspace=+1, surface=cylinder)

cells[3].addSurface(halfspace=+1, surface=left)
cells[3].addSurface(halfspace=-1, surface=right)
cells[3].addSurface(halfspace=-1, surface=top)
cells[3].addSurface(halfspace=+1, surface=bottom)

universes[0].addCells(cells[0:2])
universes[1].addCell(cells[2])
universes[2].addCell(cells[3])


###############################################################################
###########################   Meshing the Cells   #############################
###############################################################################

print('Meshing the Cells...')

'''
mesh = RadialMesh()
mesh.setNumRings(3)
mesh.setMaxRadius(1.0)
mesh.setMinRadius(0.0)
new_cells = mesh.subdivideCell(cell=cells[0], universe=universes[0])

mesh = RadialMesh()
mesh.setNumRings(3)
mesh.setMaxRadius(2.3)
mesh.setMinRadius(1.0)
mesh.setWithOuter(True)
new_cells = mesh.subdivideCell(cell=cells[1], universe=universes[0])

mesh = SectorMesh(num_sectors=8)
mesh.subdivideUniverse(universe=universes[0])
'''

###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

print('Creating Lattices...')

# Initialize Lattices
lattices = list()
lattices.append(Lattice(name='2x2 assembly'))
lattices.append(Lattice(name='2x2 core'))

# 2x2 Pin fuel assembly
lattices[0].setWidth((2.0, 2.0))
lattices[0].setDimension((2, 2))

pin_template = [[1, 1],
                [1, 1]]

for i in range(len(pin_template)):
  for j in range(len(pin_template[i])):
    pin_template[i][j] = universes[pin_template[i][j]-1]

lattices[0].setUniverses(pin_template)
cells[2].setFill(lattices[0])

# 2x2 fuel assembly core
lattices[1].setWidth((4.0,4.0))
lattices[1].setDimension((2,2))

assembly_template = [[universes[1], universes[1]],
                     [universes[1], universes[1]]]

lattices[1].setUniverses(assembly_template)
cells[3].setFill(lattices[1])

###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

print('Creating Geometry...')

geometry = Geometry()
geometry.setRootUniverse(universes[2])

geometry.initializeCellOffsets()

###############################################################################
###############################   Ray Tracing   ###############################
###############################################################################

print('Tracing Sample Rays...')

rays = geometry.generateRays(num_rays=3000)
rays = geometry.traceRays(rays)

###############################################################################
##########################   Plotting the Geometry   ##########################
###############################################################################

#print('Plotting Geometry...')

plotter.plot_segments(rays, geometry)
#plotter.plot_cells(geometry)
#plotter.plot_materials(geometry)
plotter.plot_regions(geometry)

