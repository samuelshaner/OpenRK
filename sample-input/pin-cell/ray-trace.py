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
left = XPlane(boundary='reflective', x0=-2.0)
right = XPlane(boundary='reflective', x0=2.0)
top = YPlane(boundary='reflective', y0=2.0)
bottom = YPlane(boundary='reflective', y0=-2.0)

###############################################################################
###########################   Creating Universes  #############################
###############################################################################

print('Creating Universes...')

universes = list()
universes.append(Universe(name='Pin'))
universes.append(Universe(universe_id=0, name='Root Universe'))


###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

print('Creating Cells...')

cells = list()
cells.append(Cell(name='Fuel', fill=uo2))
cells.append(Cell(name='Water', fill=water))
cells.append(Cell(name='Root Cell'))

cells[0].addSurface(halfspace=-1, surface=cylinder)
cells[1].addSurface(halfspace=+1, surface=cylinder)
cells[2].addSurface(halfspace=+1, surface=left)
cells[2].addSurface(halfspace=-1, surface=right)
cells[2].addSurface(halfspace=-1, surface=top)
cells[2].addSurface(halfspace=+1, surface=bottom)

universes[0].addCells(cells[0:2])
universes[1].addCell(cells[2])


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

lattice = Lattice(name='1x1')
lattice.setWidth((4.0, 4.0))
lattice.setDimension((1, 1))
lattice.setUniverses([[universes[0]]])

cells[2].setFill(lattice)


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

print('Creating Geometry...')

geometry = Geometry()
geometry.setRootUniverse(universes[1])

geometry.initializeCellOffsets()
geometry.setVolume(volume=16., tolerance=1e-1)

###############################################################################
###############################   Ray Tracing   ###############################
###############################################################################

print('Tracing Sample Rays...')

segments = geometry.traceSampleRays(num_rays=50000)


###############################################################################
##########################   Plotting the Geometry   ##########################
###############################################################################

#print('Plotting Geometry...')

plotter.plot_segments(segments, geometry)
#plotter.plot_cells(geometry)
#plotter.plot_materials(geometry)
#plotter.plot_regions(geometry)

