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
universes.append(Universe(name='Pin 1'))
universes.append(Universe(name='Pin 2'))
universes.append(Universe(name='Pin 3'))
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
cells.append(Cell(name='Root Cell'))

cells[0].addSurface(halfspace=-1, surface=cylinders[0])
cells[1].addSurface(halfspace=+1, surface=cylinders[0])
cells[2].addSurface(halfspace=-1, surface=cylinders[1])
cells[3].addSurface(halfspace=+1, surface=cylinders[1])
cells[4].addSurface(halfspace=-1, surface=cylinders[2])
cells[5].addSurface(halfspace=+1, surface=cylinders[2])

cells[6].addSurface(halfspace=+1, surface=planes[0])
cells[6].addSurface(halfspace=-1, surface=planes[1])
cells[6].addSurface(halfspace=+1, surface=planes[2])
cells[6].addSurface(halfspace=-1, surface=planes[3])

universes[0].addCells(cells[0:2])
universes[1].addCells(cells[2:4])
universes[2].addCells(cells[4:6])
universes[3].addCell(cells[6])


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

#mesh = SectorMesh(num_sectors=8)
#mesh.subdivideUniverse(universe=universes[0])

#mesh = SectorMesh(num_sectors=8)
#mesh.subdivideUniverse(universe=universes[1])

#mesh = SectorMesh(num_sectors=8)
#mesh.subdivideUniverse(universe=universes[2])


###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

print('Creating Lattices...')

lattices = list()
lattices.append(Lattice(name='4x4'))
lattices[0].setWidth((1.0, 1.0))
lattices[0].setLowerLeft((-2.0, -2.0))
lattices[0].setDimension((4, 4))

template = [[1, 2, 1, 2],
            [2, 3, 2, 3],
            [1, 2, 1, 2],
            [2, 3, 2, 3]]

for i in range(len(template)):
  for j in range(len(template[i])):
    template[i][j] = universes[template[i][j]-1]

lattices[0].setUniverses(template)

cells[6].setFill(lattices[0])


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

print('Creating Geometry...')

geometry = Geometry()
geometry.setRootUniverse(universes[3])

geometry.initializeCellOffsets()
geometry.setVolume(16., tolerance=1e-1)

num_regions = geometry.getNumRegions()
print('The Geometry contains %d regions' % num_regions)

tot_volume = 0.
for region in range(num_regions):
  coords = geometry.findRegion(region)
  coords = coords.getTailNode()
  cell = coords.getCell()
  volume = cell.getVolume()
  tot_volume += volume
  print('region = %d, volume = %f, tot volume = %f' %
        (region, volume, tot_volume))

print('Total volume = %f' % tot_volume)

###############################################################################
##########################   Plotting the Geometry   ##########################
###############################################################################

print('Plotting Geometry...')

plotter.plot_cells(geometry)
#plotter.plot_materials(geometry)
#plotter.plot_regions(geometry)
