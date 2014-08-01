from opencsg import *
import opencsg.plotter as plotter
import numpy as np
import matplotlib as plt
import pylab as pl

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


rays = []
bounds = geometry.getBounds()

for ray in xrange(1000):
  edge = np.random.choice([0,1,2,3])
  if edge == 0:
    x = bounds[edge] + 1e-10
    y = np.random.uniform(bounds[2], bounds[3])
    z = np.random.uniform(1e-12, 1e12)
  elif edge == 1:
    x = bounds[edge] - 1e-10
    y = np.random.uniform(bounds[2], bounds[3])
    z = np.random.uniform(1e-12, 1e12)
  elif edge == 2:
    x = np.random.uniform(bounds[0], bounds[1])
    y = bounds[edge] + 1e-10
    z = np.random.uniform(1e-12, 1e12)
  else:
    x = np.random.uniform(bounds[0], bounds[1])
    y = bounds[edge] - 1e-10
    z = np.random.uniform(1e-12, 1e12)

  u, v, w = np.random.rand(3)
  point = Point(x=x,y=y,z=z)
  direction = Direction(u=u,v=v,w=w)
  rays.append((point, direction))

colors = []
segments = []

#Code below generates rays that do not go through surfaces
#Picture is still, however, off center

for ray in rays:
  intersect = geometry.getNearestIntersection(ray[0], ray[1])
  if not intersect is None and ray[0].distanceToPoint(intersect)<10:
    segments.append([ray[0]._coords[:2], intersect._coords[:2]])
    if cylinder.evaluate(ray[0]) < 0:
      colors.append((1, 0, 1, 1))
    else:
      colors.append((1, 1, 0, 1))

#Code below attempts to generate rays that go through surfaces
#Picture is off center and goes crazy
'''
for ray in rays:
  intersect = geometry.getNearestIntersection(ray[0], ray[1])
  while not intersect is None:
    print ray[0]
    segments.append([ray[0]._coords[:2], intersect._coords[:2]])
    if cylinder.evaluate(ray[0]) < 0:
      colors.append((1, 0, 1, 1))
    else:
      colors.append((1, 1, 0, 1))
    ray[0].setCoords((intersect._coords + 1e-10*ray[1]._comps))
    print ray[0]
    intersect = geometry.getNearestIntersection(ray[0], ray[1])
'''

def startToEnd(rays, geometry):
  #Code below plots initial and intersect points
  init_x_vals = []
  init_y_vals = []
  int_x_vals = []
  int_y_vals = []
  for ray in rays:
    intersect = geometry.getNearestIntersection(ray[0], ray[1])
    if not intersect is None and ray[0].distanceToPoint(intersect)<10:
      init_x_vals.append(ray[0]._coords[0])
      init_y_vals.append(ray[0]._coords[1])
      int_x_vals.append(intersect._coords[0])
      int_y_vals.append(intersect._coords[1])

  fig = plt.pyplot.figure()
  plt.pyplot.plot(init_x_vals, init_y_vals, 'ro')
  plt.pyplot.plot(int_x_vals, int_y_vals, 'bo')
  plt.pyplot.show()
  fig.savefig('plots/start-end.png')

c = np.array(colors)

lc = plt.collections.LineCollection(segments, colors=c, linewidths=1)
fig, ax = pl.subplots()
ax.add_collection(lc)
ax.autoscale()
ax.margins(0)
fig.savefig('plots/rays-xy.png')

startToEnd(rays, geometry)

###############################################################################
##########################   Plotting the Geometry   ##########################
###############################################################################

print('Plotting Geometry...')

#plotter.plot_cells(geometry)
#plotter.plot_materials(geometry)
#plotter.plot_regions(geometry)

