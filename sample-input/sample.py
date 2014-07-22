from opencsg.point import *
from opencsg.surface import *
from opencsg.universe import *
from opencsg.geometry import *

point = Point()
direction = Direction(u=1, v=0, w=0)
plane = XPlane(surface_id=0, name='Billy', x0=1.)
universe = Universe(universe_id=0, name='Bonnie')
geometry = Geometry()
material = Material(material_id=0, name='Bo')
cell = Cell(cell_id=0, name='Bob')
cell.addSurface(plane, +1)
cell.setFill(material)
universe.addCell(cell)
geometry.setRootUniverse(universe)


intersect = geometry.getNearestIntersection(point, direction)
print intersect
print plane.evaluate(intersect)
