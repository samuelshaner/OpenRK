import openrk as rk
import openrk.compatible as compat
import openrk.mesh
import openrk.material

from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.materialize as materialize
from openmoc.options import Options

options = Options()

num_threads = options.getNumThreads()
track_spacing = options.getTrackSpacing()
num_azim = options.getNumAzimAngles()
tolerance = options.getTolerance()
max_iters = options.getMaxIterations()

log.set_log_level('NORMAL')


###############################################################################
###########################   Creating Materials   ############################
###############################################################################

log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = materialize.materialize('../OpenMOC-mitcrpg/sample-input/c5g7-materials.h5')

uo2_id = materials['UO2'].getId()
water_id = materials['Water'].getId()


###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

circles = []
planes = []
planes.append(XPlane(x=-2.0))
planes.append(XPlane(x=2.0))
planes.append(YPlane(y=-2.0))
planes.append(YPlane(y=2.0))
circles.append(Circle(x=0.0, y=0.0, radius=0.4))
circles.append(Circle(x=0.0, y=0.0, radius=0.3))
circles.append(Circle(x=0.0, y=0.0, radius=0.2))
for plane in planes: plane.setBoundaryType(REFLECTIVE)


###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

cells = []
cells.append(CellBasic(universe=1, material=uo2_id, rings=3, sectors=8))
cells.append(CellBasic(universe=1, material=water_id, sectors=8))
cells.append(CellBasic(universe=2, material=uo2_id, rings=3, sectors=8))
cells.append(CellBasic(universe=2, material=water_id, sectors=8))
cells.append(CellBasic(universe=3, material=uo2_id, rings=3, sectors=8))
cells.append(CellBasic(universe=3, material=water_id, sectors=8))
cells.append(CellFill(universe=0, universe_fill=5))

cells[0].addSurface(halfspace=-1, surface=circles[0])
cells[1].addSurface(halfspace=+1, surface=circles[0])
cells[2].addSurface(halfspace=-1, surface=circles[1])
cells[3].addSurface(halfspace=+1, surface=circles[1])
cells[4].addSurface(halfspace=-1, surface=circles[2])
cells[5].addSurface(halfspace=+1, surface=circles[2])

cells[6].addSurface(halfspace=+1, surface=planes[0])
cells[6].addSurface(halfspace=-1, surface=planes[1])
cells[6].addSurface(halfspace=+1, surface=planes[2])
cells[6].addSurface(halfspace=-1, surface=planes[3])


###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating simple 4 x 4 lattice...')

lattice = Lattice(id=5, width_x=1.0, width_y=1.0)
lattice.setLatticeCells([[1, 2, 1, 2],
                         [2, 3, 2, 3],
                         [1, 2, 1, 2],
                         [2, 3, 2, 3]])

###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()
for material in materials.values(): geometry.addMaterial(material)
for cell in cells: geometry.addCell(cell)
geometry.addLattice(lattice)

geometry.setTcmfdLattice(2,2);
geometry.setTcmfdGroupStructure([1,4,8])
geometry.initializeFlatSourceRegions()

###############################################################################
########################   Creating the TrackGenerator   ######################
###############################################################################

log.py_printf('NORMAL', 'Initializing the track generator...')

track_generator = TrackGenerator(geometry, num_azim, track_spacing)
track_generator.setNumThreads(num_threads)
track_generator.generateTracks()

moc_mesh = compat.extract_openrk_fsr_mesh(geometry)
print moc_mesh

###############################################################################
###########################   Running a Simulation   ##########################
###############################################################################

solver = TransientCPUSolver(geometry, track_generator)
solver.setNumThreads(num_threads)
solver.setSourceConvergenceThreshold(tolerance)
solver.convergeSource(max_iters)
solver.printTimerReport()

tcmfd_mesh = compat.extract_openrk_tcmfd_mesh(geometry)
print tcmfd_mesh

