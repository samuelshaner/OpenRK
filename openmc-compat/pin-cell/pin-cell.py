from openmoc import *
import openmoc.log as log
import openmoc.materialize as materialize
import openmoc.compatible.openmc as openmc
import numpy as np

###############################################################################
###########################   Creating Materials   ############################
###############################################################################

log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = materialize.materialize('../c5g7-materials.h5')

uo2_id = materials['UO2'].getId()
water_id = materials['Water'].getId()


###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

circle = Circle(x=0.0, y=0.0, radius=1.0)
left = XPlane(x=-2.0)
right = XPlane(x=2.0)
top = YPlane(y=2.0)
bottom = YPlane(y=-2.0)

left.setBoundaryType(REFLECTIVE)
right.setBoundaryType(REFLECTIVE)
top.setBoundaryType(REFLECTIVE)
bottom.setBoundaryType(REFLECTIVE)


###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

cells = []
cells.append(CellBasic(universe=1, material=uo2_id, rings=2, sectors=8))
cells.append(CellBasic(universe=1, material=water_id))
cells.append(CellFill(universe=0, universe_fill=2))

cells[0].addSurface(halfspace=-1, surface=circle)
cells[1].addSurface(halfspace=+1, surface=circle)
cells[2].addSurface(halfspace=+1, surface=left)
cells[2].addSurface(halfspace=-1, surface=right)
cells[2].addSurface(halfspace=+1, surface=bottom)
cells[2].addSurface(halfspace=-1, surface=top)


###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating simple pin cell lattice...')

lattice = Lattice(id=2, width_x=16.0, width_y=16.0)
lattice.setLatticeCells([[1,1],
                         [1,1]])


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()
for material in materials.values(): geometry.addMaterial(material)
for cell in cells: geometry.addCell(cell)
geometry.addLattice(lattice)


###############################################################################
###################   Exporting to OpenMC XML Input Files  ####################
###############################################################################


geometry_file = openmc.GeometryFile()
geometry_file.createGeometrySubelements(geometry)
geometry_file.exportToXML()

settings_file = openmc.SettingsFile()
settings_file.createEigenvalueSubelement()
settings_file.createOutputSubelement(['tallies'])
settings_file.createStatepointSubelement(batches=[1,2,3,4,5])
settings_file.createSourceSpaceSubelement(type='box', params=[-1,-1,-1,1,1,1])
settings_file.exportToXML()

plot_file = openmc.PlotsFile()
plot_file.addNewPlot(id=1, width=[4.0,4.0], origin=[0.,0.,0.], pixels=[100,100])
plot_file.exportToXML()


h1 = Isotope(isotope_id(), 'H-1')
o16 = Isotope(isotope_id(), 'O-16')
u234 = Isotope(isotope_id(), 'U-234')
u235 = Isotope(isotope_id(), 'U-235')
u238 = Isotope(isotope_id(), 'U-238')
h1.setNumEnergyGroups(7)
o16.setNumEnergyGroups(7)
u234.setNumEnergyGroups(7)
u235.setNumEnergyGroups(7)
u238.setNumEnergyGroups(7)

h2o = IsoMaterial(material_id(), 'H2O')
h2o.setDensity(0.99, 'g/cc')
h2o.addIsotope(h1, 2.0)
h2o.addIsotope(o16, 1.0)

u2o = IsoMaterial(material_id(), 'UO2')
u2o.setDensity(10., 'g/cc')
u2o.addIsotope(o16, 2.0)
u2o.addIsotope(u234, 0.001)
u2o.addIsotope(u235, 0.04)
u2o.addIsotope(u238, 0.959)

materials = [u2o, h2o]

materials_file = openmc.MaterialsFile()
materials_file.createDefaultXSSubelement()
materials_file.createMaterialsSubelements(materials)
materials_file.exportToXML()

tallies_file = openmc.TalliesFile()
tallies_file.addMeshSubelement(id=1, dimension=[100,100,1], lower_left=[-1.,-1,-1.], upper_right=[-2.,-2.,-2.])
tally = openmc.Tally()
tally.setLabel(label='first tally')
tally.addNuclide('H-1')
tally.addNuclide('U-238')
tally.addNuclide('total')
tally.addScore('flux')
tally.addScore('total')
tally.addScore('fission')
tally.addFilter(type='cell', bins=[1,3,5])
tally.addFilter(type='cellborn', bins=[1,2,3])
tally.addFilter(type='surface', bins=[1,2])
tally.addFilter(type='material', bins=[10,20])
tally.addFilter(type='universe', bins=[20,30])
tally.addFilter(type='energy', bins=np.linspace(1e-6, 2e6, 10))
tally.addFilter(type='energyout', bins=[1,2,3])
tally.addFilter(type='mesh', bins=1)
tallies_file.addTallySubelement(tally)

tally = openmc.Tally()
tally.addScore('flux')
tallies_file.addTallySubelement(tally)

tallies_file.exportToXML()