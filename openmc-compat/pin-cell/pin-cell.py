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
left = XPlane(x=-4.0)
right = XPlane(x=4.0)
top = YPlane(y=4.0)
bottom = YPlane(y=-4.0)

left.setBoundaryType(REFLECTIVE)
right.setBoundaryType(REFLECTIVE)
top.setBoundaryType(REFLECTIVE)
bottom.setBoundaryType(REFLECTIVE)


###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

cells = []
cells.append(CellBasic(universe=1, material=10008, rings=2, sectors=4))
cells.append(CellBasic(universe=1, material=10007))
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

lattice = Lattice(id=2, width_x=4.0, width_y=4.0)
lattice.setLatticeCells([[1,1],
                         [1,1]])


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()

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

uo2 = IsoMaterial(material_id(), 'UO2')
uo2.setDensity(10., 'g/cc')
uo2.addIsotope(o16, 2.0)
uo2.addIsotope(u234, 0.001)
uo2.addIsotope(u235, 0.04)
uo2.addIsotope(u238, 0.959)

geometry.addMaterial(uo2)
geometry.addMaterial(h2o)

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
settings_file.setBatches(25)
settings_file.setParticles(10000)
settings_file.createEigenvalueSubelement()
#settings_file.createOutputSubelement(['tallies'])
settings_file.createStatepointSubelement(batches=[1,2,3,4,5])
settings_file.createSourceSpaceSubelement(type='box', params=[-1,-1,-1,1,1,1])
settings_file.exportToXML()

plot_file = openmc.PlotsFile()
plot_file.addNewPlot(id=1, width=[8.0,8.0], origin=[0.,0.,0.], pixels=[1000,1000])
plot_file.exportToXML()

materials_file = openmc.MaterialsFile()
materials_file.createDefaultXSSubelement()
materials_file.createMaterialSubelement(uo2)
materials_file.createMaterialSubelement(h2o)
materials_file.exportToXML()

#tallies_file = openmc.TalliesFile()
#tally = openmc.Tally()
#tally.setLabel(label='first tally')
#tally.addNuclide('total')
#tally.addScore('flux')
#tally.addScore('total')
#tally.addScore('fission')
#tally.addFilter(type='cell', bins=[10001])
#tally.addFilter(type='energy', bins=np.linspace(1e-6, 2e6, 10))
#tallies_file.addTallySubelement(tally)
#tallies_file.exportToXML()