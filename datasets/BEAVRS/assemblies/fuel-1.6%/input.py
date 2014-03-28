from openmoc import *
import openmoc.log as log
import openmoc.compatible.openmc as openmc



###############################################################################
############################   Creating Isotopes   ############################
###############################################################################

log.py_printf('NORMAL', 'Creating Isotopes...')

# 1.6% enriched uranium fuel
u234 = openmc.Isotope(isotope_id(), 'U-234')
u235 = openmc.Isotope(isotope_id(), 'U-235')
u238 = openmc.Isotope(isotope_id(), 'U-238')
o16 = openmc.Isotope(isotope_id(), 'O-16')

# Borated water
b10 = openmc.Isotope(isotope_id(), 'B-10')
b11 = openmc.Isotope(isotope_id(), 'B-11')
h1 = openmc.Isotope(isotope_id(), 'H-1')
h2 = openmc.Isotope(isotope_id(), 'H-2')

# Helium
he4 = openmc.Isotope(isotope_id(), 'He-4')

# Zircaloy
cr50 = openmc.Isotope(isotope_id(), 'Cr-50')
cr52 = openmc.Isotope(isotope_id(), 'Cr-52')
cr53 = openmc.Isotope(isotope_id(), 'Cr-53')
cr54 = openmc.Isotope(isotope_id(), 'Cr-54')
fe54 = openmc.Isotope(isotope_id(), 'Fe-54')
fe56 = openmc.Isotope(isotope_id(), 'Fe-56')
fe57 = openmc.Isotope(isotope_id(), 'Fe-57')
fe58 = openmc.Isotope(isotope_id(), 'Fe-58')
zr90 = openmc.Isotope(isotope_id(), 'Zr-90')
zr91 = openmc.Isotope(isotope_id(), 'Zr-91')
zr92 = openmc.Isotope(isotope_id(), 'Zr-92')
zr94 = openmc.Isotope(isotope_id(), 'Zr-94')
zr96 = openmc.Isotope(isotope_id(), 'Zr-96')
sn112 = openmc.Isotope(isotope_id(), 'Sn-112')
sn114 = openmc.Isotope(isotope_id(), 'Sn-114')
sn115 = openmc.Isotope(isotope_id(), 'Sn-115')
sn116 = openmc.Isotope(isotope_id(), 'Sn-116')
sn117 = openmc.Isotope(isotope_id(), 'Sn-117')
sn118 = openmc.Isotope(isotope_id(), 'Sn-118')
sn119 = openmc.Isotope(isotope_id(), 'Sn-119')
sn120 = openmc.Isotope(isotope_id(), 'Sn-120')
sn122 = openmc.Isotope(isotope_id(), 'Sn-122')
sn124 = openmc.Isotope(isotope_id(), 'Sn-124')

u234.setNumEnergyGroups(1)
u235.setNumEnergyGroups(1)
u238.setNumEnergyGroups(1)
o16.setNumEnergyGroups(1)
b10.setNumEnergyGroups(1)
b11.setNumEnergyGroups(1)
h1.setNumEnergyGroups(1)
h2.setNumEnergyGroups(1)
he4.setNumEnergyGroups(1)
cr50.setNumEnergyGroups(1)
cr52.setNumEnergyGroups(1)
cr53.setNumEnergyGroups(1)
cr54.setNumEnergyGroups(1)
fe54.setNumEnergyGroups(1)
fe56.setNumEnergyGroups(1)
fe57.setNumEnergyGroups(1)
fe58.setNumEnergyGroups(1)
zr90.setNumEnergyGroups(1)
zr91.setNumEnergyGroups(1)
zr92.setNumEnergyGroups(1)
zr94.setNumEnergyGroups(1)
zr96.setNumEnergyGroups(1)
sn112.setNumEnergyGroups(1)
sn114.setNumEnergyGroups(1)
sn115.setNumEnergyGroups(1)
sn116.setNumEnergyGroups(1)
sn117.setNumEnergyGroups(1)
sn118.setNumEnergyGroups(1)
sn119.setNumEnergyGroups(1)
sn120.setNumEnergyGroups(1)
sn122.setNumEnergyGroups(1)
sn124.setNumEnergyGroups(1)



###############################################################################
###########################   Creating Materials   ############################
###############################################################################

log.py_printf('NORMAL', 'Creating Materials...')

fuel = openmc.IsoMaterial(material_id(), '1.6% Enr. Fuel')
fuel.setNumEnergyGroups(1)
fuel.setDensity(10.31341, 'g/cc')
fuel.addIsotope(u234, 3.0131e-6)
fuel.addIsotope(u235, 3.7503e-4)
fuel.addIsotope(u238, 2.2625e-2)
fuel.addIsotope(o16, 4.6007e-2)    # Includes O-17 and O-18 number densities

moderator = openmc.IsoMaterial(material_id(), 'Moderator')
moderator.setNumEnergyGroups(1)
moderator.setDensity(0.740582, 'g/cc')
moderator.addIsotope(b10, 8.0042e-6)
moderator.addIsotope(b11, 3.2218e-5)
moderator.addIsotope(h1, 4.9457e-2)
moderator.addIsotope(h2, 7.4196e-6)
moderator.addIsotope(o16, 2.4732e-2)  # Includes O-17 and O-18 number densities

gap = openmc.IsoMaterial(material_id(), 'Gap')
gap.setNumEnergyGroups(1)
gap.setDensity(0.001598, 'g/cc')
gap.addIsotope(he4, 2.4044e-4)

clad = openmc.IsoMaterial(material_id(), 'Clad')
clad.setNumEnergyGroups(1)
clad.setDensity(6.55, 'g/cc')
clad.addIsotope(o16, 3.0818e-4)   # Includes O-17 and O-18 number densities
clad.addIsotope(cr50, 3.2962e-6)
clad.addIsotope(cr52, 6.3564e-5)
clad.addIsotope(cr53, 7.2076e-6)
clad.addIsotope(cr54, 1.7941e-6)
clad.addIsotope(fe54, 8.6699e-6)
clad.addIsotope(fe56, 1.3610e-4)
clad.addIsotope(fe57, 3.1431e-6)
clad.addIsotope(fe58, 4.1829e-7)
clad.addIsotope(zr90, 2.1827e-2)
clad.addIsotope(zr91, 4.7600e-3)
clad.addIsotope(zr92, 7.2758e-3)
clad.addIsotope(zr94, 7.3734e-3)
clad.addIsotope(zr96, 1.1879e-3)
clad.addIsotope(sn112, 4.6735e-6)
clad.addIsotope(sn114, 3.1799e-6)
clad.addIsotope(sn115, 1.6381e-6)
clad.addIsotope(sn116, 7.0055e-5)
clad.addIsotope(sn117, 3.7003e-5)
clad.addIsotope(sn118, 1.1669e-4)
clad.addIsotope(sn119, 4.1387e-5)
clad.addIsotope(sn120, 1.5697e-4)
clad.addIsotope(sn122, 2.2308e-5)
clad.addIsotope(sn124, 2.7897e-5)



###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating Surfaces...')

pitch = 1.25984                                            # centimeters
lattice_pin_dims = 17.                                     # 17 pins
lattice_width = lattice_pin_dims * pitch                   # Lattice width [cm]
slice_height = 10.

surfaces = {}

surfaces['Fuel Radius'] = Circle(x=0.0, y=0.0, radius=0.39218)
surfaces['Gap Radius'] = Circle(x=0.0, y=0.0, radius=0.40005)
surfaces['Fuel Clad Radius'] = Circle(x=0.0, y=0.0, radius=0.45720)
surfaces['GT Radius'] = Circle(x=0.0, y=0.0, radius=0.50419)
surfaces['GT Clad Radius'] = Circle(x=0.0, y=0.0, radius=0.54610)

surfaces['Left Boundary'] = XPlane(x=-lattice_width / 2.)
surfaces['Right Boundary'] = XPlane(x=lattice_width / 2.)
surfaces['Top Boundary'] = YPlane(y=lattice_width / 2.)
surfaces['Bottom Boundary'] = YPlane(y=-lattice_width / 2.)
surfaces['Upper Boundary'] = ZPlane(z=slice_height / 2.)
surfaces['Lower Boundary'] = ZPlane(z=-slice_height / 2.)

surfaces['Left Boundary'].setBoundaryType(REFLECTIVE)
surfaces['Right Boundary'].setBoundaryType(REFLECTIVE)
surfaces['Top Boundary'].setBoundaryType(REFLECTIVE)
surfaces['Bottom Boundary'].setBoundaryType(REFLECTIVE)
surfaces['Upper Boundary'].setBoundaryType(REFLECTIVE)
surfaces['Lower Boundary'].setBoundaryType(REFLECTIVE)



###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

log.py_printf('NORMAL', 'Creating Cells...')

universes = {}

# Fuel pin
universes[1] = []

fuel_region = CellBasic(universe=1, material=fuel.getId(), rings=3, sectors=4)
fuel_region.addSurface(halfspace=-1, surface=surfaces['Fuel Radius'])

gap_region = CellBasic(universe=1, material=gap.getId(), sectors=4)
gap_region.addSurface(halfspace=+1, surface=surfaces['Fuel Radius'])
gap_region.addSurface(halfspace=-1, surface=surfaces['Gap Radius'])

clad_region = CellBasic(universe=1, material=clad.getId(), sectors=4)
clad_region.addSurface(halfspace=+1, surface=surfaces['Gap Radius'])
clad_region.addSurface(halfspace=-1, surface=surfaces['Fuel Clad Radius'])

moderator_region = CellBasic(universe=1, material=moderator.getId(), sectors=4)
moderator_region.addSurface(halfspace=+1, surface=surfaces['Fuel Clad Radius'])

universes[1].append(fuel_region)
universes[1].append(gap_region)
universes[1].append(clad_region)
universes[1].append(moderator_region)


# Empty guide tube
universes[2] = []

inner_water_region = CellBasic(universe=2, material=moderator.getId(),
                               rings=3, sectors=4)
inner_water_region.addSurface(halfspace=-1, surface=surfaces['GT Radius'])

clad_region = CellBasic(universe=2, material=clad.getId(), sectors=4)
clad_region.addSurface(halfspace=+1, surface=surfaces['GT Radius'])
clad_region.addSurface(halfspace=-1, surface=surfaces['GT Clad Radius'])

outer_water_region = CellBasic(universe=2, material=moderator.getId(), sectors=4)
outer_water_region.addSurface(halfspace=+1, surface=surfaces['GT Clad Radius'])

universes[2].append(inner_water_region)
universes[2].append(clad_region)
universes[2].append(outer_water_region)


# Root cell - full geometry
universes[0] = []
root = CellFill(universe=0, universe_fill=10)
root.addSurface(halfspace=+1, surface=surfaces['Left Boundary'])
root.addSurface(halfspace=-1, surface=surfaces['Right Boundary'])
root.addSurface(halfspace=+1, surface=surfaces['Bottom Boundary'])
root.addSurface(halfspace=-1, surface=surfaces['Top Boundary'])
root.addSurface(halfspace=+1, surface=surfaces['Lower Boundary'])
root.addSurface(halfspace=-1, surface=surfaces['Upper Boundary'])
universes[0].append(root)



###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating Lattices...')

lattice = Lattice(id=10, width_x=pitch, width_y=pitch)
lattice.setLatticeCells([[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                         [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                         [1,1,1,1,1,2,1,1,2,1,1,2,1,1,1,1,1],
                         [1,1,1,2,1,1,1,1,1,1,1,1,1,2,1,1,1],
                         [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                         [1,1,2,1,1,2,1,1,2,1,1,2,1,1,2,1,1],
                         [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                         [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                         [1,1,2,1,1,2,1,1,2,1,1,2,1,1,2,1,1],
                         [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                         [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                         [1,1,2,1,1,2,1,1,2,1,1,2,1,1,2,1,1],
                         [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                         [1,1,1,2,1,1,1,1,1,1,1,1,1,2,1,1,1],
                         [1,1,1,1,1,2,1,1,2,1,1,2,1,1,1,1,1],
                         [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                         [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]])



###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()
geometry.addMaterial(fuel)
geometry.addMaterial(moderator)
geometry.addMaterial(gap)
geometry.addMaterial(clad)

for cells in universes.values():
  for cell in cells:
    geometry.addCell(cell)

geometry.addLattice(lattice)



###############################################################################
###################   Exporting to OpenMC XML Input Files  ####################
###############################################################################

log.py_printf('NORMAL', 'Exporting to OpenMC XML Files...')

# settings.xml
settings_file = openmc.SettingsFile()
settings_file.setBatches(25)
settings_file.setParticles(10000)
settings_file.createEigenvalueSubelement()
settings_file.createSourceSpaceSubelement(type='box',
                                          params=[-pitch/2., -pitch/2.,
                                                  -slice_height/2., pitch/2.,
                                                  pitch/2., slice_height/2.])
settings_file.exportToXML()

# plots.xml
plot_file = openmc.PlotsFile()
plot_file.addNewPlot(id=1, width=[lattice_width, lattice_width],
                     origin=[0.,0.,0.], pixels=[1000,1000])
plot_file.exportToXML()

# materials.xml
materials_file = openmc.MaterialsFile()
materials_file.createDefaultXSSubelement()
materials_file.createMaterialSubelement(fuel)
materials_file.createMaterialSubelement(moderator)
materials_file.createMaterialSubelement(gap)
materials_file.createMaterialSubelement(clad)
materials_file.exportToXML()

# geometry.xml
geometry_file = openmc.GeometryFile()
geometry_file.createGeometrySubelements(geometry)
geometry_file.exportToXML()