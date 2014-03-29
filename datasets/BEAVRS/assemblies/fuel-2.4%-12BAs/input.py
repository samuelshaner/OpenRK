import openmoc.log as log
import openmoc.compatible.openmc as openmc
from lattices import *
import materials
import surfaces
import pincells


###############################################################################
###################   Simulation Input File Parameters   ######################
###############################################################################

# Get the appropriate lattice from the lattices module
lattice = lattices['2.4% Fuel - 12BA']
lattice_id = lattice.getId()

# Height of the axial slice
slice_height = 10.

# OpenMC simulation parameters
batches = 25
inactive = 10
particles = 10000


###############################################################################
######################   Creating Bounding Surfaces   #########################
###############################################################################

log.py_printf('NORMAL', 'Creating the bounding Surfaces...')

boundaries = dict()

boundaries['X-Min'] = XPlane(x=-lattice_width / 2.)
boundaries['X-Max'] = XPlane(x=lattice_width / 2.)
boundaries['Y-Min'] = YPlane(y=-lattice_width / 2.)
boundaries['Y-Max'] = YPlane(y=lattice_width / 2.)
boundaries['Z-Min'] = ZPlane(z=-slice_height / 2.)
boundaries['Z-Max'] = ZPlane(z=slice_height / 2.)

for index in boundaries.keys():
  boundaries[index].setBoundaryType(REFLECTIVE)
  surfaces.surfaces[index] = boundaries[index]


###############################################################################
######################   Creating Root Universe   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating the root Universe...')

# Root cell encapsulates the full geometry
root = CellFill(universe=0, universe_fill=lattice_id)
root.addSurface(halfspace=+1, surface=boundaries['X-Min'])
root.addSurface(halfspace=-1, surface=boundaries['X-Max'])
root.addSurface(halfspace=+1, surface=boundaries['Y-Min'])
root.addSurface(halfspace=-1, surface=boundaries['Y-Max'])
root.addSurface(halfspace=+1, surface=boundaries['Z-Min'])
root.addSurface(halfspace=-1, surface=boundaries['Z-Max'])

# Append the root to the list of all Cells
pincells.cells.append(root)


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating the Geometry...')

geometry = Geometry()

# Add all Materials to the Geometry
for material in materials.materials.values():
  geometry.addMaterial(material)

# Add all Surfaces to the Geometry
for surface in surfaces.surfaces.values():
  geometry.addSurface(surface)

# Add all Cells to the Geometry
for cell in pincells.cells:
  geometry.addCell(cell)

# Add lattice to
geometry.addLattice(lattice)



###############################################################################
###################   Exporting to OpenMC XML Input Files  ####################
###############################################################################

log.py_printf('NORMAL', 'Exporting to OpenMC XML Files...')

# settings.xml
settings_file = openmc.SettingsFile()
settings_file.setBatches(batches)
settings_file.setInactive(inactive)
settings_file.setParticles(particles)
settings_file.createEigenvalueSubelement()

source = [-pin_pitch/2., -pin_pitch/2., -slice_height/2.,
          pin_pitch/2., pin_pitch/2., slice_height/2.]
settings_file.createSourceSpaceSubelement(type='box', params=source)

settings_file.exportToXML()

# plots.xml
plot_file = openmc.PlotsFile()
plot_file.addNewPlot(id=1, width=[lattice_width, lattice_width],
                     origin=[0.,0.,0.], pixels=[1000,1000])
plot_file.exportToXML()

# materials.xml
materials_file = openmc.MaterialsFile()
materials_file.createDefaultXSSubelement()
materials_file.createMaterialsSubelements(materials.materials.values())
materials_file.exportToXML()

# geometry.xml
geometry_file = openmc.GeometryFile()
geometry_file.createGeometrySubelements(geometry)
geometry_file.exportToXML()