import openmc
from infermc.build import MicroXSTallyFactory
from datasets.energy_groups import group_structures


###############################################################################
#                      Simulation Input File Parameters
###############################################################################

# OpenMC simulation parameters
batches = 15
inactive = 5
particles = 10000


###############################################################################
#                 Exporting to OpenMC materials.xml File
###############################################################################

# Instantiate some Nuclides
h1 = openmc.Nuclide('H-1')
o16 = openmc.Nuclide('O-16')
u235 = openmc.Nuclide('U-235')
u238 = openmc.Nuclide('U-238')
zr90 = openmc.Nuclide('Zr-90')

# Instantiate a Material and register the Nuclides
inf_medium = openmc.Material(name='moderator')
inf_medium.set_density('g/cc', 5.)
inf_medium.add_nuclide(h1,  0.028999667)
inf_medium.add_nuclide(o16, 0.01450188)
inf_medium.add_nuclide(u235, 0.000114142)
inf_medium.add_nuclide(u238, 0.006886019)
inf_medium.add_nuclide(zr90, 0.002116053)

# Instantiate a MaterialsFile, register all Materials, and export to XML
materials_file = openmc.MaterialsFile()
materials_file.set_default_xs('71c')
materials_file.add_material(inf_medium)
materials_file.export_to_xml()


###############################################################################
#                 Exporting to OpenMC geometry.xml File
###############################################################################

# Instantiate boundary Planes
left = openmc.XPlane(bc_type='reflective', x0=-0.63)
right = openmc.XPlane(bc_type='reflective', x0=0.63)
bottom = openmc.YPlane(bc_type='reflective', y0=-0.63)
top = openmc.YPlane(bc_type='reflective', y0=0.63)

# Instantiate Cells
cell = openmc.Cell(cell_id=1, name='cell')

# Register Surfaces with Cells
cell.add_surface(surface=left, halfspace=+1)
cell.add_surface(surface=right, halfspace=-1)
cell.add_surface(surface=bottom, halfspace=+1)
cell.add_surface(surface=top, halfspace=-1)

# Register Material with Cells
cell.set_fill(inf_medium)

# Instantiate Universe
root = openmc.Universe(universe_id=0, name='root universe')

# Register Cell with Universe
root.add_cell(cell)

# Instantiate a Geometry and register the root Universe
geometry = openmc.Geometry()
geometry.set_root_universe(root)

# Instantiate a GeometryFile, register Geometry, and export to XML
geometry_file = openmc.GeometryFile()
geometry_file.set_geometry(geometry)
geometry_file.export_to_xml()


###############################################################################
#                   Exporting to OpenMC settings.xml File
###############################################################################

# Instantiate a SettingsFile, set all runtime parameters, and export to XML
settings_file = openmc.SettingsFile()
settings_file.set_batches(batches)
settings_file.set_inactive(inactive)
settings_file.set_particles(particles)
settings_file.set_source_space('box', [-0.63, -0.63, -0.63, 0.63, 0.63, 0.63])
settings_file.export_to_xml()


###############################################################################
##################   Exporting to OpenMC tallies.xml File  ####################
###############################################################################

groups = group_structures['CASMO']['2-group']

tally_factory = MicroXSTallyFactory(geometry)
tally_factory.createAllMultiGroupXS(groups, domain_type='material')
tally_factory.createTalliesFile()