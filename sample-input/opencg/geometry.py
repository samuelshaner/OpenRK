import opencg as cg
import openmc as mc

###############################################################################
######################   Creating Bounding Surfaces   #########################
###############################################################################

print('Creating the bounding Surfaces...')

boundaries = dict()

boundaries['Box'] = cg.ZSquarePrism(boundary='vacuum', R=50.)
boundaries['Z-Min'] = cg.ZPlane(z0=-50., boundary='vacuum')
boundaries['Z-Max'] = cg.ZPlane(z0=50., boundary='vacuum')


###############################################################################
#########################   Discretize Pin Cells  #############################
###############################################################################

print('Discretizing the Cells...')

# create material
mat = mc.Material(name='3.1% Fuel')
mat.set_density('g/cm3', 10.30166)
mat.add_nuclide(mc.Nuclide('U-234'), 5.7987e-6)
mat.add_nuclide(mc.Nuclide('U-235'), 7.2175e-4)
mat.add_nuclide(mc.Nuclide('U-238'), 2.2253e-2)
mat.add_nuclide(mc.Nuclide('O-16'), 4.5940e-2)

universe = cg.Universe(universe_id=0, name='Root Universe')
opencg_mat = mc.opencg_compatible.get_opencg_material(mat)
cell = cg.Cell(name='Root Cell', fill=opencg_mat)
cell.addSurface(surface=boundaries['Box'], halfspace=-1)
cell.addSurface(surface=boundaries['Z-Min'], halfspace=+1)
cell.addSurface(surface=boundaries['Z-Max'], halfspace=-1)
universe.addCell(cell)

###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

print('Creating the Geometry...')

geometry = cg.Geometry()
geometry.setRootUniverse(universe)

geometry.initializeCellOffsets()

num_regions = geometry._num_regions

print('# regions = %d' % num_regions)
