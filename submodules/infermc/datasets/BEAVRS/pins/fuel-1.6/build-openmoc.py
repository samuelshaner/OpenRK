from openmoc.compatible.opencg_compatible import *
from geometry import geometry
import openmoc
import h5py, numpy


###############################################################################
########################   Creating OpenMOC Geometry  #########################
###############################################################################

water_mesh = opencg.RadialMesh()
water_mesh.setNumRings(3)
water_mesh.setMaxRadius(0.7)
water_mesh.setMinRadius(0.4572)
water_mesh.setWithOuter(True)

fuel_mesh = opencg.RadialMesh()
fuel_mesh.setNumRings(5)
fuel_mesh.setMaxRadius(0.39218)
fuel_mesh.setMinRadius(0.)
fuel_mesh.setWithOuter(False)

mesh = opencg.SectorMesh(num_sectors=8)

universes = geometry.getAllMaterialUniverses()

cells = geometry._root_universe.getAllCells()
for cell_id, cell in cells.items():
  if cell._type == 'material':
    if cell._fill._id == 10003:
      new_cells = water_mesh.subdivideCell(cell=cell, universe=universes[10000])
    if cell._fill._id == 10000:
      new_cells = fuel_mesh.subdivideCell(cell=cell, universe=universes[10000])

mesh.subdivideUniverse(universe=universes[10000])

openmoc_geometry = get_openmoc_geometry(geometry)

cells = openmoc_geometry.getRootUniverse().getAllCells()

num_groups = 70

f = h5py.File('multigroupxs/multigroupxs-70-group.h5', 'r')

mats = f['material']

for mat_key in mats.keys():
  mat_id = int(mat_key.split(' ')[1])
  mat_group = mats[mat_key]

  openmoc_material = openmoc.Material(mat_id)
  openmoc_material.setNumEnergyGroups(num_groups)
  openmoc_material.thisown = 0   # FIXME: Can SWIG do this on its own??

  macro_xs = dict()
  macro_xs['total'] = numpy.zeros(num_groups)
  macro_xs['transport'] = numpy.zeros(num_groups)
  macro_xs['scatter matrix'] = numpy.zeros((num_groups, num_groups))
  macro_xs['absorption'] = numpy.zeros(num_groups)
  macro_xs['fission'] = numpy.zeros(num_groups)
  macro_xs['nu-fission'] = numpy.zeros(num_groups)
  macro_xs['chi'] = numpy.zeros(num_groups)

  for nuclide in mat_group.keys():
    nuclide_group = mat_group[nuclide]
    density = nuclide_group['density'][0]

    if nuclide == 'total':
      continue

    for rxn_type in nuclide_group.keys():
      if rxn_type in macro_xs.keys():
        micro_xs = nuclide_group[rxn_type]['average'][...]
        macro_xs[rxn_type] += micro_xs * density

  if mat_id != 10004:
    openmoc_material.setSigmaT(macro_xs['transport'])
    openmoc_material.setSigmaA(macro_xs['absorption'])
    openmoc_material.setSigmaF(macro_xs['fission'])
    openmoc_material.setNuSigmaF(macro_xs['nu-fission'])
    openmoc_material.setSigmaS(macro_xs['scatter matrix'].ravel())
    openmoc_material.setChi(macro_xs['chi'])
  else:
    openmoc_material.setSigmaT(macro_xs['total'])
    openmoc_material.setSigmaA(macro_xs['absorption'])
    openmoc_material.setSigmaF(macro_xs['fission'])
    openmoc_material.setNuSigmaF(macro_xs['nu-fission'])
    openmoc_material.setSigmaS(macro_xs['scatter matrix'].ravel())
    openmoc_material.setChi(macro_xs['chi'])

  cells = openmoc_geometry.getRootUniverse().getAllCells()

  for cell_id, cell in cells.items():
    if cell.getType() == openmoc.MATERIAL:
      cell = openmoc.castCellToCellBasic(cell)
      if cell.getMaterial().getId() == mat_id:
        cell.setMaterial(openmoc_material)


# FIXME: The ring/sectors in OpenMOC are not working
#cells = openmoc_geometry.getRootUniverse().getAllCells()
#for cell_id, cell in cells.items():
#  if cell.getType() == openmoc.MATERIAL:
#    cell = openmoc.castCellToCellBasic(cell)
#    if cell.getMaterial().getId() == 10000:
#      cell.setNumRings(3)
    #  cell.setNumSectors(8)


# Must initialize FSRs before track generation
openmoc_geometry.initializeFlatSourceRegions()
track_generator = openmoc.TrackGenerator(openmoc_geometry, 128, 0.05)
track_generator.generateTracks()

solver = openmoc.CPUSolver(openmoc_geometry, track_generator)
solver.setSourceConvergenceThreshold(1E-5)
solver.convergeSource(1000)
solver.printTimerReport()