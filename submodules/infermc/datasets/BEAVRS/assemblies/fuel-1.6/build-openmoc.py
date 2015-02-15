from openmoc.compatible.opencg_compatible import *
from geometry import geometry
import opencg, openmoc
import openmoc.plotter
import opencg.plotter
import h5py, numpy


###############################################################################
########################   Creating OpenMOC Geometry  #########################
###############################################################################

openmoc_geometry = get_openmoc_geometry(geometry)
openmoc_geometry.initializeFlatSourceRegions()

cells = openmoc_geometry.getRootUniverse().getAllCells()

num_groups = 2

f = h5py.File('multigroupxs/multigroupxs.h5', 'r')

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

openmoc.plotter.plot_cells(openmoc_geometry)
openmoc.plotter.plot_materials(openmoc_geometry)

# Must initialize FSRs before track generation
openmoc_geometry.initializeFlatSourceRegions()
track_generator = openmoc.TrackGenerator(openmoc_geometry, 32, 0.1)
track_generator.generateTracks()

solver = openmoc.CPUSolver(openmoc_geometry, track_generator)
solver.setSourceConvergenceThreshold(1E-5)
solver.convergeSource(1000)
solver.printTimerReport()
