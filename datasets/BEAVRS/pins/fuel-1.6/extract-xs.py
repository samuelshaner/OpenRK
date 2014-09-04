from datasets.energy_groups import group_structures
import openmc
from openmc.statepoint import StatePoint
from infermc.process import XSTallyExtractor, MicroXSTallyExtractor
from infermc.multigroupxs import xs_types
import infermc.plotter as plotter
import infermc


#batches = range(25, 105, 5)
batches = [30]

groups = group_structures['CASMO']['2-group']

for batch in batches:

  print batch

  filename = 'statepoint.{0:02d}.h5'.format(batch)

  # Initialize a handle on the OpenMC statepoint file
  statepoint = openmc.statepoint.StatePoint(filename)

  ## MICROS
  micro_extractor = MicroXSTallyExtractor(statepoint)
  micro_extractor.extractAllMultiGroupXS(groups, 'material')
  micro_extractor.extractAllMultiGroupXS(groups, 'distribcell')
  micro_extractor.checkXS()

  nuclides = micro_extractor._openmc_geometry.getAllNuclides()

  for xs_type in xs_types:

    if xs_type != 'scatter matrix':

      for nuclide_name, nuclide_tuple in nuclides.items():
        plotter.scatter_micro_xs(micro_extractor, xs_type, nuclide_tuple[0],
                              domain_types=['distribcell', 'material'],
                              filename='{0}-{1}-{2}-batches'.format(nuclide_name, xs_type, batch))

  materials = micro_extractor._openmc_geometry.getAllMaterials()

  # DUMP-TO-FILE and PRINT XS
  for material in materials:
    for xs_type in xs_types:
      xs = micro_extractor._multigroup_xs['material'][material._id][xs_type]
      xs.printPDF(directory='micro', filename='material-{0}-{1}'.format(material._id, xs_type))

  openmc.reset_auto_ids()
  del micro_extractor, statepoint

  '''
  ## MACROS
  extractor = XSTallyExtractor(statepoint)
  extractor.extractAllMultiGroupXS(groups, 'material')
  extractor.extractAllMultiGroupXS(groups, 'distribcell')
  extractor.checkXS()

  for xs_type in xs_types:

    if xs_type != 'scatter matrix':
      plotter.scatter_multigroup_xs(extractor, xs_type,
                            domain_types=['distribcell', 'material'],
                            colors=['cell', 'material'],
                            filename='{0}-{1}-batches'.format(xs_type,batch))

  materials = extractor._openmc_geometry.getAllMaterials()

  # DUMP-TO-FILE and PRINT XS
  for material in materials:
    for xs_type in xs_types:
      xs = extractor._multigroup_xs['material'][material._id][xs_type]
      xs.dumpToFile(directory='macro', filename='material-{0}-{1}'.format(material._id, xs_type))
      xs.printXS()
      xs.exportResults()
      xs.printPDF(directory='macro', filename='material-{0}-{1}'.format(material._id, xs_type))

  openmc.reset_auto_ids()
  del extractor, statepoint
  '''