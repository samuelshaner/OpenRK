import openmc
import openmc.statepoint
from datasets.energy_groups import group_structures
from infermc.process import XSTallyExtractor, MicroXSTallyExtractor
from infermc.multigroupxs import xs_types
import infermc.plotter as plotter


#batches = range(10, 35, 5)
batches = [30]

groups = group_structures['CASMO']['2-group']

for batch in batches:

  print batch

  filename = 'statepoint.{0}.h5'.format(batch)

  # Initialize a handle on the OpenMC statepoint file
  statepoint = openmc.statepoint.StatePoint(filename)

  ## MICROS
  micro_extractor = MicroXSTallyExtractor(statepoint)
  micro_extractor.extractAllMultiGroupXS(groups, 'material', corr=True)
  micro_extractor.extractAllMultiGroupXS(groups, 'distribcell', corr=True)
  micro_extractor.checkXS()

#  plotter.scatter_micro_xs(micro_extractor,
#                           domain_types=['distribcell', 'material'],
#                           colors=['cell', 'material'],
#                           filename='{0}-batch'.format(batch))

  materials = micro_extractor._openmc_geometry.getAllMaterials()

  # DUMP-TO-FILE and PRINT XS
  for material in materials:
    for xs_type in xs_types:
      xs = micro_extractor._multigroup_xs['material'][material._id][xs_type]
#      xs.dumpToFile(directory='micro', filename='material-{0}-{1}'.format(material._id, xs_type))
#      xs.exportResults()
      xs.printPDF(directory='micro', filename='material-{0}-{1}'.format(material._id, xs_type))

  openmc.reset_auto_ids()
  del micro_extractor, statepoint

  '''
  ## MACROS
  extractor = XSTallyExtractor(statepoint)
  extractor.extractAllMultiGroupXS(groups, 'material')
  extractor.extractAllMultiGroupXS(groups, 'distribcell')
  extractor.checkXS()

  plotter.scatter_micro_xs(extractor,
                           domain_types=['distribcell', 'material'],
                           filename='{0}-batch'.format(batch))

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
