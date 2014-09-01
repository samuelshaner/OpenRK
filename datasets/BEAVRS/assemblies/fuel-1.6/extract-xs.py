from datasets.energy_groups import group_structures
import openmc
from openmc.statepoint import StatePoint
from infermc.process import XSTallyExtractor, MicroXSTallyExtractor
from infermc.multigroupxs import xs_types
from infermc.plotter import scatter_multigroup_xs
import infermc


batches = range(10, 35, 5)
#batches = [10]

groups = group_structures['CASMO']['2-group']

for batch in batches:

  print batch

  filename = 'statepoint.{0}.h5'.format(batch)

  # Initialize a handle on the OpenMC statepoint file
  statepoint = openmc.statepoint.StatePoint(filename)

  ## MICROS
  micro_extractor = MicroXSTallyExtractor(statepoint)
  micro_extractor.extractAllMultiGroupXS(groups, 'material')
  micro_extractor.extractAllMultiGroupXS(groups, 'distribcell')
#  micro_extractor.checkXS()

  materials = micro_extractor._openmc_geometry.getAllMaterials()

  # DUMP-TO-FILE and PRINT XS
  for material in materials:
    for xs_type in xs_types:
      xs = micro_extractor._multigroup_xs['material'][material._id][xs_type]
      xs.dumpToFile(directory='micro', filename='material-{0}-{1}'.format(material._id, xs_type))
#      xs.printXS()
      xs.exportResults()
      xs.printPDF(directory='micro', filename='material-{0}-{1}'.format(material._id, xs_type))
      xs.checkXS()

  # RESTORE-FROM-FILE
  for material in materials:
    for xs_type in xs_types:

      if xs_type == 'total':
        xs = infermc.MicroTotalXS(material, 'material')
        xs.restoreFromFile(directory='micro', filename='material-{0}-{1}'.format(material._id, xs_type))
      elif xs_type == 'chi':
        xs = infermc.MicroChi(material, 'material')
        xs.restoreFromFile(directory='micro', filename='material-{0}-{1}'.format(material._id, xs_type))
      elif xs_type == 'transport':
        xs = infermc.MicroTransportXS(material, 'material')
        xs.restoreFromFile(directory='micro', filename='material-{0}-{1}'.format(material._id, xs_type))
      elif xs_type == 'scatter-matrix':
        xs = infermc.MicroScatterMatrixXS(material, 'material')
        xs.restoreFromFile(directory='micro', filename='material-{0}-{1}'.format(material._id, xs_type))

  '''
  openmc.reset_auto_ids()

  ## MACROS
  extractor = XSTallyExtractor(statepoint)
  extractor.extractAllMultiGroupXS(groups, 'material')
  extractor.extractAllMultiGroupXS(groups, 'distribcell')
  extractor.checkXS()

  for xs_type in xs_types:

    if xs_type != 'scatter matrix':
      scatter_multigroup_xs(extractor, xs_type,
                            domain_types=['distribcell', 'material'],
                            colors=['neighbors', 'material'],
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
  '''

  openmc.reset_auto_ids()
  del micro_extractor, statepoint