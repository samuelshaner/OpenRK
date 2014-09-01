from datasets.energy_groups import group_structures
import openmc
from openmc.statepoint import StatePoint
from infermc.process import XSTallyExtractor
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

  # Initialize an InferMC XSTallyExtractor object to compute cross-sections
  extractor = XSTallyExtractor(statepoint)

  ## MICROS
  extractor.extractAllMicroXS(groups, 'material')
  extractor.extractAllMicroXS(groups, 'distribcell')
  extractor.checkXS()

  materials = extractor._openmc_geometry.getAllMaterials()

  # DUMP-TO-FILE and PRINT XS
  for material in materials:
    for xs_type in xs_types:
      xs = extractor._multigroup_xs['material'][material._id][xs_type]
      xs.dumpToFile(filename='material-{0}-{1}'.format(material._id, xs_type))
      xs.printXS()
      xs.exportResults()
      xs.printPDF(filename='material-{0}-{1}'.format(material._id, xs_type))
      xs.checkXS()



  # RESTORE-FROM-FILE
  for material in materials:
    for xs_type in xs_types:

      if xs_type == 'total':
        xs = infermc.MicroTotalXS(material, 'material')
        xs.restoreFromFile(filename='material-{0}-{1}'.format(material._id, xs_type))
      elif xs_type == 'chi':
        xs = infermc.MicroChi(material, 'material')
        xs.restoreFromFile(filename='material-{0}-{1}'.format(material._id, xs_type))
      elif xs_type == 'transport':
        xs = infermc.MicroTransportXS(material, 'material')
        xs.restoreFromFile(filename='material-{0}-{1}'.format(material._id, xs_type))
      elif xs_type == 'scatter-matrix':
        xs = infermc.MicroScatterMatrixXS(material, 'material')
        xs.restoreFromFile(filename='material-{0}-{1}'.format(material._id, xs_type))

  '''
  ## MACROS
  extractor.extractAllMultiGroupXS(groups, 'material')
  extractor.extractAllMultiGroupXS(groups, 'distribcell')

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
      xs.dumpToFile(filename='material-{0}-{1}'.format(material._id, xs_type))
      xs.printXS()
      xs.exportResults()
      xs.printPDF(filename='material-{0}-{1}'.format(material._id, xs_type))
  '''

  openmc.reset_auto_ids()
  del extractor, statepoint