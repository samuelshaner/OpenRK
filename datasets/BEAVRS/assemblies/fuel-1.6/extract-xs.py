from datasets.energy_groups import group_structures
import openmc
from openmc.statepoint import StatePoint
from infermc.process import XSTallyExtractor
from infermc.multigroupxs import xs_types
from infermc.plotter import scatter_multigroup_xs
import infermc


#batches = range(10, 35, 5)
batches = [10]

groups = group_structures['CASMO']['2-group']

for batch in batches:

  print batch

  filename = 'statepoint.{0}.h5'.format(batch)

  # Initialize a handle on the OpenMC statepoint file
  statepoint = openmc.statepoint.StatePoint(filename)

  # Initialize an InferMC XSTallyExtractor object to compute cross-sections
  extractor = XSTallyExtractor(statepoint)

  extractor.extractAllMicroXS(groups, 'material')
  extractor.extractAllMicroXS(groups, 'distribcell')

  domains = extractor._openmc_geometry.getAllMaterialCells()
  domain_type = 'distribcell'


  # DUMP-TO-FILE and PRINT XS
  for domain in domains:
    for xs_type in xs_types:
      xs = extractor._multigroup_xs[domain_type][domain._id][xs_type]
      xs.dumpToFile(filename='material-{0}-{1}'.format(domain._id, xs_type))
      xs.printXS(nuclide='all')


  # RESTORE-FROM-FILE
  for domain in domains:
    for xs_type in xs_types:

      if xs_type == 'total':
        xs = infermc.MicroTotalXS(domain, domain_type)
        xs.restoreFromFile(filename='material-{0}-{1}'.format(domain._id, xs_type))
      elif xs_type == 'chi':
        xs = infermc.MicroChi(domain, domain_type)
        xs.restoreFromFile(filename='material-{0}-{1}'.format(domain._id, xs_type))
      elif xs_type == 'transport':
        xs = infermc.MicroTransportXS(domain, domain_type)
        xs.restoreFromFile(filename='material-{0}-{1}'.format(domain._id, xs_type))
      elif xs_type == 'scatter-matrix':
        xs = infermc.MicroScatterMatrixXS(domain, domain_type)
        xs.restoreFromFile(filename='material-{0}-{1}'.format(domain._id, xs_type))


  '''
  extractor.extractAllMultiGroupXS(groups, 'material')
  extractor.extractAllMultiGroupXS(groups, 'distribcell')
  '''

  '''
  for xs_type in xs_types:

    cells = extractor._opencsg_geometry.getAllMaterialCells()

    if xs_type != 'scatter matrix' and xs_type != 'transport':
      print xs_type
      scatter_multigroup_xs(extractor, xs_type,
                            domain_types=['distribcell', 'material'],
                            energy_groups=(1,2),
                            colors=['neighbors', 'material'], extension='png',
                            filename='{0}-{1}-batches'.format(xs_type, batch))
  '''

  openmc.reset_auto_ids()
