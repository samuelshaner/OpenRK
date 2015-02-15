import openmc
from openmc.statepoint import StatePoint
from openmc.summary import Summary
from datasets.energy_groups import group_structures
from infermc.process import XSTallyExtractor, MicroXSTallyExtractor
from infermc.multigroupxs import xs_types
import infermc.plotter as plotter

batches = range(10, 35, 5)
batches = [10]

groups = group_structures['CASMO']['2-group']

summary = Summary('summary.h5')

for batch in batches:

  print batch

  filename = 'statepoint.{0:02}.h5'.format(batch)

  # Initialize a handle on the OpenMC statepoint file
  statepoint = StatePoint(filename)

  ## MICROS
  micro_extractor = MicroXSTallyExtractor(statepoint, summary)
  micro_extractor.extractAllMultiGroupXS(groups, 'cell')
  micro_extractor.extractAllMultiGroupXS(groups, 'distribcell')
#  micro_extractor.checkXS()

#  plotter.scatter_micro_xs(micro_extractor,
#                           domain_types=['distribcell', 'cell'],
#                           colors=['cell', 'cell'],
#                           filename='{0}-batch'.format(batch))

#  micro_extractor.buildNeighborMaps(unique=True)
#  plotter.scatter_all_neighbors(micro_extractor,
#                                filename='{0}-batch'.format(batch))

  cells = micro_extractor._openmc_geometry.get_all_material_cells()

  # DUMP-TO-FILE and PRINT XS
  for cell in cells:
    for xs_type in xs_types:
      xs = micro_extractor._multigroup_xs['cell'][cell._id][xs_type]
      xs.exportResults()
#      xs.printPDF(filename='cell-{0}-{1}'.format(cell._id, xs_type))

  '''
  micro_extractor.buildNeighborMaps(unique=True, first_level=1)

  test_xs = micro_extractor._multigroup_xs['distribcell'][10000]['fission']
  avg_xs = test_xs.getDomainAveragedXS()
  avg_xs.printPDF(filename='distribcell-10000-fission')

  all_neighbor_xs = test_xs.getAllNeighborAveragedXS()

  for neighbor in all_neighbor_xs:
    neighbor_xs = all_neighbor_xs[neighbor]
    neighbor_xs.printPDF(filename='neighbor-{0}-fission'.format(neighbor))


  # Plotting data colored by neighbors
  plotter.scatter_all_neighbors(micro_extractor, uncertainties=False,
                              filename='{0}-batch'.format(batch))

  '''
  openmc.reset_auto_ids()
  del micro_extractor, statepoint