from datasets.energy_groups import group_structures
import openmc
from openmc.statepoint import StatePoint
from infermc.process import XSTallyExtractor
from infermc.multigroupxs import xs_types
from infermc.plotter import scatter_multigroup_xs
import numpy as np


batches = range(10, 55, 5)

groups = group_structures['CASMO']['2-group']

for batch in batches:

  print batch

  filename = 'statepoint.{0}.h5'.format(batch)

  # Initialize a handle on the OpenMC statepoint file
  statepoint = openmc.statepoint.StatePoint(filename)

  # Initialize an InferMC XSTallyExtractor object to compute cross-sections
  extractor = XSTallyExtractor(statepoint)

  extractor.extractAllMultiGroupXS(groups, 'material')
  extractor.extractAllMultiGroupXS(groups, 'distribcell')

  for xs_type in xs_types:

    cells = extractor._opencsg_geometry.getAllMaterialCells()

#    for cell in cells:

#      xs = extractor.getMultiGroupXS(xs_type, cell._id, 'distribcell')
#      filename = 'cell-{0}-{1}'.format(cell._id, xs_type)
#      xs.dumpToFile(filename)

    if xs_type != 'scatter matrix':
      scatter_multigroup_xs(extractor, xs_type,
                            domain_types=['distribcell', 'material'],
                            colors=['unique neighbors', 'material'], extension='png',
                            filename='{0}-{1}-batches'.format(xs_type,batch))

  openmc.reset_auto_ids()