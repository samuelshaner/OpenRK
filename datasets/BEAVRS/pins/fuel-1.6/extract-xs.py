from datasets.energy_groups import group_structures
import openmc
from openmc.statepoint import StatePoint
from infermc.process import XSTallyExtractor
from infermc.plotter import scatter_multigroup_xs
from infermc.multigroupxs import xs_types


groups = group_structures['CASMO']['2-group']

batches = [10, 15, 20, 25, 30]

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

    if xs_type != 'scatter matrix':
      scatter_multigroup_xs(extractor, xs_type,
                            domain_types=['distribcell', 'material'],
                            colors=['neighbors', 'material'],
                            filename='{0}-{1}-batches'.format(xs_type,batch))

  openmc.reset_auto_ids()
  del extractor, statepoint