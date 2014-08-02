from datasets.energy_groups import group_structures
import openmc
from openmc.statepoint import StatePoint
from infermc.process import XSTallyExtractor
from infermc.plotter import scatter_multigroup_xs
from infermc.multigroupxs import xs_types


groups = group_structures['CASMO']['2-group']

batches = range(10, 75, 5)
#batches = [10]

for batch in batches:

  filename = 'statepoint.{0}.h5'.format(batch)

  # Initialize a handle on the OpenMC statepoint file
  statepoint = openmc.statepoint.StatePoint(filename)

  # Initialize an InferMC XSTallyExtractor object to compute cross-sections
  extractor = XSTallyExtractor(statepoint)

  extractor.extractAllMultiGroupXS(groups, 'material')
  extractor.extractAllMultiGroupXS(groups, 'distribcell')

  for xs_type in xs_types:

    print batch, xs_type

    if xs_type != 'scatter matrix':
      scatter_multigroup_xs(extractor, xs_type,
                            domain_types=['distribcell', 'material'],
                            colors=['material', 'material'],
                            filename='{0}-{1}-batches'.format(xs_type,batch))

  openmc.reset_auto_ids()
