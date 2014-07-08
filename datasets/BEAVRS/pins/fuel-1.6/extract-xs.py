from datasets.energy_groups import group_structures
import openmc
from openmc.statepoint import StatePoint
from infermc.process import XSTallyExtractor
from infermc.plotter import scatter_multigroup_xs
from infermc.multigroupxs import xs_types


groups = group_structures['CASMO']['2-group']

batches = [50]

for batch in batches:

  print batch

  filename = 'statepoint.{0}.h5'.format(batch)

  # Initialize a handle on the OpenMC statepoint file
  statepoint = openmc.statepoint.StatePoint(filename)

  # Initialize an InferMC XSTallyExtractor object to compute cross-sections
  extractor = XSTallyExtractor(statepoint)

  extractor.extractAllMultiGroupXS(groups, 'material')
  extractor.extractAllMultiGroupXS(groups, 'distribcell')

  multigroup_xs = extractor._multigroup_xs['material']

  for material in multigroup_xs.keys():
    for xs_type in xs_types:
      xs = multigroup_xs[material][xs_type]
      xs.printPDF(filename='material-{0}-{1}'.format(material, xs_type),
                  directory='validate', uncertainties=True)

  openmc.reset_auto_ids()
  del extractor, statepoint