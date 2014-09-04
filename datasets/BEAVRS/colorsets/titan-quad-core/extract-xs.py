import openmc
import openmc.statepoint
from datasets.energy_groups import group_structures
from infermc.process import XSTallyExtractor, MicroXSTallyExtractor
from infermc.multigroupxs import xs_types
import infermc.plotter as plotter




#batches = range(55, 255, 5)
min_batch = 255
max_batch = 2000
batch_interval = 1


groups = group_structures['CASMO']['2-group']

curr_batch =


for batch in batches:

  print batch

  filename = 'statepoint.{0:03d}.h5'.format(batch)

  # Initialize a handle on the OpenMC statepoint file
  statepoint = openmc.statepoint.StatePoint(filename)

  ## MICROS
  micro_extractor = MicroXSTallyExtractor(statepoint)
  micro_extractor.extractAllMultiGroupXS(groups, 'material')
  micro_extractor.extractAllMultiGroupXS(groups, 'distribcell')
  micro_extractor.checkXS()

#  plotter.scatter_micro_xs(micro_extractor,
#                           domain_types=['distribcell', 'material'],
#                           colors=['cell', 'material'],
#                           filename='{0}-batch'.format(batch))

  materials = micro_extractor._openmc_geometry.getAllMaterials()

  openmc.reset_auto_ids()
  del micro_extractor, statepoint