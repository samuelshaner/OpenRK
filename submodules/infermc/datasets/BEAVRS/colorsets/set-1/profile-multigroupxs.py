import openmc
import openmc.statepoint
from datasets.energy_groups import group_structures
from infermc.process import MicroXSTallyExtractor
import infermc.plotter as plotter


groups = group_structures['CASMO']['2-group']


def profile():

  batch = 50

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


import cProfile
cProfile.run('profile()', 'stats')

import pstats
p = pstats.Stats('stats')
p.sort_stats('time').print_stats(40)
