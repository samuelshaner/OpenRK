from datasets.energy_groups import group_structures
import openmc
from openmc.statepoint import StatePoint
from infermc.process import MicroXSTallyExtractor
import infermc
import infermc.plotter as plotter


def profile():

  #batches = range(10, 55, 5)
  batches = [30]

  groups = group_structures['CASMO']['2-group']

  for batch in batches:

    print(batch)

    filename = 'statepoint.{0:03d}.h5'.format(batch)

    # Initialize a handle on the OpenMC statepoint file
    statepoint = openmc.statepoint.StatePoint(filename)

    ## MICROS
    micro_extractor = MicroXSTallyExtractor(statepoint)
    micro_extractor.extractAllMultiGroupXS(groups, 'material')
    micro_extractor.extractAllMultiGroupXS(groups, 'distribcell')
    micro_extractor.checkXS()

    nuclides = micro_extractor._openmc_geometry.getAllNuclides()

    for xs_type in infermc.xs_types:

      if xs_type != 'scatter matrix':

        for nuclide_name, nuclide_tuple in nuclides.items():
          plotter.scatter_micro_xs(micro_extractor, xs_type, nuclide_tuple[0],
                                domain_types=['distribcell', 'material'],
                                filename='{0}-{1}-{2}-batches'.format(nuclide_name, xs_type, batch))


import cProfile
cProfile.run('profile()', 'stats')

import pstats
p = pstats.Stats('stats')
p.sort_stats('time').print_stats(40)
