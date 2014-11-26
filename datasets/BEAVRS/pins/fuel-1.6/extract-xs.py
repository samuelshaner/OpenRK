import openmc
from openmc.statepoint import StatePoint
from openmc.summary import Summary
from datasets.energy_groups import group_structures
from infermc.process import XSTallyExtractor, MicroXSTallyExtractor
from infermc.multigroupxs import xs_types
import infermc.plotter as plotter


#batches = range(15, 35, 5)
batches = [30]

groups = group_structures['CASMO']['2-group']

summary = Summary('summary.h5')

for batch in batches:

  print batch

  filename = 'statepoint.{0:02}.h5'.format(batch)

  # Initialize a handle on the OpenMC statepoint file
  statepoint = StatePoint(filename)

  ## MICROS
  micro_extractor = MicroXSTallyExtractor(statepoint, summary)
  micro_extractor.extractAllMultiGroupXS(groups, 'material')
  micro_extractor.extractAllMultiGroupXS(groups, 'distribcell')
  micro_extractor.checkXS()


#  plotter.scatter_micro_xs(micro_extractor,
#                           domain_types=['distribcell', 'material'],
#                           colors=['cell', 'material'],
#                           filename='{0}-batch'.format(batch))

  materials = micro_extractor._openmc_geometry.get_all_materials()

    # DUMP-TO-FILE and PRINT XS
  for material in materials:
    for xs_type in xs_types:
      xs = micro_extractor._multigroup_xs['material'][material._id][xs_type]
      xs.exportResults()
#      xs.printPDF(directory='micro', filename='material-{0}-{1}'.format(material._id, xs_type))

  openmc.reset_auto_ids()
  del micro_extractor, statepoint
