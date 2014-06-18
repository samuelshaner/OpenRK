from openmc.statepoint import StatePoint
from openmc import TalliesFile, Filter
from openmc.opencsg_compatible import get_opencsg_geometry
from infermc.process import get_path

# Retrieve statepoint files
files = ['statepoint.10.h5']

sp = StatePoint('statepoint.10.h5')

sp.read_results()
sp.read_source()
sp.compute_ci()

tallies = sp._tallies

for tally_id in tallies:
  tally = sp._tallies[tally_id]
  tally = sp.get_tally(tally._scores[0], tally._filters,
                       tally._label, tally._estimator)
#  print tally

#TODO: Check to see if we can get a value for a tally given it's path!!!
#TODO: Check to see if we can get a value for a tally (cell, universe or material)

openmc_geometry = sp._geometry
opencsg_geometry = get_opencsg_geometry(openmc_geometry)
num_regions = opencsg_geometry._num_regions


for region in range(num_regions):

  # Get the path to the region in the geometry
  coord = opencsg_geometry.findRegion(region)
  path = get_path(coord)

  filters = list()

  cell_id = path[-1]
  path_filter = Filter(type='distribcell', bin_edges=cell_id)
  filters.append(path_filter)

  energy_filter = Filter(type='energy', bin_edges=[0.0, 0.625e-6, 10.])
  filters.append(energy_filter)

  tally = sp.get_tally('flux', filters, label='2 groups')

  true_path_filter = tally.findFilter('distribcell', [cell_id])
  offset = openmc_geometry.getOffset(path, true_path_filter._offset)

  test_value = tally.getValue('flux', filters, [offset, 0])
  print test_value


tallies_file = TalliesFile()

for tally_id in tallies:
  tallies_file.addTally(tallies[tally_id])

# tallies_file.exportToXML()