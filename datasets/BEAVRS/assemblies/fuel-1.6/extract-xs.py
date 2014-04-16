from statepoint import StatePoint
import glob
from geometry import geometry
import infermc.plotter as plotter

# Get statepoint files
files = glob.glob('statepoint.*.h5')

for file in files:
  sp = StatePoint(file)
  sp.read_results()
  plotter.plot_fluxes(geometry, sp, energies=[0, 1], gridsize=200)

  #  num_regions = geometry.getNumRegions()
#  print('The Geometry contains %d regions' % num_regions)
#  cells_to_tallies = plotter.get_cells_to_tallies(sp)

#  fluxes = np.zeros((num_regions, 2))

#  for i in range(num_regions):
#    coords = geometry.findRegion(i)
#    path = plotter.get_path(coords)

#    cell_id = path[-1]
#    tally_id = cells_to_tallies[cell_id]
#    filters = [('distribcell', path)]

#    try:
#      flux = sp.get_value(tally_id, filters, 0)
#    except:
#      print i, cell_id, tally_id, path
#      flux = [0., 0.]

#    fluxes[i,:] = flux
#    print fluxes