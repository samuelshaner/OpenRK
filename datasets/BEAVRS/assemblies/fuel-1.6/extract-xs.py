from statepoint import StatePoint
import glob
from geometry import geometry
import infermc.plotter as plotter
from infermc.process import TallyExtractor
import numpy as np
import matplotlib.pyplot as plt

# Get statepoint files
files = glob.glob('statepoint.*.h5')

for file in files:
  sp = StatePoint(file)
  sp.read_results()
#  plotter.plot_fluxes(geometry, sp, energies=[0, 1], gridsize=200)

  tally_extractor = TallyExtractor(statepoint=sp, geometry=geometry)

  num_regions = geometry.getNumRegions()
  flux = np.zeros((num_regions, 2))
  absorb = np.zeros((num_regions, 2))

  for region in range(num_regions):

    flux_data = tally_extractor.getDistribcellTallyData(region, 'flux')
    flux[region, :] = flux_data

    absorb_data = tally_extractor.getDistribcellTallyData(region, 'absorption')
    absorb[region, :] = absorb_data

  print flux
  print absorb

  absorb_xs = absorb / flux
  print absorb_xs

  fig = plt.figure()
  plt.scatter(absorb_xs[:,0], absorb_xs[:,1])
  plt.title('Absorption XS')
  plt.xlabel('Thermal Group $\Sigma_{abs}$ [cm$^{-1}$]')
  plt.ylabel('Fast Group $\Sigma_{abs}$ [cm$^{-1}$]')
  plt.grid()
  plt.savefig('abs-xs.png', bbox_inches='tight')


