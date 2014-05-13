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

  plotter.plot_fluxes(geometry, sp, energies=[0, 1], gridsize=200)

  tally_extractor = TallyExtractor(statepoint=sp, geometry=geometry)

  num_regions = geometry.getNumRegions()
  flux = np.zeros((num_regions, 2))
  total = np.zeros((num_regions, 2))
  fiss = np.zeros((num_regions, 2))
  absorb = np.zeros((num_regions, 2))
  scatter = np.zeros((num_regions, 2))

  for region in range(num_regions):

    flux_data = tally_extractor.getDistribcellTallyData(region, 'flux')
    flux[region, :] = flux_data

    total_data = tally_extractor.getDistribcellTallyData(region, 'total')
    total[region, :] = total_data

    fiss_data = tally_extractor.getDistribcellTallyData(region, 'fission')
    fiss[region, :] = fiss_data

    absorb_data = tally_extractor.getDistribcellTallyData(region, 'absorption')
    absorb[region, :] = absorb_data

    scatter_data = tally_extractor.getDistribcellTallyData(region, 'scatter')
    scatter[region, :] = scatter_data

  total_xs = total / flux
  fiss_xs = fiss / flux
  absorb_xs = absorb / flux
  scatter_xs = scatter / flux

  # Create scatter plot of the total cross-sections in each region, group
  fig = plt.figure()
  ax = fig.add_subplot(111)
  plt.scatter(total_xs[:,0], total_xs[:,1])
  plt.title('Total XS')
  plt.xlabel('Thermal Group $\Sigma_{abs}$ [cm$^{-1}$]')
  plt.ylabel('Fast Group $\Sigma_{abs}$ [cm$^{-1}$]')
  plt.grid()
  filename = 'tot-xs-' + str(sp.current_batch) + '.png'
  plt.savefig(filename, bbox_inches='tight')
  ax.text(3, 8, 'boxed italics text in data coords', style='italic',
          bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})

  # Create scatter plot of the fission cross-sections in each region, group
  fig = plt.figure()
  plt.scatter(fiss_xs[:,0], fiss_xs[:,1])
  plt.title('Fission XS')
  plt.xlabel('Thermal Group $\Sigma_{abs}$ [cm$^{-1}$]')
  plt.ylabel('Fast Group $\Sigma_{abs}$ [cm$^{-1}$]')
  plt.grid()
  filename = 'fiss-xs-' + str(sp.current_batch) + '.png'
  plt.savefig(filename, bbox_inches='tight')

  # Create scatter plot of the scatter cross-sections in each region, group
  fig = plt.figure()
  plt.scatter(scatter_xs[:,0], scatter_xs[:,1])
  plt.title('Scatter XS')
  plt.xlabel('Thermal Group $\Sigma_{abs}$ [cm$^{-1}$]')
  plt.ylabel('Fast Group $\Sigma_{abs}$ [cm$^{-1}$]')
  plt.grid()
  filename = 'scatter-xs-' + str(sp.current_batch) + '.png'
  plt.savefig(filename, bbox_inches='tight')

  # Create scatter plot of the absorption cross-sections in each region, group
  fig = plt.figure()
  plt.scatter(absorb_xs[:,0], absorb_xs[:,1])
  plt.title('Absorption XS')
  plt.xlabel('Thermal Group $\Sigma_{abs}$ [cm$^{-1}$]')
  plt.ylabel('Fast Group $\Sigma_{abs}$ [cm$^{-1}$]')
  plt.grid()
  filename = 'abs-xs-' + str(sp.current_batch) + '.png'
  plt.savefig(filename, bbox_inches='tight')