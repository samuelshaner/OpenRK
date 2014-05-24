from statepoint import StatePoint
from geometry import geometry
from infermc.process import XSTallyExtractor
from infermc.build import EnergyGroups
import numpy as np
import matplotlib.pyplot as plt


# Get statepoint files
files = ['statepoint.10.h5', 'statepoint.15.h5']

energy_groups = EnergyGroups()
energy_groups.setGroupEdges([0.0, 0.625e-6, 20.])
edges = energy_groups.getGroupEdges()


for file in files:

  statepoint = StatePoint(file)
  statepoint.read_results()

  extractor = XSTallyExtractor(statepoint=statepoint, geometry=geometry)

  num_regions = geometry.getNumRegions()

  # Initialize empty arrays of
  total_xs = np.zeros((num_regions, 2))
  transport_xs = np.zeros((num_regions, 2))
  absorb_xs = np.zeros((num_regions, 2))
  fission_xs = np.zeros((num_regions, 2))
  nufission_xs = np.zeros((num_regions, 2))
  scatter_xs = np.zeros((num_regions, 2))
  scatter_matrix = np.zeros((num_regions, 2, 2))
  chi = np.zeros((num_regions, 2))

  for region in range(num_regions):

    total_xs[region, :] = extractor.getXS('total', edges, region)
    transport_xs[region, :] = extractor.getXS('transport', edges, region)
    absorb_xs[region, :] = extractor.getXS('absorption', edges, region)
    fission_xs[region, :] = extractor.getXS('fission', edges, region)
    nufission_xs[region, :] = extractor.getXS('nu-fission', edges, region)
    scatter_xs[region, :] = extractor.getXS('nu-scatter', edges, region)
    scatter_matrix[region, :, :] = extractor.getXS('scatter matrix', edges, region)
    chi[region, :] = extractor.getXS('chi', edges, region)


  # Create scatter plot of the total cross-sections in each region, group
  fig = plt.figure()
  ax = fig.add_subplot(111)
  plt.scatter(total_xs[:,0], total_xs[:,1])
  plt.title('Total XS')
  plt.xlabel('Thermal Group $\Sigma_{tot}$ [cm$^{-1}$]')
  plt.ylabel('Fast Group $\Sigma_{tot}$ [cm$^{-1}$]')
  plt.grid()
  filename = 'tot-xs-' + str(statepoint.current_batch) + '.png'
  plt.savefig(filename, bbox_inches='tight')
  ax.text(3, 8, 'boxed italics text in data coords', style='italic',
          bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
  plt.close(fig)

  fig = plt.figure()
  ax = fig.add_subplot(111)
  plt.scatter(transport_xs[:,0], transport_xs[:,1])
  plt.title('Transport XS')
  plt.xlabel('Thermal Group $\Sigma_{tr}$ [cm$^{-1}$]')
  plt.ylabel('Fast Group $\Sigma_{tr}$ [cm$^{-1}$]')
  plt.grid()
  filename = 'trans-xs-' + str(statepoint.current_batch) + '.png'
  plt.savefig(filename, bbox_inches='tight')
  ax.text(3, 8, 'boxed italics text in data coords', style='italic',
          bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
  plt.close(fig)

  fig = plt.figure()
  ax = fig.add_subplot(111)
  plt.scatter(absorb_xs[:,0], absorb_xs[:,1])
  plt.title('Absorption XS')
  plt.xlabel('Thermal Group $\Sigma_{abs}$ [cm$^{-1}$]')
  plt.ylabel('Fast Group $\Sigma_{abs}$ [cm$^{-1}$]')
  plt.grid()
  filename = 'abs-xs-' + str(statepoint.current_batch) + '.png'
  plt.savefig(filename, bbox_inches='tight')
  ax.text(3, 8, 'boxed italics text in data coords', style='italic',
          bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
  plt.close(fig)

  fig = plt.figure()
  ax = fig.add_subplot(111)
  plt.scatter(fission_xs[:,0], fission_xs[:,1])
  plt.title('Fission XS')
  plt.xlabel('Thermal Group $\Sigma_{fiss}$ [cm$^{-1}$]')
  plt.ylabel('Fast Group $\Sigma_{fis}$ [cm$^{-1}$]')
  plt.grid()
  filename = 'fiss-xs-' + str(statepoint.current_batch) + '.png'
  plt.savefig(filename, bbox_inches='tight')
  ax.text(3, 8, 'boxed italics text in data coords', style='italic',
          bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
  plt.close(fig)

  fig = plt.figure()
  ax = fig.add_subplot(111)
  plt.scatter(nufission_xs[:,0], nufission_xs[:,1])
  plt.title('NuFission XS')
  plt.xlabel('Thermal Group $\nu\Sigma_{fiss}$ [cm$^{-1}$]')
  plt.ylabel('Fast Group $\nu\Sigma_{fiss}$ [cm$^{-1}$]')
  plt.grid()
  filename = 'nufiss-xs-' + str(statepoint.current_batch) + '.png'
  plt.savefig(filename, bbox_inches='tight')
  ax.text(3, 8, 'boxed italics text in data coords', style='italic',
          bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
  plt.close(fig)

  fig = plt.figure()
  ax = fig.add_subplot(111)
  plt.scatter(scatter_xs[:,0], scatter_xs[:,1])
  plt.title('Scatter XS')
  plt.xlabel('Thermal Group $\Sigma_{scat}$ [cm$^{-1}$]')
  plt.ylabel('Fast Group $\Sigma_{scat}$ [cm$^{-1}$]')
  plt.grid()
  filename = 'scat-xs-' + str(statepoint.current_batch) + '.png'
  plt.savefig(filename, bbox_inches='tight')
  ax.text(3, 8, 'boxed italics text in data coords', style='italic',
          bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
  plt.close(fig)

  fig = plt.figure()
  ax = fig.add_subplot(111)
  plt.scatter(chi[:,0], chi[:,1])
  plt.title('Chi')
  plt.xlabel('Thermal Group $\Xi$ [cm$^{-1}$]')
  plt.ylabel('Fast Group $\Xi$ [cm$^{-1}$]')
  plt.grid()
  filename = 'chi-' + str(statepoint.current_batch) + '.png'
  plt.savefig(filename, bbox_inches='tight')
  ax.text(3, 8, 'boxed italics text in data coords', style='italic',
          bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
  plt.close(fig)