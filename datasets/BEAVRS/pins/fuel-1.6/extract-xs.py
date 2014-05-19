from statepoint import StatePoint
from geometry import geometry
from infermc.process import XSTallyExtractor
from infermc.build import EnergyGroups
import numpy as np
import matplotlib.pyplot as plt


# Get statepoint files
files = ['statepoint.10.h5', 'statepoint.15.h5']

energy_groups = EnergyGroups()
energy_groups.setGroupEdges([0.0, 0.625, 1e7])
edges = energy_groups.getGroupEdges()


for file in files:

  statepoint = StatePoint(file)
  statepoint.read_results()

  extractor = XSTallyExtractor(statepoint=statepoint, geometry=geometry)

  num_regions = geometry.getNumRegions()

  # Initialize empty arrays of
  tot_xs = np.zeros((num_regions, 2))
  abs_xs = np.zeros((num_regions, 2))
  fiss_xs = np.zeros((num_regions, 2))
  nufiss_xs = np.zeros((num_regions, 2))
  scat_xs = np.zeros((num_regions, 2))
  chi = np.zeros((num_regions, 2))


  print chi

  for region in range(num_regions):

    tot_xs[region, :] = extractor.getXS('total', edges, domain_id=region)
    abs_xs[region, :] = extractor.getXS('absorption', edges, domain_id=region)
    fiss_xs[region, :] = extractor.getXS('fission', edges, domain_id=region)
    nufiss_xs[region, :] = extractor.getXS('nu-fission', edges, domain_id=region)
    scat_xs[region, :] = extractor.getXS('scatter', edges, domain_id=region)
    chi[region, :] = extractor.getXS('chi', edges, domain_id=region)


  # Create scatter plot of the total cross-sections in each region, group
  fig = plt.figure()
  ax = fig.add_subplot(111)
  plt.scatter(tot_xs[:,0], tot_xs[:,1])
  plt.title('Total XS')
  plt.xlabel('Thermal Group $\Sigma_{tot}$ [cm$^{-1}$]')
  plt.ylabel('Fast Group $\Sigma_{tot}$ [cm$^{-1}$]')
  plt.grid()
  filename = 'tot-xs-' + str(statepoint.current_batch) + '.png'
  plt.savefig(filename, bbox_inches='tight')
  ax.text(3, 8, 'boxed italics text in data coords', style='italic',
          bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})

  fig = plt.figure()
  ax = fig.add_subplot(111)
  plt.scatter(abs_xs[:,0], abs_xs[:,1])
  plt.title('Absorption XS')
  plt.xlabel('Thermal Group $\Sigma_{abs}$ [cm$^{-1}$]')
  plt.ylabel('Fast Group $\Sigma_{abs}$ [cm$^{-1}$]')
  plt.grid()
  filename = 'abs-xs-' + str(statepoint.current_batch) + '.png'
  plt.savefig(filename, bbox_inches='tight')
  ax.text(3, 8, 'boxed italics text in data coords', style='italic',
          bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})

  fig = plt.figure()
  ax = fig.add_subplot(111)
  plt.scatter(fiss_xs[:,0], fiss_xs[:,1])
  plt.title('Fission XS')
  plt.xlabel('Thermal Group $\Sigma_{fiss}$ [cm$^{-1}$]')
  plt.ylabel('Fast Group $\Sigma_{fis}$ [cm$^{-1}$]')
  plt.grid()
  filename = 'fiss-xs-' + str(statepoint.current_batch) + '.png'
  plt.savefig(filename, bbox_inches='tight')
  ax.text(3, 8, 'boxed italics text in data coords', style='italic',
          bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})

  fig = plt.figure()
  ax = fig.add_subplot(111)
  plt.scatter(nufiss_xs[:,0], nufiss_xs[:,1])
  plt.title('NuFission XS')
  plt.xlabel('Thermal Group $\nu\Sigma_{fiss}$ [cm$^{-1}$]')
  plt.ylabel('Fast Group $\nu\Sigma_{fiss}$ [cm$^{-1}$]')
  plt.grid()
  filename = 'nufiss-xs-' + str(statepoint.current_batch) + '.png'
  plt.savefig(filename, bbox_inches='tight')
  ax.text(3, 8, 'boxed italics text in data coords', style='italic',
          bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})

  fig = plt.figure()
  ax = fig.add_subplot(111)
  plt.scatter(scat_xs[:,0], scat_xs[:,1])
  plt.title('Scatter XS')
  plt.xlabel('Thermal Group $\Sigma_{scat}$ [cm$^{-1}$]')
  plt.ylabel('Fast Group $\Sigma_{scat}$ [cm$^{-1}$]')
  plt.grid()
  filename = 'scat-xs-' + str(statepoint.current_batch) + '.png'
  plt.savefig(filename, bbox_inches='tight')
  ax.text(3, 8, 'boxed italics text in data coords', style='italic',
          bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})

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
