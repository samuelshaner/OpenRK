from statepoint import StatePoint
from geometry import geometry
from infermc.process import XSTallyExtractor
from infermc.build import EnergyGroups
import numpy as np
import matplotlib.pyplot as plt
from openmoc import Material, material_id


# Get statepoint files
files = ['statepoint.10.h5']

energy_groups = EnergyGroups()
energy_groups.setGroupEdges([0.0, 0.625e-6, 20.])

edges = energy_groups.getGroupEdges()
num_groups = energy_groups.getNumGroups()
num_regions = geometry.getNumRegions()


for file in files:

  statepoint = StatePoint(file)
  statepoint.read_results()
  statepoint.generate_stdev()

  extractor = XSTallyExtractor(statepoint=statepoint, geometry=geometry)

  # Initialize empty arrays of cross-sections for each region, group
  sigma_t = np.zeros((num_regions, num_groups))
  sigma_a = np.zeros((num_regions, num_groups))
  sigma_f = np.zeros((num_regions, num_groups))
  nusigma_f = np.zeros((num_regions, num_groups))
  sigma_s = np.zeros((num_regions, num_groups, num_groups))
  chi = np.zeros((num_regions, num_groups))

  # Initialize empty array of OpenMOC Materials
  materials = np.empty(num_regions, dtype=Material)


  for region in range(num_regions):

    sigma_t[region,:] = extractor.getXS('transport', edges, region)
    sigma_a[region,:] = extractor.getXS('absorption', edges, region)
    sigma_f[region,:] = extractor.getXS('fission', edges, region)
    nusigma_f[region,:] = extractor.getXS('nu-fission', edges, region)
    sigma_s[region,:,:] = extractor.getXS('scatter matrix', edges, region)
    chi[region,:] = extractor.getXS('chi', edges, region)

    #FIXME: Must reorder arrays in energy groups!!!
    #FIXME: sigma_s must be raveled!!
    materials[region] = Material(material_id())
    materials[region].setNumEnergyGroups(num_groups)
    materials[region].setSigmaT(sigma_t[region,:])
    materials[region].setSigmaA(sigma_a[region,:])
    materials[region].setSigmaF(sigma_f[region,:])
    materials[region].setNuSigmaF(nusigma_f[region,:])
    materials[region].setSigmaS(np.ravel(sigma_s[region,:,:]))
    materials[region].setChi(chi[region,:])


  fig = plt.figure()
  ax = fig.add_subplot(111)
  plt.scatter(sigma_t[:,0], sigma_t[:,1])
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
  plt.scatter(sigma_a[:,0], sigma_a[:,1])
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
  plt.scatter(sigma_f[:,0], sigma_f[:,1])
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
  plt.scatter(nusigma_f[:,0], nusigma_f[:,1])
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
  plt.scatter(sigma_s[:,0], sigma_s[:,1])
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