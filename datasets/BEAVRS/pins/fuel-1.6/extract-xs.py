from statepoint import StatePoint
from geometry import geometry
from infermc.process import XSTallyExtractor
from infermc.build import EnergyGroups
import numpy as np
import matplotlib.pyplot as plt
from openmoc import Material, material_id
from openmoc.materialize import export_materials


# Get statepoint files
files = ['statepoint.10.h5']

energy_groups = EnergyGroups()
energy_groups.setGroupEdges([0.0, 0.625e-6, 20.])

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
  sigma_s = np.zeros((num_regions, num_groups**2))
  chi = np.zeros((num_regions, num_groups))

  # Initialize empty dictionary of OpenMOC Materials
  materials = dict()


  for region in range(num_regions):

    sigma_t[region,:] = extractor.getXS('transport', energy_groups, region)
    sigma_a[region,:] = extractor.getXS('absorption', energy_groups, region)
    sigma_f[region,:] = extractor.getXS('fission', energy_groups, region)
    nusigma_f[region,:] = extractor.getXS('nu-fission', energy_groups, region)
    sigma_s[region,:] = extractor.getXS('scatter matrix', energy_groups, region)
    chi[region,:] = extractor.getXS('chi', energy_groups, region)

    #FIXME: Must reorder arrays in energy groups!!!
    #FIXME: sigma_s must be raveled!!
    name = 'FSR-%d' % region
    materials[name] = Material(material_id())
    materials[name].setNumEnergyGroups(num_groups)
    materials[name].setSigmaT(sigma_t[region,:])
    materials[name].setSigmaA(sigma_a[region,:])
    materials[name].setSigmaF(sigma_f[region,:])
    materials[name].setNuSigmaF(nusigma_f[region,:])
    materials[name].setSigmaS(sigma_s[region,:])
    materials[name].setChi(chi[region,:])

  print sigma_s[0,:]

  export_materials(materials, use_hdf5=True)


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