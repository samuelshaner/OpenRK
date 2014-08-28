from infermc.checktype import *
from infermc.uncorr_math import *
import openmc
import numpy as np
import os

# Type-checking support
from typecheck import accepts, Or, Exact, Self


xs_types = ['total',
            'transport',
            'absorption',
            'scatter',
            'nu-scatter',
            'scatter matrix',
            'fission',
            'nu-fission',
            'chi']
#            'diffusion']

domain_types = ['cell',
                'distribcell',
                'universe',
                'material',
                'mesh']

domains = [openmc.Cell, openmc.Universe, openmc.Material, openmc.Mesh]

# For use with typecheck module
xs_types_check = Or(*[Exact(xs_type) for xs_type in xs_types])
domain_types_check = Or(*[Exact(domain_type) for domain_type in domain_types])
domains_check = Or(*[domain for domain in domains])


class EnergyGroups(object):

  def __init__(self):

    self._group_edges = None
    self._num_groups = None


  @property
  def group_edges(self):
    return self._group_edges


  @group_edges.setter
  @accepts(Self(), Or(list, tuple, np.ndarray))
  def group_edges(self, edges):
    self._group_edges = np.array(edges)
    self._num_groups = len(edges)-1


  @accepts(Self(), Or(int, float), Or(int, float), int,
           Or(Exact('linear'), Exact('logarithmic')))
  def generateBinEdges(self, start, stop, num_groups, type='linear'):

    self._num_groups = num_groups

    if type == 'linear':
      self._group_edges = np.linspace(start, stop, num_groups+1)
    else:
      self._group_edges = np.logspace(np.log(start), np.log(stop), num_groups+1)


  @accepts(Self(), Or(int, float))
  def getGroup(self, energy):

    # Assumes energy is in eV

    if self._group_edges is None:
      msg = 'Unable to get energy group for energy {0} eV since ' \
            'the group edges have not yet been set'.format(energy)
      raise ValueError(msg)

    index = np.where(self._group_edges > energy)[0]
    group = self._num_groups - index
    return group


  @accepts(Self(), int)
  def getGroupBounds(self, group):

    if self._group_edges is None:
      msg = 'Unable to get energy group bounds for group {0} since ' \
            'the group edges have not yet been set'.format(group)
      raise ValueError(msg)

    lower = self._group_edges[self._num_groups-group]
    upper = self._group_edges[self._num_groups-group+1]
    return (lower, upper)



class MultiGroupXS(object):

  def __init__(self, domain=None, domain_type=None, energy_groups=None):

    self._xs_type = None
    self._domain = None
    self._domain_type = None
    self._energy_groups = None
    self._num_groups = None
    self._tallies = dict()
    self._xs = None

    self._offset = None

    # A dictionary used to compute indices into the xs array
    # Keys   - Domain ID (ie, Material ID, Region ID for districell, etc)
    # Values - Offset/stride into xs array
    self._subdomain_offsets = dict()

    if not domain_type is None:
      self.domain_type = domain_type

    if not domain is None:
      self.domain = domain

    if not energy_groups is None:
      self.energy_groups = energy_groups


  @property
  def energy_groups(self):
    return self._energy_groups


  @energy_groups.setter
  @accepts(Self(), EnergyGroups)
  def energy_groups(self, energy_groups):
    self._energy_groups = energy_groups
    self._num_groups = energy_groups._num_groups


  @property
  def num_groups(self):
    return self._num_groups


  @property
  def domain(self):
    return self._domain


  @domain.setter
  @accepts(Self(), domains_check)
  def domain(self, domain):

    self._domain = domain

    if self._domain_type in ['material', 'cell', 'universe']:
      self._subdomain_offsets[domain._id] = 0


  @property
  def domain_type(self):
    return self._domain_type


  @domain_type.setter
  @accepts(Self(), domain_types_check)
  def domain_type(self, domain_type):
    self._domain_type = domain_type


  def findDomainOffset(self):

    if self._domain_type in ['material', 'cell', 'universe']:
      self._offset = 0

    # Distribcell tallies
    else:
      first_tally = self._tallies[self._tallies.keys()[0]]
      domain_filter = first_tally.findFilter('distribcell', [self._domain._id])
      self._offset = domain_filter._offset


  @accepts(Self(), int, int)
  def setSubDomainOffset(self, domain_id, offset):
    self._subdomain_offsets[domain_id] = offset


  def createTallies(self):

    if self._energy_groups is None:
      msg = 'Unable to create Tallies since the energy groups has not been set'
      raise ValueError(msg)

    elif self._domain is None:
      msg = 'Unable to create Tallies since the domain has not been set'
      raise ValueError(msg)

    elif self._domain_type is None:
      msg = 'Unable to create Tallies since the domain type has not been set'
      raise ValueError(msg)

    return


  @accepts(Self(), int, Or(Exact(None), int), str)
  def getXS(self, group, subdomain=None, metric='mean'):

    if self._xs is None:
      msg = 'Unable to get cross-section since it has not been computed'
      raise ValueError(msg)

    # Compute the energy group index into the array
    # Subtract one from group number
    group_index = group - 1

    # Compute the subdomain index into the array
    if self._domain_type == 'distribcell':
      subdomain_index = self._subdomain_offsets[subdomain]
    else:
      subdomain_index = 0

    # Batch mean
    if metric == 'mean':
      return self._xs[0, subdomain_index, group_index]

    # Batch standard deviation
    else:
      return self._xs[1, subdomain_index, group_index]


  @accepts(Self(), str, str)
  def dumpToFile(self, filename='multigroupxs', directory='multigroupxs'):

    import pickle

    # Make directory if it does not exist
    if not os.path.exists(directory):
      os.makedirs(directory)

    # Create an empty dictionary to store the data
    xs_results = dict()

    # Store all of this MultiGroupXS' class attributes in the dictionary
    xs_results['xs type'] = self._xs_type
    xs_results['domain type'] = self._domain_type
    xs_results['domain'] = self._domain
    xs_results['energy groups'] = self._energy_groups
    xs_results['tallies'] = self._tallies
    xs_results['xs'] = self._xs
    xs_results['offset'] = self._offset
    xs_results['subdomain offsets'] = self._subdomain_offsets

    # Pickle the MultiGroupXS results to a file
    filename = directory + '/' + filename + '.pkl'
    filename = filename.replace(' ', '-')
    pickle.dump(xs_results, open(filename, 'wb'))


  @accepts(Self(), str, str)
  def restoreFromFile(self, filename, directory='.'):

    import pickle
    import os.path

    filename = directory + '/' + filename + '.pkl'
    filename = filename.replace(' ', '-')

    # Check that the file exists
    if not os.path.exists(filename):
      msg = 'Unable to import cross-section from filename={0} ' \
            'which does not exist'.format(filename)
      raise ValueError(msg)

    # Load the pickle file into a dictionary
    xs_results = pickle.load(open(filename, 'rb'))

    # Extract the MultiGroupXS' class attributes in the dictionary
    xs_type = xs_results['xs type']
    domain_type = xs_results['domain type']
    domain = xs_results['domain']
    energy_groups = xs_results['energy groups']
    tallies = xs_results['tallies']
    xs = xs_results['xs']
    offset = xs_results['offset']
    subdomain_offsets = xs_results['subdomain offsets']

    # Store the MultiGroupXS class attributes
    self._xs_type = xs_type
    self.domain_type = domain_type
    self.domain = domain
    self.energy_groups = energy_groups
    self._tallies = tallies
    self._xs = xs
    self._offset = offset
    self._subdomain_offsets = subdomain_offsets



  #FIXME: Modularize!!!
  @accepts(Self(), Or(Exact(None), int), str, str, str, bool, bool)
  def exportSubdomainResults(self, subdomain=None, filename='multigroupxs',
                             directory='multigroupxs', format='hdf5',
                             append=True, uncertainties=False):

    if self._xs is None:
      msg = 'Unable to export cross-section {0} since it has not yet ' \
            'been computed'.format(self._xs_type)
      raise ValueError(msg)

    # Make directory if it does not exist
    if not os.path.exists(directory):
      os.makedirs(directory)

    # HDF5 binary file
    if format == 'hdf5':

      import h5py

      filename = directory + '/' + filename + '.h5'
      filename = filename.replace(' ', '-')

      if append:
        xs_results = h5py.File(filename, 'a')
      else:
        xs_results = h5py.File(filename, 'w')

      # Create an HDF5 group within the file for this particular MultiGroupXS
      domain_type_group = xs_results.require_group(self._domain_type)
      domain_group = domain_type_group.require_group('Cell {0}'.format(self._domain._id))
      xs_group = domain_group.require_group(self._xs_type)

      # Compute the subdomain index into the array
      if self._domain_type == 'distribcell':
        subdomain_index = self._subdomain_offsets[subdomain]
      else:
        subdomain_index = 0

      # Add MultiGroupXS results data to the HDF5 group
      average = self._xs[0, ...]
      shape = average[subdomain_index,:].shape
      xs_group.require_dataset('average', dtype=np.float64, shape=shape,
                               data=average[subdomain_index,:])

      if uncertainties:
        std_dev = self._xs[1, ...]
        xs_group.require_dataset('std. dev.', dtype=np.float64, shape=shape,
                                 data=std_dev[subdomain_index,:])

      # Close the MultiGroup results HDF5 file
      xs_results.close()


    # Python pickle binary file
    elif format == 'pkl':

      import pickle

      # Load the dictionary from the Pickle file
      filename = directory + '/' + filename + '.pkl'
      filename = filename.replace(' ', '-')

      if os.path.exists(filename) and append:
        xs_results = pickle.load(open(filename, 'rb'))
      else:
        xs_results = dict()

      # Create a nested dictionary within the file for this particular MultiGroupXS
      if not xs_results.has_key(self._domain_type):
        domain_type_group = xs_results[self._domain_type] = dict()
      else:
        domain_type_group = xs_results[self._domain_type]

      if not domain_type_group.has_key('Cell {0}'.format(self._domain._id)):
        domain_group = domain_type_group['Cell {0}'.format(self._domain._id)] = dict()
      else:
        domain_group = domain_type_group['Cell {0}'.format(self._domain._id)]

      xs_group = domain_group[self._xs_type] = dict()

      # Compute the subdomain index into the array
      if self._domain_type == 'distribcell':
        subdomain_index = self._subdomain_offsets[subdomain]
      else:
        subdomain_index = 0

      # Add MultiGroupXS results data to the dictionary
      average = self._xs[0, ...]
      xs_group['average'] = average[subdomain_index,:]

      if uncertainties:
        std_dev = self._xs[1, ...]
        xs_group['std. dev.'] = std_dev[subdomain_index,:]

      # Pickle the MultiGroupXS results to a file
      pickle.dump(xs_results, open(filename, 'wb'))


    # LaTeX script
    elif format == 'latex':

      import tabulate

      # Load the LaTeX file
      filename = directory + '/' + filename + '.tex'
      filename = filename.replace(' ', '-')

      if os.path.exists(filename) and append:

        xs_results = open(filename, 'r')
        lines = xs_results.readlines()
        xs_results.close()
        xs_results = open(filename, 'w')

        for line in lines:
          if line != '\\end{document}\n':
            xs_results.write(line)

      else:
        xs_results = open(filename, 'w')
        xs_results.write('\\documentclass[preview, 12pt, border=1mm]{standalone}\n')
        xs_results.write('\\usepackage{caption}\n')
        xs_results.write('\\begin{document}\n')

      # Compute the subdomain index into the array
      if self._domain_type == 'distribcell':
        subdomain_index = self._subdomain_offsets[subdomain]
      else:
        subdomain_index = 0

      # Add MultiGroupXS results data to the table list
      table = list()
      average = self._xs[0, ...]

      if uncertainties:
        std_dev = self._xs[1, ...]

      if self._xs_type != 'scatter matrix':

        headers = list()
        headers.append('Group')
        headers.append('Average XS')

        if uncertainties:
          headers.append('Std. Dev.')

        for group in range(self._num_groups):
          subtable = list()
          subtable.append(group+1)
          subtable.append(average[subdomain_index,group])

          if uncertainties:
            subtable.append(std_dev[subdomain_index,group])

          table.append(subtable)

      # Scattering matrix
      else:

        headers = list()
        headers.append('Group In')
        headers.append('Group Out')
        headers.append('Average XS')

        if uncertainties:
          headers.append('Std. Dev.')

        for in_group in range(self._num_groups):
          for out_group in range(self._num_groups):

            subtable = list()
            subtable.append(in_group+1)
            subtable.append(out_group+1)
            subtable.append(average[subdomain_index, in_group, out_group])

            if uncertainties:
              subtable.append(std_dev[subdomain_index, in_group, out_group])

            table.append(subtable)

      greek = dict()
      greek['total'] = '$\\Sigma_{t}$'
      greek['transport'] = '$\\Sigma_{tr}$'
      greek['absorption'] = '$\\Sigma_{a}$'
      greek['scatter'] = '$\\Sigma_{s}$'
      greek['nu-scatter'] = '$\\nu\\Sigma_{s}$'
      greek['scatter matrix'] = '$\\Sigma_{s}$'
      greek['chi'] = '$\\chi$'
      greek['fission'] = '$\\Sigma_{f}$'
      greek['nu-fission'] = '$\\nu\\Sigma_{f}$'

      if subdomain is None:
        caption = '\\caption{{{0} {1} {2}}}'.format(
          self._domain_type.capitalize(), self._domain._id, greek[self._xs_type])
      else:
        caption = '\\caption{{{0} {1} (Region {2}) {3}}}'.format(
          self._domain_type.capitalize(), self._domain._id,
          subdomain, greek[self._xs_type])

      # Write the MultiGroupXS results to a file
      xs_results.write('\\begin{table}\n')
      xs_results.write('\\begin{center}\n')
      xs_results.write(tabulate.tabulate(table, headers, tablefmt='latex'))
      xs_results.write(caption)
      xs_results.write('\\end{center}\n')
      xs_results.write('\\end{table}\n')
      xs_results.write('\\end{document}\n')

      xs_results.close()


  @accepts(Self(), Or(Exact(None), int), str, str, bool)
  def printPDF(self, subdomain=None, filename='multigroupxs', directory='.',
               uncertainties=False):

    if self._xs is None:
      msg = 'Unable to print PDF for cross-section {0} since it has not yet ' \
            'been computed'.format(self._xs_type)
      raise ValueError(msg)

    # Make directory if it does not exist
    if not os.path.exists(directory):
      os.makedirs(directory)

    filename = filename.replace(' ', '-')

    # Generate LaTeX file
    self.exportSubdomainResults(subdomain, filename, '.',
                                'latex', False, uncertainties)

    # Compile LaTeX to PDF
    import subprocess
    proc = subprocess.Popen(['pdflatex', '{0}.tex'.format(filename)])
    proc.communicate()

    # Move PDF to requested directory and cleanup temporary LaTeX files
    if directory != '.':
      os.system('mv {0}.pdf {1}'.format(filename, directory))

    os.system('rm {0}.tex {0}.aux {0}.log'.format(filename))


  @accepts(Self(), Or(Exact(None), int))
  def printXS(self, subdomain=None):

    if self._xs is None:
      msg = 'Unable to print cross-section {0} since it has not yet ' \
            'been computed'.format(self._xs_type)
      raise ValueError(msg)

    string = 'Multi-Group XS\n'
    string += '{0: <16}{1}{2}\n'.format('\tType', '=\t', self._xs_type)
    string += '{0: <16}{1}{2}\n'.format('\tDomain Type', '=\t', self._domain_type)
    string += '{0: <16}{1}{2}\n'.format('\tDomain ID', '=\t', self._domain._id)

    # Compute the subdomain index into the array
    if self._domain_type == 'distribcell':
      subdomain_index = self._subdomain_offsets[subdomain]
    else:
      subdomain_index = 0

    string += '{0: <16}\n'.format('\tCross-Sections [cm^-1]:')

    # Loop over energy groups ranges
    for group in range(1,self._num_groups+1):
      bounds = self._energy_groups.getGroupBounds(group)
      string += '{0: <12}Group {1} [{2: <10} - ' \
                '{3: <10}MeV]:\t'.format('', group, bounds[0], bounds[1])
      average = self._xs[0, subdomain_index, group-1]
      std_dev = self._xs[1, subdomain_index, group-1]
      string += '{:.2e}+/-{:.2e}'.format(average, std_dev)
      string += '\n'

    print(string)


class TotalXS(MultiGroupXS):

  def __init__(self, domain=None, domain_type=None, energy_groups=None):
    super(TotalXS, self).__init__(domain, domain_type, energy_groups)
    self._xs_type = 'total'

  def createTallies(self):

    super(TotalXS, self).createTallies()

    # Create Filter objects
    group_edges = self._energy_groups._group_edges
    energy_filter = openmc.Filter('energy', group_edges)
    domain_filter = openmc.Filter(self._domain_type, self._domain._id)

    # Build a Tally label string
    label = '{0} groups'.format(self._num_groups)

    flux = openmc.Tally()
    flux.addScore('flux')
    flux.addFilter(domain_filter)
    flux.addFilter(energy_filter)
    flux.setLabel(label)
    self._tallies['flux'] = flux

    total = openmc.Tally()
    total.addScore('total')
    total.addFilter(domain_filter)
    total.addFilter(energy_filter)
    total.setLabel(label)
    self._tallies['total'] = total


  def computeXS(self):

    # TODO: Error checking...

    # Find the flux, total Tallies
    flux_tally = self._tallies['flux']
    total_tally = self._tallies['total']

    flux_mean = flux_tally._mean
    total_mean = total_tally._mean
    flux_std_dev = flux_tally._std_dev
    total_std_dev = total_tally._std_dev

    # Get the number of domains
    num_subdomains = flux_tally._mean.shape[0] / self._num_groups
    new_shape = (num_subdomains, self._num_groups)

#    flux_mean = np.reshape(flux_mean, new_shape)
#    total_mean = np.reshape(total_mean, new_shape)
#    flux_std_dev = np.reshape(flux_std_dev, new_shape)
#    total_std_dev = np.reshape(total_std_dev, new_shape)

    flux = np.array([flux_mean, flux_std_dev])
    total = np.array([total_mean, total_std_dev])

    # Set any zero fluxes and reaction rates to a negative value
    flux_indices = flux[0, ...] == 0.
    total_indices = total[0, ...] == 0
    flux[:, flux_indices] = -1.
    total[:, total_indices] = -1.

    # Compute the xs with uncertainty propagation
    self._xs = uncorr_divide_by_array(total, flux)

    # For any region without flux or reaction rate, convert xs to zero
    self._xs[:, flux_indices] = 0.
    self._xs[:, total_indices] = 0.

    # Correct -0.0 to +0.0
    self._xs += 0.

    # Reverse array so that it is ordered intuitively from high to low energy
    self._xs = self._xs[..., ::-1]



class TransportXS(MultiGroupXS):

  def __init__(self, domain=None, domain_type=None, energy_groups=None):
    super(TransportXS, self).__init__(domain, domain_type, energy_groups)
    self._xs_type = 'transport'

  def createTallies(self):

    super(TransportXS, self).createTallies()

    # Create Filter objects
    group_edges = self._energy_groups._group_edges
    energy_filter = openmc.Filter('energy', group_edges)
    domain_filter = openmc.Filter(self._domain_type, self._domain._id)

    # Build a Tally label string
    label = '{0} groups'.format(self._num_groups)

    flux = openmc.Tally()
    flux.addScore('flux')
    flux.addFilter(domain_filter)
    flux.addFilter(energy_filter)
    flux.setLabel(label)
    flux.setEstimator('analog')
    self._tallies['flux'] = flux

    total = openmc.Tally()
    total.addScore('total')
    total.addFilter(domain_filter)
    total.addFilter(energy_filter)
    total.setLabel(label)
    total.setEstimator('analog')
    self._tallies['total'] = total

    scatter1 = openmc.Tally()
    scatter1.addScore('scatter-1')
    scatter1.addFilter(domain_filter)
    scatter1.addFilter(energy_filter)
    scatter1.setLabel(label)
    scatter1.setEstimator('analog')
    self._tallies['scatter-1'] = scatter1


  def computeXS(self):

    # TODO: Error checking...

    # Find the flux, total, and scatter-1 Tallies
    flux_tally = self._tallies['flux']
    total_tally = self._tallies['total']
    scatter1_tally = self._tallies['scatter-1']

    flux_mean = flux_tally._mean
    total_mean = total_tally._mean
    scatter1_mean = scatter1_tally._mean
    flux_std_dev = flux_tally._std_dev
    total_std_dev = total_tally._std_dev
    scatter1_std_dev = scatter1_tally._std_dev

    # Get the number of domains and create a new array shape
    num_subdomains = flux_tally._mean.shape[0] / self._num_groups
    new_shape = (num_subdomains, self._num_groups)

#    flux_mean = np.reshape(flux_mean, new_shape)
#    total_mean = np.reshape(total_mean, new_shape)
#    scatter1_mean = np.reshape(scatter1_mean, new_shape)
#    flux_std_dev = np.reshape(flux_std_dev, new_shape)
#    total_std_dev = np.reshape(total_std_dev, new_shape)
#    scatter1_std_dev = np.reshape(scatter1_std_dev, new_shape)

    flux = np.array([flux_mean, flux_std_dev])
    total = np.array([total_mean, total_std_dev])
    scatter1 = np.array([scatter1_mean, scatter1_std_dev])

    # Set any zero fluxes to a negative value
    flux_indices = flux[0, ...] == 0.
    flux[:, flux_indices] = -1.

    delta = uncorr_sub(total, scatter1)

    # Set any zero reaction rates to a negative value
    delta_indices = delta[0, ...] == 0.
    delta[:, delta_indices] = -1.

    # Compute the xs with uncertainty propagation
    self._xs = uncorr_divide_by_array(delta, flux)

    # For any region without flux or reaction rate, convert xs to zero
    self._xs[:, flux_indices] = 0.
    self._xs[:, delta_indices] = 0.

    # Correct -0.0 to +0.0
    self._xs += 0.

    # Reverse array so that it is ordered intuitively from high to low energy
    self._xs = self._xs[..., ::-1]



class AbsorptionXS(MultiGroupXS):

  def __init__(self, domain=None, domain_type=None, energy_groups=None):
    super(AbsorptionXS, self).__init__(domain, domain_type, energy_groups)
    self._xs_type = 'absorption'

  def createTallies(self):

    super(AbsorptionXS, self).createTallies()

    # Create Filter objects
    group_edges = self._energy_groups._group_edges
    energy_filter = openmc.Filter('energy', group_edges)
    domain_filter = openmc.Filter(self._domain_type, self._domain._id)

    # Build a Tally label string
    label = '{0} groups'.format(self._num_groups)

    flux = openmc.Tally()
    flux.addScore('flux')
    flux.addFilter(domain_filter)
    flux.addFilter(energy_filter)
    flux.setLabel(label)
    self._tallies['flux'] = flux

    absorption = openmc.Tally()
    absorption.addScore('absorption')
    absorption.addFilter(domain_filter)
    absorption.addFilter(energy_filter)
    absorption.setLabel(label)
    self._tallies['absorption'] = absorption


  def computeXS(self):

    # TODO: Error checking...

    # Find the flux, total Tallies
    flux_tally = self._tallies['flux']
    absorption_tally = self._tallies['absorption']

    flux_mean = flux_tally._mean
    absorption_mean = absorption_tally._mean
    flux_std_dev = flux_tally._std_dev
    absorption_std_dev = absorption_tally._std_dev

    # Get the number of domains
    num_subdomains = flux_tally._mean.shape[0] / self._num_groups
    new_shape = (num_subdomains, self._num_groups)

#    flux_mean = np.reshape(flux_mean, new_shape)
#    absorption_mean = np.reshape(absorption_mean, new_shape)
#    flux_std_dev = np.reshape(flux_std_dev, new_shape)
#    absorption_std_dev = np.reshape(absorption_std_dev, new_shape)

    flux = np.array([flux_mean, flux_std_dev])
    absorption = np.array([absorption_mean, absorption_std_dev])

    # Set any zero fluxes and reaction rates to a negative value
    flux_indices = flux[0, ...] == 0.
    absorption_indices = absorption[0, ...] == 0.
    flux[:, flux_indices] = -1.
    absorption[:, absorption_indices] = -1.

    # Compute the xs with uncertainty propagation
    self._xs = uncorr_divide_by_array(absorption, flux)

    # For any region without flux, convert xs to zero
    self._xs[:, flux_indices] = 0.
    self._xs[:, absorption_indices] = 0.

    # Correct -0.0 to +0.0
    self._xs += 0.

    # Reverse array so that it is ordered intuitively from high to low energy
    self._xs = self._xs[..., ::-1]



class FissionXS(MultiGroupXS):

  def __init__(self, domain=None, domain_type=None, energy_groups=None):
    super(FissionXS, self).__init__(domain, domain_type, energy_groups)
    self._xs_type = 'fission'

  def createTallies(self):

    super(FissionXS, self).createTallies()

    # Create Filter objects
    group_edges = self._energy_groups._group_edges
    energy_filter = openmc.Filter('energy', group_edges)
    domain_filter = openmc.Filter(self._domain_type, self._domain._id)

    # Build a Tally label string
    label = '{0} groups'.format(self._num_groups)

    flux = openmc.Tally()
    flux.addScore('flux')
    flux.addFilter(domain_filter)
    flux.addFilter(energy_filter)
    flux.setLabel(label)
    self._tallies['flux'] = flux

    fission = openmc.Tally()
    fission.addScore('fission')
    fission.addFilter(domain_filter)
    fission.addFilter(energy_filter)
    fission.setLabel(label)
    self._tallies['fission'] = fission


  def computeXS(self):

    # TODO: Error checking...

    # Find the flux, fission Tallies
    flux_tally = self._tallies['flux']
    fission_tally = self._tallies['fission']

    flux_mean = flux_tally._mean
    fission_mean = fission_tally._mean
    flux_std_dev = flux_tally._std_dev
    fission_std_dev = fission_tally._std_dev

    # Get the number of domains
    num_subdomains = flux_tally._mean.shape[0] / self._num_groups
    new_shape = (num_subdomains, self._num_groups)

 #   flux_mean = np.reshape(flux_mean, new_shape)
 #   fission_mean = np.reshape(fission_mean, new_shape)
 #   flux_std_dev = np.reshape(flux_std_dev, new_shape)
 #   fission_std_dev = np.reshape(fission_std_dev, new_shape)

    flux = np.array([flux_mean, flux_std_dev])
    fission = np.array([fission_mean, fission_std_dev])

    # Set any zero fluxes and reaction rates to a negative value
    flux_indices = flux[0, ...] == 0.
    fission_indices = fission[0, ...] == 0.
    flux[:, flux_indices] = -1.
    fission[:, fission_indices] = -1.

    # Compute the xs with uncertainty propagation
    self._xs = uncorr_divide_by_array(fission, flux)

    # For any region without flux or reaction rate, convert xs to zero
    self._xs[:, flux_indices] = 0.
    self._xs[:, fission_indices] = 0.

    # Correct -0.0 to +0.0
    self._xs += 0.

    # Reverse array so that it is ordered intuitively from high to low energy
    self._xs = self._xs[..., ::-1]



class NuFissionXS(MultiGroupXS):

  def __init__(self, domain=None, domain_type=None, energy_groups=None):
    super(NuFissionXS, self).__init__(domain, domain_type, energy_groups)
    self._xs_type = 'nu-fission'

  def createTallies(self):

    super(NuFissionXS, self).createTallies()

    # Create Filter objects
    group_edges = self._energy_groups._group_edges
    energy_filter = openmc.Filter('energy', group_edges)
    domain_filter = openmc.Filter(self._domain_type, self._domain._id)

    # Build a Tally label string
    label = '{0} groups'.format(self._num_groups)

    flux = openmc.Tally()
    flux.addScore('flux')
    flux.addFilter(domain_filter)
    flux.addFilter(energy_filter)
    flux.setLabel(label)
    self._tallies['flux'] = flux

    nu_fission = openmc.Tally()
    nu_fission.addScore('nu-fission')
    nu_fission.addFilter(domain_filter)
    nu_fission.addFilter(energy_filter)
    nu_fission.setLabel(label)
    self._tallies['nu-fission'] = nu_fission


  def computeXS(self):

    # TODO: Error checking...

    # Find the flux, nu-fission Tallies
    flux_tally = self._tallies['flux']
    nu_fission_tally = self._tallies['nu-fission']

    flux_mean = flux_tally._mean
    nu_fission_mean = nu_fission_tally._mean
    flux_std_dev = flux_tally._std_dev
    nu_fission_std_dev = nu_fission_tally._std_dev

    # Get the number of domains
    num_subdomains = flux_tally._mean.shape[0] / self._num_groups
    new_shape = (num_subdomains, self._num_groups)

#    flux_mean = np.reshape(flux_mean, new_shape)
#    nu_fission_mean = np.reshape(nu_fission_mean, new_shape)
#    flux_std_dev = np.reshape(flux_std_dev, new_shape)
#    nu_fission_std_dev = np.reshape(nu_fission_std_dev, new_shape)

    flux = np.array([flux_mean, flux_std_dev])
    nu_fission = np.array([nu_fission_mean, nu_fission_std_dev])

    # Set any zero fluxes and reaction rates to a negative value
    flux_indices = flux[0, ...] == 0.
    nu_fission_indices = nu_fission[0, ...] == 0.
    flux[:, flux_indices] = -1.
    nu_fission[:, nu_fission_indices] = -1.

    # Compute the xs with uncertainty propagation
    self._xs = uncorr_divide_by_array(nu_fission, flux)

    # For any region without flux or reaction rate, convert xs to zero
    self._xs[:, flux_indices] = 0.
    self._xs[:, nu_fission_indices] = 0.

    # Correct -0.0 to +0.0
    self._xs += 0.

    # Reverse array so that it is ordered intuitively from high to low energy
    self._xs = self._xs[..., ::-1]



class ScatterXS(MultiGroupXS):

  def __init__(self, domain=None, domain_type=None, energy_groups=None):
    super(ScatterXS, self).__init__(domain, domain_type, energy_groups)
    self._xs_type = 'scatter'

  def createTallies(self):

    super(ScatterXS, self).createTallies()

    # Create Filter objects
    group_edges = self._energy_groups._group_edges
    energy_filter = openmc.Filter('energy', group_edges)
    domain_filter = openmc.Filter(self._domain_type, self._domain._id)

    # Build a Tally label string
    label = '{0} groups'.format(self._num_groups)

    flux = openmc.Tally()
    flux.addScore('flux')
    flux.addFilter(domain_filter)
    flux.addFilter(energy_filter)
    flux.setLabel(label)
    self._tallies['flux'] = flux

    scatter = openmc.Tally()
    scatter.addScore('scatter')
    scatter.addFilter(domain_filter)
    scatter.addFilter(energy_filter)
    scatter.setLabel(label)
    self._tallies['scatter'] = scatter


  def computeXS(self):

    # TODO: Error checking...

    # Find the flux, scatter Tallies
    flux_tally = self._tallies['flux']
    scatter_tally = self._tallies['scatter']

    flux_mean = flux_tally._mean
    scatter_mean = scatter_tally._mean
    flux_std_dev = flux_tally._std_dev
    scatter_std_dev = scatter_tally._std_dev

    # Get the number of domains
    num_subdomains = flux_tally._mean.shape[0] / self._num_groups
    new_shape = (num_subdomains, self._num_groups)

#    flux_mean = np.reshape(flux_mean, new_shape)
#    scatter_mean = np.reshape(scatter_mean, new_shape)
#    flux_std_dev = np.reshape(flux_std_dev, new_shape)
#    scatter_std_dev = np.reshape(scatter_std_dev, new_shape)

    flux = np.array([flux_mean, flux_std_dev])
    scatter = np.array([scatter_mean, scatter_std_dev])

    # Set any zero fluxes and reaction rates to a negative value
    flux_indices = flux[0, ...] == 0.
    scatter_indices = scatter[0, ...] == 0.
    flux[:, flux_indices] = -1.
    scatter[:, scatter_indices] = -1.

    # Compute the xs with uncertainty propagation
    self._xs = uncorr_divide_by_array(scatter, flux)

    # For any region without flux or reaction rate convert xs to zero
    self._xs[:, flux_indices] = 0.
    self._xs[:, scatter_indices] = 0.

    # Correct -0.0 to +0.0
    self._xs += 0.

    # Reverse array so that it is ordered intuitively from high to low energy
    self._xs = self._xs[..., ::-1]



class NuScatterXS(MultiGroupXS):

  def __init__(self, domain=None, domain_type=None, energy_groups=None):
    super(NuScatterXS, self).__init__(domain, domain_type, energy_groups)
    self._xs_type = 'nu-scatter'

  def createTallies(self):

    super(NuScatterXS, self).createTallies()

    # Create Filter objects
    group_edges = self._energy_groups._group_edges
    energy_filter = openmc.Filter('energy', group_edges)
    domain_filter = openmc.Filter(self._domain_type, self._domain._id)

    # Build a Tally label string
    label = '{0} groups'.format(self._num_groups)

    flux = openmc.Tally()
    flux.addScore('flux')
    flux.addFilter(domain_filter)
    flux.addFilter(energy_filter)
    flux.setLabel(label)
    flux.setEstimator('analog')
    self._tallies['flux'] = flux

    nu_scatter = openmc.Tally()
    nu_scatter.addScore('nu-scatter')
    nu_scatter.addFilter(domain_filter)
    nu_scatter.addFilter(energy_filter)
    nu_scatter.setLabel(label)
    nu_scatter.setEstimator('analog')
    self._tallies['nu-scatter'] = nu_scatter


  def computeXS(self):

    # TODO: Error checking...

    # Find the flux, nu-scatter Tallies
    flux_tally = self._tallies['flux']
    nu_scatter_tally = self._tallies['nu-scatter']

    flux_mean = flux_tally._mean
    nu_scatter_mean = nu_scatter_tally._mean
    flux_std_dev = flux_tally._std_dev
    nu_scatter_std_dev = nu_scatter_tally._std_dev

    # Get the number of domains
    num_subdomains = flux_tally._mean.shape[0] / self._num_groups
    new_shape = (num_subdomains, self._num_groups)

#    flux_mean = np.reshape(flux_mean, new_shape)
#    nu_scatter_mean = np.reshape(nu_scatter_mean, new_shape)
#    flux_std_dev = np.reshape(flux_std_dev, new_shape)
#    nu_scatter_std_dev = np.reshape(nu_scatter_std_dev, new_shape)

    flux = np.array([flux_mean, flux_std_dev])
    nu_scatter = np.array([nu_scatter_mean, nu_scatter_std_dev])

    # Set any zero fluxes and reaction rates to a negative value
    flux_indices = flux[0, ...] == 0.
    nu_scatter_indices = nu_scatter[0, ...] == 0.
    flux[:, flux_indices] = -1.
    nu_scatter[:, nu_scatter_indices] = -1.

    # Compute the xs with uncertainty propagation
    self._xs = uncorr_divide_by_array(nu_scatter, flux)

    # For any region without flux or reaction rate, convert xs to zero
    self._xs[:, flux_indices] = 0.
    self._xs[:, nu_scatter_indices] = 0.

    # Correct -0.0 to +0.0
    self._xs += 0.

    # Reverse array so that it is ordered intuitively from high to low energy
    self._xs = self._xs[..., ::-1]



class ScatterMatrixXS(MultiGroupXS):

  def __init__(self, domain=None, domain_type=None, energy_groups=None):
    super(ScatterMatrixXS, self).__init__(domain, domain_type, energy_groups)
    self._xs_type = 'scatter matrix'

  def createTallies(self):

    super(ScatterMatrixXS, self).createTallies()

    # Create Filter objects
    group_edges = self._energy_groups._group_edges
    energy_filter = openmc.Filter('energy', group_edges)
    energyout_filter = openmc.Filter('energyout', group_edges)
    domain_filter = openmc.Filter(self._domain_type, self._domain._id)

    # Build a Tally label string
    label = '{0} groups'.format(self._num_groups)

    flux = openmc.Tally()
    flux.addScore('flux')
    flux.addFilter(domain_filter)
    flux.addFilter(energy_filter)
    flux.setLabel(label)
    flux.setEstimator('analog')
    self._tallies['flux'] = flux

    nu_scatter = openmc.Tally()
    nu_scatter.addScore('nu-scatter')
    nu_scatter.addFilter(domain_filter)
    nu_scatter.addFilter(energy_filter)
    nu_scatter.addFilter(energyout_filter)
    nu_scatter.setLabel(label)
    nu_scatter.setEstimator('analog')
    self._tallies['nu-scatter'] = nu_scatter


  def computeXS(self):

    # TODO: Error checking...

    # Find the flux, nu-scatter Tallies
    flux_tally = self._tallies['flux']
    nu_scatter_tally = self._tallies['nu-scatter']

    flux_mean = flux_tally._mean
    nu_scatter_mean = nu_scatter_tally._mean
    flux_std_dev = flux_tally._std_dev
    nu_scatter_std_dev = nu_scatter_tally._std_dev

    # Get the number of domains
    num_subdomains = flux_tally._mean.shape[0] / self._num_groups
    new_shape = (num_subdomains, self._num_groups, self._num_groups)

#    flux_mean = np.reshape(flux_mean, new_shape[0:2])
#    nu_scatter_mean = np.reshape(nu_scatter_mean, new_shape)
#    flux_std_dev = np.reshape(flux_std_dev, new_shape[0:2])
#    nu_scatter_std_dev = np.reshape(nu_scatter_std_dev, new_shape)

    flux = np.array([flux_mean, flux_std_dev])
    nu_scatter = np.array([nu_scatter_mean, nu_scatter_std_dev])

    # Set any zero fluxes or reaction rates to a negative value
    flux_indices = flux[0, ...] == 0.
    nu_scatter_indices = nu_scatter[0, ...] == 0.
    flux[:, flux_indices] = -1.
    nu_scatter[:, nu_scatter_indices] = -1.

    # Tile the flux into a 3D array corresponding to the nu-scatter array
#    flux = np.reshape(flux, (2, num_subdomains, self._num_groups, 1))
    # FIXME: should axis=2 or axis=3?
#    flux = np.repeat(flux, self._num_groups, axis=3)

    # Compute the xs with uncertainty propagation
#    self._xs = uncorr_divide_by_array(nu_scatter, flux)

    # For any region without flux or reaction rate, convert xs to zero
#    self._xs[:, flux_indices] = 0.
#    self._xs[:, nu_scatter_indices] = 0.

    # Correct -0.0 to +0.0
#    self._xs += 0.

    # Reverse array so that it is ordered intuitively from high to low energy
#    self._xs = self._xs[..., ::-1, ::-1]


  def printXS(self, subdomain=None):

    if not subdomain is None and not subdomain in self._subdomain_offsets.keys():
      msg = 'Unable to print cross-section for domain {0} since it is ' \
            'not one of the domain offsets'.format(subdomain)
      raise ValueError(msg)

    elif subdomain is None and self._domain_type == 'distribcell':
      msg = 'Unable to print cross-section for a distribcell ' \
            'since no subdomain was provided'
      raise ValueError(msg)

    elif self._xs is None:
      msg = 'Unable to print cross-section {0} since it has not yet ' \
            'been computed'.format(self._xs_type)
      raise ValueError(msg)

    string = 'Multi-Group XS\n'
    string += '{0: <16}{1}{2}\n'.format('\tType', '=\t', self._xs_type)
    string += '{0: <16}{1}{2}\n'.format('\tDomain Type', '=\t', self._domain_type)
    string += '{0: <16}{1}{2}\n'.format('\tDomain ID', '=\t', self._domain._id)

    # Compute the subdomain index into the array
    if self._domain_type == 'distribcell':
      subdomain_index = self._subdomain_offsets[subdomain]
    else:
      subdomain_index = 0

    string += '{0: <16}\n'.format('\tEnergy Groups:')

    # Loop over energy groups ranges
    for group in range(1,self._num_groups+1):
      bounds = self._energy_groups.getGroupBounds(group)
      string += '{0: <12}Group {1} [{2: <10} - ' \
                '{3: <10}MeV]\n'.format('', group, bounds[0], bounds[1])

    string += '{0: <16}\n'.format('\tCross-Sections [cm^-1]:')

    # Loop over energy groups ranges
    for in_group in range(1,self._num_groups+1):
      for out_group in range(1,self._num_groups+1):
        string += '{0: <12}Group {1} -> Group {2}:\t\t'.format('', in_group, out_group)
        average = self._xs[0, subdomain_index, in_group-1, out_group-1]
        std_dev = self._xs[1, subdomain_index, in_group-1, out_group-1]
        string += '{:.2e}+/-{:.2e}'.format(average, std_dev)
        string += '\n'

    print(string)


  def getXS(self, in_group, out_group, subdomain=None, metric='mean'):

    if not is_integer(in_group) or not is_integer(out_group):
      msg = 'Unable to get cross-section for non-integer in-scatter' \
            'group {0} and out-scatter group {1}'.format(in_group, out_group)
      raise ValueError(msg)

    elif in_group < 1 or in_group > self._num_groups:
      msg = 'Unable to get cross-section for non-integer in-scatter group ' \
            '{0} which is oustide the energy group bound'.format(in_group)
      raise ValueError(msg)

    elif out_group < 1 or out_group > self._num_groups:
      msg = 'Unable to get cross-section for non-integer out-scatter group ' \
            '{0} which is oustide the energy group bound'.format(out_group)
      raise ValueError(msg)

    if not subdomain is None and not subdomain in self._subdomain_offsets.keys():
      msg = 'Unable to get cross-section for domain {0} since it is ' \
            'not one of the domain offsets'.format(subdomain)
      raise ValueError(msg)

    elif subdomain is None and self._domain_type == 'distribcell':
      msg = 'Unable to get cross-section for a distribcell ' \
            'since no subdomain was provided'
      raise ValueError(msg)

    elif not metric in ['mean', 'std_dev']:
      msg = 'Unable to get cross-section for metric {0} which is not ' \
            '\'mean\' or \'std_dev\''.format(metric)
      raise ValueError(msg)

    elif self._xs is None:
      msg = 'Unable to get cross-section since it has not been computed'
      raise ValueError(msg)

    # Compute the in/out-scatter energy group indices into the array
    # Subtract one from group number
    in_group_index = in_group - 1
    out_group_index = out_group - 1

    # Compute the subdomain index into the array
    if self._domain_type == 'distribcell':
      subdomain_index = self._subdomain_offsets[subdomain]
    else:
      subdomain_index = 0

    if metric == 'mean':
      return self._xs[0, subdomain_index, in_group_index, out_group_index]
    else:
      return self._xs[1, subdomain_index, in_group_index, out_group_index]



class DiffusionCoeff(MultiGroupXS):

  def __init__(self, domain=None, domain_type=None, energy_groups=None):
    super(DiffusionCoeff, self).__init__(domain, domain_type, energy_groups)
    self._xs_type = 'diffusion'

  def createTallies(self):

    super(DiffusionCoeff, self).createTallies()

    msg = 'DiffusionCoeff is not yet able to build tallies'
    raise ValueError(msg)

  def computeXS(self):

    # TODO: Error checking...
    msg = 'Unable to compute diffusion coeff'
    raise ValueError(msg)




class Chi(MultiGroupXS):

  def __init__(self, domain=None, domain_type=None, energy_groups=None):
    super(Chi, self).__init__(domain, domain_type, energy_groups)
    self._xs_type = 'chi'

  def createTallies(self):

    super(Chi, self).createTallies()

    # Create Filter objects
    group_edges = self._energy_groups._group_edges
    energy_filter = openmc.Filter('energy', group_edges)
    energyout_filter = openmc.Filter('energyout', group_edges)
    domain_filter = openmc.Filter(self._domain_type, self._domain._id)

    # Build a Tally label string
    label = '{0} groups'.format(self._num_groups)

    nu_fission_in = openmc.Tally()
    nu_fission_in.addScore('nu-fission')
    nu_fission_in.addFilter(domain_filter)
    nu_fission_in.addFilter(energy_filter)
    nu_fission_in.setEstimator('analog')
    nu_fission_in.setLabel(label)
    self._tallies['nu-fission-in'] = nu_fission_in

    nu_fission_out = openmc.Tally()
    nu_fission_out.addScore('nu-fission')
    nu_fission_out.addFilter(domain_filter)
    nu_fission_out.addFilter(energyout_filter)
    nu_fission_out.setEstimator('analog')
    nu_fission_out.setLabel(label)
    self._tallies['nu-fission-out'] = nu_fission_out


  def computeXS(self):

    # TODO: Error checking...

    # Find the nu-fission Tallies
    nu_fission_in_tally = self._tallies['nu-fission-in']
    nu_fission_out_tally = self._tallies['nu-fission-out']

    nu_fission_in_mean = nu_fission_in_tally._mean
    nu_fission_out_mean = nu_fission_out_tally._mean
    nu_fission_in_std_dev = nu_fission_in_tally._std_dev
    nu_fission_out_std_dev = nu_fission_out_tally._std_dev

    # Get the number of domains
    num_subdomains = nu_fission_in_tally._mean.shape[0] / self._num_groups
    new_shape = (num_subdomains, self._num_groups)

#    nu_fission_in_mean = np.reshape(nu_fission_in_mean, new_shape)
#    nu_fission_out_mean = np.reshape(nu_fission_out_mean, new_shape)
#    nu_fission_in_std_dev = np.reshape(nu_fission_in_std_dev, new_shape)
#    nu_fission_out_std_dev = np.reshape(nu_fission_out_std_dev, new_shape)

    nu_fission_in = np.array([nu_fission_in_mean, nu_fission_in_std_dev])
    nu_fission_out = np.array([nu_fission_out_mean, nu_fission_out_std_dev])

    # Compute the xs with uncertainty propagation
    norm = nu_fission_in[:].sum()

    # If the xs are from a non-fissionable domain
    if norm == 0.:
      self._xs = np.zeros(nu_fission_in.shape)

    # If the xs are from a fissionable domain
    else:

      # Set any zero reaction rates to a negative value
      nu_fission_in_indices = nu_fission_in[0, :, :] == 0.
      nu_fission_in[0, nu_fission_in_indices] = -1.

      # FIXME - uncertainty propagation
      self._xs = uncorr_divide_by_scalar(nu_fission_out,
                                         nu_fission_in.sum(2)[0, :, np.newaxis])

      # Compute the total across all groups per subdomain
      norm = self._xs.sum(2)[0, :, np.newaxis]

      # Set any norms to a negative value
      norm_indices = norm == 0.
      norm[norm_indices] = -1.

      # Normalize chi to 1.0
      # FIXME - uncertainty propagation
      self._xs = uncorr_divide_by_scalar(self._xs, norm)

      # For any region without flux or reaction rate, convert xs to zero
      self._xs[..., norm_indices] = 0.

      # Correct -0.0 to +0.0
      self._xs += 0.

      # Reverse array so that it is ordered intuitively from high to low energy
      self._xs = self._xs[..., ::-1]