import openmc
import infermc
import numpy as np
import os, abc

# Type-checking support
from typecheck import accepts, Or, Exact, Self

# Supported cross-section types
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

# Supported domain types
domain_types = ['cell',
                'distribcell',
                'universe',
                'material',
                'mesh']

# Supported domain objects
domains = [openmc.Cell, openmc.Universe, openmc.Material, openmc.Mesh]

# Type collections for use with the typecheck module
xs_types_check = Or(*[Exact(xs_type) for xs_type in xs_types])
domain_types_check = Or(*[Exact(domain_type) for domain_type in domain_types])
domains_check = Or(*[domain for domain in domains])

# LaTeX Greek symbols for each cross-section type
greek = dict()
greek['total'] = '$\\Sigma_{t}$'
greek['transport'] = '$\\Sigma_{tr}$'
greek['absorption'] = '$\\Sigma_{a}$'
greek['scatter'] = '$\\Sigma_{s}$'
greek['nu-scatter'] = '$\\nu\\Sigma_{s}$'
greek['scatter matrix'] = '$\\Sigma_{s}$'
greek['fission'] = '$\\Sigma_{f}$'
greek['nu-fission'] = '$\\nu\\Sigma_{f}$'
greek['chi'] = '$\\chi$'
greek['diffusion'] = '$D$'

# Tally batch metrics
metrics = dict()
metrics['mean'] = 0
metrics['std_dev'] = 1


def flip_axis(arr, axis=0):
    ''' Flip contents of `axis` in array `arr`
    Taken verbatim from:
    https://github.com/nipy/nibabel/blob/master/nibabel/orientations.py
    '''
    arr = np.asanyarray(arr)
    arr = arr.swapaxes(0, axis)
    arr = np.flipud(arr)
    return arr.swapaxes(axis, 0)


class MultiGroupXS(object):

  # This is an abstract class which cannot be instantiated
  __metaclass__ = abc.ABCMeta

  def __init__(self, domain=None, domain_type=None, energy_groups=None):

    self._xs_type = None
    self._domain = None
    self._domain_type = None
    self._energy_groups = None
    self._num_groups = None
    self._tallies = dict()
    self._xs = None

    # A dictionary used to compute indices into the xs array
    # Keys   - Domain ID (ie, Material ID, Region ID for districell, etc)
    # Values - Offset/stride into xs array
    self._subdomain_offsets = dict()
    self._offset = None

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
  @accepts(Self(), infermc.EnergyGroups)
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

    if self._domain_type in ['material', 'cell', 'universe', 'mesh']:
      self._subdomain_offsets[domain._id] = 0


  @property
  def domain_type(self):
    return self._domain_type


  @domain_type.setter
  @accepts(Self(), domain_types_check)
  def domain_type(self, domain_type):
    self._domain_type = domain_type


  def findDomainOffset(self):

    if self._domain_type in ['material', 'cell', 'universe', 'mesh']:
      self._offset = 0

    # Distribcell tallies
    else:
      first_tally = self._tallies[self._tallies.keys()[0]]
      domain_filter = first_tally.findFilter('distribcell', [self._domain._id])
      self._offset = domain_filter._offset


  @accepts(Self(), int, int)
  def setSubDomainOffset(self, domain_id, offset):
    self._subdomain_offsets[domain_id] = offset


  @abc.abstractmethod
  @accepts(Self(), [str], [[openmc.Filter]], [str],
           Or(Exact('analog'), Exact('tracklength')), )
  def createTallies(self, scores, filters, keys, estimator):

    if self._energy_groups is None:
      msg = 'Unable to create Tallies since the energy groups has not been set'
      raise ValueError(msg)

    elif self._domain is None:
      msg = 'Unable to create Tallies since the domain has not been set'
      raise ValueError(msg)

    elif self._domain_type is None:
      msg = 'Unable to create Tallies since the domain type has not been set'
      raise ValueError(msg)

    # Create a domain Filter object
    domain_filter = openmc.Filter(self._domain_type, self._domain._id)

    # Create a label for the Tallies
    label = '{0} groups'.format(self._num_groups)

    for i, score in enumerate(scores):
      key = keys[i]
      self._tallies[key] = openmc.Tally(label=label)
      self._tallies[key].addScore(score)
      self._tallies[key].setEstimator(estimator)
      self._tallies[key].addFilter(domain_filter)

      # Add all non-domain specific Filters (ie, energy) to the Tally
      for filter in filters[i]:
        self._tallies[key].addFilter(filter)


  @abc.abstractmethod
  def computeXS(self):

    if self._tallies is None:
      msg = 'Unable to compute cross-section without any Tallies'
      raise ValueError(msg)

    tally_data = dict()
    zero_indices = dict()

    for key, tally in self._tallies.items():

      # Get the Tally batch mean and std. dev.
      mean = tally._mean
      std_dev = tally._std_dev

      # Determine shape of the Tally data from its Filters, Nuclides
      new_shape = tuple()
      energy_axes = list()

      for i, filter in enumerate(tally._filters):
        new_shape += (filter.getNumBins(), )

        if 'energy' in filter._type:
          energy_axes.append(i)

      new_shape += (tally.getNumNuclides(), )

      # Reshape the array
      mean = np.reshape(mean, new_shape)
      std_dev = np.reshape(std_dev, new_shape)

      # Reverse arrays so they are ordered from high to low energy
      for energy_axis in energy_axes:
        mean = flip_axis(mean, axis=energy_axis)
        std_dev = flip_axis(std_dev, axis=energy_axis)

      # Get the array of indices of zero elements
      zero_indices[key] = mean[...] == 0.

      # Store the tally data to the dictionary
      tally_data[key] = np.array([mean, std_dev])

    return tally_data, zero_indices


  @accepts(Self(), Or(str, tuple, list, np.ndarray))
  def getSubdomainOffsets(self, subdomains='all'):

    if subdomains == 'all':
      offsets = self._subdomain_offsets.values()

    else:
      offsets = np.zeros(len(subdomains), dtype=np.int64)

      for i, subdomain in enumerate(subdomains):
        if self._subdomain_offsets.has_key(subdomain):
          offsets[i] = self._subdomain_offsets[subdomain]
        else:
          msg = 'Unable to get subdomain index for subdomain {0} since it is ' \
                'not one of the subdomains in the cross-section'.format(subdomain)
          raise ValueError(msg)

    return offsets


  @accepts(Self(), Or(str, tuple, list, np.ndarray))
  def getSubdomains(self, offsets='all'):

    if offsets == 'all':
      subdomains = self._subdomain_offsets.keys()

    else:
      subdomains = np.zeros(len(offsets), dtype=np.int64)

      keys = self._subdomain_offsets.keys()
      values = self._subdomain_offsets.values()

      for i, offset in enumerate(offsets):
        if offset in values:
          subdomains[i] = keys[values.index(offset)]
        else:
          msg = 'Unable to get subdomain for offset {0} since it is ' \
                'not one of the offsets in the cross-section'.format(offset)
          raise ValueError(msg)

    return subdomains


  @accepts(Self(), Or(str, tuple, list, np.ndarray),
           Or(str, tuple, list, np.ndarray), str)
  def getXS(self, groups='all', subdomains='all', metric='mean'):

    if self._xs is None:
      msg = 'Unable to get cross-section since it has not been computed'
      raise ValueError(msg)

    global metrics
    groups = self._energy_groups.getGroupIndices(groups)
    offsets = self.getSubdomainOffsets(subdomains)

    xs = self._xs[metrics[metric], offsets, ...]
    xs = xs[..., groups, :]
    return xs


  @accepts(Self(), Or(str, tuple, list, np.ndarray))
  def printXS(self, subdomains='all'):

    string = 'Multi-Group XS\n'
    string += '{0: <16}{1}{2}\n'.format('\tType', '=\t', self._xs_type)
    string += '{0: <16}{1}{2}\n'.format('\tDomain Type', '=\t', self._domain_type)
    string += '{0: <16}{1}{2}\n'.format('\tDomain ID', '=\t', self._domain._id)

    if subdomains == 'all':
      subdomains = self._subdomain_offsets.keys()

    # Loop over all subdomains
    for subdomain in subdomains:

      if self._domain_type == 'distribcell':
        string += '{0: <16}{1}{2}\n'.format('\tSubdomain', '=\t', subdomain)

      string += '{0: <16}\n'.format('\tCross-Sections [cm^-1]:')

      # Loop over energy groups ranges
      for group in range(1,self._num_groups+1):
        bounds = self._energy_groups.getGroupBounds(group)
        string += '{0: <12}Group {1} [{2: <10} - ' \
                  '{3: <10}MeV]:\t'.format('', group, bounds[0], bounds[1])
        average = self.getXS([group], [subdomain], 'mean')
        std_dev = self.getXS([group], [subdomain], 'std_dev')
        string += '{:.2e}+/-{:.2e}'.format(average[0,0], std_dev[0,0])
        string += '\n'

      string += '\n'

    print(string)


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
  def restoreFromFile(self, filename='multigroupxs', directory='multigroupxs'):

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


  @accepts(Self(), Or(str, tuple, list, np.ndarray), str, str, str, bool)
  def exportResults(self, subdomains='all', filename='multigroupxs',
                    directory='multigroupxs', format='hdf5', append=True):

    # Make directory if it does not exist
    if not os.path.exists(directory):
      os.makedirs(directory)

    global metrics
    offsets = self.getSubdomainOffsets(subdomains)
    subdomains = self.getSubdomains(offsets)
    average = self._xs[metrics['mean'], offsets, ...]
    std_dev = self._xs[metrics['std_dev'], offsets, ...]

    # HDF5 binary file
    if format == 'hdf5':

      import h5py

      filename = directory + '/' + filename + '.h5'
      filename = filename.replace(' ', '-')

      if append:
        xs_results = h5py.File(filename, 'a')
      else:
        xs_results = h5py.File(filename, 'w')

      # Create an HDF5 group within the file for the domain
      domain_type_group = xs_results.require_group(self._domain_type)
      group_name = '{0} {1}'.format(self._domain_type, self._domain._id)
      domain_group = domain_type_group.require_group(group_name)

      # Loop over all subdomains
      for i, subdomain in enumerate(subdomains):

        # Create an HDF5 group for the subdomain and xs type
        if self._domain_type == 'distribcell':
          group_name = '{0}'.format(subdomain)
          subdomain_group = domain_group.require_group(group_name)
        else:
          subdomain_group = domain_group

        xs_group = subdomain_group.require_group(self._xs_type)

        # Add MultiGroupXS results data to the HDF5 group
        xs_group.require_dataset('average', dtype=np.float64,
                                 shape=average[i, ...].shape,
                                 data=average[i, ...])
        xs_group.require_dataset('std. dev.', dtype=np.float64,
                                 shape=std_dev[i, ...].shape,
                                 data=std_dev[i, ...])

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

      group_name = '{0} {1}'.format(self._domain_type, self._domain._id)

      if not xs_results.has_key(group_name):
        domain_group = xs_results[group_name] = dict()
      else:
        domain_group = xs_results[group_name]

      # Loop over all subdomains
      for i, subdomain in enumerate(subdomains):

        # Create an HDF5 group for the subdomain and xs type
        if self._domain_type == 'distribcell':
          group_name = '{0}'.format(subdomain)

          if not domain_group.has_key(group_name):
            subdomain_group = domain_group[group_name] = dict()
          else:
            subdomain_group = domain_group[group_name]

        else:
          subdomain_group = domain_group

        xs_group = subdomain_group[self._xs_type] = dict()

        # Add MultiGroupXS results data to the dictionary
        xs_group['average'] = average[i, ...]
        xs_group['std. dev.'] = std_dev[i, ...]

      # Pickle the MultiGroupXS results to a file
      pickle.dump(xs_results, open(filename, 'wb'))


    # LaTeX script
    elif format == 'latex':

      import tabulate
      global greek

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

      # Loop over all subdomains
      for i, subdomain in enumerate(subdomains):

        # Add MultiGroupXS results data to the table list
        table = list()

        if self._xs_type != 'scatter matrix':
          headers = list()
          headers.append('Group')
          headers.append('Average XS')
          headers.append('Std. Dev.')

          for group in range(self._num_groups):
            subtable = list()
            subtable.append(group+1)
            subtable.append(average[i, group, ...])
            subtable.append(std_dev[i, group, ...])
            table.append(subtable)

        # Scattering matrix
        else:
          headers = list()
          headers.append('Group In')
          headers.append('Group Out')
          headers.append('Average XS')
          headers.append('Std. Dev.')

          for in_group in range(self._num_groups):
            for out_group in range(self._num_groups):
              subtable = list()
              subtable.append(in_group+1)
              subtable.append(out_group+1)
              subtable.append(average[i, in_group, out_group, ...])
              subtable.append(std_dev[i, in_group, out_group, ...])
              table.append(subtable)

        if self._domain_type == 'distribcell':
          caption = '\\caption{{{0} {1}, (subdomain {2}) {3} [cm$^-1]}}'.format(
            self._domain_type.capitalize(), self._domain._id,
            subdomain, greek[self._xs_type])
        else:
          caption = '\\caption{{{0} {1} {2} [cm^-1]}}'.format(
            self._domain_type.capitalize(), subdomain, greek[self._xs_type])

        # Write the MultiGroupXS results to a file
        xs_results.write('\\begin{table}\n')
        xs_results.write('\\begin{center}\n')
        xs_results.write(tabulate.tabulate(table, headers, tablefmt='latex'))
        xs_results.write(caption)
        xs_results.write('\\end{center}\n')
        xs_results.write('\\end{table}\n')
        xs_results.write('\\\n')

      xs_results.write('\\end{document}\n')
      xs_results.close()


  @accepts(Self(), Or(str, tuple, list, np.ndarray), str, str)
  def printPDF(self, subdomains='all', filename='multigroupxs',
               directory='multigroupxs'):

    import subprocess

    # Make directory if it does not exist
    if not os.path.exists(directory):
      os.makedirs(directory)

    filename = filename.replace(' ', '-')

    # Generate LaTeX file
    self.exportResults(subdomains, filename, '.', 'latex', False)

    # Compile LaTeX to PDF
    FNULL = open(os.devnull, 'w')
    subprocess.check_call('pdflatex {0}.tex'.format(filename),
                          shell=True, stdout=FNULL)

    # Move PDF to requested directory and cleanup temporary LaTeX files
    if directory != '.':
      os.system('mv {0}.pdf {1}'.format(filename, directory))

    os.system('rm {0}.tex {0}.aux {0}.log'.format(filename))


class TotalXS(MultiGroupXS):

  def __init__(self, domain=None, domain_type=None, energy_groups=None):
    super(TotalXS, self).__init__(domain, domain_type, energy_groups)
    self._xs_type = 'total'


  def createTallies(self):

    # Create a list of scores for each Tally to be created
    scores = ['flux', 'total']
    estimator = 'tracklength'
    keys = scores

    # Create the non-domain specific Filters for the Tallies
    group_edges = self._energy_groups._group_edges
    energy_filter = openmc.Filter('energy', group_edges)
    filters = [[energy_filter], [energy_filter]]

    # Intialize the Tallies
    super(TotalXS, self).createTallies(scores, filters, keys, estimator)


  def computeXS(self):

    # Extract and clean the Tally data
    tally_data, zero_indices = super(TotalXS, self).computeXS()
    total = tally_data['total']
    flux = tally_data['flux']

    # Set any subdomain's zero fluxes and reaction rates to a negative value
    flux[:, zero_indices['flux']] = -1.
    total[:, zero_indices['total']] = -1.

    # Compute the xs with uncertainty propagation
    self._xs = infermc.uncorr_math.divide_by_array(total, flux)

    # For any region without flux or reaction rate, convert xs to zero
    all_zero_indices = np.logical_or(zero_indices['flux'],
                                     zero_indices['total'])
    self._xs[:, all_zero_indices] = 0.

    # Correct -0.0 to +0.0
    self._xs += 0.


class TransportXS(MultiGroupXS):

  def __init__(self, domain=None, domain_type=None, energy_groups=None):
    super(TransportXS, self).__init__(domain, domain_type, energy_groups)
    self._xs_type = 'transport'


  def createTallies(self):

    # Create a list of scores for each Tally to be created
    scores = ['flux', 'total', 'scatter-1']
    estimator = 'analog'
    keys = scores

    # Create the non-domain specific Filters for the Tallies
    group_edges = self._energy_groups._group_edges
    energy_filter = openmc.Filter('energy', group_edges)
    filters = [[energy_filter], [energy_filter], [energy_filter]]

    # Intialize the Tallies
    super(TransportXS, self).createTallies(scores, filters, keys, estimator)


  def computeXS(self):

    # Extract and clean the Tally data
    tally_data, zero_indices = super(TransportXS, self).computeXS()
    total = tally_data['total']
    scatter1 = tally_data['scatter-1']
    flux = tally_data['flux']

    # Set any zero fluxes to a negative value
    flux[:, zero_indices['flux']] = -1.

    delta = infermc.uncorr_math.sub(total, scatter1)

    # Set any subdomain's zero reaction rates to a negative value
    delta_indices = delta[0, ...] == 0.
    delta[:, delta_indices] = -1.

    # Compute the xs with uncertainty propagation
    self._xs = infermc.uncorr_math.divide_by_array(delta, flux)

    # For any region without flux or reaction rate, convert xs to zero
    all_zero_indices = np.logical_or(zero_indices['flux'], delta_indices)
    self._xs[:, all_zero_indices] = 0.

    # Correct -0.0 to +0.0
    self._xs += 0.


class AbsorptionXS(MultiGroupXS):

  def __init__(self, domain=None, domain_type=None, energy_groups=None):
    super(AbsorptionXS, self).__init__(domain, domain_type, energy_groups)
    self._xs_type = 'absorption'


  def createTallies(self):

    # Create a list of scores for each Tally to be created
    scores = ['flux', 'absorption']
    estimator = 'tracklength'
    keys = scores

    # Create the non-domain specific Filters for the Tallies
    group_edges = self._energy_groups._group_edges
    energy_filter = openmc.Filter('energy', group_edges)
    filters = [[energy_filter], [energy_filter]]

    # Intialize the Tallies
    super(AbsorptionXS, self).createTallies(scores, filters, keys, estimator)


  def computeXS(self):

    # Extract and clean the Tally data
    tally_data, zero_indices = super(AbsorptionXS, self).computeXS()
    absorption = tally_data['absorption']
    flux = tally_data['flux']

    # Set any subdomain's zero fluxes and reaction rates to a negative value
    flux[:, zero_indices['flux']] = -1.
    absorption[:, zero_indices['absorption']] = -1.

    # Compute the xs with uncertainty propagation
    self._xs = infermc.uncorr_math.divide_by_array(absorption, flux)

    # For any region without flux or reaction rate, convert xs to zero
    all_zero_indices = np.logical_or(zero_indices['flux'],
                                     zero_indices['absorption'])
    self._xs[:, all_zero_indices] = 0.

    # Correct -0.0 to +0.0
    self._xs += 0.


class FissionXS(MultiGroupXS):

  def __init__(self, domain=None, domain_type=None, energy_groups=None):
    super(FissionXS, self).__init__(domain, domain_type, energy_groups)
    self._xs_type = 'fission'


  def createTallies(self):

    # Create a list of scores for each Tally to be created
    scores = ['flux', 'fission']
    estimator = 'tracklength'
    keys = scores

    # Create the non-domain specific Filters for the Tallies
    group_edges = self._energy_groups._group_edges
    energy_filter = openmc.Filter('energy', group_edges)
    filters = [[energy_filter], [energy_filter]]

    # Intialize the Tallies
    super(FissionXS, self).createTallies(scores, filters, keys, estimator)


  def computeXS(self):

    # Extract and clean the Tally data
    tally_data, zero_indices = super(FissionXS, self).computeXS()
    fission = tally_data['fission']
    flux = tally_data['flux']

    # Set any subdomain's zero fluxes and reaction rates to a negative value
    flux[:, zero_indices['flux']] = -1.
    fission[:, zero_indices['fission']] = -1.

    # Compute the xs with uncertainty propagation
    self._xs = infermc.uncorr_math.divide_by_array(fission, flux)

    # For any region without flux or reaction rate, convert xs to zero
    all_zero_indices = np.logical_or(zero_indices['flux'],
                                     zero_indices['fission'])
    self._xs[:, all_zero_indices] = 0.

    # Correct -0.0 to +0.0
    self._xs += 0.


class NuFissionXS(MultiGroupXS):

  def __init__(self, domain=None, domain_type=None, energy_groups=None):
    super(NuFissionXS, self).__init__(domain, domain_type, energy_groups)
    self._xs_type = 'nu-fission'

  def createTallies(self):

    # Create a list of scores for each Tally to be created
    scores = ['flux', 'nu-fission']
    estimator = 'tracklength'
    keys = scores

    # Create the non-domain specific Filters for the Tallies
    group_edges = self._energy_groups._group_edges
    energy_filter = openmc.Filter('energy', group_edges)
    filters = [[energy_filter], [energy_filter]]

    # Intialize the Tallies
    super(NuFissionXS, self).createTallies(scores, filters, keys, estimator)


  def computeXS(self):

    # Extract and clean the Tally data
    tally_data, zero_indices = super(NuFissionXS, self).computeXS()
    nu_fission = tally_data['nu-fission']
    flux = tally_data['flux']

    # Set any subdomain's zero fluxes and reaction rates to a negative value
    flux[:, zero_indices['flux']] = -1.
    nu_fission[:, zero_indices['nu-fission']] = -1.

    # Compute the xs with uncertainty propagation
    self._xs = infermc.uncorr_math.divide_by_array(nu_fission, flux)

    # For any region without flux or reaction rate, convert xs to zero
    all_zero_indices = np.logical_or(zero_indices['flux'],
                                     zero_indices['nu-fission'])
    self._xs[:, all_zero_indices] = 0.

    # Correct -0.0 to +0.0
    self._xs += 0.


class ScatterXS(MultiGroupXS):

  def __init__(self, domain=None, domain_type=None, energy_groups=None):
    super(ScatterXS, self).__init__(domain, domain_type, energy_groups)
    self._xs_type = 'scatter'

  def createTallies(self):

    # Create a list of scores for each Tally to be created
    scores = ['flux', 'scatter']
    estimator = 'tracklength'
    keys = scores

    # Create the non-domain specific Filters for the Tallies
    group_edges = self._energy_groups._group_edges
    energy_filter = openmc.Filter('energy', group_edges)
    filters = [[energy_filter], [energy_filter]]

    # Intialize the Tallies
    super(ScatterXS, self).createTallies(scores, filters, keys, estimator)


  def computeXS(self):

    # Extract and clean the Tally data
    tally_data, zero_indices = super(ScatterXS, self).computeXS()
    scatter = tally_data['scatter']
    flux = tally_data['flux']

    # Set any subdomain's zero fluxes and reaction rates to a negative value
    flux[:, zero_indices['flux']] = -1.
    scatter[:, zero_indices['scatter']] = -1.

    # Compute the xs with uncertainty propagation
    self._xs = infermc.uncorr_math.divide_by_array(scatter, flux)

    # For any region without flux or reaction rate, convert xs to zero
    all_zero_indices = np.logical_or(zero_indices['flux'],
                                     zero_indices['scatter'])
    self._xs[:, all_zero_indices] = 0.

    # Correct -0.0 to +0.0
    self._xs += 0.


class NuScatterXS(MultiGroupXS):

  def __init__(self, domain=None, domain_type=None, energy_groups=None):
    super(NuScatterXS, self).__init__(domain, domain_type, energy_groups)
    self._xs_type = 'nu-scatter'

  def createTallies(self):

    # Create a list of scores for each Tally to be created
    scores = ['flux', 'nu-scatter']
    estimator = 'analog'
    keys = scores

    # Create the non-domain specific Filters for the Tallies
    group_edges = self._energy_groups._group_edges
    energy_filter = openmc.Filter('energy', group_edges)
    filters = [[energy_filter], [energy_filter]]

    # Intialize the Tallies
    super(NuScatterXS, self).createTallies(scores, filters, keys, estimator)


  def computeXS(self):

    # Extract and clean the Tally data
    tally_data, zero_indices = super(NuScatterXS, self).computeXS()
    nu_scatter = tally_data['nu-scatter']
    flux = tally_data['flux']

    # Set any subdomain's zero fluxes and reaction rates to a negative value
    flux[:, zero_indices['flux']] = -1.
    nu_scatter[:, zero_indices['nu-scatter']] = -1.

    # Compute the xs with uncertainty propagation
    self._xs = infermc.uncorr_math.divide_by_array(nu_scatter, flux)

    # For any region without flux or reaction rate, convert xs to zero
    all_zero_indices = np.logical_or(zero_indices['flux'],
                                     zero_indices['nu-scatter'])
    self._xs[:, all_zero_indices] = 0.

    # Correct -0.0 to +0.0
    self._xs += 0.


class ScatterMatrixXS(MultiGroupXS):

  def __init__(self, domain=None, domain_type=None, energy_groups=None):
    super(ScatterMatrixXS, self).__init__(domain, domain_type, energy_groups)
    self._xs_type = 'scatter matrix'


  def createTallies(self):

    # Create a list of scores for each Tally to be created
    scores = ['flux', 'nu-scatter']
    estimator = 'analog'
    keys = scores

    # Create the non-domain specific Filters for the Tallies
    group_edges = self._energy_groups._group_edges
    energy_filter = openmc.Filter('energy', group_edges)
    energyout_filter = openmc.Filter('energyout', group_edges)
    filters = [[energy_filter], [energy_filter, energyout_filter]]

    # Intialize the Tallies
    super(ScatterMatrixXS, self).createTallies(scores, filters, keys, estimator)


  def computeXS(self):

    # Extract and clean the Tally data
    tally_data, zero_indices = super(ScatterMatrixXS, self).computeXS()
    nu_scatter = tally_data['nu-scatter']
    flux = tally_data['flux']

    # Set any subdomain's zero fluxes and reaction rates to a negative value
    flux[:, zero_indices['flux']] = -1.
    nu_scatter[:, zero_indices['nu-scatter']] = -1.

    # Tile the flux to correspond to the nu-scatter array
    flux = np.repeat(flux[:,:,np.newaxis,:,:], self._num_groups, axis=2)

    # Compute the xs with uncertainty propagation
    self._xs = infermc.uncorr_math.divide_by_array(nu_scatter, flux)

    # For any region without flux or reaction rate, convert xs to zero
    all_zero_indices = np.logical_or(zero_indices['flux'][...,np.newaxis],
                                     zero_indices['nu-scatter'])
    self._xs[:, all_zero_indices] = 0.

    # Correct -0.0 to +0.0
    self._xs += 0.


  @accepts(Self(), Or(str, tuple, list, np.ndarray),
           Or(str, tuple, list, np.ndarray),
           Or(str, tuple, list, np.ndarray), str)
  def getXS(self, in_groups='all', out_groups='all',
            subdomains='all', metric='mean'):

    if self._xs is None:
      msg = 'Unable to get cross-section since it has not been computed'
      raise ValueError(msg)

    global metrics
    in_groups = self._energy_groups.getGroupIndices(in_groups)
    out_groups = self._energy_groups.getGroupIndices(out_groups)
    offsets = self.getSubdomainOffsets(subdomains)

    xs = self._xs[metrics[metric], offsets, ...]
    xs = xs[..., in_groups, out_groups, :]
    return xs


  @accepts(Self(), Or(str, tuple, list, np.ndarray))
  def printXS(self, subdomains='all'):

    string = 'Multi-Group XS\n'
    string += '{0: <16}{1}{2}\n'.format('\tType', '=\t', self._xs_type)
    string += '{0: <16}{1}{2}\n'.format('\tDomain Type', '=\t', self._domain_type)
    string += '{0: <16}{1}{2}\n'.format('\tDomain ID', '=\t', self._domain._id)

    string += '{0: <16}\n'.format('\tEnergy Groups:')

    # Loop over energy groups ranges
    for group in range(1,self._num_groups+1):
      bounds = self._energy_groups.getGroupBounds(group)
      string += '{0: <12}Group {1} [{2: <10} - ' \
                '{3: <10}MeV]\n'.format('', group, bounds[0], bounds[1])

    if subdomains == 'all':
      subdomains = self._subdomain_offsets.keys()

    for subdomain in subdomains:

      if self._domain_type == 'distribcell':
        string += '{0: <16}{1}{2}\n'.format('\tSubdomain', '=\t', subdomain)

      string += '{0: <16}\n'.format('\tCross-Sections [cm^-1]:')

      # Loop over energy groups ranges
      for in_group in range(1,self._num_groups+1):
        for out_group in range(1,self._num_groups+1):
          string += '{0: <12}Group {1} -> Group {2}:\t\t'.format('', in_group, out_group)
          average = self.getXS([in_group], [out_group], [subdomain], 'mean')
          std_dev = self.getXS([in_group], [out_group], [subdomain], 'std_dev')
          string += '{:.2e}+/-{:.2e}'.format(average[0,0], std_dev[0,0])
          string += '\n'

      string += '\n'

    print(string)


class DiffusionCoeff(MultiGroupXS):

  def __init__(self, domain=None, domain_type=None, energy_groups=None):
    super(DiffusionCoeff, self).__init__(domain, domain_type, energy_groups)
    self._xs_type = 'diffusion'


  def createTallies(self):
    msg = 'DiffusionCoeff is not yet able to build tallies'
    raise NotImplementedError(msg)


  def computeXS(self):
    msg = 'Unable to compute diffusion coeff'
    raise NotImplementedError(msg)


class Chi(MultiGroupXS):

  def __init__(self, domain=None, domain_type=None, energy_groups=None):
    super(Chi, self).__init__(domain, domain_type, energy_groups)
    self._xs_type = 'chi'


  def createTallies(self):

    # Create a list of scores for each Tally to be created
    scores = ['nu-fission', 'nu-fission']
    estimator = 'analog'
    keys = ['nu-fission-in', 'nu-fission-out']

    # Create the non-domain specific Filters for the Tallies
    group_edges = self._energy_groups._group_edges
    energy_filter = openmc.Filter('energy', group_edges)
    energyout_filter = openmc.Filter('energyout', group_edges)
    filters = [[energy_filter], [energyout_filter]]

    # Intialize the Tallies
    super(Chi, self).createTallies(scores, filters, keys, estimator)


  def computeXS(self):

    # Extract and clean the Tally data
    tally_data, zero_indices = super(Chi, self).computeXS()
    nu_fission_in = tally_data['nu-fission-in']
    nu_fission_out = tally_data['nu-fission-out']

    # Set any zero reaction rates to -1
    nu_fission_in[0, zero_indices['nu-fission-in']] = -1.

    # FIXME - uncertainty propagation
    self._xs = infermc.uncorr_math.divide_by_scalar(nu_fission_out,
                                   nu_fission_in.sum(2)[0, :, np.newaxis, ...])

    # Compute the total across all groups per subdomain
    norm = self._xs.sum(2)[0, :, np.newaxis, ...]

    # Set any zero norms (in non-fissionable domains) to -1
    norm_indices = norm == 0.
    norm[norm_indices] = -1.

    # Normalize chi to 1.0
    # FIXME - uncertainty propagation
    self._xs = infermc.uncorr_math.divide_by_scalar(self._xs, norm)

    # For any region without flux or reaction rate, convert xs to zero
    self._xs[:, norm_indices] = 0.

    # FIXME - uncertainty propagation - this is just a temporary fix
    self._xs[1, ...] = 0.

    # Correct -0.0 to +0.0
    self._xs += 0.