import openmc
from infermc.checkvalue import *
from infermc.uncorr_math import *
from uncertainties import unumpy
import os


xs_types = ['total',
            'transport',
            'absorption',
            'scatter',
            'nu-scatter',
            'scatter matrix',
            'fission',
            'nu-fission',
            'chi']
#            'diffusion'

domain_types = ['cell',
                'distribcell',
                'universe',
                'material']
#               'mesh'


class EnergyGroups(object):

  def __init__(self):

    self._group_edges = None
    self._num_groups = None


  def setGroupEdges(self, edges):

    if not isinstance(edges, (tuple, list, np.ndarray)):
      msg = 'Unable to set energy group edges from {0} which ' \
            'is not a Python tuple/list or NumPy array'.format(edges)
      raise ValueError(msg)

    for i, edge in enumerate(edges):

      if not is_integer(edge) and not is_float(edge):
        msg = 'Unable to set energy group edge = {0} which is ' \
              'not an integer or floating point value'.format(edge)
        raise ValueError(msg)

      if i > 0 and edges[i] <= edges[i-1]:
        msg = 'Unable to set energy group edges with {0} for edge {1} ' \
              'and {2} for edge {3}'.format(edges[i], i, edges[i-1], i-1)
        raise ValueError(msg)

    self._group_edges = np.array(edges)
    self._num_groups = len(edges)-1


  def generateBinEdges(self, start, stop, num_groups, type='linear'):

    if not is_integer(start) and not is_float(start):
      msg = 'Unable to generate energy group edges with start = {0} ' \
            'which is not an integer or floating point value'.format(start)
      raise ValueError(msg)

    if not is_integer(stop) and not is_float(stop):
      msg = 'Unable to generate energy group edges with stop = {0} ' \
            'which is not an integer or floating point value'.format(stop)
      raise ValueError(msg)

    if not is_integer(num_groups):
      msg = 'Unable to generate energy group edges with num groups = {0} ' \
            'which is not an integer value'.format(num_groups)
      raise ValueError(msg)

    if not type in ['linear', 'logarithmic']:
      msg = 'Unable to generate energy group edges with type = {0} which is ' \
            'neither linear nor logarithmic'.format(type)
      raise ValueError(msg)

    if start < 0.:
      msg = 'Unable to generate energy group edges with start = {0} which is ' \
           'a negative value'.format(start)
      raise ValueError(msg)

    if stop < 0.:
      msg = 'Unable to generate energy group edges with stop = {0} which is ' \
            'a negative value'.format(stop)
      raise ValueError(msg)

    if start >= stop:
      msg = 'Unable to generate energy group edges with start = {0} which is ' \
           'greater than stop = {1}'.format(start, stop)
      raise ValueError(msg)

    if num_groups < 0.:
      msg = 'Unable to generate energy group edges with num groups {0} ' \
            'which is a negative value'.format(num_groups)
      raise ValueError(msg)

    self._num_groups = num_groups

    if type == 'linear':
      self._group_edges = np.linspace(start, stop, num_groups+1)

    else:
      self._group_edges = np.logspace(np.log(start), np.log(stop), num_groups+1)


  def getGroup(self, energy):

    # Assumes energy is in eV

    if not is_integer(energy) and not is_float(energy):
      msg = 'Unable to get energy group for energy {0} since it is ' \
            'neither an integer or floating point value'.format(energy)
      raise ValueError(msg)

    if self._group_edges is None:
      msg = 'Unable to get energy group for energy {0} eV since ' \
            'the group edges have not yet been set'.format(energy)
      raise ValueError(msg)

    # Loop over all edges and search for the group for this energy
    for i, edge in enumerate(self._group_edges):

      if energy <= edge:
        return self._num_groups - i

    msg = 'Unable to find energy group for energy {0} eV since it does not ' \
          'lie within the bounds of the energy group edges'.format(energy)
    raise ValueError(msg)


  def getGroupBounds(self, group):

    if not is_integer(group):
      msg = 'Unable to get energy group bounds for group {0} since it ' \
            'is not an integer value'.format(group)
      raise ValueError(msg)

    elif group < 1 or group > self._num_groups:
      msg = 'Unable to get energy group bounds for group {0} since it ' \
            'is outside the range of energy groups'.format(group)
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
      self.setDomainType(domain_type)

    if not domain is None:
      self.setDomain(domain)

    if not energy_groups is None:
      self.setEnergyGroups(energy_groups)


  def setEnergyGroups(self, energy_groups):

    if not isinstance(energy_groups, EnergyGroups):
      msg = 'Unable to set the energy groups to {0} which is not an ' \
            'EnergyGroups object'.format(energy_groups)
      raise ValueError(msg)

    self._energy_groups = energy_groups
    self._num_groups = energy_groups._num_groups


  def setDomain(self, domain):

    if not isinstance(domain, (openmc.Material, openmc.Cell,
                               openmc.Universe, openmc.Mesh)):
      msg = 'Unable to set the domain to {0} which is not an OpenMC ' \
            'Material, Cell, Universe or Mesh object'.format(domain)
      raise ValueError(msg)

    self._domain = domain

    if self._domain_type in ['material', 'cell', 'universe']:
      self._subdomain_offsets[domain._id] = 0


  def setDomainType(self, domain_type):

    if not domain_type in domain_types:
      msg = 'Unable to set the domain type to {0} which is ' \
            'not a supported domain type'.format(domain_type)
      raise ValueError(msg)

    self._domain_type = domain_type


  def findDomainOffset(self):

    if self._domain_type in ['material', 'cell', 'universe']:
      self._offset = 0

    # Distribcell tallies
    else:
      first_tally = self._tallies[self._tallies.keys()[0]]
      domain_filter = first_tally.findFilter('distribcell', [self._domain._id])
      self._offset = domain_filter._offset


  def setSubDomainOffset(self, domain_id, offset):

    if not is_integer(domain_id):
      msg = 'Unable to set the domain offset for domain {0} which is not ' \
            'an integer value'.format(domain_id)
      raise ValueError(msg)

    elif not is_integer(offset):
      msg = 'Unable to set the domain offset for domain {0} to offset {1} ' \
            'which is not an integer'.format(domain_id, offset)
      raise ValueError(msg)

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


  def getXS(self, group, subdomain=None, metric='mean'):

    if not is_integer(group):
      msg = 'Unable to get cross-section for non-integer group {0}'.format(group)
      raise ValueError(msg)

    elif group < 1 or group > self._num_groups:
      msg = 'Unable to get cross-section for non-integer group {0} ' \
            'which is oustide the energy group bound'.format(group)
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

    # Compute the energy group index into the array
    # Subtract one from group number
    group_index = group - 1

    # Compute the subdomain index into the array
    if self._domain_type == 'distribcell':
      subdomain_index = self._subdomain_offsets[subdomain]
    else:
      subdomain_index = 0

    xs = self._xs[subdomain_index, group_index]

    if metric == 'mean':
      return unumpy.nominal_values(xs)
    else:
      return unumpy.std_devs(xs)


  def dumpToFile(self, filename='multigroupxs', directory='multigroupxs'):

    if not is_string(filename):
      msg = 'Unable to dump cross-section to filename={0} ' \
            'since it is not a string'.format(filename)
      raise ValueError(msg)

    elif not is_string(directory):
      msg = 'Unable to dump cross-section to directory={0} ' \
            'since it is not a string'.format(directory)
      raise ValueError(msg)

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



  def restoreFromFile(self, filename, directory='.'):

    if not is_string(filename):
      msg = 'Unable to import cross-section from filename={0} ' \
            'since it is not a string'.format(filename)
      raise ValueError(msg)

    elif not is_string(directory):
      msg = 'Unable to import cross-section from directory={0} ' \
            'since it is not a string'.format(directory)
      raise ValueError(msg)

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
    self.setDomainType(domain_type)
    self.setDomain(domain)
    self.setEnergyGroups(energy_groups)
    self._tallies = tallies
    self._xs = xs
    self._offset = offset
    self._subdomain_offsets = subdomain_offsets



  #FIXME: Modularize!!!
  def exportSubdomainResults(self, subdomain=None, filename='multigroupxs',
                             directory='multigroupxs', format='hdf5',
                             append=True, uncertainties=False):

    if not subdomain is None and not subdomain in self._subdomain_offsets.keys():
      msg = 'Unable to export cross-section for domain {0} since it is ' \
            'not one of the domain offsets'.format(subdomain)
      raise ValueError(msg)

    elif subdomain is None and self._domain_type == 'distribcell':
      msg = 'Unable to export cross-section for a distribcell ' \
            'since no subdomain was provided'
      raise ValueError(msg)

    elif not is_string(filename):
      msg = 'Unable to export cross-section to filename={0} ' \
            'since it is not a string'.format(filename)
      raise ValueError(msg)

    elif not is_string(directory):
      msg = 'Unable to export cross-section to directory={0} ' \
            'since it is not a string'.format(directory)
      raise ValueError(msg)

    elif not format in ['hdf5', 'pkl', 'latex']:
      msg = 'Unable to export cross-section to format {0} ' \
            'since it is not supported'.format(format)
      raise ValueError(msg)

    elif not isinstance(append, (bool, np.bool)):
      msg = 'Unable to export cross-section since the append ' \
            'parameters is not True/False'.format(append)
      raise ValueError(msg)

    elif not isinstance(uncertainties, (bool, np.bool)):
      msg = 'Unable to export cross-section since the uncertainties ' \
            'parameters is not True/False'.format(uncertainties)
      raise ValueError(msg)

    elif self._xs is None:
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
      average = unumpy.nominal_values(self._xs)
      shape = average[subdomain_index,:].shape
      xs_group.require_dataset('average', dtype=np.float64, shape=shape,
                               data=average[subdomain_index,:])

      if uncertainties:
        std_dev = unumpy.std_devs(self._xs)
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
      average = unumpy.nominal_values(self._xs)
      xs_group['average'] = average[subdomain_index,:]

      if uncertainties:
        std_dev = unumpy.std_devs(self._xs)
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
      average = unumpy.nominal_values(self._xs)

      if uncertainties:
        std_dev = unumpy.std_devs(self._xs)

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


  def printPDF(self, subdomain=None, filename='multigroupxs', directory='.',
               uncertainties=False):

    if not subdomain is None and not subdomain in self._subdomain_offsets.keys():
      msg = 'Unable to print PDF for cross-section for domain {0} since it is ' \
            'not one of the domain offsets'.format(subdomain)
      raise ValueError(msg)

    elif subdomain is None and self._domain_type == 'distribcell':
      msg = 'Unable to print PDF for cross-section for a distribcell ' \
            'since no subdomain was provided'
      raise ValueError(msg)

    elif not is_string(filename):
      msg = 'Unable to print PDF for cross-section to filename={0} ' \
            'since it is not a string'.format(filename)
      raise ValueError(msg)

    elif not is_string(directory):
      msg = 'Unable to print PDF for cross-section to directory={0} ' \
            'since it is not a string'.format(directory)
      raise ValueError(msg)

    elif not isinstance(uncertainties, (bool, np.bool)):
      msg = 'Unable to print PDF for cross-section since the uncertainties ' \
            'parameters is not True/False'.format(uncertainties)
      raise ValueError(msg)

    elif self._xs is None:
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

    string += '{0: <16}\n'.format('\tCross-Sections [cm^-1]:')

    # Loop over energy groups ranges
    for group in range(1,self._num_groups+1):
      bounds = self._energy_groups.getGroupBounds(group)
      string += '{0: <12}Group {1} [{2: <10} - ' \
                '{3: <10}MeV]:\t'.format('', group, bounds[0], bounds[1])
      string += '{:.2e}'.format(self._xs[subdomain_index, group-1])
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
    newshape = (num_subdomains, self._num_groups)

    flux_mean = np.reshape(flux_mean, newshape)
    total_mean = np.reshape(total_mean, newshape)
    flux_std_dev = np.reshape(flux_std_dev, newshape)
    total_std_dev = np.reshape(total_std_dev, newshape)

    flux = unumpy.uarray(flux_mean, flux_std_dev)
    total = unumpy.uarray(total_mean, total_std_dev)

    # Set any zero fluxes to a negative value
    flux[flux == 0.] = unumpy.uarray([-1.], [0.])

    # Compute the xs with uncertainty propagation
    self._xs = unumpy_uncorr_divide(total, flux)

    # For any region without flux, convert xs to zero
    self._xs[flux == 0.] = unumpy.uarray([0.], [0.])

    # Correct -0.0 to +0.0
    self._xs += 0.

    # Reverse array so that it is ordered intuitively from high to low energy
    # Codes such as OpenMOC expect cross-sections ordered in this way
    self._xs = self._xs[:,::-1]



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
    newshape = (num_subdomains, self._num_groups)

    flux_mean = np.reshape(flux_mean, newshape)
    total_mean = np.reshape(total_mean, newshape)
    scatter1_mean = np.reshape(scatter1_mean, newshape)
    flux_std_dev = np.reshape(flux_std_dev, newshape)
    total_std_dev = np.reshape(total_std_dev, newshape)
    scatter1_std_dev = np.reshape(scatter1_std_dev, newshape)

    flux = unumpy.uarray(flux_mean, flux_std_dev)
    total = unumpy.uarray(total_mean, total_std_dev)
    scatter1 = unumpy.uarray(scatter1_mean, scatter1_std_dev)

    # Set any zero fluxes to a negative value
    flux[flux == 0.] = unumpy.uarray([-1.], [0.])

    # Compute the xs with uncertainty propagation
    self._xs = unumpy_uncorr_divide(unumpy_uncorr_sub(total, scatter1), flux)

    # For any region without flux, convert xs to zero
    self._xs[flux == 0.] = unumpy.uarray([0.], [0.])

    # Correct -0.0 to +0.0
    self._xs += 0.

    # Reverse array so that it is ordered intuitively from high to low energy
    # Codes such as OpenMOC expect cross-sections ordered in this way
    self._xs = self._xs[:,::-1]



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
    newshape = (num_subdomains, self._num_groups)

    flux_mean = np.reshape(flux_mean, newshape)
    absorption_mean = np.reshape(absorption_mean, newshape)
    flux_std_dev = np.reshape(flux_std_dev, newshape)
    absorption_std_dev = np.reshape(absorption_std_dev, newshape)

    flux = unumpy.uarray(flux_mean, flux_std_dev)
    absorption = unumpy.uarray(absorption_mean, absorption_std_dev)

    # Set any zero fluxes to a negative value
    flux[flux == 0.] = unumpy.uarray([-1.], [0.])

    # Compute the xs with uncertainty propagation
    self._xs = unumpy_uncorr_divide(absorption, flux)

    # For any region without flux, convert xs to zero
    self._xs[flux == 0.] = unumpy.uarray([0.], [0.])

    # Correct -0.0 to +0.0
    self._xs += 0.

    # Reverse array so that it is ordered intuitively from high to low energy
    # Codes such as OpenMOC expect cross-sections ordered in this way
    self._xs = self._xs[:,::-1]



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
    newshape = (num_subdomains, self._num_groups)

    flux_mean = np.reshape(flux_mean, newshape)
    fission_mean = np.reshape(fission_mean, newshape)
    flux_std_dev = np.reshape(flux_std_dev, newshape)
    fission_std_dev = np.reshape(fission_std_dev, newshape)

    flux = unumpy.uarray(flux_mean, flux_std_dev)
    fission = unumpy.uarray(fission_mean, fission_std_dev)

    # Set any zero fluxes to a negative value
    flux[flux == 0.] = unumpy.uarray([-1.], [0.])

    # Compute the xs with uncertainty propagation
    self._xs = unumpy_uncorr_divide(fission, flux)

    # For any region without flux, convert xs to zero
    self._xs[flux == 0.] = unumpy.uarray([0.], [0.])

    # Correct -0.0 to +0.0
    self._xs += 0.

    # Reverse array so that it is ordered intuitively from high to low energy
    # Codes such as OpenMOC expect cross-sections ordered in this way
    self._xs = self._xs[:,::-1]



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
    newshape = (num_subdomains, self._num_groups)

    flux_mean = np.reshape(flux_mean, newshape)
    nu_fission_mean = np.reshape(nu_fission_mean, newshape)
    flux_std_dev = np.reshape(flux_std_dev, newshape)
    nu_fission_std_dev = np.reshape(nu_fission_std_dev, newshape)

    flux = unumpy.uarray(flux_mean, flux_std_dev)
    nu_fission = unumpy.uarray(nu_fission_mean, nu_fission_std_dev)

    # Set any zero fluxes to a negative value
    flux[flux == 0.] = unumpy.uarray([-1.], [0.])

    # Compute the xs with uncertainty propagation
    self._xs = unumpy_uncorr_divide(nu_fission, flux)

    # For any region without flux, convert xs to zero
    self._xs[flux == 0.] = unumpy.uarray([0.], [0.])

    # Correct -0.0 to +0.0
    self._xs += 0.

    # Reverse array so that it is ordered intuitively from high to low energy
    # Codes such as OpenMOC expect cross-sections ordered in this way
    self._xs = self._xs[:,::-1]



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
    newshape = (num_subdomains, self._num_groups)

    flux_mean = np.reshape(flux_mean, newshape)
    scatter_mean = np.reshape(scatter_mean, newshape)
    flux_std_dev = np.reshape(flux_std_dev, newshape)
    scatter_std_dev = np.reshape(scatter_std_dev, newshape)

    flux = unumpy.uarray(flux_mean, flux_std_dev)
    scatter = unumpy.uarray(scatter_mean, scatter_std_dev)

    # Set any zero fluxes to a negative value
    flux[flux == 0.] = unumpy.uarray([-1.], [0.])

    # Compute the xs with uncertainty propagation
    self._xs = unumpy_uncorr_divide(scatter, flux)

    # For any region without flux, convert xs to zero
    self._xs[flux == 0.] = unumpy.uarray([0.], [0.])

    # Correct -0.0 to +0.0
    self._xs += 0.

    # Reverse array so that it is ordered intuitively from high to low energy
    # Codes such as OpenMOC expect cross-sections ordered in this way
    self._xs = self._xs[:,::-1]



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
    newshape = (num_subdomains, self._num_groups)

    flux_mean = np.reshape(flux_mean, newshape)
    nu_scatter_mean = np.reshape(nu_scatter_mean, newshape)
    flux_std_dev = np.reshape(flux_std_dev, newshape)
    nu_scatter_std_dev = np.reshape(nu_scatter_std_dev, newshape)

    flux = unumpy.uarray(flux_mean, flux_std_dev)
    nu_scatter = unumpy.uarray(nu_scatter_mean, nu_scatter_std_dev)

    # Set any zero fluxes to a negative value
    flux[flux == 0.] = unumpy.uarray([-1.], [0.])

    # Compute the xs with uncertainty propagation
    self._xs = unumpy_uncorr_divide(nu_scatter, flux)

    # For any region without flux, convert xs to zero
    self._xs[flux == 0.] = unumpy.uarray([0.], [0.])

    # Reverse array so that it is ordered intuitively from high to low energy
    # Correct -0.0 to +0.0
    self._xs += 0.

    # Codes such as OpenMOC expect cross-sections ordered in this way
    self._xs = self._xs[:,::-1]



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
    newshape = (num_subdomains, self._num_groups, self._num_groups)

    flux_mean = np.reshape(flux_mean, newshape[0:2])
    nu_scatter_mean = np.reshape(nu_scatter_mean, newshape)
    flux_std_dev = np.reshape(flux_std_dev, newshape[0:2])
    nu_scatter_std_dev = np.reshape(nu_scatter_std_dev, newshape)

    flux = unumpy.uarray(flux_mean, flux_std_dev)
    nu_scatter = unumpy.uarray(nu_scatter_mean, nu_scatter_std_dev)

    # Set any zero fluxes to a negative value
    flux[flux == 0.] = unumpy.uarray([-1.], [0.])

    # Tile the flux into a 3D array corresponding to the nu-scatter array
    flux = np.reshape(flux, (num_subdomains, self._num_groups, 1))
    flux = np.repeat(flux, self._num_groups, axis=2)

    # Compute the xs with uncertainty propagation
    self._xs = unumpy_uncorr_divide(nu_scatter, flux)

    # For any region without flux, convert xs to zero
    self._xs[flux == 0.] = unumpy.uarray([0.], [0.])

    # Correct -0.0 to +0.0
    self._xs += 0.

    # Reverse array so that it is ordered intuitively from high to low energy
    # Codes such as OpenMOC expect cross-sections ordered in this way
    self._xs = self._xs[:,::-1,::-1]


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
        string += '{:.2e}'.format(self._xs[subdomain_index, in_group-1, out_group-1])
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

    xs = self._xs[subdomain_index, in_group_index, out_group_index]

    if metric == 'mean':
      return unumpy.nominal_values(xs)
    else:
      return unumpy.std_devs(xs)



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
    newshape = (num_subdomains, self._num_groups)

    nu_fission_in_mean = np.reshape(nu_fission_in_mean, newshape)
    nu_fission_out_mean = np.reshape(nu_fission_out_mean, newshape)
    nu_fission_in_std_dev = np.reshape(nu_fission_in_std_dev, newshape)
    nu_fission_out_std_dev = np.reshape(nu_fission_out_std_dev, newshape)

    nu_fission_in = unumpy.uarray(nu_fission_in_mean, nu_fission_in_std_dev)
    nu_fission_out = unumpy.uarray(nu_fission_out_mean, nu_fission_out_std_dev)

    # Compute the xs with uncertainty propagation
    norm = nu_fission_in.sum()

    # If the xs are from a non-fissionable domain
    if norm.nominal_value == 0.:
      self._xs = unumpy.uarray(np.zeros(nu_fission_in.shape),
                               np.zeros(nu_fission_in.shape))

    # If the xs are from a fissionable domain
    else:

      self._xs = nu_fission_out / nu_fission_in.sum(1)[:, np.newaxis]

      # Normalize chi to 1.0
      self._xs = self._xs / self._xs.sum(1)[:, np.newaxis]

      # Correct -0.0 to +0.0
      self._xs += 0.

      # Reverse array so that it is ordered intuitively from high to low energy
      # Codes such as OpenMOC expect cross-sections ordered in this way
      self._xs = self._xs[:,::-1]