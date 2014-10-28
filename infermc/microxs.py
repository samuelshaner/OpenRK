import numpy as np
import os, copy
import openmc
import infermc

# LaTeX Greek symbols for each cross-section type
greek = dict()
greek['total'] = '$\\sigma_{t}$'
greek['transport'] = '$\\sigma_{tr}$'
greek['absorption'] = '$\\sigma_{a}$'
greek['scatter'] = '$\\sigma_{s}$'
greek['nu-scatter'] = '$\\nu\\sigma_{s}$'
greek['scatter matrix'] = '$\\sigma_{s}$'
greek['fission'] = '$\\sigma_{f}$'
greek['nu-fission'] = '$\\nu\\sigma_{f}$'
greek['chi'] = '$\\chi$'
greek['diffusion'] = '$D$'


class MicroXS(infermc.MultiGroupXS):

  def __init__(self, domain=None, domain_type=None,
               energy_groups=None, nuclides=None, densities=None):

    super(MicroXS, self).__init__(domain, domain_type, energy_groups)

    # Initialize empty lists of Nuclides and number densities (at/b-cm)
    self._nuclides = list()
    self._densities = np.zeros(0)
    self._num_nuclides = 0

    if not nuclides is None:
      self.addNuclides(nuclides, densities)


  def __deepcopy__(self, memo):

    existing = memo.get(id(self))

    # If this is the first time we have tried to copy this object, create a copy
    if existing is None:

      clone = type(self).__new__(type(self))
      clone._xs_type = self._xs_type
      clone._domain = self._domain
      clone._domain_type = self._domain_type
      clone._energy_groups = copy.deepcopy(self._energy_groups, memo)
      clone._num_groups = self._num_groups
      clone._xs = copy.deepcopy(self._xs, memo)
      clone._colors = copy.deepcopy(self._colors, memo)
      clone._subdomain_offsets = copy.deepcopy(self._subdomain_offsets, memo)
      clone._offset = copy.deepcopy(self._offset, memo)
      clone._densities = np.zeros(0)
      clone._num_nuclides = 0

      clone._tallies = dict()
      for tally_type, tally in self._tallies.items():
        clone._tallies[tally_type] = copy.deepcopy(tally, memo)

      clone._nuclides = list()
      for i, nuclide in enumerate(self._nuclides):
        clone.addNuclide(copy.deepcopy(nuclide), self._densities[i])

      memo[id(self)] = clone

      return clone

    # If this object has been copied before, return the first copy made
    else:
      return existing


  def addNuclide(self, nuclide, density=None):

    self._nuclides.append(nuclide)

    if not density is None:
      self._densities = np.append(self._densities, density)

    # FIXME: Add a dummy density if one was not specified
    else:
      self._densities = np.append(self._densities, -1)

    self._num_nuclides += 1


  def addNuclides(self, nuclides, densities=None):

    for i, nuclide in enumerate(nuclides):
      if densities is None:
        self.addNuclide(nuclide)
      else:
        self.addNuclide(nuclide, densities[i])


  def addNuclidesToTallies(self):

    for tally_id, tally in self._tallies.items():
      if 'flux' in tally._scores:
        continue
      for nuclide in self._nuclides:
        tally.addNuclide(nuclide)


  def containsNuclide(self, nuclide):
    return (nuclide in self._nuclides)


  def getNuclideIndices(self, nuclides='all'):

    if nuclides == 'all':
      indices = np.arange(self._num_nuclides)

    else:
      indices = np.zeros(len(nuclides), dtype=np.int64)

      for i, nuclide in enumerate(nuclides):
        try:
          indices[i] = self._nuclides.index(nuclide)
        except ValueError:
          msg = 'Unable to get index for Nuclide {0} since it is not one ' \
                'of the Nuclides in the cross-section'.format(nuclide)
          raise ValueError(msg)

    return indices


  def computeXS(self, corr=True):

    super(MicroXS, self).computeXS(corr)

    # Divide out the densities to convert xs to barns
    if self._xs_type != 'chi':
      self._xs = infermc.error_prop.arithmetic.divide_by_scalar(self._xs, self._densities, corr)


  def getXS(self, groups='all', nuclides='all', subdomains='all', metric='mean'):
    xs = super(MicroXS, self).getXS(groups, subdomains, metric)
    nuclides = self.getNuclideIndices(nuclides)
    return xs[..., nuclides]


  def printXS(self, nuclides='all', subdomains='all'):

    string = 'Micro XS\n'
    string += '{0: <16}{1}{2}\n'.format('\tType', '=\t', self._xs_type)
    string += '{0: <16}{1}{2}\n'.format('\tDomain Type', '=\t', self._domain_type)
    string += '{0: <16}{1}{2}\n'.format('\tDomain ID', '=\t', self._domain._id)

    if subdomains == 'all':
      subdomains = self._subdomain_offsets.keys()

    if nuclides == 'all':
      nuclides = self._nuclides

    # Loop over all subdomains
    for subdomain in subdomains:

      if self._domain_type == 'distribcell':
        string += '{0: <16}{1}{2}\n'.format('\tSubDomain', '=\t', subdomain)

      # Loop over all Nuclides
      for nuclide in nuclides:
        string += '{0: <16}{1}{2}\n'.format('\tNuclide', '=\t', nuclide._name)
        string += '{0: <16}\n'.format('\tCross-Sections [barns]:')

        # Loop over energy group ranges
        for group in range(1,self._num_groups+1):
          bounds = self._energy_groups.getGroupBounds(group)
          string += '{0: <12}Group {1} [{2: <10} - ' \
                    '{3: <10}MeV]:\t'.format('', group, bounds[0], bounds[1])
          average = self.getXS([group], [nuclide], [subdomain], 'mean')
          std_dev = self.getXS([group], [nuclide], [subdomain], 'std_dev')

          std_dev[average == 0.] = 0.
          average[average == 0.] = 1.
          rel_err = (std_dev / average) * 100.

          string += '{:.2e}+/-{:1.2e}%'.format(average[0,0,0], std_dev[0,0,0])
          string += '\n'

        string += '\n'
      string += '\n'

    print(string)


  def dumpToFile(self, filename='multigroupxs', directory='multigroupxs'):

    # Export all data to the file except for the Nuclides
    super(MicroXS, self).dumpToFile(filename, directory)

    import pickle

    filename = directory + '/' + filename + '.pkl'
    filename = filename.replace(' ', '-')

    # Load pickle file, append the Nuclides, and save it to the file
    xs_results = pickle.load(open(filename, 'rb'))
    xs_results['nuclides'] = self._nuclides
    xs_results['densities'] = self._densities
    pickle.dump(xs_results, open(filename, 'wb'))


  def restoreFromFile(self, filename='multigroupxs', directory='multigroupxs'):

    # Import all data from the file except for the Nuclides
    super(MicroXS, self).restoreFromFile(filename, directory)

    import pickle

    filename = directory + '/' + filename + '.pkl'
    filename = filename.replace(' ', '-')

    # Load the pickle file into a dictionary
    xs_results = pickle.load(open(filename, 'rb'))

    # Extract the Nuclides and store them to the class attribute
    nuclides = xs_results['nuclides']
    densities = xs_results['densities']
    num_nuclides = len(nuclides)
    self.addNuclides(nuclides, densities)


  def exportResults(self, nuclides='all', subdomains='all',
                    filename='multigroupxs', directory='multigroupxs',
                    format='hdf5', append=True):

    # Make directory if it does not exist
    if not os.path.exists(directory):
      os.makedirs(directory)

    offsets = self.getSubDomainOffsets(subdomains)
    subdomains = self.getSubDomains(offsets)
    nuclides = self.getNuclideIndices(nuclides)

    average = self._xs[infermc.metrics['mean'], offsets, ...]
    average = average[..., nuclides]
    std_dev = self._xs[infermc.metrics['std_dev'], offsets, ...]
    std_dev = std_dev[..., nuclides]

    std_dev[average == 0.] = 0.
    average[average == 0.] = 1.
    rel_err = (std_dev / average) * 100.

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

        # Create an HDF5 group for the subdomain
        if self._domain_type == 'distribcell':
          group_name = '{0}'.format(subdomain)
          subdomain_group = domain_group.require_group(group_name)
        else:
          subdomain_group = domain_group

        # Loop over all Nuclides
        for nuclide in nuclides:

          # Create an HDF5 group for the Nuclide and xs type
          group_name = self._nuclides[nuclide]._name
          nuclide_group = subdomain_group.require_group(group_name)
          xs_group = nuclide_group.require_group(self._xs_type)

          # Add MultiGroupXS results data to the HDF5 group
          xs_group.require_dataset('average', dtype=np.float64,
                                   shape=average[i, ..., nuclide].shape,
                                   data=average[i, ..., nuclide])
          xs_group.require_dataset('rel. err.', dtype=np.float64,
                                   shape=rel_err[i, ..., nuclide].shape,
                                   data=rel_err[i, ..., nuclide])

      # Close the results HDF5 file
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
        domain_group = dict()
      else:
        domain_group = xs_results[group_name]

      # Loop over all subdomains
      for i, subdomain in enumerate(subdomains):

        # Create an HDF5 group for the subdomain
        if self._domain_type == 'distribcell':
          group_name = '{0}'.format(subdomain)

          if not domain_group.has_key(group_name):
            subdomain_group = domain_group[group_name] = dict()
          else:
            subdomain_group = domain_group[group_name]

        else:
          subdomain_group = domain_group

        # Loop over all Nuclides
        for nuclide in nuclides:

          # Create an HDF5 group for the Nuclide and xs type
          group_name = self._nuclides[nuclide]._name

          if not subdomain_group.has_key(group_name):
            nuclide_group = subdomain_group[group_name] = dict()
          else:
            nuclide_group = subdomain_group[group_name]

          xs_group = nuclide_group[self._xs_type] = dict()

          # Add MultiGroupXS results data to the dictionary
          xs_group['average'] = average[i, ..., nuclide]
          xs_group['rel. err.'] = rel_err[i, ..., nuclide]

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

      # Loop over all subdomains and Nuclides
      for i, subdomain in enumerate(subdomains):
        for nuclide in nuclides:

          # Add MultiGroupXS results data to the table list
          table = list()

          if self._xs_type != 'scatter matrix':
            headers = list()
            headers.append('Group')
            headers.append('Average XS')
            headers.append('Rel. Err. (\%)')

            for group in range(self._num_groups):
              subtable = list()
              subtable.append(group+1)
              subtable.append(average[i, group, ..., nuclide])
              subtable.append(rel_err[i, group, ..., nuclide])
              table.append(subtable)

          # Scattering matrix
          else:
            headers = list()
            headers.append('Group In')
            headers.append('Group Out')
            headers.append('Average XS')
            headers.append('Rel. Err. (\%)')

            for in_group in range(self._num_groups):
              for out_group in range(self._num_groups):
                subtable = list()
                subtable.append(in_group+1)
                subtable.append(out_group+1)
                subtable.append(average[i, in_group, out_group, ..., nuclide])
                subtable.append(rel_err[i, in_group, out_group, ..., nuclide])
                table.append(subtable)

          if self._domain_type == 'distribcell':
            caption = '\\caption{{{0} {1}, (subdomain {2}) {3} {4} [barns]}}'.format(
              self._domain_type.capitalize(), self._domain._id, subdomain,
              self._nuclides[nuclide]._name, greek[self._xs_type])
          else:
            caption = '\\caption{{{0} {1} {2} {3} [barns]}}'.format(
              self._domain_type.capitalize(), subdomain,
              self._nuclides[nuclide]._name, greek[self._xs_type])

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


  def printPDF(self, nuclides='all', subdomains='all',
               filename='multigroupxs', directory='multigroupxs'):

    import subprocess

    # Make directory if it does not exist
    if not os.path.exists(directory):
      os.makedirs(directory)

    filename = filename.replace(' ', '-')

    # Generate LaTeX file
    self.exportResults(nuclides, subdomains, filename, '.', 'latex', False)

    # Compile LaTeX to PDF
    FNULL = open(os.devnull, 'w')
    subprocess.check_call('pdflatex {0}.tex'.format(filename),
                          shell=True, stdout=FNULL)

    # Move PDF to requested directory and cleanup temporary LaTeX files
    if directory != '.':
      os.system('mv {0}.pdf {1}'.format(filename, directory))

    os.system('rm {0}.tex {0}.aux {0}.log'.format(filename))


  def checkXS(self):

    xs = super(MicroXS, self).getXS()

    if self._xs_type == 'chi':
      return

    total_index = self.getNuclideIndices([openmc.Nuclide('total')])
    nuclide_indices = self.getNuclideIndices()
    nuclide_indices = nuclide_indices[nuclide_indices != total_index]

    # Get the macroscopic cross-sections
    total_xs = xs[..., total_index]
    all_xs = xs[..., nuclide_indices] * self._densities[:-1]
    all_xs = all_xs.sum(axis=-1)

    if not np.allclose(total_xs.ravel(), all_xs.ravel()):
      print('The nuclide micro xs {0} in {1} {2} is not equal to the total macro ' \
            'macro xs'.format(self._xs_type, self._domain_type, self._domain._id))


class MicroTotalXS(MicroXS, infermc.TotalXS):

  def createTallies(self):
    super(MicroTotalXS, self).createTallies()
    self.addNuclidesToTallies()


class MicroTransportXS(MicroXS, infermc.TransportXS):

  def createTallies(self):
    super(MicroTransportXS, self).createTallies()
    self.addNuclidesToTallies()


class MicroAbsorptionXS(MicroXS, infermc.AbsorptionXS):

  def createTallies(self):
    super(MicroAbsorptionXS, self).createTallies()
    self.addNuclidesToTallies()


class MicroCaptureXS(MicroXS, infermc.CaptureXS):

  def createTallies(self):
    super(MicroCaptureXS, self).createTallies()
    self.addNuclidesToTallies()


class MicroFissionXS(MicroXS, infermc.FissionXS):

  def createTallies(self):
    super(MicroFissionXS, self).createTallies()
    self.addNuclidesToTallies()


class MicroNuFissionXS(MicroXS, infermc.NuFissionXS):

  def createTallies(self):
    super(MicroNuFissionXS, self).createTallies()
    self.addNuclidesToTallies()


class MicroScatterXS(MicroXS, infermc.ScatterXS):

  def createTallies(self):
    super(MicroScatterXS, self).createTallies()
    self.addNuclidesToTallies()


class MicroNuScatterXS(MicroXS, infermc.NuScatterXS):

  def createTallies(self):
    super(MicroNuScatterXS, self).createTallies()
    self.addNuclidesToTallies()


class MicroScatterMatrixXS(MicroXS, infermc.ScatterMatrixXS):

  def createTallies(self):
    super(MicroScatterMatrixXS, self).createTallies()
    self.addNuclidesToTallies()


  def getXS(self, in_groups='all', out_groups='all', nuclides='all',
            subdomains='all', metric='mean'):

    xs = super(MicroXS, self).getXS(in_groups, out_groups, subdomains, metric)
    nuclides = self.getNuclideIndices(nuclides)
    return xs[..., nuclides]


  def printXS(self, nuclides='all', subdomains='all'):

    string = 'Micro XS\n'
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

    if nuclides == 'all':
      nuclides = self._nuclides

    for subdomain in subdomains:

      if self._domain_type == 'distribcell':
        string += '{0: <16}{1}{2}\n'.format('\tSubDomain', '=\t', subdomain)

      # Loop over energy groups ranges
      for nuclide in nuclides:
        string += '{0: <16}{1}{2}\n'.format('\tNuclide', '=\t', nuclide._name)
        string += '{0: <16}\n'.format('\tCross-Sections [barns]:')

        # Loop over energy groups ranges
        for in_group in range(1,self._num_groups+1):
          for out_group in range(1,self._num_groups+1):
            string += '{0: <12}Group {1} -> Group {2}:\t\t'.format('', in_group, out_group)
            average = self.getXS([in_group], [out_group], [nuclide], [subdomain], 'mean')
            std_dev = self.getXS([in_group], [out_group], [nuclide], [subdomain], 'std_dev')

            std_dev[average == 0.] = 0.
            average[average == 0.] = 1.
            rel_err = (std_dev / average) * 100.

            string += '{:.2e}+/-{:1.2e}%'.format(average[0,0,0], rel_err[0,0,0])
            string += '\n'

        string += '\n'
      string += '\n'

    print(string)


class MicroDiffusionCoeff(MicroXS, infermc.DiffusionCoeff):

  def createTallies(self):
    super(MicroDiffusionCoeff, self).createTallies()
    self.addNuclidesToTallies()


class MicroChi(MicroXS, infermc.Chi):

  def createTallies(self):
    super(MicroChi, self).createTallies()
    self.addNuclidesToTallies()
