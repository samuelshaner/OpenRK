import openmc
import infermc
import numpy as np

# Type-checking support
from typecheck import accepts, Or, Exact, Self


class MicroXS(infermc.MultiGroupXS):

  def __init__(self, domain=None, domain_type=None,
               energy_groups=None, nuclides=None):

    super(MicroXS, self).__init__(domain, domain_type, energy_groups)

    # Initialize an empty list for OpenMC Nuclide objects
    self._nuclides = list()

    if not nuclides is None:
      self.addNuclides(nuclides)


  @accepts(Self(), int, openmc.Nuclide, Or(Exact(None), int), str)
  def getXS(self, group, nuclide, subdomain=None, metric='mean'):

    # Get the cross-sections for all Nuclides
    xs = super(MicroXS, self).getXS(group, subdomain, metric)

    # Get index of the Nuclide in array and return the corresponding value
    nuclide_index = self._nuclides.index(nuclide)
    return xs[nuclide_index]


  @accepts(Self(), Or(str, openmc.Nuclide), Or(Exact(None), int))
  def printXS(self, nuclide='all', subdomain=None):

    if self._xs is None:
      msg = 'Unable to print cross-section {0} since it has not yet ' \
            'been computed'.format(self._xs_type)
      raise ValueError(msg)

    if isinstance(nuclide, openmc.Nuclide):
      nuclides = [nuclide]
    else:
      nuclides = self._nuclides

    string = 'Micro XS\n'
    string += '{0: <16}{1}{2}\n'.format('\tType', '=\t', self._xs_type)
    string += '{0: <16}{1}{2}\n'.format('\tDomain Type', '=\t', self._domain_type)
    string += '{0: <16}{1}{2}\n'.format('\tDomain ID', '=\t', self._domain._id)

    string += '{0: <16}\n'.format('\tCross-Sections [barns]:')

    # Loop over energy groups ranges
    for nuclide in nuclides:
      string += '{0: <16}{1}{2}\n'.format('\tNuclide', '=\t', nuclide._name)

      for group in range(1,self._num_groups+1):
        bounds = self._energy_groups.getGroupBounds(group)
        string += '{0: <12}Group {1} [{2: <10} - ' \
                  '{3: <10}MeV]:\t'.format('', group, bounds[0], bounds[1])
        average = self.getXS(group, nuclide, subdomain, 'mean')
        std_dev = self.getXS(group, nuclide, subdomain, 'std_dev')
        string += '{:.2e}+/-{:.2e}'.format(average, std_dev)
        string += '\n'

      string += '\n'

    print(string)


  @accepts(Self(), str, str)
  def dumpToFile(self, filename='multigroupxs', directory='multigroupxs'):

    # Export all data to the file except for the Nuclides
    super(MicroXS, self).dumpToFile(filename, directory)

    import pickle

    filename = directory + '/' + filename + '.pkl'
    filename = filename.replace(' ', '-')

    # Load pickle file, append the Nuclides, and save it to the file
    xs_results = pickle.load(open(filename, 'rb'))
    xs_results['nuclides'] = self._nuclides
    pickle.dump(xs_results, open(filename, 'wb'))


  @accepts(Self(), str, str)
  def restoreFromFile(self, filename, directory='multigroupxs'):

    # Import all data from the file except for the Nuclides
    super(MicroXS, self).restoreFromFile(filename, directory)

    import pickle

    filename = directory + '/' + filename + '.pkl'
    filename = filename.replace(' ', '-')

    # Load the pickle file into a dictionary
    xs_results = pickle.load(open(filename, 'rb'))

    # Extract the Nuclides and store them to the class attribute
    nuclides = xs_results['nuclides']
    self._nuclides = nuclides


  #FIXME:
#  def exportSubdomainResults(self, subdomain=None, filename='multigroupxs',
#                             directory='multigroupxs', format='hdf5',
#                             append=True, uncertainties=False):

  # FIXME
#  def printPDF(self, subdomain=None, filename='multigroupxs', directory='.',
#               uncertainties=False):


  @accepts(Self(), openmc.Nuclide)
  def addNuclide(self, nuclide):
    self._nuclides.append(nuclide)


  @accepts(Self(), Or(list, tuple, np.ndarray))
  def addNuclides(self, nuclides):
    for nuclide in nuclides:
      self.addNuclide(nuclide)


  def addNuclidesToTallies(self):

    if len(self._tallies) == 0:
      msg = 'Unable to add Nuclides to MicroXS since its ' \
            'Tallies have not yet been created'
      raise ValueError(msg)

    # Loop over all tallies
    for tally_id, tally in self._tallies.items():

      # If a flux Tally, do not add nuclides to it
      if 'flux' in tally._scores:
        continue

      # Add all nuclides to this Tally
      for nuclide in self._nuclides:
        tally.addNuclide(nuclide)


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


  @accepts(Self(), int, int, openmc.Nuclide, Or(Exact(None), int), str)
  def getXS(self, in_group, out_group, nuclide, subdomain=None, metric='mean'):

    # Get the cross-sections for all Nuclides
    xs = super(MicroXS, self).getXS(in_group, out_group, subdomain, metric)

    # Get index of the Nuclide in array and return the corresponding value
    nuclide_index = self._nuclides.index(nuclide)
    return xs[nuclide_index]


  @accepts(Self(), Or(str, openmc.Nuclide), Or(Exact(None), int))
  def printXS(self, nuclide='all', subdomain=None):

    if self._xs is None:
      msg = 'Unable to print cross-section {0} since it has not yet ' \
            'been computed'.format(self._xs_type)
      raise ValueError(msg)

    if isinstance(nuclide, openmc.Nuclide):
      nuclides = [nuclide]
    else:
      nuclides = self._nuclides

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

    string += '{0: <16}\n'.format('\tCross-Sections [barns]:')

    # Loop over energy groups ranges
    for nuclide in nuclides:
      string += '{0: <16}{1}{2}\n'.format('\tNuclide', '=\t', nuclide._name)

      # Loop over energy groups ranges
      for in_group in range(1,self._num_groups+1):
        for out_group in range(1,self._num_groups+1):
          string += '{0: <12}Group {1} -> Group {2}:\t\t'.format('', in_group, out_group)
          average = self.getXS(in_group, out_group, nuclide, subdomain, 'mean')
          std_dev = self.getXS(in_group, out_group, nuclide, subdomain, 'std_dev')
          string += '{:.2e}+/-{:.2e}'.format(average, std_dev)
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