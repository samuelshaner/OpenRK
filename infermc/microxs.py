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


  #FIXME
  def computeXS(self):

    # Define a routine to convert 0 to 1
    nonzero = lambda val: 1 if not val else val

    num_nuclides = self.getNumNuclides()
    num_subdomains = self.getNumSubdomains()

    # Reshape xs array to be indexed by (domain, nuclide, group)
    new_shape = (nonzero(num_subdomains),
                 nonzero(num_nuclides),
                 nonzero(self._num_groups))

    self._xs = np.reshape(self._xs, new_shape)


  #FIXME
#  def getXS(self, group, nuclide, subdomain=None, metric='mean'):

  # FIXME
#  def dumpToFile(self, filename='multigroupxs', directory='multigroupxs'):

  # FIXME
#  def restoreFromFile(self, filename, directory='.'):

  #FIXME:
#  def exportSubdomainResults(self, subdomain=None, filename='multigroupxs',
#                             directory='multigroupxs', format='hdf5',
#                             append=True, uncertainties=False):

  # FIXME
#  def printPDF(self, subdomain=None, filename='multigroupxs', directory='.',
#               uncertainties=False):

  # FIXME
#  def printXS(self, subdomain=None):


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

    # Create flux and total reaction rate Tallies for all Nuclides
    super(MicroTotalXS, self).createTallies()
    self.addNuclidesToTallies()


class MicroTransportXS(MicroXS, infermc.TransportXS):

  def createTallies(self):

    # Create flux, total and scatter-p1 reaction rate Tallies for all Nuclides
    super(MicroTransportXS, self).createTallies()
    self.addNuclidesToTallies()


class MicroAbsorptionXS(MicroXS, infermc.AbsorptionXS):

  def createTallies(self):

    # Create flux and absorption reaction rate Tallies for all Nuclides
    super(MicroAbsorptionXS, self).createTallies()
    self.addNuclidesToTallies()


class MicroFissionXS(MicroXS, infermc.FissionXS):

  def createTallies(self):

    # Create flux and fission reaction rate Tallies for all Nuclides
    super(MicroFissionXS, self).createTallies()
    self.addNuclidesToTallies()


class MicroNuFissionXS(MicroXS, infermc.NuFissionXS):

  def createTallies(self):

    # Create flux and nu-fission reaction rate Tallies for all Nuclides
    super(MicroNuFissionXS, self).createTallies()
    self.addNuclidesToTallies()


class MicroScatterXS(MicroXS, infermc.ScatterXS):

  def createTallies(self):

    # Create flux and scatter reaction rate Tallies for all Nuclides
    super(MicroScatterXS, self).createTallies()
    self.addNuclidesToTallies()


class MicroNuScatterXS(MicroXS, infermc.NuScatterXS):

  def createTallies(self):

    # Create flux and nu-scatter reaction rate Tallies for all Nuclides
    super(MicroNuScatterXS, self).createTallies()
    self.addNuclidesToTallies()


class MicroScatterMatrixXS(MicroXS, infermc.ScatterMatrixXS):

  def createTallies(self):

    # Create flux and scatter reaction rate Tallies for all Nuclides
    super(MicroScatterMatrixXS, self).createTallies()
    self.addNuclidesToTallies()


class MicroDiffusionCoeff(MicroXS, infermc.DiffusionCoeff):

  def createTallies(self):

    # Create flux and reaction rate Tallies for all Nuclides
    super(MicroDiffusionCoeff, self).createTallies()
    self.addNuclidesToTallies()


class MicroChi(MicroXS, infermc.Chi):

  def createTallies(self):

    # Create flux and nu-fission reaction rate Tallies for all Nuclides
    super(MicroChi, self).createTallies()
    self.addNuclidesToTallies()