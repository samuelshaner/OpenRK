from infermc.multigroupxs import *
import openmc


class MicroXS(MultiGroupXS):

  def __init__(self, domain=None, domain_type=None, \
               energy_groups=None, nuclides=None):

    super(MicroXS, self).__init__(domain, domain_type, energy_groups)

    # Initialize an empty list for OpenMC Nuclide objects
    self._nuclides = list()

    if not nuclides is None:
      self.addNuclides(nuclides)


  def addNuclide(self, nuclide):

    if not isinstance(nuclide, openmc.Nuclide):
      msg = 'Unable to add Nuclide {0} to MicroXS since ' \
            'it is not an OpenMC Nuclide'.format(nuclide)
      raise ValueError(msg)

    self._nuclides.append(nuclide)


  def addNuclides(self, nuclides):

    if not isinstance(nuclides, (np.ndarray, list, tuple)):
      msg = 'Unable to add Nuclides {0} to MicroXS since it is not a ' \
            'NumPy array or Python tuple list of Nuclides'.format(nuclides)
      raise ValueError(msg)

    for nuclide in nuclides:
      self.addNuclide(nuclide)


  def extractNuclidesFromDomain(self):

    if self._domain is None:

    nuclides = self._domain.getAllNuclides()

    if len(self._tallies) == 0:
      msg = 'Unable to add Nuclides to MicroXS since its ' \
            'Tallies have not yet been created'
      raise ValueError(msg)



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


class MicroTotalXS(MicroXS, TotalXS):

  def createTallies(self):

    # Create the flux and total reaction rate Tallies
    super(MicroTotalXS, self).createTallies()

    # Add all nuclides to the Tallies
    self.addNuclidesToTallies()


class MicroTransportXS(MicroXS, TransportXS):

  def createTallies(self):

    # Create the flux and total and scatter-p1 reaction rate Tallies
    super(MicroTransportXS, self).createTallies()

    # Add all nuclides to the Tallies
    self.addNuclidesToTallies()


class MicroAbsorptionXS(MicroXS, AbsorptionXS):

  def createTallies(self):

    # Create the flux and absorption reaction rate Tallies
    super(MicroAbsorptionXS, self).createTallies()

    # Add all nuclides to the Tallies
    self.addNuclidesToTallies()


class MicroFissionXS(MicroXS, FissionXS):

  def createTallies(self):

    # Create the flux and fission reaction rate Tallies
    super(MicroFissionXS, self).createTallies()

    # Add all nuclides to the Tallies
    self.addNuclidesToTallies()


class MicroNuFissionXS(MicroXS, NuFissionXS):

  def createTallies(self):

    # Create the flux and nu-fission reaction rate Tallies
    super(MicroNuFissionXS, self).createTallies()

    # Add all nuclides to the Tallies
    self.addNuclidesToTallies()


class MicroScatterXS(MicroXS, ScatterXS):

  def createTallies(self):

    # Create the flux and scatter reaction rate Tallies
    super(MicroScatterXS, self).createTallies()

    # Add all nuclides to the Tallies
    self.addNuclidesToTallies()


class MicroNuScatterXS(MicroXS, NuScatterXS):

  def createTallies(self):

    # Create the flux and nu-scatter reaction rate Tallies
    super(MicroNuScatterXS, self).createTallies()

    # Add all nuclides to the Tallies
    self.addNuclidesToTallies()


class MicroScatterMatrixXS(MicroXS, ScatterMatrixXS):

  def createTallies(self):

    # Create the flux and scatter reaction rate Tallies
    super(MicroScatterMatrixXS, self).createTallies()

    # Add all nuclides to the Tallies
    self.addNuclidesToTallies()


class MicroDiffusionCoeff(MicroXS, DiffusionCoeff):

  def createTallies(self):

    # Create the flux and  reaction rate Tallies
    super(MicroDiffusionCoeff, self).createTallies()

    # Add all nuclides to the Tallies
    self.addNuclidesToTallies()


class MicroChi(MicroXS, Chi):

  def createTallies(self):

    # Create the flux and nu-fission reaction rate Tallies
    super(MicroChi, self).createTallies()

    # Add all nuclides to the Tallies
    self.addNuclidesToTallies()
