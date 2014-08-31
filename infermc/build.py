from openmc import Nuclide, Geometry
from openmc.tallies import TalliesFile
import infermc
from infermc.microxs import *

# Type-checking support
from typecheck import accepts, Or, Exact, Self


class XSTallyFactory(object):

  def __init__(self, geometry=None):

    self._all_xs = list()
    self._tallies_file = TalliesFile()
    self._geometry = None

    if not geometry is None:
      self.geometry = geometry


  @property
  def geometry(self):
    return self._geometry


  @geometry.setter
  @accepts(Self(), Geometry)
  def geometry(self, geometry):
    self._geometry = geometry


  @accepts(Self(), str, infermc.EnergyGroups, infermc.domains_check, infermc.domain_types_check)
  def createXS(self, xs_type, energy_groups, domain, domain_type='distribcell'):

    if xs_type == 'total':
      xs = infermc.TotalXS(domain, domain_type, energy_groups)
    elif xs_type == 'transport':
      xs = infermc.TransportXS(domain, domain_type, energy_groups)
    elif xs_type == 'absorption':
      xs = infermc.AbsorptionXS(domain, domain_type, energy_groups)
    elif xs_type == 'fission':
      xs = infermc.FissionXS(domain, domain_type, energy_groups)
    elif xs_type == 'nu-fission':
      xs = infermc.NuFissionXS(domain, domain_type, energy_groups)
    elif xs_type == 'scatter':
      xs = infermc.ScatterXS(domain, domain_type, energy_groups)
    elif xs_type == 'nu-scatter':
      xs = infermc.NuScatterXS(domain, domain_type, energy_groups)
    elif xs_type == 'scatter matrix':
      xs = infermc.ScatterMatrixXS(domain, domain_type, energy_groups)
    elif xs_type == 'diffusion':
      xs = infermc.DiffusionCoeff(domain, domain_type, energy_groups)
    elif xs_type == 'chi':
      xs = infermc.Chi(domain, domain_type, energy_groups)
    else:
      msg = 'The XSTallyBuilder is unable to create cross-section type ' \
            '{0} which is not one of the supported types'.format(xs_type)
      raise ValueError(msg)

    xs.createTallies()

    self._all_xs.append(xs)


  @accepts(Self(), str, infermc.EnergyGroups, infermc.domains_check,
           infermc.domain_types_check, Or(Exact(None), [Nuclide]))
  def createMicroXS(self, xs_type, energy_groups, domain,
                    domain_type='distribcell', nuclides=None):

    # Add all of the Nuclides in the domain to the MicroXS
    if nuclides is None:
      nuclides = domain.getAllNuclides()
      nuclides = nuclides.values()
      nuclides.append(openmc.Nuclide('total'))

    if xs_type == 'total':
      xs = MicroTotalXS(domain, domain_type, energy_groups, nuclides)
    elif xs_type == 'transport':
      xs = MicroTransportXS(domain, domain_type, energy_groups, nuclides)
    elif xs_type == 'absorption':
      xs = MicroAbsorptionXS(domain, domain_type, energy_groups, nuclides)
    elif xs_type == 'fission':
      xs = MicroFissionXS(domain, domain_type, energy_groups, nuclides)
    elif xs_type == 'nu-fission':
      xs = MicroNuFissionXS(domain, domain_type, energy_groups, nuclides)
    elif xs_type == 'scatter':
      xs = MicroScatterXS(domain, domain_type, energy_groups, nuclides)
    elif xs_type == 'nu-scatter':
      xs = MicroNuScatterXS(domain, domain_type, energy_groups, nuclides)
    elif xs_type == 'scatter matrix':
      xs = MicroScatterMatrixXS(domain, domain_type, energy_groups, nuclides)
    elif xs_type == 'diffusion':
      xs = MicroDiffusionCoeff(domain, domain_type, energy_groups, nuclides)
    elif xs_type == 'chi':
      xs = MicroChi(domain, domain_type, energy_groups, nuclides)
    else:
      msg = 'The XSTallyBuilder is unable to create micro cross-section ' \
            'type {0} which is not one of the supported types'.format(xs_type)
      raise ValueError(msg)

    xs.createTallies()

    self._all_xs.append(xs)


  @accepts(Self(), infermc.EnergyGroups, infermc.domain_types_check)
  def createAllXS(self, energy_groups, domain_type='distribcell'):

    if domain_type == 'distribcell' or domain_type == 'cell':
      domains = self._geometry.getAllMaterialCells()
    elif domain_type == 'universe':
      domains = self._geometry.getAllMaterialUniverses()
    elif domain_type == 'material':
      domains = self._geometry.getAllMaterials()

    for domain in domains:
      for xs_type in infermc.xs_types:
        self.createXS(xs_type, energy_groups, domain, domain_type)


  @accepts(Self(), infermc.EnergyGroups, infermc.domain_types_check, Or(Exact(None), [Nuclide]))
  def createAllMicroXS(self, energy_groups, domain_type='distribcell',
                       nuclides=None):

    if domain_type == 'distribcell' or domain_type == 'cell':
      domains = self._geometry.getAllMaterialCells()
    elif domain_type == 'universe':
      domains = self._geometry.getAllMaterialUniverses()
    elif domain_type == 'material':
      domains = self._geometry.getAllMaterials()

    for domain in domains:
      for xs_type in infermc.xs_types:
        self.createMicroXS(xs_type, energy_groups, domain, domain_type, nuclides)


  def createTalliesFile(self):

    tallies = set()

    for xs in self._all_xs:
      tallies = tallies.union(set(xs._tallies.values()))

    for tally in tallies:
      self._tallies_file.addTally(tally)

    self._tallies_file.exportToXML()