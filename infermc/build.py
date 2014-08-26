from openmc import Material, Cell, Universe, Geometry
from openmc.tallies import TalliesFile
from infermc.multigroupxs import *
from infermc.microxs import *


class XSTallyFactory(object):

  def __init__(self, geometry=None):

    self._all_xs = list()
    self._tallies_file = TalliesFile()
    self._geometry = None

    if not geometry is None:
      self.setGeometry(geometry)


  def setGeometry(self, geometry):

    if not isinstance(geometry, Geometry):
      msg = 'Unable to set the Geometry for XSTallyBuilder to {0} since ' \
            'it is not an OpenMC Geometry object'.format(geometry)
      raise ValueError(msg)

    self._geometry = geometry


  def createXS(self, xs_type, energy_groups, domain,
               domain_type='distribcell'):

    if not is_string(xs_type):
      msg = 'The XSTallyBuilder is unable to create cross-section type {0} ' \
            'which is not a string value'.format(xs_type)
      raise ValueError(msg)

    elif not isinstance(energy_groups, EnergyGroups):
      msg = 'The XSTallyBuilder is unable to create cross-section ' \
            'type {0} with energy groups {1} which is not an EnergyGroups ' \
            'object'.format(xs_type, energy_groups)
      raise ValueError(msg)

    elif not isinstance(domain, (Material, Cell, Universe)):
      msg = 'The XSTallyBuilder is unable to create cross-section ' \
            'type {0} in {1} {2} since the domain is not a Material, ' \
            'Cell or Universe'.format(xs_type, domain_type, domain)
      raise ValueError(msg)

    elif not domain_type in domain_types:
      msg = 'The XSTallyBuilder is unable to create cross-section ' \
            'type {0} in {1} {2} since it is not a supported domain ' \
            'type'.format(xs_type, domain_type, domain._id)
      raise ValueError(msg)


    if xs_type == 'total':
      xs = TotalXS(domain, domain_type, energy_groups)
    elif xs_type == 'transport':
      xs = TransportXS(domain, domain_type, energy_groups)
    elif xs_type == 'absorption':
      xs = AbsorptionXS(domain, domain_type, energy_groups)
    elif xs_type == 'fission':
      xs = FissionXS(domain, domain_type, energy_groups)
    elif xs_type == 'nu-fission':
      xs = NuFissionXS(domain, domain_type, energy_groups)
    elif xs_type == 'scatter':
      xs = ScatterXS(domain, domain_type, energy_groups)
    elif xs_type == 'nu-scatter':
      xs = NuScatterXS(domain, domain_type, energy_groups)
    elif xs_type == 'scatter matrix':
      xs = ScatterMatrixXS(domain, domain_type, energy_groups)
    elif xs_type == 'diffusion':
      xs = DiffusionCoeff(domain, domain_type, energy_groups)
    elif xs_type == 'chi':
      xs = Chi(domain, domain_type, energy_groups)
    else:
      msg = 'The XSTallyBuilder is unable to create cross-section type ' \
            '{0} which is not one of the supported types'.format(xs_type)
      raise ValueError(msg)

    xs.createTallies()

    self._all_xs.append(xs)


  def createMicroXS(self, xs_type, energy_groups, domain,
               domain_type='distribcell', nuclides=None):

    if not is_string(xs_type):
      msg = 'The XSTallyBuilder is unable to create cross-section type {0} ' \
            'which is not a string value'.format(xs_type)
      raise ValueError(msg)

    elif not isinstance(energy_groups, EnergyGroups):
      msg = 'The XSTallyBuilder is unable to create cross-section ' \
            'type {0} with energy groups {1} which is not an EnergyGroups ' \
            'object'.format(xs_type, energy_groups)
      raise ValueError(msg)

    elif not isinstance(domain, (Material, Cell, Universe)):
      msg = 'The XSTallyBuilder is unable to create cross-section ' \
            'type {0} in {1} {2} since the domain is not a Material, ' \
            'Cell or Universe'.format(xs_type, domain_type, domain)
      raise ValueError(msg)

    elif not domain_type in domain_types:
      msg = 'The XSTallyBuilder is unable to create cross-section ' \
            'type {0} in {1} {2} since it is not a supported domain ' \
            'type'.format(xs_type, domain_type, domain._id)
      raise ValueError(msg)

    # Add all of the Nuclides in the domain to the MicroXS
    if nuclides is None:
      nuclides = domain.getAllNuclides()
      nuclides = nuclides.values()

    elif not isinstance(nuclides, (np.ndarray, list, tuple)):
      msg = 'The XSTallyBuilder is unable to create cross-section ' \
            'type {0} in {1} {2} for nuclides {3} since it is not ' \
            'a NumPy array or Python list/tuple of nuclides'.format(
            xs_type, domain_type, domain._id, nuclides)
      raise ValueError(msg)

    for nuclide in nuclides:
      if not isinstance(nuclide, openmc.Nuclide):
        msg = 'The XSTallyBuilder is unable to create micro cross-section ' \
              'type {0} in {1} {2} for nuclide {3} since it is not an OpenMC ' \
              'nuclide'.format(xs_type, domain_type, domain._id, nuclide)
        raise ValueError(msg)

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


  def createAllXS(self, energy_groups, domain_type='distribcell'):

    if not domain_type in domain_types:
      msg = 'The XSTallyBuilder is unable to create all cross-sections for ' \
            'domain {0} since it is not a supported type'.format(domain_type)
      raise ValueError(msg)

    if domain_type == 'distribcell' or domain_type == 'cell':
      domains = self._geometry.getAllMaterialCells()
    elif domain_type == 'universe':
      domains = self._geometry.getAllMaterialUniverses()
    elif domain_type == 'material':
      domains = self._geometry.getAllMaterials()

    for domain in domains:
      for xs_type in xs_types:
        self.createXS(xs_type, energy_groups, domain, domain_type)


  def createAllMicroXS(self, energy_groups, domain_type='distribcell',
                       nuclides=None):

    if not domain_type in domain_types:
      msg = 'The XSTallyBuilder is unable to create all micro cross-sections ' \
            'for domain {0} since it is not a supported type'.format(domain_type)
      raise ValueError(msg)

    if domain_type == 'distribcell' or domain_type == 'cell':
      domains = self._geometry.getAllMaterialCells()
    elif domain_type == 'universe':
      domains = self._geometry.getAllMaterialUniverses()
    elif domain_type == 'material':
      domains = self._geometry.getAllMaterials()

    for domain in domains:
      for xs_type in xs_types:
        self.createMicroXS(xs_type, energy_groups, domain, domain_type, nuclides)


  def createTalliesFile(self):

    tallies = set()

    for xs in self._all_xs:
      tallies = tallies.union(set(xs._tallies.values()))

    for tally in tallies:
      self._tallies_file.addTally(tally)

    self._tallies_file.exportToXML()