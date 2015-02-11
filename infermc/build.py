from openmc import Nuclide
from openmc.tallies import TalliesFile
from infermc.multigroupxs import *
from infermc.microxs import *


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
  def geometry(self, geometry):
    self._geometry = geometry


  def createMultiGroupXS(self, xs_type, energy_groups, domain,
                         domain_type='distribcell'):

    if xs_type == 'total':
      xs = TotalXS(domain, domain_type, energy_groups)
    elif xs_type == 'transport':
      xs = TransportXS(domain, domain_type, energy_groups)
    elif xs_type == 'absorption':
      xs = AbsorptionXS(domain, domain_type, energy_groups)
    elif xs_type == 'capture':
      xs = CaptureXS(domain, domain_type, energy_groups)
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
    elif xs_type == 'nu-scatter matrix':
      xs = NuScatterMatrixXS(domain, domain_type, energy_groups)
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


  def createAllMultiGroupXS(self, energy_groups, domain_type='distribcell'):

    if domain_type == 'distribcell' or domain_type == 'cell':
      domains = self._geometry.get_all_material_cells()
    elif domain_type == 'universe':
      domains = self._geometry.get_all_material_universes()
    elif domain_type == 'material':
      domains = self._geometry.get_all_materials()

    for domain in domains:
      for xs_type in infermc.xs_types:
        self.createMultiGroupXS(xs_type, energy_groups, domain, domain_type)


  def createTalliesFile(self):

    tallies = set()

    for xs in self._all_xs:
      tallies = tallies.union(set(xs._tallies.values()))

    for tally in tallies:
      self._tallies_file.add_tally(tally)

    self._tallies_file.export_to_xml()


class MicroXSTallyFactory(XSTallyFactory):


  def createMultiGroupXS(self, xs_type, energy_groups, domain,
                    domain_type='distribcell', nuclides=None):

    # Add all of the Nuclides in the domain to the MicroXS
    if nuclides is None:

      nuclides_densities = domain.get_all_nuclides().values()

      nuclides = list()

      for nuclide_density in nuclides_densities:
        nuclides.append(nuclide_density[0])

      nuclides.append(Nuclide('total'))

    if xs_type == 'total':
      xs = MicroTotalXS(domain, domain_type, energy_groups, nuclides)
    elif xs_type == 'transport':
      xs = MicroTransportXS(domain, domain_type, energy_groups, nuclides)
    elif xs_type == 'absorption':
      xs = MicroAbsorptionXS(domain, domain_type, energy_groups, nuclides)
    elif xs_type == 'capture':
      xs = MicroCaptureXS(domain, domain_type, energy_groups, nuclides)
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
    elif xs_type == 'nu-scatter matrix':
      xs = MicroNuScatterMatrixXS(domain, domain_type, energy_groups, nuclides)
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


  def createAllMultiGroupXS(self, energy_groups, domain_type='distribcell',
                            nuclides=None):

    if domain_type == 'distribcell' or domain_type == 'cell':
      domains = self._geometry.get_all_material_cells()
    elif domain_type == 'universe':
      domains = self._geometry.get_all_material_universes()
    elif domain_type == 'material':
      domains = self._geometry.get_all_materials()

    for domain in domains:
      for xs_type in infermc.xs_types:
        self.createMultiGroupXS(xs_type, energy_groups, domain, domain_type, nuclides)