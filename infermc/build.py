from opencsg import Material, Cell, Universe, Geometry
from openmc.tallies import TalliesFile
from infermc.multigroupxs import *


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
            'it is not an OpenCSG Geometry object'.format(geometry)
      raise ValueError(msg)

    self._geometry = geometry


  def createXS(self, xs_type, energy_groups, domain, domain_type='distribcell'):

    global xs_types, domain_types

    if not is_string(xs_type):
      msg = 'The XSTallyBuilder is unable to create cross-section type {0} ' \
            'which is not a string value'.format(xs_type)
      raise ValueError(msg)

    if not isinstance(energy_groups, EnergyGroups):
      msg = 'The XSTallyBuilder is unable to create cross-section ' \
            'type {0} with energy groups {1} which is not an EnergyGroups ' \
            'object'.format(xs_type, energy_groups)
      raise ValueError(msg)

    if not isinstance(domain, (Material, Cell, Universe)):
      msg = 'The XSTallyBuilder is unable to create cross-section ' \
            'type {0} in {1} {2} since the domain is not a Material, ' \
            'Cell or Universe'.format(xs_type, domain_type, domain)
      raise ValueError(msg)

    if not domain_type in domain_types:
      msg = 'The XSTallyBuilder is unable to create cross-section ' \
            'type {0} in {1} {2} since it is not a supported domain ' \
            'type'.format(xs_type, domain_type, domain._id)
      raise ValueError(msg)

    if xs_type == 'total':
      xs = TotalXS()
    elif xs_type == 'transport':
      xs = TransportXS()
    elif xs_type == 'absorption':
      xs = AbsorptionXS()
    elif xs_type == 'fission':
      xs = FissionXS()
    elif xs_type == 'nu-fission':
      xs = NuFissionXS()
    elif xs_type == 'scatter':
      xs = ScatterXS()
    elif xs_type == 'nu-scatter':
      xs = NuScatterXS()
    elif xs_type == 'scatter matrix':
      xs = ScatterMatrixXS()
    elif xs_type == 'diffusion':
      xs = DiffusionCoeff()
    elif xs_type == 'chi':
      xs = Chi()
    else:
      msg = 'The XSTallyBuilder is unable to create cross-section type {0} ' \
            'which is not one of the supported types'.format(xs_type)
      raise ValueError(msg)

    xs.setEnergyGroups(energy_groups)
    xs.createTallies()

    # Create a Filter for the Tally domain
    filter = Filter(domain_type, domain._id)

    # Add the Tally domain Filter to all Tallies in the cross-section
    for tally in xs._tallies:
      tally.addFilter(filter)

    self._all_xs.append(xs)


  def createAllXS(self, energy_groups, domain_type='distribcell'):

    global xs_types, domain_types

    if not domain_type in domain_types:
      msg = 'The XSTallyBuilder is unable to create all cross-sections for ' \
            'domain {0} since it is not a supported types'.format(domain_type)
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


  def createTalliesFile(self):

    tallies = set()

    for xs in self._all_xs:
      tallies = tallies.union(set(xs._tallies))

    for tally in tallies:
      self._tallies_file.addTally(tally)

    self._tallies_file.exportToXML()