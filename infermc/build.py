from opencsg import Material, Cell, Universe, Geometry
from openmc.tallies import TalliesFile
from infermc.multigroupxs import *


class XSTallyBuilder(object):

  def __init__(self, geometry=None):

    self._all_xs = list()
    self._tallies_file = TalliesFile()
    self._geometry = None

    if not geometry is None:
      self.setGeometry(geometry)


  def getTalliesFile(self):
    return self._tallies_file


  def getGeometry(self):
    return self._geometry


  def setGeometry(self, geometry):

    if not isinstance(geometry, Geometry):
      exit('Unable to set the Geometry for XSTallyBuilder to %s since it is '
           'not an OpenCSG Geometry object' % str(geometry))

    self._geometry = geometry


  def createXS(self, xs_type, energy_groups, domain, domain_type='distribcell'):

    global xs_types, domain_types

    if not is_string(xs_type):
      exit('The XSTallyBuilder is unable to create cross-section type %s '
           'which is not a string value' % str(xs_type))

    if not isinstance(energy_groups, EnergyGroups):
      exit('The XSTallyBuilder is unable to create cross-section type %s '
           'with energy groups %s which is not an EnergyGroups object' %
           (xs_type, str(energy_groups)))

    if not isinstance(domain, (Material, Cell, Universe)):
      exit('The XSTallyBuilder is unable to create cross-section type %s '
           'in %s %s since the domain is not a Material, Cell or Universe' %
           (xs_type, str(domain_type), str(domain)))

    if not domain_type in domain_types:
      exit('The XSTallyBuilder is unable to create cross-section type %s '
           'in %s %d since it is not a supported domain type' %
           (xs_type, str(domain_type), domain._id))

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
      exit('The XSTallyBuilder is unable to create cross-section type %s '
           'which is not one of the supported types' % xs_type)

    xs.setEnergyGroups(energy_groups)
    xs.createTallies()

    for tally in xs._tallies:
      tally.addFilter(type=domain_type, bins=domain._id)

    self._all_xs.append(xs)


  def createAllXS(self, energy_groups, domain_type='distribcell'):

    global xs_types, domain_types

    if not domain_type in domain_types:
      exit('The XSTallyBuilder is unable to create all cross-sections for '
           'domain %s since it is not a supported types' % str(domain_types))

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
#      tallies.union_update(set(xs._tallies))

    for tally in tallies:
      self._tallies_file.addTally(tally)

    self._tallies_file.exportToXML()