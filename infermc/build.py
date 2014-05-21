from opencsg import Cell
from opencsg.checkvalue import *
from openmc.input.tallies import Tally, TalliesFile
from sets import Set


xs_types = ['total',
            'absorption',
            'scatter',
            'scatter matrix',
            'fission',
            'nu-fission',
            'chi']

#            'diffusion',
#            'transport']


class EnergyGroups(object):

  def __init__(self):
    self._groups_edges = None
    self._num_groups = None


  def setGroupEdges(self, edges):

    if not isinstance(edges, (tuple, list, np.ndarray)):
      exit('Unable to set energy group edges from %s which '
           'is not a Python tuple/list or NumPy array' % str(edges))

    for i, edge in enumerate(edges):

      if not is_integer(edge) and not is_float(edge):
        exit('Unable to set energy group edge = %s which is '
             'not an integer or floating point value' % str(edge))

      if i > 0 and edges[i] <= edges[i-1]:
        exit('Unable to set energy group edges with %s for edge %d '
             'and %s for edge %d' % (str(edges[i], i, edges[i-1], i-1)))

    self._groups_edges = np.array(edges)
    self._num_groups = len(edges)-1


  def generateBinEdges(self, start, stop, num_groups, type='linear'):

    if not is_integer(start) and not is_float(start):
      exit('Unable to generate energy group edges with start = %s '
           'which is not an integer or floating point value' % str(start))

    if not is_integer(stop) and not is_float(stop):
      exit('Unable to generate energy group edges with stop = %s '
           'which is not an integer or floating point value' % str(stop))

    if not is_integer(num_groups):
      exit('Unable to generate energy group edges with num groups = %s '
           'which is not an integer value' % str(num_groups))

    if not type in ['linear', 'logarithmic']:
      exit('Unable to generate energy group edges with type = %s which is '
           'neither linear nor logarithmic' % str(type))

    if start < 0.:
      exit('Unable to generate energy group edges with start = %s which is '
           'a negative value' % str(start))

    if stop < 0.:
      exit('Unable to generate energy group edges with stop = %s which is '
           'a negative value' % str(stop))

    if start >= stop:
      exit('Unable to generate energy group edges with start = %s which is '
           'greater than stop = %s' % (str(start), str(stop)))

    if num_groups < 0.:
      exit('Unable to generate energy group edges with num groups %s '
           'which is a negative value' % str(num_groups))

    self._num_groups = num_groups

    if type == 'linear':
      self._groups_edges = np.linspace(start, stop, num_groups+1)

    else:
      self._group_edges = np.linspace(np.log(start), np.log(stop), num_groups+1)


  def getNumGroups(self):

    if self._groups_edges is None:
      exit('Unable to get the number of energy groups since '
           'the group edges have not yet been set')

    return self._num_groups


  def getGroupEdges(self):
    return self._groups_edges


  def getGroup(self, energy):

    # Assumes energy is in eV

    if not is_integer(energy) and not is_float(energy):
      exit('Unable to get energy group for energy %s since it is '
           'neither an integer or floating point value' % str(energy))

    if self._groups_edges is None:
      exit('Unable to get energy group for energy %s eV since '
           'the group edges have not yet been set' % str(energy))

    # Loop over all edges and search for the group for this energy
    for i, edge in enumerate(self._groups_edges):

      if energy <= edge:
        return self._num_groups - i

    exit('Unable to find energy group for energy %s eV since it does not '
         'lie within the bounds of the energy group edges' % str(energy))


class MultiGroupXS(object):

  def __init__(self):
    self._type = None
    self._energy_groups = None
    self._tallies = list()


  def getType(self):
    self._type


  def getEnergyGroups(self):
    return self._energy_groups


  def getNumGroups(self):

    if self._energy_groups is None:
      exit('Unable to get the number of energy groups for multi group xs'
           'since the energy groups has not yet been set')

    return self._energy_groups.getNumGroups()


  def getTallies(self):
    return self._tallies


  def setEnergyGroups(self, energy_groups):

    if not isinstance(energy_groups, EnergyGroups):
      exit('Unable to set the energy groups to %s which is not an '
           'EnergyGroups object' % str(energy_groups))

    self._energy_groups = energy_groups


  def createTallies(self):

    if self._energy_groups is None:
      exit('Unable to create tallies for a MultiGroupXS since the '
           'energy groups has not been set')

    group_edges = self._energy_groups.getGroupEdges()

    flux = Tally()
    flux.addScore(score='flux')
    flux.addFilter(type='energy', bins=group_edges)
    flux.setLabel(label='%d groups' % self._energy_groups.getNumGroups())
    self._tallies.append(flux)

    return


class TotalXS(MultiGroupXS):

  def __init__(self):
    self._type = 'total'
    super(TotalXS, self).__init__()


  def createTallies(self):

    super(TotalXS, self).createTallies()

    group_edges = self._energy_groups.getGroupEdges()

    rxn_rate = Tally()
    rxn_rate.addScore(score='total')
    rxn_rate.addFilter(type='energy', bins=group_edges)
    rxn_rate.setLabel(label='%d groups' % self._energy_groups.getNumGroups())
    self._tallies.append(rxn_rate)


class AbsorptionXS(MultiGroupXS):

  def __init__(self):
    self._type = 'absorption'
    super(AbsorptionXS, self).__init__()


  def createTallies(self):

    super(AbsorptionXS, self).createTallies()

    group_edges = self._energy_groups.getGroupEdges()

    rxn_rate = Tally()
    rxn_rate.addScore(score='absorption')
    rxn_rate.addFilter(type='energy', bins=group_edges)
    rxn_rate.setLabel(label='%d groups' % self._energy_groups.getNumGroups())
    self._tallies.append(rxn_rate)


class FissionXS(MultiGroupXS):

  def __init__(self):
    self._type = 'fission'
    super(FissionXS, self).__init__()


  def createTallies(self):

    super(FissionXS, self).createTallies()

    group_edges = self._energy_groups.getGroupEdges()

    rxn_rate = Tally()
    rxn_rate.addScore(score='fission')
    rxn_rate.addFilter(type='energy', bins=group_edges)
    rxn_rate.setLabel(label='%d groups' % self._energy_groups.getNumGroups())
    self._tallies.append(rxn_rate)


class NuFissionXS(MultiGroupXS):

  def __init__(self):
    self._type = 'nu-fission'
    super(NuFissionXS, self).__init__()


  def createTallies(self):

    super(NuFissionXS, self).createTallies()

    group_edges = self._energy_groups.getGroupEdges()

    rxn_rate = Tally()
    rxn_rate.addScore(score='nu-fission')
    rxn_rate.addFilter(type='energy', bins=group_edges)
    rxn_rate.setLabel(label='%d groups' % self._energy_groups.getNumGroups())
    self._tallies.append(rxn_rate)


class ScatterXS(MultiGroupXS):

  def __init__(self):
    self._type = 'scatter'
    super(ScatterXS, self).__init__()

  def createTallies(self):

    super(ScatterXS, self).createTallies()

    group_edges = self._energy_groups.getGroupEdges()

    rxn_rate = Tally()
    rxn_rate.addScore(score='scatter')
    rxn_rate.addFilter(type='energy', bins=group_edges)
    rxn_rate.setLabel(label='%d groups' % self._energy_groups.getNumGroups())
    self._tallies.append(rxn_rate)


class NuScatterXS(MultiGroupXS):

  def __init__(self):
    self._type = 'nu-scatter'
    super(ScatterXS, self).__init__()

  def createTallies(self):

    super(ScatterXS, self).createTallies()

    group_edges = self._energy_groups.getGroupEdges()

    rxn_rate = Tally()
    rxn_rate.addScore(score='nu-scatter')
    rxn_rate.addFilter(type='energy', bins=group_edges)
    rxn_rate.setLabel(label='%d groups' % self._energy_groups.getNumGroups())
    self._tallies.append(rxn_rate)


class ScatterMatrixXS(MultiGroupXS):

  def __init__(self):
    self._type = 'scatter matrix'
    super(ScatterMatrixXS, self).__init__()

  def createTallies(self):

    super(ScatterMatrixXS, self).createTallies()

    group_edges = self._energy_groups.getGroupEdges()

    rxn_rate = Tally()
    rxn_rate.addScore(score='nu-scatter')
    rxn_rate.addFilter(type='energy', bins=group_edges)
    rxn_rate.addFilter(type='energyout', bins=group_edges)
    rxn_rate.setLabel(label='%d groups' % self._energy_groups.getNumGroups())
    self._tallies.append(rxn_rate)


class TransportXS(MultiGroupXS):

  def __init__(self):
    self._type = 'transport'
    super(TransportXS, self).__init__()

  def createTallies(self):

    super(TransportXS, self).createTallies()

    exit('TransportXS is not yet able to build tallies')


class DiffusionCoeff(MultiGroupXS):

  def __init__(self):
    self._type = 'diffusion'
    super(DiffusionCoeff, self).__init__()

  def createTallies(self):

    super(DiffusionCoeff, self).createTallies()

    exit('DiffusionCoeff is not yet able to build tallies')


class Chi(MultiGroupXS):

  def __init__(self):
    self._type = 'chi'
    super(Chi, self).__init__()

  def createTallies(self):

    super(Chi, self).createTallies()

    group_edges = self._energy_groups.getGroupEdges()

    rxn_rate = Tally()
    rxn_rate.addScore(score='nu-fission')
    rxn_rate.addFilter(type='energy', bins=group_edges)
    rxn_rate.setLabel(label='%d groups' % self._energy_groups.getNumGroups())

    emit_rate = Tally()
    emit_rate.addScore(score='nu-fission')
    emit_rate.addFilter(type='energyout', bins=group_edges)
    emit_rate.setLabel(label='%d groups' % self._energy_groups.getNumGroups())

    self._tallies.append(rxn_rate)
    self._tallies.append(emit_rate)


class XSTallyBuilder(object):

  def __init__(self):

    self._all_xs = list()
    self._tallies_file = TalliesFile()


  def getTalliesFile(self):
    return self._tallies_file


  def createXS(self, xs_type, energy_groups, cell):

    global xs_types

    if not is_string(xs_type):
      exit('The XSTallyBuilder is unable to create cross-section type %s '
           'which is not a string value' % str(xs_type))

    if not isinstance(energy_groups, EnergyGroups):
      exit('The XSTallyBuilder is unable to create cross-section type %s '
           'with energy groups %s which is not an EnergyGroups object' %
           (xs_type, str(energy_groups)))

    if not isinstance(cell, Cell):
      exit('The XSTallyBuilder is unable to create cross-section type %s '
           'with cell %s which is not an OpenCSG Cell object' %
           (xs_type, str(cell)))

    if xs_type == 'total':
      xs = TotalXS()
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
    elif xs_type == 'transport':
      xs = TransportXS()
    elif xs_type == 'chi':
      xs = Chi()
    else:
      exit('The XSTallyBuilder is unable to create cross-section type %s '
           'which is not one of the supported types' % xs_type)

    xs.setEnergyGroups(energy_groups)
    xs.createTallies()

    for tally in xs.getTallies():
      tally.addFilter(type='distribcell', bins=cell.getId())

    self._all_xs.append(xs)


  def createAllXS(self, energy_groups, cells):

    global xs_types

    if not isinstance(cells, (tuple, list, np.ndarray)):
      cells = [cells]

    for cell in cells:
      for xs_type in xs_types:
        self.createXS(xs_type, energy_groups, cell)


  def createTalliesFile(self):

    tallies = Set()

    for xs in self._all_xs:
      tallies.union_update(Set(xs.getTallies()))

    for tally in tallies:
      self._tallies_file.addTally(tally)

    self._tallies_file.exportToXML()