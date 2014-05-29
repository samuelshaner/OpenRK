from openmc.tallies import Tally
from infermc.checkvalue import *


xs_types = ['total',
            'transport',
            'absorption',
            'scatter',
            'nu-scatter',
            'scatter matrix',
            'fission',
            'nu-fission',
            'chi']
#            'diffusion'

domain_types = ['cell',
                'distribcell',
                'universe',
                'material']


class EnergyGroups(object):

  def __init__(self):
    self._group_edges = None
    self._num_groups = None


  def setGroupEdges(self, edges):

    if not isinstance(edges, (tuple, list, np.ndarray)):
      msg = 'Unable to set energy group edges from {0} which ' \
            'is not a Python tuple/list or NumPy array'.format(edges)
      raise ValueError(msg)

    for i, edge in enumerate(edges):

      if not is_integer(edge) and not is_float(edge):
        msg = 'Unable to set energy group edge = {0} which is ' \
              'not an integer or floating point value'.format(edge)
        raise ValueError(msg)

      if i > 0 and edges[i] <= edges[i-1]:
        msg = 'Unable to set energy group edges with {0} for edge {1} ' \
              'and {2} for edge {3}'.format(edges[i], i, edges[i-1], i-1)
        raise ValueError(msg)

    self._group_edges = np.array(edges)
    self._num_groups = len(edges)-1


  def generateBinEdges(self, start, stop, num_groups, type='linear'):

    if not is_integer(start) and not is_float(start):
      msg = 'Unable to generate energy group edges with start = {0} ' \
            'which is not an integer or floating point value'.format(start)
      raise ValueError(msg)

    if not is_integer(stop) and not is_float(stop):
      msg = 'Unable to generate energy group edges with stop = {0} ' \
            'which is not an integer or floating point value'.format(stop)
      raise ValueError(msg)

    if not is_integer(num_groups):
      msg = 'Unable to generate energy group edges with num groups = {0} ' \
            'which is not an integer value'.format(num_groups)
      raise ValueError(msg)

    if not type in ['linear', 'logarithmic']:
      msg = 'Unable to generate energy group edges with type = {0} which is ' \
            'neither linear nor logarithmic'.format(type)
      raise ValueError(msg)

    if start < 0.:
      msg = 'Unable to generate energy group edges with start = {0} which is ' \
           'a negative value'.format(start)
      raise ValueError(msg)

    if stop < 0.:
      msg = 'Unable to generate energy group edges with stop = {0} which is ' \
            'a negative value'.format(stop)
      raise ValueError(msg)

    if start >= stop:
      msg = 'Unable to generate energy group edges with start = {0} which is ' \
           'greater than stop = {1}'.format(start, stop)
      raise ValueError(msg)

    if num_groups < 0.:
      msg = 'Unable to generate energy group edges with num groups {0} ' \
            'which is a negative value'.format(num_groups)
      raise ValueError(msg)

    self._num_groups = num_groups

    if type == 'linear':
      self._group_edges = np.linspace(start, stop, num_groups+1)

    else:
      self._group_edges = np.linspace(np.log(start), np.log(stop), num_groups+1)


  def getGroup(self, energy):

    # Assumes energy is in eV

    if not is_integer(energy) and not is_float(energy):
      msg = 'Unable to get energy group for energy {0} since it is ' \
            'neither an integer or floating point value'.format(energy)
      raise ValueError(msg)

    if self._group_edges is None:
      msg = 'Unable to get energy group for energy {0} eV since ' \
            'the group edges have not yet been set'.format(energy)
      raise ValueError(msg)

    # Loop over all edges and search for the group for this energy
    for i, edge in enumerate(self._group_edges):

      if energy <= edge:
        return self._num_groups - i

    msg = 'Unable to find energy group for energy {0} eV since it does not ' \
          'lie within the bounds of the energy group edges'.format(energy)
    raise ValueError(msg)


class MultiGroupXS(object):

  def __init__(self):
    self._type = None
    self._energy_groups = None
    self._tallies = list()
    self._data = dict()
    self._xs_estimate = None
    self._xs_uncertainty = None


  def setEnergyGroups(self, energy_groups):

    if not isinstance(energy_groups, EnergyGroups):
      msg = 'Unable to set the energy groups to {0} which is not an ' \
           'EnergyGroups object'.format(energy_groups)
      raise ValueError(msg)

    self._energy_groups = energy_groups


  def createTallies(self):

    if self._energy_groups is None:
      msg = 'Unable to create tallies for a MultiGroupXS since the ' \
            'energy groups has not been set'
      raise ValueError(msg)

    group_edges = self._energy_groups._group_edges

    flux = Tally()
    flux.addScore(score='flux')
    flux.addFilter(type='energy', bins=group_edges)
    flux.setLabel(label='%d groups' % self._energy_groups._num_groups)
    self._tallies.append(flux)

    return


class TotalXS(MultiGroupXS):

  def __init__(self):
    self._type = 'total'
    super(TotalXS, self).__init__()


  def createTallies(self):

    super(TotalXS, self).createTallies()

    group_edges = self._energy_groups._group_edges

    rxn_rate = Tally()
    rxn_rate.addScore(score='total')
    rxn_rate.addFilter(type='energy', bins=group_edges)
    rxn_rate.setLabel(label='%d groups' % self._energy_groups._num_groups)
    self._tallies.append(rxn_rate)


class TransportXS(MultiGroupXS):

  def __init__(self):
    self._type = 'transport'
    super(TransportXS, self).__init__()

  def createTallies(self):

    super(TransportXS, self).createTallies()

    group_edges = self._energy_groups._group_edges

    flux = Tally()
    flux.addScore(score='flux')
    flux.addFilter(type='energy', bins=group_edges)
    flux.setLabel(label='%d groups' % self._energy_groups._num_groups)
    flux.setEstimator(estimator='analog')
    self._tallies.append(flux)

    rxn_rate = Tally()
    rxn_rate.addScore(score='total')
    rxn_rate.addFilter(type='energy', bins=group_edges)
    rxn_rate.setLabel(label='%d groups' % self._energy_groups._num_groups)
    rxn_rate.setEstimator(estimator='analog')
    self._tallies.append(rxn_rate)

    rxn_rate = Tally()
    rxn_rate.addScore(score='scatter-P1')
    rxn_rate.addFilter(type='energy', bins=group_edges)
    rxn_rate.setLabel(label='%d groups' % self._energy_groups._num_groups)
    rxn_rate.setEstimator(estimator='analog')
    self._tallies.append(rxn_rate)


class AbsorptionXS(MultiGroupXS):

  def __init__(self):
    self._type = 'absorption'
    super(AbsorptionXS, self).__init__()


  def createTallies(self):

    super(AbsorptionXS, self).createTallies()

    group_edges = self._energy_groups._group_edges

    rxn_rate = Tally()
    rxn_rate.addScore(score='absorption')
    rxn_rate.addFilter(type='energy', bins=group_edges)
    rxn_rate.setLabel(label='%d groups' % self._energy_groups._num_groups)
    self._tallies.append(rxn_rate)


class FissionXS(MultiGroupXS):

  def __init__(self):
    self._type = 'fission'
    super(FissionXS, self).__init__()


  def createTallies(self):

    super(FissionXS, self).createTallies()

    group_edges = self._energy_groups._group_edges

    rxn_rate = Tally()
    rxn_rate.addScore(score='fission')
    rxn_rate.addFilter(type='energy', bins=group_edges)
    rxn_rate.setLabel(label='%d groups' % self._energy_groups._num_groups)
    self._tallies.append(rxn_rate)


class NuFissionXS(MultiGroupXS):

  def __init__(self):
    self._type = 'nu-fission'
    super(NuFissionXS, self).__init__()


  def createTallies(self):

    super(NuFissionXS, self).createTallies()

    group_edges = self._energy_groups._group_edges

    rxn_rate = Tally()
    rxn_rate.addScore(score='nu-fission')
    rxn_rate.addFilter(type='energy', bins=group_edges)
    rxn_rate.setLabel(label='%d groups' % self._energy_groups._num_groups)
    self._tallies.append(rxn_rate)


class ScatterXS(MultiGroupXS):

  def __init__(self):
    self._type = 'scatter'
    super(ScatterXS, self).__init__()

  def createTallies(self):

    super(ScatterXS, self).createTallies()

    group_edges = self._energy_groups._group_edges

    rxn_rate = Tally()
    rxn_rate.addScore(score='scatter')
    rxn_rate.addFilter(type='energy', bins=group_edges)
    rxn_rate.setLabel(label='%d groups' % self._energy_groups._num_groups)
    self._tallies.append(rxn_rate)


class NuScatterXS(MultiGroupXS):

  def __init__(self):
    self._type = 'nu-scatter'
    super(NuScatterXS, self).__init__()

  def createTallies(self):

    super(NuScatterXS, self).createTallies()

    group_edges = self._energy_groups._group_edges

    flux = Tally()
    flux.addScore(score='flux')
    flux.addFilter(type='energy', bins=group_edges)
    flux.setLabel(label='%d groups' % self._energy_groups._num_groups)
    flux.setEstimator(estimator='analog')
    self._tallies.append(flux)

    rxn_rate = Tally()
    rxn_rate.addScore(score='nu-scatter')
    rxn_rate.addFilter(type='energy', bins=group_edges)
    rxn_rate.setLabel(label='%d groups' % self._energy_groups._num_groups)
    rxn_rate.setEstimator(estimator='analog')
    self._tallies.append(rxn_rate)


class ScatterMatrixXS(MultiGroupXS):

  def __init__(self):
    self._type = 'scatter matrix'
    super(ScatterMatrixXS, self).__init__()

  def createTallies(self):

    super(ScatterMatrixXS, self).createTallies()

    group_edges = self._energy_groups._group_edges

    flux = Tally()
    flux.addScore(score='flux')
    flux.addFilter(type='energy', bins=group_edges)
    flux.setLabel(label='%d groups' % self._energy_groups._num_groups)
    flux.setEstimator(estimator='analog')
    self._tallies.append(flux)

    rxn_rate = Tally()
    rxn_rate.addScore(score='nu-scatter')
    rxn_rate.addFilter(type='energy', bins=group_edges)
    rxn_rate.addFilter(type='energyout', bins=group_edges)
    rxn_rate.setLabel(label='%d groups' % self._energy_groups._num_groups)
    rxn_rate.setEstimator(estimator='analog')
    self._tallies.append(rxn_rate)


class DiffusionCoeff(MultiGroupXS):

  def __init__(self):
    self._type = 'diffusion'
    super(DiffusionCoeff, self).__init__()

  def createTallies(self):

    super(DiffusionCoeff, self).createTallies()

    msg = 'DiffusionCoeff is not yet able to build tallies'
    raise ValueError(msg)


class Chi(MultiGroupXS):

  def __init__(self):
    self._type = 'chi'
    super(Chi, self).__init__()

  def createTallies(self):

    super(Chi, self).createTallies()

    group_edges = self._energy_groups._group_edges

    rxn_rate = Tally()
    rxn_rate.addScore(score='nu-fission')
    rxn_rate.addFilter(type='energy', bins=group_edges)
    rxn_rate.setEstimator(estimator='analog')
    rxn_rate.setLabel(label='%d groups' % self._energy_groups._num_groups)

    emit_rate = Tally()
    emit_rate.addScore(score='nu-fission')
    emit_rate.addFilter(type='energyout', bins=group_edges)
    emit_rate.setEstimator(estimator='analog')
    emit_rate.setLabel(label='%d groups' % self._energy_groups._num_groups)

    self._tallies.append(rxn_rate)
    self._tallies.append(emit_rate)
