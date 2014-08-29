from typecheck import accepts, Or, Exact, Self
import numpy as np


class EnergyGroups(object):

  def __init__(self):

    self._group_edges = None
    self._num_groups = None


  @property
  def group_edges(self):
    return self._group_edges


  @group_edges.setter
  @accepts(Self(), Or(list, tuple, np.ndarray))
  def group_edges(self, edges):
    self._group_edges = np.array(edges)
    self._num_groups = len(edges)-1


  @accepts(Self(), Or(int, float), Or(int, float), int,
           Or(Exact('linear'), Exact('logarithmic')))
  def generateBinEdges(self, start, stop, num_groups, type='linear'):

    self._num_groups = num_groups

    if type == 'linear':
      self._group_edges = np.linspace(start, stop, num_groups+1)
    else:
      self._group_edges = np.logspace(np.log(start), np.log(stop), num_groups+1)


  @accepts(Self(), Or(int, float))
  def getGroup(self, energy):

    # Assumes energy is in eV

    if self._group_edges is None:
      msg = 'Unable to get energy group for energy {0} eV since ' \
            'the group edges have not yet been set'.format(energy)
      raise ValueError(msg)

    index = np.where(self._group_edges > energy)[0]
    group = self._num_groups - index
    return group


  @accepts(Self(), int)
  def getGroupBounds(self, group):

    if self._group_edges is None:
      msg = 'Unable to get energy group bounds for group {0} since ' \
            'the group edges have not yet been set'.format(group)
      raise ValueError(msg)

    lower = self._group_edges[self._num_groups-group]
    upper = self._group_edges[self._num_groups-group+1]
    return (lower, upper)
