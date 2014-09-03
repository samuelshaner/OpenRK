import numpy as np


class EnergyGroups(object):

  def __init__(self):

    self._group_edges = None
    self._num_groups = None


  @property
  def group_edges(self):
    return self._group_edges


  @group_edges.setter
  def group_edges(self, edges):
    self._group_edges = np.array(edges)
    self._num_groups = len(edges)-1


  def generateBinEdges(self, start, stop, num_groups, type='linear'):

    self._num_groups = num_groups

    if type == 'linear':
      self._group_edges = np.linspace(start, stop, num_groups+1)
    else:
      self._group_edges = np.logspace(np.log(start), np.log(stop), num_groups+1)


  def getGroup(self, energy):

    # Assumes energy is in eV

    if self._group_edges is None:
      msg = 'Unable to get energy group for energy {0} eV since ' \
            'the group edges have not yet been set'.format(energy)
      raise ValueError(msg)

    index = np.where(self._group_edges > energy)[0]
    group = self._num_groups - index
    return group


  def getGroupBounds(self, group):

    if self._group_edges is None:
      msg = 'Unable to get energy group bounds for group {0} since ' \
            'the group edges have not yet been set'.format(group)
      raise ValueError(msg)

    lower = self._group_edges[self._num_groups-group]
    upper = self._group_edges[self._num_groups-group+1]
    return (lower, upper)


  def getGroupIndices(self, groups='all'):

    if self._group_edges is None:
      msg = 'Unable to get energy group indices for groups {0} since ' \
            'the group edges have not yet been set'.format(groups)
      raise ValueError(msg)

    if groups == 'all':
      indices = np.arange(self._num_groups)

    else:

      indices = np.zeros(len(groups), dtype=np.int64)

      for i, group in enumerate(groups):

        if group > 0 and group <= self._num_groups:
          indices[i] = group - 1
        else:
          msg = 'Unable to get energy group index for group {0} since ' \
                'it is outside the group bounds'.format(group)
          raise ValueError(msg)

    return indices
