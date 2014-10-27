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


  def __eq__(self, other):
    if not isinstance(other, EnergyGroups):
      return False
    elif self._group_edges != other._group_edges:
      return False


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


  def getCondensedGroups(self, coarse_groups):
    '''This routine takes in a collection of 2-tuples of energy groups'''

    if not isinstance(coarse_groups, (tuple, list, np.array)):
      msg = 'Unable to condense MultiGroupXS with group_bounds {0} which ' \
            'is not a Python tuple/list or NumPy array'.format(coarse_groups)
      raise ValueError(msg)

    for group in coarse_groups:
      if group[0] < 1 or group[0] > self._num_groups:
        msg = 'Unable to condense MultiGroupXS with group bound {0}'.format(group)
        raise ValueError(msg)
      elif group[1] < 1 or group[1] > self._num_groups:
        msg = 'Unable to condense MultiGroupXS with group bound {0}'.format(group)
        raise ValueError(msg)
      elif group[0] >= group[1]:
        msg = 'Unable to condense MultiGroupXS with groups {0} which ' \
              'is not monotonically increasing'.format(group)
        raise ValueError(msg)

    # Compute the group indices into the coarse group
    group_bounds = list()
    for group in coarse_groups:
      group_bounds.append(group[0])
    group_bounds.append(coarse_groups[-1][1])

    # Determine the indices mapping the fine-to-coarse energy groups
    group_bounds = np.asarray(group_bounds)
    group_indices = np.flipud(self._num_groups - group_bounds)
    group_indices[-1] += 1

    # Determine the edges between coarse energy groups and sort
    # in increasing order in case the user passed in unordered groups
    group_edges = self._group_edges[group_indices]
    group_edges = np.sort(group_edges)

    # Create a new condensed EnergyGroups object
    condensed_groups = EnergyGroups()
    condensed_groups.group_edges = group_edges
    return condensed_groups