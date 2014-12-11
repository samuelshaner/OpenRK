from opencg.checkvalue import *

# A static variable for auto-generated Material UIDs
AUTO_MATERIAL_UID = 1


class Material(object):

  def __init__(self, material_id=None, name=''):

    # Initialize class attributes
    global AUTO_MATERIAL_UID
    self._uid = AUTO_MATERIAL_UID
    AUTO_MATERIAL_UID += 1

    self._id = None
    self._set_id = False
    self._name = ''

    # Set the Material class attributes
    if not material_id is None:
      self.setId(material_id)
    else:
      self.setId(self._uid)
      self._set_id = False

    self.setName(name)


  def __deepcopy__(self, memo):

    existing = memo.get(id(self))

    # If this is the first time we have tried to copy this object, create a copy
    if existing is None:

      clone = type(self).__new__(type(self))
      clone._uid = self._uid
      clone._id = self._id
      clone._set_id = self._set_id
      clone._name = self._name

      memo[id(self)] = clone

      return clone

    # If this object has been copied before, return the first copy made
    else:
      return existing


  def __gt__(self, other):
    return (id(self) > id(other))


  def __ge__(self, other):
    return (id(self) >= id(other))


  def __lt__(self, other):
    return (id(self) < id(other))


  def __le__(self, other):
    return (id(self) <= id(other))


  def setId(self, material_id):

    # Check that the ID is a non-negative integer
    if is_integer(material_id):

      if material_id >= 0:
        self._id = material_id
        self._set_id = True
      else:
        msg = 'Unable to set Material ID to {0} since it must be a ' \
              'non-negative integer'.format(material_id)
        raise ValueError(msg)

    else:
      msg = 'Unable to set Material ID to non-integer {0}'.format(material_id)
      raise ValueError(msg)


  def setName(self, name):

    if not is_string(name):
      msg = 'Unable to set name for Material ID={0} with a non-string ' \
            'value {1}'.format(self._id, (name))
      raise ValueError(msg)

    else:
      self._name = name


  def __repr__(self):

    string = 'Material\n'
    string += '{0: <16}{1}{2}\n'.format('\tID', '=\t', self._id)
    string += '{0: <16}{1}{2}\n'.format('\tName', '=\t', self._name)

    return string
