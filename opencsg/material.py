__author__ = 'Will Boyd'
__email__ = 'wboyd@mit.edu'


from opencsg.checkvalue import *


# A list of all IDs for all Materials created
MATERIAL_IDS = list()

# A static variable for auto-generated Material IDs
AUTO_MATERIAL_ID = 10000


class Material(object):

  def __init__(self, material_id=None, name=''):

    # Initialize class attributes
    self._id = None
    self._name = ''

    # Set the Material class attributes
    self.setId(material_id)
    self.setName(name)


  def __deepcopy__(self, memo):

    existing = memo.get(self)

    # If this is the first time we have tried to copy this object, create a copy
    if existing is None:

      clone = type(self).__new__(type(self))
      clone._id = self._id
      clone._name = self._name

      return clone

    # If this object has been copied before, return the first copy made
    else:
      return existing


  def setId(self, material_id=None):

    global MATERIAL_IDS

    if material_id is None:
      global AUTO_MATERIAL_ID
      self._id = AUTO_MATERIAL_ID
      MATERIAL_IDS.append(AUTO_MATERIAL_ID)
      AUTO_MATERIAL_ID += 1

    # Check that the ID is an integer and wasn't already used
    elif is_integer(material_id):

      # If the Material already has an ID, remove it from global list
      if not self._id is None:
        MATERIAL_IDS.remove(self._id)

      if material_id in MATERIAL_IDS:
        msg = 'Unable to set Material ID to {0} since a Material ' \
              'with this ID was already initialized'.format(material_id)
        raise ValueError(msg)

      if material_id < 0:
        msg = 'Unable to set Material ID to {0} since it must be a ' \
              'non-negative integer'.format(material_id)
        raise ValueError(msg)

      else:
        self._id = material_id
        MATERIAL_IDS.append(material_id)

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
