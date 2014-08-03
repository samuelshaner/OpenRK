__author__ = 'Will Boyd'
__email__ = 'wboyd@mit.edu'


from opencsg.checkvalue import *


# A list of all IDs for all Materials created
material_ids = list()

# A static variable for auto-generated Material IDs
auto_material_id = 10000


class Material(object):

  def __init__(self, material_id=None, name=''):

    # Initialize class attributes
    self._id = None
    self._name = ''

    # Set the Material class attributes
    self.setId(material_id)
    self.setName(name)


  def setId(self, material_id=None):

    global material_ids

    if material_id is None:
      global auto_material_id
      self._id = auto_material_id
      material_ids.append(auto_material_id)
      auto_material_id += 1

    # Check that the ID is an integer and wasn't already used
    elif is_integer(material_id):

      # If the Material already has an ID, remove it from global list
      if not self._id is None:
        material_ids.remove(self._id)

      if material_id in material_ids:
        msg = 'Unable to set Material ID to {0} since a Material ' \
              'with this ID was already initialized'.format(material_id)
        raise ValueError(msg)

      if material_id < 0:
        msg = 'Unable to set Material ID to {0} since it must be a ' \
              'non-negative integer'.format(material_id)
        raise ValueError(msg)

      else:
        self._id = material_id
        material_ids.append(material_id)

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
