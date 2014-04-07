import numpy as np


# A list of all IDs for all Materials created
material_ids = list()

# A static variable for auto-generated Material IDs
auto_material_id = 10000

# Units for density supported by OpenMC
density_units = ['g/cm3', 'kg/cm3', 'at/b-cm', 'at/cm3', 'sum']


class Material(object):

    def __init__(self, material_id=None, name=''):

        # Initialize class attributes
        self._id = None
        self._name = ''

        # Set the Material class attributes
        self.setId(material_id)
        self.setName(name)


    def getId(self):
        return self._id


    def getName(self):
        return self._name


    def setId(self, material_id=None):

        global material_ids

        if material_id is None:
            global auto_material_id
            self._id = auto_material_id
            material_ids.append(auto_material_id)
            auto_material_id += 1

        # Check that the ID is an integer and wasn't already used
        elif isinstance(material_id, (int, np.int32, np.int64)):

            # If the Material already has an ID, remove it from global list
            if not self._id is None:
                material_ids.remove(self._id)

            if material_id in material_ids:
                exit('Unable to set Material ID to %s since a Material '
                      'with this ID was already initialized.', str(material_id))

            if material_id < 0:
                exit('Unable to set Material ID to %d since it must be a '
                     'non-negative integer', material_id)

            else:
                self._id = material_id
                material_ids.append(material_id)

        else:
            exit('Unable to set a non-integer Material ID %s', str(material_id))


    def setName(self, name):

        if not isinstance(name, str):
            exit('Unable to set name for Material ID=%d with a non-string '
                 'value %s', self._id, str(name))

        else:
            self._name = name


    def toString(self):

        string = ''

        string += 'Material\n'

        material_id = '{0: <16}'.format('\tID') + '=\t' + str(self._id)
        string += material_id + '\n'

        name = '{0: <16}'.format('\tName') + '=\t' + self._name
        string += name + '\n'

        return string


    def printString(self):
        print(self.toString())