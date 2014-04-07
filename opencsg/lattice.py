__author__ = 'Will Boyd'
__email__ = 'wboyd@mit.edu'


import numpy as np
from cell import universe_ids, auto_universe_id


# Lists of all IDs for all Lattices (Universes) created
global universe_ids


# A static variable for auto-generated Lattice (Universe) IDs
global auto_universe_id


class Lattice(object):

    def __init__(self, lattice_id=None, name='', type='rectangular'):

        # Initialize Lattice class attributes
        self._id = None
        self._name = None
        self._type = ''
        self._dimension = None
        self._lower_left = None
        self._width = None
        self._outside = None
        self._universes = None

        self.setId(lattice_id)
        self.setName(name)
        self.setType(type)

    def getId(self):
        return self._id

    def getName(self):
        return self._name

    def getType(self):
        return self._type

    def getLowerLeft(self):
        return self._lower_left

    def getWidth(self):
        return self._width

    def getOutside(self):
        return self._outside

    def getUniverses(self):
        return self._universes

    def setId(self, lattice_id=None):

        global universe_ids

        if lattice_id is None:
            global auto_universe_id
            self._id = auto_universe_id
            universe_ids.append(auto_universe_id)
            auto_universe_id += 1

        # Check that the ID is an integer and wasn't already used
        elif isinstance(lattice_id, (int, np.int32, np.int64)):

            # If the Lattice already has an ID, remove it from global list
            if not self._id is None:
                universe_ids.remove(self._id)

            if lattice_id in universe_ids:
                exit('Unable to set Lattice ID to %s since a Lattice '
                      'with this ID was already initialized.', str(lattice_id))

            if lattice_id < 0:
                exit('Unable to set Lattice ID to %d since it must be a '
                     'non-negative integer', lattice_id)

            else:
                self._id = lattice_id
                universe_ids.append(lattice_id)

        else:
            exit('Unable to set a non-integer Lattice ID %s', str(lattice_id))


    def setName(self, name):

        if not isinstance(name, str):
            exit('Unable to set name for Lattice ID=%d with a non-string '
                 'value %s', self._id, str(name))

        else:
            self._name = name


    def setType(self, type):

        if not isinstance(type, str):
            exit('Unable to set type for Lattice ID=%d with a non-string '
                 'value %s', self._id, str(type))

        if not type in ['rectangular', 'hexagonal']:
            exit('Unable to set type for Lattice ID=%d as %s since it is not '
                 'rectangular or hexagonal', self._id, type)

        self._type = type


    def setDimension(self, dimension):

        if not isinstance(dimension, (tuple, list, np.ndarray)):
            exit('Unable to set Lattice ID=%d dimension to %s since it is not '
                 'a Python tuple/list or NumPy array', self._id, str(dimension))

        if len(dimension) != 2 and len(dimension) != 3:
            exit('Unable to set Lattice ID=%d dimension to %s since it does '
                 'not contain 2 or 3 coordinates', self._id, str(dimension))

        for dim in dimension:

            if not isinstance(dim, (int, np.int32, np.int64)) \
              and not isinstance(dim, (float, np.float32, np.float64)):
                exit('Unable to set Lattice ID=%d dimension to %s since it '
                     'is not an integer or floating point value',
                     self._id, str(dim))

            if dim < 0:
                exit('Unable to set Lattice ID=%d dimension to %s since it '
                     'is a negative value', self._id, dim)

        self._dimension = dimension


    def setLowerLeft(self, lower_left):

        if not isinstance(lower_left, (tuple, list, np.ndarray)):
            exit('Unable to set Lattice ID=%d lower_left to %s since '
                 'it is not a Python tuple/list or NumPy array',
                 self._id, str(lower_left))

        if len(lower_left) != 2 and len(lower_left) != 3:
            exit('Unable to set Lattice ID=%d lower_left to %s since it does '
                 'not contain 2 or 3 coordinates', self._id, str(lower_left))

        for dim in lower_left:

            if not isinstance(dim, (int, np.int32, np.int64)) \
              and not isinstance(dim, (float, np.float32, np.float64)):
                exit('Unable to set Lattice ID=%d lower_left to %s since it '
                     'is not an integer or floating point value',
                     self._id, str(dim))

        self._lower_left = lower_left


    def setWidth(self, width):

        if not isinstance(width, (tuple, list, np.ndarray)):
            exit('Unable to set Lattice ID=%d width to %s since '
                 'it is not a Python tuple/list or NumPy array',
                 self._id, str(width))

        if len(width) != 2 and len(width) != 3:
            exit('Unable to set Lattice ID=%d width to %s since it does '
                 'not contain 2 or 3 coordinates', self._id, str(width))

        for dim in width:

            if not isinstance(dim, (int, np.int32, np.int64)) \
              and not isinstance(dim, (float, np.float32, np.float64)):
                exit('Unable to set Lattice ID=%d width to %s since it '
                     'is not an integer or floating point value',
                     self._id, str(dim))

            if dim < 0:
                exit('Unable to set Lattice ID=%d width to %s since it '
                     'is a negative value', self._id, dim)

        self._width = width


    def setOutside(self, outside):

        if not isinstance(outside, (int, np.int32, np.int64)):
            exit('Unable to set Lattice ID=%d outside universe to a '
                 'non-integer value %s', self._id, str(outside))

        if outside < 0:
            exit('Unable to set Lattice ID=%d outside universe to a '
                 'negative value %d', self._id, outside)

        self._outside = outside


    def setUniverses(self, universes):

        if not isinstance(universes, (tuple, list, np.ndarray)):
            exit('Unable to set Lattice ID=%d universes to %s since '
                 'it is not a Python tuple/list or NumPy array',
                 self._id, str(universes))

        self._universes = universes


    def toString(self):

        string = ''

        string += 'Lattice\n'

        cell_id = '{0: <16}'.format('\tID') + '=\t' + str(self._id)
        string += cell_id + '\n'

        name = '{0: <16}'.format('\tName') + '=\t' + self._name
        string += name + '\n'

        type = '{0: <16}'.format('\tType') + '=\t'
        type += str(self._type)
        string += type + '\n'

        dimension = '{0: <16}'.format('\tDimension') + '=\t'
        dimension += str(self._dimension)
        string += dimension + '\n'

        lower_left = '{0: <16}'.format('\tLower Left') + '=\t'
        lower_left += str(self._lower_left)
        string += lower_left + '\n'

        width = '{0: <16}'.format('\tWidth') + '=\t'
        width += str(self._width)
        string += width + '\n'

        outside = '{0: <16}'.format('\tOutside') + '=\t'
        outside += str(self._outside)
        string += outside + '\n'

        universes = '{0: <16}'.format('\tUniverses') + '\n'
        for i in range(len(self._universes)):
            universes += '\t'
            for j in range(len(self._universes[0])):
                universes += '%s ' % str(int(self._universes[i][j]))
            universes += '\n'
        string += universes.rstrip('\n')

        return string


    def printString(self):
        print(self.toString())