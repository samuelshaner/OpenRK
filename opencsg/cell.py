import numpy as np


# Lists of all IDs for all Cells and Universes created
cell_ids = list()
universe_ids = list()

# A static variable for auto-generated Cell and Universe IDs
auto_cell_id = 10000
auto_universe_id = 10000


class Cell(object):

    def __init__(self, cell_id=None, name='', universe=0,
                 fill=None, material=None, surfaces=None):

        # Initialize Cell class attributes
        self._id = None
        self._name = None
        self._universe = None
        self._fill = None
        self._material = None
        self._surfaces = list()
        self._rotation = None
        self._translation = None

        self.setId(cell_id)
        self.setName(name)
        self.setUniverse(universe)

        if not fill is None:
            self.setFill(fill)

        if not material is None:
            self.setMaterial(material)

        if not surfaces is None:
            self.setSurfaces(surfaces)

    def getId(self):
        return self._id

    def getName(self):
        return self._name

    def getUniverse(self):
        return self._universe

    def getFill(self):
        return self._fill

    def getMaterial(self):
        return self._material

    def getSurfaces(self):
        return self._surfaces

    def getRotation(self):
        return self._rotation

    def getTranslation(self):
        return self._translation

    def setId(self, cell_id=None):

        global cell_ids

        if cell_id is None:
            global auto_cell_id
            self._id = auto_cell_id
            cell_ids.append(auto_cell_id)
            auto_cell_id += 1

        # Check that the ID is an integer and wasn't already used
        elif isinstance(cell_id, (int, np.int32, np.int64)):

            # If the Cell already has an ID, remove it from global list
            if not self._id is None:
                cell_ids.remove(self._id)

            if cell_id in cell_ids:
                exit('Unable to set Cell ID to %s since a Material '
                      'with this ID was already initialized.', str(cell_id))

            if cell_id < 0:
                exit('Unable to set Cell ID to %d since it must be a '
                     'non-negative integer', cell_id)

            else:
                self._id = cell_id
                cell_ids.append(cell_id)

        else:
            exit('Unable to set a non-integer Cell ID %s', str(cell_id))


    def setName(self, name):

        if not isinstance(name, str):
            exit('Unable to set name for Cell ID=%d with a non-string '
                 'value %s', self._id, str(name))

        else:
            self._name = name


    def setUniverse(self, universe):

        global universe_ids

        if universe is None:
            global auto_universe_id
            self._universe = auto_universe_id
            auto_cell_id += 1

        elif not isinstance(universe, (int, np.int32, np.int64)):
            exit('Unable to set Cell ID=%d to use a non-integer universe %s',
                 self._id, str(universe))

        if universe < 0:
            exit('Unable to set Cell ID=%d to use a negative universe %d',
                 universe)

        else:
            self._universe = universe

        if not self._universe in universe_ids:
            universe_ids.append(auto_universe_id)


    def setFill(self, fill):

        if not isinstance(fill, (int, np.int32, np.int64)):
            exit('Unable to set Cell ID=%d to use a non-integer fill %s',
                 self._id, str(fill))

        if fill < 0:
            exit('Unable to set Cell ID=%d to use a negative fill %d',
                 fill)

        self._fill = fill


    def setMaterial(self, material):

        if not isinstance(material, (int, np.int32, np.int64)):
            exit('Unable to set Cell ID=%d to use a non-integer material '
                 '%s', self._id, str(material))

        if material < 0:
            exit('Unable to set Cell ID=%d to use a negative material %d',
                 material)

        self._material = material


    def addSurface(self, surface):

        if not isinstance(surface, (int, np.int32, np.int64)):
            exit('Unable to set Cell ID=%d to use a non-integer surface %s',
                  self._id, str(surface))

        if not surface in self._surfaces:
            self._surfaces.append(surface)


    def removeSurface(self, surface):

        if surface in self._surfaces:
            self._surfaces.remove(surface)


    def setSurfaces(self, surfaces):

        if not isinstance(surfaces, (tuple, list, np.ndarray)):
            exit('Unable to set Cell ID=%d with surfaces %s since it is not a '
                 'Python tuple/list or NumPy array', self._id, str(surfaces))

        for surface in surfaces:
            self.addSurface(surface)


    def setRotation(self, rotation):

        if not isinstance(rotation, (tuple, list, np.ndarray)):
            exit('Unable to add rotation %s to Cell ID=%d since it is not a '
                 'Python tuple/list or NumPy array', str(rotation), self._id)

        if len(rotation) != 3:
            exit('Unable to add rotation %s to Cell ID=%d since it does not '
                 'contain 3 values', str(rotation), self._id)

        for axis in rotation:
            if not isinstance(axis, (int, np.int32, np.int64)) \
              and not isinstance(axis, (float, np.float32, np.float64)):
                exit('Unable to add rotation %s to Cell ID=%d since it is '
                     'not an integer or floating point value',
                     str(axis), self._id)

        self._rotation = rotation


    def setTranslation(self, translation):

        if not isinstance(translation, (tuple, list, np.ndarray)):
            exit('Unable to add translation %s to Cell ID=%d since it is not a '
                 'Python tuple/list or NumPy array', str(translation), self._id)

        if len(translation) != 3:
            exit('Unable to add translation %s to Cell ID=%d since it does not '
                 'contain 3 values', str(translation), self._id)

        for axis in translation:
            if not isinstance(axis, (int, np.int32, np.int64)) \
              and not isinstance(axis, (float, np.float32, np.float64)):
                exit('Unable to add translation %s to Cell ID=%d since it is '
                     'not an integer or floating point value',
                     str(axis), self._id)

        self._translation = translation


    def toString(self):

        string = ''

        string += 'Cell\n'

        cell_id = '{0: <16}'.format('\tID') + '=\t' + str(self._id)
        string += cell_id + '\n'

        name = '{0: <16}'.format('\tName') + '=\t' + self._name
        string += name + '\n'

        universe = '{0: <16}'.format('\tUniverse') + '=\t'
        universe += str(self._universe)
        string += universe + '\n'

        fill = '{0: <16}'.format('\tFill') + '=\t'
        fill += str(self._fill)
        string += fill + '\n'

        material = '{0: <16}'.format('\tMaterial') + '=\t'
        material += str(self._material)
        string += material + '\n'

        surfaces = '{0: <16}'.format('\tSurfaces') + '=\t'
        surfaces += str(self._surfaces)
        string += surfaces + '\n'

        rotation = '{0: <16}'.format('\tRotation') + '=\t'
        rotation += str(self._rotation)
        string += rotation + '\n'

        translation = '{0: <16}'.format('\tTranslation') + '=\t'
        translation += str(self._translation)
        string += translation

        return string


    def printString(self):
        print(self.toString())