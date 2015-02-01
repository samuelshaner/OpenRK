__author__ = 'Samuel Shaner'
__email__ = 'shaner@mit.edu'

# Import modules
import numpy as np
import checkvalue as cv
from material import Material
from surface import Surface

# A static variable for auto-generated Material UIDs
AUTO_CELL_UID = 1


class Cell(object):
    def __init__(self):

        # Initialize class attributes
        global AUTO_CELL_UID
        self._id = AUTO_CELL_UID
        AUTO_CELL_UID += 1

        # Initialize class attributes
        self._material = None
        self._id = None

    def set_material(self, material):

        if not isinstance(material, Material):
            msg = 'Unable to set material for Cell ID={0} with a ' \
                  'non-material value {1}'.format(self._id, material)

            raise ValueError(msg)

        else:
            self._material = material

    def get_material(self):

        if self._material is None:
            msg = 'Unable to get material for Cell ID={0} since the ' \
                  'material has not been set'.format(self._id)
            raise ValueError(msg)

        else:
            return self._material


class CmfdCell(Cell):
    def __init__(self):

        # initialize FunctionalMaterial class attributes
        super(CmfdCell, self).__init__()

        # Initialize class attributes
        self._surfaces = np.empty(4, dtype=object)

    def get_surface(self, side):

        # Check input values
        cv.check_is_int(side, 'Cell ID={0} get surface'.format(self._id), 'side')

        if side < 0 or side > 3:
            msg = 'Cannot get surface for cell ID={0} for invalid side: {1}' \
                .format(self._id, side)
            raise ValueError(msg)

        elif self._surfaces[side] is None:
            msg = 'Cannot get surface for cell ID={0} since the surfaces have not been initialized!' \
                .format(self._id, side)
            raise ValueError(msg)

        else:

            return self._surfaces[side]
