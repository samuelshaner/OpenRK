__author__ = 'Will Boyd'
__email__ = 'wboyd@mit.edu'

from surface import *
from cell import *
from lattice import *


class GeometryFile(object):

    def __init__(self):

        # Initialize Geometry class attributes
        self._geometry_file = ET.Element("geometry")
        self._surfaces = list()
        self._cells = list()
        self._lattices = list()


    def addSurface(self, surface):

        if not isinstance(surface, Surface):
            exit('Unable to add %s to GeometryFile since it is '
                 'not a Surface', str(surface))

        if not surface in self._cells:
            self._surfaces.append(surface)


    def removeSurface(self, surface):

        if surface in self._surfaces:
            self._surfaces.remove(surface)


    def addCell(self, cell):

        if not isinstance(cell, Cell):
            exit('Unable to add %s to GeometryFile since it is '
                 'not a Cell', str(cell))

        if not cell in self._cells:
            self._cells.append(cell)


    def removeCell(self, cell):

        if cell in self._cells:
            self._cells.remove(cell)


    def addLattice(self, lattice):

        if not isinstance(lattice, Lattice):
            exit('Unable to add %s to GeometryFile since it is '
                 'not a Lattice', str(lattice))

        if not lattice in self._lattices:
            self._lattices.append(lattice)


    def removeLattice(self, lattice):

        if lattice in self._lattices:
            self._lattices.remove(lattice)


    def toString(self):

      string = ''

      for cell in self._cells:
          string += cell.toString()

      for lattice in self._lattices:
          string += lattice.toString()

      return string


    def printString(self):
        print(self.toString())