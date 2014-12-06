
from math import *
import numpy as np
from Cell import *
from Surface import *


class Mesh(object):

  def __init__(self, width, height):
        
    # Initialize class attributes
    self._cells = None
    self._width = None
    self._height = None
    self._num_cells = None
    self._boundaries = np.zeros(4)
    
    # Set Mesh properties
    self.setWidth(width)
    self.setHeight(height)


  def setWidth(self, width):

    if not is_float(width):
      msg = 'Unable to set width to non-integer {0}'\
          .format(width)
      raise ValueError(msg)

    elif width <= 0.0:
      msg = 'Unable to set width to non-positive {0}'\
          .format(width)
      raise ValueError(msg)
        
    else:
      self._width = width


  def setHeight(self, height):

    if not is_float(height):
      msg = 'Unable to set height to non-integer {0}'\
          .format(height)
      raise ValueError(msg)

    elif height <= 0.0:
      msg = 'Unable to set height to non-positive {0}'\
          .format(height)
      raise ValueError(msg)
        
    else:
      self._height = height




class StructuredMesh(Mesh):

  def __init__(self, width, height, cells_x, cells_y):

    # initialize FunctionalMaterial class attributes
    super(Mesh, self).__init__(width, height)
    
    # Initialize class attributes
    self._cells_x = None
    self._cells_y = None
    
    # Set Mesh properties
    self.setNumCellsX(cells_x)
    self.setNumCellsY(cells_y)


  def initializeCells(self):
                
    self._cells = np.empty(shape=[self._cells_y, self._cells_x], dtype=object)
   
    # create cell objects
    for y in range(self._cells_y):
      for x in range(self._cells_x):
        self._cells[y][x] = Cell()
                                
    # set neighbor cells and allocate for each cell
    for y in xrange(self._cells_y):
      for x in xrange(self._cells_x):
         
        cell = self._cells[y][x]
         
        if x != 0:
          cell._neighbor_cells[0] = self._cells[y][x - 1]
        
        if y != 0:
          cell._neighbor_cells[1] = self._cells[y-1][x]
                    
        if x != self._cells_x - 1:
          cell._neighbor_cells[2] = self._cells[y][x + 1]

        if y != self._cells_y - 1:
          cell._neighbor_cells[3] = self._cells[y+1][x]

    for y in xrange(self._cells_y):
      for x in xrange(self.cells_x):
         
        cell = self.cells[y*self.cells_x+x]
        
        # Surface 0
        if x == 0:
          cell.surfaces[0] = Surface(cell.material.num_groups)
        else:
          cell.surfaces[0] = self._cells[y][x-1].surfaces[2]
                    
        # Surface 1
        if y == 0:
          cell.surfaces[1] = Surface(cell.material.num_groups)
        else:
          cell.surfaces[1] = self._cells[y-1][x].surfaces[3]
      
        # Surface 2 
        cell.surfaces[2] = Surface(cell.material.num_groups)

        # Surface 3
        cell.surfaces[3] = Surface(cell.material.num_groups)






  def setNumCellsX(self, cells_x):

    if not is_integer(cells_x):
      msg = 'Unable to set num cells x to non-integer {0}'\
          .format(cells_x)
      raise ValueError(msg)

    elif cells_x < 1:
      msg = 'Unable to set num cells x to non-positive {0}'\
          .format(cells_x)
      raise ValueError(msg)
        
    else:
      self._cells_x = cells_x

    if self._cells_y is not None:
      self.initializeCells()


  def setNumCellsY(self, cells_y):

    if not is_integer(cells_y):
      msg = 'Unable to set num cells y to non-integer {0}'\
          .format(cells_y)
      raise ValueError(msg)

    elif cells_y < 1:
      msg = 'Unable to set num cells x to non-positive {0}'\
          .format(cells_y)
      raise ValueError(msg)
        
    else:
      self._cells_y = cells_y

    if self._cells_x is not None:
      self.initializeCells()
