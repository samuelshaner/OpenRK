import numpy as np
from math import *
from cell import *
from material import *
from mesh import *
import openrk

class Solver(object):

  def __init__(self, mesh):

    self._mesh = mesh
    self._num_energy_groups = mesh._cells[0][0]._material._num_energy_groups
    nc = mesh._cells_x * mesh._cells_y
    ng = self._num_energy_groups
    self._A = np.zeros((nc, ng*(ng+4)))
    self._M = np.zeros((nc, ng*ng))
    self._phi = np.ones(nc*ng)
    self._b = np.zeros(nc*ng)
    self._num_threads = 1

    self._k_eff = 1.0
    self._k_eff_old = 0.0

  def computeDs(self):

    cw = self._mesh._cells_x
    ch = self._mesh._cells_y
    ng = self._num_energy_groups
    width = self._mesh._cell_width
    height = self._mesh._cell_height

    for y in xrange(ch):
      for x in xrange(cw):
        
        cell = self._mesh._cells[y][x]
        
        for side in xrange(4):
          
          surface = cell._surfaces[side]
          cell_next = cell._neighbor_cells[side]
          
          for e in xrange(ng):
            
            d = cell._material._dif_coef[e]
            
            # set the length of the surface parallel to and perpendicular from surface
            if side == 0 or side == 2:
              length_perpen = width
            elif side == 1 or side == 3:
              length_perpen = height

            if cell_next is None:
              dif_hat = 2 * d / length_perpen / (1 + 4 * d / length_perpen)
              dif_hat *= self._mesh._boundaries[side]

            else:

              if side == 0 or side == 2:
                next_length_perpen = width
              elif side == 1 or side == 3:
                next_length_perpen = height

              d_next = cell_next._material._dif_coef[e]
              dif_hat = 2 * d * d_next / (length_perpen * d + next_length_perpen * d_next)

            surface._dif_hat[e] = dif_hat


  def makeAM(self):
      
    cw = self._mesh._cells_x
    ch = self._mesh._cells_y
    ng = self._num_energy_groups
    height = self._mesh._cell_height
    width = self._mesh._cell_width
    
    # reinitialize matrices to zero
    self._A.fill(0.0)
    self._M.fill(0.0)

    for y in xrange(ch): 
      for x in xrange(cw):

        for e in range(ng):

          cell = self._mesh._cells[y][x]

          # absorption term on diagonal
          self._A[y*cw+x][e*(ng+4)+e+2] += cell._material._sigma_a[e] * cell._volume

          # out scattering term on diagonal
          for g in xrange(ng):
            if e != g:
              self._A[y*cw+x][e*(ng+4)+e+2] += cell._material._sigma_s[e,g] * cell._volume
              
          # fission terms on diagonal
          for g in xrange(ng):
            self._M[y*cw+x][e*ng+g] = cell._material._chi[e] * cell._material._nu_sigma_f[g] * cell._volume

          # in scattering terms on off diagonals
          for g in xrange(ng):
            if e != g:
              self._A[y*cw+x][e*(ng+4)+g+2] -= cell._material._sigma_s[g,e] * cell._volume

          # RIGHT SURFACE
                    
          # transport term on diagonal
          self._A[y*cw+x][e*(ng+4)+e+2] += cell._surfaces[2]._dif_hat[e] * height
                        
          # transport terms on off diagonals
          if x != cw - 1:
            self._A[y*cw+x][e*(ng+4)+ng+2] -= cell._surfaces[2]._dif_hat[e] * height

          # LEFT SURFACE
                    
          # transport term on diagonal
          self._A[y*cw+x][e*(ng+4)+e+2] += cell._surfaces[0]._dif_hat[e] * height
                        
          # transport terms on off diagonals
          if x != 0:
            self._A[y*cw+x][e*(ng+4)] -= cell._surfaces[0]._dif_hat[e] * height

          # BOTTOM SURFACE

          # transport term on diagonal
          self._A[y*cw+x][e*(ng+4)+e+2] += cell._surfaces[1]._dif_hat[e] * width

          # transport terms on off diagonals
          if y != 0:
            self._A[y*cw+x][e*(ng+4)+1] -= cell._surfaces[1]._dif_hat[e] * width

          # TOP SURFACE

          # transport term on diagonal
          self._A[y*cw+x][e*(ng+4)+e+2] += cell._surfaces[3]._dif_hat[e] * width

          # transport terms on off diagonals
          if y != ch - 1:
            self._A[y*cw+x][e*(ng+4)+ng+3] -= cell._surfaces[3]._dif_hat[e] * width


  def solve(self, tol, iterations):

    self.computeDs()
    self.makeAM()
    nc = self._mesh._cells_x * self._mesh._cells_y
    ng = self._num_energy_groups
    old_source = np.zeros(nc*ng)
    flux_temp = np.zeros(nc*ng)
    self._k_eff = openrk.eigenvalueSolve(self._A, self._M, self._phi, self._b, old_source, flux_temp,
                                         ng, self._mesh._cells_x, self._mesh._cells_y, tol)

    print 'DIFFUSION: --- k_eff = ' + str(self._k_eff)[0:10]

  def setNumThreads(self, num_threads):

    self._num_threads = num_threads
    openrk.setNumThreads(num_threads)
