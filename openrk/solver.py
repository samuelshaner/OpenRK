import numpy as np
from math import *
from cell import *
from material import *
from mesh import *
import openrk
from checkvalue import *

class Solver(object):

  def __init__(self):

    self._A = None
    self._AM = None
    self._M = None
    self._b = None

    openrk.setNumThreads(1)


  def setNumThreads(self, num_threads):

    # Check input values
    check_is_int(num_threads, 'Solver number of threads', 'num threads')

    self._num_threads = num_threads


class CmfdSolver(self, cmfd_mesh):

  def __init__(self):

    super(CmfdSolver, self).__init__()

    self._cmfd_mesh = cmfd_mesh
    self.setNumCMFDEnergyGroups(cmfd_mesh.getNumCmfdEnergyGroups())
    nc = cmfd_mesh.getCellsX() * cmfd_mesh.getCellsY()
    ng = self._num_energy_groups
    self._A   = np.zeros((nc, ng*(ng+4)))
    self._AM  = np.zeros((nc, ng*(ng+4)))
    self._M   = np.zeros((nc, ng*ng))
    self._phi = np.ones(nc*ng)
    self._b   = np.zeros(nc*ng)


  def computeSurfaceDifCoefs(self):

    cw = self._cmfd_mesh.getCellsX()
    ch = self._mesh.getCellsY()
    ng = self._num_energy_groups
    width = self._mesh.getCellWidth()
    height = self._mesh.getCellHeight()

    for y in xrange(ch):
      for x in xrange(cw):
        
        cell = self._mesh.getCell(x,y)
        
        for side in xrange(4):
          
          surface = cell.getSurface(side)
          cell_next = self._mesh.getNeighborCell(x, y, side)
          material = cell.getMaterial()
          
          for e in xrange(ng):
            
            d = material.getDifCoefByGroup(e)
            
            # set the length of the surface parallel to and perpendicular from surface
            if side == 0 or side == 2:
              length_perpen = width
            elif side == 1 or side == 3:
              length_perpen = height

            if cell_next is None:
              dif_hat = 2 * d / length_perpen / (1 + 4 * d / length_perpen)
              dif_hat *= self._mesh.getBoundary(side)

            else:

              if side == 0 or side == 2:
                next_length_perpen = width
              elif side == 1 or side == 3:
                next_length_perpen = height

              d_next = cell_next.getMaterial().getDifCoefByGroup(e)
              dif_hat = 2 * d * d_next / (length_perpen * d + next_length_perpen * d_next)

            surface.setDifCoefLinearByGroup(dif_hat, group)


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

          cell = self._mesh.getCell(x, y)
          material = cell.getMaterial()

          # absorption term on diagonal
          self._A[y*cw+x][e*(ng+4)+e+2] += material.getSigmaAByGroup(e) * cell.getVolume()

          # out scattering term on diagonal
          for g in xrange(ng):
            if e != g:
              self._A[y*cw+x][e*(ng+4)+e+2] += material.getSigmaSByGroup[e,g] * cell.getVolume()
              
          # fission terms on diagonal
          for g in xrange(ng):
            self._M[y*cw+x][e*ng+g] = material.getChiByGroup(e) * material.getNuSigmaFByGroup(g) * cell.getVolume()

          # in scattering terms on off diagonals
          for g in xrange(ng):
            if e != g:
              self._A[y*cw+x][e*(ng+4)+g+2] -= material.getSigmaSByGroup(g,e) * cell.getVolume()

          # RIGHT SURFACE
          surface = cell.getSurface(2)

          # transport term on diagonal
          self._A[y*cw+x][e*(ng+4)+e+2] += surface.getDifCoefLinearByGroup(e) * height
                        
          # transport terms on off diagonals
          if x != cw - 1:
            self._A[y*cw+x][e*(ng+4)+ng+2] -= surface.getDifCoefLinearByGroup(e) * height

          # LEFT SURFACE
          surface = cell.getSurface(0)
                    
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

    # pass flux to cells
    cw = self._mesh._cells_x
    ch = self._mesh._cells_y
    for y in xrange(ch): 
      for x in xrange(cw):
        cell = self._mesh._cells[y][x]
        for e in range(ng):          
          cell._flux[e] = self._phi[(y*cw+x)*ng + e]


class NonlinearCmfdSolver(self, cmfd_mesh):

  def __init__(self):

    super(NonlinearCmfdSolver, self).__init__()


  def computeSurfaceDifCoefs(self):

    cw = self._cmfd_mesh.getCellsX()
    ch = self._mesh.getCellsY()
    ng = self._num_energy_groups
    width = self._mesh.getCellWidth()
    height = self._mesh.getCellHeight()

    for y in xrange(ch):
      for x in xrange(cw):
        
        cell = self._mesh.getCell(x,y)
        
        for side in xrange(4):
          
          surface = cell.getSurface(side)
          cell_next = self._mesh.getNeighborCell(x, y, side)
          material = cell.getMaterial()
          
          for e in xrange(ng):
            
            d = material.getDifCoefByGroup(e)
            
            # set the length of the surface parallel to and perpendicular from surface
            if side == 0 or side == 2:
              length_perpen = width
            elif side == 1 or side == 3:
              length_perpen = height

            if cell_next is None:
              dif_hat = 2 * d / length_perpen / (1 + 4 * d / length_perpen)
              dif_hat *= self._mesh.getBoundary(side)

            else:

              if side == 0 or side == 2:
                next_length_perpen = width
              elif side == 1 or side == 3:
                next_length_perpen = height

              d_next = cell_next.getMaterial().getDifCoefByGroup(e)
              dif_hat = 2 * d * d_next / (length_perpen * d + next_length_perpen * d_next)

            surface.setDifCoefLinearByGroup(dif_hat, group)


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
          self._A[y*cw+x][e*(ng+4)+e+2] += cell._material._sigma_a[e] * cell.getVolume()

          # out scattering term on diagonal
          for g in xrange(ng):
            if e != g:
              self._A[y*cw+x][e*(ng+4)+e+2] += cell._material._sigma_s[e,g] * cell.getVolume()
              
          # fission terms on diagonal
          for g in xrange(ng):
            self._M[y*cw+x][e*ng+g] = cell._material._chi[e] * cell._material._nu_sigma_f[g] * cell.getVolume()

          # in scattering terms on off diagonals
          for g in xrange(ng):
            if e != g:
              self._A[y*cw+x][e*(ng+4)+g+2] -= cell._material._sigma_s[g,e] * cell.getVolume()

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






class TransientMOCSolver(Solver):

  def __init__(self, moc_mesh, tcmfd_mesh, cmfd_mesh=None):

    super(CmfdSolver, self).__init__()

    self._moc_mesh = moc_mesh
    self._cmfd_mesh = cmfd_mesh
    self._tcmfd_mesh = cmfd_mesh

    self._num_moc_energy_groups = mov_mesh._num_moc_energy_groups
    self._num_cmfd_energy_groups = cmfd_mesh._num_cmfd_energy_groups
    self._num_tcmfd_energy_groups = tcmfd_mesh._num_tcmfd_energy_groups

    self._A = None
    self._M = None
    self._AM = None
    self._b = None













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
          self._A[y*cw+x][e*(ng+4)+e+2] += cell._material._sigma_a[e] * cell.getVolume()

          # out scattering term on diagonal
          for g in xrange(ng):
            if e != g:
              self._A[y*cw+x][e*(ng+4)+e+2] += cell._material._sigma_s[e,g] * cell.getVolume()
              
          # fission terms on diagonal
          for g in xrange(ng):
            self._M[y*cw+x][e*ng+g] = cell._material._chi[e] * cell._material._nu_sigma_f[g] * cell.getVolume()

          # in scattering terms on off diagonals
          for g in xrange(ng):
            if e != g:
              self._A[y*cw+x][e*(ng+4)+g+2] -= cell._material._sigma_s[g,e] * cell.getVolume()

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

    # pass flux to cells
    cw = self._mesh._cells_x
    ch = self._mesh._cells_y
    for y in xrange(ch): 
      for x in xrange(cw):
        cell = self._mesh._cells[y][x]
        for e in range(ng):          
          cell._flux[e] = self._phi[(y*cw+x)*ng + e]

