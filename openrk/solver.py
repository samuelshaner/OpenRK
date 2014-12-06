import numpy as np
import matplotlib.pyplot as plt
from math import *
from Cell import *
from Material import *
from Mesh import *
import os
import sys
import getopt
import time
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve


class Solver(object):

  def __init__(self, mesh):

    self._mesh = mesh
    self._num_energy_groups = mesh._cells[0]._material._num_energy_groups
    nc = mesh._cells_x * mesh._cells_y
    ng = self._num_energy_groups
    self._A = np.zeros((nc, ng*(ng+4)))
    self._M = np.zeros((nc, ng*ng))
    self._phi = np.ones(nc*ng)
    self._b = np.zeros(nc*ng)

    self._k_eff = 1.0
    self._k_eff_old = 0.0

  def computeDs(self):

    cw = self._mesh._cells_x
    ch = self._mesh._cells_y
    ng = self._num_energy_groups

    for y in xrange(ch):
      for x in xrange(cw):
        
        cell = self._mesh._cells[y][x]
        
        for side in xrange(4):
          
          surface = cell._surfaces[side]
          cell_next = cell._neighbor_cells[side]
          
          for e in xrange(ng):
            
            d = cell._material._dif_coef[e]
            
            # set the sense of the surface
            if side == 0 or side == 1:
              sense = -1.0
            else:
              sense = 1.0
              
            # set the length of the surface parallel to and perpendicular from surface
            if side == 0 or side == 2:
              length_perpen = cell.width
            elif side == 1 or side == 3:
              length_perpen = cell.height

            if cell_next is None:
              if surface._boundary == 'reflective':
                dif_hat = 0.0                    
              else:
                dif_hat = 2 * d / length_perpen / (1 + 4 * d / length_perpen)

            else:

              if side == 0 or side == 2:
                next_length_perpen = cellNext.width
              elif side == 1 or side == 3:
                next_length_perpen = cellNext.height

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

          cell = self.mesh.cells[y*cw+x]

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
            self._A[y*cw+x][e*(ng+4)+ng+2] -= cell._surfaces[2]._dif_coef[e] * height

          # LEFT SURFACE
                    
          # transport term on diagonal
          self._A[y*cw+x][e*(ng+4)+e+2] += cell._surfaces[0]._dif_coef[e] * .height
                        
          # transport terms on off diagonals
          if x != 0:
            self._A[y*cw+x][e*(ng+4)] -= cell._surfaces[0].D_dif_coefDif[e] * cell.height

                    # BOTTOM SURFACE
                    
                    # transport term on diagonal
                    if self.method == 'NEM4' or self.method == 'NEM2':
                        self.A[(y*cw+x)*ng + e, (y*cw+x)*ng + e] += (cell.surfaces[1].DHat[e] - cell.surfaces[1].DTilde[e]) * cell.width
                    elif self.method == 'diffusion':
                        self.A[(y*cw+x)*ng + e, (y*cw+x)*ng + e] += cell.surfaces[1].DDif[e] * cell.width
                        
                    # transport terms on off diagonals
                    if y != ch - 1:
                        if self.method == 'NEM4' or self.method == 'NEM2':
                            self.A[(y*cw+x)*ng + e, ((y+1)*cw+x)*ng + e] -= (cell.surfaces[1].DHat[e] + cell.surfaces[1].DTilde[e]) * cell.width
                        elif self.method == 'diffusion':
                            self.A[(y*cw+x)*ng + e, ((y+1)*cw+x)*ng + e] -= cell.surfaces[1].DDif[e] * cell.width

                    # TOP SURFACE
                    
                    # transport term on diagonal
                    if self.method == 'NEM4' or self.method == 'NEM2':
                        self.A[(y*cw+x)*ng + e, (y*cw+x)*ng + e] += (cell.surfaces[3].DHat[e] + cell.surfaces[3].DTilde[e]) * cell.width
                    elif self.method == 'diffusion':
                        self.A[(y*cw+x)*ng + e, (y*cw+x)*ng + e] += cell.surfaces[3].DDif[e] * cell.width
                        
                    # transport terms on off diagonals
                    if y != 0:
                        if self.method == 'NEM4' or self.method == 'NEM2':
                            self.A[(y*cw+x)*ng + e, ((y-1)*cw+x)*ng + e] -= (cell.surfaces[3].DHat[e] - cell.surfaces[3].DTilde[e]) * cell.width
                        elif self.method == 'diffusion':
                            self.A[(y*cw+x)*ng + e, ((y-1)*cw+x)*ng + e] -= cell.surfaces[3].DDif[e] * cell.width

        self.A = self.A.tocsr()
        self.M = self.M.tocsr()
                
           
    def makeN(self):
        
        cw = self.mesh.cells_x
        ch = self.mesh.cells_y
        ng = self.ng
        
        if self.method == 'NEM4':
            order = 4
            TL = np.zeros((3,3))
            TL[0,0] = 1.0
            TL[1,0] = 1.0
            TL[2,0] = 1.0
            TL[0,1] = -2.0
            TL[1,1] = 0.0
            TL[2,1] = 2.0
            TL[0,2] = -6.0
            TL[1,2] = 0.0
            TL[2,2] = -6.0
            f = np.zeros(3)
            current = np.zeros(3)
        elif self.method == 'NEM2':
            order = 2
        
        self.N = self.N.tolil()
        self.N = self.N * 0.0
        
        for y in range(ch): 
            for x in range(cw):
                
                cell = self.mesh.cells[y*cw+x]
                
                for e in range(ng):
                    
                    cn = 2*order*e + 2*order*ng*(y*cw+x)
                    cr = 2*order*e + 2*order*ng*(y*cw+x+1)
                    cl = 2*order*e + 2*order*ng*(y*cw+x-1)
                    cd = 2*order*e + 2*order*ng*((y+1)*cw+x)
                    cu = 2*order*e + 2*order*ng*((y-1)*cw+x)                    
                    
                    # SURFACE 0
                    if x == 0:
                        # CURRENT BOUNDARY
                        if order == 4:
                            self.N[cn,cn]   = self.dP1(0.0)
                            self.N[cn,cn+1] = self.dP2(0.0)
                            self.N[cn,cn+2] = self.dP3(0.0)
                            self.N[cn,cn+3] = self.dP4(0.0)
                        else:
                            self.N[cn,cn]   = self.dP1(0.0)
                            self.N[cn,cn+1] = self.dP2(0.0)
                    else:
                        # FLUX INTERFACE
                        self.N[cn,cn]   = - self.P1(0.0)
                        self.N[cn,cn+1] = - self.P2(0.0)
                        self.N[cn,cl]   = self.P1(1.0)
                        self.N[cn,cl+1] = self.P2(1.0)    
                        
                        self.b[cn] = - self.mesh.cells[y*cw+x-1].flux[e] + cell.flux[e]                        

                    # SURFACE 1
                    if y == ch - 1:
                        # CURRENT BOUNDARY
                        if order == 4:
                            self.N[cn+1,cn+4] = self.dP1(1.0)
                            self.N[cn+1,cn+5] = self.dP2(1.0)
                            self.N[cn+1,cn+6] = self.dP3(1.0)
                            self.N[cn+1,cn+7] = self.dP4(1.0)
                        else:
                            self.N[cn+1,cn+2] = self.dP1(1.0)
                            self.N[cn+1,cn+3] = self.dP2(1.0)                        
                    else:
                        cell_next = self.mesh.cells[(y+1)*cw+x]

                        # CURRENT INTERFACE                        
                        if order == 4:
                            self.N[cn+1,cn+4] = - self.dP1(1.0) * cell.material.D[e] / cell.height
                            self.N[cn+1,cn+5] = - self.dP2(1.0) * cell.material.D[e] / cell.height
                            self.N[cn+1,cn+6] = - self.dP3(1.0) * cell.material.D[e] / cell.height
                            self.N[cn+1,cn+7] = - self.dP4(1.0) * cell.material.D[e] / cell.height
                            self.N[cn+1,cd+4] = self.dP1(0.0) * cell_next.material.D[e] / cell_next.height
                            self.N[cn+1,cd+5] = self.dP2(0.0) * cell_next.material.D[e] / cell_next.height
                            self.N[cn+1,cd+6] = self.dP3(0.0) * cell_next.material.D[e] / cell_next.height
                            self.N[cn+1,cd+7] = self.dP4(0.0) * cell_next.material.D[e] / cell_next.height                            
                        else:
                            self.N[cn+1,cn+2] = - self.dP1(1.0) * cell.material.D[e] / cell.height
                            self.N[cn+1,cn+3] = - self.dP2(1.0) * cell.material.D[e] / cell.height
                            self.N[cn+1,cd+2] = self.dP1(0.0) * cell_next.material.D[e] / cell_next.height
                            self.N[cn+1,cd+3] = self.dP2(0.0) * cell_next.material.D[e] / cell_next.height
                        
                    # SURFACE 2
                    if x == cw - 1:
                        # CURRENT BOUNDARY 
                        if order == 4:
                            self.N[cn+2,cn]   = self.dP1(1.0)
                            self.N[cn+2,cn+1] = self.dP2(1.0)
                            self.N[cn+2,cn+2] = self.dP3(1.0)
                            self.N[cn+2,cn+3] = self.dP4(1.0)
                        else:
                            self.N[cn+2,cn]   = self.dP1(1.0)
                            self.N[cn+2,cn+1] = self.dP2(1.0)                        
                    else:  
                        # CURRENT INTERFACE
                        cell_next = self.mesh.cells[y*cw+x+1]
                        
                        if order == 4:
                            self.N[cn+2,cn]   = - self.dP1(1.0) * cell.material.D[e] / cell.width
                            self.N[cn+2,cn+1] = - self.dP2(1.0) * cell.material.D[e] / cell.width
                            self.N[cn+2,cn+2] = - self.dP3(1.0) * cell.material.D[e] / cell.width
                            self.N[cn+2,cn+3] = - self.dP4(1.0) * cell.material.D[e] / cell.width
                            self.N[cn+2,cr]   = self.dP1(0.0) * cell_next.material.D[e] / cell_next.width
                            self.N[cn+2,cr+1] = self.dP2(0.0) * cell_next.material.D[e] / cell_next.width
                            self.N[cn+2,cr+2] = self.dP3(0.0) * cell_next.material.D[e] / cell_next.width
                            self.N[cn+2,cr+3] = self.dP4(0.0) * cell_next.material.D[e] / cell_next.width                            
                        else:
                            self.N[cn+2,cn]   = - self.dP1(1.0) * cell.material.D[e] / cell.width
                            self.N[cn+2,cn+1] = - self.dP2(1.0) * cell.material.D[e] / cell.width
                            self.N[cn+2,cr]   = self.dP1(0.0) * cell_next.material.D[e] / cell_next.width
                            self.N[cn+2,cr+1] = self.dP2(0.0) * cell_next.material.D[e] / cell_next.width                        
                    
                    # SURFACE 3
                    if y == 0:
                        # CURRENT BOUNDARY
                        if order == 4:
                            self.N[cn+3,cn+4] = self.dP1(0.0)
                            self.N[cn+3,cn+5] = self.dP2(0.0)
                            self.N[cn+3,cn+6] = self.dP3(0.0)
                            self.N[cn+3,cn+7] = self.dP4(0.0)
                        else:
                            self.N[cn+3,cn+2] = self.dP1(0.0)
                            self.N[cn+3,cn+3] = self.dP2(0.0)
                      
                    else:
                        # FLUX INTERFACE
                        if order == 4:
                            self.N[cn+3,cn+4] = - self.P1(0.0)
                            self.N[cn+3,cn+5] = - self.P2(0.0)
                            self.N[cn+3,cu+4] = self.P1(1.0)
                            self.N[cn+3,cu+5] = self.P2(1.0)  
                        else:
                            self.N[cn+3,cn+2] = - self.P1(0.0)
                            self.N[cn+3,cn+3] = - self.P2(0.0)
                            self.N[cn+3,cu+2] = self.P1(1.0)
                            self.N[cn+3,cu+3] = self.P2(1.0)    
                            
                        self.b[cn+3] = - self.mesh.cells[(y-1)*cw+x].flux[e] + cell.flux[e]                      
                        
                        
                    # WEIGHTED RESIDUALS
                    if order == 4:
                        
                        # RESIDUAL X1
                        self.N[cn+4,cn]   = cell.material.sigma_r[e] / 3.0
                        self.N[cn+4,cn+1] = 0.0
                        self.N[cn+4,cn+2] = cell.material.sigma_r[e] / 5.0 + 12 * cell.material.D[e] / cell.width**2
                        self.N[cn+4,cn+3] = 0.0
                        
                        # transverse leakage
                        f[:] = 0.0
                        current[:] = 0.0
                        current[1] = cell.surfaces[1].current[e] - cell.surfaces[3].current[e]
                        if x != 0:
                            current[0] = self.mesh.cells[y*cw+x - 1].surfaces[1].current[e] - self.mesh.cells[y*cw+x - 1].surfaces[3].current[e]
                        if x != cw - 1:
                            current[2] = self.mesh.cells[y*cw+x + 1].surfaces[1].current[e] - self.mesh.cells[y*cw+x + 1].surfaces[3].current[e]

                        current = current * self.relax
                        f = np.linalg.solve(TL, current)
#                         self.b[cn+4] = -f[1]/3.0
                        self.b[cn+4] = 0.25/3.0*(current[0] + current[2])

                        cnp = 2*order*ng*(y*cw+x)
                            
                        # in scattering and fission
                        for g in range(ng):
                            # fission
                            self.N[cn+4,cnp+2*order*g]   -= cell.material.chi[e] * cell.material.nu_sigma_f[g] / self.keff / 3.0
                            self.N[cn+4,cnp+2*order*g+2] -= cell.material.chi[e] * cell.material.nu_sigma_f[g] / self.keff / 5.0
                                
                            # in scattering
                            if g != e:
                                self.N[cn+4,cnp+2*order*g]   -= cell.material.sigma_s[g,e] / 3.0
                                self.N[cn+4,cnp+2*order*g+2] -= cell.material.sigma_s[g,e] / 5.0

                            
                        # RESIDUAL X2
                        self.N[cn+5,cn]   = 0.0
                        self.N[cn+5,cn+1] = cell.material.sigma_r[e] / 5.0
                        self.N[cn+5,cn+2] = 0.0
                        self.N[cn+5,cn+3] = - 3.0 * cell.material.sigma_r[e] / 35.0 - 12.0 * cell.material.D[e] / cell.width**2
 
                        # transverse leakage
#                         self.b[cn+5] = -f[2]/5.0
                        self.b[cn+5] = 1.0/60.0*(current[0] - 2*current[1] + current[2])
                            
                        # in scattering and fission
                        for g in range(ng):
                            # fission
                            self.N[cn+5,cnp+2*order*g+1] -= cell.material.chi[e] * cell.material.nu_sigma_f[g] / self.keff / 5.0
                            self.N[cn+5,cnp+2*order*g+3] -= cell.material.chi[e] * cell.material.nu_sigma_f[g] / self.keff * (-3.0) / 35.0
                                
                            # in scattering
                            if g != e:
                                self.N[cn+5,cnp+2*order*g+1] -= cell.material.sigma_s[g,e] / 5.0
                                self.N[cn+5,cnp+2*order*g+3] -= cell.material.sigma_s[g,e] * (-3.0) / 35.0

                        # RESIDUAL Y1
                        self.N[cn+6,cn+4] = cell.material.sigma_r[e] / 3.0
                        self.N[cn+6,cn+5] = 0.0
                        self.N[cn+6,cn+6] = cell.material.sigma_r[e] / 5.0 + 12 * cell.material.D[e] / cell.height**2
                        self.N[cn+6,cn+7] = 0.0
                        
                        # transverse leakage
                        f[:] = 0.0
                        current[:] = 0.0
                        current[1] = cell.surfaces[2].current[e] - cell.surfaces[0].current[e]
                        if y != 0:
                            current[0] = self.mesh.cells[(y-1)*cw+x].surfaces[2].current[e] - self.mesh.cells[(y-1)*cw+x].surfaces[0].current[e]
                        if y != ch - 1:
                            current[2] = self.mesh.cells[(y+1)*cw+x].surfaces[2].current[e] - self.mesh.cells[(y+1)*cw+x].surfaces[0].current[e]

                        current = current * self.relax
                        f = np.linalg.solve(TL, current)
#                         self.b[cn+6] = -f[1]/3.0
                        self.b[cn+6] = 0.25/3.0*(current[0] + current[2])

                        # in scattering and fission
                        for g in range(ng):
                            # fission
                            self.N[cn+6,cnp+2*order*g+4] -= cell.material.chi[e] * cell.material.nu_sigma_f[g] / self.keff / 3.0
                            self.N[cn+6,cnp+2*order*g+6] -= cell.material.chi[e] * cell.material.nu_sigma_f[g] / self.keff / 5.0
                                
                            # in scattering
                            if g != e:
                                self.N[cn+6,cnp+2*order*g+4] -= cell.material.sigma_s[g,e] / 3.0
                                self.N[cn+6,cnp+2*order*g+6] -= cell.material.sigma_s[g,e] / 5.0
                          
                        # RESIDUAL Y2
                        self.N[cn+7,cn+4] = 0.0
                        self.N[cn+7,cn+5] = cell.material.sigma_r[e] / 5.0
                        self.N[cn+7,cn+6] = 0.0
                        self.N[cn+7,cn+7] = - 3.0 * cell.material.sigma_r[e] / 35.0 - 12.0 * cell.material.D[e] / cell.height**2
                            
                        # transverse leakage
#                         self.b[cn+7] = -f[2]/5.0
                        self.b[cn+7] = 1.0/60.0*(current[0] - 2*current[1] + current[2])
                            
                        # in scattering and fission
                        for g in range(ng):
                            # fission
                            self.N[cn+7,cnp+2*order*g+5] -= cell.material.chi[e] * cell.material.nu_sigma_f[g] / self.keff / 5.0
                            self.N[cn+7,cnp+2*order*g+7] -= cell.material.chi[e] * cell.material.nu_sigma_f[g] / self.keff * (-3.0) / 35.0
                                
                            # in scattering
                            if g != e:
                                self.N[cn+7,cnp+2*order*g+5] -= cell.material.sigma_s[g,e] / 5.0
                                self.N[cn+7,cnp+2*order*g+7] -= cell.material.sigma_s[g,e] * (-3.0) / 35.0
                        
                        
        self.N = self.N.tocsr()

    
    def computeTransLeak(self):
                 
        cw = self.mesh.cells_x
        ch = self.mesh.cells_y
        ng = self.ng
                 
        for y in range(ch): 
            for x in range(cw):
                 
                cell = self.mesh.cells[y*cw+x]
                 
                for e in range(ng):
                     
                    # leakage in x and y directions
                    cell.leak[e]    = cell.surfaces[2].current[e] - cell.surfaces[0].current[e]
                    cell.leak[ng+e] = cell.surfaces[1].current[e] - cell.surfaces[3].current[e]                    
    
        
    def computeCoeffs(self):
        
        self.coeffs = spsolve(self.N,self.b)
        self.N = self.N.tolil()
        
    # compute the currents        
    def computeCurrents(self):
        
        cw = self.mesh.cells_x
        ch = self.mesh.cells_y
        ng = self.ng
        
        if self.method == 'NEM4':
            order = 4
        elif self.method == 'NEM2':
            order = 2
        
        for y in range(ch): 
            for x in range(cw):
                
                cell = self.mesh.cells[y*cw+x]
                
                for e in range(ng):
                    
                    for side in range(4):
                        surface = cell.surfaces[side]
                        
                        if cell.neighborCells[side] is not None:
                            cn = 2*order*ng*(y*cw+x) + 2*order*e
                            
                            if order == 4:
                                if side == 0:
                                    surface.current[e] = - cell.material.D[e] / cell.width  * (self.coeffs[cn]   * self.dP1(0.0) + self.coeffs[cn+1] * self.dP2(0.0) + self.coeffs[cn+2] * self.dP3(0.0) + self.coeffs[cn+3] * self.dP4(0.0))
                                elif side == 1:
                                    surface.current[e] = - cell.material.D[e] / cell.height * (self.coeffs[cn+4] * self.dP1(1.0) + self.coeffs[cn+5] * self.dP2(1.0) + self.coeffs[cn+6] * self.dP3(1.0) + self.coeffs[cn+7] * self.dP4(1.0))
                                elif side == 2:
                                    surface.current[e] = - cell.material.D[e] / cell.width  * (self.coeffs[cn]   * self.dP1(1.0) + self.coeffs[cn+1] * self.dP2(1.0) + self.coeffs[cn+2] * self.dP3(1.0) + self.coeffs[cn+3] * self.dP4(1.0))
                                else:
                                    surface.current[e] = - cell.material.D[e] / cell.height * (self.coeffs[cn+4] * self.dP1(0.0) + self.coeffs[cn+5] * self.dP2(0.0) + self.coeffs[cn+6] * self.dP3(0.0) + self.coeffs[cn+7] * self.dP4(0.0))
                                                                        
                            elif order == 2:
                                if side == 0:
                                    surface.current[e] = - cell.material.D[e] / cell.width  * (self.coeffs[cn]   * self.dP1(0.0) + self.coeffs[cn+1] * self.dP2(0.0))
                                elif side == 1:
                                    surface.current[e] = - cell.material.D[e] / cell.height * (self.coeffs[cn+2] * self.dP1(1.0) + self.coeffs[cn+3] * self.dP2(1.0))
                                elif side == 2:
                                    surface.current[e] = - cell.material.D[e] / cell.width * (self.coeffs[cn]    * self.dP1(1.0) + self.coeffs[cn+1] * self.dP2(1.0))
                                else:
                                    surface.current[e] = - cell.material.D[e] / cell.height * (self.coeffs[cn+2] * self.dP1(0.0) + self.coeffs[cn+3] * self.dP2(0.0))
                                    
    def P1(self, val):
        return (2 * val - 1)
    
    def dP1(self, val):
        return 2
    
    def P2(self, val):
        return (6.0*val*(1-val)-1)
    
    def dP2(self, val):
        return (6 - 12.0*val)
    
    def P3(self, val):
        return (6.0*val*(1-val)*(2*val-1))

    def dP3(self, val):
        return (-6.0*(1-6*val+6*val**2))
    
    def P4(self, val):
        return (6.0*val*(1-val)*(5*val**2-5*val+1))
    
    def dP4(self, val):
        return (6.0-72*val+180*val**2-120*val**3)     
        

    def solve(self, tol, iterations):
        
        if self.mesh.num_cells > 10:
            self.relax = min(1.0/self.mesh.cells[0].width, 0.66)
        else:
            self.relax = 0.66
            
        print 'relaxation factor - ' + str(self.relax)[0:6]

        
        if self.method == 'NEM4' or self.method == 'NEM2':
            for iteration in range(iterations):
                print 'CMFD outer iteration ' + str(iteration) + ' --- k_eff = ' + str(self.keff)[0:10]
                self.computeDs()
                self.makeAM()
                self.computeFlux(tol)
                self.makeN()
                self.computeTransLeak()
                self.computeCoeffs()
                self.computeCurrents()
        
                if abs(self.keff_old - self.keff) < tol:
                    print self.method + ': Converged in ' + str(iteration+1) + ' iterations --- k_eff = ' + str(self.keff)[0:10]
                    break
                                
        elif self.method == 'diffusion':
            self.computeDs()
            self.makeAM()
            self.computeFlux(tol)
            print 'DIFFUSION: --- k_eff = ' + str(self.keff)[0:10]


    def computeFlux(self, tol):
                 
        max_iter = 1000
        ng = self.mesh.cells[0].material.num_groups
        cw = self.mesh.cells_x
        ch = self.mesh.cells_y
        sold = np.zeros(cw*ch*ng)
        snew = np.zeros(cw*ch*ng)
        res = np.zeros(cw*ch*ng)
        self.phi = np.ones(cw*ch*ng)
        self.keff_old = self.keff

        # get initial source and find initial keff
        snew = self.M * self.phi
        sumnew = np.sum(snew)
        sold = self.A * self.phi
        sumold = np.sum(sold)
        self.keff = sumnew / sumold
        
        # recompute and initialize the initial source
        sold = self.M * self.phi
        sumold = sum(sold)
        sold = sold * (cw*ch*ng / sumold)
        sumold = cw*ch*ng
        
        for i in range(max_iter):
            self.phi = spsolve(self.A, sold)
            snew = self.M * self.phi
            sumnew = np.sum(snew)
            
            # compute and set keff
            self.keff = sumnew / sumold
            sold = sold * self.keff
             
            # compute L2 norm of source error
            snew += 1.e-10
            res = np.divide(sold, snew)
            res = res - 1.0
            l2_norm = np.linalg.norm(res)
         
            # compute error
            l2_norm = l2_norm / (cw*ch*ng)
            snew = snew * (cw*ch*ng) / sumnew
             
            sold = snew
            
#             print 'iteration: ' + str(i) + '  --- keff: ' + str(self.keff)[0:10]
             
            if l2_norm < tol:
                
                for i in range(cw*ch):
                    for e in range(ng):
                        self.mesh.cells[i].flux[e] = self.phi[i*ng+e]
                break            


def main():
    
    start = time.time()

    # set default values
    tol = 1.e-8
    cell_size    = 10.0
    solve_method = 'NEM4'
    iterations = 100

    # create mesh
    mesh = Mesh([cell_size,cell_size], [cell_size])
    
    # create fuel
    fuel = Material(2, 'fuel')
    fuel.setSigmaA([0.005, 0.10])
    fuel.setD([1.5, 0.40])
    fuel.setNuSigmaF([0.005, 0.15])
    fuel.setChi([1.0, 0.0])
    fuel.setSigmaS(np.array([[0.0, 0.02],[0.0, 0.0]]))
    
    # create fuel
    moderator = Material(2, 'moderator')
    moderator.setSigmaA([0.0, 0.01])
    moderator.setD([1.5, 0.20])
    moderator.setSigmaS(np.array([[0.0, 0.025],[0.0, 0.0]]))
    
    if solve_method == 'NEM4':
        order = 4
    else:
        order = 2
  
    # add materials to cells
    mesh.cells[0].setMaterial(fuel, order)
    mesh.cells[1].setMaterial(moderator, order)
#     mesh = mesh.refineMesh(.1)
    mesh.makeSurfaces()
    
    # plot the mesh
    pttr.plotMesh(mesh)
    
    # create solver
    solver = Solver(mesh, solve_method)   

    # solve the matrix problem to get flux profile and keff
    solver.solve(tol, iterations)
         
    # plot the flux    
    pttr.plotFlux(solver)
    pttr.plotCellFlux(solver)
    pttr.plotCurrent(solver)

    stop = time.time()
    
    print 'Ran time ' + str(stop-start)[0:5] + ' seconds'

    print '----------------------------------------------------------------------'

if __name__ == '__main__':

    main()



