__author__ = 'Samuel Shaner'
__email__ = 'shaner@mit.edu'

# Import modules
import openmc
import copy
import numpy as np
import openrk as rk


def get_openrk_currents(tally, dimension, num_groups):

  num_x = dimension[0]
  num_y = dimension[1]
  num_z = dimension[2]

  rk_currents = np.zeros(num_x*num_y*num_z*num_groups*6)

  # Get the array of currents
  tally.compute_std_dev()
  mc_currents = tally._mean

  # Compute the net current across each surface
  for z in xrange(num_z+1):
    for y in xrange(num_y+1):
      for x in xrange(num_x+1):
        rk_x = x-1 
        rk_y = y-1 
        rk_z = z-1 
        rk_cell = (rk_z)*num_x*num_y + (rk_y)*num_x + (rk_x)
        mc_cell = x*(num_y+1)*(num_z+1) + y*(num_z+1) + z
        
        for g in xrange(num_groups):
          
          # RIGHT SURFACE
          current = -mc_currents[mc_cell*6*num_groups + g*6][0][0] + mc_currents[mc_cell*6*num_groups + g*6 + 1][0][0]
          rk_cell_next = (z-1)*num_x*num_y + (y-1)*num_x + (x)
          if rk_x >= 0 and rk_x < num_x and rk_y >= 0 and rk_y < num_y and rk_z >= 0 and rk_z < num_z:
            rk_currents[(rk_cell*6 + 3)*num_groups + g] = current
          if rk_x >= -1 and rk_x < num_x-1 and rk_y >= 0 and rk_y < num_y and rk_z >= 0 and rk_z < num_z:
            rk_currents[(rk_cell_next*6 + 0)*num_groups + g] = current

          # FRONT SURFACE
          current = -mc_currents[mc_cell*6*num_groups + g*6 + 2][0][0] + mc_currents[mc_cell*6*num_groups + g*6 + 3][0][0]
          rk_cell_next = (z-1)*num_x*num_y + (y)*num_x + (x-1)
          if rk_x >= 0 and rk_x < num_x and rk_y >= 0 and rk_y < num_y and rk_z >= 0 and rk_z < num_z:
            rk_currents[(rk_cell*6 + 4)*num_groups + g] = current
          if rk_x >= 0 and rk_x < num_x and rk_y >= -1 and rk_y < num_y-1 and rk_z >= 0 and rk_z < num_z:
            rk_currents[(rk_cell_next*6 + 1)*num_groups + g] = current

          # TOP SURFACE
          current = -mc_currents[mc_cell*6*num_groups + g*6 + 4][0][0] + mc_currents[mc_cell*6*num_groups + g*6 + 5][0][0]
          rk_cell_next = (z)*num_x*num_y + (y-1)*num_x + (x-1)
          if rk_x >= 0 and rk_x < num_x and rk_y >= 0 and rk_y < num_y and rk_z >= 0 and rk_z < num_z:
            rk_currents[(rk_cell*6 + 5)*num_groups + g] = current
          if rk_x >= 0 and rk_x < num_x and rk_y >= 0 and rk_y < num_y and rk_z >= -1 and rk_z < num_z-1:
            rk_currents[(rk_cell_next*6 + 2)*num_groups + g] = current
  

  return rk_currents
