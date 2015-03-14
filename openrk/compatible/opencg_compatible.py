__author__ = 'Samuel Shaner'
__email__ = 'shaner@mit.edu'

# Import modules
import opencg
import copy
import numpy as np
import openrk as rk


def get_openrk_amp_mesh(opencg_mesh):

  if not isinstance(opencg_mesh, opencg.Mesh):
    msg = 'Unable to get OpenRK amp mesh from {0} which is ' \
        'not an OpenCG Mesh object'.format(opencg_mesh)
    raise ValueError(msg)

  # Get parameters from OpenCG mesh
  width_x = opencg_mesh._width_x
  width_y = opencg_mesh._width_y
  width_z = opencg_mesh._width_z
  num_x = opencg_mesh._num_x
  num_y = opencg_mesh._num_y
  num_z = opencg_mesh._num_z

  # Create openmc mesh
  openrk_amp_mesh = rk.AmpMesh(width_x, width_y, width_z, num_x, num_y, num_z)
  openrk_amp_mesh.setShapeMap(opencg_mesh._region_ids)
  
  return openrk_amp_mesh


def get_openrk_shape_mesh(opencg_mesh):

  if not isinstance(opencg_mesh, opencg.Mesh):
    msg = 'Unable to get OpenRK shape mesh from {0} which is ' \
        'not an OpenCG Mesh object'.format(opencg_mesh)
    raise ValueError(msg)

  # Get parameters from OpenCG mesh
  width_x = opencg_mesh._width_x
  width_y = opencg_mesh._width_y
  width_z = opencg_mesh._width_z
  num_regions = opencg_mesh._num_regions

  # Create openmc mesh
  openrk_shape_mesh = rk.UnstructuredShapeMesh(width_x, width_y, width_z, num_regions)

  return openrk_amp_mesh
