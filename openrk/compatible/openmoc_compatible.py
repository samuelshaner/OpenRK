__author__ = 'Samuel Shaner'
__email__ = 'shaner@mit.edu'

# Import modules
import openmoc
import copy
import numpy as np
import openrk as rk


def extract_openmoc_fsr_fluxes(tcmfd, moc_mesh, fluxes=['MOC new flux', 'MOC old flux']):

  for flux in fluxes:
    tcmfd.exportFSRScalarFluxes(moc_mesh._flux[flux])


def extract_openmoc_cmfd_fluxes(cmfd, cmfd_mesh, fluxes=['AMP new flux', 'AMP old flux']):

  if 'AMP new flux' in fluxes:
    cmfd.exportCmfdFluxNew(cmfd_mesh._flux['AMP new flux'])

  if 'AMP old flux' in fluxes:
    cmfd.exportCmfdFluxOld(cmfd_mesh._flux['AMP old flux'])


def extract_openrk_fsr_mesh(geometry):

  if not isinstance(geometry, openmoc.Geometry):
    msg = 'Unable to create an OpenRK Mesh from {0} ' \
          'which is not an OpenMOC Geometry'.format(openmoc_geometry)
    raise ValueError(msg)

  # Get general Geometry info (num_fsrs, width, height, num_moc_groups)
  num_fsrs = geometry.getNumFSRs()
  x_min = geometry.getXMin()
  x_max = geometry.getXMax()
  y_min = geometry.getYMin()
  y_max = geometry.getYMax()
  width = x_max - x_min
  height = y_max - y_min
  num_energy_groups = geometry.getNumEnergyGroups()
  fsr_mesh = rk.MOCMesh(name='MOC mesh', width=width, height=height)
  fsr_mesh.setXMin(x_min)
  fsr_mesh.setXMax(x_max)
  fsr_mesh.setYMin(y_min)
  fsr_mesh.setYMax(y_max)
  fsr_mesh.setNumCells(num_fsrs)
  fsr_mesh.setNumShapeEnergyGroups(num_energy_groups)

  # Initialize FSR mesh
  fsr_mesh.initializeFlux()
  fsr_mesh.initializeCells()

  return fsr_mesh
  

def create_openmoc_material(openrk_material, time=None, temp=None, clock_position='CURRENT'):

  # Check input value
  if not isinstance(openrk_material, rk.Material):
    msg = 'Unable to create an OpenMOC Material from {0} ' \
          'which is not an OpenRK Material'.format(openrk_material)
    raise ValueError(msg)

  # Create OpenMOC material
  openmoc_material = openmoc.TransientMaterial(openmoc.material_id())
  copy_openrk_material_to_openmoc(openrk_material, openmoc_material, time, temp, clock_position)

  return openmoc_material


def copy_openrk_material_to_openmoc(openrk_material, openmoc_material, time=None, temp=None, clock_position='CURRENT'):

  # Check input value
  if not isinstance(openrk_material, rk.Material):
    msg = 'Unable to copy to an OpenMOC Material from {0} ' \
          'which is not an OpenRK Material'.format(openrk_material)
    raise ValueError(msg)

  # Create OpenMOC material
  openmoc_material.setNumEnergyGroups(openrk_material._num_energy_groups)

  if isinstance(openrk_material, rk.FunctionalMaterial):

    # Check input require to create a transient material

    check_is_float_or_int(time, 'create openmoc material', 'time')
    check_is_float_or_int(temp, 'create openmoc material', 'temp')
    check_clock_position(clock_position, 'create openmoc material')

    for e in xrange(openrk_material._num_energy_groups):
      openmoc_material.setSigmaTByGroup(openrk_material.getSigmaTByGroup(e, time, temp), e)
      openmoc_material.setSigmaAByGroup(openrk_material.getSigmaAByGroup(e, time, temp), e)
      openmoc_material.setSigmaFByGroup(openrk_material.getSigmaFByGroup(e, time, temp), e)
      openmoc_material.setNuSigmaFByGroup(openrk_material.getNuSigmaFByGroup(e, time, temp), e)
      openmoc_material.setChiByGroup(openrk_material.getChiByGroup(e, time, temp), e)
      openmoc_material.setVelocityByGroup(openrk_material.getVelocityByGroup(e), e)

      for g in xrange(openrk_material._num_energy_groups):
        openmoc_material.setSigmaSByGroup(openrk_material.getSigmaSByGroup(e, g, time, temp), e, g)

    openmoc_material.setNumDelayedGroups(openrk_material._num_delayed_groups)
    for d in xrange(openrk_material._num_delayed_groups):
      openmoc_material.setDelayedFractionByGroup(openrk_material.getDelayedFractionByGroup(d), d)
      openmoc_material.setDecayConstantByGroup(openrk_material.getDecayConstantByGroup(d), d)    
      openmoc_material.setPrecursorConcByGroup(openrk_material.getPrecursorConcByGroup(d, clock_position), d)

  else:
    openmoc_material.setSigmaT(openrk_material._sigma_t)
    openmoc_material.setSigmaA(openrk_material._sigma_a)
    openmoc_material.setSigmaF(openrk_material._sigma_f)
    openmoc_material.setNuSigmaF(openrk_material._nu_sigma_f)
    openmoc_material.setSigmaS(openrk_material._sigma_s)
    openmoc_material.setChi(openrk_material._chi)


def extract_openrk_tcmfd_mesh(geometry):

  if not isinstance(geometry, openmoc.TransientGeometry):
    msg = 'Unable to create an OpenRK Mesh from {0} ' \
          'which is not an OpenMOC TransientGeometry'.format(openmoc_geometry)
    raise ValueError(msg)

  # Get general Geometry info (num_fsrs, width, height, num_moc_groups)
  tcmfd = geometry.getTcmfd()
  mesh = tcmfd.getMesh()
  width = mesh.getWidth()
  height = mesh.getHeight()
  num_x = mesh.getNumX()
  num_y = mesh.getNumY()
  num_moc_energy_groups = mesh.getNumMOCGroups()
  num_tcmfd_energy_groups = mesh.getNumMeshGroups()

  # Create and initialize CMFD Mesh
  tcmfd_mesh = rk.AmpMesh(name='Amp Mesh', width=width, height=height, num_x=num_x, num_y=num_y)
  tcmfd_mesh.setNumShapeEnergyGroups(num_moc_energy_groups)
  tcmfd_mesh.setNumAmpEnergyGroups(num_tcmfd_energy_groups)
  tcmfd_mesh.setXMin(-width/2.0)
  tcmfd_mesh.setXMax(width/2.0)
  tcmfd_mesh.setYMin(-height/2.0)
  tcmfd_mesh.setYMax(height/2.0)

  # Initialize tcmfd cells and flux arrays
  tcmfd_mesh.initializeCells()
  tcmfd_mesh.initializeFlux()

  # Extract fsr to tcmfd cells map
  fsrs_per_cell = np.zeros(num_x*num_y, dtype=np.int32)
  mesh.getCellToFSRIdsSizes(fsrs_per_cell)

  print 'fsrs per cell'
  print fsrs_per_cell

  for index,num_fsrs in enumerate(fsrs_per_cell):
    cell = tcmfd_mesh._cells[index/num_x][index%num_x]
    #cell.setNumFSRs(fsrs_per_cell[index])
    #geometry.getFSRsInTcmfdCell(cell.getFSRsArray(), index)

    #cmfd.getMOCToTCMFDGroups(cell._moc_to_tcmfd_groups, index)
    #tcmfd_mesh.generateAmpToShapeGroupsMap()

  return tcmfd_mesh


def copy_openrk_xs_to_openmoc(moc_mesh, geometry):

  # get current time
  current_time = moc_mesh._clock.getCurrentTime()
  
  # loop over FSRs / moc_mesh cells
  for fsr_id,cell in enumerate(moc_mesh._cells):
    
    # get openrk and openmoc materials
    openrk_material = cell._material
    openmoc_material = geometry.findFSRMaterial(fsr_id)
    temp = cell.getCurrentTemperature()
    copy_openrk_material_to_openmoc(openrk_material, openmoc_material, current_time, temp, 'CURRENT')

  
def extract_openrk_cmfd_mesh(geometry):

  if not isinstance(geometry, openmoc.TransientGeometry):
    msg = 'Unable to create an OpenRK Mesh from {0} ' \
          'which is not an OpenMOC TransientGeometry'.format(openmoc_geometry)
    raise ValueError(msg)

  # Get general Geometry info (num_fsrs, width, height, num_moc_groups)
  cmfd = geometry.getCmfd()
  mesh = cmfd.getMesh()
  width = mesh.getWidth()
  height = mesh.getHeight()
  num_x = mesh.getNumX()
  num_y = mesh.getNumY()
  num_moc_energy_groups = mesh.getNumMOCGroups()
  num_tcmfd_energy_groups = mesh.getNumMeshGroups()

  # Create and initialize CMFD Mesh
  cmfd_mesh = rk.AmpMesh(name='Amp Mesh', width=width, height=height, num_x=num_x, num_y=num_y)
  cmfd_mesh.setNumShapeEnergyGroups(num_moc_energy_groups)
  cmfd_mesh.setNumAmpEnergyGroups(num_tcmfd_energy_groups)
  cmfd_mesh.setXMin(-width/2.0)
  cmfd_mesh.setXMax(width/2.0)
  cmfd_mesh.setYMin(-height/2.0)
  cmfd_mesh.setYMax(height/2.0)

  cmfd_mesh.initializeCells()
  cmfd_mesh.initializeFlux()

  # Extract fsr to tcmfd cells map
  fsrs_per_cell = np.zeros(num_x*num_y, dtype=int)
  #geometry.getFSRsToTcmfdCellsArray(fsrs_per_cell)

  for index,num_fsrs in enumerate(fsrs_per_cell):
    cell = cmfd_mesh._cells[index/num_x][index%num_x]
    #cell.setNumFSRs(fsrs_per_cell[index])
    #geometry.getFSRsInTcmfdCell(cell.getFSRsArray(), index)

    #cmfd.getMOCToTCMFDGroups(cell._moc_to_tcmfd_groups, index)
    #tcmfd_mesh.generateAmpToShapeGroupsMap()

  return cmfd_mesh
