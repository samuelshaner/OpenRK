#include "UnstructuredShapeMesh.h"


/**
 * @brief Constructor sets the ID and unique ID for the Material.
 * @param id the user-defined ID for the material
 */
UnstructuredShapeMesh::UnstructuredShapeMesh(double width, double height, double depth, int num_cells) : ShapeMesh(width, height, depth){

  _volume = NULL;
  _mesh_type = UNSTRUCTURED_SHAPE_MESH;
  
  return;
}


/**
 * @brief Destructor deletes all cross-section data structures from memory.
 */
UnstructuredShapeMesh::~UnstructuredShapeMesh() {

  if (_volume != NULL)
    delete [] _volume;
}


UntructuredShapeMesh* UnstructuredShapeMesh::clone(){

  if (_num_shape_energy_groups == 0)
    log_printf(ERROR, "Unable to clone the UnstructuredShapeMesh since the number "
               "of shape energy groups has not been set");

  if (_num_amp_energy_groups == 0)
    log_printf(ERROR, "Unable to clone the UnstructuredShapeMesh since the number "
               "of amp energy groups has not been set");

  if (_num_delayed_groups == 0)
    log_printf(ERROR, "Unable to clone the UnstructuredShapeMesh since the number "
               "of delayed groups has not been set");
  
  UnstructuredShapeMesh* mesh = new UnstructuredShapeMesh(getWidth(), getHeight(), getDepth(), 
                                                          _num_cells);

  mesh->setNumShapeEnergyGroups(_num_shape_energy_groups);
  mesh->setNumAmpEnergyGroups(_num_amp_energy_groups);
  mesh->setNumDelayedGroups(_num_delayed_groups);

  if (_k_eff_0 != 0.0)
    mesh->setKeff0(_k_eff_0);
 
  for (int s=0; s < 6; s++)
    mesh->setBoundary(s, _boundaries[s]);

  if (_clock != NULL)
    mesh->setClock(_clock);

  if (_decay_constants != NULL)
    mesh->setDecayConstants(_decay_constants, _num_delayed_groups);

  if (_delayed_fractions != NULL)
    mesh->setDelayedFractions(_delayed_fractions, _num_delayed_groups);

  return mesh;  
}


void UnstructuredShapeMesh::setVolume(double* volume, int num_cells){

  if (num_cells != _num_cells)
    log_printf(ERROR, "Unable to set the volume for an Unstructured Shape Mesh "
               "with an array of length %i that is not equal to the number of"
               " cells: %i", num_cells, _num_cells);

  if (_volume != NULL)
    delete [] _volume;

  _volume = new double[_num_cells];

  std::copy(volume, volume + _num_cells, _volumes);
}


void UnstructuredShapeMesh::setTemperature(double* temperature, int num_cells){

  if (num_cells != _num_cells)
    log_printf(ERROR, "Unable to set the temperature for an Unstructured Shape "
               "Mesh with an array of length %i that is not equal to the number of"
               " cells: %i", num_cells, _num_cells);

  for (int c=0; c < 8; c++)    
    std::copy(temperature, temperature + _num_cells, _temperature[c]);
}


void UnstructuredShapeMesh::getNumCells(){
  return _num_cells;
}


double UnstructuredShapeMesh::getCellVolume(int cell){

  if (_volume == NULL)
    log_printf(ERROR, "Unable to get volume for an Unstructured Shape Mesh since "
               "the mesh volumes have not been set");

  if (cell >= _num_cells || cell < 0)
    log_printf(ERROR, "Unable to get volume for cell %i since there are only "
               "%i cells", cell, _num_cells);
  
  return _volume[i];
}
