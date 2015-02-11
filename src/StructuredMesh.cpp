#include "StructuredMesh.h"


/**
 * @brief Constructor sets the ID and unique ID for the Material.
 * @param id the user-defined ID for the material
 */
StructuredMesh::StructuredMesh(double width, double height, int num_x, int num_y) :
  Mesh(width, height){

  /* Set mesh properties */
  setNumX(num_x);
  setNumY(num_y);
  _fuel_volume = 0.0;
  
  for (int c=0; c < 8; c++){
    _current[c] = NULL;
    _dif_linear[c] = NULL;
    _dif_nonlinear[c] = NULL;
  }

  return;
}


/**
 * @brief Destructor deletes all cross-section data structures from memory.
 */
StructuredMesh::~StructuredMesh() {

  for (int c=0; c < 8; c++){
    if (_current[c] != NULL)
      delete [] _current[c];

    if (_dif_linear[c] != NULL)
      delete [] _dif_linear[c];

    if (_dif_nonlinear[c] != NULL)
      delete [] _dif_nonlinear[c];

  }

  _current.clear();
  _dif_linear.clear();
  _dif_nonlinear.clear();
}


void StructuredMesh::computeFuelVolume(){

  _fuel_volume = 0.0;
  for (int i=0; i < _num_x*_num_y; i++){
    if (_materials[i]->isFissionable())
      _fuel_volume += getCellVolume();
  }
}


int StructuredMesh::getNeighborCell(int x, int y, int side){

  if (side < 0 || side > 3)
    log_printf(ERROR, "Unable to get neighbor cell for side %d as there are only"
               " 4 geometry sides", side);

  int neighbor_cell = -1;

  if (side == 0){
    if (x != 0)
      neighbor_cell = y*_num_x + x - 1;
  }
  else if (side == 1){
    if (y != 0)
      neighbor_cell = (y-1)*_num_x + x;
  }
  else if (side == 2){
    if (x != _num_x-1)
      neighbor_cell = y*_num_x + x + 1;
  }
  else if (side == 3){
    if (y != _num_y-1)
      neighbor_cell = (y+1)*_num_x + x;
  }

  return neighbor_cell;
}


Material* StructuredMesh::getNeighborMaterial(int x, int y, int side){

  int neighbor_cell = getNeighborCell(x, y, side);

  if (neighbor_cell == -1)
    return NULL;
  else
    return _materials[neighbor_cell];
  
}


int StructuredMesh::getNumX(){
  return _num_x;
}


int StructuredMesh::getNumY(){
  return _num_y;
}


double StructuredMesh::getCellWidth(){
  return _cell_width;
}


double StructuredMesh::getCellHeight(){
  return _cell_height;
}


double StructuredMesh::getCellVolume(){
  return _cell_width * _cell_height;
}


void StructuredMesh::setNumX(int num_x){
  _num_x = num_x;
  _cell_width = (_x_max - _x_min) / num_x;
}


void StructuredMesh::setNumY(int num_y){
  _num_y = num_y;
  _cell_height = (_y_max - _y_min) / num_y;
}


int StructuredMesh::findCell(double x, double y){

  int i = floor((x - _x_min) / _cell_width);
  int j = floor((y - _y_min) / _cell_height);

  return j * _num_x + i;  
}


void StructuredMesh::uniquifyMaterials(){

  Material* material;
  
  for (int i=0; i < _num_x*_num_y; i++){
    material = _materials[i]->clone();
    _materials[i] = material;
  }
}


double StructuredMesh::getMaxTemperature(int position){

  double temp_max = _temperature[position][0];
  for (int i=0; i < _num_x*_num_y; i++){
    if (_temperature[position][i] > temp_max)
      temp_max = _temperature[position][i];
  }

  return temp_max;
}


void StructuredMesh::copyPower(int position_from, int position_to){
  std::copy(_power[position_from], _power[position_from] + _num_x*_num_y, _power[position_to]);
}


void StructuredMesh::copyTemperature(int position_from, int position_to){
  std::copy(_temperature[position_from], _temperature[position_from] + _num_x*_num_y, _temperature[position_to]);
}


double* StructuredMesh::getCurrent(int position){
  return _current[position];
}


double* StructuredMesh::getDifLinear(int position){
  return _dif_linear[position];
}


double* StructuredMesh::getDifNonlinear(int position){
  return _dif_nonlinear[position];
}


void StructuredMesh::setTemperature(double temperature){

  for (int c=0; c < 8; c++){
    for (int i=0; i < _num_x*_num_y; i++)
      _temperature[c][i] = temperature;
  }
}


void StructuredMesh::setClock(Clock* clock){

  _clock = clock;

  for (int i=0; i < _num_x*_num_y; i++)
    _materials[i]->setClock(clock);

}
