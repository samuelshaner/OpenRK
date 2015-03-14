#include "StructuredMesh.h"


/**
 * @brief Constructor sets the ID and unique ID for the Material.
 * @param id the user-defined ID for the material
 */
StructuredMesh::StructuredMesh(double width_x, double width_y, double width_z, \
                               int num_x, int num_y, int num_z) :
  Mesh(width_x, width_y, width_z){

  /* Set mesh properties */
  setNumX(num_x);
  setNumY(num_y);
  setNumZ(num_z);
  
  return;
}


/**
 * @brief Destructor deletes all cross-section data structures from memory.
 */
StructuredMesh::~StructuredMesh() {
}


int StructuredMesh::getNeighborCell(int x, int y, int z, int side){

  if (side < 0 || side > 5)
    log_printf(ERROR, "Unable to get neighbor cell for side %d as there are only"
               " 6 geometry sides", side);

  int neighbor_cell = -1;

  if (side == X_MIN_SURFACE){
    if (x != 0)
      neighbor_cell = z*_num_x*_num_y + y*_num_x + x - 1;
  }
  else if (side == Y_MIN_SURFACE){
    if (y != 0)
      neighbor_cell = z*_num_x*_num_y + (y-1)*_num_x + x;
  }
  else if (side == Z_MIN_SURFACE){
    if (z != 0)
      neighbor_cell = (z-1)*_num_x*_num_y + y*_num_x + x;
  }
  else if (side == X_MAX_SURFACE){
    if (x != _num_x-1)
      neighbor_cell = z*_num_x*_num_y + y*_num_x + x + 1;
  }
  else if (side == Y_MAX_SURFACE){
    if (y != _num_y-1)
      neighbor_cell = z*_num_x*_num_y + (y+1)*_num_x + x;
  }
  else if (side == Z_MAX_SURFACE){
    if (z != _num_z-1)
      neighbor_cell = (z+1)*_num_x*_num_y + y*_num_x + x;
  }

  return neighbor_cell;
}


int StructuredMesh::getNumX(){
  return _num_x;
}


int StructuredMesh::getNumY(){
  return _num_y;
}


int StructuredMesh::getNumZ(){
  return _num_z;
}


double StructuredMesh::getCellWidthX(){
  return _cell_width_x;
}


double StructuredMesh::getCellWidthY(){
  return _cell_width_y;
}

double StructuredMesh::getCellWidthZ(){
  return _cell_width_z;
}


double StructuredMesh::getCellVolume(){
  return _cell_width_x * _cell_width_y * _cell_width_z;
}


void StructuredMesh::setNumX(int num_x){
  _num_x = num_x;
  _cell_width_x = (_x_max - _x_min) / num_x;
}


void StructuredMesh::setNumY(int num_y){
  _num_y = num_y;
  _cell_width_y = (_y_max - _y_min) / num_y;
}


void StructuredMesh::setNumZ(int num_z){
  _num_z = num_z;
  _cell_width_z = (_z_max - _z_min) / num_z;
}


int StructuredMesh::findCell(double x, double y, double z){

  int i = floor((x - _x_min) / _cell_width_x);
  int j = floor((y - _y_min) / _cell_width_y);
  int k = floor((z - _z_min) / _cell_width_z);

  return k * _num_x * _num_y + j * _num_x + i;  
}


int StructuredMesh::getNumCells(){
  return _num_x * _num_y * _num_z;
}
