#include "GeometryDiffusion.h"


/**
 * @brief Constructor sets the ID and unique ID for the Material.
 * @param id the user-defined ID for the material
 */
GeometryDiffusion::GeometryDiffusion(double width, double height, double depth) : Geometry(width, height, depth){

  /* Set fine mesh properties */
  setShapeMeshDimensions();
  
  return;
}


/**
 * @brief Destructor deletes all cross-section data structures from memory.
 */
GeometryDiffusion::~GeometryDiffusion() {
}


void GeometryDiffusion::setShapeMeshDimensions(int num_x, int num_y, int num_z){

  _num_x_shape = num_x;
  _num_y_shape = num_y;
  _num_z_shape = num_z;

  setNumShapeCells(num_x * num_y * num_z);
  
  for (int i=0; i < _num_shape_cells; i++){
    _volumes[i] = (getWidth()/num_x) * (getHeight()/num_y) * (getDepth()/num_z);
  }
}


int GeometryDiffusion::getNumXShape(){
  return _num_x_shape;
}


int GeometryDiffusion::getNumYShape(){
  return _num_y_shape;
}


int GeometryDiffusion::getNumZShape(){
  return _num_z_shape;
}


int GeometryDiffusion::getNeighborShapeCell(int x, int y, int z, int side){

  if (side < 0 || side >= NUM_SURFACES)
    log_printf(ERROR, "Unable to get neighbor cell for side %d as there are only"
               " %d geometry sides", side, NUM_SURFACES);

  int neighbor_cell = -1;

  if (side == SURFACE_X_MIN){
    if (x != 0)
      neighbor_cell = z*_num_x_shape*_num_y_shape + y*_num_x_shape + x - 1;
  }
  else if (side == SURFACE_Y_MIN){
    if (y != 0)
      neighbor_cell = z*_num_x_shape*_num_y_shape + (y-1)*_num_x_shape + x;
  }
  else if (side == SURFACE_Z_MIN){
    if (z != 0)
      neighbor_cell = (z-1)*_num_x_shape*_num_y_shape + y*_num_x_shape + x;
  }
  else if (side == SURFACE_X_MAX){
    if (x != _num_x_shape-1)
      neighbor_cell = z*_num_x_shape*_num_y_shape + y*_num_x_shape + x + 1;
  }
  else if (side == SURFACE_Y_MAX){
    if (y != _num_y_shape-1)
      neighbor_cell = z*_num_x_shape*_num_y_shape + (y+1)*_num_x_shape + x;
  }
  else if (side == SURFACE_Z_MAX){
    if (z != _num_z_shape-1)
      neighbor_cell = (z+1)*_num_x_shape*_num_y_shape + y*_num_x_shape + x;
  }

  return neighbor_cell;
}


GeometryDiffusion* GeometryDiffusion::uniformRefine(int refine_x, int refine_y, int refine_z){

  GeometryDiffusion* geometry = clone();

  int nx = _num_x_shape * refine_x;
  int ny = _num_y_shape * refine_y;
  int nz = _num_z_shape * refine_z;
  geometry->setShapeMeshDimensions(nx, ny, nz);

  for (int z=0; z < nz; z++){
    int zz = z / refine_z;
    for (int y=0; y < ny; y++){
      int yy = y / refine_y;
      for (int x=0; x < nx; x++){    
        int xx = x / refine_x;
        geometry->setMaterial(getMaterial(zz*_num_x_shape*_num_y_shape + yy*_num_x_shape + xx), z*nx*ny + y*nx + x);
      }
    }
  }

  return geometry;
}


GeometryDiffusion* GeometryDiffusion::clone(){

  GeometryDiffusion* geometry = new GeometryDiffusion(getWidth(), getHeight(), getDepth());
  geometry->setAmpMeshDimensions(_num_x_amp, _num_y_amp, _num_z_amp);
  geometry->setShapeMeshDimensions(_num_x_shape, _num_y_shape, _num_z_shape);
  geometry->setNumShapeCells(_num_shape_cells);

  for (int i=0; i < NUM_SURFACES; i++)
    geometry->setBoundary(i, _boundaries[i]);

  for (int i=0; i < _num_shape_cells; i++)
    geometry->setMaterial(_materials[i]->clone(), i);

  return geometry;
}


void GeometryDiffusion::generateCellMap(){

  int num_refines_x = _num_x_shape / _num_x_amp;
  int num_refines_y = _num_y_shape / _num_y_amp;
  int num_refines_z = _num_z_shape / _num_z_amp;

  for (int z=0; z < _num_z_shape; z++){
    int zz = z / num_refines_z;
    for (int y=0; y < _num_y_shape; y++){
      int yy = y / num_refines_y;
      for (int x=0; x < _num_x_shape; x++){
        int xx = x / num_refines_x;
        int amp_cell = zz*_num_x_amp*_num_y_amp + yy*_num_x_amp + xx;
        int shape_cell = z*_num_x_shape*_num_y_shape + y*_num_x_shape + x;
        addShapeCellToAmpCell(shape_cell, amp_cell);
      }
    }
  }  
}


int GeometryDiffusion::findCell(double x, double y, double z){
  
  double width = getWidth() / _num_x_shape;
  double height = getHeight() / _num_y_shape;
  double depth = getDepth() / _num_z_shape;
  
  int xc = floor( (x - getXMin()) / width);
  int yc = floor( (y - getYMin()) / height);
  int zc = floor( (z - getZMin()) / depth);

  return zc*(_num_x_shape*_num_y_shape) + yc*_num_x_shape + xc;
}
