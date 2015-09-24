#include "Geometry.h"


/**
 * @brief Constructor sets the ID and unique ID for the Material.
 * @param id the user-defined ID for the material
 */
Geometry::Geometry(double width, double height, double depth) {

  _materials = NULL;
  _shape_to_amp = NULL;
  _volumes = NULL;
  _num_energy_groups = 0;
  _num_delayed_groups = 0;
  
  /* Set amp mesh properties */
  setAmpMeshDimensions();
  _num_amp_cells = 1;
  _num_shape_cells = 0;
  
  /* Set mesh properties */
  setXMin(-width/2.0);
  setXMax(width/2.0);
  setYMin(-height/2.0);
  setYMax(height/2.0);
  setZMin(-depth/2.0);
  setZMax(depth/2.0);

  _boundaries = new int[NUM_SURFACES];
  for (int i=0; i < NUM_SURFACES; i++)
    _boundaries[i] = REFLECTIVE;
  
  return;
}


/**
 * @brief Destructor deletes all cross-section data structures from memory.
 */
Geometry::~Geometry() {
}


void Geometry::setAmpMeshDimensions(int num_x, int num_y, int num_z){

  if (_amp_to_shape.empty() == false)
    _amp_to_shape.clear();

  _num_x_amp = num_x;
  _num_y_amp = num_y;
  _num_z_amp = num_z;
  _num_amp_cells = num_x * num_y * num_z;

  /* Allocate memory for mesh cell FSR vectors */
  for (int z = 0; z < num_z; z++){
    for (int y = 0; y < num_y; y++){
      for (int x = 0; x < num_x; x++){
        std::vector<int> *cells = new std::vector<int>;
        _amp_to_shape.push_back(*cells);
      }
    }
  }  
}


int Geometry::getNumXAmp(){
  return _num_x_amp;
}


int Geometry::getNumYAmp(){
  return _num_y_amp;
}


int Geometry::getNumZAmp(){
  return _num_z_amp;
}


void Geometry::setXMin(double x_min){
  _x_min = x_min;
}


void Geometry::setXMax(double x_max){
  _x_max = x_max;
}


void Geometry::setYMin(double y_min){
  _y_min = y_min;
}


void Geometry::setYMax(double y_max){
  _y_max = y_max;
}

void Geometry::setZMin(double z_min){
  _z_min = z_min;
}


void Geometry::setZMax(double z_max){
  _z_max = z_max;
}


double Geometry::getXMin(){
  return _x_min;
}


double Geometry::getXMax(){
  return _x_max;
}


double Geometry::getYMin(){
  return _y_min;
}


double Geometry::getYMax(){
  return _y_max;
}


double Geometry::getZMin(){
  return _z_min;
}


double Geometry::getZMax(){
  return _z_max;
}


double Geometry::getWidth(){
  return _x_max - _x_min;
}


double Geometry::getHeight(){
  return _y_max - _y_min;
}


double Geometry::getDepth(){
  return _z_max - _z_min;
}


void Geometry::setBoundary(int side, int boundary){

  if (side < 0 || side >= NUM_SURFACES)
    log_printf(ERROR, "Unable to set boundary for side %d as there are only"
               " %d geometry sides", side, NUM_SURFACES);

  _boundaries[side] = boundary;
}


int Geometry::getBoundary(int side){

  if (side < 0 || side >= NUM_SURFACES)
    log_printf(ERROR, "Unable to get boundary for side %d as there are only"
               " %d geometry sides", side, NUM_SURFACES);
  
  return _boundaries[side];
}


Material* Geometry::getMaterial(int cell){

  if (_materials == NULL)
    log_printf(ERROR, "Cannot get material since the materials "
               "array has not been initialized");

  return _materials[cell];
}


void Geometry::setMaterial(Material* material, int cell){

  if (_materials == NULL)
    log_printf(ERROR, "Cannot set material since the materials "
               "array has not been initialized");

  _materials[cell] = material;
  _num_energy_groups = material->getNumEnergyGroups();

  if (_num_delayed_groups == 0)
    _num_delayed_groups = material->getNumDelayedGroups();
}


void Geometry::uniquifyMaterials(){

  if (_materials == NULL)
    log_printf(ERROR, "Cannot uniquify materials since the materials "
               "array has not been initialized");

  Material* material;
  
  for (int i=0; i < _num_shape_cells; i++){
    material = _materials[i]->clone();
    setMaterial(material, i);
  }
}


void Geometry::setNumShapeCells(int num_shape_cells){

  if (num_shape_cells <= 0)
    log_printf(ERROR, "Unable to set num shape cells to non-positive "
               "number: %i", num_shape_cells);

  _num_shape_cells = num_shape_cells;

  if (_shape_to_amp != NULL)
    delete [] _shape_to_amp;

  if (_materials != NULL)
    delete [] _materials;

  if (_volumes != NULL)
    delete [] _volumes;

  _shape_to_amp = new int[num_shape_cells];
  _materials = new Material*[num_shape_cells];
  _volumes = new double[num_shape_cells];
}


int Geometry::getNumShapeCells(){
  return _num_shape_cells;
}


int Geometry::getNumAmpCells(){
  return _num_amp_cells;
}


void Geometry::addShapeCellToAmpCell(int shape_cell, int amp_cell){
  _amp_to_shape[amp_cell].push_back(shape_cell);
  _shape_to_amp[shape_cell] = amp_cell;
}


int Geometry::findAmpCellContainingShapeCell(int shape_cell){
  return _shape_to_amp[shape_cell];
}


std::vector< std::vector<int> >* Geometry::getAmpToShapeMap(){
  return &_amp_to_shape;
}


int* Geometry::getShapeToAmpMap(){
  return _shape_to_amp;
}


int Geometry::getNeighborAmpCell(int x, int y, int z, int side){

  if (side < 0 || side >= NUM_SURFACES)
    log_printf(ERROR, "Unable to get neighbor cell for side %d as there are only"
               " %d geometry sides", side, NUM_SURFACES);

  int neighbor_cell = -1;

  if (side == SURFACE_X_MIN){
    if (x != 0)
      neighbor_cell = z*_num_x_amp*_num_y_amp + y*_num_x_amp + x - 1;
  }
  else if (side == SURFACE_Y_MIN){
    if (y != 0)
      neighbor_cell = z*_num_x_amp*_num_y_amp + (y-1)*_num_x_amp + x;
  }
  else if (side == SURFACE_Z_MIN){
    if (z != 0)
      neighbor_cell = (z-1)*_num_x_amp*_num_y_amp + y*_num_x_amp + x;
  }
  else if (side == SURFACE_X_MAX){
    if (x != _num_x_amp-1)
      neighbor_cell = z*_num_x_amp*_num_y_amp + y*_num_x_amp + x + 1;
  }
  else if (side == SURFACE_X_MAX){
    if (y != _num_y_amp-1)
      neighbor_cell = z*_num_x_amp*_num_y_amp + (y+1)*_num_x_amp + x;
  }
  else if (side == SURFACE_Z_MAX){
    if (z != _num_z_amp-1)
      neighbor_cell = (z+1)*_num_x_amp*_num_y_amp + y*_num_x_amp + x;
  }

  return neighbor_cell;
}


Geometry* Geometry::clone(){

  Geometry* geometry = new Geometry(getWidth(), getHeight(), getDepth());
  geometry->setAmpMeshDimensions(_num_x_amp, _num_y_amp, _num_z_amp);
  geometry->setNumShapeCells(_num_shape_cells);

  for (int i=0; i < NUM_SURFACES; i++)
    geometry->setBoundary(i, _boundaries[i]);

  for (int i=0; i < _num_shape_cells; i++)
    geometry->setMaterial(_materials[i]->clone(), i);

  return geometry;
}


int Geometry::getNumEnergyGroups(){
  return _num_energy_groups;
}


int Geometry::getNumDelayedGroups(){
  return _num_delayed_groups;
}


double Geometry::getVolume(int cell){
  return _volumes[cell];
}