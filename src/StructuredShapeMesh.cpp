#include "StructuredShapeMesh.h"


/**
 * @brief Constructor sets the ID and unique ID for the Material.
 * @param id the user-defined ID for the material
 */
StructuredShapeMesh::StructuredShapeMesh(double width, double height, double depth, int num_x, int num_y, int num_z) : ShapeMesh(width, height, depth){

  setNumX(num_x);
  setNumY(num_y);
  setNumZ(num_z);

  _mesh_type = STRUCTURED_SHAPE_MESH;
  
  for (int c=0; c < 8; c++)
    _dif_linear[c] = NULL;
  
  return;
}


/**
 * @brief Destructor deletes all cross-section data structures from memory.
 */
StructuredShapeMesh::~StructuredShapeMesh() {

  for (int c=0; c < 8; c++){
    if (_dif_linear[c] != NULL)
      delete [] _dif_linear[c];
  }

  _dif_linear.clear();
}


void StructuredShapeMesh::setDifLinearByValue(double dif_linear, int cell, int group, int side, int position){

  if (cell >= getNumCells() || cell < 0)
    log_printf(ERROR, "Unable to set dif linear by value for Structured Shape Mesh "
               "cell %i since there are only %i cells", cell, getNumCells());

  if (group >= _num_shape_energy_groups || group < 0)
    log_printf(ERROR, "Unable to set dif linear by value for Structured Shape Mesh group "
               "%i since there are only %i groups", group, _num_shape_energy_groups);

  if (position >= 8 || position < 0)
    log_printf(ERROR, "Unable to set dif linear by value for Structured Shape Mesh clock "
               "position %i since there are only 8 clock positions", position);

  if (side < 0 || side >= 6)
    log_printf(NORMAL, "Unable to set dif linear by value for Structured Shape Mesh "
               "for side %i since there are only 6 sides", side);

  _dif_linear[position][(cell*6 + side)*_num_shape_energy_groups + group] = dif_linear;
}


double StructuredShapeMesh::getDifLinearByValue(int cell, int group, int side, int position){

  if (cell >= getNumCells() || cell < 0)
    log_printf(ERROR, "Unable to get dif linear by value for Structured Shape Mesh "
               "cell %i since there are only %i cells", cell, getNumCells());

  if (group >= _num_shape_energy_groups || group < 0)
    log_printf(ERROR, "Unable to get dif linear by value for Structured Shape Mesh group "
               "%i since there are only %i groups", group, _num_shape_energy_groups);

  if (position >= 8 || position < 0)
    log_printf(ERROR, "Unable to get dif linear by value for Structured Shape Mesh clock "
               "position %i since there are only 8 clock positions", position);

  if (side < 0 || side >= 6)
    log_printf(NORMAL, "Unable to get dif linear by value for Structured Shape Mesh "
               "for side %i since there are only 6 sides", side);

  return _dif_linear[position][(cell*6 + side)*_num_shape_energy_groups + group];
}


void StructuredShapeMesh::copyDifLinear(int position_from, int position_to){

  if (position_from >= 8 || position_from < 0)
    log_printf(ERROR, "Unable to copy dif linear from position %i since "
               "there are only 8 clock positions", position_from);

  if (position_to >= 8 || position_to < 0)
    log_printf(ERROR, "Unable to copy dif linear to position %i since "
               "there are only 8 clock positions", position_to);
    
  if (_dif_linear[position_from] == NULL)
    log_printf(ERROR, "Unable to copy dif linear from position %i since "
               "the dif linear have not been initialized yet", position_from);

  if (_dif_linear[position_to] == NULL)
    log_printf(ERROR, "Unable to copy dif linear to position %i since "
               "the dif linear have not been initialized yet", position_to);

  std::copy(_dif_linear[position_from], 
            _dif_linear[position_from] + getNumCells()*_num_shape_energy_groups*6, 
            _dif_linear[position_to]);
}


void StructuredShapeMesh::initialize(){

  ShapeMesh::initialize();

  for(int c=0; c < 8; c++){
    _dif_linear[c] = new double[getNumCells() * _num_shape_energy_groups * 6];
    _flux[c] = new double[getNumCells() * _num_shape_energy_groups];

    memset(_dif_linear[c], 0.0, sizeof(double) * getNumCells() * 
           _num_shape_energy_groups * 6);
    memset(_flux[c], 1.0, sizeof(double) * getNumCells() * _num_shape_energy_groups);
  }  
}


void StructuredShapeMesh::computeDifCoefs(int position){

  if (_num_shape_energy_groups == 0)
    log_printf(ERROR, "Unable to compute the dif coefs since the number "
               "of shape energy groups has not been set");

  int nx = _num_x;
  int ny = _num_y;
  int nz = _num_z;
  int ng = _num_shape_energy_groups;
  double width = getCellWidth();
  double height = getCellHeight();
  double depth = getCellDepth();
  double* temps = _temperature[position];
  
  for (int z=0; z < nz; z++){
    #pragma omp parallel for
    for (int y=0; y < ny; y++){
      for (int x=0; x < nx; x++){
        
        
        int cell = z*nx*ny + y*nx + x;
        double temp = temps[cell];
        int sense;
        double length, length_perpen;
        double dif_coef, current, flux, dif_linear, flux_next, dif_coef_next;
        
        for (int s=0; s < 6; s++){
          
          int cell_next = getNeighborCell(x, y, z, s);
          
          if (s == 0 || s == 1 || sense == 2)
            sense = -1;
          else
            sense = 1;
          
          if (s == 0 || s == 3){
            length = height * depth;
            length_perpen = width;
          }
          else if (s == 1 || s == 4){
            length = width * depth;
            length_perpen = height;
          }
          else{
            length = width*height;
            length_perpen = depth;            
          }
          
          for (int g=0; g < ng; g++){
            
            dif_coef = _materials[cell]->getDifCoefByGroup(g, position, temp);
            current = getCurrentByValue(cell, g, s, position);
            flux = getFluxByValue(cell, g, position);
            
            if (cell_next == -1){
              dif_linear = 2 * dif_coef / length_perpen / (1 + 4 * dif_coef / length_perpen);
              dif_linear *= _boundaries[s];
            }
            else{
              flux_next = getFluxByValue(cell_next, g, position);
              dif_coef_next = _materials[cell_next]->getDifCoefByGroup(g, position, temps[cell_next]);
              dif_linear = 2 * dif_coef * dif_coef_next / (length_perpen * dif_coef +
                                                           length_perpen * dif_coef_next);
            }
            
            setDifLinearByValue(dif_linear, cell, g, s, position);
          }
        }
      }
    }
  }
}


StructuredShapeMesh* StructuredShapeMesh::uniformRefine(int refine_x, int refine_y, int refine_z){

  StructuredShapeMesh* mesh = clone();

  int nx = _num_x * refine_x;
  int ny = _num_y * refine_y;
  int nz = _num_z * refine_z;
  mesh->setNumX(nx);
  mesh->setNumY(ny);
  mesh->setNumZ(nz);
  mesh->initialize();

  for (int z=0; z < nz; z++){
    int zz = z / refine_z;
    for (int y=0; y < ny; y++){
      int yy = y / refine_y;
      for (int x=0; x < nx; x++){    
        int xx = x / refine_x;
        mesh->setMaterial(getMaterial(zz*_num_x*_num_y + yy*_num_x+xx), z*nx*ny + y*nx+x);
        mesh->getTemperature(CURRENT)[z*nx*ny+y*nx+x] = _temperature[CURRENT][zz*_num_x*_num_y + yy*_num_x+xx];
      }
    }
  }

  return mesh;
}


StructuredShapeMesh* StructuredShapeMesh::clone(){

  if (_num_shape_energy_groups == 0)
    log_printf(ERROR, "Unable to clone the StructuredShapeMesh since the number "
               "of shape energy groups has not been set");

  if (_num_amp_energy_groups == 0)
    log_printf(ERROR, "Unable to clone the StructuredShapeMesh since the number "
               "of amp energy groups has not been set");

  if (_num_delayed_groups == 0)
    log_printf(ERROR, "Unable to clone the StructuredShapeMesh since the number "
               "of delayed groups has not been set");
  
  StructuredShapeMesh* mesh = new StructuredShapeMesh(getWidth(), getHeight(), getDepth(), 
                                                      _num_x, _num_y, _num_z);

  mesh->setNumShapeEnergyGroups(_num_shape_energy_groups);
  mesh->setNumAmpEnergyGroups(_num_amp_energy_groups);
  mesh->setNumDelayedGroups(_num_delayed_groups);
  mesh->setBuckling(_buckling);

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


int StructuredShapeMesh::getNumCells(){
  return _num_x*_num_y*_num_z;
}


double StructuredShapeMesh::getCellVolume(int cell){

  if (cell >= _num_cells || cell < 0)
    log_printf(ERROR, "Unable to get volume for cell %i since there are only "
               "%i cells", cell, _num_cells);
  
  return _cell_width*_cell_height*_cell_depth;
}


int StructuredShapeMesh::getNumX(){
  return _num_x;
}


int StructuredShapeMesh::getNumY(){
  return _num_y;
}


int StructuredShapeMesh::getNumZ(){
  return _num_z;
}


double StructuredShapeMesh::getCellWidth(){
  return _cell_width;
}


double StructuredShapeMesh::getCellHeight(){
  return _cell_height;
}


double StructuredShapeMesh::getCellDepth(){
  return _cell_depth;
}


void StructuredShapeMesh::setNumX(int num_x){

  if (num_x < 1)
    log_printf(ERROR, "Unable to set num x for Structured Shape Mesh to non"
               " positive number: %i", num_x);

  _num_x = num_x;
  _cell_width = (_x_max - _x_min) / num_x;
}


void StructuredShapeMesh::setNumY(int num_y){

  if (num_y < 1)
    log_printf(ERROR, "Unable to set num y for Structured Shape Mesh to non"
               " positive number: %i", num_y);

  _num_y = num_y;
  _cell_height = (_y_max - _y_min) / num_y;
}


void StructuredShapeMesh::setNumZ(int num_z){

  if (num_z < 1)
    log_printf(ERROR, "Unable to set num z for Structured Shape Mesh to non"
               " positive number: %i", num_z);

  _num_z = num_z;
  _cell_depth = (_z_max - _z_min) / num_z;
}


int StructuredShapeMesh::findCell(double x, double y, double z){

  if (x > _x_max || x < _x_min)
    log_printf(ERROR, "Unable to find cell with x value %f outside "
               "the mesh x bounds (%f, %f)", x, _x_min, _x_max);

  if (y > _y_max || y < _y_min)
    log_printf(ERROR, "Unable to find cell with y value %f outside "
               "the mesh y bounds (%f, %f)", y, _y_min, _y_max);

  if (z > _z_max || z < _z_min)
    log_printf(ERROR, "Unable to find cell with z value %f outside "
               "the mesh z bounds (%f, %f)", z, _z_min, _z_max);

  int i = floor((x - _x_min) / _cell_width);
  int j = floor((y - _y_min) / _cell_height);
  int k = floor((z - _z_min) / _cell_depth);

  return k * _num_x * _num_y + j * _num_x + i;  
}
