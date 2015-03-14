#include "StructuredShapeMesh.h"


/**
 * @brief Constructor sets the ID and unique ID for the Material.
 * @param id the user-defined ID for the material
 */
StructuredShapeMesh::StructuredShapeMesh(double width, double height, double depth, int num_x, int num_y, int num_z) : StructuredMesh(width, height, depth, num_x, num_y, num_z){

  _amp_mesh = NULL;
  _group_indices = NULL;
  _amp_map = NULL;
  _buckling = 0.0;
  
  return;
}


/**
 * @brief Destructor deletes all cross-section data structures from memory.
 */
StructuredShapeMesh::~StructuredShapeMesh() {

  if (_group_indices != NULL)
    delete [] _group_indices;

  if (_amp_map != NULL)
    delete [] _amp_map;
}


void StructuredShapeMesh::setAmpMesh(AmpMesh* mesh){
  _amp_mesh = mesh;

  int nx = mesh->getNumX();
  int ny = mesh->getNumY();
  int nz = mesh->getNumZ();
  int num_refines_x = _num_x / nx;
  int num_refines_y = _num_y / ny;
  int num_refines_z = _num_z / nz;
  _amp_map = new int[_num_x * _num_y * _num_z];
  
  for (int z=0; z < _num_z; z++){
    int zz = z / num_refines_z;
    for (int y=0; y < _num_y; y++){
      int yy = y / num_refines_y;
      for (int x=0; x < _num_x; x++){
        int xx = x / num_refines_x;
        int amp_cell = zz*nx*ny + yy*nx + xx;
        int shape_cell = z*_num_x*_num_y + y*_num_x + x;
        _amp_map[shape_cell] = amp_cell;
      }
    }
  }
}


void StructuredShapeMesh::setFluxByValue(double flux, int cell, int group, int position){
  _flux[position][cell*_num_shape_energy_groups + group] = flux;
}


void StructuredShapeMesh::setCurrentByValue(double current, int cell, int group, int side, int position){
  _current[position][(cell*6 + side)*_num_shape_energy_groups + group] = current;
}


void StructuredShapeMesh::setDifLinearByValue(double dif_linear, int cell, int group, int side, int position){
  _dif_linear[position][(cell*6 + side)*_num_shape_energy_groups + group] = dif_linear;
}


void StructuredShapeMesh::setDifNonlinearByValue(double dif_nonlinear, int cell, int group, int side, int position){
  _dif_nonlinear[position][(cell*6 + side)*_num_shape_energy_groups + group] = dif_nonlinear;
}


double StructuredShapeMesh::getFluxByValue(int cell, int group, int position){
  return _flux[position][cell*_num_shape_energy_groups + group];
}


double StructuredShapeMesh::getCurrentByValue(int cell, int group, int side, int position){
  return _current[position][(cell*6 + side)*_num_shape_energy_groups + group];
}


double StructuredShapeMesh::getDifLinearByValue(int cell, int group, int side, int position){
  return _dif_linear[position][(cell*6 + side)*_num_shape_energy_groups + group];
}


double StructuredShapeMesh::getDifNonlinearByValue(int cell, int group, int side, int position){
  return _dif_nonlinear[position][(cell*6 + side)*_num_shape_energy_groups + group];
}


void StructuredShapeMesh::copyFlux(int position_from, int position_to){
  std::copy(_flux[position_from], _flux[position_from] + _num_x*_num_y*_num_z*_num_shape_energy_groups, _flux[position_to]);
}


void StructuredShapeMesh::copyCurrent(int position_from, int position_to){
  std::copy(_current[position_from], _current[position_from] + _num_x*_num_y*_num_z*_num_shape_energy_groups*6, _current[position_to]);
}


void StructuredShapeMesh::copyDifLinear(int position_from, int position_to){
  std::copy(_dif_linear[position_from], _dif_linear[position_from] + _num_x*_num_y*_num_z*_num_shape_energy_groups*6, _dif_linear[position_to]);
}


void StructuredShapeMesh::copyDifNonlinear(int position_from, int position_to){
  std::copy(_dif_nonlinear[position_from], _dif_nonlinear[position_from] + _num_x*_num_y*_num_z*_num_shape_energy_groups*6, _dif_nonlinear[position_to]);
}


void StructuredShapeMesh::initialize(){

  _materials = new Material*[_num_x*_num_y*_num_z];
    
  for(int c=0; c < 8; c++){
    _dif_linear[c] = new double[_num_x * _num_y * _num_z * _num_shape_energy_groups * 6];
    _dif_nonlinear[c] = new double[_num_x * _num_y * _num_z * _num_shape_energy_groups * 6];
    _current[c] = new double[_num_x * _num_y * _num_z * _num_shape_energy_groups * 6];
    _flux[c] = new double[_num_x * _num_y * _num_z * _num_shape_energy_groups];
    _temperature[c] = new double[_num_x * _num_y * _num_z];
    _power[c] = new double[_num_x * _num_y * _num_z];

    memset(_dif_linear[c], 0.0, sizeof(double) * _num_x * _num_y * _num_z * _num_shape_energy_groups * 6);
    memset(_dif_nonlinear[c], 0.0, sizeof(double) * _num_x * _num_y * _num_z * _num_shape_energy_groups * 6);
    memset(_current[c], 0.0, sizeof(double) * _num_x * _num_y * _num_z * _num_shape_energy_groups * 6);
    memset(_flux[c], 1.0, sizeof(double) * _num_x * _num_y * _num_z * _num_shape_energy_groups);
    memset(_temperature[c], 300.0, sizeof(double) * _num_x * _num_y * _num_z);
    memset(_power[c], 0.0, sizeof(double) * _num_x * _num_y * _num_z);
  }  
}


void StructuredShapeMesh::synthesizeFlux(int position){

  int ngs = _num_shape_energy_groups;
  int nga = _num_amp_energy_groups;
  double dt = _clock->getTime(FORWARD_OUT) - _clock->getTime(PREVIOUS_OUT);
  double wt_begin = (_clock->getTime(FORWARD_OUT) - _clock->getTime(position)) / dt;
  double wt_end = (_clock->getTime(position) - _clock->getTime(PREVIOUS_OUT)) / dt;

  #pragma omp parallel for
  for (int i=0; i < _num_x * _num_y * _num_z; i++){
    double shape_previous, shape_forward, shape_current;
    for (int g=0; g < ngs; g++){
      shape_previous = _flux[PREVIOUS_OUT][i*ngs+g]
        / _amp_mesh->getFlux(PREVIOUS_OUT)[_amp_map[i]*nga+_group_indices[g]];
      shape_forward = _flux[FORWARD_OUT][i*ngs+g]
        / _amp_mesh->getFlux(FORWARD_OUT)[_amp_map[i]*nga+_group_indices[g]];

      shape_current = wt_begin * shape_previous + wt_end * shape_forward;
      _flux[position][i*ngs+g] = shape_current *
        _amp_mesh->getFlux(position)[_amp_map[i]*nga+_group_indices[g]];

    }
  }
}


void StructuredShapeMesh::reconstructFlux(int position, int position_shape, int position_amp){

  int ngs = _num_shape_energy_groups;
  int nga = _num_amp_energy_groups;

  #pragma omp parallel for
  for (int i=0; i < _num_x * _num_y * _num_z; i++){
    double shape;
    for (int g=0; g < ngs; g++){
      shape = _flux[position_shape][i*ngs+g]
        / _amp_mesh->getFlux(position_shape)[_amp_map[i]*nga+_group_indices[g]];
      _flux[position][i*ngs+g] = shape *
        _amp_mesh->getFlux(position_amp)[_amp_map[i]*nga+_group_indices[g]];
    }
  }
}


void StructuredShapeMesh::computePower(int position){

  #pragma omp parallel for
  for (int i=0; i < _num_x * _num_y * _num_z; i++){
    double fission_rate = 0.0;
    double temp = _temperature[position][i];
    
    for (int g=0; g < _num_shape_energy_groups; g++){
      fission_rate += _materials[i]->getSigmaFByGroup(g, position, temp) * getFluxByValue(i, g, position);
    }
    
    _power[position][i] = fission_rate * _materials[i]->getEnergyPerFission();
  }
}


void StructuredShapeMesh::computeDifCoefs(int position){

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
  
  StructuredShapeMesh* mesh = new StructuredShapeMesh(getWidth(), getHeight(), getDepth(), _num_x, _num_y, _num_z);

  mesh->setNumShapeEnergyGroups(_num_shape_energy_groups);
  mesh->setNumAmpEnergyGroups(_num_amp_energy_groups);
  mesh->setNumDelayedGroups(_num_delayed_groups);
  mesh->setBuckling(_buckling);
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


void StructuredShapeMesh::computeInitialPrecursorConc(int position){

  double* temps = _temperature[position];

  #pragma omp parallel for
  for (int i=0; i < _num_x * _num_y * _num_z; i++){
    
    double fission_rate = 0.0;
    
    if (_materials[i]->isFissionable()){
      for (int g=0; g < _num_shape_energy_groups; g++){
        fission_rate += _materials[i]->getNuSigmaFByGroup(g, position, temps[i]) *
          getFluxByValue(i, g, position);
      }
      
      for (int d=0; d < _num_delayed_groups; d++){
        double conc = fission_rate * _delayed_fractions[d] / _k_eff_0 / _decay_constants[d];
        _materials[i]->setPrecursorConcByGroup(conc, d, position);
        
      }
    }
  }
}


void StructuredShapeMesh::integratePrecursorConc(int position_from, int position_to){

  double* temps_from = _temperature[position_from];
  double* temps_to = _temperature[position_to];
  double dt = _clock->getTime(position_to) - _clock->getTime(position_from);

  #pragma omp parallel for
  for (int i=0; i < _num_x * _num_y * _num_z; i++){

    double fission_rate_from = 0.0;
    double fission_rate_to = 0.0;
    
    if (_materials[i]->isFissionable()){
      for (int g=0; g < _num_shape_energy_groups; g++){
        fission_rate_from += _materials[i]->getNuSigmaFByGroup(g, position_from, temps_from[i]) *
          getFluxByValue(i, g, position_from);
        fission_rate_to += _materials[i]->getNuSigmaFByGroup(g, position_to, temps_to[i]) *
          getFluxByValue(i, g, position_to);
      }
      
      for (int d=0; d < _num_delayed_groups; d++){
        double k1 = exp(-_decay_constants[d] * dt);
        double k2 = 1.0 - (1.0 - k1) / (_decay_constants[d] * dt);
        double k3 = k1 + (k2 - 1.0);
        double conc = _materials[i]->getPrecursorConcByGroup(d, position_from);
        double new_conc = k1 * conc + k2 * _delayed_fractions[d] / _decay_constants[d] / _k_eff_0
          * fission_rate_to - k3 * _delayed_fractions[d] / _decay_constants[d] / _k_eff_0 * fission_rate_from;

        
        _materials[i]->setPrecursorConcByGroup(new_conc, d, position_to);
      }
    }
  }
}


void StructuredShapeMesh::integrateTemperature(int position_from, int position_to){

  double* temps_from = _temperature[position_from];
  double* temps_to = _temperature[position_to];
  double dt = _clock->getTime(position_to) - _clock->getTime(position_from);

  #pragma omp parallel for
  for (int i=0; i < _num_x * _num_y * _num_z; i++){

    double fission_rate_from = 0.0;
    double fission_rate_to = 0.0;
    
    if (_materials[i]->isFissionable()){
      for (int g=0; g < _num_shape_energy_groups; g++){
        fission_rate_from += _materials[i]->getNuSigmaFByGroup(g, position_from, temps_from[i]) *
          getFluxByValue(i, g, position_from);
        fission_rate_to += _materials[i]->getNuSigmaFByGroup(g, position_to, temps_to[i]) *
          getFluxByValue(i, g, position_to);
      }
      
      _temperature[position_to][i] = _temperature[position_from][i] + dt * 0.5 *
        (fission_rate_from + fission_rate_to) * _materials[i]->getTemperatureConversionFactor();
    }
  }
}


void StructuredShapeMesh::setGroupStructure(int* group_indices, int length_group_indices){

  /* Allocate memory */
  if (_group_indices != NULL)
    delete [] _group_indices;

  _group_indices = new int[_num_shape_energy_groups];

  if (group_indices[0] != 0)
    log_printf(ERROR, "The first value in group indices must be 0!");

  /* Set first group indice to 0 */
  int index = 0;
  
  /* Set MOC group bounds for rest of CMFD energy groups */
  for (int g=0; g < _num_shape_energy_groups; g++){
    if (g == group_indices[index])
      index++;
    
    _group_indices[g] = index-1;
  }      
}


void StructuredShapeMesh::scaleFlux(int position, double scale_val){

  #pragma omp parallel for
  for (int i=0; i < _num_x*_num_y*_num_z*_num_shape_energy_groups; i++)
    _flux[position][i] *= scale_val;
}


double StructuredShapeMesh::computeAveragePower(int position){

  if (_fuel_volume == 0.0)
    computeFuelVolume();   

  computePower(position);

  return pairwise_sum(_power[position], _num_x*_num_y*_num_z) * getCellVolume() / _fuel_volume;
}  


double StructuredShapeMesh::computePowerL2Norm(int position_1, int position_2){

  computePower(position_1);
  computePower(position_2);

  double* power_residual = new double[_num_x * _num_y * _num_z];
  memset(power_residual, 0.0, sizeof(double) * _num_x * _num_y * _num_z);

  #pragma omp parallel for
  for (int i=0; i < _num_x * _num_y * _num_z; i++){
    if (_power[position_1][i] > 0.0){
      power_residual[i] = pow((_power[position_1][i] - _power[position_2][i]) / _power[position_1][i], 2);
    }
  }
  
  double residual = sqrt(pairwise_sum(power_residual, _num_x*_num_y*_num_z));

  delete [] power_residual;
  return residual;
}


int StructuredShapeMesh::getAmpGroup(int shape_group){
  return _group_indices[shape_group];
}


void StructuredShapeMesh::uniquifyMaterials(){

  Material* material;
  
  for (int i=0; i < _num_x*_num_y*_num_z; i++){
    material = _materials[i]->clone();
    _materials[i] = material;
  }
}


double StructuredShapeMesh::getMaxTemperature(int position){

  double temp_max = _temperature[position][0];
  for (int i=0; i < _num_x*_num_y*_num_z; i++){
    if (_temperature[position][i] > temp_max)
      temp_max = _temperature[position][i];
  }

  return temp_max;
}
