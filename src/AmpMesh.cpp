#include "AmpMesh.h"


/**
 * @brief Constructor sets the ID and unique ID for the Material.
 * @param id the user-defined ID for the material
 */
AmpMesh::AmpMesh(double width_x, double width_y, double width_z, int num_x, int num_y, int num_z) :
  StructuredMesh(width_x, width_y, width_z, num_x, num_y, num_z){

  _shape_mesh = NULL;
  
  return;
}


/**
 * @brief Destructor deletes all cross-section data structures from memory.
 */
AmpMesh::~AmpMesh() {
}


void AmpMesh::setShapeMesh(StructuredShapeMesh* mesh){
  _shape_mesh = mesh;

  for (int i=0; i < _num_x*_num_y*_num_z; i++){
    std::vector<int> *shape_cells = new std::vector<int>;
    _shape_map.push_back(*shape_cells);
  }
  
  int nx = mesh->getNumX();
  int ny = mesh->getNumY();
  int nz = mesh->getNumZ();
  int num_refines_x = nx / _num_x;
  int num_refines_y = ny / _num_y;
  int num_refines_z = nz / _num_z;

  for (int z=0; z < nz; z++){
    int zz = z / num_refines_z;
    for (int y=0; y < ny; y++){
      int yy = y / num_refines_y;
      for (int x=0; x < nx; x++){
        int xx = x / num_refines_x;
        int amp_cell = zz*_num_x*_num_y + yy*_num_x + xx;
        int shape_cell = z*nx*ny + y*nx + x;
        _shape_map[amp_cell].push_back(shape_cell);
      }
    }
  }
}


void AmpMesh::computePower(int position){

  double* temperature = _shape_mesh->getFieldVariable(TEMPERATURE, position);
  Material* matrial;
  
  for (int i=0; i < getNumCells(); i++){
    double fission_rate = 0.0;
    double temp = _field_variables[TEMPERATURE][position][i];
    std::vector<int>::iterator iter;

    for (iter=_shape_map[i].begin(); iter != _shape_map[i].end(); iter++){
      material = _shape_mesh->getMaterial(*iter);
      for (int g=0; g < _num_energy_groups; g++)
        fission_rate += material->getSigmaFByGroup(g, position, temperature[*iter])\
          * getFieldVariableByValue(AMPLITUDE, position, i, g)
          * _shape_mesh->getFieldVariableByValue(SHAPE, position, *iter, g);
    }
    
    setFieldVariableByValue(POWER, position, fission_rate, i);
  }
}


void AmpMesh::computeCurrent(int position){

  int sm_nx = _shape_mesh->getNumX();
  int sm_ny = _shape_mesh->getNumY();
  int sm_nz = _shape_mesh->getNumZ();
  int sm_cw = _shape_mesh->getCellWidth();
  int sm_ch = _shape_mesh->getCellHeight();
  int sm_cd = _shape_mesh->getCellDepth();
  int num_refines_x = sm_nx / _num_x;
  int num_refines_y = sm_ny / _num_y;
  int num_refines_z = sm_nz / _num_z;
  double* temps = _shape_mesh->getTemperature(position);

  #pragma omp parallel for 
  for (int i=0; i < _num_x * _num_y * _num_z; i++){

    std::vector<int>::iterator iter;
    double flux, dif_linear, d, d_next, temp, temp_next, flux_next;
    int cell_next;
    double length_perpen, length;
    Material *mat, *mat_next;

    for (int g=0; g < _num_energy_groups; g++){
      for (int s=0; s < 6; s++){
        double current = 0.0;
        for (iter = _shape_map[i].begin(); iter != _shape_map[i].end(); ++iter){

          int x = (*iter % (sm_nx*sm_ny)) % sm_nx;
          int y = (*iter % (sm_nx*sm_ny)) / sm_nx;
          int z = (*iter / (sm_nx*sm_ny));
          mat = _shape_mesh->getMaterial(*iter);
          temp = temps[*iter];
          cell_next = _shape_mesh->getNeighborCell(x, y, z, s);
          temp_next = temps[cell_next];
          mat_next = _shape_mesh->getNeighborMaterial(x, y, z, s);

          if (s == 0 && x % num_refines_x == 0){
            flux = _shape_mesh->getFluxByValue(*iter, g, position);
            d = mat->getDifCoefByGroup(g, position, temp);
            length_perpen = sm_cw;
            length = sm_ch * sm_cd;
            
            if (mat_next == NULL){
              dif_linear = 2 * d / length_perpen / (1 + 4 * d / length_perpen);
              dif_linear *= _shape_mesh->getBoundary(s);
              current += - dif_linear * flux * length;
            }
            else{
              d_next = mat_next->getDifCoefByGroup(g, position, temp_next);
              dif_linear = 2 * d * d_next / (length_perpen * d + length_perpen * d_next);
              flux_next = _shape_mesh->getFluxByValue(cell_next, gg, position);
              current += - dif_linear * (flux - flux_next) * length;
            }
          }

          else if (s == 1 && y % num_refines_y == 0){
            flux = _shape_mesh->getFluxByValue(*iter, g, position);
            d = mat->getDifCoefByGroup(g, position, temp);
            length_perpen = sm_ch;
            length = sm_cw * sm_cd;
            
            if (mat_next == NULL){
              dif_linear = 2 * d / length_perpen / (1 + 4 * d / length_perpen);
              dif_linear *= _shape_mesh->getBoundary(s);
              current += - dif_linear * flux * length;
            }
            else{
              d_next = mat_next->getDifCoefByGroup(g, position, temp_next);
              dif_linear = 2 * d * d_next / (length_perpen * d + length_perpen * d_next);
              flux_next = _shape_mesh->getFluxByValue(cell_next, g, position);
              current += - dif_linear * (flux - flux_next) * length;
            }
          }

          else if (s == 2 && z % num_refines_z == 0){
            flux = _shape_mesh->getFluxByValue(*iter, g, position);
            d = mat->getDifCoefByGroup(g, position, temp);
            length_perpen = sm_cd;
            length = sm_cw * sm_ch;
            
            if (mat_next == NULL){
              dif_linear = 2 * d / length_perpen / (1 + 4 * d / length_perpen);
              dif_linear *= _shape_mesh->getBoundary(s);
              current += - dif_linear * flux * length;
            }
            else{
              d_next = mat_next->getDifCoefByGroup(g, position, temp_next);
              dif_linear = 2 * d * d_next / (length_perpen * d + length_perpen * d_next);
              flux_next = _shape_mesh->getFluxByValue(cell_next, g, position);
              current += - dif_linear * (flux - flux_next) * length;
            }
          }

          else if (s == 3 && x % num_refines_x == num_refines_x - 1){
            flux = _shape_mesh->getFluxByValue(*iter, g, position);
            d = mat->getDifCoefByGroup(g, position, temp);
            length_perpen = sm_cw;
            length = sm_ch * sm_cd;
            
            if (mat_next == NULL){
              dif_linear = 2 * d / length_perpen / (1 + 4 * d / length_perpen);
              dif_linear *= _shape_mesh->getBoundary(s);
              current += dif_linear * flux * length;
            }
            else{
              d_next = mat_next->getDifCoefByGroup(g, position, temp_next);
              dif_linear = 2 * d * d_next / (length_perpen * d + length_perpen * d_next);
              flux_next = _shape_mesh->getFluxByValue(cell_next, g, position);
              current += - dif_linear * (flux_next - flux) * length;
            }
          }            

          else if (s == 4 && y % num_refines_y == num_refines_y - 1){
            flux = _shape_mesh->getFluxByValue(*iter, g, position);
            d = mat->getDifCoefByGroup(g, position, temp);
            length_perpen = sm_ch;
            length = sm_cw * sm_cd;
            
            if (mat_next == NULL){
              dif_linear = 2 * d / length_perpen / (1 + 4 * d / length_perpen);
              dif_linear *= _shape_mesh->getBoundary(s);
              current += dif_linear * flux * length;
            }
            else{
              d_next = mat_next->getDifCoefByGroup(g, position, temp_next);
              dif_linear = 2 * d * d_next / (length_perpen * d + length_perpen * d_next);
              flux_next = _shape_mesh->getFluxByValue(cell_next, g, position);
              current += - dif_linear * (flux_next - flux) * length;
            }
          }

          else if (s == 5 && z % num_refines_z == num_refines_z - 1){
            flux = _shape_mesh->getFluxByValue(*iter, g, position);
            d = mat->getDifCoefByGroup(g, position, temp);
            length_perpen = sm_cd;
            length = sm_cw * sm_ch;
            
            if (mat_next == NULL){
              dif_linear = 2 * d / length_perpen / (1 + 4 * d / length_perpen);
              dif_linear *= _shape_mesh->getBoundary(s);
              current += dif_linear * flux * length;
            }
            else{
              d_next = mat_next->getDifCoefByGroup(g, position, temp_next);
              dif_linear = 2 * d * d_next / (length_perpen * d + length_perpen * d_next);
              flux_next = _shape_mesh->getFluxByValue(cell_next, g, position);
              current += - dif_linear * (flux_next - flux) * length;
            }            
          }
        }
        
        _current[position][(i*6+s) * _num_energy_groups + g] = current;
      }
    }
  }
}


double AmpMesh::computeDifCorrect(double dif_coef, double length){

  double f = 0.0;
  
  if (_optically_thick){
    double mu = cos(asin(0.798184));
    double expon = exp(- length / (3 * dif_coef * mu));
    double alpha = (1 + expon) / (1 - expon) - 2 * 3 * dif_coef * mu / length;
    double rho = mu * alpha;
    f = 1 + length * rho / (2 * dif_coef);
  }
  else{
    f = 1.0;
  }
  
  return f;
}


void AmpMesh::computeDifCoefs(int position){

  int nx = _num_x;
  int ny = _num_y;
  int nz = _num_z;
  int ng = _num_energy_groups;
  double width = getCellWidth();
  double height = getCellHeight();
  double depth = getCellDepth();
  double* temps = _temperature[position];
  

  for (int z=0; z < nz; z++){
    #pragma omp parallel for  
    for (int y=0; y < ny; y++){
      for (int x=0; x < nx; x++){

        int sense;
        double length, length_perpen;
        double dif_coef, current, flux, f, dif_linear, dif_nonlinear, flux_next, dif_coef_next, f_next;
        int cell = z*nx*ny+y*nx+x;
        double temp = temps[cell];
        
        for (int s=0; s < 6; s++){
          
          int cell_next = getNeighborCell(x, y, z, s);
          
          if (s == 0 || s == 1 || s == 2)
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
           length = width * height;
           length_perpen = depth;
          }
          
          for (int g=0; g < ng; g++){
            
            dif_coef = _materials[cell]->getDifCoefByGroup(g, position, temp);
            current = getCurrentByValue(cell, g, s, position);
            flux = getFluxByValue(cell, g, position);
            f = computeDifCorrect(dif_coef, length_perpen);
            
            if (cell_next == -1){
              dif_linear = 2 * dif_coef * f / length_perpen / (1 + 4 * dif_coef * f / length_perpen);
              dif_nonlinear = (sense * dif_linear * flux - current / length) / flux;
              dif_linear *= _boundaries[s];
              dif_nonlinear *= _boundaries[s];
            }
            else{
              flux_next = getFluxByValue(cell_next, g, position);
              dif_coef_next = _materials[cell_next]->getDifCoefByGroup(g, position, temps[cell_next]);
              f_next = computeDifCorrect(dif_coef_next, length_perpen);
              dif_linear = 2 * dif_coef * f * dif_coef_next * f_next / (length_perpen * dif_coef * f +
                                                                        length_perpen * dif_coef_next * f_next);
              dif_nonlinear = - (sense * dif_linear * (flux_next - flux) + current / length)
                / (flux_next + flux);
              
              if (dif_nonlinear > dif_linear){
                if (sense == -1){
                  if (dif_nonlinear > 0.0){
                    dif_linear = - current / (2 * flux);
                    dif_nonlinear = - current / (2 * flux);
                  }
                  else{
                    dif_linear = current / (2 * flux_next);
                    dif_nonlinear = - current / (2 * flux_next);                  
                  }
                }
                else{
                  if (dif_nonlinear > 0.0){
                    dif_linear = - current / (2 * flux_next);
                    dif_nonlinear = - current / (2 * flux_next);
                  }
                  else{
                    dif_linear = current / (2 * flux);
                    dif_nonlinear = - current / (2 * flux);
                  }
                }
              }
            }
            
            setDifLinearByValue(dif_linear, cell, g, s, position);
            setDifNonlinearByValue(dif_nonlinear, cell, g, s, position);          
          }
        }
      }
    }
  }
}


AmpMesh* AmpMesh::clone(){
  
  AmpMesh* mesh = new AmpMesh(getWidth(), getHeight(), getDepth(), _num_x, _num_y, _num_z);

  mesh->setNumShapeEnergyGroups(_num_shape_energy_groups);
  mesh->setNumAmpEnergyGroups(_num_energy_groups);
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
}


void AmpMesh::setGroupStructure(int* group_indices, int length_group_indices){

  /* Allocate memory */
  if (_group_indices == NULL){
    _group_indices = new int[length_group_indices];
  }

  if (group_indices == NULL){
    for (int i = 0; i < length_group_indices; i++){
      _group_indices[i] = i;
    }
  }
  else{
    if (group_indices[0] != 0)
      log_printf(ERROR, "The first value in group indices must be 1!");

    /* Set first group indice to 0 */
    _group_indices[0] = 0;

    /* Set MOC group bounds for rest of CMFD energy groups */
    for (int i = 1; i < length_group_indices; i++){
      /* Check that the group indices are always increasing */
      if (group_indices[i] <= group_indices[i-1])
        log_printf(ERROR, "The group indices must be increasing!");

      _group_indices[i] = group_indices[i];
    }
  }
}


double AmpMesh::computeAveragePower(int position){

  if (_fuel_volume == 0.0)
    computeFuelVolume();   
  
  computePower(position);

  return pairwise_sum(_power[position], _num_x*_num_y*_num_z) * getCellVolume() / _fuel_volume;
}  


double AmpMesh::computePowerL2Norm(int position_1, int position_2){

  computePower(position_1);
  computePower(position_2);

  double* power_residual = new double[_num_x * _num_y * _num_z];
  memset(power_residual, 0.0, sizeof(double) * _num_x * _num_y * _num_z);

  #pragma omp parallel for
  for (int i=0; i < _num_x * _num_y * _num_z; i++){
    if (_power[position_1][i] > 0.0)
      power_residual[i] = pow((_power[position_1][i] - _power[position_2][i]) / _power[position_1][i], 2);
  }

  double residual = sqrt(pairwise_sum(power_residual, _num_x*_num_y*_num_z));
  delete [] power_residual;

  return residual;
}


void AmpMesh::interpolateDifNonlinear(int position_begin, int position_end, int position){

  double dt = _clock->getTime(position_end) - _clock->getTime(position_begin);
  double wt_begin = (_clock->getTime(position_end) - _clock->getTime(position)) / dt;
  double wt_end = (_clock->getTime(position) - _clock->getTime(position_begin)) / dt;

  #pragma omp parallel for
  for (int i=0; i < _num_x*_num_y*_num_z*_num_energy_groups*6; i++){
    _dif_nonlinear[position][i] = _dif_nonlinear[position_begin][i] * wt_begin
      + _dif_nonlinear[position_end][i] * wt_end;
  }
}
