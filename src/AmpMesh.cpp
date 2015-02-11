#include "AmpMesh.h"


/**
 * @brief Constructor sets the ID and unique ID for the Material.
 * @param id the user-defined ID for the material
 */
AmpMesh::AmpMesh(double width, double height, int num_x, int num_y) :
  StructuredMesh(width, height, num_x, num_y){

  _shape_mesh = NULL;
  _optically_thick = false;
  _group_indices = NULL;
  _energy_per_fission = 0.0;
  
  return;
}


/**
 * @brief Destructor deletes all cross-section data structures from memory.
 */
AmpMesh::~AmpMesh() {

  if (_group_indices != NULL)
    delete [] _group_indices;
}


void AmpMesh::setEnergyPerFission(double energy_per_fission){
  _energy_per_fission = energy_per_fission;
}

void AmpMesh::setOpticallyThick(bool optically_thick){
  _optically_thick = optically_thick;
}


void AmpMesh::setShapeMesh(StructuredShapeMesh* mesh){
  _shape_mesh = mesh;

  for (int i=0; i < _num_x*_num_y; i++){
    std::vector<int> *shape_cells = new std::vector<int>;
    _shape_map.push_back(*shape_cells);
  }
  
  int nx = mesh->getNumX();
  int ny = mesh->getNumY();
  int num_refines_x = nx / _num_x;
  int num_refines_y = ny / _num_y;
  
  for (int y=0; y < ny; y++){
    int yy = y / num_refines_y;
    for (int x=0; x < nx; x++){
      int xx = x / num_refines_x;
      int amp_cell = yy*_num_x + xx;
      int shape_cell = y*nx + x;
      _shape_map[amp_cell].push_back(shape_cell);
    }
  }
}


void AmpMesh::setFluxByValue(double flux, int cell, int group, int position){
  _flux[position][cell*_num_amp_energy_groups + group] = flux;
}


void AmpMesh::setCurrentByValue(double current, int cell, int group, int side, int position){
  _current[position][(cell*4 + side)*_num_amp_energy_groups + group] = current;
}


void AmpMesh::setDifLinearByValue(double dif_linear, int cell, int group, int side, int position){
  _dif_linear[position][(cell*4 + side)*_num_amp_energy_groups + group] = dif_linear;
}


void AmpMesh::setDifNonlinearByValue(double dif_nonlinear, int cell, int group, int side, int position){
  _dif_nonlinear[position][(cell*4 + side)*_num_amp_energy_groups + group] = dif_nonlinear;
}


double AmpMesh::getFluxByValue(int cell, int group, int position){
  return _flux[position][cell*_num_amp_energy_groups + group];
}


double AmpMesh::getCurrentByValue(int cell, int group, int side, int position){
  return _current[position][(cell*4 + side)*_num_amp_energy_groups + group];
}


double AmpMesh::getDifLinearByValue(int cell, int group, int side, int position){
  return _dif_linear[position][(cell*4 + side)*_num_amp_energy_groups + group];
}


double AmpMesh::getDifNonlinearByValue(int cell, int group, int side, int position){
  return _dif_nonlinear[position][(cell*4 + side)*_num_amp_energy_groups + group];
}


void AmpMesh::copyFlux(int position_from, int position_to){
  std::copy(_flux[position_from], _flux[position_from] + _num_x*_num_y*_num_amp_energy_groups, _flux[position_to]);
}


void AmpMesh::copyCurrent(int position_from, int position_to){
  std::copy(_current[position_from], _current[position_from] + _num_x*_num_y*_num_amp_energy_groups*4, _current[position_to]);
}


void AmpMesh::copyDifLinear(int position_from, int position_to){
  std::copy(_dif_linear[position_from], _dif_linear[position_from] + _num_x*_num_y*_num_amp_energy_groups*4, _dif_linear[position_to]);
}


void AmpMesh::copyDifNonlinear(int position_from, int position_to){
  std::copy(_dif_nonlinear[position_from], _dif_nonlinear[position_from] + _num_x*_num_y*_num_amp_energy_groups*4, _dif_nonlinear[position_to]);
}


void AmpMesh::condenseMaterials(int position, bool save_flux){

  double shape_cell_volume = _shape_mesh->getCellVolume();
  double amp_cell_volume = getCellVolume();
  double* temps = _shape_mesh->getTemperature(position);
  double sigma_a, sigma_t, sigma_f, nu_sigma_f, dif_coef, rxn, production, velocity, chi_tally;
  double* chi = new double[_num_amp_energy_groups];
  double* sigma_s = new double[_num_amp_energy_groups];
  double sigma_t_group, rxn_group, volume;
  Material* amp_mat;
  Material* shape_mat;
  double shape;
  double precursor_conc, temp;
  
  for (int i=0; i < _num_x*_num_y; i++){
    
    amp_mat = getMaterial(i);
    std::vector<int>::iterator iter;
    
    for (int g=0; g < _num_amp_energy_groups; g++){

      sigma_a = 0.0;
      sigma_t = 0.0;
      sigma_f = 0.0;
      nu_sigma_f = 0.0;
      dif_coef = 0.0;
      rxn = 0.0;
      production = 0.0;
      velocity = 0.0;

      for (int h=0; h < _num_amp_energy_groups; h++){
        sigma_s[h] = 0.0;
        chi[h] = 0.0;
      }
      
      for (iter = _shape_map[i].begin(); iter != _shape_map[i].end(); ++iter){

        shape_mat = _shape_mesh->getMaterial(*iter);
        temp = temps[*iter];

        for (int e=0; e < _num_amp_energy_groups; e++){
          chi_tally = 0.0;

          for (int ee=_group_indices[e]; ee < _group_indices[e+1]; ee++)
            chi_tally += shape_mat->getChiByGroup(ee, position, temp);
          
          for (int h=0; h < _num_shape_energy_groups; h++){
            chi[e] += chi_tally * shape_mat->getNuSigmaFByGroup(h, position, temp) *
              _shape_mesh->getFluxByValue(*iter, h, position) * shape_cell_volume;
            production += chi_tally * shape_mat->getNuSigmaFByGroup(h, position, temp) *
              _shape_mesh->getFluxByValue(*iter, h, position) * shape_cell_volume;
          }
        }
      }

      for (int gg=_group_indices[g]; gg < _group_indices[g+1]; gg++){

        rxn_group = 0.0;
        volume = 0.0;
        sigma_t_group = 0.0;

      for (iter = _shape_map[i].begin(); iter != _shape_map[i].end(); ++iter){

          shape_mat = _shape_mesh->getMaterial(*iter);
          temp = temps[*iter];

          shape = _shape_mesh->getFluxByValue(*iter, gg, position);
          sigma_a += shape_mat->getSigmaAByGroup(gg, position, temp) * shape * shape_cell_volume;
          sigma_t += shape_mat->getSigmaTByGroup(gg, position, temp) * shape * shape_cell_volume;
          sigma_f += shape_mat->getSigmaFByGroup(gg, position, temp) * shape * shape_cell_volume;
          nu_sigma_f += shape_mat->getNuSigmaFByGroup(gg, position, temp) * shape * shape_cell_volume;
          rxn += shape * shape_cell_volume;
          rxn_group += shape * shape_cell_volume;
          volume += shape_cell_volume;
          sigma_t_group += shape_mat->getSigmaTByGroup(gg, position, temp) * shape * shape_cell_volume;
          velocity += 1.0 / shape_mat->getVelocityByGroup(gg, position, temp) * shape * shape_cell_volume;

          for (int h=0; h < _num_shape_energy_groups; h++){
            sigma_s[_shape_mesh->getAmpGroup(h)] += shape_mat->getSigmaSByGroup(gg, h, position, temp) *
              shape * shape_cell_volume;
          }
        }

        dif_coef += rxn_group / (3.0 * sigma_t_group / rxn_group);
      }

      amp_mat->setSigmaAByGroup(sigma_a / rxn, g, position);
      amp_mat->setSigmaTByGroup(sigma_t / rxn, g, position);
      amp_mat->setSigmaFByGroup(sigma_f / rxn, g, position);
      amp_mat->setNuSigmaFByGroup(nu_sigma_f / rxn, g, position);
      amp_mat->setDifCoefByGroup(dif_coef / rxn, g, position);
      amp_mat->setVelocityByGroup(1.0 / (velocity / rxn), g, position);

      if (save_flux)
        _flux[position][i*_num_amp_energy_groups + g] = rxn / volume;

      if (production != 0.0)
        amp_mat->setChiByGroup(chi[g] / production, g, position);
      else
        amp_mat->setChiByGroup(0.0, g, position);      

      for (int e=0; e < _num_amp_energy_groups; e++)
        amp_mat->setSigmaSByGroup(sigma_s[e] / rxn, g, e, position);
    }

    for (int d=0; d < _num_delayed_groups; d++){
      precursor_conc = 0.0;

      for (iter = _shape_map[i].begin(); iter != _shape_map[i].end(); ++iter){
        if (_shape_mesh->getMaterial(*iter)->isFissionable()){
          precursor_conc += _shape_mesh->getMaterial(*iter)->getPrecursorConcByGroup(d, position) * shape_cell_volume;
        }
      }
        
      amp_mat->setPrecursorConcByGroup(precursor_conc / amp_cell_volume, d, position);
    }
  }

  delete [] chi;
  delete [] sigma_s;
}


void AmpMesh::initialize(){

  _materials = new Material*[_num_x*_num_y];
  
  for (int i=0; i < _num_x*_num_y; i++){
    Material* material = new Material();
    material->setNumEnergyGroups(_num_amp_energy_groups);
    material->setNumDelayedGroups(_num_delayed_groups);
    _materials[i] = material;
  }

  for(int c=0; c < 8; c++){
    _dif_linear[c] = new double[_num_x * _num_y * _num_amp_energy_groups * 4];
    _dif_nonlinear[c] = new double[_num_x * _num_y * _num_amp_energy_groups * 4];
    _current[c] = new double[_num_x * _num_y * _num_amp_energy_groups * 4];
    _flux[c] = new double[_num_x * _num_y * _num_amp_energy_groups];
    _temperature[c] = new double[_num_x * _num_y];
    _power[c] = new double[_num_x * _num_y];

    memset(_dif_linear[c], 0.0, sizeof(double) * _num_x * _num_y * _num_amp_energy_groups * 4);
    memset(_dif_nonlinear[c], 0.0, sizeof(double) * _num_x * _num_y * _num_amp_energy_groups * 4);
    memset(_current[c], 0.0, sizeof(double) * _num_x * _num_y * _num_amp_energy_groups * 4);
    memset(_flux[c], 1.0, sizeof(double) * _num_x * _num_y * _num_amp_energy_groups);
    memset(_temperature[c], 300.0, sizeof(double) * _num_x * _num_y);
    memset(_power[c], 0.0, sizeof(double) * _num_x * _num_y);
  }  
}


void AmpMesh::computePower(int position){

  for (int i=0; i < _num_x * _num_y; i++){
    double fission_rate = 0.0;
    double temp = _temperature[position][i];
    
    for (int g=0; g < _num_amp_energy_groups; g++)
      fission_rate += _materials[i]->getSigmaFByGroup(g, position, temp) * getFluxByValue(i, g, position);
    
    _power[position][i] = fission_rate * _energy_per_fission;
  }
}


void AmpMesh::computeCurrent(int position){

  int sm_nx = _shape_mesh->getNumX();
  int sm_cw = _shape_mesh->getCellWidth();
  int sm_ch = _shape_mesh->getCellHeight();
  int num_refines = sm_nx / _num_x;
  double* temps = _shape_mesh->getTemperature(position);
  double current, flux, dif_linear, d, d_next, temp, temp_next, flux_next;
  Material *mat, *mat_next;
  int cell_next;
  double length_perpen, length;
  
  for (int i=0; i < _num_x * _num_y; i++){
    std::vector<int>::iterator iter;
    for (int g=0; g < _num_amp_energy_groups; g++){
      for (int s=0; s < 4; s++){
        current = 0.0;
        for (iter = _shape_map[i].begin(); iter != _shape_map[i].end(); ++iter){

          mat = _shape_mesh->getMaterial(*iter);
          temp = temps[*iter];
          cell_next = _shape_mesh->getNeighborCell(*iter % sm_nx, *iter / sm_nx, s);
          temp_next = temps[cell_next];
          mat_next = _shape_mesh->getNeighborMaterial(*iter % sm_nx, *iter / sm_nx, s);

          if (s == 0 && (*iter % sm_nx) % num_refines == 0){
            for (int gg=_group_indices[g]; gg < _group_indices[g+1]; gg++){
              flux = _shape_mesh->getFluxByValue(*iter, gg, position);
              d = mat->getDifCoefByGroup(gg, position, temp);
              length_perpen = sm_cw;
              length = sm_ch;

              if (mat_next == NULL){
                dif_linear = 2 * d / length_perpen / (1 + 4 * d / length_perpen);
                dif_linear *= _shape_mesh->getBoundary(s);
                current += - dif_linear * flux * length;
              }
              else{
                d_next = mat_next->getDifCoefByGroup(gg, position, temp_next);
                dif_linear = 2 * d * d_next / (length_perpen * d + length_perpen * d_next);
                flux_next = _shape_mesh->getFluxByValue(cell_next, gg, position);
                current += - dif_linear * (flux - flux_next) * length;
              }
            }
          }
          else if (s == 1 && (*iter / sm_nx) % num_refines == 0){
            for (int gg=_group_indices[g]; gg < _group_indices[g+1]; gg++){
              flux = _shape_mesh->getFluxByValue(*iter, gg, position);
              d = mat->getDifCoefByGroup(gg, position, temp);
              length_perpen = sm_ch;
              length = sm_cw;

              if (mat_next == NULL){
                dif_linear = 2 * d / length_perpen / (1 + 4 * d / length_perpen);
                dif_linear *= _shape_mesh->getBoundary(s);
                current += - dif_linear * flux * length;
              }
              else{
                d_next = mat_next->getDifCoefByGroup(gg, position, temp_next);
                dif_linear = 2 * d * d_next / (length_perpen * d + length_perpen * d_next);
                flux_next = _shape_mesh->getFluxByValue(cell_next, gg, position);
                current += - dif_linear * (flux - flux_next) * length;
              }
            }            
          }
          else if (s == 2 && (*iter % sm_nx) % num_refines == num_refines - 1){
            for (int gg=_group_indices[g]; gg < _group_indices[g+1]; gg++){
              flux = _shape_mesh->getFluxByValue(*iter, gg, position);
              d = mat->getDifCoefByGroup(gg, position, temp);
              length_perpen = sm_cw;
              length = sm_ch;

              if (mat_next == NULL){
                dif_linear = 2 * d / length_perpen / (1 + 4 * d / length_perpen);
                dif_linear *= _shape_mesh->getBoundary(s);
                current += dif_linear * flux * length;
              }
              else{
                d_next = mat_next->getDifCoefByGroup(gg, position, temp_next);
                dif_linear = 2 * d * d_next / (length_perpen * d + length_perpen * d_next);
                flux_next = _shape_mesh->getFluxByValue(cell_next, gg, position);
                current += - dif_linear * (flux_next - flux) * length;
              }
            }            
          }
          else if (s == 3 && (*iter / sm_nx) % num_refines == num_refines - 1){
            for (int gg=_group_indices[g]; gg < _group_indices[g+1]; gg++){
              flux = _shape_mesh->getFluxByValue(*iter, gg, position);
              d = mat->getDifCoefByGroup(gg, position, temp);
              length_perpen = sm_ch;
              length = sm_cw;

              if (mat_next == NULL){
                dif_linear = 2 * d / length_perpen / (1 + 4 * d / length_perpen);
                dif_linear *= _shape_mesh->getBoundary(s);
                current += dif_linear * flux * length;
              }
              else{
                d_next = mat_next->getDifCoefByGroup(gg, position, temp_next);
                dif_linear = 2 * d * d_next / (length_perpen * d + length_perpen * d_next);
                flux_next = _shape_mesh->getFluxByValue(cell_next, gg, position);
                current += - dif_linear * (flux_next - flux) * length;
              }
            }            
          }
        }

        _current[position][(i*4+s) * _num_amp_energy_groups + g] = current;
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
  int ng = _num_amp_energy_groups;
  double width = getCellWidth();
  double height = getCellHeight();
  double* temps = _temperature[position];
  int sense;
  double length, length_perpen;
  double dif_coef, current, flux, f, dif_linear, dif_nonlinear, flux_next, dif_coef_next, f_next;
  
  
  for (int x=0; x < nx; x++){
    for (int y=0; y < ny; y++){
      int cell = y*nx+x;
      double temp = temps[cell];

      for (int s=0; s < 4; s++){

        int cell_next = getNeighborCell(x, y, s);

        if (s == 0 || s == 1)
          sense = -1;
        else
          sense = 1;

        if (s == 0 || s == 2){
          length = height;
          length_perpen = width;
        }
        else{
          length = width;
          length_perpen = height;
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


AmpMesh* AmpMesh::clone(){
  
  AmpMesh* mesh = new AmpMesh(getWidth(), getHeight(), _num_x, _num_y);

  mesh->setNumShapeEnergyGroups(_num_shape_energy_groups);
  mesh->setNumAmpEnergyGroups(_num_amp_energy_groups);
  mesh->setNumDelayedGroups(_num_delayed_groups);
  mesh->setBuckling(_buckling);
  mesh->setKeff0(_k_eff_0);
 
  for (int s=0; s < 4; s++)
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

  return pairwise_sum(_power[position], _num_x*_num_y) * getCellVolume() / _fuel_volume;
}  


double AmpMesh::computePowerL2Norm(int position_1, int position_2){

  computePower(position_1);
  computePower(position_2);

  double* power_residual = new double[_num_x * _num_y];
  memset(power_residual, 0.0, sizeof(double) * _num_x * _num_y);
  
  for (int i=0; i < _num_x * _num_y; i++){
    if (_power[position_1][i] > 0.0)
      power_residual[i] = pow((_power[position_1][i] - _power[position_2][i]) / _power[position_1][i], 2);
  }


  double residual = sqrt(pairwise_sum(power_residual, _num_x*_num_y));
  delete [] power_residual;

  return residual;
}


void AmpMesh::interpolateDifNonlinear(int position_begin, int position_end, int position){

  double dt = _clock->getTime(position_end) - _clock->getTime(position_begin);
  double wt_begin = (_clock->getTime(position_end) - _clock->getTime(position)) / dt;
  double wt_end = (_clock->getTime(position) - _clock->getTime(position_begin)) / dt;

  for (int i=0; i < _num_x*_num_y*_num_amp_energy_groups*4; i++){
    _dif_nonlinear[position][i] = _dif_nonlinear[position_begin][i] * wt_begin
      + _dif_nonlinear[position_end][i] * wt_end;
  }
}
