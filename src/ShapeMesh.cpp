#include "ShapeMesh.h"


/**
 * @brief Constructor sets the ID and unique ID for the Material.
 * @param id the user-defined ID for the material
 */
ShapeMesh::ShapeMesh(double width, double height, double depth) : 
  Mesh(width, height, depth){

  _amp_mesh = NULL;
  _group_indices = NULL;
  _amp_map = NULL;
  
  return;
}


/**
 * @brief Destructor deletes all cross-section data structures from memory.
 */
ShapeMesh::~ShapeMesh() {

  if (_group_indices != NULL)
    delete [] _group_indices;

  if (_amp_map != NULL)
    delete [] _amp_map;
}


void ShapeMesh::setAmpMesh(AmpMesh* mesh){
  _amp_mesh = mesh;
}


void ShapeMesh::setFluxByValue(double flux, int cell, int group, int position){

  if (cell >= getNumCells() || cell < 0)
    log_printf(ERROR, "Unable to set flux by value for Shape Mesh "
               "cell %i since there are only %i cells", cell, getNumCells());

  if (group >= _num_shape_energy_groups || group < 0)
    log_printf(ERROR, "Unable to set flux by value for Shape Mesh group "
               "%i since there are only %i groups", group, _num_shape_energy_groups);

  if (position >= 8 || position < 0)
    log_printf(ERROR, "Unable to set flux by value for Shape Mesh clock "
               "position %i since there are only 8 clock positions", position);

  if (flux < 0.0)
    log_printf(NORMAL, "Unable to set flux by value for Shape Mesh "
               "for cell %i and group %i to a negative value %f", cell, group, flux);

  _flux[position][cell*_num_shape_energy_groups + group] = flux;
}


double ShapeMesh::getFluxByValue(int cell, int group, int position){

  if (cell >= getNumCells() || cell < 0)
    log_printf(ERROR, "Unable to get flux by value for Shape Mesh "
               "cell %i since there are only %i cells", cell, getNumCells());

  if (group >= _num_shape_energy_groups || group < 0)
    log_printf(ERROR, "Unable to get flux by value for Shape Mesh group "
               "%i since there are only %i groups", group, _num_shape_energy_groups);

  if (position >= 8 || position < 0)
    log_printf(ERROR, "Unable to get flux by value for Shape Mesh clock "
               "position %i since there are only 8 clock positions", position);

  if (_flux[position] == NULL)
    log_printf(NORMAL, "Unable to get flux by value for Shape Mesh "
               "since the flux has not been initialized");

  return _flux[position][cell*_num_shape_energy_groups + group];
}


void ShapeMesh::copyFlux(int position_from, int position_to){

  if (position_from >= 8 || position_from < 0)
    log_printf(ERROR, "Unable to copy flux from position %i since "
               "there are only 8 clock positions", position_from);

  if (position_to >= 8 || position_to < 0)
    log_printf(ERROR, "Unable to copy flux to position %i since "
               "there are only 8 clock positions", position_to);
    
  if (_flux[position_from] == NULL)
    log_printf(ERROR, "Unable to copy flux from position %i since "
               "the flux have not been initialized yet", position_from);

  if (_flux[position_to] == NULL)
    log_printf(ERROR, "Unable to copy flux to position %i since "
               "the flux have not been initialized yet", position_to);

  std::copy(_flux[position_from], 
            _flux[position_from] + getNumCells()*_num_shape_energy_groups, 
            _flux[position_to]);
}


void ShapeMesh::initialize(){

  Mesh::initialize();

  for (int i=0; i < getNumCells(); i++){
    Material* material = new Material();
    material->setNumEnergyGroups(_num_amp_energy_groups);
    material->setNumDelayedGroups(_num_delayed_groups);
    _materials[i] = material;
  }
}


void ShapeMesh::synthesizeFlux(int position){

  if (position >= 8 || position < 0)
    log_printf(ERROR, "Unable to synthesize flux at position %i since "
               "there are only 8 clock positions", position);

  if (_num_shape_energy_groups == 0)
    log_printf(ERROR, "Unable to synthesize the flux since the number "
               "of shape energy groups has not been set");

  if (_num_amp_energy_groups == 0)
    log_printf(ERROR, "Unable to synthesize the flux since the number "
               "of amp energy groups has not been set");

  if (_num_delayed_groups == 0)
    log_printf(ERROR, "Unable to synthesize the flux since the number "
               "of delayed groups has not been set");

  int ngs = _num_shape_energy_groups;
  int nga = _num_amp_energy_groups;
  double dt = _clock->getTime(FORWARD_OUT) - _clock->getTime(PREVIOUS_OUT);
  double wt_begin = (_clock->getTime(FORWARD_OUT) - _clock->getTime(position)) / dt;
  double wt_end = (_clock->getTime(position) - _clock->getTime(PREVIOUS_OUT)) / dt;

  #pragma omp parallel for
  for (int i=0; i < getNumCells(); i++){
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


void ShapeMesh::reconstructFlux(int position, int position_shape, int position_amp){

  if (position >= 8 || position < 0)
    log_printf(ERROR, "Unable to reconstruct flux at position %i since "
               "there are only 8 clock positions", position);

  if (_num_shape_energy_groups == 0)
    log_printf(ERROR, "Unable to reconstruct the flux since the number "
               "of shape energy groups has not been set");

  if (_num_amp_energy_groups == 0)
    log_printf(ERROR, "Unable to reconstruct the flux since the number "
               "of amp energy groups has not been set");

  if (_num_delayed_groups == 0)
    log_printf(ERROR, "Unable to reconstruct the flux since the number "
               "of delayed groups has not been set");

  int ngs = _num_shape_energy_groups;
  int nga = _num_amp_energy_groups;

  #pragma omp parallel for
  for (int i=0; i < getNumCells(); i++){
    double shape;
    for (int g=0; g < ngs; g++){
      shape = _flux[position_shape][i*ngs+g]
        / _amp_mesh->getFlux(position_shape)[_amp_map[i]*nga+_group_indices[g]];
      _flux[position][i*ngs+g] = shape *
        _amp_mesh->getFlux(position_amp)[_amp_map[i]*nga+_group_indices[g]];
    }
  }
}


void ShapeMesh::computePower(int position){

  if (position >= 8 || position < 0)
    log_printf(ERROR, "Unable to compute the power at position %i since "
               "there are only 8 clock positions", position);

  if (_num_shape_energy_groups == 0)
    log_printf(ERROR, "Unable to compute the power since the number "
               "of shape energy groups has not been set");

  #pragma omp parallel for
  for (int i=0; i < getNumCells(); i++){
    double fission_rate = 0.0;
    double temp = _temperature[position][i];
    
    for (int g=0; g < _num_shape_energy_groups; g++){
      fission_rate += _materials[i]->getSigmaFByGroup(g, position, temp) * getFluxByValue(i, g, position);
    }
    
    _power[position][i] = fission_rate * _materials[i]->getEnergyPerFission();
  }
}


void ShapeMesh::computeInitialPrecursorConc(int position){

  if (position >= 8 || position < 0)
    log_printf(ERROR, "Unable to compute the initial precursor conc at position %i since "
               "there are only 8 clock positions", position);

  if (_num_shape_energy_groups == 0)
    log_printf(ERROR, "Unable to compute the initial precursor conc since the number "
               "of shape energy groups has not been set");

  if (_num_delayed_groups == 0)
    log_printf(ERROR, "Unable to compute the initial precursor conc since the number "
               "of delayed groups has not been set");

  double* temps = _temperature[position];

  #pragma omp parallel for
  for (int i=0; i < getNumCells(); i++){
    
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


void ShapeMesh::integratePrecursorConc(int position_from, int position_to){

  if (position_from >= 8 || position_from < 0)
    log_printf(ERROR, "Unable to integrate the precursor conc from position %i since "
               "there are only 8 clock positions", position_from);

  if (position_to >= 8 || position_to < 0)
    log_printf(ERROR, "Unable to integrate the precursor conc to position %i since "
               "there are only 8 clock positions", position_to);

  if (_num_shape_energy_groups == 0)
    log_printf(ERROR, "Unable to integrate the precursor conc since the number "
               "of shape energy groups has not been set");

  if (_num_delayed_groups == 0)
    log_printf(ERROR, "Unable to integrate the precursor conc since the number "
               "of delayed groups has not been set");

  double* temps_from = _temperature[position_from];
  double* temps_to = _temperature[position_to];
  double dt = _clock->getTime(position_to) - _clock->getTime(position_from);

  #pragma omp parallel for
  for (int i=0; i < getNumCells(); i++){

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
          * fission_rate_to - k3 * _delayed_fractions[d] / _decay_constants[d] / _k_eff_0 
          * fission_rate_from;

        
        _materials[i]->setPrecursorConcByGroup(new_conc, d, position_to);
      }
    }
  }
}


void ShapeMesh::integrateTemperature(int position_from, int position_to){

  if (position_from >= 8 || position_from < 0)
    log_printf(ERROR, "Unable to integrate the temperature from position %i since "
               "there are only 8 clock positions", position_from);

  if (position_to >= 8 || position_to < 0)
    log_printf(ERROR, "Unable to integrate the temperature to position %i since "
               "there are only 8 clock positions", position_to);

  if (_num_shape_energy_groups == 0)
    log_printf(ERROR, "Unable to integrate the temperature since the number "
               "of shape energy groups has not been set");

  if (_num_delayed_groups == 0)
    log_printf(ERROR, "Unable to integrate the temperature since the number "
               "of delayed groups has not been set");

  double* temps_from = _temperature[position_from];
  double* temps_to = _temperature[position_to];
  double dt = _clock->getTime(position_to) - _clock->getTime(position_from);

  #pragma omp parallel for
  for (int i=0; i < getNumCells(); i++){

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


void ShapeMesh::setGroupStructure(int* group_indices, int length_group_indices){

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


void ShapeMesh::scaleFlux(int position, double scale_val){

  if (position >= 8 || position < 0)
    log_printf(ERROR, "Unable to scale the flux at position %i since "
               "there are only 8 clock positions", position);

  if (_num_shape_energy_groups == 0)
    log_printf(ERROR, "Unable to scale the flux since the number "
               "of shape energy groups has not been set");

  #pragma omp parallel for
  for (int i=0; i < getNumCells()*_num_shape_energy_groups; i++)
    _flux[position][i] *= scale_val;
}


double ShapeMesh::computeAveragePower(int position){

  if (_fuel_volume == 0.0)
    computeFuelVolume();   

  computePower(position);

  double* power = new double[getNumCells()];
  
  for (int i=0; i < getNumCells(); i++){
      power[i] = _power[position][i] * getCellVolume(i);
  }

  double average_power = pairwise_sum(power, getNumCells()) / _fuel_volume;
  delete [] power;

  return average_power;
}  


double ShapeMesh::computePowerL2Norm(int position_1, int position_2){

  computePower(position_1);
  computePower(position_2);

  double* power_residual = new double[getNumCells()];
  memset(power_residual, 0.0, sizeof(double) * getNumCells());

  #pragma omp parallel for
  for (int i=0; i < getNumCells(); i++){
    if (_power[position_1][i] > 0.0){
      power_residual[i] = pow((_power[position_1][i] - _power[position_2][i]) / _power[position_1][i], 2);
    }
  }
  
  double residual = sqrt(pairwise_sum(power_residual, getNumCells()));
  delete [] power_residual;

  return residual;
}


int ShapeMesh::getAmpGroup(int shape_group){
  
  if (shape_group < 0 || shape_group >= _num_shape_energy_groups)
    log_printf(ERROR, "Unable to get amp group for shape group %i since there "
               "are only %i shape energy groups", 
               shape_group, _num_shape_energy_groups);

  return _group_indices[shape_group];
}


void ShapeMesh::setFlux(double* flux, int num_cells_times_groups){

  if (num_cells_times_groups != getNumCells() * _num_shape_energy_groups)
    log_printf(ERROR, "Unable to set flux with %i values when the number of "
               "cells is %i and shape energy groups is %i", _num_cells,
               _num_shape_energy_groups);

  for (int c=0; c < 8; c++){
    if (_flux[c] == NULL)
      log_printf(ERROR, "Unable to set current for Shape Mesh since the"
                 " flux has not been initialized");
    
    std::copy(flux, flux + num_cells_times_groups, _flux[c]);
  }  
}


void ShapeMesh::setAmpCellContainingShapeCell(int shape_cell, int amp_cell){

  /* Check if shape cell is valid */
  if (shape_cell >= _num_cells)
    log_printf(ERROR, "Unable to set amp cell containing shape cell for shape cell %i"
               " since there are only %i shape cells", shape_cell, _num_cells);

  /* Allocate memory for amp map if it has not been initialized */
  if (_amp_map == NULL)
    _amp_map = new int[_num_cells];
    
  /* Set amp cell containing shape cell */
  _amp_map[shape_cell] = amp_cell;
}
