#include "Solver.h"

Solver::Solver(Geometry* geometry){

  /* Initialize variables */
  _k_eff_0 = 0.0;
  _geometry = geometry;
  _num_energy_groups = _geometry->getNumEnergyGroups();
  _num_delayed_groups = _geometry->getNumDelayedGroups();
  _num_cells = _geometry->getNumShapeCells();
  _clock = new Clock();
  _buckling = 0.0;
  _initial_power = 1.e-6;
  _flux_solve_tolerance = FLUX_SOLVE_TOLERANCE;

  /* Initialize coarse mesh matrices */
  int num_x = _geometry->getNumXAmp();
  int num_y = _geometry->getNumYAmp();
  int num_z = _geometry->getNumZAmp();
  _AM = new Matrix(num_x, num_y, num_z, _num_energy_groups);
  _M  = new Matrix(num_x, num_y, num_z, _num_energy_groups);
  _A  = new Matrix(num_x, num_y, num_z, _num_energy_groups);
  _b  = new Vector(num_x, num_y, num_z, _num_energy_groups);
}


Solver::~Solver(){

  if (_AM != NULL)
    delete _AM;

  if (_M != NULL)
    delete _M;

  if (_A != NULL)
    delete _A;

  if (_b != NULL)
    delete _b;
}


double Solver::getKeff0(){
  return _k_eff_0;
}


double Solver::getBuckling(){
  return _buckling;
}


Geometry* Solver::getGeometry(){
  return _geometry;
}


Matrix* Solver::getAM(){
  return _AM;
}


Matrix* Solver::getA(){
  return _A;
}


Matrix* Solver::getM(){
  return _M;
}


Vector* Solver::getb(){
  return _b;
}


Vector* Solver::getTemperature(int state){
  return _temperature[state];
}


Vector* Solver::getFlux(int state){
  return _flux[state];
}


Vector* Solver::getPower(int state){
  return _power[state];
}


Vector* Solver::getWeight(int state){
  return _weight[state];
}


double Solver::getTemperatureByValue(int cell, int state){
  return _temperature[state]->getValue(cell, 0);
}


double Solver::getPowerByValue(int cell, int state){
  return _power[state]->getValue(cell, 0);
}


double Solver::getFluxByValue(int cell, int group, int state){
  return _flux[state]->getValue(cell, group);
}


double Solver::getWeightByValue(int cell, int group, int state){
  return _weight[state]->getValue(cell, group);
}


void Solver::setBuckling(double buckling){
  _buckling = buckling;
}


void Solver::setInitialPower(double power){
  _initial_power = power;
}


void Solver::setFluxSolveTolerance(double tolerance){
  _flux_solve_tolerance = tolerance;
}


void Solver::setEndTime(double time){
  _clock->setEndTime(time);
}

void Solver::setTemperatureByValue(double value, int cell, int state){
  _temperature[state]->setValue(cell, 0,  value);
}


void Solver::setFluxByValue(double value, int cell, int group, int state){
  _flux[state]->setValue(cell, group, value);
}


void Solver::setPowerByValue(double value, int cell, int state){
  _power[state]->setValue(cell, 0, value);
}


void Solver::setWeightByValue(double value, int cell, int group, int state){
  _weight[state]->setValue(cell, group, value);
}


void Solver::integratePrecursorConcentrations(int state_from, int state_to){

  double dt = _clock->getTime(state_to) - _clock->getTime(state_from);

  #pragma omp parallel for
  for (int i=0; i < _num_cells; i++){

    double fission_rate_from = 0.0;
    double fission_rate_to   = 0.0;
    double temp_from         = getTemperatureByValue(i, state_from);
    double temp_to           = getTemperatureByValue(i, state_to);
    Material* material       = _geometry->getMaterial(i);
    
    if (material->isFissionable()){
      for (int g=0; g < _num_energy_groups; g++){
        fission_rate_from += material->getNuSigmaFByGroup
          (g, state_from,temp_from) * getFluxByValue(i, g, state_from);
        
        fission_rate_to   += material->getNuSigmaFByGroup
          (g, state_to, temp_to) * getFluxByValue(i, g, state_to);
      }
      
      for (int d=0; d < _num_delayed_groups; d++){
        double delayed_fraction = material->getDelayedFractionByGroup
          (d, state_from);
        double decay_constant   = material->getDecayConstantByGroup(d);
        double k1 = exp(- decay_constant * dt);
        double k2 = 1.0 - (1.0 - k1) / (decay_constant * dt);
        double k3 = k1 + (k2 - 1.0);
        double conc = material->getPrecursorConcByGroup(d, state_from);
        double new_conc = k1 * conc + k2 * delayed_fraction / decay_constant
          / _k_eff_0 * fission_rate_to - k3 * delayed_fraction / decay_constant
          / _k_eff_0 * fission_rate_from;

        material->setPrecursorConcByGroup(new_conc, d, state_to);
      }
    }
  }
}


void Solver::integrateTemperature(int state_from, int state_to){
  
  double dt = _clock->getTime(state_to) - _clock->getTime(state_from);

  #pragma omp parallel for
  for (int i=0; i < _num_cells; i++){

    double fission_rate_from = 0.0;
    double fission_rate_to   = 0.0;
    double temp_from         = getTemperatureByValue(i, state_from);
    double temp_to           = getTemperatureByValue(i, state_to);    
    Material* material       = _geometry->getMaterial(i);
    
    if (material->isFissionable()){
      for (int g=0; g < _num_energy_groups; g++){
        fission_rate_from += material->getSigmaFByGroup
          (g, state_from, temp_from) * getFluxByValue(i, g, state_from);
        fission_rate_to   += material->getSigmaFByGroup
          (g, state_to  , temp_to  ) * getFluxByValue(i, g, state_to);
      }
      
      temp_to = temp_from + dt * 0.5 * (fission_rate_from + fission_rate_to) *
        material->getTemperatureConversionFactor();
      setTemperatureByValue(temp_to, i, state_to);
    }
  }
}


void Solver::computePrecursorConcentrations(int state){

  #pragma omp parallel for
  for (int i=0; i < _num_cells; i++){

    double fission_rate = 0.0;
    Material* material  = _geometry->getMaterial(i);
    double temp         = getTemperatureByValue(i, state);
    
    if (material->isFissionable()){
      for (int g=0; g < _num_energy_groups; g++){
        fission_rate += material->getNuSigmaFByGroup(g, state, temp) * 
          getFluxByValue(i, g, state);
      }

      for (int d=0; d < _num_delayed_groups; d++){
        double delayed_fraction = material->getDelayedFractionByGroup(d, state);
        double decay_constant = material->getDecayConstantByGroup(d);
        double conc = delayed_fraction / decay_constant / _k_eff_0 * fission_rate;
        material->setPrecursorConcByGroup(conc, d, state);
      }
    }
  }
}


void Solver::computePower(int state){

  #pragma omp parallel for
  for (int i=0; i < _num_cells; i++){

    double cell_power  = 0.0;
    Material* material = _geometry->getMaterial(i);
    double temp        = getTemperatureByValue(i, state);
    
    if (material->isFissionable()){
      for (int g=0; g < _num_energy_groups; g++){
        cell_power += material->getSigmaFByGroup(g, state, temp) * 
          getFluxByValue(i, g, state) * material->getEnergyPerFission()
          * _geometry->getVolume(i);
      }
    }

    setPowerByValue(cell_power, i, state);
  }
}


double Solver::computeAveragePower(int state){
  
  computePower(state);
  double average_power = 0.0;
  double fuel_volume = 0.0;

  for (int i=0; i < _num_cells; i++){

    Material* material = _geometry->getMaterial(i);
    average_power += getPowerByValue(i, state);

    if (material->isFissionable())
      fuel_volume += _geometry->getVolume(i);
  }

  /* Divide cumulative power by fuel volume */
  average_power /= fuel_volume;

  return average_power;
}


double Solver::computePowerRMSError(int state_1, int state_2){
  
  computePower(state_1);
  computePower(state_2);
  Vector* power_1 = getPower(state_1);
  Vector* power_2 = getPower(state_2);
  double power_rmse = computeRMSE(power_1, power_2, false);
  return power_rmse;
}


void Solver::normalizeFlux(double value, int state){

  /* Compute the average power */
  double average_power = computeAveragePower(state);

  for (int i=0; i < _num_cells; i++){
    for (int e=0; e < _num_energy_groups; e++){
      double old_flux = getFluxByValue(i, e, state);
      setFluxByValue(old_flux * value / average_power, i, e, state);
    }
  }
}


void Solver::initializeClock(){
  
  for (int i=0; i < _num_cells; i++)
    _geometry->getMaterial(i)->setClock(_clock);
}


void Solver::copyPrecursors(int state_from, int state_to){

  for (int i=0; i < _num_cells; i++){
    Material* material = _geometry->getMaterial(i);
    if (material->isFissionable()){
      for (int d=0; d < _num_delayed_groups; d++){
        double conc = material->getPrecursorConcByGroup(d, state_from);
        material->setPrecursorConcByGroup(conc, d, state_to);
      }
    }
  }
}


void Solver::copyVariables(int state_from, int state_to){

  copyPrecursors(state_from, state_to);

  _temperature[state_from]->copyTo(_temperature[state_to]);
  _flux       [state_from]->copyTo(_flux       [state_to]);
  _power      [state_from]->copyTo(_power      [state_to]);
  _weight     [state_from]->copyTo(_weight     [state_to]);
}


void Solver::copyVariablesToAll(int state_from){

  for (int t=0; t < NUM_STATES; t++)
    copyVariables(state_from, t);
}
