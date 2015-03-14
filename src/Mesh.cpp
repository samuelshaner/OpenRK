#include "Mesh.h"


/**
 * @brief Constructor sets the ID and unique ID for the Material.
 * @param id the user-defined ID for the material
 */
Mesh::Mesh(double width_x, double width_y, double width_z) {

  _materials = NULL;
  _clock = NULL;
  _decay_constants = NULL;
  _delayed_fractions = NULL;
  _k_eff_0 = 0.0;
  
  /* Set mesh properties */
  setXMin(-width_x/2.0);
  setXMax(width_x/2.0);
  setYMin(-width_y/2.0);
  setYMax(width_y/2.0);
  setZMin(-width_z/2.0);
  setZMax(width_z/2.0);
  
  _boundaries = new boundaryType[6];
  for (int i=0; i < 6; i++)
    _boundaries[i] = REFLECTIVE;
  
  return;
}


/**
 * @brief Destructor deletes all cross-section data structures from memory.
 */
Mesh::~Mesh() {

  std::map<int, std::map<int, double*> >::iterator iter;
  
  for (iter=_field_variables.begin(); iter != _field_variables.end(); ++iter){
    for (int c=0; c < NUM_CLOCK_POSITIONS; c++){
      if (iter->second[c] != NULL)
        delete [] iter->second[c];
    }

    iter->second.clear();
  }
  
  _field_variables.clear();
}


void Mesh::setKeff0(double k_eff_0){
  _k_eff_0 = k_eff_0;
}


double Mesh::getKeff0(){
  return _k_eff_0;
}


Clock* Mesh::getClock(){
  return _clock;
}


void Mesh::setXMin(double x_min){
  _x_min = x_min;
}


void Mesh::setXMax(double x_max){
  _x_max = x_max;
}


void Mesh::setYMin(double y_min){
  _y_min = y_min;
}


void Mesh::setYMax(double y_max){
  _y_max = y_max;
}

void Mesh::setZMin(double z_min){
  _z_min = z_min;
}


void Mesh::setZMax(double z_max){
  _z_max = z_max;
}


double Mesh::getXMin(){
  return _x_min;
}


double Mesh::getXMax(){
  return _x_max;
}


double Mesh::getYMin(){
  return _y_min;
}


double Mesh::getYMax(){
  return _y_max;
}


double Mesh::getZMin(){
  return _z_min;
}


double Mesh::getZMax(){
  return _z_max;
}


double Mesh::getWidthX(){
  return _x_max - _x_min;
}


double Mesh::getWidthY(){
  return _y_max - _y_min;
}


double Mesh::getWidthZ(){
  return _z_max - _z_min;
}


void Mesh::setBoundary(surfaceLocation side, boundaryType boundary){

  if (side < 0 || side > 5)
    log_printf(ERROR, "Unable to set boundary for side %d as there are only"
               " 6 geometry sides", side);

  _boundaries[side] = boundary;
}


boundaryType Mesh::getBoundary(surfaceLocation side){

  if (side < 0 || side > 5)
    log_printf(ERROR, "Unable to get boundary for side %d as there are only"
               " 6 geometry sides", side);
  
  return _boundaries[side];
}


void Mesh::setNumEnergyGroups(int num_groups){
  _num_energy_groups = num_groups;
}


void Mesh::setNumDelayedGroups(int num_groups){
  _num_delayed_groups = num_groups;

  if (_decay_constants != NULL)
    delete [] _decay_constants;

  if (_delayed_fractions != NULL)
    delete [] _delayed_fractions;

  _decay_constants = new double[_num_delayed_groups];
  _delayed_fractions = new double[_num_delayed_groups];  

  memset(_decay_constants, 0.0, sizeof(double) * _num_delayed_groups);
  memset(_delayed_fractions, 0.0, sizeof(double) * _num_delayed_groups);

}


int Mesh::getNumEnergyGroups(){
  return _num_energy_groups;
}


int Mesh::getNumDelayedGroups(){
  return _num_delayed_groups;
}


double* Mesh::getFieldVariable(fieldVariable name, int position){
  return _field_variables[name][position];
}


double Mesh::getFieldVariableByValue(fieldVariable name, int position, int cell, int group){
  int num_groups = getNumFieldVariableGroups(name);
  return _power[name][position][cell*num_groups + group];
}

double Mesh::copyFieldVariable(fieldVariable name, int position_from, int position_to){
  int num_groups = getNumFieldVariableGroups(name);
  std::copy(_field_variables[name][position_from], _field_variables[name][position_from]\
            + getNumCells() * num_groups, _field_variables[name][position_to]);
}


void Mesh::setFieldVariableByValue(fieldVariable name, double value, int cell, \
                                   int position, int group){
  int num_groups = getNumFieldVariableGroups(name);
  _field_variables[name][position][cell*num_groups + group] = value;
}


void Mesh::setBuckling(double buckling){
  _buckling = buckling;
}


double Mesh::getBuckling(){
  return _buckling;
}


void Mesh::setDecayConstants(double* xs, int num_groups){

  for (int i=0; i < num_groups; i++)
    _decay_constants[i] = xs[i];
}


void Mesh::setDelayedFractions(double* xs, int num_groups){

  for (int i=0; i < num_groups; i++)
    _delayed_fractions[i] = xs[i];  
}


double* Mesh::getDecayConstants(){
  return _decay_constants;
}


double* Mesh::getDelayedFractions(){
  return _delayed_fractions;
}


double Mesh::getDecayConstantByGroup(int group){
  return _decay_constants[group];
}


double Mesh::getDelayedFractionByGroup(int group){
  return _delayed_fractions[group];
}


void Mesh::setClock(Clock* clock){
  _clock = clock;
}
