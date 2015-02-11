#include "Mesh.h"


/**
 * @brief Constructor sets the ID and unique ID for the Material.
 * @param id the user-defined ID for the material
 */
Mesh::Mesh(double width, double height) {

  _materials = NULL;
  _offset = NULL;
  _clock = NULL;
  _decay_constants = NULL;
  _delayed_fractions = NULL;
  _k_eff_0 = 0.0;
  
  /* Set mesh properties */
  setXMin(-width/2.0);
  setXMax(width/2.0);
  setYMin(-height/2.0);
  setYMax(height/2.0);

  for (int c=0; c < 8; c++){
    _flux[c] = NULL;
    _power[c] = NULL;
    _temperature[c] = NULL;
  }

  _boundaries = new boundaryType[4];
  for (int i=0; i < 4; i++)
    _boundaries[i] = REFLECTIVE;
  
  return;
}


/**
 * @brief Destructor deletes all cross-section data structures from memory.
 */
Mesh::~Mesh() {

  for (int c=0; c < 8; c++){
    if (_flux[c] != NULL)
      delete [] _flux[c];

    if (_power[c] != NULL)
      delete [] _power[c];

    if (_temperature[c] != NULL)
      delete [] _temperature[c];

  }

  _flux.clear();
  _power.clear();
  _temperature.clear();

  if (_materials != NULL)
    delete [] _materials;
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


double Mesh::getWidth(){
  return _x_max - _x_min;
}


double Mesh::getHeight(){
  return _y_max - _y_min;
}


void Mesh::setBoundary(int side, boundaryType boundary){

  if (side < 0 || side > 3)
    log_printf(ERROR, "Unable to set boundary for side %d as there are only"
               " 4 geometry sides", side);

  _boundaries[side] = boundary;
}


boundaryType Mesh::getBoundary(int side){

  if (side < 0 || side > 3)
    log_printf(ERROR, "Unable to get boundary for side %d as there are only"
               " 4 geometry sides", side);
  
  return _boundaries[side];
}


void Mesh::setNumShapeEnergyGroups(int num_groups){
  _num_shape_energy_groups = num_groups;
}


void Mesh::setNumAmpEnergyGroups(int num_groups){
  _num_amp_energy_groups = num_groups;
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


int Mesh::getNumShapeEnergyGroups(){
  return _num_shape_energy_groups;
}


int Mesh::getNumAmpEnergyGroups(){
  return _num_amp_energy_groups;
}


int Mesh::getNumDelayedGroups(){
  return _num_delayed_groups;
}


double* Mesh::getFlux(int position){
  return _flux[position];
}


double* Mesh::getPower(int position){
  return _power[position];
}


double* Mesh::getTemperature(int position){
  return _temperature[position];
}


double Mesh::getPowerByValue(int cell, int position){
  return _power[position][cell];
}


double Mesh::getTemperatureByValue(int cell, int position){
  return _temperature[position][cell];
}


void Mesh::setMaterial(Material* material, int cell){
  _materials[cell] = material;
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


Material* Mesh::getMaterial(int cell){
  return _materials[cell];
}