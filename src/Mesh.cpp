#include "Mesh.h"


/**
 * @brief Constructor sets the ID and unique ID for the Material.
 * @param id the user-defined ID for the material
 */
Mesh::Mesh(double width, double height, double depth) {

  _materials = NULL;
  _clock = NULL;
  _decay_constants = NULL;
  _delayed_fractions = NULL;
  _k_eff_0 = 0.0;
  _fuel_volume = 0.0;
  _buckling = 0.0;
  _num_shape_energy_groups = 0;
  _num_amp_energy_groups = 0;
  _num_delayed_groups = 0;
  
  /* Set mesh properties */
  setXMin(-width/2.0);
  setXMax(width/2.0);
  setYMin(-height/2.0);
  setYMax(height/2.0);
  setZMin(-depth/2.0);
  setZMax(depth/2.0);

  for (int c=0; c < 8; c++){
    _flux[c] = NULL;
    _power[c] = NULL;
    _temperature[c] = NULL;
  }

  _boundaries = new boundaryType[6];
  for (int i=0; i < 6; i++)
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

  if (_k_eff_0 <= 0.0)
    log_printf(ERROR, "Unable to set the initial k-eff to a non-positive "
               "number: %f", _k_eff_0);

  _k_eff_0 = k_eff_0;
}


double Mesh::getKeff0(){

  if (_k_eff_0 == 0.0)
    log_printf(ERROR, "Unable to get the initial k-eff as it has not "
               "been set");

  return _k_eff_0;
}


Clock* Mesh::getClock(){
  
  if (_clock == NULL)
    log_printf(ERROR, "Unable to get the Mesh Clock since it has not been set");

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


double Mesh::getWidth(){
  return _x_max - _x_min;
}


double Mesh::getHeight(){
  return _y_max - _y_min;
}


double Mesh::getDepth(){
  return _z_max - _z_min;
}


void Mesh::setBoundary(int side, boundaryType boundary){

  if (side < 0 || side > 5)
    log_printf(ERROR, "Unable to set boundary for side %d as there are only"
               " 6 geometry sides", side);

  _boundaries[side] = boundary;
}


boundaryType Mesh::getBoundary(int side){

  if (side < 0 || side > 5)
    log_printf(ERROR, "Unable to get boundary for side %d as there are only"
               " 6 geometry sides", side);
  
  return _boundaries[side];
}


void Mesh::setNumShapeEnergyGroups(int num_groups){
  
  if (num_groups <= 0)
    log_printf(ERROR, "Unable to set the number of shape energy groups to a "
               "non-positive number: %i", num_groups);

  _num_shape_energy_groups = num_groups;
}


void Mesh::setNumAmpEnergyGroups(int num_groups){

  if (num_groups <= 0)
    log_printf(ERROR, "Unable to set the number of amp energy groups to a "
               "non-positive number: %i", num_groups);

  _num_amp_energy_groups = num_groups;
}


void Mesh::setNumDelayedGroups(int num_groups){

  if (num_groups <= 0)
    log_printf(ERROR, "Unable to set the number of delayed groups to a "
               "non-positive number: %i", num_groups);

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

  if (_num_shape_energy_groups == 0)
    log_printf(ERROR, "Unable to get the number of shape energy groups "
               "since it has not been set");

  return _num_shape_energy_groups;
}


int Mesh::getNumAmpEnergyGroups(){

  if (_num_amp_energy_groups == 0)
    log_printf(ERROR, "Unable to get the number of amp energy groups "
               "since it has not been set");

  return _num_amp_energy_groups;
}


int Mesh::getNumDelayedGroups(){

  if (_num_delayed_groups == 0)
    log_printf(ERROR, "Unable to get the number of delayed groups "
               "since it has not been set");

  return _num_delayed_groups;
}


double* Mesh::getFlux(int position){

  if (position < 0 || position >= 8)
    log_printf(ERROR, "Unable to get the Mesh flux for position %i since there"
               " only 8 groups", position);

  if (_flux[c] == NULL)
    log_printf(ERROR, "Unable to get the Mesh flux for position %i since the "
               "flux has not been initialized", position);

  return _flux[position];
}


double* Mesh::getPower(int position){

  if (position < 0 || position >= 8)
    log_printf(ERROR, "Unable to get the Mesh power for position %i since there"
               " only 8 groups", position);

  if (_power[c] == NULL)
    log_printf(ERROR, "Unable to get the Mesh power for position %i since the "
               "power has not been initialized", position);

  return _power[position];
}


double* Mesh::getTemperature(int position){

  if (position < 0 || position >= 8)
    log_printf(ERROR, "Unable to get the Mesh temperature for position %i since there"
               " only 8 groups", position);

  if (_temperature[c] == NULL)
    log_printf(ERROR, "Unable to get the Mesh temperature for position %i since the "
               "temperature has not been initialized", position);

  return _temperature[position];
}


double Mesh::getPowerByValue(int cell, int position){

  if (cell < 0 || cell >= getNumCells())
    log_printf(ERROR, "Unable to get the Mesh power for cell %i since there"
               " only %i cells", cell, getNumCells());

  if (position < 0 || position >= 8)
    log_printf(ERROR, "Unable to get the Mesh power for position %i since there"
               " only 8 groups", position);

  if (_power[c] == NULL)
    log_printf(ERROR, "Unable to get the Mesh power for position %i since the "
               "power has not been initialized", position);

  return _power[position][cell];
}


double Mesh::getTemperatureByValue(int cell, int position){

  if (cell < 0 || cell >= getNumCells())
    log_printf(ERROR, "Unable to get the Mesh temperature for cell %i since there"
               " only %i cells", cell, getNumCells());

  if (position < 0 || position >= 8)
    log_printf(ERROR, "Unable to get the Mesh temperature for position %i since there"
               " only 8 groups", position);

  if (_temperature[c] == NULL)
    log_printf(ERROR, "Unable to get the Mesh temperature for position %i since the "
               "temperature has not been initialized", position);

  return _temperature[position][cell];
}


void Mesh::setMaterial(Material* material, int cell){

  if (cell < 0 || cell >= getNumCells())
    log_printf(ERROR, "Unable to set the Mesh Material for cell %i since there"
               " only %i cells", cell, getNumCells());
  
  _materials[cell] = material;
}


void Mesh::setBuckling(double buckling){
  _buckling = buckling;
}


double Mesh::getBuckling(){
  return _buckling;
}


void Mesh::setDecayConstants(double* xs, int num_groups){

  if (num_groups != _num_delayed_groups)
    log_printf(ERROR, "Unable to set the decay constants for %i groups since "
               "there are %i delayed groups", num_groups, _num_delayed_groups);

  for (int i=0; i < num_groups; i++){

    if (xs[i] <= 0.0)
      log_printf(ERROR, "Unable to set the decay constant for %i group to a "
                 "non-positive value: %f", i, xs[i]);

    _decay_constants[i] = xs[i];
  }
}


void Mesh::setDelayedFractions(double* xs, int num_groups){

  if (num_groups != _num_delayed_groups)
    log_printf(ERROR, "Unable to set the delayed fraction for %i groups since "
               "there are %i delayed groups", num_groups, _num_delayed_groups);

  for (int i=0; i < num_groups; i++){

    if (xs[i] <= 0.0)
      log_printf(ERROR, "Unable to set the delayed fraction for %i group to a "
                 "non-positive value: %f", i, xs[i]);

    _delayed_fractions[i] = xs[i];  
  }
}


double* Mesh::getDecayConstants(){

  if (_decay_constants == NULL)
    log_printf(ERROR, "Unable to get the decay constants since they have "
               "not been initialilzed");

  return _decay_constants;
}


double* Mesh::getDelayedFractions(){

  if (_delayed_fractions == NULL)
    log_printf(ERROR, "Unable to get the delayed fractions since they have "
               "not been initialilzed");

  return _delayed_fractions;
}


double Mesh::getDecayConstantByGroup(int group){

  if (group < 0 || group >= _num_delayed_groups)
    log_printf(ERROR, "Unable to get the decay constant for group %i since "
               "there are only %i delayed groups", group, _num_delayed_groups);

  if (_decay_constants == NULL)
    log_printf(ERROR, "Unable to get the decay constant for group %i since "
               "they have not been initialilzed", group);

  return _decay_constants[group];
}


double Mesh::getDelayedFractionByGroup(int group){

  if (group < 0 || group >= _num_delayed_groups)
    log_printf(ERROR, "Unable to get the delayed fraction for group %i since "
               "there are only %i delayed groups", group, _num_delayed_groups);

  if (_delayed_fractions == NULL)
    log_printf(ERROR, "Unable to get the delayed fraction for group %i since "
               "they have not been initialilzed", group);

  return _delayed_fractions[group];
}


Material* Mesh::getMaterial(int cell){

  if (cell < 0 || cell >= getNumCells())
    log_printf(ERROR, "Unable to get the Mesh Material for cell %i since there"
               " only %i cells", cell, getNumCells());

  return _materials[cell];
}


void Mesh::computeFuelVolume(){

  if (_materials == NULL)
    log_printf(ERROR, "Unable to compute the fuel volume since the materials"
               " have not been initialized");

  _fuel_volume = 0.0;
  for (int i=0; i < getNumCells(); i++){
    if (_materials[i]->isFissionable())
      _fuel_volume += getCellVolume(i);
  }
}


void Mesh::uniquifyMaterials(){

  if (_materials == NULL)
    log_printf(ERROR, "Unable to uniquify the materials since the materials"
               " have not been initialized");

  Material* material;
  
  for (int i=0; i < getNumCells(); i++){
    material = _materials[i]->clone();
    _materials[i] = material;
  }
}


void Mesh::copyPower(int position_from, int position_to){

  if (position_from >= 8 || position_from < 0)
    log_printf(ERROR, "Unable to copy power from position %i since "
               "there are only 8 clock positions", position_from);

  if (position_to >= 8 || position_to < 0)
    log_printf(ERROR, "Unable to copy power to position %i since "
               "there are only 8 clock positions", position_to);
    
  if (_power[position_from] == NULL)
    log_printf(ERROR, "Unable to copy power from position %i since "
               "the powers have not been initialized yet", position_from);

  if (_power[position_to] == NULL)
    log_printf(ERROR, "Unable to copy power to position %i since "
               "the powers have not been initialized yet", position_to);

  std::copy(_power[position_from], 
            _power[position_from] + getNumCells(), 
            _power[position_to]);
}


void Mesh::copyTemperature(int position_from, int position_to){

  if (position_from >= 8 || position_from < 0)
    log_printf(ERROR, "Unable to copy temperature from position %i since "
               "there are only 8 clock positions", position_from);

  if (position_to >= 8 || position_to < 0)
    log_printf(ERROR, "Unable to copy temperature to position %i since "
               "there are only 8 clock positions", position_to);
    
  if (_temperature[position_from] == NULL)
    log_printf(ERROR, "Unable to copy temperature from position %i since "
               "the temperatures have not been initialized yet", position_from);

  if (_temperature[position_to] == NULL)
    log_printf(ERROR, "Unable to copy temperature to position %i since "
               "the temperatures have not been initialized yet", position_to);

  std::copy(_temperature[position_from], 
            _temperature[position_from] + getNumCells(), 
            _temperature[position_to]);
}


meshType Mesh::getMeshType(){
  return _mesh_type;
}


void Mesh::initialize(){
  
  if (_num_shape_energy_groups == 0)
    log_printf(ERROR, "Unable to initialize the Mesh since the number "
               "of shape energy groups has not been set");

  if (_num_amp_energy_groups == 0)
    log_printf(ERROR, "Unable to initialize the Mesh since the number "
               "of amp energy groups has not been set");

  if (_num_delayed_groups == 0)
    log_printf(ERROR, "Unable to initialize the Mesh since the number "
               "of delayed groups has not been set");

  for (int c=0; c < 8; c++){
    if (_power[c] != NULL)
      delete [] _power[c];

    if (_temperature[c] != NULL)
      delete [] _temperature[c];
  }

  for(int c=0; c < 8; c++){
    _temperature[c] = new double[getNumCells()];
    _power[c] = new double[getNumCells()];

    memset(_temperature[c], 300.0, sizeof(double) * getNumCells());
    memset(_power[c], 0.0, sizeof(double) * getNumCells());
  }  

  if (_materials != NULL)
    delete [] _materials;

  _materials = new Material*[getNumCells()];

}


void Mesh::setClock(Clock* clock){
  _clock = clock;

  if (_materials != NULL){
    for (int i=0; i < getNumCells(); i++)
      _materials[i]->getClock(clock);
  }
}
