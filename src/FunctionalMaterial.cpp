#include "FunctionalMaterial.h"

/**
 * @brief Constructor sets the ID and unique ID for the Material.
 * @param id the user-defined ID for the material
 */
FunctionalMaterial::FunctionalMaterial(int id) : Material(id) {

  _doppler_coefficients = NULL;
  _time_steps = NULL;
  _num_time_steps = 0;
  
  return;
}


/**
 * @brief Destructor deletes all cross-section data structures from memory.
 */
FunctionalMaterial::~FunctionalMaterial() {

  if (_doppler_coefficients != NULL)
    delete [] _doppler_coefficients;
  
  if (_time_steps != NULL)
    delete [] _time_steps;  
}


/**
 * @brief Set the number of energy groups for this Material.
 * @param num_groups the number of energy groups.
 */
void FunctionalMaterial::setNumEnergyGroups(const int num_groups) {

  if (num_groups < 0)
    log_printf(ERROR, "Unable to set the number of energy groups for "
               "material %d to %d", _id, num_groups);

  if (_num_time_steps == 0)
    log_printf(ERROR, "Unable to set the number of energy groups for "
               "material %d since the number of time steps has not been set", _id);
      
  _num_energy_groups = num_groups;

  /* Free old data arrays if they were allocated for a previous simulation */
  if (_sigma_t != NULL)
    delete [] _sigma_t;
  
  if (_sigma_a != NULL)
    delete [] _sigma_a;
  
  if (_sigma_s != NULL)
    delete [] _sigma_s;
  
  if (_sigma_f != NULL)
    delete [] _sigma_f;
  
  if (_nu_sigma_f != NULL)
    delete [] _nu_sigma_f;
  
  if (_chi != NULL)
    delete [] _chi;
  
  if (_dif_coef != NULL)
    delete [] _dif_coef;
  
  if (_velocity != NULL)
    delete [] _velocity;

  if (_doppler_coefficients != NULL)
    delete [] _doppler_coefficients;

  /* Allocate memory for data arrays */
  _sigma_t = new double[_num_time_steps*_num_energy_groups];
  _sigma_a = new double[_num_time_steps*_num_energy_groups];
  _sigma_f = new double[_num_time_steps*_num_energy_groups];
  _nu_sigma_f = new double[_num_time_steps*_num_energy_groups];
  _chi = new double[_num_time_steps*_num_energy_groups];
  _sigma_s = new double[_num_time_steps*_num_energy_groups*_num_energy_groups];
  _velocity = new double[_num_time_steps*_num_energy_groups];
  _dif_coef = new double[_num_time_steps*_num_energy_groups];
  _doppler_coefficients = new double[_num_energy_groups];
  
  /* Assign the null vector to each data array */
  memset(_sigma_t, 0.0, _num_time_steps * sizeof(double) * _num_energy_groups);
  memset(_sigma_a, 0.0, _num_time_steps * sizeof(double) * _num_energy_groups);
  memset(_sigma_f, 0.0, _num_time_steps * sizeof(double) * _num_energy_groups);
  memset(_nu_sigma_f, 0.0, _num_time_steps * sizeof(double) * _num_energy_groups);
  memset(_chi, 0.0, _num_time_steps * sizeof(double) * _num_energy_groups);
  memset(_sigma_s, 0.0, _num_time_steps * sizeof(double) * _num_energy_groups * _num_energy_groups);
  memset(_velocity, 0.0, _num_time_steps * sizeof(double) * _num_energy_groups);
  memset(_dif_coef, 0.0, _num_time_steps * sizeof(double) * _num_energy_groups);
  memset(_doppler_coefficients, 0.0, sizeof(double) * _num_energy_groups);
}


void FunctionalMaterial::setTimeSteps(double* time_steps, int num_steps) {

  if (num_steps < 0)
    log_printf(ERROR, "Unable to set the number of time steps for "
               "material %d to %d", _id, num_steps);

  if (_time_steps != NULL)
    delete [] _time_steps;  
  
  _num_time_steps = num_steps;
  _time_steps = new double[num_steps];
  
  for (int i=0; i < _num_time_steps; i++)
    _time_steps[i] = time_steps[i];
}



/**
 * @brief Set the Material's array of total cross-sections.
 * @details This method is a helper function to allow OpenMOC users to assign
 *          the Material's nuclear data in Python. A user must initialize a
 *          NumPy array of the correct size (i.e., a float64 array the length
 *          of the number of energy groups) as input to this function. This
 *          function then fills the NumPy array with the data values for the
 *          Material's total cross-sections. An example of how this function
 *          might be called in Python is as follows:
 *
 * @code
 *          sigma_t = numpy.array([0.05, 0.1, 0.15, ... ])
 *          material = openmoc.Material(openmoc.material_id())
 *          material.setSigmaT(sigma_t)
 * @endcode
 *
 * @param xs the array of total cross-sections
 * @param num_groups the number of energy groups
 */
void FunctionalMaterial::setSigmaT(double* xs, int num_time_steps, int num_groups) {

  if (_num_time_steps == 0)
    log_printf(ERROR, "Unable to set the sigma t for material %d since"
               " the number of time steps has not been set", _id);

  if (_num_energy_groups != num_groups)
    log_printf(ERROR, "Unable to set sigma_t with %d groups for Material "
               "%d which contains %d energy groups",
               num_groups, _id, _num_energy_groups);

  if (_num_time_steps != num_time_steps)
    log_printf(ERROR, "Unable to set sigma_t with %d time steps for Material "
               "%d which contains %d time steps",
               num_time_steps, _id, _num_time_steps);

  for (int c=0; c < _num_time_steps; c++){
    for (int i=0; i < _num_energy_groups; i++)
      _sigma_t[c*_num_energy_groups+i] = double(xs[c*_num_energy_groups+i]);
  }
}


/**
 * @brief Set the Material's total cross-section for some energy group.
 * @param xs the total cross-section
 * @param group the energy group
 */
void FunctionalMaterial::setSigmaTByGroup(double xs, int group, int position) {

  log_printf(ERROR, "Unable to set the sigma t by group for FunctionalMaterial %d"
             " since this Material subclass doesn't support settig xs by group", _id);
}


/**
 * @brief Set the Material's total cross-section for some energy group.
 * @param xs the total cross-section
 * @param group the energy group
 */
double FunctionalMaterial::getSigmaTByGroup(int group, int position, double temp) {

  if (_num_time_steps == 0)
    log_printf(ERROR, "Unable to get the sigma t for material %d and "
               "group %d since the number of time steps has not been set",
               _id, group);
  
  if (group < 0 || group >= _num_energy_groups)
    log_printf(ERROR, "Unable to get sigma_t for group %d for Material "
               "%d which contains %d energy groups", group, _id, _num_energy_groups);

  int time_step = getTimeStep(position);
  int ng = _num_energy_groups;
  
  double sigma_t_left = _sigma_t[time_step*ng + group];
  double sigma_t_right = _sigma_t[(time_step+1)*ng + group];
  double time = _clock->getTime(position);
  double time_left = _time_steps[time_step];
  double time_right = _time_steps[time_step+1];

  /* Interpolate sigma t */
  double sigma_t = sigma_t_left + (time - time_left) * (sigma_t_right - sigma_t_left)
    / (time_right - time_left);
  
  return sigma_t;
}


/**
 * @brief Set the Material's array of absorption scattering cross-sections.
 * @details This method is a helper function to allow OpenMOC users to assign
 *          the Material's nuclear data in Python. A user must initialize a
 *          NumPy array of the correct size (i.e., a float64 array the length
 *          of the number of energy groups) as input to this function. This
 *          function then fills the NumPy array with the data values for the
 *          Material's absorption cross-sections. An example of how this
 *          function might be called in Python is as follows:
 *
 * @code
 *          sigma_a = numpy.array([0.05, 0.1, 0.15, ... ])
 *          material = openmoc.Material(openmoc.material_id())
 *          material.setSigmaA(sigma_a)
 * @endcode
 *
 * @param xs the array of absorption scattering cross-sections
 * @param num_groups the number of energy groups
 */
void FunctionalMaterial::setSigmaA(double* xs, int num_time_steps, int num_groups) {

  if (_num_time_steps == 0)
    log_printf(ERROR, "Unable to set the sigma_a for material %d since"
               " the number of time steps has not been set", _id);
  
  if (_num_energy_groups != num_groups)
    log_printf(ERROR, "Unable to set sigma_a with %d groups for Material "
               "%d which contains %d energy groups",
               num_groups, _id, _num_energy_groups);

  if (_num_time_steps != num_time_steps)
    log_printf(ERROR, "Unable to set sigma_a with %d time steps for Material "
               "%d which contains %d time steps",
               num_time_steps, _id, _num_time_steps);

  for (int c=0; c < _num_time_steps; c++){
    for (int i=0; i < _num_energy_groups; i++)
      _sigma_a[c*_num_energy_groups+i] = double(xs[c*_num_energy_groups+i]);
  }
}


/**
 * @brief Set the Material's absorption cross-section for some energy group.
 * @param xs the absorption cross-section
 * @param group the energy group
 */
void FunctionalMaterial::setSigmaAByGroup(double xs, int group, int position) {

  log_printf(ERROR, "Unable to set the sigma t by group for FunctionalMaterial %d"
             " since this Material subclass doesn't support settig xs by group", _id);
}


/**
 * @brief Set the Material's absorption cross-section for some energy group.
 * @param xs the absorption cross-section
 * @param group the energy group
 */
double FunctionalMaterial::getSigmaAByGroup(int group, int position, double temp) {

  if (_num_time_steps == 0)
    log_printf(ERROR, "Unable to get the sigma a for material %d and "
               "group %d since the number of time steps has not been set",
               _id, group);

  if (group < 0 || group >= _num_energy_groups)
    log_printf(ERROR, "Unable to get sigma_a for group %d for material "
               "%d which contains %d energy groups", group, _id, _num_energy_groups);

  if (temp == 0.0)
    log_printf(ERROR, "Unable to get sigma_a by group for FunctionalMaterial %d"
               " with temperature %.6f", _id, temp);
  
  if (_doppler_coefficients == NULL)
    log_printf(ERROR, "Unable to get sigma_a by group for FunctionalMaterial %d"
               " since the doppler coefficients have not been set", _id);

  int time_step = getTimeStep(position);
  int ng = _num_energy_groups;

  double sigma_a_left = _sigma_a[time_step*ng + group] *
    (1.0 + _doppler_coefficients[group] *
     (sqrt(temp) - sqrt(300.0)));
  double sigma_a_right = _sigma_a[(time_step+1)*ng + group]  *
    (1.0 + _doppler_coefficients[group] *
     (sqrt(temp) - sqrt(300.0)));
  double time = _clock->getTime(position);
  double time_left = _time_steps[time_step];
  double time_right = _time_steps[time_step+1];
  
  /* Interpolate sigma t */
  double sigma_a = sigma_a_left + (time - time_left) * (sigma_a_right - sigma_a_left)
    / (time_right - time_left);

  //log_printf(NORMAL, "time: %f, time left: %f, time right: %f, sigma a: %f", time, time_left, time_right, sigma_a);
  
  return sigma_a;
}


/**
 * @brief Set the Material's 2D array of scattering cross-sections.
 * @details This assumes that the scattering matrix passed in has the standard
 *          notation: the (i,j) element is for scattering from group i to j. For
 *          efficient caching of the elements of this matrix during fixed source
 *          iteration, the matrix transpose is what is actually stored in the
 *          Material.
 *
 *          This method is a helper function to allow OpenMOC users to assign
 *          the Material's nuclear data in Python. A user must initialize a
 *          NumPy array of the correct size (i.e., a float64 array the length
 *          of the number of energy groups) as input to this function. This
 *          function then fills the NumPy array with the data values for the
 *          Material's scattering cross-sections. An example of how this
 *          function might be called in Python is as follows:
 *
 * @code
 *          sigma_s = numpy.array([[0.05, 0.1, 0.15, ... ],
 *                                 [...      ...      ...],
 *                                           ...
 *                                 [...      ...      ...]])
 *          material = openmoc.Material(openmoc.material_id())
 *          material.setSigmaS(sigma_s)
 * @endcode
 *
 * @param xs the array of scattering cross-sections
 * @param num_groups_squared the number of energy groups squared
 */
void FunctionalMaterial::setSigmaS(double* xs, int num_time_steps, int num_groups_squared) {
  
  if (_num_time_steps == 0)
    log_printf(ERROR, "Unable to set the sigma s for material %d since"
               " the number of time steps has not been set", _id);

  int ng = _num_energy_groups;
  int nt = _num_time_steps;
  
  if (ng*ng != num_groups_squared)
    log_printf(ERROR, "Unable to set sigma_s with %f groups for Material %d "
               "which contains %d energy groups",
               float(sqrt(num_groups_squared)), _id, _num_energy_groups, _num_time_steps);
  
  if (_num_time_steps != num_time_steps)
    log_printf(ERROR, "Unable to set sigma_s with %d time steps for Material "
               "%d which contains %d time steps",
               num_time_steps, _id, _num_time_steps);

  for (int c=0; c < _num_time_steps; c++){
    for (int i=0; i < ng; i++) {
      for (int j=0; j < ng; j++)
        _sigma_s[c*ng*ng+j*ng+i] = xs[c*ng*ng+j*ng+i];
    }
  }
}


/**
 * @brief Set the Material's scattering cross-section for some energy group.
 * @param xs the scattering cross-section
 * @param group1 the row index in the scattering matrix
 * @param group2 the column index in the scattering matrix
 */
void FunctionalMaterial::setSigmaSByGroup(double xs, int group_from, int group_to, int position) {

  log_printf(ERROR, "Unable to set the sigma t by group for FunctionalMaterial %d"
             " since this Material subclass doesn't support settig xs by group", _id);
}


/**
 * @brief Set the Material's scattering cross-section for some energy group.
 * @param xs the scattering cross-section
 * @param group1 the row index in the scattering matrix
 * @param group2 the column index in the scattering matrix
 */
double FunctionalMaterial::getSigmaSByGroup(int group_from, int group_to, int position, double temp) {

  if (_num_time_steps == 0)
    log_printf(ERROR, "Unable to get the sigma s for material %d and "
               "group %d to %d since the number of time steps has not been set",
               _id, group_from, group_to);
  
  if (group_from < 0 || group_to < 0 || group_from >= _num_energy_groups || group_to >= _num_energy_groups)
    log_printf(ERROR, "Unable to get sigma_s for group %d to %d for Material %d "
               "which contains %d energy groups",
               group_from, group_to, _id, _num_energy_groups);

  int time_step = getTimeStep(position);
  int ng = _num_energy_groups;
  
  double sigma_s_left = _sigma_s[time_step*ng*ng + group_from*ng + group_to];
  double sigma_s_right = _sigma_s[(time_step+1)*ng*ng + group_from*ng + group_to];
  double time = _clock->getTime(position);
  double time_left = _time_steps[time_step];
  double time_right = _time_steps[time_step+1];

  /* Interpolate sigma t */
  double sigma_s = sigma_s_left + (time - time_left) * (sigma_s_right - sigma_s_left)
    / (time_right - time_left);
  
  return sigma_s;
}


/**
 * @brief Set the Material's array of fission cross-sections.
 * @details This method is a helper function to allow OpenMOC users to assign
 *          the Material's nuclear data in Python. A user must initialize a
 *          NumPy array of the correct size (i.e., a float64 array the length
 *          of the number of energy groups) as input to this function. This
 *          function then fills the NumPy array with the data values for the
 *          Material's fission cross-sections. An example of how this
 *          function might be called in Python is as follows:
 *
 * @code
 *          sigma_f = numpy.array([0.05, 0.1, 0.15, ... ])
 *          material = openmoc.Material(openmoc.material_id())
 *          material.setSigmaF(sigma_f)
 * @endcode
 *
 * @param xs the array of fission cross-sections
 * @param num_groups the number of energy groups
 */
void FunctionalMaterial::setSigmaF(double* xs, int num_time_steps, int num_groups) {

  if (_num_time_steps == 0)
    log_printf(ERROR, "Unable to set the sigma_f for material %d since "
               "the number of time steps has not been set", _id);

  if (_num_energy_groups != num_groups)
    log_printf(ERROR, "Unable to set sigma_f with %d groups for Material "
               "%d which contains %d energy groups",
               num_groups, _id, _num_energy_groups);

  if (_num_time_steps != num_time_steps)
    log_printf(ERROR, "Unable to set sigma_f with %d time steps for Material "
               "%d which contains %d time steps",
               num_time_steps, _id, _num_time_steps);

  for (int c=0; c < _num_time_steps; c++){
    for (int i=0; i < _num_energy_groups; i++)
      _sigma_f[c*_num_energy_groups+i] = xs[c*_num_energy_groups+i];
  }
  
  /* Determine whether or not this Material is fissionable */
  _fissionable = false;

  for (int i=0; i < _num_energy_groups*_num_time_steps; i++) {
    if (_sigma_f[i] > 0.0) {
      _fissionable = true;
      return;
    }
  }
}


/**
 * @brief Set the Material's fission cross-section for some energy group.
 * @param xs the fission cross-section
 * @param group the energy group
 */
void FunctionalMaterial::setSigmaFByGroup(double xs, int group, int position) {

  log_printf(ERROR, "Unable to set the sigma t by group for FunctionalMaterial %d"
             " since this Material subclass doesn't support settig xs by group", _id);
}


/**
 * @brief Set the Material's fission cross-section for some energy group.
 * @param xs the fission cross-section
 * @param group the energy group
 */
double FunctionalMaterial::getSigmaFByGroup(int group, int position, double temp) {
  
  if (_num_time_steps == 0)
    log_printf(ERROR, "Unable to get the sigma f for material %d and "
               "group %d since the number of time steps has not been set",
               _id, group);
  
  if (group < 0 || group >= _num_energy_groups)
    log_printf(ERROR, "Unable to get sigma_f for group %d for Material "
               "%d which contains %d energy groups", group, _id, _num_energy_groups);

  int time_step = getTimeStep(position);
  int ng = _num_energy_groups;

  double sigma_f_left = _sigma_f[time_step*ng + group];
  double sigma_f_right = _sigma_f[(time_step+1)*ng + group];
  double time = _clock->getTime(position);
  double time_left = _time_steps[time_step];
  double time_right = _time_steps[time_step+1];

  /* Interpolate sigma f */
  double sigma_f = sigma_f_left + (time - time_left) * (sigma_f_right - sigma_f_left)
    / (time_right - time_left);

  return sigma_f;
}


/**
 * @brief Set the Material's array of fission cross-sections multiplied by
 *         \f$ \nu \f$
 * @param xs the array of fission cross-sections multiplied by nu
 *        \f$ \nu \f$
 * @param num_groups the number of energy groups
*/
void FunctionalMaterial::setNuSigmaF(double* xs, int num_time_steps, int num_groups) {

  if (_num_time_steps == 0)
    log_printf(ERROR, "Unable to set the nu_sigma_f for material %d since "
               "the number of time steps has not been set", _id);

  if (_num_energy_groups != num_groups)
    log_printf(ERROR, "Unable to set nu_sigma_f with %d groups for Material "
               "%d which contains %d energy groups",
               num_groups, _id, _num_energy_groups);

  if (_num_time_steps != num_time_steps)
    log_printf(ERROR, "Unable to set nu_sigma_f with %d time steps for Material "
               "%d which contains %d time steps",
               num_time_steps, _id, _num_time_steps);

  for (int c=0; c < _num_time_steps; c++){
    for (int i=0; i < _num_energy_groups; i++)
      _nu_sigma_f[c*_num_energy_groups+i] = xs[c*_num_energy_groups+i];
  }
}


/**
 * @brief Set the Material's fission cross-section multiplied by \f$ \nu \f$
 *        for some energy group.
 * @details This method is a helper function to allow OpenMOC users to assign
 *          the Material's nuclear data in Python. A user must initialize a
 *          NumPy array of the correct size (i.e., a float64 array the length
 *          of the number of energy groups) as input to this function. This
 *          function then fills the NumPy array with the data values for the
 *          Material's nu*fission cross-sections. An example of how this
 *          function might be called in Python is as follows:
 *
 * @code
 *          nu_sigma_f = numpy.array([0.05, 0.1, 0.15, ... ])
 *          material = openmoc.Material(openmoc.material_id())
 *          material.setNuSigmaF(nu_sigma_f)
 * @endcode
 *
 * @param xs the fission cross-section multiplied by nu \f$ \nu \f$
 * @param group the energy group
 */
void FunctionalMaterial::setNuSigmaFByGroup(double xs, int group, int position) {

  log_printf(ERROR, "Unable to set the sigma t by group for FunctionalMaterial %d"
             " since this Material subclass doesn't support settig xs by group", _id);
}


/**
 * @brief Set the Material's fission cross-section multiplied by \f$ \nu \f$
 *        for some energy group.
 * @details This method is a helper function to allow OpenMOC users to assign
 *          the Material's nuclear data in Python. A user must initialize a
 *          NumPy array of the correct size (i.e., a float64 array the length
 *          of the number of energy groups) as input to this function. This
 *          function then fills the NumPy array with the data values for the
 *          Material's nu*fission cross-sections. An example of how this
 *          function might be called in Python is as follows:
 *
 * @code
 *          nu_sigma_f = numpy.array([0.05, 0.1, 0.15, ... ])
 *          material = openmoc.Material(openmoc.material_id())
 *          material.setNuSigmaF(nu_sigma_f)
 * @endcode
 *
 * @param xs the fission cross-section multiplied by nu \f$ \nu \f$
 * @param group the energy group
 */
double FunctionalMaterial::getNuSigmaFByGroup(int group, int position, double temp) {

  if (_num_time_steps == 0)
    log_printf(ERROR, "Unable to get the nu sigma f for material %d and "
               "group %d since the number of time steps has not been set",
               _id, group);
  
  if (group < 0 || group >= _num_energy_groups)
    log_printf(ERROR, "Unable to get nu_sigma_f for group %d for Material "
               "%d which contains %d energy groups", group, _id, _num_energy_groups);

  int time_step = getTimeStep(position);
  int ng = _num_energy_groups;

  double nu_sigma_f_left = _nu_sigma_f[time_step*ng + group];
  double nu_sigma_f_right = _nu_sigma_f[(time_step+1)*ng + group];
  double time = _clock->getTime(position);

  double time_left = _time_steps[time_step];
  double time_right = _time_steps[time_step+1];

  /* Interpolate sigma f */
  double nu_sigma_f = nu_sigma_f_left + (time - time_left) * (nu_sigma_f_right - nu_sigma_f_left)
    / (time_right - time_left);
  
  return nu_sigma_f;
}



/**
 * @brief Set the Material's array of chi \f$ \chi \f$ values.
 * @details This method is a helper function to allow OpenMOC users to assign
 *          the Material's nuclear data in Python. A user must initialize a
 *          NumPy array of the correct size (i.e., a float64 array the length
 *          of the number of energy groups) as input to this function. This
 *          function then fills the NumPy array with the data values for the
 *          Material's chi distribution. An example of how this function might
 *          be called in Python is as follows:
 *
 * @code
 *          chi = numpy.array([0.05, 0.1, 0.15, ... ])
 *          material = openmoc.Material(openmoc.material_id())
 *          material.setChi(chi)
 * @endcode
 *
 * @param xs the array of chi \f$ \chi \f$ values
 * @param num_groups the number of energy groups
 */
void FunctionalMaterial::setChi(double* xs, int num_time_steps, int num_groups) {

  if (_num_time_steps == 0)
    log_printf(ERROR, "Unable to set the chi for material %d since "
               "the number of time steps has not been set", _id);

  if (_num_energy_groups != num_groups)
    log_printf(ERROR, "Unable to set chi with %d groups for Material "
               "%d which contains %d energy groups",
               num_groups, _id, _num_energy_groups);

  if (_num_time_steps != num_time_steps)
    log_printf(ERROR, "Unable to set chi with %d time steps for Material "
               "%d which contains %d time steps",
               num_time_steps, _id, _num_time_steps);  
  
  for (int c=0; c < _num_time_steps; c++){
    for (int i=0; i < _num_energy_groups; i++)
      _chi[c*_num_energy_groups+i] = xs[c*_num_energy_groups+i];
  }
}


/**
 * @brief Set the Material's chi value for some energy group.
 * @param xs the chi value (\f$ \Chi \f$)
 * @param group the energy group
 */
void FunctionalMaterial::setChiByGroup(double xs, int group, int position) {

  log_printf(ERROR, "Unable to set the sigma t by group for FunctionalMaterial %d"
             " since this Material subclass doesn't support settig xs by group", _id);
}


/**
 * @brief Set the Material's chi value for some energy group.
 * @param xs the chi value (\f$ \Chi \f$)
 * @param group the energy group
 */
double FunctionalMaterial::getChiByGroup(int group, int position, double temp) {

  if (_num_time_steps == 0)
    log_printf(ERROR, "Unable to get the chi for material %d and "
               "group %d since the number of time steps has not been set",
               _id, group);
  
  if (group < 0 || group >= _num_energy_groups)
    log_printf(ERROR, "Unable to get chi for group %d for Material "
               "%d which contains %d energy groups", group, _id, _num_energy_groups);

  int time_step = getTimeStep(position);
  int ng = _num_energy_groups;
  
  double chi_left = _chi[time_step*ng + group];
  double chi_right = _chi[(time_step+1)*ng + group];
  double time = _clock->getTime(position);
  double time_left = _time_steps[time_step];
  double time_right = _time_steps[time_step+1];

  /* Interpolate sigma f */
  double chi = chi_left + (time - time_left) * (chi_right - chi_left)
    / (time_right - time_left);
  
  return chi;
}


/**
 * @brief Set the Material's array of diffusion coefficients.
 * @details This method is a helper function to allow OpenMOC users to assign
 *          the Material's nuclear data in Python. A user must initialize a
 *          NumPy array of the correct size (i.e., a float64 array the length
 *          of the number of energy groups) as input to this function. This
 *          function then fills the NumPy array with the data values for the
 *          Material's diffusion coefficients. An example of how this
 *          function might be called in Python is as follows:
 *
 * @code
 *          dif_coef = numpy.array([0.05, 0.1, 0.15, ... ])
 *          material = openmoc.Material(openmoc.material_id())
 *          material.setDifCoef(dif_coef)
 * @endcode
 *
 * @param xs the array of diffusion coefficents
 * @param num_groups the number of energy groups
 */
void FunctionalMaterial::setDifCoef(double* xs, int num_time_steps, int num_groups) {

  if (_num_time_steps == 0)
    log_printf(ERROR, "Unable to set the dif_coef for material %d since "
               "the number of time steps has not been set", _id);

  if (_num_energy_groups != num_groups)
    log_printf(ERROR, "Unable to set dif_coef with %d groups for Material "
               "%d which contains %d energy groups",
               num_groups, _id, _num_energy_groups);

  if (_num_time_steps != num_time_steps)
    log_printf(ERROR, "Unable to set dif_coef with %d time steps for Material "
               "%d which contains %d time steps",
               num_time_steps, _id, _num_time_steps);    

  for (int c=0; c < _num_time_steps; c++){
    for (int i=0; i < _num_energy_groups; i++)
      _dif_coef[c*_num_energy_groups+i] = xs[c*_num_energy_groups+i];
  }
}


/**
 * @brief Set the Material's diffusion coefficient for some energy group.
 * @param xs the diffusion coefficient
 * @param group the energy group
 */
void FunctionalMaterial::setDifCoefByGroup(double xs, int group, int position) {

  log_printf(ERROR, "Unable to set the sigma t by group for FunctionalMaterial %d"
             " since this Material subclass doesn't support settig xs by group", _id);
}


/**
 * @brief Set the Material's diffusion coefficient for some energy group.
 * @param xs the diffusion coefficient
 * @param group the energy group
 */
double FunctionalMaterial::getDifCoefByGroup(int group, int position, double temp) {

  if (_num_time_steps == 0)
    log_printf(ERROR, "Unable to get the dif coef for material %d and "
               "group %d since the number of time steps has not been set",
               _id, group);

  if (group < 0 || group >= _num_energy_groups)
    log_printf(ERROR, "Unable to get diffusion coefficient for group %d for "
               "Material %d which contains %d energy groups",
               group, _id, _num_energy_groups);

  int time_step = getTimeStep(position);
  int ng = _num_energy_groups;
  
  double dif_coef_left = _dif_coef[time_step*ng + group];
  double dif_coef_right = _dif_coef[(time_step+1)*ng + group];
  double time = _clock->getTime(position);
  double time_left = _time_steps[time_step];
  double time_right = _time_steps[time_step+1];

  /* Interpolate sigma f */
  double dif_coef = dif_coef_left + (time - time_left) * (dif_coef_right - dif_coef_left)
    / (time_right - time_left);
  
  return dif_coef;
}


/**
 * @brief Set the Material's array of diffusion coefficients.
 * @details This method is a helper function to allow OpenMOC users to assign
 *          the Material's nuclear data in Python. A user must initialize a
 *          NumPy array of the correct size (i.e., a float64 array the length
 *          of the number of energy groups) as input to this function. This
 *          function then fills the NumPy array with the data values for the
 *          Material's buckling coefficients. An example of how this function
 *          might be called in Python is as follows:
 *
 * @code
 *          buckling = numpy.array([0.05, 0.1, 0.15, ... ])
 *          material = openmoc.Material(openmoc.material_id())
 *          material.setBuckling(buckling)
 * @endcode
 *
 * @param xs the array of diffusion coefficents
 * @param num_groups the number of energy groups
 */
void FunctionalMaterial::setVelocity(double* velocity, int num_time_steps, int num_groups) {

  if (_num_time_steps == 0)
    log_printf(ERROR, "Unable to set the velocity for material %d since "
               "the number of time steps has not been set", _id);

  if (_num_energy_groups != num_groups)
    log_printf(ERROR, "Unable to set velocity with %d groups for Material "
               "%d which contains %d energy groups",
               num_groups, _id, _num_energy_groups);

  if (_num_time_steps != num_time_steps)
    log_printf(ERROR, "Unable to set velocity with %d time steps for Material "
               "%d which contains %d time steps",
               num_time_steps, _id, _num_time_steps);    

  for (int c=0; c < _num_time_steps; c++){
    for (int i=0; i < _num_energy_groups; i++)
      _velocity[c*_num_energy_groups+i] = velocity[c*_num_energy_groups+i];
  }
}


/**
 * @brief Set the Material's diffusion coefficient for some energy group.
 * @param xs the diffusion coefficient
 * @param group the energy group
 */
void FunctionalMaterial::setVelocityByGroup(double velocity, int group, int position) {

  log_printf(ERROR, "Unable to set the sigma t by group for FunctionalMaterial %d"
             " since this Material subclass doesn't support settig xs by group", _id);
}


/**
 * @brief Set the Material's diffusion coefficient for some energy group.
 * @param xs the diffusion coefficient
 * @param group the energy group
 */
double FunctionalMaterial::getVelocityByGroup(int group, int position, double temp) {

  if (_num_time_steps == 0)
    log_printf(ERROR, "Unable to get the velocity for material %d and "
               "group %d since the number of time steps has not been set",
               _id, group);

  if (group < 0 || group >= _num_energy_groups)
    log_printf(ERROR, "Unable to set velocity for group %d for "
               "Material %d which contains %d energy groups",
               group, _id, _num_energy_groups);

  int time_step = getTimeStep(position);
  int ng = _num_energy_groups;
  
  double velocity_left = _velocity[time_step*ng + group];
  double velocity_right = _velocity[(time_step+1)*ng + group];
  double time = _clock->getTime(position);
  double time_left = _time_steps[time_step];
  double time_right = _time_steps[time_step+1];

  /* Interpolate sigma f */
  double velocity = velocity_left + (time - time_left) * (velocity_right - velocity_left)
    / (time_right - time_left);
  
  return velocity;
}


void FunctionalMaterial::setDopplerCoefficients(double* xs, int num_groups) {

  if (_num_energy_groups != num_groups)
    log_printf(ERROR, "Unable to set doppler coefficients with %d groups for Material "
               "%d which contains %d energy groups",
               num_groups, _num_energy_groups);

  for (int i=0; i < _num_energy_groups; i++)
    _doppler_coefficients[i] = xs[i];
}


double FunctionalMaterial::getDopplerCoefficientByGroup(int group){

  if (group < 0 || group >= _num_energy_groups)
    log_printf(ERROR, "Unable to get the doppler coefficient for group %d for "
               "Material %d which contains %d energy groups",
               group, _id, _num_energy_groups);

  return _doppler_coefficients[group];  
}


int FunctionalMaterial::getTimeStep(int position){

  if (_clock == NULL){
    log_printf(ERROR, "Unable to get time step for FunctionalMaterial %d"
               " since Clock has not been set", _id);
  }

  double time = _clock->getTime(position);

  int step = 1;
  while (time > _time_steps[step]){
    step++;
  }

  return (step-1);
}


/**
 * @brief Converts this Material's attributes to a character array
 *        representation.
 * @details The character array returned includes the user-defined ID, and each
 *          of the absorption, total, fission, nu multiplied by fission and
 *          scattering cross-sections and chi for all energy groups.
 * @return character array of this Material's attributes
 */
std::string FunctionalMaterial::toString() {

  std::stringstream string;
  int ng = _num_energy_groups;
  
  string << "FunctionalMaterial id = " << _id;
  string << "\n\t Energy per fission = " << _energy_per_fission;
  string << "\n\t Temperature Conversion Factor = " << _temperature_conversion_factor;
  
  for (int step=0; step < _num_time_steps; step++){

    string << "\n\t Time Step = " << _time_steps[step];
    
    if (_sigma_a != NULL) {
      string << "\n\t\tSigma_a = ";
      for (int e = 0; e < _num_energy_groups; e++)
        string << _sigma_a[step*ng+e] << ", ";
    }
    
    if (_sigma_t != NULL) {
      string << "\n\t\tSigma_t = ";
      for (int e = 0; e < _num_energy_groups; e++)
        string << _sigma_t[step*ng+e] << ", ";
    }
    
    if (_sigma_f != NULL) {
      string << "\n\t\tSigma_f = ";
      for (int e = 0; e < _num_energy_groups; e++)
        string << _sigma_f[step*ng+e] << ", ";
    }
    
    if (_nu_sigma_f != NULL) {
      string << "\n\t\tnu_sigma_f = ";
      for (int e = 0; e < _num_energy_groups; e++)
        string << _nu_sigma_f[step*ng+e] << ", ";
    }
    
    if (_sigma_s != NULL) {
      string << "\n\t\tSigma_s = \n\t\t";
      for (int G = 0; G < _num_energy_groups; G++) {
        for (int g = 0; g < _num_energy_groups; g++)
          string << _sigma_s[step*ng*ng+G*ng+g] << "\t\t ";
        string << "\n\t\t";
      }
    }
    
    if (_chi != NULL) {
      string << "\n\t\tChi = ";
      for (int e = 0; e < _num_energy_groups; e++)
        string << _chi[step*ng+e] << ", ";
    }
    
    if (_dif_coef != NULL) {
      string << "\n\t\tDiffusion Coefficient = ";
      for (int e = 0; e < _num_energy_groups; e++)
        string << _dif_coef[step*ng+e] << ", ";
    }
    
    if (_velocity != NULL) {
      string << "\n\t\tVelocity = ";
      for (int e = 0; e < _num_energy_groups; e++)
        string << _velocity[step*ng+e] << ", ";
    }

  }
  
  for (int position=0; position < 8; position++){

    string << "\n\t Clock Position = " << _clock->getPositionName(position);

    if (_precursor_conc != NULL) {
      string << "\n\t\tPrecursor Conc = ";
      for (int e = 0; e < _num_delayed_groups; e++)
        string << _precursor_conc[position*_num_delayed_groups+e] << ", ";
    }
  }

  string << std::endl;
  
  return string.str();
}


/**
 * @brief Create a duplicate of the Material.
 * @return a pointer to the clone
 */
FunctionalMaterial* FunctionalMaterial::clone(){

  FunctionalMaterial* clone = new FunctionalMaterial(getId());

  clone->setTimeSteps(_time_steps, _num_time_steps);
  clone->setNumEnergyGroups(_num_energy_groups);
  clone->setNumDelayedGroups(_num_delayed_groups);
  clone->setEnergyPerFission(_energy_per_fission);
  clone->setTemperatureConversionFactor(_temperature_conversion_factor);
  
  int ng = _num_energy_groups;
  int nt = _num_time_steps;

  clone->setSigmaT(_sigma_t, nt, ng);
  clone->setSigmaA(_sigma_a, nt, ng);
  clone->setSigmaF(_sigma_f, nt, ng);
  clone->setNuSigmaF(_nu_sigma_f, nt, ng);
  clone->setChi(_chi, nt, ng);
  clone->setDifCoef(_dif_coef, nt, ng);
  clone->setVelocity(_velocity, nt, ng);  
  clone->setSigmaS(_sigma_s, nt, ng*ng);
  clone->setDopplerCoefficients(_doppler_coefficients, ng);

  for (int c=0; c < 8; c++){
    for (int d=0; d < _num_delayed_groups; d++)
      clone->setPrecursorConcByGroup(_precursor_conc[c*_num_delayed_groups+d], d, c);
  }

  return clone;
}


void FunctionalMaterial::copy(int position_from, int position_to){

  for (int d=0; d < _num_delayed_groups; d++)
    setPrecursorConcByGroup(_precursor_conc[position_from*_num_delayed_groups+d], d, position_to);
}

