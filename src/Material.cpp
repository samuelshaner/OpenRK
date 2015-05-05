#include "Material.h"

static int auto_id = 10000;


/**
 * @brief Returns an auto-generated unique Material ID.
 * @details This method is intended as a utility method for user's writing
 *          OpenMOC input files. The method makes use of a static Material
 *          ID which is incremented each time the method is called to enable
 *          unique generation of monotonically increasing IDs. The method's
 *          first ID begins at 10000. Hence, user-defined material IDs greater
 *          than or equal to 10000 is prohibited.
 */
int material_id() {
  int id = auto_id;
  auto_id++;
  return id;
}


/**
 * @brief Constructor sets the ID and unique ID for the Material.
 * @param id the user-defined ID for the material
 */
Material::Material(int id) {

  /* If the user did not define an optional ID, create one */
  if (id == 0)
    _id = material_id();
  
  /* Use the user-defined ID */
  else
    _id = id;

  _sigma_t = NULL;
  _sigma_a = NULL;
  _sigma_s = NULL;
  _sigma_f = NULL;
  _nu_sigma_f = NULL;
  _chi = NULL;
  _dif_coef = NULL;
  _velocity = NULL;
  _precursor_conc = NULL;
  _decay_constant = NULL;
  _delayed_fraction = NULL;
  _energy_per_fission = 0.0;
  _clock = NULL;
  _temperature_conversion_factor = 0.0;
  _num_energy_groups = 0;
  _num_delayed_groups = 0;
  
  _fissionable = false;

  return;
}


/**
 * @brief Destructor deletes all cross-section data structures from memory.
 */
Material::~Material() {

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

  if (_precursor_conc != NULL)
    delete [] _precursor_conc;

  if (_decay_constant != NULL)
    delete [] _decay_constant;

  if (_delayed_fraction != NULL)
    delete [] _delayed_fraction;
}


/**
 * @brief Return the Material's user-defined ID
 * @return the Material's user-defined ID
 */
int Material::getId() const {
  return _id;
}


/**
 * @brief Returns the number of energy groups for this Material's nuclear data.
 * @return the number of energy groups
 */
int Material::getNumEnergyGroups() const {
  return _num_energy_groups;
}


/**
 * @brief Returns the number of delayed groups for this Material's nuclear data.
 * @return the number of delayed groups
 */
int Material::getNumDelayedGroups() const {
  return _num_delayed_groups;
}


/**
 * @brief Return the array of the Material's total cross-sections.
 * @return the pointer to the Material's array of total cross-sections
 */
double* Material::getSigmaT() {
  if (_sigma_t == NULL)
    log_printf(ERROR, "Unable to return Material %d's total "
               "cross-section since it has not yet been set", _id);

  return _sigma_t;
}


/**
 * @brief Return the array of the Material's absorption cross-sections.
 * @return the pointer to the Material's array of absorption cross-sections
 */
double* Material::getSigmaA() {
  if (_sigma_a == NULL)
      log_printf(ERROR, "Unable to return Material %d's absorption "
                 "cross-section since it has not yet been set", _id);

  return _sigma_a;
}


/**
 * @brief Return the array of the Material's scattering cross-section matrix.
 * @return the pointer to the Material's array of scattering cross-sections
 */
double* Material::getSigmaS() {
  if (_sigma_s == NULL)
    log_printf(ERROR, "Unable to return Material %d's scattering "
               "cross-section since it has not yet been set", _id);

  return _sigma_s;
}


/**
 * @brief Return the array of the Material's fission cross-sections.
 * @return the pointer to the Material's array of fission cross-sections
 */
double* Material::getSigmaF() {
  if (_sigma_f == NULL)
    log_printf(ERROR, "Unable to return material %d's fission "
               "cross-section since it has not yet been set", _id);

  return _sigma_f;
}


/**
 * @brief Return the array of the Material's fission cross-sections
 *        multiplied by nu \f$ \nu \f$.
 * @return the pointer to the Material's array of fission cross-sections
 *         multiplied by nu \f$ \nu \f$
 */
double* Material::getNuSigmaF() {
  if (_nu_sigma_f == NULL)
    log_printf(ERROR, "Unable to return Material %d's nu times fission "
               "cross-section since it has not yet been set", _id);

  return _nu_sigma_f;
}


/**
 * @brief Return the array of the Material's chi \f$ \chi \f$.
 * @return the pointer to the Material's array of chi \f$ \chi \f$ values
 */
double* Material::getChi() {
  if (_chi == NULL)
    log_printf(ERROR, "Unable to return Material %d's chi spectrum "
               "since it has not yet been set", _id);

  return _chi;
}


/**
 * @brief Return the array of the Material's diffusion coefficients.
 * @return the pointer to the Material's array of diffusion coefficients
 */
double* Material::getDifCoef() {

  if (_dif_coef == NULL)
    log_printf(ERROR, "Unable to return Material %d's diffusion coefficients "
               "since it has not yet been set", _id);

  return _dif_coef;
}


/**
 * @brief Return the array of the Material's velocities
 * @return the pointer to the Material's array of velocities
 */
double* Material::getVelocity() {

  if (_velocity == NULL)
    log_printf(ERROR, "Unable to return Material %d's velocities "
               "since it has not yet been set", _id);

  return _velocity;
}


/**
 * @brief Return the array of the Material's precursor concentrations.
 * @return the pointer to the Material's array of precursor concentrations
 */
double* Material::getPrecursorConc() {

  if (_precursor_conc == NULL)
    log_printf(ERROR, "Unable to return Material %d's precursor concentrations "
               "since it has not yet been set", _id);

  return _precursor_conc;
}


/**
 * @brief Return the array of the Material's decay constants.
 * @return the pointer to the Material's array of decay constants
 */
double* Material::getDecayConstant() {

  if (_decay_constant == NULL)
    log_printf(ERROR, "Unable to return Material %d's decay constants "
               "since it has not yet been set", _id);

  return _decay_constant;
}


/**
 * @brief Return the array of the Material's delayed fractions.
 * @return the pointer to the Material's array of delayed fractions
 */
double* Material::getDelayedFractions() {

  if (_delayed_fraction == NULL)
    log_printf(ERROR, "Unable to return Material %d's delayed fractions "
               "since it has not yet been set", _id);

  return _delayed_fractions;
}


/**
 * @brief Returns whether or not the Material contains a fissionable (non-zero)
 *        fission cross-section.
 * @return true if fissionable, false otherwise
 */
bool Material::isFissionable() {
  return _fissionable;
}


/**
 * @brief Set the number of energy groups for this Material.
 * @param num_groups the number of energy groups.
 */
void Material::setNumEnergyGroups(const int num_groups) {

  if (num_groups < 0)
    log_printf(ERROR, "Unable to set the number of energy groups for "
               "material %d to %d", _id, num_groups);

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

  /* Allocate memory for data arrays */
  _sigma_t = new double[8*_num_energy_groups];
  _sigma_a = new double[8*_num_energy_groups];
  _sigma_f = new double[8*_num_energy_groups];
  _nu_sigma_f = new double[8*_num_energy_groups];
  _chi = new double[8*_num_energy_groups];
  _sigma_s = new double[8*_num_energy_groups*_num_energy_groups];
  _velocity = new double[8*_num_energy_groups];
  _dif_coef = new double[8*_num_energy_groups];
  
  /* Assign the null vector to each data array */
  memset(_sigma_t, 0.0, 8 * sizeof(double) * _num_energy_groups);
  memset(_sigma_a, 0.0, 8 * sizeof(double) * _num_energy_groups);
  memset(_sigma_f, 0.0, 8 * sizeof(double) * _num_energy_groups);
  memset(_nu_sigma_f, 0.0, 8 * sizeof(double) * _num_energy_groups);
  memset(_chi, 0.0, 8 * sizeof(double) * _num_energy_groups);
  memset(_sigma_s, 0.0, 8 * sizeof(double) * _num_energy_groups * _num_energy_groups);
  memset(_velocity, 0.0, 8 * sizeof(double) * _num_energy_groups);
  memset(_dif_coef, 0.0, 8 * sizeof(double) * _num_energy_groups);
}


/**
 * @brief Set the number of delayed groups for this Material.
 * @param num_groups the number of delayed groups.
 */
void Material::setNumDelayedGroups(const int num_groups) {

  if (num_groups < 0)
    log_printf(ERROR, "Unable to set the number of delayed groups for "
               "material %d to %d", _id, num_groups);

  _num_delayed_groups = num_groups;

  /* Free old data arrays if they were allocated for a previous simulation */
  if (_precursor_conc != NULL)
    delete [] _precursor_conc;

  /* Free old data arrays if they were allocated for a previous simulation */
  if (_decay_constant != NULL)
    delete [] _decay_constant;

  /* Free old data arrays if they were allocated for a previous simulation */
  if (_delayed_fraction != NULL)
    delete [] _delayed_fraction;

  /* Allocate memory for data arrays */
  _precursor_conc = new double[8*_num_delayed_groups];
  _decay_constant = new double[8*_num_delayed_groups];
  _delayed_fraction = new double[8*_num_delayed_groups];
  
  /* Assign the null vector to each data array */
  memset(_precursor_conc, 0.0, 8 * sizeof(double) * _num_delayed_groups);
  memset(_decay_constant, 0.0, 8 * sizeof(double) * _num_delayed_groups);
  memset(_delayed_fraction, 0.0, 8 * sizeof(double) * _num_delayed_groups);
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
void Material::setSigmaT(double* xs, int num_groups) {

  if (_num_energy_groups != num_groups)
    log_printf(ERROR, "Unable to set sigma_t with %d groups for Material "
               "%d which contains %d energy groups", num_groups,  _num_energy_groups);

  for (int c=0; c < 8; c++){
    for (int i=0; i < _num_energy_groups; i++)
      _sigma_t[c*_num_energy_groups+i] = double(xs[i]);
  }
}


/**
 * @brief Set the Material's total cross-section for some energy group.
 * @param xs the total cross-section
 * @param group the energy group
 */
void Material::setSigmaTByGroup(double xs, int group, int position) {

  if (group < 0 || group >= _num_energy_groups)
    log_printf(ERROR, "Unable to set sigma_t for group %d for Material "
               "%d which contains %d energy groups", group, _id, _num_energy_groups);

  _sigma_t[position*_num_energy_groups+group] = xs;
}


/**
 * @brief Set the Material's total cross-section for some energy group.
 * @param xs the total cross-section
 * @param group the energy group
 */
double Material::getSigmaTByGroup(int group, int position, double temp) {

  if (group < 0 || group >= _num_energy_groups)
    log_printf(ERROR, "Unable to get sigma_t for group %d for Material "
               "%d which contains %d energy groups", group, _id, _num_energy_groups);

  return _sigma_t[position*_num_energy_groups+group];
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
void Material::setSigmaA(double* xs, int num_groups) {

  if (_num_energy_groups != num_groups)
    log_printf(ERROR, "Unable to set sigma_a with %d groups for Material "
               "%d which contains %d energy groups", num_groups, _num_energy_groups);

  for (int c=0; c < 8; c++){
    for (int i=0; i < _num_energy_groups; i++)
      _sigma_a[c*_num_energy_groups+i] = double(xs[i]);
  }
}


/**
 * @brief Set the Material's absorption cross-section for some energy group.
 * @param xs the absorption cross-section
 * @param group the energy group
 */
void Material::setSigmaAByGroup(double xs, int group, int position) {

  if (group < 0 || group >= _num_energy_groups)
    log_printf(ERROR, "Unable to set sigma_a for group %d for material "
               "%d which contains %d energy groups", group, _id, _num_energy_groups);

  _sigma_a[position*_num_energy_groups+group] = xs;
}


/**
 * @brief Set the Material's absorption cross-section for some energy group.
 * @param xs the absorption cross-section
 * @param group the energy group
 */
double Material::getSigmaAByGroup(int group, int position, double temp) {

  if (group < 0 || group >= _num_energy_groups)
    log_printf(ERROR, "Unable to get sigma_a for group %d for material "
               "%d which contains %d energy groups", group, _id, _num_energy_groups);

  return _sigma_a[position*_num_energy_groups+group];
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
void Material::setSigmaS(double* xs, int num_groups_squared) {

  if (_num_energy_groups*_num_energy_groups != num_groups_squared)
    log_printf(ERROR, "Unable to set sigma_s with %f groups for Material %d "
               "which contains %d energy groups",
                float(sqrt(num_groups_squared)), _num_energy_groups);

  int ng = _num_energy_groups;
  
  for (int c=0; c < 8; c++){
    for (int i=0; i < ng; i++) {
      for (int j=0; j < ng; j++)
        _sigma_s[c*ng*ng+j*ng+i] = xs[j*ng+i];
    }
  }
}


/**
 * @brief Set the Material's scattering cross-section for some energy group.
 * @param xs the scattering cross-section
 * @param group1 the row index in the scattering matrix
 * @param group2 the column index in the scattering matrix
 */
void Material::setSigmaSByGroup(double xs, int group_from, int group_to, int position) {

  if (group_from < 0 || group_to < 0 || group_from >= _num_energy_groups || group_to >= _num_energy_groups)
    log_printf(ERROR, "Unable to set sigma_s for group %d to %d for Material %d "
               "which contains %d energy groups",
               group_from, group_to, _id, _num_energy_groups);

  int ng = _num_energy_groups;
  _sigma_s[position*ng*ng + ng*(group_from) + (group_to)] = xs;
}


/**
 * @brief Set the Material's scattering cross-section for some energy group.
 * @param xs the scattering cross-section
 * @param group1 the row index in the scattering matrix
 * @param group2 the column index in the scattering matrix
 */
double Material::getSigmaSByGroup(int group_from, int group_to, int position, double temp) {

  if (group_from < 0 || group_to < 0 || group_from >= _num_energy_groups || group_to >= _num_energy_groups)
    log_printf(ERROR, "Unable to get sigma_s for group %d to %d for Material %d "
               "which contains %d energy groups",
               group_from, group_to, _id, _num_energy_groups);

  int ng = _num_energy_groups;
  return _sigma_s[position*ng*ng + ng*(group_from) + (group_to)];
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
void Material::setSigmaF(double* xs, int num_groups) {

  if (_num_energy_groups != num_groups)
    log_printf(ERROR, "Unable to set sigma_f with %d groups for Material "
               "%d which contains %d energy groups", num_groups, _num_energy_groups);

  for (int c=0; c < 8; c++){
    for (int i=0; i < _num_energy_groups; i++)
      _sigma_f[c*_num_energy_groups+i] = xs[i];
  }
  
  /* Determine whether or not this Material is fissionable */
  _fissionable = false;

  for (int i=0; i < _num_energy_groups; i++) {
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
void Material::setSigmaFByGroup(double xs, int group, int position) {

  if (group < 0 || group >= _num_energy_groups)
    log_printf(ERROR, "Unable to set sigma_f for group %d for Material "
               "%d which contains %d energy groups", group, _id, _num_energy_groups);

  _sigma_f[position*_num_energy_groups+group] = xs;

  /* Determine whether or not this Material is fissionable */
  _fissionable = false;

  for (int i=0; i < _num_energy_groups; i++) {
    if (_sigma_f[position*_num_energy_groups+i] > 0.0) {
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
double Material::getSigmaFByGroup(int group, int position, double temp) {

  if (group < 0 || group >= _num_energy_groups)
    log_printf(ERROR, "Unable to get sigma_f for group %d for Material "
               "%d which contains %d energy groups", group, _id, _num_energy_groups);

  return _sigma_f[position*_num_energy_groups+group];
}


/**
 * @brief Set the Material's array of fission cross-sections multiplied by
 *         \f$ \nu \f$
 * @param xs the array of fission cross-sections multiplied by nu
 *        \f$ \nu \f$
 * @param num_groups the number of energy groups
*/
void Material::setNuSigmaF(double* xs, int num_groups) {

  if (_num_energy_groups != num_groups)
    log_printf(ERROR, "Unable to set nu_sigma_f with %d groups for Material %d "
              "which contains %d energy groups", num_groups, _id, _num_energy_groups);

  for (int c=0; c < 8; c++){
    for (int i=0; i < _num_energy_groups; i++)
      _nu_sigma_f[c*_num_energy_groups+i] = xs[i];
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
void Material::setNuSigmaFByGroup(double xs, int group, int position) {

  if (group < 0 || group >= _num_energy_groups)
    log_printf(ERROR, "Unable to set nu_sigma_f for group %d for Material "
               "%d which contains %d energy groups", group, _id, _num_energy_groups);

  _nu_sigma_f[position*_num_energy_groups+group] = xs;
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
double Material::getNuSigmaFByGroup(int group, int position, double temp) {

  if (group < 0 || group >= _num_energy_groups)
    log_printf(ERROR, "Unable to get nu_sigma_f for group %d for Material "
               "%d which contains %d energy groups", group, _id, _num_energy_groups);

  return _nu_sigma_f[position*_num_energy_groups+group];
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
void Material::setChi(double* xs, int num_groups) {

  if (_num_energy_groups != num_groups)
    log_printf(ERROR, "Unable to set chi with %d groups for Material "
               "%d which contains %d energy groups", num_groups, _num_energy_groups);

  for (int c=0; c < 8; c++){
    for (int i=0; i < _num_energy_groups; i++)
      _chi[c*_num_energy_groups+i] = xs[i];
  }
}


/**
 * @brief Set the Material's chi value for some energy group.
 * @param xs the chi value (\f$ \Chi \f$)
 * @param group the energy group
 */
void Material::setChiByGroup(double xs, int group, int position) {

  if (group < 0 || group >= _num_energy_groups)
    log_printf(ERROR, "Unable to set chi for group %d for Material "
              "%d which contains %d energy groups", group, _num_energy_groups, _id);

  _chi[position*_num_energy_groups+group] = xs;
}


/**
 * @brief Set the Material's chi value for some energy group.
 * @param xs the chi value (\f$ \Chi \f$)
 * @param group the energy group
 */
double Material::getChiByGroup(int group, int position, double temp) {

  if (group < 0 || group >= _num_energy_groups)
    log_printf(ERROR, "Unable to get chi for group %d for Material "
               "%d which contains %d energy groups", group, _id, _num_energy_groups);

  return _chi[position*_num_energy_groups+group];
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
void Material::setDifCoef(double* xs, int num_groups) {

  if (_num_energy_groups != num_groups)
    log_printf(ERROR, "Unable to set diffusion coefficient with %d groups for "
               "Material %d which contains %d energy groups", num_groups,
               _num_energy_groups);

  for (int c=0; c < 8; c++){
    for (int i=0; i < _num_energy_groups; i++)
      _dif_coef[c*_num_energy_groups+i] = xs[i];
  }
}


/**
 * @brief Set the Material's diffusion coefficient for some energy group.
 * @param xs the diffusion coefficient
 * @param group the energy group
 */
void Material::setDifCoefByGroup(double xs, int group, int position) {

  if (group < 0 || group >= _num_energy_groups)
    log_printf(ERROR, "Unable to set diffusion coefficient for group %d for "
               "Material %d which contains %d energy groups",
               group, _id, _num_energy_groups);

  _dif_coef[position*_num_energy_groups+group] = xs;
}


/**
 * @brief Set the Material's diffusion coefficient for some energy group.
 * @param xs the diffusion coefficient
 * @param group the energy group
 */
double Material::getDifCoefByGroup(int group, int position, double temp) {

  if (group < 0 || group >= _num_energy_groups)
    log_printf(ERROR, "Unable to get diffusion coefficient for group %d for "
               "Material %d which contains %d energy groups",
               group, _id, _num_energy_groups);

  return _dif_coef[position*_num_energy_groups+group];
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
void Material::setVelocity(double* velocity, int num_groups) {

  if (_num_energy_groups != num_groups)
    log_printf(ERROR, "Unable to set velocity with %d groups for "
               "Material %d which contains %d energy groups", num_groups,
               _num_energy_groups);

  for (int c=0; c < 8; c++){
    for (int i=0; i < _num_energy_groups; i++)
      _velocity[c*_num_energy_groups+i] = velocity[i];
  }
}


/**
 * @brief Set the Material's diffusion coefficient for some energy group.
 * @param xs the diffusion coefficient
 * @param group the energy group
 */
void Material::setVelocityByGroup(double velocity, int group, int position) {

  if (group < 0 || group >= _num_energy_groups)
    log_printf(ERROR, "Unable to set velocity for group %d for "
               "Material %d which contains %d energy groups",
               group, _id, _num_energy_groups);

  _velocity[position*_num_energy_groups+group] = velocity;
}


/**
 * @brief Set the Material's diffusion coefficient for some energy group.
 * @param xs the diffusion coefficient
 * @param group the energy group
 */
double Material::getVelocityByGroup(int group, int position, double temp) {

  if (group < 0 || group >= _num_energy_groups)
    log_printf(ERROR, "Unable to set velocity for group %d for "
               "Material %d which contains %d energy groups",
               group, _id, _num_energy_groups);

  return _velocity[position*_num_energy_groups+group];
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
void Material::setPrecursorConc(double* precursor_conc, int num_groups) {

  if (_num_delayed_groups != num_groups)
    log_printf(ERROR, "Unable to set precursor conc with %d groups for "
               "Material %d which contains %d delayed groups", num_groups,
               _num_delayed_groups);

  for (int c=0; c < 8; c++){
    for (int i=0; i < _num_delayed_groups; i++)
      _precursor_conc[c*_num_delayed_groups+i] = precursor_conc[i];
  }
}


/**
 * @brief Set the Material's diffusion coefficient for some energy group.
 * @param xs the diffusion coefficient
 * @param group the energy group
 */
void Material::setPrecursorConcByGroup(double conc, int group, int position) {

  if (group < 0 || group >= _num_delayed_groups)
    log_printf(ERROR, "Unable to set precursor conc for group %d for "
               "Material %d which contains %d delayed groups",
               group, _id, _num_delayed_groups);

  _precursor_conc[position*_num_delayed_groups+group] = conc;
}


/**
 * @brief Set the Material's diffusion coefficient for some energy group.
 * @param xs the diffusion coefficient
 * @param group the energy group
 */
double Material::getPrecursorConcByGroup(int group, int position, double temp) {

  if (group < 0 || group >= _num_delayed_groups)
    log_printf(ERROR, "Unable to get precursor conc for group %d for "
               "Material %d which contains %d delayed groups",
               group, _id, _num_delayed_groups);

  return _precursor_conc[position*_num_delayed_groups+group];
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
void Material::setDecayConstant(double* decay_constant, int num_groups) {

  if (_num_delayed_groups != num_groups)
    log_printf(ERROR, "Unable to set decay constant with %d groups for "
               "Material %d which contains %d delayed groups", num_groups,
               _num_delayed_groups);

  for (int c=0; c < 8; c++){
    for (int i=0; i < _num_delayed_groups; i++)
      _decay_constant[c*_num_delayed_groups+i] = decay_constant[i];
  }
}


/**
 * @brief Set the Material's diffusion coefficient for some energy group.
 * @param xs the diffusion coefficient
 * @param group the energy group
 */
void Material::setDecayConstantByGroup(double decay_constant, int group, int position) {

  if (group < 0 || group >= _num_delayed_groups)
    log_printf(ERROR, "Unable to set decay constant for group %d for "
               "Material %d which contains %d delayed groups",
               group, _id, _num_delayed_groups);

  _decay_constant[position*_num_delayed_groups+group] = decay_constant;
}


/**
 * @brief Set the Material's diffusion coefficient for some energy group.
 * @param xs the diffusion coefficient
 * @param group the energy group
 */
double Material::getDecayConstantByGroup(int group, int position, double temp) {

  if (group < 0 || group >= _num_delayed_groups)
    log_printf(ERROR, "Unable to get decay constant for group %d for "
               "Material %d which contains %d delayed groups",
               group, _id, _num_delayed_groups);

  return _decay_constant[position*_num_delayed_groups+group];
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
void Material::setDelayedFraction(double* delayed_fraction, int num_groups) {

  if (_num_delayed_groups != num_groups)
    log_printf(ERROR, "Unable to set delayed fractions with %d groups for "
               "Material %d which contains %d delayed groups", num_groups,
               _num_delayed_groups);

  for (int c=0; c < 8; c++){
    for (int i=0; i < _num_delayed_groups; i++)
      _delayed_fraction[c*_num_delayed_groups+i] = delayed_fraction[i];
  }
}


/**
 * @brief Set the Material's diffusion coefficient for some energy group.
 * @param xs the diffusion coefficient
 * @param group the energy group
 */
void Material::setDelayedFractionByGroup(double delayed_fraction, int group, int position) {

  if (group < 0 || group >= _num_delayed_groups)
    log_printf(ERROR, "Unable to set delayed fraction for group %d for "
               "Material %d which contains %d delayed groups",
               group, _id, _num_delayed_groups);

  _delayed_fraction[position*_num_delayed_groups+group] = delayed_fraction;
}


/**
 * @brief Set the Material's diffusion coefficient for some energy group.
 * @param xs the diffusion coefficient
 * @param group the energy group
 */
double Material::getDelayedFractionByGroup(int group, int position, double temp) {

  if (group < 0 || group >= _num_delayed_groups)
    log_printf(ERROR, "Unable to get delayed fraction for group %d for "
               "Material %d which contains %d delayed groups",
               group, _id, _num_delayed_groups);

  return _delayed_fraction[position*_num_delayed_groups+group];
}


/**
 * @brief Converts this Material's attributes to a character array
 *        representation.
 * @details The character array returned includes the user-defined ID, and each
 *          of the absorption, total, fission, nu multiplied by fission and
 *          scattering cross-sections and chi for all energy groups.
 * @return character array of this Material's attributes
 */
std::string Material::toString() {

  std::stringstream string;
  int ng = _num_energy_groups;
  
  string << "Material id = " << _id;
  string << "\n\t Energy per fission = " << _energy_per_fission;
  string << "\n\t Temperature Conversion Factor = " << _temperature_conversion_factor;
  
  for (int position=0; position < 8; position++){

    string << "\n\t Clock Position = " << _clock->getPositionName(position);
    
    if (_sigma_a != NULL) {
      string << "\n\t\tSigma_a = ";
      for (int e = 0; e < _num_energy_groups; e++)
        string << _sigma_a[position*ng+e] << ", ";
    }
    
    if (_sigma_t != NULL) {
      string << "\n\t\tSigma_t = ";
      for (int e = 0; e < _num_energy_groups; e++)
        string << _sigma_t[position*ng+e] << ", ";
    }
    
    if (_sigma_f != NULL) {
      string << "\n\t\tSigma_f = ";
      for (int e = 0; e < _num_energy_groups; e++)
        string << _sigma_f[position*ng+e] << ", ";
    }
    
    if (_nu_sigma_f != NULL) {
      string << "\n\t\tnu_sigma_f = ";
      for (int e = 0; e < _num_energy_groups; e++)
        string << _nu_sigma_f[position*ng+e] << ", ";
    }
    
    if (_sigma_s != NULL) {
      string << "\n\t\tSigma_s = \n\t\t";
      for (int G = 0; G < _num_energy_groups; G++) {
        for (int g = 0; g < _num_energy_groups; g++)
          string << _sigma_s[position*ng*ng+G*ng+g] << "\t\t ";
        string << "\n\t\t";
      }
    }
    
    if (_chi != NULL) {
      string << "\n\t\tChi = ";
      for (int e = 0; e < _num_energy_groups; e++)
        string << _chi[position*ng+e] << ", ";
    }
    
    if (_dif_coef != NULL) {
      string << "\n\t\tDiffusion Coefficient = ";
      for (int e = 0; e < _num_energy_groups; e++)
        string << _dif_coef[position*ng+e] << ", ";
    }
    
    if (_velocity != NULL) {
      string << "\n\t\tVelocity = ";
      for (int e = 0; e < _num_energy_groups; e++)
        string << _velocity[position*ng+e] << ", ";
    }
    
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
 * @brief Prints a string representation of all of the Material's attributes to
 *        the console.
 */
void Material::printString() {
  std::cout << toString().c_str();
}


/**
 * @brief Create a duplicate of the Material.
 * @return a pointer to the clone
 */
Material* Material::clone(){

  Material* clone = new Material(getId());

  clone->setNumEnergyGroups(_num_energy_groups);
  clone->setNumDelayedGroups(_num_delayed_groups);
  clone->setEnergyPerFission(_energy_per_fission);
  clone->setTemperatureConversionFactor(_temperature_conversion_factor);
  
  int ng = _num_energy_groups;
  
  for (int c=0; c < 8; c++){
    for (int i=0; i < ng; i++) {
      clone->setSigmaTByGroup(_sigma_t[c*ng+i], i, c);
      clone->setSigmaAByGroup(_sigma_a[c*ng+i], i, c);
      clone->setSigmaFByGroup(_sigma_f[c*ng+i], i, c);
      clone->setNuSigmaFByGroup(_nu_sigma_f[c*ng+i], i, c);
      clone->setChiByGroup(_chi[c*ng+i], i, c);
      clone->setDifCoefByGroup(_dif_coef[c*ng+i], i, c);
      clone->setVelocityByGroup(_velocity[c*ng+i], i, c);
      
      for (int j=0; j < _num_energy_groups; j++)
        clone->setSigmaSByGroup(_sigma_s[c*ng*ng + i*ng + j], i, j, c);
    }

    for (int d=0; d < _num_delayed_groups; d++)
      clone->setPrecursorConcByGroup(_precursor_conc[c*_num_delayed_groups+d], d, c);
  }

  return clone;
}


void Material::copy(int position_from, int position_to){

  int ng = _num_energy_groups;
  
  for (int i=0; i < ng; i++) {
    setSigmaTByGroup(_sigma_t[position_from*ng+i], i, position_to);
    setSigmaAByGroup(_sigma_a[position_from*ng+i], i, position_to);
    setSigmaFByGroup(_sigma_f[position_from*ng+i], i, position_to);
    setNuSigmaFByGroup(_nu_sigma_f[position_from*ng+i], i, position_to);
    setChiByGroup(_chi[position_from*ng+i], i, position_to);
    setDifCoefByGroup(_dif_coef[position_from*ng+i], i, position_to);
    setVelocityByGroup(_velocity[position_from*ng+i], i, position_to);
    
    for (int j=0; j < _num_energy_groups; j++)
      setSigmaSByGroup(_sigma_s[position_from*ng*ng + i*ng + j], i, j, position_to);
  }
  
  for (int d=0; d < _num_delayed_groups; d++)
    setPrecursorConcByGroup(_precursor_conc[position_from*_num_delayed_groups+d], d, position_to);
}


void Material::setEnergyPerFission(double energy_per_fission){
  _energy_per_fission = energy_per_fission;
}


double Material::getEnergyPerFission(){
  return _energy_per_fission;
}


void Material::setClock(Clock* clock){
  _clock = clock;
}


void Material::setTemperatureConversionFactor(double conversion_factor){
  _temperature_conversion_factor = conversion_factor;
}


double Material::getTemperatureConversionFactor(){
  return _temperature_conversion_factor;
}
