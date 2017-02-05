#include "Array.h"

/**
 * @brief Constructor initializes Array object as a floating point array
 *        and sets the array dimensions.
 * @detail The array has arbitrary dimensions that are set in the constructor.
 * @param dimensions An integer array of dimensions.
 * @param num_dimensions The number of dimensions.
 */
Array::Array(long int* dimensions, int num_dimensions) {

  _values = NULL;
  _dimensions = NULL;

  //srand (static_cast <unsigned> (time(0)));

  setDimensions(dimensions, num_dimensions);
}


/**
 * @brief Destructor deletes the arrays used to represent the array.
 */
Array::~Array() {
  clearObjects();
}


void Array::clearObjects() {

  if (_values != NULL)
    delete [] _values;

  if (_dimensions != NULL)
    delete [] _dimensions;
}


/**
 * @brief Print the array object to the log file.
 */
void Array::printString() {

  std::stringstream string;
  string << std::setprecision(6);
  string << "Array" << std::endl;
  string << " Dimensions: ( ";

  for (int d=0; d < _num_dimensions; d++) {
    string << _dimensions[d];
    if (d < _num_dimensions-1)
      string << ", ";
  }

  string << ")" << std::endl;

  long int indices[_num_dimensions];
  for (long int i=0; i < getSize(); i++) {
    getIndices(i, indices);

    /* Print the indices */
    string << "( ";
    for (int d=0; d < _num_dimensions; d++) {
      string << indices[d];
      if (d < _num_dimensions-1)
        string << ", ";
    }

    string << "): ";

    /* Print the value */
    string << _values[i] << std::endl;
  }

  string << "End Array" << std::endl;
  std::cout << string.str();
}


void Array::setDimensions(long int* dimensions, int num_dimensions) {

  /* Clear objects */
  clearObjects();

  /* Create dimensions array */
  _num_dimensions = num_dimensions;
  _dimensions = new long int[num_dimensions];
  std::copy(dimensions, dimensions + num_dimensions, _dimensions);

  /* Create new array and initialize to zero */
  _values = new double[getSize()];
  setAll(0.);
}


void Array::setAll(double val) {
  std::fill_n(_values, getSize(), val);
}


long int Array::getSize() {

  long int size = 1;
  for (int d=0; d < _num_dimensions; d++)
    size *= _dimensions[d];

  return size;
}


void Array::getIndices(long int index, long int* indices) {

  long int size = 1;
  long int offset;
  for (int d=0; d < _num_dimensions; d++) {
    offset = 1;
    for (int e=d+1; e < _num_dimensions; e++)
      offset *= _dimensions[e];
    indices[d] = index / offset;
    index = index % offset;
  }
}


long int Array::getIndex(long int* dimensions, int num_dimensions) {

  long int index = 0;
  long int offset;
  for (int d=0; d < _num_dimensions; d++) {
    offset = 1;
    for (int e=d+1; e < _num_dimensions; e++)
      offset *= _dimensions[e];
    index += offset * dimensions[d];
  }

  return index;
}


void Array::clear() {
  setAll(0.);
}


void Array::scaleByValue(double val) {
  for (long int i=0; i < getSize(); i++)
    _values[i] *= val;
}


int Array::getNumDimensions() {
  return _num_dimensions;
}


long int Array::getShape(int dimension) {
  return _dimensions[dimension];
}


void Array::multiply(Array* array, Array* result) {

  if (getSize() != array->getSize() ||
      getSize() != result->getSize())
    log_printf(ERROR, "Cannot multiply arrays with different sizes");

  double* array_values = array->getValues();
  double* result_values = result->getValues();

  for (long int i=0; i < getSize(); i++)
    result_values[i] = _values[i] * array_values[i];
}


void Array::divide(Array* array, Array* result) {

  if (getSize() != array->getSize() ||
      getSize() != result->getSize())
    log_printf(ERROR, "Cannot divide arrays with different sizes");

  double* array_values = array->getValues();
  double* result_values = result->getValues();

  for (long int i=0; i < getSize(); i++) {
    if (array_values[i] == 0.)
      result_values[i] = _values[i] / 1.e-12;
    else
      result_values[i] = _values[i] / array_values[i];
  }
}


void Array::add(Array* array, Array* result) {

  if (getSize() != array->getSize() ||
      getSize() != result->getSize())
    log_printf(ERROR, "Cannot add arrays with different sizes");

  double* array_values = array->getValues();
  double* result_values = result->getValues();

  for (long int i=0; i < getSize(); i++)
    result_values[i] = _values[i] + array_values[i];
}


void Array::subtract(Array* array, Array* result) {

  if (getSize() != array->getSize() ||
      getSize() != result->getSize())
    log_printf(ERROR, "Cannot multiply arrays with different sizes");

  double* array_values = array->getValues();
  double* result_values = result->getValues();

  for (long int i=0; i < getSize(); i++)
    result_values[i] = _values[i] - array_values[i];
}


double* Array::getValues() {
  return _values;
}


void Array::outputValues(double* np_array, long int num_values) {
  std::copy(_values, _values + num_values, np_array);
}


/**
 * @brief Get the sum of all the values in the array.
 * @return The sum of all the values in the array.
 */
double Array::getSum() {
  return pairwise_sum(_values, getSize());
}


double Array::getValue(long int* dimensions, int num_dimensions) {
  return _values[getIndex(dimensions, num_dimensions)];
}


void Array::setValue(long int index, double value) {
  _values[index] = value;
}


void Array::incrementValue(long int index, double value) {
  _values[index] += value;
}


Array* Array::sumAxis(int axis) {

  /* Create new dimensions array */
  int num_dimensions = _num_dimensions - 1;
  long int dimensions[num_dimensions];

  for (int d=0; d < _num_dimensions; d++) {
    if (d < axis)
      dimensions[d] = _dimensions[d];
    else if (d > axis)
      dimensions[d-1] = _dimensions[d];
  }

  /* Create a new Array object */
  Array* new_array = new Array(dimensions, num_dimensions);

  /* Sum over axis */
  long int indices[_num_dimensions];
  long int new_indices[num_dimensions];
  long int new_index;
  for (long int i=0; i < getSize(); i++) {
    getIndices(i, indices);

    /* Get the new indices */
    for (int d=0; d < _num_dimensions; d++) {
      if (d < axis)
        new_indices[d] = indices[d];
      else if (d > axis)
        new_indices[d-1] = indices[d];
    }

    /* Increment the new array values */
    new_index = new_array->getIndex(new_indices, num_dimensions);
    new_array->incrementValue(new_index, _values[i]);
  }

  return new_array;
}


Array* Array::tile(int factor) {

  /* Create new dimensions array */
  long int dimensions[_num_dimensions];
  std::copy(_dimensions, _dimensions + _num_dimensions, dimensions);
  dimensions[_num_dimensions-1] = _dimensions[_num_dimensions-1] * factor;

  /* Create a new Array object */
  Array* new_array = new Array(dimensions, _num_dimensions);

  /* Sum over axis */
  long int indices[_num_dimensions];
  long int new_index;
  for (long int i=0; i < getSize(); i++) {
    getIndices(i, indices);

    for (int f=0; f<factor; f++) {
      new_index = new_array->getIndex(indices, _num_dimensions);
      new_array->setValue(new_index, _values[i]);
      indices[_num_dimensions-1] += _dimensions[_num_dimensions-1];
    }
  }

  return new_array;
}


Array* Array::repeat(int factor) {

  /* Create new dimensions array */
  long int dimensions[_num_dimensions];
  std::copy(_dimensions, _dimensions + _num_dimensions, dimensions);
  dimensions[_num_dimensions-1] = _dimensions[_num_dimensions-1] * factor;

  /* Create a new Array object */
  Array* new_array = new Array(dimensions, _num_dimensions);

  /* Sum over axis */
  long int indices[_num_dimensions];
  long int new_index;
  for (long int i=0; i < getSize(); i++) {
    getIndices(i, indices);
    indices[_num_dimensions-1] *= factor;

    for (int f=0; f<factor; f++) {
      new_index = new_array->getIndex(indices, _num_dimensions);
      new_array->setValue(new_index, _values[i]);
      indices[_num_dimensions-1]++;
    }
  }

  return new_array;
}


void Array::reshape(long int* dimensions, int num_dimensions) {

  long int old_size = getSize();
  long int new_size = 1;
  for (int d=0; d < num_dimensions; d++)
    new_size *= dimensions[d];

  if (old_size != new_size)
    log_printf(NORMAL, "Cannot reshape array since new size of %d is not equal"
               "to old size of %d", new_size, old_size);

  /* Create new dimensions array */
  delete [] _dimensions;
  _num_dimensions = num_dimensions;
  _dimensions = new long int[_num_dimensions];
  std::copy(dimensions, dimensions + num_dimensions, _dimensions);
}


void Array::flatten() {

  delete [] _dimensions;
  _num_dimensions = 1;
  _dimensions = new long int[1];
  _dimensions[0] = getSize();
}


void Array::copyTo(Array* array) {

  if (_num_dimensions != array->getNumDimensions())
    log_printf(NORMAL, "Cannot copy array since the copy has a different number"
               " of dimensions (%d, %d)", _num_dimensions, array->getNumDimensions());

  for (int d=0; d < _num_dimensions; d++) {
    if (_dimensions[d] != array->getShape(d))
    log_printf(NORMAL, "Cannot copy array since dimension %d does not agree:"
               " (%d, %d)", d, _dimensions[d], array->getShape(d));
  }

  /* Copy values */
  std::copy(_values, _values + getSize(), array->getValues());
}


void Array::setValues(double* values, long int num_values) {

  if (num_values != getSize())
    log_printf(NORMAL, "Cannot set array values since the values have a different number"
               " of values (%d, %d)", num_values, getSize());

  std::copy(values, values + num_values, _values);
}


void Array::fillWithRandom() {
  for (long int i=0; i < getSize(); i++)
    _values[i] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
}
