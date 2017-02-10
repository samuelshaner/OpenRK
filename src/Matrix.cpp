#include "Matrix.h"

/**
 * @brief Constructor initializes Matrix as a list of lists
 *        and sets the matrix dimensions.
 * @detail The matrix object uses a "lists of lists" structure (implemented as
 *         a map of lists) to allow for easy setting and incrementing of the
 *         values in the object. When the matrix is needed to perform linear
 *         algebra operations, it is converted to compressed row storage (CSR)
 *         form. The matrix is ordered by cell (as opposed to by group) on the
 *         outside. Locks are used to make the matrix thread-safe against
 *         concurrent writes the same value. One lock locks out multiple rows of
 *         the matrix at a time reprsenting multiple groups in the same cell.
 * @param num_x The number of cells in the x direction.
 * @param num_y The number of cells in the y direction.
 * @param num_z The number of cells in the z direction.
 * @param num_groups The number of energy groups in each cell.
 */
Matrix::Matrix(long int num_cells) {

  setNumCells(num_cells);

  /* Initialize variables */
  for (long int i=0; i < _num_cells; i++)
    _LIL.push_back(std::map<long int, double>());

  _A = NULL;
  _IA = NULL;
  _JA = NULL;
  _DIAG = NULL;
  _modified = true;
  _diags = NULL;
  _num_diags = 0;
}


/**
 * @brief Destructor clears list of lists and deletes the arrays
 *        used to represent the matrix in CSR form.
 */
Matrix::~Matrix() {

  if (_A != NULL)
    delete [] _A;

  if (_IA != NULL)
    delete [] _IA;

  if (_JA != NULL)
    delete [] _JA;

  if (_DIAG != NULL)
    delete [] _DIAG;

  if (_diags != NULL)
    delete [] _diags;

  for (long int i=0; i < _num_cells; i++)
    _LIL[i].clear();
  _LIL.clear();
}


/**
 * @brief Increment a value in the matrix.
 * @detail This method takes a cell and group of origin (cell/group from)
 *         and cell and group of destination (cell/group to) and floating
 *         point value. The origin and destination are used to compute the
 *         row and column in the matrix. If a value exists for the row/column,
 *         the value is incremented by val; otherwise, it is set to val.
 * @param cell_from The origin cell.
 * @param group_from The origin group.
 * @param cell_to The destination cell.
 * @param group_from The destination group.
 * @param val The value used to increment the row/column location.
 */
void Matrix::incrementValue(long int col, long int row, double value) {

  if (col >= _num_cells || col < 0)
    log_printf(ERROR, "Unable to increment Matrix value for col %d"
               " which is not between 0 and %d", col, _num_cells-1);
  else if (row >= _num_cells || row < 0)
    log_printf(ERROR, "Unable to increment Matrix value for row %d"
               " which is not between 0 and %d", row, _num_cells-1);

  _LIL[row][col] += value;

  /* Set global modified flag to true */
  _modified = true;
}


/**
 * @brief Set a value in the matrix.
 * @detail This method takes a cell and group of origin (cell/group from)
 *         and cell and group of destination (cell/group to) and floating
 *         point value. The origin and destination are used to compute the
 *         row and column in the matrix. The location specified by the
 *         row/column is set to val.
 * @param cell_from The origin cell.
 * @param group_from The origin group.
 * @param cell_to The destination cell.
 * @param group_from The destination group.
 * @param val The value used to set the row/column location.
 */
void Matrix::setValue(long int col, long int row, double value) {

  if (col >= _num_cells || col < 0)
    log_printf(ERROR, "Unable to set Matrix value for cell_from %d"
               " which is not between 0 and %d", col, _num_cells-1);
  else if (row >= _num_cells || row < 0)
    log_printf(ERROR, "Unable to set Matrix value for cell_to %d"
               " which is not between 0 and %d", row, _num_cells-1);

  _LIL[row][col] = value;

  /* Set global modified flag to true */
  _modified = true;
}


/**
 * @brief Clear all values in the matrix list of lists.
 */
void Matrix::clear() {
  for (long int i=0; i < _num_cells; i++)
    _LIL[i].clear();

  _modified = true;
}


/**
 * @brief Convert the matrix lists of lists to compressed row (CSR) storage
 *        form.
 */
void Matrix::convertToCSR() {

  /* Get number of nonzero values */
  long int NNZ = getNNZ();

  /* Deallocate memory for arrays if previously allocated */
  if (_A != NULL)
    delete [] _A;

  if (_IA != NULL)
    delete [] _IA;

  if (_JA != NULL)
    delete [] _JA;

  if (_DIAG != NULL)
    delete [] _DIAG;

  /* Allocate memory for arrays */
  _A = new double[NNZ];
  _IA = new long int[_num_cells+1];
  _JA = new long int[NNZ];
  _DIAG = new double[_num_cells];
  std::fill_n(_DIAG, _num_cells, 0.0);

  /* Form arrays */
  long int j = 0;
  std::map<long int, double>::iterator iter;
  for (long int row=0; row < _num_cells; row++) {
    _IA[row] = j;
    for (iter = _LIL[row].begin(); iter != _LIL[row].end(); ++iter) {
      if (iter->second != 0.0) {
        _JA[j] = iter->first;
        _A[j] = iter->second;

        if (row == iter->first)
          _DIAG[row] = iter->second;

        j++;
      }
    }
  }

  _IA[_num_cells] = NNZ;

  /* Reset flat indicating the CSR objects have the same values as the
   * LIL object */
  _modified = false;
}



/**
 * @brief Print the matrix object to the log file.
 */
void Matrix::printString() {

  /* Convert to CSR form */
  convertToCSR();

  std::stringstream string;
  string << std::setprecision(6);
  string << " Matrix Object " << std::endl;
  string << " Num cells: " << _num_cells << std::endl;
  string << " NNZ      : " << getNNZ() << std::endl;

  for (long int row=0; row < _num_cells; row++) {
    for (long int i = _IA[row]; i < _IA[row+1]; i++)
      string << " ( " << row << ", " << _JA[i] << "): " << _A[i] << std::endl;
  }

  string << "End Matrix " << std::endl;
  std::cout << string.str();
}


/**
 * @brief Get a value in the matrix.
 * @detail This method takes a cell and group of origin (cell/group from)
 *         and cell and group of destination (cell/group to).
 *         The origin and destination are used to compute the
 *         row and column in the matrix. The value at the location specified
 *         by the row/column is returned.
 * @param cell_from The origin cell.
 * @param group_from The origin group.
 * @param cell_to The destination cell.
 * @param group_from The destination group.
 * @return The value at the corresponding row/column location.
 */
double Matrix::getValue(long int col, long int row) {

  if (col >= _num_cells || col < 0)
    log_printf(ERROR, "Unable to get Matrix value for col %d"
               " which is not between 0 and %d", col, _num_cells-1);
  else if (row >= _num_cells || row < 0)
    log_printf(ERROR, "Unable to get Matrix value for row %d"
               " which is not between 0 and %d", row, _num_cells-1);

  return _LIL[row][col];
}


/**
 * @brief Get the A component of the CSR form of the matrix object.
 * @return A pointer to the A component of the CSR form matrix object.
 */
double* Matrix::getA() {

  if (_modified)
    convertToCSR();

  return _A;
}


/**
 * @brief Get the IA component of the CSR form of the matrix object.
 * @return A pointer to the IA component of the CSR form matrix object.
 */
long int* Matrix::getIA() {

  if (_modified)
    convertToCSR();

  return _IA;
}


/**
 * @brief Get the JA component of the CSR form of the matrix object.
 * @return A pointer to the JA component of the CSR form matrix object.
 */
long int* Matrix::getJA() {

  if (_modified)
    convertToCSR();

  return _JA;
}


/**
 * @brief Get the diagonal component of the matrix object.
 * @return A pointer to the diagonal component of the matrix object.
 */
double* Matrix::getDiag() {

  if (_modified)
    convertToCSR();

  return _DIAG;
}


/**
 * @brief Get the number of cells in the x dimension.
 * @return The number of cells in the x dimension.
 */
long int Matrix::getNumCells() {
  return _num_cells;
}


/**
 * @brief Get the number of non-zero values in the matrix.
 * @return The number of non-zero values in the matrix.
 */
long int Matrix::getNNZ() {

  long int NNZ = 0;
  std::map<long int, double>::iterator iter;
  for (long int row=0; row < _num_cells; row++) {
    for (iter = _LIL[row].begin(); iter != _LIL[row].end(); ++iter) {
      if (iter->second != 0.0)
        NNZ++;
    }
  }

  return NNZ;
}


/**
 * @brief Set the number of cells in the x dimension.
 * @param num_x The number of cells in the x dimension.
 */
void Matrix::setNumCells(long int num_cells) {

  if (num_cells < 1)
    log_printf(ERROR, "Unable to set Matrix num x to non-positive value %d",
               num_cells);

  _num_cells = num_cells;
}


/**
 * @brief Transpose the matrix in place.
 */
void Matrix::transpose() {

  Matrix temp(_num_cells);
  convertToCSR();
  long int col;
  double val;

  /* Transpose matrix to temp */
  for (long int row=0; row < _num_cells; row++) {
    for (long int i = _IA[row]; i < _IA[row+1]; i++) {
      col = _JA[i];
      val = _A[i];
      temp.setValue(row, col, val);
    }
  }

  /* Copy temp to current matrix */
  clear();
  temp.convertToCSR();
  long int* IA = temp.getIA();
  long int* JA = temp.getJA();
  double* A = temp.getA();

  for (long int row=0; row < _num_cells; row++) {
    for (long int i = IA[row]; i < IA[row+1]; i++) {
      col = JA[i];
      val = A[i];
      setValue(col, row, val);
    }
  }
}


void Matrix::setDiags(int* diags, int num_diags) {

  if (_diags == NULL)
    delete [] _diags;

  /* Create dimensions array */
  _num_diags = num_diags;
  _diags = new int[num_diags];
  std::copy(diags, diags + num_diags, _diags);
}


Matrix* Matrix::diags(Array* array) {

  long int shape[2];
  shape[0] = _num_diags;
  shape[1] = array->getSize() / _num_diags;
  array->reshape(shape, 2);
  double* values = array->getValues();
  long int num_values = array->getShape(1);
  long int row, col;
  for (int d=0; d < _num_diags; d++) {
    if (_diags[d] >= 0) {
      for (row = 0; row < _num_cells - _diags[d]; row++) {
        col = row + _diags[d];
        _LIL[row][col] = values[d * num_values + row];
      }
    }
    else {
      for (row = -_diags[d]; row < _num_cells; row++) {
        col = row + _diags[d];
        _LIL[row][col] = values[d * num_values + col];
      }
    }
  }
}


void Matrix::blockDiags(Array* array, int block_size) {

  int bs2 = block_size * block_size;

  if (array->getShape(0) != _num_cells / block_size)
    log_printf(ERROR, "number of block diags not equal to the 1st dimension "
               "of the array: (%d, %d)", array->getShape(0), _num_cells / block_size);
  if (array->getSize() != _num_cells * block_size)
    log_printf(ERROR, "array size not equal to number of block values"
               " : (%d, %d)", array->getSize(), _num_cells*block_size);

  long int shape[3];
  shape[0] = _num_cells / block_size;
  shape[1] = block_size;
  shape[2] = block_size;
  array->reshape(shape, 3);
  double* values = array->getValues();
  long int row, col;
  for (long int i=0; i < _num_cells / block_size; i++) {
    for (int r=0; r < block_size; r++) {
      row = i*block_size + r;
      for (int c=0; c < block_size; c++) {
        col = i*block_size + c;
        _LIL[row][col] = values[i * bs2 + r * block_size + c];
      }
    }
  }
}


void Matrix::fillWithRandom() {
  for (long int row=0; row < _num_cells; row++) {
    for (long int col=0; col < _num_cells; col++) {
      _LIL[row][col] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
    }
  }
}


/**
 * @brief Performs a matrix vector multiplication.
 * @details This function takes in a Matrix (A), a variable Vector (X),
 *          and a solution Vector (B) and computes the matrix vector product.
 *          The solution Vector is modified in place.
 * @param A a Matrix object
 * @param X the variable Vector object
 * @param B the solution Vector object
 */
void Matrix::scaleByValue(double val) {

  std::map<long int, double>::iterator iter;
  for (long int row=0; row < _num_cells; row++) {
    for (iter = _LIL[row].begin(); iter != _LIL[row].end(); ++iter)
      iter->second *= val;
  }
}


void Matrix::add(Matrix* matrix) {

  if (getNumCells() != matrix->getNumCells())
    log_printf(ERROR, "Cannot add matrices with different sizes");

  long int* IA = matrix->getIA();
  long int* JA = matrix->getJA();
  double* a = matrix->getA();

  for (long int row = 0; row < size; row++) {
    for (long int i = IA[row]; i < IA[row+1]; i++)
      _LIL[row][JA[i]] += a[i];
  }
}


void Matrix::subtract(Matrix* matrix) {

  if (getNumCells() != matrix->getNumCells())
    log_printf(ERROR, "Cannot add matrices with different sizes");

  long int* IA = matrix->getIA();
  long int* JA = matrix->getJA();
  double* a = matrix->getA();

  for (long int row = 0; row < size; row++) {
    for (long int i = IA[row]; i < IA[row+1]; i++)
      _LIL[row][JA[i]] -= a[i];
  }
}
