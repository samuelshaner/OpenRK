#include "linalg.h"

/**
 * @brief Solves a generalized eigenvalue problem using the Power method.
 * @details This function takes in a loss + streaming Matrix (A),
 *          a fission gain Matrix (M), a flux Vector (X), a tolerance used
 *          for both the power method and linear solve convergence (tol), and
 *          a successive over-relaxation factor (SOR_factor) and computes the
 *          dominant eigenvalue and eigenvector using the Power method. The
 *          eigenvalue is returned and the input X Vector is modified in
 *          place to be the corresponding eigenvector.
 * @param A the loss + streaming Matrix object
 * @param M the fission gain Matrix object
 * @param X the flux Vector object
 * @param tol the power method and linear solve source convergence threshold
 * @param SOR_factor the successive over-relaxation factor
 * @return k_eff the dominant eigenvalue
 */
double eigenvalueSolve(Matrix* A, Matrix* M, Array* X, double tol,
                       double SOR_factor) {

  log_printf(INFO, "Computing the Matrix-Vector eigenvalue...");

  long int size = X->getSize();

  /* Check for consistency of matrix and vector dimensions */
  if (A->getNumCells() != M->getNumCells() || A->getNumCells() != size)
    log_printf(ERROR, "Cannot compute the Matrix-Vector eigenvalue with "
               "different x dimensions for the A matrix, M matrix, and X vector"
               ": (%d, %d, %d)", A->getNumCells(), M->getNumCells(), size);

  /* Initialize variables */
  long int dimensions[1];
  dimensions[0] = size;
  Array* old_source = new Array(dimensions, 1);
  Array* new_source = new Array(dimensions, 1);
  double residual, _k_eff;
  int iter;

  /* Compute and normalize the initial source */
  old_source = matrixMultiplication(M, X, old_source);
  old_source->scaleByValue(size / old_source->sum());
  X->scaleByValue(size / old_source->sum());

  /* Power iteration Matrix-Vector solver */
  for (iter = 0; iter < MAX_LINALG_POWER_ITERATIONS; iter++) {

    /* Solve X = A^-1 * old_source */
    linearSolve(A, M, old_source, X, tol*1e1, SOR_factor);

    /* Compute the new source */
    new_source = matrixMultiplication(M, X, new_source);

    /* Compute and set keff */
    _k_eff = new_source->sum() / size;

    /* Scale the old source by keff */
    old_source->scaleByValue(_k_eff);

    /* Compute the residual */
    residual = computeRMSE(new_source, old_source);

    /* Normalize the new source to have an average value of 1.0 */
    new_source->scaleByValue(size / new_source->sum());
    new_source->copyTo(old_source);

    log_printf(NORMAL, "Eigen iter: %d, keff: %f, res: %f",
               iter, _k_eff, residual);

    /* Check for convergence */
    if (residual < tol && iter > 10)
      break;
  }

  log_printf(INFO, "Matrix-Vector eigenvalue solve iterations: %d", iter);

  delete new_source;
  delete old_source;

  return _k_eff;
}


/**
 * @brief Solves a linear system using Red-Black Gauss Seidel with
 *        successive over-relaxation.
 * @details This function takes in a loss + streaming Matrix (A),
 *          a fission gain Matrix (M), a flux Vector (X), a source Vector (B),
 *          a source convergence tolerance (tol) and a successive
 *          over-relaxation factor (SOR_factor) and computes the
 *          solution to the linear system. The input X Vector is modified in
 *          place to be the solution vector.
 * @param A the loss + streaming Matrix object
 * @param M the fission gain Matrix object
 * @param X the flux Vector object
 * @param B the source Vector object
 * @param tol the power method and linear solve source convergence threshold
 * @param SOR_factor the successive over-relaxation factor
 */
Array* linearSolve(Matrix* A, Matrix* M, Array* B, Array* X, double tol,
                 double SOR_factor) {

  if (X == NULL) {
    X = B->copy();
    X->setAll(1.0);
  }

  long int size = X->getSize();

  /* Check for consistency of matrix and vector dimensions */
  if (A->getNumCells() != B->getSize() || A->getNumCells() != size ||
      A->getNumCells() != M->getNumCells())
    log_printf(ERROR, "Cannot perform linear solve with different x dimensions"
               " for the A matrix, M matrix, B vector, and X vector: "
               "(%d, %d, %d, %d)", A->getNumCells(), M->getNumCells(),
               B->getSize(), size);

  /* Initialize variables */
  double residual;
  int iter = 0;
  long int dimensions[1];
  dimensions[0] = size;
  Array* X_old = new Array(dimensions, 1);
  double* x_old = X_old->getValues();
  long int* IA = A->getIA();
  long int* JA = A->getJA();
  double* DIAG = A->getDiag();
  double* a = A->getA();
  double* x = X->getValues();
  double* b = B->getValues();
  long int col;
  Array* old_source = new Array(dimensions, 1);
  Array* new_source = new Array(dimensions, 1);

  /* Compute initial source */
  old_source = matrixMultiplication(M, X, old_source);

  while (iter < MAX_LINEAR_SOLVE_ITERATIONS) {

    /* Pass new flux to old flux */
    X->copyTo(X_old);

    /* Iteration over red/black cells */
    for (long int row = 0; row < size; row++) {

      /* Over-relax the x array */
      x[row] = (1.0 - SOR_factor) * x[row];

      for (long int i = IA[row]; i < IA[row+1]; i++) {

        /* Get the column index */
        col = JA[i];

        if (row == col)
          x[row] += SOR_factor * b[row] / DIAG[row];
        else
          x[row] -= SOR_factor * a[i] * x[col] / DIAG[row];
      }
    }

    /* Compute the new source */
    new_source = matrixMultiplication(M, X, new_source);

    /* Compute the residual */
    residual = computeRMSE(new_source, old_source);

    /* Copy the new source to the old source */
    new_source->copyTo(old_source);

    /* Increment the interations counter */
    iter++;

    log_printf(INFO, "SOR iter: %d, residual: %g", iter, residual);

    if (residual < tol && iter > 10)
      break;
  }

  delete old_source;
  delete new_source;
  delete X_old;

  log_printf(INFO, "linear solve iterations: %d", iter);

  return X;
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
Array* matrixMultiplication(Matrix* A, Array* X, Array* B) {

  if (B == NULL)
    B = X->copy();

  /* Check for consistency of matrix and vector dimensions */
  if (A->getNumCells() != B->getSize() || A->getNumCells() != X->getSize())
    log_printf(ERROR, "Cannot perform matrix multiplication  with different x "
               "dimensions for the A matrix, B vector, and X vector: "
               "(%d, %d, %d)", A->getNumCells(), B->getSize(), X->getSize());

  B->setAll(0.0);
  long int* IA = A->getIA();
  long int* JA = A->getJA();
  double* a = A->getA();
  double* x = X->getValues();
  double* b = B->getValues();
  long int size = X->getSize();

  #pragma omp parallel for
  for (long int row = 0; row < size; row++) {
    for (long int i = IA[row]; i < IA[row+1]; i++)
      b[row] += a[i] * x[JA[i]];
  }

  return B;
}


/**
 * @brief Computes the Root Mean Square Error of two Vectors.
 * @details This function takes in two vectors (X and Y) and computes the
 *          Root Mean Square Error of the Vector Y with respect to Vector X.
 *          The boolean integrated must also be given to indicate whether the
 *          operation on the vector should be group-wise integrated before
 *          performing the RMSE operation.
 * @param X a Vector object
 * @param Y a second Vector object
 * @param integrated a boolean indicating whether to group-wise integrate.
 */
double computeRMSE(Array* X, Array* Y) {

  /* Check for consistency of vector dimensions */
  if (X->getSize() != Y->getSize())
    log_printf(ERROR, "Cannot compute RMSE with different vector dimensions: "
               "(%d) and (%d)", X->getSize(), Y->getSize());

  double rmse;
  long int size = X->getSize();

  long int dimensions[1];
  dimensions[0] = size;
  Array residual(dimensions, 1);
  double* x_values = X->getValues();
  double* y_values = Y->getValues();

  /* Compute the RMSE */
#pragma omp parallel for
  for (long int i = 0; i < size; i++) {
    if (x_values[i] != 0.0)
      residual.setValue(i, pow((x_values[i] - y_values[i]) / x_values[i], 2));
  }

  rmse = sqrt(residual.sum() / size);

  return rmse;
}


void setNumThreads(int num_threads) {
  omp_set_num_threads(num_threads);
}
