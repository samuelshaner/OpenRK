/**
 * @file linalg.h
 * @details This file contains a library of functions for performing linear
 *          algebra operations.
 * @date July 3, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

/* File: linalg.h */

#ifndef LINALG_H_
#define LINALG_H_

#ifdef __cplusplus
#include <math.h>
#include "log.h"
#include <omp.h>
#include "Matrix.h"
#include "Vector.h"
#include "constants.h"
#endif

double eigenvalueSolve(Matrix* A, Matrix* M, Vector* X, double tol,
                             double SOR_factor=1.5);
void linearSolve(Matrix* A, Matrix* M, Vector* X, Vector* B, double tol,
                 double SOR_factor=1.5);
void matrixMultiplication(Matrix* A, Vector* X, Vector* B);
double computeRMSE(Vector* x, Vector* y, bool integrated);
void setNumThreads(int num_threads);

#endif /* LINALG_H_ */
