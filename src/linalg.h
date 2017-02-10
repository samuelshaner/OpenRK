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
#include <omp.h>
#include "Matrix.h"
#include "Array.h"
#include "constants.h"
#endif

double eigenvalueSolve(Matrix* A, Matrix* M, Array* X, double tol,
                             double SOR_factor=1.5);
Array* linearSolve(Matrix* A, Matrix* M, Array* B, Array* X=NULL, double tol=1.e-6,
                 double SOR_factor=1.5);
Array* matrixMultiplication(Matrix* A, Array* X, Array* B=NULL);
double computeRMSE(Array* x, Array* y);
void setNumThreads(int num_threads);

#endif /* LINALG_H_ */
