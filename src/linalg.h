
/* File: linalg.h */

#ifndef LINALG_H_
#define LINALG_H_

#ifdef __cplusplus
#include<math.h>
#include "log.h"
#include "pairwise_sum.h"
#include<omp.h>
#include "Matrix.h"
#include "Vector.h"
#endif

double eigenvalueSolve(Matrix* A, Matrix* M, Vector* X, double tol);
void linearSolve(Matrix* A, Vector* X, Vector* B, double tol);
void matrixMultiplication(Matrix* A, Vector* X, Vector* B);
void setNumThreads(int num_threads);

#endif /* LINALG_H_ */
