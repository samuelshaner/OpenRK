
/* File: linalg.h */

#ifndef LINALG_H_
#define LINALG_H_

#ifdef __cplusplus
#include<math.h>
#include "log.h"
#include "pairwise_sum.h"
#include<omp.h>
#include "gpu_linalg.h"
#endif

/** Indexing macro for the scalar flux in each FSR and energy group */
#define A(r,e) (A[(r)*ng*(ng+6) + (e)])

double eigenvalueSolve(double *A, double *M, double* flux, double* old_source, double* new_source,
                       double* flux_temp, int ng, int cx, int cy, int cz, double tol);

void linearSolve(double *A, double* flux, double* source, double* flux_temp,
                 int cx, int cy, int cz,int ng, double tol);

void matrix_multiplication(double *matrix, double* vector_x, 
                           double* vector_y, int num_blocks, 
                           int block_width);

void vector_scale(double* vector, double scale_value, int length);

void setNumThreads(int num_threads);

#endif /* LINALG_H_ */
