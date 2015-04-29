
/* File: linalg.h */

#ifndef LINALG_H_
#define LINALG_H_

#ifdef __cplusplus
#include<math.h>
#include "log.h"
#include "pairwise_sum.h"
#include<omp.h>
#endif

double eigenvalueSolve(double *A, int A1, int A2, 
                       double *M, int M1, int M2,
                       double* flux, int flux1,
                       double* old_source, int old_source1,
                       double* new_source, int new_source1,
                       double* flux_temp, int flux_temp1,
                       int ng, int cx, int cy, double tol);

double eigenvalueSolve2d(double **A, int A1, int A2, 
                         double **M, int M1, int M2,
                         double* flux, int flux1,
                         double* old_source, int old_source1,
                         double* new_source, int new_source1,
                         double* flux_temp, int flux_temp1,
                         int ng, int cx, int cy, int cz, double tol);

void linearSolve(double *A, int A1, int A2, 
                 double* flux, int flux1, 
                 double* source, int source1, 
                 double* flux_temp, int flux_temp1, 
                 int cx, int cy, int ng, double tol);

void linearSolve2d(double **A, int A1, int A2, 
                   double* flux, int flux1, 
                   double* source, int source1, 
                   double* flux_temp, int flux_temp1, 
                   int cx, int cy, int cz, int ng, double tol);

void matrix_multiplication(double *matrix, double* vector_x, 
                           double* vector_y, int num_blocks, 
                           int block_width);
void matrix_multiplication2d(double **matrix, double* vector_x, 
                             double* vector_y, int num_blocks, 
                             int block_width);
void vector_scale(double* vector, double scale_value, int length);
void setNumThreads(int num_threads);
void matMultA(double** A, double* flux, double* vec_y, int nx, int ny, int nz, int ng);

#endif /* LINALG_H_ */
