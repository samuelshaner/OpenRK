
/* File: linalg.h */

#include<math.h>
#include "log.h"
#include "pairwise_sum.h"
#include<omp.h>

int fact(int n);
double rms(double* seq, int n);

double eigenvalueSolve(double *A, int A1, int A2, 
                       double *M, int M1, int M2,
                       double* flux, int flux1,
                       double* old_source, int old_source1,
                       double* new_source, int new_source1,
                       double* flux_temp, int flux_temp1,
                       int ng, int cx, int cy, double tol);

void linearSolve(double *A, int A1, int A2, 
                 double* flux, int flux1, 
                 double* source, int source1, 
                 double* flux_temp, int flux_temp1, 
                 int cx, int cy, int ng, double tol);

void vector_copy(double* vector_from, int length_from,
                 double* vector_to, int length_to);
void matrix_zero(double *matrix, int width, int length);
void vector_zero(double* vector, int length);
void matrix_multiplication(double *matrix, double* vector_x, 
                           double* vector_y, int num_blocks, 
                           int block_width);
void vector_scale(double* vector, double scale_value, int length);
void setNumThreads(int num_threads);
