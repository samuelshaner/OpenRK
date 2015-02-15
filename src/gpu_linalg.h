
/* File: gpu_linalg.h */

#ifndef GPU_LINALG_H_
#define GPU_LINALG_H_

#ifdef __cplusplus
#include<math.h>
#include "pairwise_sum.h"
#endif

/** Indexing macro for the scalar flux in each FSR and energy group */
#define dev_A(r,e) (dev_A[(r)*ng*(ng+6) + (e)])


void GPUlinearSolve(double *A, double* flux, double* source, double* flux_temp,
                   int cx, int cy, int cz, int ng, double tol);

#endif /* GPU_LINALG_H_ */
