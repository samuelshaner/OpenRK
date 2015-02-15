
/* File: gpu_linalg.cu */

#include "gpu_linalg.h"


/**
 * @brief Solve the linear system Ax=b using Gauss Seidel with SOR.
 * @param pointer to A matrix
 * @param pointer to x vector
 * @param pointer to b vector
 * @param flux convergence criteria
 * @param the maximum number of iterations
 */
void GPUlinearSolve(double *A, double* flux, double* source, double* flux_temp, 
                   int cx, int cy, int cz, int ng, double tol){
  
  /* Initialize variable pointers for GPU */
  double *dev_A, *dev_flux, *dev_source, *dev_flux_temp;
  
  /* Allocate and copy memory to GPU */
  HANDLE_ERROR( cudaMalloc((void**)&dev_A, sizeof(double) * cx*cy*cz*ng*(ng+6)) );
  HANDLE_ERROR( cudaMalloc((void**)&dev_flux, sizeof(double) * cx*cy*cz*ng) );
  HANDLE_ERROR( cudaMalloc((void**)&dev_source, sizeof(double) * cx*cy*cz*ng) );
  HANDLE_ERROR( cudaMalloc((void**)&dev_flux_temp, sizeof(double) * cx*cy*cz*ng) );

  HANDLE_ERROR( cudaMemcpy( dev_A, A, sizeof(double) * cx*cy*cz*ng*(ng+6), cudaMemcpyHostToDevice) );
  HANDLE_ERROR( cudaMemcpy( dev_flux, flux, sizeof(double) * cx*cy*cz*ng, cudaMemcpyHostToDevice) );
  HANDLE_ERROR( cudaMemcpy( dev_source, source, sizeof(double) * cx*cy*cz*ng, cudaMemcpyHostToDevice));

  int iter = 0;
  double _SOR_factor = 1.5;
  double residual;

  while (iter < 1000){

   /* Pass new flux to old flux */
    HANDLE_ERROR(cudaMemcpy( flux_temp, flux, sizeof(double) * cx*cy*cz*ng, cudaMemcpyDeviceToDevice));

    gaussSeidel<<<256, 256>>>(dev_A, dev_flux, dev_source, dev_flux_temp, 
                              cx, cy, cz, ng, _SOR_factor, 0);

    gaussSeidel<<<256, 256>>>(dev_A, dev_flux, dev_source, dev_flux_temp, 
                              cx, cy, cz, ng, _SOR_factor, 1);

    HANDLE_ERROR( cudaMemcpy( &flux_temp, dev_flux_temp, sizeof(double) * cx*cy*cz*ng, cudaMemcpyDeviceToHost) );
    residual = pow(pairwise_sum(flux_temp, cx*cy*cz*ng), 0.5) / (cx*cy*cz*ng);

    iter++;

    log_printf(NORMAL, "GS iter: %i, res: %f", iter, residual);
    
    if (residual < tol && iter > 10)
      break;
  }

  /* Retrieve flux from GPU */
  HANDLE_ERROR( cudaMemcpy( &flux, dev_flux, sizeof(double) * cx*cy*cz*ng, cudaMemcpyDeviceToHost) );

  /* Deallocate memory */
  HANDLE_ERROR( cudaFree(dev_A) );
  HANDLE_ERROR( cudaFree(dev_flux) );
  HANDLE_ERROR( cudaFree(dev_source) );
  HANDLE_ERROR( cudaFree(dev_flux_temp) );
  HANDLE_ERROR( cudaFree(dev_residual) );

  return 0;
}
 

__global__ void gaussSeidel(double *dev_A, double* dev_flux, double* dev_source, double* dev_flux_temp,
                               int cx, int cy, int cz, int ng, double _SOR_factor, int color){

  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  double val;
  int x, y, z;
  int row, cell;

  /* Iteration over red cells */
  while (tid < (cx*cy*cz*ng+(1-color))/2){

    x = (2*(tid+color) % (cx*cy)) % cx;
    y = (2*(tid+color) % (cx*cy)) / cx;
    z = (2*(tid+color) / (cx*cy));
    cell = z*cx*cy+y*cx+x;
    
    for (int g = 0; g < ng; g++){
      
      row = cell*ng + g;
      val = 0.0;
      
      /* Previous flux term */
      val += (1.0 - _SOR_factor) * dev_flux[row];
      
      /* Source term */
      val += _SOR_factor * dev_source[row] / dev_A(cell, g*(ng+6)+g+3);
      
      /* Left surface */
      if (x != 0)
        val -= _SOR_factor * dev_flux[row - ng] *
          dev_A(cell, g*(ng+6)) /
          dev_A(cell, g*(ng+6)+g+3);
      
      /* Back surface */
      if (y != 0)
        val -= _SOR_factor * dev_flux[row - cx * ng] *
          dev_A(cell, g*(ng+6)+1) /
          dev_A(cell, g*(ng+6)+g+3);
      
      /* Bottom surface */
      if (z != 0)
        val -= _SOR_factor * dev_flux[row - cx * cy * ng] *
          dev_A(cell, g*(ng+6)+2) /
          dev_A(cell, g*(ng+6)+g+3);
      
      /* Group-to-group */
      for (int e = 0; e < ng; e++){
        if (e != g)
          val -= _SOR_factor * dev_flux[cell*ng+e] *
            dev_A(cell, g*(ng+6)+3+e) /
            dev_A(cell, g*(ng+6)+g+3);
      }
      
      /* Right surface */
      if (x != cx - 1)
        val -= _SOR_factor * dev_flux[row + ng] *
          dev_A(cell, g*(ng+6)+ng+3) /
          dev_A(cell, g*(ng+6)+g+3);
      
      /* Front surface */
      if (y != cy - 1)
        val -= _SOR_factor * dev_flux[row + ng*cx] *
          dev_A(cell, g*(ng+6)+ng+4) /
          dev_A(cell, g*(ng+6)+g+3);
      
      /* Front surface */
      if (z != cz - 1)
        val -= _SOR_factor * dev_flux[row + ng*cx*cy] *
          dev_A(cell, g*(ng+6)+ng+5) /
          dev_A(cell, g*(ng+6)+g+3);
      
      dev_flux[row] = val;

      /* Store the square residual */
      if (dev_flux[row] != 0.0)
        dev_flux_temp[row] = (dev_flux[row] - dev_flux_temp[row]) / dev_flux[row] * (dev_flux[row] - dev_flux_temp[row]) / dev_flux[row];
      }
    
    tid += blockDim.x * gridDim.x;
  }
}
