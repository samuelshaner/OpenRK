
/* File: linalg.cpp */

#include "linalg.h"


double eigenvalueSolve(double *A, int A1, int A2,
                       double *M, int M1, int M2,
                       double* flux, int flux1,
                       double* old_source, int old_source1,
                       double* new_source, int new_source1,
                       double* flux_temp, int flux_temp1,
                       int ng, int cx, int cy, double tol){

  log_printf(NORMAL, "computing eigenvalue...");

  /* Initialize variables */
  double sum_new, sum_old, val, residual, scale_val, _k_eff;
  int row;

  /* Compute and normalize the initial source */
  matrix_multiplication(M, flux, old_source, cx*cy, ng);
  sum_old = pairwise_sum(old_source, cx*cy*ng);
  scale_val = (cx*cy*ng) / sum_old;
  vector_scale(old_source, scale_val, cx*cy*ng);
  vector_scale(flux, scale_val, cx*cy*ng);
  sum_old = cx*cy*ng;

  /* Power iteration diffusion solver */
  for (int iter = 0; iter < 25000; iter++){

    /* Solve phi = A^-1 * old_source */
    linearSolve(A, A1, A2, flux, flux1, old_source, old_source1,
                flux_temp, flux_temp1, cx, cy, ng, tol);

    /* Compute the new source */
    matrix_multiplication(M, flux, new_source, cx*cy, ng);
    sum_new = pairwise_sum(new_source, cx*cy*ng);

    /* Compute and set keff */
    _k_eff = sum_new / sum_old;

    /* Scale the old source by keff */
    vector_scale(old_source, _k_eff, cx*cy*ng);

    /* Compute the L2 norm of source error */
    residual = 0.0;
    for (int i = 0; i < cx*cy*ng; i++){
      if (new_source[i] != 0.0)
        residual += pow((new_source[i] - old_source[i]) / new_source[i], 2);
    }

    /* Compute the source RMS error */
    residual = sqrt(residual / (cx*cy*ng));

    /* Normalize the new source to have an average value of 1.0 */
    scale_val = (cx*cy*ng) / sum_new;
    vector_scale(new_source, scale_val, cx*cy*ng);
    std::copy(new_source, new_source + cx*cy*ng, old_source);

    log_printf(NORMAL, "CMFD iter: %i, keff: %f, error: %f",
               iter, _k_eff, residual);

    /* Check for convergence */
    if (residual < tol && iter > 10)
      break;
  }

  return _k_eff;
}


double eigenvalueSolve2d(double **A, int A1, int A2,
                         double **M, int M1, int M2,
                         double* flux, int flux1,
                         double* old_source, int old_source1,
                         double* new_source, int new_source1,
                         double* flux_temp, int flux_temp1,
                         int ng, int cx, int cy, int cz, double tol){

  log_printf(NORMAL, "computing eigenvalue...");

  /* Initialize variables */
  double sum_new, sum_old, val, residual, scale_val, _k_eff;
  int row;

  /* Compute and normalize the initial source */
  matrix_multiplication2d(M, flux, old_source, cx*cy*cz, ng);
  sum_old = pairwise_sum(old_source, cx*cy*cz*ng);
  scale_val = (cx*cy*cz*ng) / sum_old;
  vector_scale(old_source, scale_val, cx*cy*cz*ng);
  vector_scale(flux, scale_val, cx*cy*cz*ng);
  sum_old = cx*cy*cz*ng;

  /* Power iteration diffusion solver */
  for (int iter = 0; iter < 25000; iter++){

    /* Solve phi = A^-1 * old_source */
    linearSolve2d(A, A1, A2, flux, flux1, old_source, old_source1,
                  flux_temp, flux_temp1, cx, cy, cz, ng, tol);

    /* Compute the new source */
    matrix_multiplication2d(M, flux, new_source, cx*cy*cz, ng);
    sum_new = pairwise_sum(new_source, cx*cy*cz*ng);

    /* Compute and set keff */
    _k_eff = sum_new / sum_old;

    /* Scale the old source by keff */
    vector_scale(old_source, _k_eff, cx*cy*cz*ng);

    /* Compute the L2 norm of source error */
    #pragma omp parallel for
    for (int i = 0; i < cx*cy*cz*ng; i++){
      if (new_source[i] != 0.0)
        old_source[i] = pow((new_source[i] - old_source[i]) / new_source[i], 2);
    }

    /* Compute the source RMS error */
    residual = sqrt(pairwise_sum(old_source, cx*cy*cz*ng) / (cx*cy*cz*ng));

    /* Normalize the new source to have an average value of 1.0 */
    scale_val = (cx*cy*cz*ng) / sum_new;
    vector_scale(new_source, scale_val, cx*cy*cz*ng);
    std::copy(new_source, new_source + cx*cy*cz*ng, old_source);

    log_printf(NORMAL, "CMFD iter: %i, keff: %f, error: %f",
               iter, _k_eff, residual);

    /* Check for convergence */
    if (residual < tol && iter > 10)
      break;
  }

  return _k_eff;
}


/**
 * @brief Solve the linear system Ax=b using Gauss Seidel with SOR.
 * @param pointer to A matrix
 * @param pointer to x vector
 * @param pointer to b vector
 * @param flux convergence criteria
 * @param the maximum number of iterations
 */
void linearSolve(double *A, int A1, int A2, 
                 double* flux, int flux1, 
                 double* source, int source1, 
                 double* flux_temp, int flux_temp1, 
                 int cx, int cy, int ng, double tol){


  double residual = 1E10;
  int row, cell;
  double val;
  int iter = 0;
  double _SOR_factor = 1.5;

  while (iter < 1000){

    /* Pass new flux to old flux */
    std::copy(flux, flux + cx*cy*ng, flux_temp);

    /* Iteration over red cells */
    #pragma omp parallel for private(row, val, cell)
    for (int y = 0; y < cy; y++){
      for (int x = y % 2; x < cx; x += 2){

        cell = y*cx+x;

        for (int g = 0; g < ng; g++){

          row = cell*ng + g;
          val = 0.0;

          /* Previous flux term */
          val += (1.0 - _SOR_factor) * flux[row];

          /* Source term */
          val += _SOR_factor * source[row] / A[cell*ng*(ng+4) + g*(ng+4)+g+2];

          /* Left surface */
          if (x != 0)
            val -= _SOR_factor * flux[row - ng] *
                A[cell*ng*(ng+4) + g*(ng+4)] /
                A[cell*ng*(ng+4) + g*(ng+4)+g+2];

          /* Bottom surface */
          if (y != 0)
            val -= _SOR_factor * flux[row - cx * ng] *
                A[cell*ng*(ng+4) + g*(ng+4)+1] /
                A[cell*ng*(ng+4) + g*(ng+4)+g+2];

          /* Group-to-group */
          for (int e = 0; e < ng; e++){
            if (e != g)
              val -= _SOR_factor * flux[cell*ng+e] *
                  A[cell*ng*(ng+4) + g*(ng+4)+2+e] /
                  A[cell*ng*(ng+4) + g*(ng+4)+g+2];
          }

          /* Right surface */
          if (x != cx - 1)
            val -= _SOR_factor * flux[row + ng] *
                A[cell*ng*(ng+4) + g*(ng+4)+ng+2] /
                A[cell*ng*(ng+4) + g*(ng+4)+g+2];

          /* Top surface */
          if (y != cy - 1)
            val -= _SOR_factor * flux[row + ng*cx] *
                A[cell*ng*(ng+4) + g*(ng+4)+ng+3] /
                A[cell*ng*(ng+4) + g*(ng+4)+g+2];

          flux[row] = val;
        }
      }
    }

    /* Iteration over black cells */
    #pragma omp parallel for private(row, val, cell)
    for (int y = 0; y < cy; y++){
      for (int x = 1 - y % 2; x < cx; x += 2){

        cell = y*cx+x;

        for (int g = 0; g < ng; g++){

          row = cell*ng + g;
          val = 0.0;

          /* Previous flux term */
          val += (1.0 - _SOR_factor) * flux[row];

          /* Source term */
          val += _SOR_factor*source[row] / A[cell*ng*(ng+4) + g*(ng+4)+g+2];

          /* Left surface */
          if (x != 0)
            val -= _SOR_factor * flux[row - ng] *
                A[cell*ng*(ng+4) + g*(ng+4)] /
                A[cell*ng*(ng+4) + g*(ng+4)+g+2];

          /* Bottom surface */
          if (y != 0)
            val -= _SOR_factor * flux[row - cx * ng] *
                A[cell*ng*(ng+4) + g*(ng+4)+1] /
                A[cell*ng*(ng+4) + g*(ng+4)+g+2];

          /* Group-to-group */
          for (int e = 0; e < ng; e++){
            if (e != g)
              val -= _SOR_factor * flux[cell*ng+e] *
                  A[cell*ng*(ng+4) + g*(ng+4)+2+e] /
                  A[cell*ng*(ng+4) + g*(ng+4)+g+2];
          }

          /* Right surface */
          if (x != cx - 1)
            val -= _SOR_factor * flux[row + ng] *
                A[cell*ng*(ng+4) + g*(ng+4)+ng+2] /
                A[cell*ng*(ng+4) + g*(ng+4)+g+2];

          /* Top surface */
          if (y != cy - 1)
            val -= _SOR_factor * flux[row + ng*cx] *
                A[cell*ng*(ng+4) + g*(ng+4)+ng+3] /
                A[cell*ng*(ng+4) + g*(ng+4)+g+2];

          flux[row] = val;
        }
      }
    }

    /* Compute the average residual */
    residual = 0.0;
    for (int i = 0; i < cx*cy*ng; i++){
      if (flux[i] != 0.0)
        residual += pow((flux[i] - flux_temp[i]) / flux[i], 2);
    }

    residual = pow(residual, 0.5) / (cx*cy*ng);

    /* Increment the interations counter */
    iter++;

    log_printf(DEBUG, "GS iter: %i, res: %f", iter, residual);

    if (residual < tol && iter > 10)
      break;
  }

  log_printf(DEBUG, "linear solver iterations: %i", iter);
}


/**
 * @brief Solve the linear system Ax=b using Gauss Seidel with SOR.
 * @param pointer to A matrix
 * @param pointer to x vector
 * @param pointer to b vector
 * @param flux convergence criteria
 * @param the maximum number of iterations
 */
void linearSolve2d(double **A, int A1, int A2, 
                   double* flux, int flux1, 
                   double* source, int source1, 
                   double* flux_temp, int flux_temp1, 
                   int cx, int cy, int cz, int ng, double tol){
  

  double residual = 1E10;
  int row, cell;
  double val;
  int iter = 0;
  double _SOR_factor = 1.5;
  
  while (iter < 1000){

    /* Pass new flux to old flux */
    std::copy(flux, flux + cx*cy*cz*ng, flux_temp);

    /* Iteration over red cells */
    for (int z = 0; z < cz; z++){
      #pragma omp parallel for private(row, val, cell)
      for (int y = 0; y < cy; y++){
        for (int x = (z+y) % 2; x < cx; x += 2){
          
          cell = z*cx*cy+y*cx+x;

          for (int g = 0; g < ng; g++){
            
            row = cell*ng + g;
            val = 0.0;
            
            /* Previous flux term */
            val += (1.0 - _SOR_factor) * flux[row];
            
            /* Source term */
            val += _SOR_factor * source[row] / A[cell][g*(ng+6)+g+3];
            
            /* Left surface */
            if (x != 0)
              val -= _SOR_factor * flux[row - ng] *
                A[cell][g*(ng+6)] /
                A[cell][g*(ng+6)+g+3];
            
            /* Back surface */
            if (y != 0)
              val -= _SOR_factor * flux[row - cx * ng] *
                A[cell][g*(ng+6)+1] /
                A[cell][g*(ng+6)+g+3];

            /* Bottom surface */
            if (z != 0)
              val -= _SOR_factor * flux[row - cx * cy * ng] *
                A[cell][g*(ng+6)+2] /
                A[cell][g*(ng+6)+g+3];
            
            /* Group-to-group */
            for (int e = 0; e < ng; e++){
              if (e != g)
                val -= _SOR_factor * flux[cell*ng+e] *
                  A[cell][g*(ng+6)+3+e] /
                  A[cell][g*(ng+6)+g+3];
            }
            
            /* Right surface */
            if (x != cx - 1)
              val -= _SOR_factor * flux[row + ng] *
                A[cell][g*(ng+6)+ng+3] /
                A[cell][g*(ng+6)+g+3];
            
            /* Front surface */
            if (y != cy - 1)
              val -= _SOR_factor * flux[row + ng*cx] *
                A[cell][g*(ng+6)+ng+4] /
                A[cell][g*(ng+6)+g+3];

            /* Front surface */
            if (z != cz - 1)
              val -= _SOR_factor * flux[row + ng*cx*cy] *
                A[cell][g*(ng+6)+ng+5] /
                A[cell][g*(ng+6)+g+3];
            
            flux[row] = val;
          }
        }
      }
    }

    /* Iteration over black cells */
    for (int z = 0; z < cz; z++){
      #pragma omp parallel for private(row, val, cell)
      for (int y = 0; y < cy; y++){
        for (int x = 1 - (z+y) % 2; x < cx; x += 2){
          
          cell = z*cx*cy+y*cx+x;

          for (int g = 0; g < ng; g++){
            
            row = cell*ng + g;
            val = 0.0;
            
            /* Previous flux term */
            val += (1.0 - _SOR_factor) * flux[row];
            
            /* Source term */
            val += _SOR_factor * source[row] / A[cell][g*(ng+6)+g+3];
            
            /* Left surface */
            if (x != 0)
              val -= _SOR_factor * flux[row - ng] *
                A[cell][g*(ng+6)] /
                A[cell][g*(ng+6)+g+3];
            
            /* Back surface */
            if (y != 0)
              val -= _SOR_factor * flux[row - cx * ng] *
                A[cell][g*(ng+6)+1] /
                A[cell][g*(ng+6)+g+3];

            /* Bottom surface */
            if (z != 0)
              val -= _SOR_factor * flux[row - cx * cy * ng] *
                A[cell][g*(ng+6)+2] /
                A[cell][g*(ng+6)+g+3];
            
            /* Group-to-group */
            for (int e = 0; e < ng; e++){
              if (e != g)
                val -= _SOR_factor * flux[cell*ng+e] *
                  A[cell][g*(ng+6)+3+e] /
                  A[cell][g*(ng+6)+g+3];
            }
            
            /* Right surface */
            if (x != cx - 1)
              val -= _SOR_factor * flux[row + ng] *
                A[cell][g*(ng+6)+ng+3] /
                A[cell][g*(ng+6)+g+3];
            
            /* Front surface */
            if (y != cy - 1)
              val -= _SOR_factor * flux[row + ng*cx] *
                A[cell][g*(ng+6)+ng+4] /
                A[cell][g*(ng+6)+g+3];

            /* Front surface */
            if (z != cz - 1)
              val -= _SOR_factor * flux[row + ng*cx*cy] *
                A[cell][g*(ng+6)+ng+5] /
                A[cell][g*(ng+6)+g+3];
            
            flux[row] = val;
          }
        }
      }
    }

    /* Compute the average residual */
    residual = 0.0;
    #pragma omp parallel for
    for (int i = 0; i < cx*cy*cz*ng; i++){
      if (flux[i] != 0.0)
        flux_temp[i] = pow((flux[i] - flux_temp[i]) / flux[i], 2);
    }

    residual = pow(pairwise_sum(flux_temp, cx*cy*cz*ng), 0.5) / (cx*cy*cz*ng);

    /* Increment the interations counter */
    iter++;

    if (tol == 1.e-9)
      log_printf(DEBUG, "GS iter: %i, res: %f", iter, residual);

    if (residual < tol && iter > 10)
      break;
  }

  log_printf(DEBUG, "linear solver iterations: %i", iter);
}


/**
 * @brief Multiply matrix by vector (i.e., y = M *x).
 * @param matrix source matrix
 * @param vector_x x vector
 * @param vector_y y vector
 * @param num_blocks number of cell blocks in M matrix.
 * @param block_width number of elements in cell blocks in M matrix.
 */
void matrix_multiplication(double *matrix, double* vector_x, 
                           double* vector_y, int num_blocks, 
                           int block_width){

  memset(vector_y, 0.0, sizeof(double) * num_blocks*block_width);
  
  #pragma omp parallel for
  for (int i = 0; i < num_blocks; i++){
    for (int g = 0; g < block_width; g++){
      for (int e = 0; e < block_width; e++){
        vector_y[i*block_width+g] += matrix[i*(block_width*block_width) + g*block_width+e] 
            * vector_x[i*block_width+e];
      }
    }
  }
}


/**
 * @brief Multiply matrix by vector (i.e., y = M *x).
 * @param matrix source matrix
 * @param vector_x x vector
 * @param vector_y y vector
 * @param num_blocks number of cell blocks in M matrix.
 * @param block_width number of elements in cell blocks in M matrix.
 */
void matrix_multiplication2d(double **matrix, double* vector_x, 
                           double* vector_y, int num_blocks, 
                           int block_width){

  memset(vector_y, 0.0, sizeof(double) * num_blocks*block_width);

  #pragma omp parallel for
  for (int i = 0; i < num_blocks; i++){
    for (int g = 0; g < block_width; g++){
      for (int e = 0; e < block_width; e++){
        vector_y[i*block_width+g] += matrix[i][g*block_width+e] 
            * vector_x[i*block_width+e];
      }
    }
  }
}


/**
 * @brief Scale vectgor by some scalar value.
 * @param vector vector to be scaled
 * @param scale_value value to scale vector
 * @param length vector length
 */
void vector_scale(double* vector, double scale_value, int length){

  #pragma omp parallel for
  for (int i = 0; i < length; i++)
    vector[i] *= scale_value;
}


/**                                                                                                                                                                                                                                        
 * @brief Sets the number of shared memory OpenMP threads to use (>0).                                                                                                                                                                     
 * @param num_threads the number of threads                                                                                                                                                                                                
 */
void setNumThreads(int num_threads) {

  if (num_threads <= 0)
      log_printf(ERROR, "Unable to set the number of threads for the Solver "
                 "to %d since it is less than or equal to 0", num_threads);

  /* Set the number of threads for OpenMP */
  omp_set_num_threads(num_threads);
}


void matMultA(double** A, double* flux, double* vec_y, int nx, int ny, int nz, int ng){

  vector_scale(vec_y, 0.0, nx*ny*nz*ng);
  int row;
  int cell;
  double val;
  
  for (int z = 0; z < nz; z++){
    for (int y = 0; y < ny; y++){
      for (int x = 0; x < nx; x++){
        for (int g = 0; g < ng; g++){
          row = (y*nx+x)*ng + g;
          
          cell = z*nx*ny+y*nx+x;
          
          for (int g = 0; g < ng; g++){
            
            row = cell*ng + g;
            val = 0.0;

            /* Left surface */
            if (x != 0)
              val += flux[row - ng] * A[cell][g*(ng+6)];
            
            /* Back surface */
            if (y != 0)
              val += flux[row - nx * ng] * A[cell][g*(ng+6)+1];
            
            /* Bottom surface */
            if (z != 0)
              val += flux[row - nx * ny * ng] * A[cell][g*(ng+6)+2];
            
            /* Group-to-group */
            for (int e = 0; e < ng; e++)
              val += flux[cell*ng+e] * A[cell][g*(ng+6)+3+e];
            
            /* Right surface */
            if (x != nx - 1)
              val += flux[row + ng] * A[cell][g*(ng+6)+ng+3];
            
            /* Front surface */
            if (y != ny - 1)
              val += flux[row + ng*nx] * A[cell][g*(ng+6)+ng+4];
            
            /* Front surface */
            if (z != nz - 1)
              val += flux[row + ng*nx*ny] * A[cell][g*(ng+6)+ng+5];
            
            vec_y[row] = -val;
          }
        }
      }
    } 
  }
}
