
/* File: linalg.cpp */

#include "linalg.h"

int fact(int n){

  if (n < 0)
    return 0;
  if (n == 0)
    return 1;
  else{
    return n * fact(n-1); 
  }
}


double rms(double* seq, int n){

  double avg = 0.0;

  for (int i = 0; i < n; i++){
    avg += seq[i];
  }

  avg = avg/n;
  
  double rms = 0.0;

  for (int i = 0; i < n; i++){
    rms += (seq[i] - avg) * (seq[i] - avg);
  }
  
  rms = sqrt(rms);

  seq[0] = 5.0;

  return rms;
}


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
    vector_copy(new_source, cx*cy*ng, old_source, cx*cy*ng);
    
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
    vector_copy(flux, cx*cy*ng, flux_temp, cx*cy*ng);

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
 * @brief Copy a vector to another vector.
 * @param vector_from vector to be copied
 * @param vector_to vector to receive copied data
 */
void vector_copy(double* vector_from, int length_from, double* vector_to, int length_to){

  for (int i = 0; i < length_to; i++)
    vector_to[i] = vector_from[i];
}


/**
 * @brief Assign all elements in a matrix to zero.
 * @param mat matrix to be zeroed
 * @param width width of matrix row
 * @param height height of matrix copy
 */
void matrix_zero(double *matrix, int width, int length){

  for (int i = 0; i < length; i++){
    for (int g = 0; g < width; g++)
      matrix[i*width + g] = 0.0;
  }
}


/**
 * @brief Assign all elements in a matrix to zero.
 * @param vecotr vector to be zeroed
 * @param length length of vector
 */
void vector_zero(double* vector, int length){

  for (int i = 0; i < length; i++)
    vector[i] = 0.0;
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

  vector_zero(vector_y, num_blocks*block_width); 

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
 * @brief Scale vectgor by some scalar value.
 * @param vector vector to be scaled
 * @param scale_value value to scale vector
 * @param length vector length
 */
void vector_scale(double* vector, double scale_value, int length){

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
