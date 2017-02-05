/* File: openrk.i */
%module openrk

%{
  #define SWIG_FILE_WITH_INIT
  #include "../src/log.h"
  #include "../src/Array.h"
  #include "../src/Matrix.h"
  #include "../src/linalg.h"

  #define printf PySys_WriteStdout

  /* Exception helpers */
  static int swig_c_error_num = 0;
  static char swig_c_err_msg[1024];

  const char* err_occurred(void) {
    if (swig_c_error_num) {
      swig_c_error_num = 0;
      return (const char*)swig_c_err_msg;
    }
    return NULL;
  }

  void set_err(const char *msg) {
    swig_c_error_num = 1;
    strncpy(swig_c_err_msg, msg, 1024);
  }

%}

%warnfilter(506) log_printf(logLevel level, const char *format, ...);
%warnfilter(511) swig::SwigPyIterator;


%exception {
  try {
    $function
  } catch (const std::exception &e) {
      SWIG_exception(SWIG_RuntimeError, e.what());
  }
}


%include "numpy.i"

%init %{
  import_array();
%}


%apply ( double* INPLACE_ARRAY1, int DIM1 ) {(double* seq, int n), (double* flux, int flux1),
         (double* old_source, int old_source1), (double* new_source, int new_source1),
         (double* flux_temp, int flux_temp1), (double* source, int source1)};
%apply ( double* INPLACE_ARRAY2, int DIM1, int DIM2 ) {(double *A, int A1, int A2),
         (double *M, int M1, int M2)};

%apply (double* IN_ARRAY1, int DIM1) {(double* xs, int num_groups),
     (double* time_steps, int num_steps)};


%apply (long int* IN_ARRAY1, int DIM1) {(long int* dimensions, int num_dimensions)};
%apply (double* IN_ARRAY1, int DIM1) {(double* values, long int num_values)};
%apply (int* IN_ARRAY1, int DIM1) {(int* diags, int num_diags)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* np_array, long int num_values)};



%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double* xs, int num_time_steps, int num_groups),
     (double* x, int num_time_steps, int num_groups_squared)};

%include <exception.i>
%include ../src/log.h
%include ../src/Array.h
%include ../src/Matrix.h
%include ../src/linalg.h

#define printf PySys_WriteStdout
