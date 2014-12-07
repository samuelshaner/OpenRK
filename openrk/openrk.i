
/* File: openrk.i */
%module openrk

%{
  #define SWIG_FILE_WITH_INIT
  #include "../src/linalg.h"
  #include "../src/log.h"


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
%include <exception.i>
%include ../src/linalg.h
%include ../src/log.h

#define printf PySys_WriteStdout
