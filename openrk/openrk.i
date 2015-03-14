/* File: openrk.i */
%module openrk

%{
  #define SWIG_FILE_WITH_INIT
  #include "../src/linalg.h"
  #include "../src/log.h"
  #include "../src/Material.h"
  #include "../src/FunctionalMaterial.h"
  #include "../src/Clock.h"
  #include "../src/Mesh.h"
  #include "../src/ShapeMesh.h"
  #include "../src/StructuredShapeMesh.h"
  #include "../src/UnstructuredShapeMesh.h"
  #include "../src/AmpMesh.h"
  #include "../src/Solver.h"
  #include "../src/Transient.h"
  
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

%apply (double* IN_ARRAY1, int DIM1) {(double* xs, int num_groups),
     (double* time_steps, int num_steps), (double* current, int num_cells_times_groups),
     (double* flux, int num_cells_times_groups};

%apply (int* IN_ARRAY1, int DIM1) {(int* group_indices, int length_group_indices)};

%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double* xs, int num_time_steps, int num_groups),
     (double* x, int num_time_steps, int num_groups_squared)};

/* Typemap for Lattice::setUniverses(int num_x, int num_y, Universe** universes)
 * method - allows users to pass in Python lists of Universes for each
 * lattice cell */
%include <std_map.i>
%typemap(in) (int** regions_to_cells, int* regions_per_cell, int num_cells){

  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError,"Expected a Python list of lists of integers"
                    "for regions to cells map");
    return NULL;
  }
  
  $3 = PySequence_Length($input);  // num_cells
  $2 = (int*) malloc(($3) * sizeof(int)); // regions_per_cell
  $1 = (int**) malloc(($3) * sizeof(int*)); // regions_to_cells

  /* Loop over cells */
  for (int i = 0; i < $3; i++) {

    /* Get the list of regions in this cell */
    PyObject* regions_in_cell = PyList_GetItem($input,i);
    $2[i]= PySequence_Length(regions_in_cell); // regions in cell i
    $3[i] = (int*) malloc(($2[i]) * sizeof(int)); // allocate regions in this cell

    for (int r=0; r < $2[i]; r++)
      $3[i][r] = PyList_GetItem(regions_in_cell, r);
  }
}


%include <exception.i>
%include ../src/linalg.h
%include ../src/log.h
%include ../src/Material.h
%include ../src/FunctionalMaterial.h
%include ../src/Clock.h
%include ../src/Mesh.h
%include ../src/ShapeMesh.h
%include ../src/StructuredShapeMesh.h
%include ../src/UnstructuredShapeMesh.h
%include ../src/AmpMesh.h
%include ../src/Solver.h
%include ../src/Transient.h

#define printf PySys_WriteStdout
