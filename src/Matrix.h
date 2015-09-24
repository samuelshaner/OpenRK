/**
 * @file Matrix.h
 * @brief A matrix object
 * @date May 5, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef MATRIX_H_
#define MATRIX_H_


#ifdef __cplusplus
#include <math.h>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <iomanip>
#include <omp.h>
#include "log.h"
#endif


class Matrix {

private:

  /** A list of lists representing the matrix */
  std::vector< std::map<int, double> > _LIL;

  /** The CSR matrix variables */
  double* _A;
  int* _IA;
  int* _JA;
  double* _DIAG;
  
  bool _modified;
  int _num_x;
  int _num_y;
  int _num_z;
  int _num_groups;
  int _num_rows;

  /** OpenMP mutual exclusion locks for atomic cell updates */
  omp_lock_t* _cell_locks;
  
  void convertToCSR();
  void setNumX(int num_x);
  void setNumY(int num_y);
  void setNumZ(int num_z);
  void setNumGroups(int num_groups);
  
public:
  Matrix(int num_x=1, int num_y=1, int num_z=1, int num_groups=1);
  virtual ~Matrix();

  /* Worker functions */
  void incrementValue(int cell_from, int group_from, int cell_to, int group_to,
                      double val);
  void clear();
  void printString();
  void transpose();

  /* Getter functions */
  double getValue(int cell_from, int group_from, int cell_to,
                        int group_to);
  double* getA();
  int* getIA();
  int* getJA();
  double* getDiag();
  int getNumX();
  int getNumY();
  int getNumZ();
  int getNumGroups();
  int getNumRows();
  int getNNZ();

  /* Setter functions */
  void setValue(int cell_from, int group_from, int cell_to, int group_to,
                double val);
};

#endif /* MATRIX_H_ */
