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
#include "Array.h"
#endif


class Matrix {

private:

  /** A list of lists representing the matrix */
  std::vector< std::map<long int, double> > _LIL;

  /** The CSR matrix variables */
  double* _A;
  long int* _IA;
  long int* _JA;
  double* _DIAG;

  int* _diags;
  int _num_diags;
  bool _modified;
  long int _num_cells;

  void convertToCSR();
  void setNumCells(long int num_cells);

public:
  Matrix(long int num_cells=1);
  virtual ~Matrix();

  /* Worker functions */
  void incrementValue(long int col, long int row, double value);
  void clear();
  void printString();
  void transpose();
  void diags(Array* array);
  void blockDiags(Array* array, int block_size);
  void fillWithRandom();

  /* Getter functions */
  double getValue(long int col, long int row);
  double* getA();
  long int* getIA();
  long int* getJA();
  double* getDiag();
  long int getNumCells();
  long int getNNZ();

  /* Setter functions */
  void setValue(long int col, long int row, double value);
  void setDiags(int* diags, int num_diags);
};

#endif /* MATRIX_H_ */
