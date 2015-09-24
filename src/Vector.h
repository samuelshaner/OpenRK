/**
 * @file Vector.h
 * @brief A vector object
 * @date May 5, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef VECTOR_H_
#define VECTOR_H_


#ifdef __cplusplus
#include <math.h>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <omp.h>
#include "log.h"
#include "pairwise_sum.h"
#endif


class Vector {

private:

  /** A list of lists representing the vector */
  double* _array;
  int _num_rows;
  int _num_x;
  int _num_y;
  int _num_z;
  int _num_groups;

  /** OpenMP mutual exclusion locks for atomic cell updates */
  omp_lock_t* _cell_locks;

  void setNumX(int num_x);
  void setNumY(int num_y);
  void setNumZ(int num_z);
  void setNumGroups(int num_groups);
  
public:
  Vector(int num_x=1, int num_y=1, int num_z=1, int num_groups=1);
  virtual ~Vector();

  /* Worker functions */
  void incrementValue(int cell, int group, double val);
  void incrementValues(int cell, int group_start, int group_end, double* vals);
  void clear();
  void scaleByValue(double val);  
  void printString();
  void copyTo(Vector* vector);
  
  /* Getter functions */
  double getValue(int cell, int group);
  double* getArray();
  int getNumX();
  int getNumY();
  int getNumZ();
  int getNumGroups();
  int getNumRows();
  double getSum();
  
  /* Setter functions */
  void setValue(int cell, int group, double val);
  void setValues(int cell, int group_start, int group_end, double* vals);
  void setAll(double val);
};

#endif /* VECTOR_H_ */
