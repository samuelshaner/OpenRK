/**
 * @file Array.h
 * @brief An array object
 * @date May 5, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef ARRAY_H_
#define ARRAY_H_


#ifdef __cplusplus
#include <math.h>
#include <map>
#include <array>
#include <string>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <omp.h>
#include <cstdlib>
#include <ctime>
#include "pairwise_sum.h"
#include "log.h"
#endif


class Array {

private:

  /** A list of lists representing the array */
  double* _values;
  int _num_dimensions;
  long int* _dimensions;

  void clearObjects();

public:
  Array(long int* dimensions, int num_dimensions);
  virtual ~Array();

  /* Worker functions */
  void clear();
  void scaleByValue(double val);
  void printString();
  void copyTo(Array* array);
  Array* tile(int factor);
  Array* repeat(int factor);
  Array* sumAxis(int axis);
  Array* copy();
  Array* multiply(Array* array, Array* result=NULL);
  Array* divide(Array* array, Array* result=NULL);
  Array* add(Array* array, Array* result=NULL);
  Array* subtract(Array* array, Array* result=NULL);
  void incrementValue(long int index, double value);
  void reshape(long int* dimensions, int num_dimensions);
  Array* flatten();
  void fillWithRandom();
  void outputValues(double* np_array, long int num_values);
  double sum();

  /* Getter functions */
  double getValue(long int* dimensions, int num_dimensions);
  double* getValues();
  int getNumDimensions();
  long int getShape(int dimension);
  long int getSize();
  long int getIndex(long int* dimensions, int num_dimensions);
  void getIndices(long int index, long int* indices);

  /* Setter functions */
  void setDimensions(long int* dimensions, int num_dimensions);
  void setValue(long int index, double value);
  void setValues(double* values, long int num_values);
  void setAll(double val);
};

#endif /* ARRAY_H_ */
