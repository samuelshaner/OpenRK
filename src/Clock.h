/**
 * @file Clock.h
 * @brief A clock object to keep track of time.
 * @date February 7, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef CLOCK_H_
#define CLOCK_H_


#ifdef __cplusplus
#include <math.h>
#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <iomanip>
#include "log.h"
#endif


/**
 * @enum clockPosition
 * @brief The positions of the clock
 */
enum clockPosition {
  START,
  PREVIOUS_OUT,
  PREVIOUS_IN,
  CURRENT,
  FORWARD_IN_OLD,
  FORWARD_OUT,
  FORWARD_OUT_OLD,
  END
};

class Clock {

private:

  /** A map of time positions to time values */
  std::map<int, double> _times;

  /** The inner time step size */
  double _dt_inner;

  /* The outer time step size */
  double _dt_outer;

public:
  Clock(double start_time=0.0, double end_time=3.0, double dt_outer=1.e-1, double dt_inner=1.e-2);
  virtual ~Clock();

  /* Getter functions */
  double getTime(int position);
  double getDtInner();
  double getDtOuter();
  std::string getPositionName(int position);
  
  /* Setter functions */
  void setDtOuter(double dt_outer);
  void setDtInner(double dt_inner);
  void setTime(int position, double time);
  void setStartTime(double time);
  void setEndTime(double time);
  
  /* Worker functions */
  void takeInnerStep();
  void takeOuterStep();
  void resetToPreviousOuterStep();
  std::string toString();
  void printString();
};

#endif /* CLOCK_H_ */
