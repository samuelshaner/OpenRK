/**
 * @file Solver.h
 * @brief A solver object to solve for the flux.
 * @date February 8, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#ifdef __cplusplus
#include <sstream>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <map>
#include "Geometry.h"
#include "omp.h"
#include "linalg.h"
#endif

/**
 * @enum boundaryType
 * @brief The boundary types
 */
enum transientMethod {
  IQS,
  THETA
};


class Solver {

protected:

  double** _amp_matrix;
  double* _amp_source;

  double _k_eff_0;
  Geometry* _geometry;
  Clock* _clock;
  transientMethod _method;
  int _num_energy_groups;
  int _num_delayed_groups;
  int _num_shape_cells;
  int _num_amp_cells;
  double _buckling;
  double _initial_power;

  /* Fine mesh field variables */
  std::vector<int, double*> _temperature;
  std::vector<int, double*> _flux;
  std::vector<int, double*> _shape;
  std::vector<int, double*> _power;

  /* Coarse mesh field variables */
  std::map<int, double*> _amplitude;
  std::map<int, double*> _current;
  std::map<int, double*> _dif_linear;
  std::map<int, double*> _dif_nonlinear;
  std::map<int, double*> _frequency;
  
public:
  Solver(Geometry* geometry);
  virtual ~Solver();

  /* Getter functions */
  double** getAmpMatrix();
  double** getAmpSource();
  double getKeff0();
  transientMethod getMethod();
  double getBuckling();

  double* getTemperature(int time);
  double* getFlux(int time);
  double* getShape(int time);
  double* getAmplitude(int time);
  double* getPower(int time);
  double* getCurrent(int time);
  double* getDifLinear(int time);
  double* getDifNonlinear(int time);
  double* getDifFrequency(int time);

  double getTemperatureByValue(int cell, int time);
  double getFluxByValue(int cell, int group, int time);
  double getShapeByValue(int cell, int group, int time);
  double getAmplitudeByValue(int cell, int group, int time);
  double getPowerByValue(int cell, int time);
  double getCurrentByValue(int cell, int group, int side, int time);
  double getDifLinearByValue(int cell, int group, int side, int time);
  double getDifNonlinearByValue(int cell, int group, int side, int time);
  double getFrequencyByValue(int cell, int group, int time);


  /* Setter functions */
  void setMethod(transientMethod method);  
  void setEndTime(double time);
  void setInnerTimeStepSize(double time);
  void setOuterTimeStepSize(double time);
  void setBuckling(double buckling);
  void setInitialPower(double power);

  void setTemperatureByValue(int value, int cell, int time);
  void setShapeByValue(int value, int cell, int group, int time);
  void setAmplitudeByValue(int value, int cell, int group, int time);
  void setPowerByValue(int value, int cell, int time);
  void setFrequencyByValue(int value, int cell, int group, int time);
  void setCurrentByValue(double value, int cell, int group, int side, int time);
  void setFluxByValue(double value, int cell, int group, int time);
  void setDifLinearByValue(double value, int cell, int group, int time);
  void setDifNonlinearByValue(double value, int cell, int group, int time);

  /* Worker functions */
  void generateAmplitudeMatrix(double wt);
  void integratePrecursorConcentrations(int time_from, int time_to);
  void integrateTemperature(int time_from, int time_to);
  void interpolateShape(int time, int time_forward, int time_backward);
  void interpolateDifNonlinear(int time, int time_forward, int time_backward);
  void computeDiffusionCoefficients(int time);
  void reconstructFlux(int time, int time_shape, int time_amp);
  void computeShape(int time, int time_flux, int time_amp);
  void computeFrequency();
  void computeInitialPrecursorConcentrations();
  void computePower(int time);
  void computePowerRMSError(int time_1, int time_2);
  void normalizeFlux();

  virtual void takeInnerStep()=0;
  virtual void takeOuterStep()=0;
  virtual void takeOuterStepOnly()=0;
  virtual void computeInitialShape()=0;

  /* Copy functions */
  void copyTemperature(int time_from, int time_to);
  void copyFlux(int time_from, int time_to);
  void copyShape(int time_from, int time_to);
  void copyAmplitude(int time_from, int time_to);
  void copyPower(int time_from, int time_to);
  void copyCurrent(int time_from, int time_to);
  void copyDifLinear(int time_from, int time_to);
  void copyDifNonlinear(int time_from, int time_to);
  void copyFrequency(int time_from, int time_to);
  void copyPrecursors(int time_from, int time_to);
  virtual void copyFieldVariables(int time_from, int time_to);
  void broadcastToAll(int time_from);
  
};

#endif /* SOLVER_H_ */
