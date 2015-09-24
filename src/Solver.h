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
#include <vector>
#include "Geometry.h"
#include "omp.h"
#include "linalg.h"
#include "Matrix.h"
#include "Vector.h"
#include "constants.h"
#endif


class Solver {

protected:

  Matrix* _amp_AM_matrix;
  Matrix* _amp_M_matrix;
  Vector* _amp_source;

  double _k_eff_0;
  Geometry* _geometry;
  Clock* _clock;
  int _method;
  int _num_energy_groups;
  int _num_delayed_groups;
  int _num_shape_cells;
  int _num_amp_cells;
  double _buckling;
  double _initial_power;

  /* Fine mesh field variables */
  std::map<int, Vector*> _temperature;
  std::map<int, Vector*> _flux;
  std::map<int, Vector*> _shape;
  std::map<int, Vector*> _power;
  std::map<int, Vector*> _weight;
  
  /* Coarse mesh field variables */
  std::map<int, Vector*> _amplitude;
  std::map<int, Vector*> _current;
  std::map<int, Vector*> _dif_linear;
  std::map<int, Vector*> _dif_nonlinear;
  std::map<int, Vector*> _frequency;
  
public:
  Solver(Geometry* geometry);
  virtual ~Solver();

  /* Getter functions */
  Matrix* getAmpAMMatrix();
  Matrix* getAmpMMatrix();
  Vector* getAmpSource();
  double getKeff0();
  int getMethod();
  double getBuckling();
  Geometry* getGeometry();
  
  Vector* getTemperature(int state);
  Vector* getFlux(int state);
  Vector* getShape(int state);
  Vector* getAmplitude(int state);
  Vector* getPower(int state);
  Vector* getCurrent(int state);
  Vector* getDifLinear(int state);
  Vector* getDifNonlinear(int state);
  Vector* getFrequency(int state);

  double getTemperatureByValue(int cell, int state);
  double getFluxByValue(int cell, int group, int state);
  double getShapeByValue(int cell, int group, int state);
  double getAmplitudeByValue(int cell, int group, int state);
  double getPowerByValue(int cell, int state);
  double getCurrentByValue(int cell, int group, int side, int state);
  double getDifLinearByValue(int cell, int group, int side, int state);
  double getDifNonlinearByValue(int cell, int group, int side, int state);
  double getFrequencyByValue(int cell, int group, int state);

  /* Setter functions */
  void setMethod(int method);
  void setEndTime(double time);
  void setInnerTimeStepSize(double time);
  void setOuterTimeStepSize(double time);
  void setBuckling(double buckling);
  void setInitialPower(double power);
  
  void setTemperatureByValue(double value, int cell, int state);
  void setShapeByValue(double value, int cell, int group, int state);
  void setAmplitudeByValue(double value, int cell, int group, int state);
  void setPowerByValue(double value, int cell, int state);
  void setFrequencyByValue(double value, int cell, int group, int state);
  void setCurrentByValue(double value, int cell, int group, int side, int state);
  void setFluxByValue(double value, int cell, int group, int state);
  void setDifLinearByValue(double value, int cell, int group, int side, int state);
  void setDifNonlinearByValue(double value, int cell, int group, int side, int state);

  /* Worker functions */
  void generateAmplitudeMatrix();
  void integratePrecursorConcentrations(int state_from, int state_to);
  void integrateTemperature(int state_from, int state_to);
  void interpolateShape(int state, int state_forward, int state_backward);
  void interpolateDifNonlinear(int state, int state_forward, int state_backward);
  void computeDiffusionCoefficients(int state);
  void reconstructFlux(int state, int state_shape, int state_amp);
  void computeShape(int state, int state_flux, int state_amp);
  void computeFrequency();
  void computeInitialPrecursorConcentrations();
  void computePower(int state);
  double computeAveragePower(int state);
  double computePowerRMSError(int state_1, int state_2);
  void normalizeFlux();
  void initializeClock();
  
  virtual void takeInnerStep()=0;
  virtual void takeOuterStep()=0;
  virtual void takeOuterStepOnly()=0;
  virtual void computeInitialShape(double tol)=0;

  /* Copy functions */
  void copyPrecursors(int state_from, int state_to);
  virtual void copyFieldVariables(int state_from, int state_to);
  void broadcastToAll(int state_from);
  
};

#endif /* SOLVER_H_ */
