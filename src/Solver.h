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

  double _k_eff_0;
  Geometry* _geometry;
  Clock* _clock;
  int _num_energy_groups;
  int _num_delayed_groups;
  int _num_cells;
  double _buckling;
  double _initial_power;
  double _flux_solve_tolerance;
  
  /* Fine mesh field variables */
  std::map<int, Vector*> _temperature;
  std::map<int, Vector*> _flux;
  std::map<int, Vector*> _power;
  std::map<int, Vector*> _weight;  

  /* Matrix and vector objects */
  Matrix* _A;
  Matrix* _M;
  Matrix* _AM;
  Vector* _b;

  
public:
  Solver(Geometry* geometry);
  virtual ~Solver();

  /* Getter functions */
  double getKeff0();
  double getBuckling();
  Geometry* getGeometry();

  Matrix* getAM();
  Matrix* getA();
  Matrix* getM();
  Vector* getb();
  
  Vector* getTemperature(int state);
  Vector* getFlux(int state);
  Vector* getPower(int state);
  Vector* getWeight(int state);
  
  double getTemperatureByValue(int cell, int state);
  double getFluxByValue(int cell, int group, int state);
  double getPowerByValue(int cell, int state);
  double getWeightByValue(int cell, int group, int state);
  
  /* Setter functions */
  void setBuckling(double buckling);
  void setInitialPower(double power);
  void setFluxSolveTolerance(double tolerance);
  void setEndTime(double time);
  
  void setTemperatureByValue(double value, int cell, int state);
  void setPowerByValue(double value, int cell, int state);
  void setFluxByValue(double value, int cell, int group, int state);
  void setWeightByValue(double value, int cell, int group, int state);

  /* Worker functions */
  void integratePrecursorConcentrations(int state_from, int state_to);
  void integrateTemperature(int state_from, int state_to);
  void computePrecursorConcentrations(int state);
  void computePower(int state);
  double computeAveragePower(int state);
  double computePowerRMSError(int state_1, int state_2);
  void normalizeFlux(double value, int state);
  void initializeClock();
  
  virtual void computeInitialFlux(double tol)=0;

  /* Copy functions */
  void copyPrecursors(int state_from, int state_to);
  virtual void copyVariables(int state_from, int state_to);
  void copyVariablesToAll(int state_from);
  
};

#endif /* SOLVER_H_ */
