/**
 * @file SolverDiffusionTheta.h
 * @brief A solver object to transient problems using the direct theta methods.
 * @date February 8, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef SOLVERDIFFUSIONTHETA_H_
#define SOLVERDIFFUSIONTHETA_H_

#ifdef __cplusplus
#include <sstream>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <map>
#include <vector>
#include "SolverDiffusion.h"
#endif


class SolverDiffusionTheta : public SolverDiffusion {

protected:
  
public:
  SolverDiffusionTheta(GeometryDiffusion* geometry);
  virtual ~SolverDiffusionTheta();

  /* Setter functions */
  void setStepSize(double time);
  
  /* Getter functions */
  double getSurfaceDifCoef(int cell, int group, int surface, int state);
  
  /* Worker functions */
  void generateInitialMatrices();
  void generateMatrices();
  void computeInitialFlux(double tol);
  void takeStep();
  
};

#endif /* SOLVERDIFFUSIONTHETA_H_ */
