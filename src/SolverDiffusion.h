/**
 * @file SolverDiffusion.h
 * @brief A solver object to solver for the flux.
 * @date February 8, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef SOLVERDIFFUSION_H_
#define SOLVERDIFFUSION_H_

#ifdef __cplusplus
#include "Solver.h"
#include "GeometryDiffusion.h"
#endif


class SolverDiffusion : public Solver {

protected:

  /* Matrix and vector objects */
  Matrix* _shape_A;
  Matrix* _shape_M;
  Matrix* _shape_AM;
  Vector* _shape_b;

  /* Fine mesh field variables */
  GeometryDiffusion* _geometry_diffusion;

public:
  SolverDiffusion(GeometryDiffusion* geometry);
  virtual ~SolverDiffusion();

  /* Getter functions */
  Matrix* getShapeA();
  Matrix* getShapeM();
  Matrix* getShapeAM();
  Vector* getShapeb();
  GeometryDiffusion* getGeometryDiffusion();
  
  /* Setter functions */

  /* Worker functions */
  virtual void takeStep()=0;
  virtual void computeInitialFlux(double tol)=0;
};

#endif /* SOLVERDIFFUSION_H_ */
