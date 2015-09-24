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

  Matrix* _shape_A_matrix;
  Matrix* _shape_M_matrix;
  Matrix* _shape_AM_matrix;
  Vector* _shape_source;

  /* Fine mesh field variables */
  std::map<int, Vector*> _dif_linear_fine;

  GeometryDiffusion* _geometry_diffusion;

public:
  SolverDiffusion(GeometryDiffusion* geometry);
  virtual ~SolverDiffusion();

  /* Getter functions */
  Matrix* getShapeAMatrix();
  Matrix* getShapeMMatrix();
  Matrix* getShapeAMMatrix();
  Vector* getShapeSource();
  Vector* getDifLinearFine(int state);
  double getDifLinearFineByValue(int cell, int group, int side, int state);
  GeometryDiffusion* getGeometryDiffusion();
  
  /* Setter functions */
  void setDifLinearFineByValue(double value, int cell, int group, int side, int state);  

  /* Worker functions */
  virtual void takeInnerStep();
  virtual void takeOuterStep();
  virtual void takeOuterStepOnly();
  virtual void computeInitialShape(double tol);
  void generateShapeMatrices(int state, int state_prev);
  void computeDiffusionCoefficientsFine(int state);
  void generateAdiabaticShapeMatrices();
  void generateAmpCurrent(int state);

};

#endif /* SOLVERDIFFUSION_H_ */
