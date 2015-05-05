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

  double** _shape_A_matrix;
  double** _shape_M_matrix;
  double** _shape_AM_matrix;
  double* _shape_source;

  /* Fine mesh field variables */
  std::map<int, double*> _dif_linear_fine;

  GeometryDiffusion* _geometry_diffusion;

public:
  SolverDiffusion(Geometry* geometry);
  virtual ~SolverDiffusion();

  /* Getter functions */
  double** getShapeAMatrix();
  double** getShapeMMatrix();
  double** getShapeAMMatrix();
  double* getShapeSource();
  double* getDifLinearFine(int time);
  double getDifLinearFineByValue(int cell, int group, int side, int time);

  /* Setter functions */
  void setDifLinearFineByValue(double value, int cell, int group, int time);  

  /* Worker functions */
  virtual void takeInnerStep();
  virtual void takeOuterStep();
  virtual void takeOuterStepOnly();
  virtual void computeInitialShape();
  void generateShapeMatrices();
  void computeDiffusionCoefficientsFine(int time);
  void generateInitialShapeMatrices();
  void generateAmpCurrent(int time);

};

#endif /* SOLVERDIFFUSION_H_ */
