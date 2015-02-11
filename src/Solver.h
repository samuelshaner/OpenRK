/**
 * @file Solver.h
 * @brief A solver object to solver for the flux.
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
#include "StructuredShapeMesh.h"
#include "omp.h"
#include "linalg.h"
#endif


class Solver {

private:

  double** _A_shape;
  double** _AM_shape;
  double** _M_shape;
  double* _b_shape;
  double** _A_amp;
  double** _AM_amp;
  double** _M_amp;
  double* _b_amp;

  double _k_eff;

  StructuredShapeMesh* _shape_mesh;
  AmpMesh* _amp_mesh;
  
public:
  Solver(StructuredShapeMesh* shape_mesh, AmpMesh* amp_mesh);
  virtual ~Solver();

  /* Getter functions */
  double** getAMShape();
  double** getAShape();
  double** getMShape();
  double* getBShape();
  double** getAMAmp();
  double** getAAmp();
  double** getMAmp();
  double* getBAmp();
  
  /* Setter functions */
  
  /* Worker functions */
  void makeAMShapeInitial(int position);
  void computeInitialShape(double tol);
  void makeAMAmp(double wt);
  void makeAMShape(double wt);
};

#endif /* SOLVER_H_ */
