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

/** Indexing macro for the scalar flux in each FSR and energy group */
#define _AM_shape(r,e) (_AM_shape[(r)*_shape_mesh->getNumShapeEnergyGroups()*(_shape_mesh->getNumShapeEnergyGroups()+6) + (e)])

/** Indexing macro for the scalar flux in each FSR and energy group */
#define _A_shape(r,e) (_A_shape[(r)*_shape_mesh->getNumShapeEnergyGroups()*(_shape_mesh->getNumShapeEnergyGroups()+6) + (e)])

/** Indexing macro for the scalar flux in each FSR and energy group */
#define _M_shape(r,e) (_M_shape[(r)*_shape_mesh->getNumShapeEnergyGroups()*_shape_mesh->getNumShapeEnergyGroups() + (e)])

/** Indexing macro for the scalar flux in each FSR and energy group */
#define _AM_amp(r,e) (_AM_amp[(r)*_amp_mesh->getNumAmpEnergyGroups()*(_amp_mesh->getNumAmpEnergyGroups()+6) + (e)])

/** Indexing macro for the scalar flux in each FSR and energy group */
#define _A_amp(r,e) (_A_amp[(r)*_amp_mesh->getNumAmpEnergyGroups()*(_amp_mesh->getNumAmpEnergyGroups()+6) + (e)])

/** Indexing macro for the scalar flux in each FSR and energy group */
#define _M_amp(r,e) (_M_amp[(r)*_amp_mesh->getNumAmpEnergyGroups()*_amp_mesh->getNumAmpEnergyGroups() + (e)])


class Solver {

private:

  double* _A_shape;
  double* _AM_shape;
  double* _M_shape;
  double* _b_shape;
  double* _A_amp;
  double* _AM_amp;
  double* _M_amp;
  double* _b_amp;

  double _k_eff;

  StructuredShapeMesh* _shape_mesh;
  AmpMesh* _amp_mesh;
  
public:
  Solver(StructuredShapeMesh* shape_mesh, AmpMesh* amp_mesh);
  virtual ~Solver();

  /* Getter functions */
  double* getAMShape();
  double* getAShape();
  double* getMShape();
  double* getBShape();
  double* getAMAmp();
  double* getAAmp();
  double* getMAmp();
  double* getBAmp();
  
  /* Setter functions */
  
  /* Worker functions */
  void makeAMShapeInitial(int position);
  void computeInitialShape(double tol);
  void makeAMAmp(double wt);
  void makeAMShape(double wt);
};

#endif /* SOLVER_H_ */
