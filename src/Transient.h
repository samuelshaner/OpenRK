/**
 * @file Transient.h
 * @brief A solver object to solver for the flux.
 * @date February 8, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef TRANSIENT_H_
#define TRANSIENT_H_

#ifdef __cplusplus
#include "Solver.h"
#endif


class Transient {

private:

  double _initial_power;
  Clock* _clock;
  StructuredShapeMesh* _shape_mesh;
  AmpMesh* _amp_mesh;
  Solver* _solver;
  double _k_eff_0;

  
public:
  Transient();
  virtual ~Transient();

  /* Getter functions */
  
  /* Setter functions */
  void setInitialPower(double initial_power);
  void setClock(Clock* clock);
  void setSolver(Solver* solver);
  void setShapeMesh(StructuredShapeMesh* mesh);
  void setAmpMesh(AmpMesh* mesh);
    
  /* Worker functions */
  void computeInitialShape();
  void takeInnerStep();
  void takeOuterStep(double tol=1.e-6);
  void takeOuterStepOnly(double tol=1.e-6);
  void broadcastToActive(clockPosition position);
  void broadcastToAll(clockPosition position);
  void broadcastToOne(clockPosition position_from, clockPosition position_to);  
};

#endif /* TRANSIENT_H_ */
