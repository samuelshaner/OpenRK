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


/**
 * @enum transientMethod
 * @brief The supported transientMethods
 */
enum transientMethod {
  FORWARD_EULER,
  BACKWARD_EULER,
  CRANK_NICOLSON,
  CUSTOM
};

class Transient {

private:

  transientMethod _inner_method;
  transientMethod _outer_method;
  double _inner_wt;
  double _outer_wt;
  double _initial_power;
  Clock* _clock;
  StructuredShapeMesh* _shape_mesh;
  AmpMesh* _amp_mesh;
  Solver* _solver;
  double _k_eff_0;

  
public:
  Transient(transientMethod inner_method=BACKWARD_EULER, transientMethod outer_method=BACKWARD_EULER, double wt_inner=1.0, double wt_outer=1.0);
  virtual ~Transient();

  /* Getter functions */
  
  /* Setter functions */
  void setInnerMethod(transientMethod inner_method, double wt=0.0);
  void setOuterMethod(transientMethod outer_method, double wt=0.0);
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
