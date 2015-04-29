/**
 * @file StructuredShapeMesh.h
 * @brief
 * @date February 8, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef STRUCTUREDSHAPEMESH_H_
#define STRUCTUREDSHAPEMESH_H_

#ifdef __cplusplus
#include "StructuredMesh.h"
#include "AmpMesh.h"
#endif

/**
 * @class StructuredShapeMesh StructuredShapeMesh.h "src/StructuredShapeMesh.h"
 * @brief The StructuredShapeMesh class represents a unique material and its relevant
 *        nuclear data (i.e., multigroup cross-sections) for neutron transport.
 */
class AmpMesh;
class StructuredShapeMesh : public StructuredMesh {

protected:

  AmpMesh* _amp_mesh;
  int* _group_indices;
  int* _amp_map;

public:
  StructuredShapeMesh(double width=1.0, double height=1.0, double depth=1.0, int num_x=1, int num_y=1, int num_z=1);
  virtual ~StructuredShapeMesh();

  /* Setter functions */
  void setAmpMesh(AmpMesh* mesh);
  void setFluxByValue(double flux, int cell, int group, int position);
  void setCurrentByValue(double current, int cell, int group, int side, int position);
  void setDifLinearByValue(double dif_linear, int cell, int group, int side, int position);
  void setDifNonlinearByValue(double dif_nonlinear, int cell, int group, int side, int position);
  void setGroupStructure(int* group_indices, int length_group_indices);
  
  /* Getter functions */
  double getFluxByValue(int cell, int group, int position);
  double getCurrentByValue(int cell, int group, int side, int position);
  double getDifLinearByValue(int cell, int group, int side, int position);
  double getDifNonlinearByValue(int cell, int group, int side, int position);

  int getAmpGroup(int shape_group);
  
  /* Worker functions */
  StructuredShapeMesh* clone();
  StructuredShapeMesh* uniformRefine(int refine_x=1, int refine_y=1, int refine_z=1);
  void initialize();
  void synthesizeFlux(int position);
  void reconstructFlux(int position, int position_shape, int position_amp);
  void computePower(int position);
  void computeInitialPrecursorConc(int position);
  void integratePrecursorConc(int position_from, int position_to);
  void integrateTemperature(int position_from, int position_to);
  void computeDifCoefs(int position);
  void copyFlux(int position_from, int position_to);
  void copyCurrent(int position_from, int position_to);
  void copyDifLinear(int position_from, int position_to);
  void copyDifNonlinear(int position_from, int position_to);
  void scaleFlux(int position, double scale_val);
  double computeAveragePower(int position);
  double computePowerL2Norm(int position_1, int position_2);
  int findAmpCell(int shape_cell);
  void saveShape();  
};

#endif /* STRUCTUREDSHAPEMESH_H_ */
