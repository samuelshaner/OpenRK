/**
 * @file AmpMesh.h
 * @brief
 * @date February 8, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef AMPMESH_H_
#define AMPMESH_H_

#ifdef __cplusplus
#include "StructuredMesh.h"
#include "StructuredShapeMesh.h"
#include <vector>
#endif

/**
 * @class AmpMesh AmpMesh.h "src/AmpMesh.h"
 * @brief The AmpMesh class represents a unique material and its relevant
 *        nuclear data (i.e., multigroup cross-sections) for neutron transport.
 */
class StructuredShapeMesh;
class AmpMesh : public StructuredMesh {

protected:

  StructuredShapeMesh* _shape_mesh;
  bool _optically_thick;
  int* _group_indices;
  std::vector< std::vector<int> > _shape_map;
  double _energy_per_fission;
  
public:
  AmpMesh(double width=1.0, double height=1.0, double depth=1.0, int num_x=1, int num_y=1, int num_z=1);
  virtual ~AmpMesh();

  /* Setter functions */
  void setOpticallyThick(bool optically_thick);
  void setShapeMesh(StructuredShapeMesh* mesh);
  void setFluxByValue(double flux, int cell, int group, int position);
  void setCurrentByValue(double current, int cell, int group, int side, int position);
  void setDifLinearByValue(double dif_linear, int cell, int group, int side, int position);
  void setDifNonlinearByValue(double dif_nonlinear, int cell, int group, int side, int position);
  void setGroupStructure(int* group_indices, int length_group_indices);
  void setEnergyPerFission(double _energy_per_fission);
  
  /* Getter functions */
  double getFluxByValue(int cell, int group, int position);
  double getCurrentByValue(int cell, int group, int side, int position);
  double getDifLinearByValue(int cell, int group, int side, int position);
  double getDifNonlinearByValue(int cell, int group, int side, int position);
  
  /* Worker functions */
  AmpMesh* clone();
  void initialize();
  void condenseMaterials(int position, bool save_flux=false);
  void computePower(int position);
  void computeCurrent(int position);
  double computeDifCorrect(double dif_coef, double length);
  void computeDifCoefs(int position);
  void copyFlux(int position_from, int position_to);
  void copyCurrent(int position_from, int position_to);
  void copyDifLinear(int position_from, int position_to);
  void copyDifNonlinear(int position_from, int position_to);
  double computeAveragePower(int position);
  double computePowerL2Norm(int position_1, int position_2);
  void interpolateDifNonlinear(int position_begin, int position_end, int position);
  
};

#endif /* AMPMESH_H_ */
