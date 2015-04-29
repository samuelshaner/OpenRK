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
  std::vector< std::vector<int> > _shape_map;
  
public:
  AmpMesh(double width=1.0, double height=1.0, double depth=1.0, int num_x=1, int num_y=1, int num_z=1);
  virtual ~AmpMesh();

  /* Setter functions */
  void setShapeMesh(StructuredShapeMesh* mesh);
  void setFluxByValue(double flux, int cell, int group, int position);
  
  /* Getter functions */
  double getFluxByValue(int cell, int group, int position);
  
  /* Worker functions */
  AmpMesh* clone();
  void initialize();
  void computePower(int position);
  void copyFlux(int position_from, int position_to);
  double computeAveragePower(int position);
  double computePowerL2Norm(int position_1, int position_2);
  
};

#endif /* AMPMESH_H_ */
