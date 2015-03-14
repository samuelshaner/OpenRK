/**
 * @file StructuredMesh.h
 * @brief
 * @date February 8, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef STRUCTUREDMESH_H_
#define STRUCTUREDMESH_H_

#ifdef __cplusplus
#include "Mesh.h"
#include <vector>
#endif

/**
 * @class StructuredMesh StructuredMesh.h "src/StructuredMesh.h"
 * @brief The StructuredMesh class represents a unique material and its relevant
 *        nuclear data (i.e., multigroup cross-sections) for neutron transport.
 */
class StructuredMesh : public Mesh {

protected:

  int _num_x;
  int _num_y;
  int _num_z;
  double _cell_width_x;
  double _cell_width_y;
  double _cell_width_z;
  

public:
  StructuredMesh(double width_x=1.0, double width_y=1.0, double width_z=1.0, \
                 int num_x=1, int num_y=1, int num_z=1);
  virtual ~StructuredMesh();

  /* Setter functions */
  void setNumX(int num_x);
  void setNumY(int num_y);
  void setNumZ(int num_z);
  
  /* Getter functions */
  int getNeighborCell(int x, int y, int z, int side);
  int getNumX();
  int getNumY();
  int getNumZ();
  int getNumCells();
  double getCellWidthX();
  double getCellWidthY();
  double getCellWidthZ();
  double getCellVolume();
  int findCell(double x, double y, double z);
  
  /* Worker functions */
};

#endif /* STRUCTUREDMESH_H_ */
