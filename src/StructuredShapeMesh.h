/**
 * @file StructuredShapeMesh.h
 * @brief
 * @date February 8, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef STRUCTUREDSHAPEMESH_H_
#define STRUCTUREDSHAPEMESH_H_

#ifdef __cplusplus
#include "AmpMesh.h"
#include "ShapeMesh.h"
#endif

/**
 * @class StructuredShapeMesh StructuredShapeMesh.h "src/StructuredShapeMesh.h"
 * @brief The StructuredShapeMesh class represents a unique material and its relevant
 *        nuclear data (i.e., multigroup cross-sections) for neutron transport.
 */
class AmpMesh;
class StructuredShapeMesh : public ShapeMesh {

protected:

  int _num_x;
  int _num_y;
  int _num_z;
  double _cell_width;
  double _cell_height;
  double _cell_depth;

  std::map<int, double*> _dif_linear;
  std::map<int, double*> _dif_nonlinear;

public:
  StructuredShapeMesh(double width=1.0, double height=1.0, double depth=1.0, 
                      int num_x=1, int num_y=1, int num_z=1);
  virtual ~StructuredShapeMesh();

  /* Setter functions */
  void setDifLinearByValue(double dif_linear, int cell, int group, int side, 
                           int position);
  void setNumX(int num_x);
  void setNumY(int num_y);
  void setNumZ(int num_z);  
  
  /* Getter functions */
  double getDifLinearByValue(int cell, int group, int side, int position);
  int getNumX();
  int getNumY();
  int getNumZ();
  double getCellWidth();
  double getCellHeight();
  double getCellDepth();
  virtual int getNumCells();
  virtual double getCellVolume(int cell=0);
  
  /* Worker functions */
  StructuredShapeMesh* clone();
  StructuredShapeMesh* uniformRefine(int refine_x=1, int refine_y=1, int refine_z=1);
  virtual void initialize();
  void computeDifCoefs(int position);
  void copyDifLinear(int position_from, int position_to);
  int findCell(double x, double y, double z);
  
};

#endif /* STRUCTUREDSHAPEMESH_H_ */
