/**
 * @file UnstructuredShapeMesh.h
 * @brief
 * @date February 20, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef UNSTRUCTUREDSHAPEMESH_H_
#define UNSTRUCTUREDSHAPEMESH_H_

#ifdef __cplusplus
#include "ShapeMesh.h"
#include "AmpMesh.h"
#endif

/**
 * @class UnstructuredShapeMesh UnstructuredShapeMesh.h "src/UnstructuredShapeMesh.h"
 * @brief The UnstructuredShapeMesh class represents a unique material and its relevant
 *        nuclear data (i.e., multigroup cross-sections) for neutron transport.
 */
class AmpMesh;
class UnstructuredShapeMesh : public ShapeMesh {

protected:

  double* _volume;

public:
  UnstructuredShapeMesh(double width=1.0, double height=1.0, double depth=1.0, int num_cells=1);
  virtual ~UnstructuredShapeMesh();

  /* Setter functions */
  void setTemperature(double* temperature, int num_cells);
  void setVolume(double* volume, int num_cells);
  
  /* Getter functions */
  virtual int getNumCells();
  virtual double getCellVolume(int cell);  

  /* Worker functions */
  UnstructuredShapeMesh* clone();
  
};

#endif /* UNSTRUCTUREDSHAPEMESH_H_ */
