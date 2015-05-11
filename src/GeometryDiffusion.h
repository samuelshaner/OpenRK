/**
 * @file GeometryDiffusion.h
 * @brief
 * @date April 12, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef GEOMETRYDIFFUSION_H_
#define GEOMETRYDIFFUSION_H_

#ifdef __cplusplus
#include "Geometry.h"
#endif


/**
 * @class GeometryDiffusion GeometryDiffusion.h "src/GeometryDiffusion.h"
 * @brief The Geometry class represents a unique material and its relevant
 *        nuclear data (i.e., multigroup cross-sections) for neutron transport.
 */
class GeometryDiffusion : public Geometry {

protected:

  /* fine mesh dimensions */
  int _num_x_shape;
  int _num_y_shape;
  int _num_z_shape;
  
public:
  GeometryDiffusion(double width=1.0, double height=1.0, double depth=1.0);
  virtual ~GeometryDiffusion();

  /* Setter functions */
  void setShapeMeshDimensions(int num_x=1, int num_y=1, int num_z=1);
  
  /* Getter functions */
  int getNumXShape();
  int getNumYShape();
  int getNumZShape();
  int getNeighborShapeCell(int x, int y, int z, int side);

  /* Worker functions */
  GeometryDiffusion* uniformRefine(int refine_x=1, int refine_y=1, int refine_z=1);
  virtual GeometryDiffusion* clone();
  void generateCellMap();
};

#endif /* GEOMETRYDIFFUSION_H_ */
