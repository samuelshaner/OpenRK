/**
 * @file Geometry.h
 * @brief
 * @date April 12, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#ifdef __cplusplus
#include "Material.h"
#include "linalg.h"
#include <map>
#endif

/**
 * @enum boundaryType
 * @brief The boundary types
 */
enum boundaryType {
  REFLECTIVE,
  VACUUM
};


/**
 * @class Geometry Geometry.h "src/Geometry.h"
 * @brief The Geometry class represents a unique material and its relevant
 *        nuclear data (i.e., multigroup cross-sections) for neutron transport.
 */
class Geometry {

protected:

  /* geometric dimensions */
  double _x_min;
  double _x_max;
  double _y_min;
  double _y_max;
  double _z_min;
  double _z_max;

  /* amplitude mesh dimensions */
  int _num_x_amp;
  int _num_y_amp;
  int _num_z_amp;

  /* geometry properties */
  int _num_shape_cells;
  int _num_amp_cells;
  Material** _materials;
  boundaryType* _boundaries;
  int _num_energy_groups;
  double* _volumes;

  /* mapping of cells between shape and amplitude */
  std::vector< std::vector<int> > _amp_to_shape;
  int* _shape_to_amp;
  
public:
  Geometry(double width=1.0, double height=1.0, double depth=1.0);
  virtual ~Geometry();

  /* Setter functions */
  void setAmpMeshDimensions(int num_x=1, int num_y=1, int num_z=1);
  void setXMin(double x_min);
  void setXMax(double x_max);
  void setYMin(double y_min);
  void setYMax(double y_max);
  void setZMin(double z_min);
  void setZMax(double z_max);
  void setBoundary(int side, boundaryType boundary);
  void setMaterial(Material* material, int cell);
  void setNumShapeCells(int num_shape_cells);
  
  /* Getter functions */
  double getWidth();
  double getHeight();
  double getDepth();
  boundaryType getBoundary(int side);
  double getXMin();
  double getXMax();
  double getYMin();
  double getYMax();
  double getZMin();
  double getZMax();
  int getNumXAmp();
  int getNumYAmp();
  int getNumZAmp();
  Material* getMaterial(int cell);
  std::vector< std::vector<int> > getAmpToShapeMap();
  int* getShapeToAmpMap();
  int getNumEnergyGroups();
  int getNeighborAmpCell(int x, int y, int z, int side);
  double getVolume(int cell);

  /* Worker functions */
  void addShapeCellToAmpCell(int shape_cell, int amp_cell);
  int findAmpCellContainingShapeCell(int shape_cell);
  void uniquifyMaterials();
  virtual Geometry* clone();
};

#endif /* GEOMETRY_H_ */
