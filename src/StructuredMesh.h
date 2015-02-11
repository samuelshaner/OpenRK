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
  double _cell_width;
  double _cell_height;
  
  std::map<int, double*> _current;
  std::map<int, double*> _dif_linear;
  std::map<int, double*> _dif_nonlinear;
  

public:
  StructuredMesh(double width=1.0, double height=1.0, int num_x=1, int num_y=1);
  virtual ~StructuredMesh();

  /* Setter functions */
  void setNumX(int num_x);
  void setNumY(int num_y);
  void setClock(Clock* clock);
  
  /* Getter functions */
  double* getCurrent(int position);
  double* getDifLinear(int position);
  double* getDifNonlinear(int position);
  int getNeighborCell(int x, int y, int side);
  Material* getNeighborMaterial(int x, int y, int side);
  int getNumX();
  int getNumY();
  double getCellWidth();
  double getCellHeight();
  double getCellVolume();
  int findCell(double x, double y);
  
  /* Worker functions */
  void computeFuelVolume();
  void uniquifyMaterials();
  double getMaxTemperature(int position);
  void copyPower(int position_from, int position_to);
  void copyTemperature(int position_from, int position_to);
  void setTemperature(double temperature);
};

#endif /* STRUCTUREDMESH_H_ */
