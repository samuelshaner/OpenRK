/**
 * @file AmpMesh.h
 * @brief
 * @date February 8, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef AMPMESH_H_
#define AMPMESH_H_

#ifdef __cplusplus
#include "Mesh.h"
#include "ShapeMesh.h"
#include <vector>
#endif

/**
 * @class AmpMesh AmpMesh.h "src/AmpMesh.h"
 * @brief The AmpMesh class represents a unique material and its relevant
 *        nuclear data (i.e., multigroup cross-sections) for neutron transport.
 */
class StructuredShapeMesh;
class UnstructuredShapeMesh;
class AmpMesh : public Mesh {

protected:

  int _num_x;
  int _num_y;
  int _num_z;
  double _cell_width;
  double _cell_height;
  double _cell_depth;

  ShapeMesh* _shape_mesh;
  bool _optically_thick;
  int* _group_indices;
  
  std::map<int, double*> _current;
  std::map<int, double*> _dif_linear;
  std::map<int, double*> _dif_nonlinear;
  std::vector< std::vector<int> > _shape_map;
  
public:
  AmpMesh(double width=1.0, double height=1.0, double depth=1.0, 
          int num_x=1, int num_y=1, int num_z=1);
  virtual ~AmpMesh();

  /* Setter functions */
  void setNumX(int num_x);
  void setNumY(int num_y);
  void setNumZ(int num_z);
  void setOpticallyThick(bool optically_thick);
  void setShapeMesh(ShapeMesh* mesh);
  void setShapeMap(int** regions_to_cells, int* regions_per_cell, int num_cells);
  void setFluxByValue(double flux, int cell, int group, int position);
  void setCurrentByValue(double current, int cell, int group, int side, int position);
  void setDifLinearByValue(double dif_linear, int cell, int group, int side, 
                           int position);
  void setDifNonlinearByValue(double dif_nonlinear, int cell, int group, int side, 
                              int position);
  void setGroupStructure(int* group_indices, int length_group_indices);
  void setCurrent(double* current, int num_cellss_times_groups);
  
  /* Getter functions */
  double* getCurrent(int position);
  double* getDifLinear(int position);
  double* getDifNonlinear(int position);
  double getFluxByValue(int cell, int group, int position);
  double getCurrentByValue(int cell, int group, int side, int position);
  double getDifLinearByValue(int cell, int group, int side, int position);
  double getDifNonlinearByValue(int cell, int group, int side, int position);
  int getNeighborCell(int x, int y, int z, int side);
  Material* getNeighborMaterial(int x, int y, int z, int side);

  int getNumX();
  int getNumY();
  int getNumZ();
  double getCellWidth();
  double getCellHeight();
  double getCellDepth();
  double getCellVolume(int cell=0);
  int getNumCells();
  
  /* Worker functions */
  AmpMesh* clone();
  void initialize();
  void interpolateDifNonlinear(int position_begin, int position_end, int position);
  int findCell(double x, double y, double z);
 
  void condenseMaterials(int position, bool save_flux=false);
  void computePower(int position);
  void computeCurrent(int position);
  double computeDifCorrect(double dif_coef, double length);
  void computeDifCoefs(int position);
  double computeAveragePower(int position);
  double computePowerL2Norm(int position_1, int position_2);

  void copyFlux(int position_from, int position_to);
  void copyCurrent(int position_from, int position_to);
  void copyDifLinear(int position_from, int position_to);
  void copyDifNonlinear(int position_from, int position_to);
  
};

#endif /* AMPMESH_H_ */
