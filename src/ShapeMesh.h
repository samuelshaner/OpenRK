/**
 * @file ShapeMesh.h
 * @brief
 * @date February 20, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef SHAPEMESH_H_
#define SHAPEMESH_H_

#ifdef __cplusplus
#include "Mesh.h"
#include "AmpMesh.h"
#endif

/**
 * @class ShapeMesh ShapeMesh.h "src/ShapeMesh.h"
 * @brief The ShapeMesh class represents a unique material and its relevant
 *        nuclear data (i.e., multigroup cross-sections) for neutron transport.
 */
class AmpMesh;
class ShapeMesh : public Mesh {

protected:

  AmpMesh* _amp_mesh;
  int* _group_indices;
  int* _amp_map;

public:
  ShapeMesh(double width=1.0, double height=1.0, double depth=1.0);
  virtual ~ShapeMesh();

  /* Setter functions */
  void setAmpMesh(AmpMesh* mesh);
  void setFluxByValue(double flux, int cell, int group, int position);
  void setGroupStructure(int* group_indices, int length_group_indices);
  void setFlux(double* flux, int num_cells_times_groups);
  void setAmpCellContainingShapeCell(int shape_cell, int amp_cell);

  /* Getter functions */
  double getFluxByValue(int cell, int group, int position);
  int getAmpGroup(int shape_group);
  virtual int getNumCells() = 0;
  virtual double getCellVolume(int cell) = 0;
  
  /* Worker functions */
  virtual void initialize();
  void synthesizeFlux(int position);
  void reconstructFlux(int position, int position_shape, int position_amp);
  void computePower(int position);
  void computeInitialPrecursorConc(int position);
  void integratePrecursorConc(int position_from, int position_to);
  void integrateTemperature(int position_from, int position_to);
  void copyFlux(int position_from, int position_to);
  void scaleFlux(int position, double scale_val);
  double computeAveragePower(int position);
  double computePowerL2Norm(int position_1, int position_2);
  
};

#endif /* SHAPEMESH_H_ */
