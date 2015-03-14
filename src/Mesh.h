/**
 * @file Mesh.h
 * @brief
 * @date February 8, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef MESH_H_
#define MESH_H_

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
 * @enum meshType
 * @brief The mesh types
 */
enum meshType {
  AMPMESH,
  STRUCTURED_SHAPE_MESH,
  UNSTRUCTURED_SHAPE_MESH
};


/**
 * @class Mesh Mesh.h "src/Mesh.h"
 * @brief The Mesh class represents a unique material and its relevant
 *        nuclear data (i.e., multigroup cross-sections) for neutron transport.
 */
class Mesh {

protected:

  double _x_min;
  double _x_max;
  double _y_min;
  double _y_max;
  double _z_min;
  double _z_max;

  double _fuel_volume;
  double _buckling;
  double* _decay_constants;
  double* _delayed_fractions;
  double _k_eff_0;
  int _num_shape_energy_groups;
  int _num_amp_energy_groups;
  int _num_delayed_groups;
  Clock* _clock;

  std::map<int, double*> _flux;
  std::map<int, double*> _temperature;
  std::map<int, double*> _power;
  Material** _materials;
  boundaryType* _boundaries;
  meshType _mesh_type;
  
public:
  Mesh(double width=1.0, double height=1.0, double depth=1.0);
  virtual ~Mesh();

  /* Setter functions */
  void setNumShapeEnergyGroups(int num_groups);
  void setNumAmpEnergyGroups(int num_groups);
  void setNumDelayedGroups(int num_groups);

  void setXMin(double x_min);
  void setXMax(double x_max);
  void setYMin(double y_min);
  void setYMax(double y_max);
  void setZMin(double z_min);
  void setZMax(double z_max);

  void setKeff0(double k_eff_0);
  void setBoundary(int side, boundaryType boundary);
  void setMaterial(Material* material, int cell);  
  void setClock(Clock* clock);

  void setBuckling(double buckling);
  void setDecayConstants(double* xs, int num_groups);
  void setDelayedFractions(double* xs, int num_groups);
  
  /* Getter functions */
  double getKeff0();
  Clock* getClock();
  double getWidth();
  double getHeight();
  double getDepth();
  boundaryType getBoundary(int side);
  meshType getMeshType();
  int getNumShapeEnergyGroups();
  int getNumAmpEnergyGroups();
  int getNumDelayedGroups();
  double getBuckling();
  double* getDecayConstants();
  double* getDelayedFractions();
  double getDecayConstantByGroup(int group);
  double getDelayedFractionByGroup(int group);
  double getXMin();
  double getXMax();
  double getYMin();
  double getYMax();
  double getZMin();
  double getZMax();
  virtual int getNumCells() = 0;
  virtual double getCellVolume(int cell) = 0;
  
  double* getFlux(int position);
  double* getPower(int position);
  double* getTemperature(int position);
  double getPowerByValue(int cell, int position);
  double getTemperatureByValue(int cell, int position);

  Material* getMaterial(int cell);
  
  /* Worker functions */
  void computeFuelVolume();
  void uniquifyMaterials();
  void copyPower(int position_from, int position_to);
  void copyTemperature(int position_from, int position_to);
  virtual void initialize();
};

#endif /* MESH_H_ */
