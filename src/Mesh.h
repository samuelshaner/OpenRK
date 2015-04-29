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
 * @class Mesh Mesh.h "src/Mesh.h"
 * @brief The Mesh class represents a unique material and its relevant
 *        nuclear data (i.e., multigroup cross-sections) for neutron transport.
 */
class Mesh {

protected:

  /* mesh geometric dimensions */
  double _x_min;
  double _x_max;
  double _y_min;
  double _y_max;
  double _z_min;
  double _z_max;

  /* coarse mesh dimensions */
  int _num_x_coarse;
  int _num_y_coarse;
  int _num_z_coarse;

  double* _offset;
  boundaryType* _boundaries;
  int _num_energy_groups;
  int _num_delayed_groups;
  Clock* _clock;

  /* arrays of field variables */
  std::map<int, double*> _flux;
  std::map<int, double*> _shape;
  std::map<int, double*> _amplitude;
  std::map<int, double*> _power;
  std::map<int, double*> _temperature;

  double _fuel_volume;
  double _buckling;
  double* _decay_constants;
  double* _delayed_fractions;
  double _k_eff_0;
  Material** _materials;

  
public:
  Mesh(double width=1.0, double height=1.0, double depth=1.0);
  virtual ~Mesh();

  /* Setter functions */
  void setNumEnergyGroups(int num_groups);
  void setNumDelayedGroups(int num_groups);
  void setCoarseMeshDimensions(int num_x=1, int num_y=1, int num_z=1);
  void setXMin(double x_min);
  void setXMax(double x_max);
  void setYMin(double y_min);
  void setYMax(double y_max);
  void setZMin(double z_min);
  void setZMax(double z_max);
  void setKeff0(double k_eff_0);
  void setBoundary(int side, boundaryType boundary);
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
  int getNumEnergyGroups();
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
  int getNumXCoarse();
  int getNumYCoarse();
  int getNumZCoarse();
  double* getFlux(int position);
  double* getShape(int position);
  double* getAmplitude(int position);
  double* getPower(int position);
  double getFluxByValue(int cell, int group, int position);
  double getShapeByValue(int cell, int group, int position);
  double getAmplitudeByValue(int cell, int group, int position);
  double getPowerByValue(int cell, int position);
  Material* getMaterial(int cell);  

};

#endif /* MESH_H_ */
