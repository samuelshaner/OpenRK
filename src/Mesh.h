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
#include "constants.h"
#include <map>
#endif


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

  boundaryType* _boundaries;
  int _num_energy_groups;
  int _num_delayed_groups;
  Clock* _clock;
  double _buckling;
  double* _decay_constants;
  double* _delayed_fractions;
  double _k_eff_0;

  std::map<int, std::map<int, double*> > _field_variables;
  
public:
  Mesh(double width_x=1.0, double width_y=1.0, double width_z=1.0);
  virtual ~Mesh();

  /* Setter functions */
  void setNumEnergyGroups(int num_groups);
  void setNumDelayedGroups(int num_groups);
  virtual void setClock(Clock* clock);
  
  void setXMin(double x_min);
  void setXMax(double x_max);
  void setYMin(double y_min);
  void setYMax(double y_max);
  void setZMin(double z_min);
  void setZMax(double z_max);

  void setKeff0(double k_eff_0);
  void setBoundary(surfaceLocation side, boundaryType boundary);
  void setBuckling(double buckling);
  void setDecayConstants(double* xs, int num_groups);
  void setDelayedFractions(double* xs, int num_groups);
  
  /* Getter functions */
  double getKeff0();
  Clock* getClock();
  double getWidthX();
  double getWidthY();
  double getWidthZ();
  boundaryType getBoundary(surfaceLocation side);
  int getNumEnergyGroups();
  int getNumDelayedGroups();
  int getNumFieldVariableGroups();
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
  virtual int getNumCells()=0;
  
  double* getFieldVariable(fieldVariable name, int position);
  double* getFieldVariableByValue(fieldVariable name, int position, int cell, int group=0);
  
  /* Worker functions */
  void copyFieldVariable(int position_from, int position_to);
};

#endif /* MESH_H_ */
