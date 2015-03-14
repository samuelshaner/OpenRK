/**
 * @file FunctionalMaterial.h
 * @brief
 * @date February 8, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef FUNCTIONALMATERIAL_H_
#define FUNCTIONALMATERIAL_H_

#ifdef __cplusplus
#include "Material.h"
#endif


/**
 * @class FunctionalMaterial FunctionalMaterial.h "src/FunctionalMaterial.h"
 * @brief The FunctionalMaterial class represents a unique material and its relevant
 *        nuclear data (i.e., multigroup cross-sections) for neutron transport with
 *        cross sections that are an explicit function of time and temperature.
 */
class FunctionalMaterial : public Material {

protected:

  /** The number of time steps */
  int _num_time_steps;

  /** An array of the doppler coefficients for each energy group */
  double* _doppler_coefficients;

  /** An array of the time values for each step */
  double* _time_steps;

public:
  FunctionalMaterial(int id=0);
  virtual ~FunctionalMaterial();

  /* Setter functions */
  virtual void setNumEnergyGroups(const int num_groups);
  void setTimeSteps(double* time_steps, int num_steps);
  
  virtual void setSigmaT(double* xs, int num_time_steps, int num_groups);
  virtual void setSigmaA(double* xs, int num_time_steps, int num_groups);
  virtual void setSigmaS(double* xs, int num_time_steps, int num_groups);
  virtual void setSigmaF(double* xs, int num_time_steps, int num_groups);
  virtual void setNuSigmaF(double* xs, int num_time_steps, int num_groups);
  virtual void setChi(double* xs, int num_time_steps, int num_groups);
  virtual void setDifCoef(double* xs, int num_time_steps, int num_groups);
  virtual void setVelocity(double* xs, int num_time_steps, int num_groups);
  void setDopplerCoefficients(double* xs, int num_groups);

  virtual void setSigmaTByGroup(double xs, int group, int position=CURRENT);
  virtual void setSigmaAByGroup(double xs, int group, int position=CURRENT);
  virtual void setSigmaFByGroup(double xs, int group, int position=CURRENT);
  virtual void setNuSigmaFByGroup(double xs, int group, int position=CURRENT);
  virtual void setSigmaSByGroup(double xs, int group_from, int group_to, int position=CURRENT);
  virtual void setChiByGroup(double xs, int group, int position=CURRENT);
  virtual void setDifCoefByGroup(double xs, int group, int position=CURRENT);
  virtual void setVelocityByGroup(double xs, int group, int position=CURRENT);

  /* Getter functions */
  virtual double getSigmaTByGroup(int group, int position=CURRENT, double temp=0.0);
  virtual double getSigmaAByGroup(int group, int position=CURRENT, double temp=0.0);
  virtual double getSigmaSByGroup(int group_from, int group_to, int position=CURRENT, double temp=0.0);
  virtual double getSigmaFByGroup(int group, int position=CURRENT, double temp=0.0);
  virtual double getNuSigmaFByGroup(int group, int position=CURRENT, double temp=0.0);
  virtual double getChiByGroup(int group, int position=CURRENT, double temp=0.0);
  virtual double getDifCoefByGroup(int group, int position=CURRENT, double temp=0.0);
  virtual double getVelocityByGroup(int group, int position=CURRENT, double temp=0.0);
  double getDopplerCoefficientByGroup(int group);
  
  int getTimeStep(int position);
  
  /* Worker functions */
  virtual std::string toString();
  virtual FunctionalMaterial* clone();
  virtual void copy(int position_from, int position_to);
};

#endif /* FUNCTIONALMATERIAL_H_ */
