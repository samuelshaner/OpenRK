/**
 * @file Material.h
 * @brief
 * @date February 7, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef MATERIAL_H_
#define MATERIAL_H_

#ifdef __cplusplus
#include <sstream>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "log.h"
#include <map>
#include "Clock.h"
#endif

int material_id();


/**
 * @class Material Material.h "src/Material.h"
 * @brief The Material class represents a unique material and its relevant
 *        nuclear data (i.e., multigroup cross-sections) for neutron transport.
 */
class Material {

protected:

  /** A user-defined ID for each Material created */
  int _id;

  /** The number of energy groups */
  int _num_energy_groups;

  /** An array of the total cross-sections for each energy group */
  double* _sigma_t;

  /** An array of the absorption cross-sections for each energy group */
  double* _sigma_a;

  /** A 2D array of the scattering cross-section matrix. The first index is
   *  row number and second index is column number */
  double* _sigma_s;

  /** An array of the fission cross-sections for each energy group */
  double* _sigma_f;

  /** An array of the fission cross-sections multiplied by nu \f$ \nu \f$
   *  for each energy group */
  double* _nu_sigma_f;

  /** An array of the chi \f$ \chi \f$ values for each energy group */
  double* _chi;

  /** An array of the diffusion coefficients for each energy group */
  double* _dif_coef;

  /** An array of the velocities for each energy group */
  double* _velocity;

  /** The energy released per fission in J/fission */
  double _energy_per_fission;

  /** A boolean representing whether or not this Material contains a non-zero
   *  fission cross-section and is fissionable */
  bool _fissionable;

  /** The number of delayed groups */
  int _num_delayed_groups;

  /** An array of the delayed neutron precursor concentration for each delayed group */
  double* _precursor_conc;

  double* _decay_constant;
  double* _delayed_fraction;

  /** A pointer to a Clock used to keep track of current state of transient */
  Clock* _clock;

  /** A conversion factor to convert energy to change in temperature */
  double _temperature_conversion_factor;  
  
public:
  Material(int id=0);
  virtual ~Material();

  /* Setter functions */
  void setNumEnergyGroups(const int num_groups);
  void setNumDelayedGroups(const int num_groups);
  void setEnergyPerFission(double energy_per_fission);
  
  virtual void setSigmaT(double* xs, int num_groups);
  virtual void setSigmaA(double* xs, int num_groups);
  virtual void setSigmaS(double* xs, int num_groups);
  virtual void setSigmaF(double* xs, int num_groups);
  virtual void setNuSigmaF(double* xs, int num_groups);
  virtual void setChi(double* xs, int num_groups);
  virtual void setDifCoef(double* xs, int num_groups);
  virtual void setVelocity(double* xs, int num_groups);
  virtual void setPrecursorConc(double* xs, int num_groups);
  virtual void setDecayConstant(double* decay_constant, int num_groups);
  virtual void setDelayedFraction(double* delayed_fraction, int num_groups);

  virtual void setSigmaTByGroup(double xs, int group, int position=CURRENT);
  virtual void setSigmaAByGroup(double xs, int group, int position=CURRENT);
  virtual void setSigmaFByGroup(double xs, int group, int position=CURRENT);
  virtual void setNuSigmaFByGroup(double xs, int group, int position=CURRENT);
  virtual void setSigmaSByGroup(double xs, int group_from, int group_to, int position=CURRENT);
  virtual void setChiByGroup(double xs, int group, int position=CURRENT);
  virtual void setDifCoefByGroup(double xs, int group, int position=CURRENT);
  virtual void setVelocityByGroup(double xs, int group, int position=CURRENT);
  virtual void setPrecursorConcByGroup(double xs, int group, int position=CURRENT);
  void setTemperatureConversionFactor(double conversion_factor);
  virtual void setDecayConstantByGroup(double xs, int group, int position=CURRENT);
  virtual void setDelayedFractionByGroup(double xs, int group, int position=CURRENT);
  
  void setClock(Clock* clock);

  /* Getter functions */
  int getId() const;

  int getNumEnergyGroups() const;
  int getNumDelayedGroups() const;
  double getEnergyPerFission();
  
  double* getSigmaT();
  double* getSigmaA();
  double* getSigmaS();
  double* getSigmaF();
  double* getNuSigmaF();
  double* getChi();
  double* getDifCoef();
  double* getVelocity();
  double* getPrecursorConc();
  double* getDecayConstant();
  double* getDelayedFraction();

  virtual double getSigmaTByGroup(int group, int position=CURRENT, double temp=0.0);
  virtual double getSigmaAByGroup(int group, int position=CURRENT, double temp=0.0);
  virtual double getSigmaSByGroup(int group_from, int group_to, int position=CURRENT, double temp=0.0);
  virtual double getSigmaFByGroup(int group, int position=CURRENT, double temp=0.0);
  virtual double getNuSigmaFByGroup(int group, int position=CURRENT, double temp=0.0);
  virtual double getChiByGroup(int group, int position=CURRENT, double temp=0.0);
  virtual double getDifCoefByGroup(int group, int position=CURRENT, double temp=0.0);
  virtual double getVelocityByGroup(int group, int position=CURRENT, double temp=0.0);
  virtual double getPrecursorConcByGroup(int group, int position=CURRENT, double temp=0.0);
  virtual double getDecayConstantByGroup(int group, int position=CURRENT, double temp=0.0);
  virtual double getDelayedFractionByGroup(int group, int position=CURRENT, double temp=0.0);
  virtual double getDelayedFractionTotal(int position=CURRENT, double temp=0.0);
  double getTemperatureConversionFactor();
  
  bool isFissionable();

  /* Worker functions */
  virtual std::string toString();
  void printString();
  virtual Material* clone();
  virtual void copy(int position_from, int position_to);
};

#endif /* MATERIAL_H_ */
