/**
 * @file TransientMaterial.h
 * @brief
 * @date February 7, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef TRANSIENTMATERIAL_H_
#define TRANSIENTMATERIAL_H_

#ifdef __cplusplus
#include "Material.h
#endif


/**
 * @class TransientMaterial TransientMaterial.h "src/TransientMaterial.h"
 * @brief The Material class represents a unique material and its relevant
 *        nuclear data (i.e., multigroup cross-sections) for neutron transport.
 */
class TransientMaterial {

protected:


public:
  Material(int id);
  virtual ~Material();

  int getUid() const;
  int getId() const;
  int getNumEnergyGroups() const;
  double* getSigmaT();
  double* getSigmaA();
  double* getSigmaS();
  double* getSigmaF();
  double* getNuSigmaF();
  double* getChi();
  double* getDifCoef();
  double* getBuckling();
  double* getDifHat();
  double* getDifTilde();
  bool isFissionable();

  void setNumEnergyGroups(const int num_groups);

  void setSigmaT(double* xs, int num_groups);
  void setSigmaA(double* xs, int num_groups);
  void setSigmaS(double* xs, int num_groups);
  void setSigmaF(double* xs, int num_groups);
  void setNuSigmaF(double* xs, int num_groups);
  void setChi(double* xs, int num_groups);
  void setDifCoef(double* xs, int num_groups);

  void setSigmaTByGroup(double xs, int group);
  void setSigmaAByGroup(double xs, int group);
  void setSigmaFByGroup(double xs, int group);
  void setNuSigmaFByGroup(double xs, int group);
  void setSigmaSByGroup(double xs, int group_from, int group_to);
  void setChiByGroup(double xs, int group);
  void setDifCoefByGroup(double xs, int group);

  std::string toString();
  void printString();

  Material* clone();
};

#endif /* MATERIAL_H_ */
