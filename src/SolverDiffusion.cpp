#include "SolverDiffusion.h"

SolverDiffusion::SolverDiffusion(GeometryDiffusion* geometry) : Solver(geometry){

  /* Initialize variables */
  _geometry_diffusion = geometry;

  int num_x = _geometry_diffusion->getNumXShape();
  int num_y = _geometry_diffusion->getNumYShape();
  int num_z = _geometry_diffusion->getNumZShape();
  
  /* Initialize coarse mesh matrices */
  _shape_A  = new Matrix(num_x, num_y, num_z, _num_energy_groups);
  _shape_M  = new Matrix(num_x, num_y, num_z, _num_energy_groups);
  _shape_AM = new Matrix(num_x, num_y, num_z, _num_energy_groups);
  _shape_b  = new Vector(num_x, num_y, num_z, _num_energy_groups);

  /* Allocate memory for field variables */
  for (int t=0; t < NUM_STATES; t++){
    Vector* temperature = new Vector(num_x, num_y, num_z, 1);
    Vector* flux        = new Vector(num_x, num_y, num_z, _num_energy_groups);
    Vector* weight      = new Vector(num_x, num_y, num_z, _num_energy_groups);
    Vector* power       = new Vector(num_x, num_y, num_z, 1);
    
    temperature->setAll(300.0);
    flux->setAll(1.0);
    weight->setAll(1.0);
    power->setAll(1.0);

    _temperature[t] = temperature;
    _flux[t]   = flux;
    _weight[t] = weight;
    _power[t]  = power;
  }
}


SolverDiffusion::~SolverDiffusion(){
}


Matrix* SolverDiffusion::getShapeA(){
  return _shape_A;
}


Matrix* SolverDiffusion::getShapeM(){
  return _shape_M;
}


Matrix* SolverDiffusion::getShapeAM(){
  return _shape_AM;
}


Vector* SolverDiffusion::getShapeb(){
  return _shape_b;
}


GeometryDiffusion* SolverDiffusion::getGeometryDiffusion(){
  return _geometry_diffusion;
}
