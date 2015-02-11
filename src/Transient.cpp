#include "Transient.h"
Transient::Transient(transientMethod inner_method, transientMethod outer_method, double wt_inner, double wt_outer){

  /* Initialize variables */
  _inner_method = BACKWARD_EULER;
  _outer_method = BACKWARD_EULER;
  _inner_wt = 1.0;
  _outer_wt = 1.0;
  _clock = NULL;
  _shape_mesh = NULL;
  _amp_mesh = NULL;
  _solver = NULL;
    
  /* Set solver parameters */
  setInnerMethod(inner_method, wt_inner);
  setOuterMethod(outer_method, wt_outer);
}


Transient::~Transient(){
}


void Transient::setInnerMethod(transientMethod inner_method, double wt){

  if (inner_method == CUSTOM){
    if (wt < 0.0 || wt > 1.0)
      log_printf(ERROR, "Unable to use CUSTOM inner method with wt ouside "
                 "range of [0,1]. wt: %d", wt);
    else{
      _inner_method = inner_method;
      _inner_wt = wt;
    }
  }
  else if (inner_method == FORWARD_EULER){
    _inner_method = inner_method;
    _inner_wt = 0.0;
  }
  else if (inner_method == BACKWARD_EULER){
    _inner_method = inner_method;
    _inner_wt = 1.0;
  }
  else if (inner_method == CRANK_NICOLSON){
    _inner_method = inner_method;
    _inner_wt = 0.5;
  }  
}


void Transient::setOuterMethod(transientMethod outer_method, double wt){

  if (outer_method == CUSTOM){
    if (wt < 0.0 || wt > 1.0)
      log_printf(ERROR, "Unable to use CUSTOM outer method with wt ouside "
                 "range of [0,1]. wt: %d", wt);
    else{
      _outer_method = outer_method;
      _outer_wt = wt;
    }
  }
  else if (outer_method == FORWARD_EULER){
    _outer_method = outer_method;
    _outer_wt = 0.0;
  }
  else if (outer_method == BACKWARD_EULER){
    _outer_method = outer_method;
    _outer_wt = 1.0;
  }
  else if (outer_method == CRANK_NICOLSON){
    _outer_method = outer_method;
    _outer_wt = 0.5;
  }  
}


void Transient::setInitialPower(double initial_power){
  _initial_power = initial_power;
}


void Transient::setSolver(Solver* solver){
  _solver = solver;
}


void Transient::setShapeMesh(StructuredShapeMesh* mesh){
  _shape_mesh = mesh;
}


void Transient::setAmpMesh(AmpMesh* mesh){
  _amp_mesh = mesh;
}


void Transient::setClock(Clock* clock){
  _clock = clock;
}

void Transient::computeInitialShape(){

  _amp_mesh->setClock(_clock);
  _shape_mesh->setClock(_clock);

  _solver->computeInitialShape(1.e-10);

  double initial_power = _shape_mesh->computeAveragePower(CURRENT);
  _shape_mesh->scaleFlux(CURRENT, _initial_power / initial_power);

  _shape_mesh->computeInitialPrecursorConc(CURRENT);

  _amp_mesh->condenseMaterials(CURRENT, true);
  _amp_mesh->computeCurrent(CURRENT);
  _amp_mesh->computeDifCoefs(CURRENT);
  broadcastToAll(CURRENT);  
}


void Transient::takeInnerStep(){

  _clock->takeInnerStep();

  _amp_mesh->interpolateDifNonlinear(PREVIOUS_OUT, FORWARD_OUT, CURRENT);

  _shape_mesh->synthesizeFlux(CURRENT);

  _shape_mesh->integrateTemperature(PREVIOUS_IN, CURRENT);

  _shape_mesh->integratePrecursorConc(PREVIOUS_IN, CURRENT);

  _amp_mesh->condenseMaterials(CURRENT);

  int ng = _amp_mesh->getNumAmpEnergyGroups();
  int nx = _amp_mesh->getNumX();
  int ny = _amp_mesh->getNumY();
  double* flux_temp = new double[nx*ny*ng];
  for (int i=0; i < nx*ny*ng; i++)
    flux_temp[i] = 0.0;
  double tol = 1.e-4;

  while (true){

    _amp_mesh->copyFlux(CURRENT, FORWARD_IN_OLD);

    _solver->makeAMAmp(_inner_wt);
    linearSolve2d(_solver->getAMAmp(), nx*ny, ng*(ng+4), _amp_mesh->getFlux(CURRENT), nx*ny*ng,
                  _solver->getBAmp(), nx*ny*ng, flux_temp, nx*ny*ng, nx, ny, ng, 1.e-8);

    _shape_mesh->synthesizeFlux(CURRENT);

    _shape_mesh->integrateTemperature(PREVIOUS_IN, CURRENT);
    
    _shape_mesh->integratePrecursorConc(PREVIOUS_IN, CURRENT);
    
    _amp_mesh->condenseMaterials(CURRENT);

    double residual = _amp_mesh->computePowerL2Norm(CURRENT, FORWARD_IN_OLD);
    double power = _shape_mesh->computeAveragePower(CURRENT);

    log_printf(NORMAL, "TIME = %1.4f, POWER = %.6e, RESIDUAL = %.6e", _clock->getTime(CURRENT), power, residual);

    if (residual < tol)
      break;
  }

  broadcastToOne(CURRENT, PREVIOUS_IN);

  delete [] flux_temp;  
}


void Transient::takeOuterStep(){

  _clock->takeOuterStep();

  broadcastToActive(FORWARD_OUT);  
  int ng = _shape_mesh->getNumShapeEnergyGroups();
  int nx = _shape_mesh->getNumX();
  int ny = _shape_mesh->getNumY();
  double* flux_temp = new double[nx*ny*ng];
  double tol = 1.e-4;

  while (_clock->getTime(CURRENT) < _clock->getTime(FORWARD_OUT) - 1.e-6)
    takeInnerStep();
  
  _shape_mesh->reconstructFlux(CURRENT, FORWARD_OUT, CURRENT);
  broadcastToOne(CURRENT, FORWARD_OUT);
  broadcastToOne(PREVIOUS_OUT, CURRENT);
  broadcastToOne(PREVIOUS_OUT, PREVIOUS_IN);
    
  while (true){

    _shape_mesh->copyFlux(FORWARD_OUT, FORWARD_OUT_OLD);

    _shape_mesh->computeDifCoefs(FORWARD_OUT);

    _solver->makeAMShape(_outer_wt);
    linearSolve2d(_solver->getAMShape(), nx*ny, ng*(ng+4), _shape_mesh->getFlux(FORWARD_OUT), nx*ny*ng,
                  _solver->getBShape(), nx*ny*ng, flux_temp, nx*ny*ng, nx, ny, ng, 1.e-8);

    _amp_mesh->computeCurrent(FORWARD_OUT);
    _amp_mesh->condenseMaterials(FORWARD_OUT, true);
    _amp_mesh->computeDifCoefs(FORWARD_OUT);

    _clock->resetToPreviousOuterStep();

    while (_clock->getTime(CURRENT) < _clock->getTime(FORWARD_OUT) - 1.e-6)
      takeInnerStep();

    _shape_mesh->reconstructFlux(CURRENT, FORWARD_OUT, CURRENT);
    broadcastToOne(CURRENT, FORWARD_OUT);

    double residual = _shape_mesh->computePowerL2Norm(FORWARD_OUT, FORWARD_OUT_OLD);
    
    log_printf(NORMAL, "OUTER RESIDUAL = %.6e", residual);

    if (residual < tol)
      break;
    else{
      broadcastToOne(PREVIOUS_OUT, CURRENT);
      broadcastToOne(PREVIOUS_OUT, PREVIOUS_IN);
    }
  }

  delete [] flux_temp;  
}


void Transient::takeOuterStepOnly(){

  _clock->takeOuterStep();

  broadcastToActive(FORWARD_OUT);  
  int ng = _shape_mesh->getNumShapeEnergyGroups();
  int nx = _shape_mesh->getNumX();
  int ny = _shape_mesh->getNumY();
  double* flux_temp = new double[nx*ny*ng];
  double tol = 1.e-8;
    
  while (true){

    _shape_mesh->integrateTemperature(PREVIOUS_OUT, FORWARD_OUT);

    _shape_mesh->integratePrecursorConc(PREVIOUS_OUT, FORWARD_OUT);

    _shape_mesh->copyFlux(FORWARD_OUT, FORWARD_OUT_OLD);
    
    _shape_mesh->computeDifCoefs(FORWARD_OUT);

    _solver->makeAMShape(_outer_wt);
    linearSolve2d(_solver->getAMShape(), nx*ny, ng*(ng+4), _shape_mesh->getFlux(FORWARD_OUT), nx*ny*ng, 
                _solver->getBShape(), nx*ny*ng, flux_temp, nx*ny*ng, nx, ny, ng, 1.e-8);

    double residual = _shape_mesh->computePowerL2Norm(FORWARD_OUT, FORWARD_OUT_OLD);
    double power = _shape_mesh->computeAveragePower(FORWARD_OUT);
    log_printf(NORMAL, "TIME = %1.4f, POWER = %.6e, RESIDUAL = %.6e", _clock->getTime(FORWARD_OUT), power, residual);

    log_printf(NORMAL, "OUTER RESIDUAL = %.6e", residual);

    if (residual < tol){
      _clock->setTime(FORWARD_OUT, _clock->getTime(CURRENT));
      break;
    }
  }

  delete [] flux_temp;  
}


void Transient::broadcastToActive(clockPosition position){

  for (int c=1; c < 7; c++){

    _amp_mesh->copyFlux(position, c);
    _amp_mesh->copyCurrent(position, c);
    _amp_mesh->copyTemperature(position, c);
    _amp_mesh->copyDifLinear(position, c);
    _amp_mesh->copyDifNonlinear(position, c);

    _shape_mesh->copyFlux(position, c);
    _shape_mesh->copyCurrent(position, c);
    _shape_mesh->copyTemperature(position, c);
    _shape_mesh->copyDifLinear(position, c);

    for (int i=0; i < _amp_mesh->getNumX() * _amp_mesh->getNumY(); i++)
      _amp_mesh->getMaterial(i)->copy(position, c);

    for (int i=0; i < _shape_mesh->getNumX() * _shape_mesh->getNumY(); i++)
      _shape_mesh->getMaterial(i)->copy(position, c);
  }  
}


void Transient::broadcastToAll(clockPosition position){

  for (int c=0; c < 8; c++){

    _amp_mesh->copyFlux(position, c);
    _amp_mesh->copyCurrent(position, c);
    _amp_mesh->copyTemperature(position, c);
    _amp_mesh->copyDifLinear(position, c);
    _amp_mesh->copyDifNonlinear(position, c);

    _shape_mesh->copyFlux(position, c);
    _shape_mesh->copyCurrent(position, c);
    _shape_mesh->copyTemperature(position, c);
    _shape_mesh->copyDifLinear(position, c);

    for (int i=0; i < _amp_mesh->getNumX() * _amp_mesh->getNumY(); i++)
      _amp_mesh->getMaterial(i)->copy(position, c);

    for (int i=0; i < _shape_mesh->getNumX() * _shape_mesh->getNumY(); i++)
      _shape_mesh->getMaterial(i)->copy(position, c);
  }  
}


void Transient::broadcastToOne(clockPosition position_from, clockPosition position_to){

  _amp_mesh->copyFlux(position_from, position_to);
  _amp_mesh->copyCurrent(position_from, position_to);
  _amp_mesh->copyTemperature(position_from, position_to);
  _amp_mesh->copyDifLinear(position_from, position_to);
  _amp_mesh->copyDifNonlinear(position_from, position_to);
  
  _shape_mesh->copyFlux(position_from, position_to);
  _shape_mesh->copyCurrent(position_from, position_to);
  _shape_mesh->copyTemperature(position_from, position_to);
  _shape_mesh->copyDifLinear(position_from, position_to);
  
  for (int i=0; i < _amp_mesh->getNumX() * _amp_mesh->getNumY(); i++)
    _amp_mesh->getMaterial(i)->copy(position_from, position_to);
  
  for (int i=0; i < _shape_mesh->getNumX() * _shape_mesh->getNumY(); i++)
    _shape_mesh->getMaterial(i)->copy(position_from, position_to); 
}
