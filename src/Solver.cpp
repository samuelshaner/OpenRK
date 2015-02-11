#include "Solver.h"

Solver::Solver(StructuredShapeMesh* shape_mesh, AmpMesh* amp_mesh){

  /* Initialize variables */
  _AM_shape = NULL;
  _A_shape = NULL;
  _M_shape = NULL;
  _b_shape = NULL;
  _AM_amp = NULL;
  _A_amp = NULL;
  _M_amp = NULL;
  _b_amp = NULL;
  _k_eff = 0.0;
  
  /* Set solver parameters */
  setNumThreads(1);
  _shape_mesh = shape_mesh;
  _amp_mesh = amp_mesh;

  int nc = _shape_mesh->getNumX() * _shape_mesh->getNumY();
  int ng = _shape_mesh->getNumShapeEnergyGroups();
  _AM_shape = new double*[nc];
  _A_shape = new double*[nc];
  _M_shape = new double*[nc];
  _b_shape = new double[nc*ng];

  for (int i=0; i < nc; i++){
    _AM_shape[i] = new double[ng*(ng+4)];
    _A_shape[i] = new double[ng*(ng+4)];
    _M_shape[i] = new double[ng*ng];
  }

  if (_amp_mesh != NULL){
    nc = _amp_mesh->getNumX() * _amp_mesh->getNumY();
    ng = _amp_mesh->getNumAmpEnergyGroups();

    _AM_amp = new double*[nc];
    _A_amp = new double*[nc];
    _M_amp = new double*[nc];
    _b_amp = new double[nc*ng];
    
    for (int i=0; i < nc; i++){
      _AM_amp[i] = new double[ng*(ng+4)];
      _A_amp[i] = new double[ng*(ng+4)];
      _M_amp[i] = new double[ng*ng];
    }
  }
}


Solver::~Solver(){

  if (_AM_shape != NULL){
    for (int i=0; i < _shape_mesh->getNumX() * _shape_mesh->getNumY(); i++)
      delete [] _AM_shape[i];

    delete [] _AM_shape;
  }    

  if (_A_shape != NULL){
    for (int i=0; i < _shape_mesh->getNumX() * _shape_mesh->getNumY(); i++)
      delete [] _A_shape[i];

    delete [] _A_shape;
  }    

  if (_M_shape != NULL){
    for (int i=0; i < _shape_mesh->getNumX() * _shape_mesh->getNumY(); i++)
      delete [] _M_shape[i];

    delete [] _M_shape;
  }    

  if (_b_shape != NULL){
    delete [] _b_shape;
  }    

  if (_AM_amp != NULL){
    for (int i=0; i < _amp_mesh->getNumX() * _amp_mesh->getNumY(); i++)
      delete [] _AM_amp[i];

    delete [] _AM_amp;
  }    

  if (_A_amp != NULL){
    for (int i=0; i < _amp_mesh->getNumX() * _amp_mesh->getNumY(); i++)
      delete [] _A_amp[i];

    delete [] _A_amp;
  }    

  if (_M_amp != NULL){
    for (int i=0; i < _amp_mesh->getNumX() * _amp_mesh->getNumY(); i++)
      delete [] _M_amp[i];

    delete [] _M_amp;
  }    

  if (_b_amp != NULL){
    delete [] _b_amp;
  }    
}


double** Solver::getAMShape(){
  return _AM_shape;
}


double** Solver::getAShape(){
  return _A_shape;
}


double** Solver::getMShape(){
  return _M_shape;
}


double* Solver::getBShape(){
  return _b_shape;
}


double** Solver::getAMAmp(){
  return _AM_amp;
}


double** Solver::getAAmp(){
  return _A_amp;
}


double** Solver::getMAmp(){
  return _M_amp;
}


double* Solver::getBAmp(){
  return _b_amp;
}


void Solver::makeAMShapeInitial(int position){

  int nx = _shape_mesh->getNumX();
  int ny = _shape_mesh->getNumY();
  int ng = _shape_mesh->getNumShapeEnergyGroups();
  double width = _shape_mesh->getCellWidth();
  double height = _shape_mesh->getCellHeight();
  double volume = width * height;
  double* temps = _shape_mesh->getTemperature(position);

  #pragma omp parallel for
  for (int i = 0; i < nx*ny; i++){
    for (int g=0; g < ng*ng; g++)
      _M_shape[i][g] = 0.0;

    for (int g=0; g < ng*(ng+4); g++){
      _A_shape[i][g] = 0.0;
      _AM_shape[i][g] = 0.0;
    }
  }

  #pragma omp parallel for
  for (int y=0; y < ny; y++){
    for (int x=0; x < nx; x++){
      int cell = y*nx+x;
      double temp = temps[cell];

      for (int e=0; e < ng; e++){
        Material* material = _shape_mesh->getMaterial(cell);

        /* Absorption term on the diagonal */
        _A_shape[cell][e * (ng+4) + e + 2] += material->getSigmaAByGroup(e, position, temp) * volume;

        /* Buckling term on the diagonal */
        _A_shape[cell][e * (ng+4) + e + 2] += material->getDifCoefByGroup(e, position, temp) * volume *
          _shape_mesh->getBuckling();

        /* Outscattering term on diagonal */
        for (int g=0; g < ng; g++){
          if (e != g)
            _A_shape[cell][e * (ng+4) + e + 2] += material->getSigmaSByGroup(e, g, position, temp) * volume;
        }
        
        /* Fission terms */
        for (int g=0; g < ng; g++){
          _M_shape[cell][e * ng + g] += material->getChiByGroup(e, position, temp) *
            material->getNuSigmaFByGroup(g, position, temp) * volume;
        }

        /* Inscattering term on off diagonals */
        for (int g=0; g < ng; g++){
          if (e != g)
            _A_shape[cell][e * (ng+4) + g + 2] -= material->getSigmaSByGroup(g, e, position, temp) * volume;
        }

        /* RIGHT SURFACE */

        /* Transport term on diagonal */
        _A_shape[cell][e * (ng+4) + e + 2] += _shape_mesh->getDifLinearByValue(cell, e, 2, position) * height;

        /* Transport term on off diagonals */
        if (x != nx - 1)
          _A_shape[cell][e * (ng+4) + ng + 2] -= _shape_mesh->getDifLinearByValue(cell, e, 2, position) * height;

        /* LEFT SURFACE */

        /* Transport term on diagonal */
        _A_shape[cell][e * (ng+4) + e + 2] += _shape_mesh->getDifLinearByValue(cell, e, 0, position) * height;

        /* Transport term on off diagonals */
        if (x != 0)
          _A_shape[cell][e * (ng+4)] -= _shape_mesh->getDifLinearByValue(cell, e, 0, position) * height;

        /* BOTTOM SURFACE */

        /* Transport term on diagonal */
        _A_shape[cell][e * (ng+4) + e + 2] += _shape_mesh->getDifLinearByValue(cell, e, 1, position) * width;

        /* Transport term on off diagonals */
        if (y != 0)
          _A_shape[cell][e * (ng+4) + 1] -= _shape_mesh->getDifLinearByValue(cell, e, 1, position) * width;

        /* LEFT SURFACE */

        /* Transport term on diagonal */
        _A_shape[cell][e * (ng+4) + e + 2] += _shape_mesh->getDifLinearByValue(cell, e, 3, position) * width;

        /* Transport term on off diagonals */
        if (y != ny - 1)
          _A_shape[cell][e * (ng+4) + ng + 3] -= _shape_mesh->getDifLinearByValue(cell, e, 3, position) * width;
      }
    }
  }    
}


void Solver::computeInitialShape(double tol){

  _shape_mesh->computeDifCoefs(CURRENT);
  makeAMShapeInitial(CURRENT);

  int nx = _shape_mesh->getNumX();
  int ny = _shape_mesh->getNumY();
  int ng = _shape_mesh->getNumShapeEnergyGroups();
  double* old_source = new double[nx*ny*ng];
  double* flux_temp = new double[nx*ny*ng];

  _k_eff = eigenvalueSolve2d(_A_shape, nx*ny, ng*(ng+4), _M_shape, nx*ny, ng*ng, _shape_mesh->getFlux(CURRENT), nx*ny*ng, _b_shape, nx*ny*ng, old_source, nx*ny*ng, flux_temp, nx*ny*ng, ng, nx, ny, tol);

  log_printf(NORMAL, "done solving eigenvalue problem");
  
  _shape_mesh->setKeff0(_k_eff);
  _amp_mesh->setKeff0(_k_eff);
  
  delete [] old_source;
  delete [] flux_temp;
}


void Solver::makeAMAmp(double wt){

  int nx = _amp_mesh->getNumX();
  int ny = _amp_mesh->getNumY();
  int ng = _amp_mesh->getNumAmpEnergyGroups();
  double width = _amp_mesh->getCellWidth();
  double height = _amp_mesh->getCellHeight();
  double volume = width * height;
  double* flux_prev = _amp_mesh->getFlux(PREVIOUS_IN);
  double dt = _amp_mesh->getClock()->getDtInner();
  double* temps = _amp_mesh->getTemperature(CURRENT);
  double* temps_prev = _amp_mesh->getTemperature(PREVIOUS_IN);

  #pragma omp parallel for
  for (int i = 0; i < nx*ny; i++){
    for (int g=0; g < ng*ng; g++)
      _M_amp[i][g] = 0.0;

    for (int g=0; g < ng*(ng+4); g++){
      _A_amp[i][g] = 0.0;
      _AM_amp[i][g] = 0.0;
    }

    for (int g=0; g < ng; g++)
      _b_amp[i*ng+g] = 0.0;
  }

  double beta = 0.0;
  for (int d=0; d < _amp_mesh->getNumDelayedGroups(); d++)
    beta += _amp_mesh->getDelayedFractionByGroup(d);

  #pragma omp parallel for
  for (int y=0; y < ny; y++){
    for (int x=0; x < nx; x++){
      int cell = y*nx+x;
      double temp = temps[cell];
      double temp_prev = temps_prev[cell];
      Material* material = _amp_mesh->getMaterial(cell);
      
      for (int e=0; e < ng; e++){

        /* Time absorption term on the diagonal */
        _AM_amp[cell][e * (ng+4) + e + 2] += 1.0 / dt / material->getVelocityByGroup(e, CURRENT, temp) * volume;
        _b_amp[cell*ng + e] += flux_prev[cell*ng + e] / dt /
          material->getVelocityByGroup(e, PREVIOUS_IN, temp_prev) * volume;

        /* Delayed neutron precursors */
        for (int d=0; d < _amp_mesh->getNumDelayedGroups(); d++){
          _b_amp[cell*ng + e] += wt * material->getChiByGroup(e, CURRENT, temp)
            * _amp_mesh->getDecayConstantByGroup(d) * material->getPrecursorConcByGroup(d, CURRENT) * volume;
          _b_amp[cell*ng + e] += (1-wt) * material->getChiByGroup(e, PREVIOUS_IN, temp_prev)
            * _amp_mesh->getDecayConstantByGroup(d) * material->getPrecursorConcByGroup(d, PREVIOUS_IN) * volume;
        }
        
        /* Absorption term on the diagonal */
        _AM_amp[cell][e * (ng+4) + e + 2] += wt * material->getSigmaAByGroup(e, CURRENT, temp) * volume;
        _b_amp[cell*ng + e] -= (1-wt) * material->getSigmaAByGroup(e, PREVIOUS_IN, temp_prev)
          * flux_prev[cell*ng + e] * volume;

        /* Buckling term on the diagonal */
        _AM_amp[cell][e * (ng+4) + e + 2] += wt * material->getDifCoefByGroup(e, CURRENT, temp) * volume *
          _amp_mesh->getBuckling();
        _b_amp[cell*ng + e] -= (1-wt) * material->getDifCoefByGroup(e, PREVIOUS_IN, temp_prev)
          * _amp_mesh->getBuckling() * flux_prev[cell*ng + e] * volume;

        /* Outscattering term on diagonal */
        for (int g=0; g < ng; g++){
          if (e != g){
            _AM_amp[cell][e * (ng+4) + e + 2] += wt * material->getSigmaSByGroup(e, g, CURRENT, temp) * volume;
            _b_amp[cell*ng+e] -= (1-wt) * material->getSigmaSByGroup(e, g, PREVIOUS_IN, temp_prev)
              * flux_prev[cell*ng + e] * volume;
          }
        }

        /* Fission terms */
        for (int g=0; g < ng; g++){
          _AM_amp[cell][e * (ng+4) + g + 2] -= wt * (1-beta) * material->getChiByGroup(e, CURRENT, temp) *
            material->getNuSigmaFByGroup(g, CURRENT, temp) / _amp_mesh->getKeff0() * volume;
          _b_amp[cell*ng+e] += (1-wt) * (1-beta) * material->getChiByGroup(e, PREVIOUS_IN, temp_prev) *
            material->getNuSigmaFByGroup(g, PREVIOUS_IN, temp_prev) / _amp_mesh->getKeff0()
            * flux_prev[cell*ng+g] * volume;
        }

        /* Inscattering term on off diagonals */
        for (int g=0; g < ng; g++){
          if (e != g){
            _AM_amp[cell][e * (ng+4) + g + 2] -= wt * material->getSigmaSByGroup(g, e, CURRENT, temp) * volume;
            _b_amp[cell*ng+e] += (1-wt) * material->getSigmaSByGroup(g, e, PREVIOUS_IN, temp_prev)
              * flux_prev[cell*ng+g] * volume;
          }
        }

        /* RIGHT SURFACE */

        /* Transport term on diagonal */
        _AM_amp[cell][e * (ng+4) + e + 2] += wt *
          (_amp_mesh->getDifLinearByValue(cell, e, 2, CURRENT) -
           _amp_mesh->getDifNonlinearByValue(cell, e, 2, CURRENT)) * height;
        _b_amp[cell*ng+e] -= (1-wt) *
          (_amp_mesh->getDifLinearByValue(cell, e, 2, PREVIOUS_IN) -
           _amp_mesh->getDifNonlinearByValue(cell, e, 2, PREVIOUS_IN)) * height * flux_prev[cell*ng+e];

        
        /* Transport term on off diagonals */
        if (x != nx - 1){
          _AM_amp[cell][e * (ng+4) + ng + 2] -= wt *
            (_amp_mesh->getDifLinearByValue(cell, e, 2, CURRENT) +
             _amp_mesh->getDifNonlinearByValue(cell, e, 2, CURRENT)) * height;
          _b_amp[cell*ng+e] += (1-wt) *
            (_amp_mesh->getDifLinearByValue(cell, e, 2, PREVIOUS_IN) +
             _amp_mesh->getDifNonlinearByValue(cell, e, 2, PREVIOUS_IN)) * height * flux_prev[(cell+1)*ng+e];
        }          

        /* LEFT SURFACE */

        /* Transport term on diagonal */
        _AM_amp[cell][e * (ng+4) + e + 2] += wt *
          (_amp_mesh->getDifLinearByValue(cell, e, 0, CURRENT) +
           _amp_mesh->getDifNonlinearByValue(cell, e, 0, CURRENT)) * height;
        _b_amp[cell*ng+e] -= (1-wt) *
          (_amp_mesh->getDifLinearByValue(cell, e, 0, PREVIOUS_IN) +
           _amp_mesh->getDifNonlinearByValue(cell, e, 0, PREVIOUS_IN)) * height * flux_prev[cell*ng+e];

        
        /* Transport term on off diagonals */
        if (x != 0){
          _AM_amp[cell][e * (ng+4)] -= wt *
            (_amp_mesh->getDifLinearByValue(cell, e, 0, CURRENT) -
             _amp_mesh->getDifNonlinearByValue(cell, e, 0, CURRENT)) * height;
          _b_amp[cell*ng+e] += (1-wt) *
            (_amp_mesh->getDifLinearByValue(cell, e, 0, PREVIOUS_IN) -
             _amp_mesh->getDifNonlinearByValue(cell, e, 0, PREVIOUS_IN)) * height * flux_prev[(cell-1)*ng+e];
        }          

        /* BOTTOM SURFACE */

        /* Transport term on diagonal */
        _AM_amp[cell][e * (ng+4) + e + 2] += wt *
          (_amp_mesh->getDifLinearByValue(cell, e, 1, CURRENT) +
           _amp_mesh->getDifNonlinearByValue(cell, e, 1, CURRENT)) * width;
        _b_amp[cell*ng+e] -= (1-wt) *
          (_amp_mesh->getDifLinearByValue(cell, e, 1, PREVIOUS_IN) +
           _amp_mesh->getDifNonlinearByValue(cell, e, 1, PREVIOUS_IN)) * width * flux_prev[cell*ng+e];

        
        /* Transport term on off diagonals */
        if (y != 0){
          _AM_amp[cell][e * (ng+4) + 1] -= wt *
            (_amp_mesh->getDifLinearByValue(cell, e, 1, CURRENT) -
             _amp_mesh->getDifNonlinearByValue(cell, e, 1, CURRENT)) * width;
          _b_amp[cell*ng+e] += (1-wt) *
            (_amp_mesh->getDifLinearByValue(cell, e, 1, PREVIOUS_IN) -
             _amp_mesh->getDifNonlinearByValue(cell, e, 1, PREVIOUS_IN)) * width * flux_prev[(cell-nx)*ng+e];
        }          
        
        /* TOP SURFACE */

        /* Transport term on diagonal */
        _AM_amp[cell][e * (ng+4) + e + 2] += wt *
          (_amp_mesh->getDifLinearByValue(cell, e, 3, CURRENT) -
           _amp_mesh->getDifNonlinearByValue(cell, e, 3, CURRENT)) * width;
        _b_amp[cell*ng+e] -= (1-wt) *
          (_amp_mesh->getDifLinearByValue(cell, e, 3, PREVIOUS_IN) -
           _amp_mesh->getDifNonlinearByValue(cell, e, 3, PREVIOUS_IN)) * width * flux_prev[cell*ng+e];

        
        /* Transport term on off diagonals */
        if (y != ny - 1){
          _AM_amp[cell][e * (ng+4) + ng + 3] -= wt *
            (_amp_mesh->getDifLinearByValue(cell, e, 3, CURRENT) +
             _amp_mesh->getDifNonlinearByValue(cell, e, 3, CURRENT)) * width;
          _b_amp[cell*ng+e] += (1-wt) *
            (_amp_mesh->getDifLinearByValue(cell, e, 3, PREVIOUS_IN) +
             _amp_mesh->getDifNonlinearByValue(cell, e, 3, PREVIOUS_IN)) * width * flux_prev[(cell+nx)*ng+e];
        }
      }
    }
  }
}


void Solver::makeAMShape(double wt){

  int nx = _shape_mesh->getNumX();
  int ny = _shape_mesh->getNumY();
  int ng = _shape_mesh->getNumShapeEnergyGroups();
  double width = _shape_mesh->getCellWidth();
  double height = _shape_mesh->getCellHeight();
  double volume = width * height;
  double* flux_prev = _shape_mesh->getFlux(PREVIOUS_OUT);
  double dt = _shape_mesh->getClock()->getDtOuter();
  double* temps = _shape_mesh->getTemperature(FORWARD_OUT);
  double* temps_prev = _shape_mesh->getTemperature(PREVIOUS_OUT);

  #pragma omp parallel for
  for (int i = 0; i < nx*ny; i++){
    for (int g=0; g < ng*ng; g++)
      _M_shape[i][g] = 0.0;

    for (int g=0; g < ng*(ng+4); g++){
      _A_shape[i][g] = 0.0;
      _AM_shape[i][g] = 0.0;
    }

    for (int g=0; g < ng; g++)
      _b_shape[i*ng+g] = 0.0;
  }

  double beta = 0.0;
  for (int d=0; d < _shape_mesh->getNumDelayedGroups(); d++)
    beta += _shape_mesh->getDelayedFractionByGroup(d);

  #pragma omp parallel for
  for (int y=0; y < ny; y++){
    for (int x=0; x < nx; x++){
      int cell = y*nx+x;
      double temp = temps[cell];
      double temp_prev = temps_prev[cell];
      Material* material = _shape_mesh->getMaterial(cell);
      
      for (int e=0; e < ng; e++){

        /* Time absorption term on the diagonal */
        _AM_shape[cell][e * (ng+4) + e + 2] += 1.0 / dt / material->getVelocityByGroup(e, FORWARD_OUT, temp) * volume;
        _b_shape[cell*ng + e] += flux_prev[cell*ng + e] / dt /
          material->getVelocityByGroup(e, PREVIOUS_OUT, temp_prev) * volume;

        /* Delayed neutron precursors */
        if (material->isFissionable()){
          for (int d=0; d < _shape_mesh->getNumDelayedGroups(); d++){
            _b_shape[cell*ng + e] += wt * material->getChiByGroup(e, FORWARD_OUT, temp)
              * _shape_mesh->getDecayConstantByGroup(d) * material->getPrecursorConcByGroup(d, FORWARD_OUT) * volume;
            _b_shape[cell*ng + e] += (1-wt) * material->getChiByGroup(e, PREVIOUS_OUT, temp_prev)
              * _shape_mesh->getDecayConstantByGroup(d) * material->getPrecursorConcByGroup(d, PREVIOUS_OUT) * volume;
          }
        }
        
        /* Absorption term on the diagonal */
        _AM_shape[cell][e * (ng+4) + e + 2] += wt * material->getSigmaAByGroup(e, FORWARD_OUT, temp) * volume;
        _b_shape[cell*ng + e] -= (1-wt) * material->getSigmaAByGroup(e, PREVIOUS_OUT, temp_prev)
          * flux_prev[cell*ng + e] * volume;

        /* Buckling term on the diagonal */
        _AM_shape[cell][e * (ng+4) + e + 2] += wt * material->getDifCoefByGroup(e, FORWARD_OUT, temp) * volume *
          _shape_mesh->getBuckling();
        _b_shape[cell*ng + e] -= (1-wt) * material->getDifCoefByGroup(e, PREVIOUS_OUT, temp_prev)
          * _shape_mesh->getBuckling() * flux_prev[cell*ng + e] * volume;

        /* Outscattering term on diagonal */
        for (int g=0; g < ng; g++){
          if (e != g){
            _AM_shape[cell][e * (ng+4) + e + 2] += wt * material->getSigmaSByGroup(e, g, FORWARD_OUT, temp) * volume;
            _b_shape[cell*ng+e] -= (1-wt) * material->getSigmaSByGroup(e, g, PREVIOUS_OUT, temp_prev)
              * flux_prev[cell*ng + e] * volume;
          }
        }
        
        /* Fission terms */
        for (int g=0; g < ng; g++){
          _AM_shape[cell][e * (ng+4) + g + 2] -= wt * (1-beta) * material->getChiByGroup(e, FORWARD_OUT, temp) *
            material->getNuSigmaFByGroup(g, FORWARD_OUT, temp) / _shape_mesh->getKeff0() * volume;
          _b_shape[cell*ng+e] += (1-wt) * (1-beta) * material->getChiByGroup(e, PREVIOUS_OUT, temp_prev) *
            material->getNuSigmaFByGroup(g, PREVIOUS_OUT, temp_prev) / _shape_mesh->getKeff0()
            * flux_prev[cell*ng+g] * volume;
        }

        /* Inscattering term on off diagonals */
        for (int g=0; g < ng; g++){
          if (e != g){
            _AM_shape[cell][e * (ng+4) + g + 2] -= wt * material->getSigmaSByGroup(g, e, FORWARD_OUT, temp) * volume;
            _b_shape[cell*ng+e] += (1-wt) * material->getSigmaSByGroup(g, e, PREVIOUS_OUT, temp_prev)
              * flux_prev[cell*ng+g] * volume;
          }
        }

        /* RIGHT SURFACE */

        /* Transport term on diagonal */
        _AM_shape[cell][e * (ng+4) + e + 2] += wt *
          _shape_mesh->getDifLinearByValue(cell, e, 2, FORWARD_OUT) * height;
        _b_shape[cell*ng+e] -= (1-wt) *
          _shape_mesh->getDifLinearByValue(cell, e, 2, PREVIOUS_OUT) * height * flux_prev[cell*ng+e];
        
        /* Transport term on off diagonals */
        if (x != nx - 1){
          _AM_shape[cell][e * (ng+4) + ng + 2] -= wt *
            _shape_mesh->getDifLinearByValue(cell, e, 2, FORWARD_OUT) * height;
          _b_shape[cell*ng+e] += (1-wt) *
            _shape_mesh->getDifLinearByValue(cell, e, 2, PREVIOUS_OUT) * height * flux_prev[(cell+1)*ng+e];
        }

        /* LEFT SURFACE */

        /* Transport term on diagonal */
        _AM_shape[cell][e * (ng+4) + e + 2] += wt *
          _shape_mesh->getDifLinearByValue(cell, e, 0, FORWARD_OUT) * height;
        _b_shape[cell*ng+e] -= (1-wt) *
          _shape_mesh->getDifLinearByValue(cell, e, 0, PREVIOUS_OUT) * height * flux_prev[cell*ng+e];

        /* Transport term on off diagonals */
        if (x != 0){
          _AM_shape[cell][e * (ng+4)] -= wt *
            _shape_mesh->getDifLinearByValue(cell, e, 0, FORWARD_OUT) * height;
          _b_shape[cell*ng+e] += (1-wt) *
            _shape_mesh->getDifLinearByValue(cell, e, 0, PREVIOUS_OUT) * height * flux_prev[(cell-1)*ng+e];
        }          

        /* BOTTOM SURFACE */

        /* Transport term on diagonal */
        _AM_shape[cell][e * (ng+4) + e + 2] += wt *
          _shape_mesh->getDifLinearByValue(cell, e, 1, FORWARD_OUT) * width;
        _b_shape[cell*ng+e] -= (1-wt) *
          _shape_mesh->getDifLinearByValue(cell, e, 1, PREVIOUS_OUT) * width * flux_prev[cell*ng+e];

        
        /* Transport term on off diagonals */
        if (y != 0){
          _AM_shape[cell][e * (ng+4) + 1] -= wt *
            _shape_mesh->getDifLinearByValue(cell, e, 1, FORWARD_OUT) * width;
          _b_shape[cell*ng+e] += (1-wt) *
            _shape_mesh->getDifLinearByValue(cell, e, 1, PREVIOUS_OUT) * width * flux_prev[(cell-nx)*ng+e];
        }          
        
        /* TOP SURFACE */

        /* Transport term on diagonal */
        _AM_shape[cell][e * (ng+4) + e + 2] += wt *
          _shape_mesh->getDifLinearByValue(cell, e, 3, FORWARD_OUT) * width;
        _b_shape[cell*ng+e] -= (1-wt) *
          _shape_mesh->getDifLinearByValue(cell, e, 3, PREVIOUS_OUT) * width * flux_prev[cell*ng+e];

        
        /* Transport term on off diagonals */
        if (y != ny - 1){
          _AM_shape[cell][e * (ng+4) + ng + 3] -= wt *
            _shape_mesh->getDifLinearByValue(cell, e, 3, FORWARD_OUT) * width;
          _b_shape[cell*ng+e] += (1-wt) *
            _shape_mesh->getDifLinearByValue(cell, e, 3, PREVIOUS_OUT) * width * flux_prev[(cell+nx)*ng+e];
        }          
      }
    }
  }    
}


