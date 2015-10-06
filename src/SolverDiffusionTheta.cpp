#include "SolverDiffusionTheta.h"

SolverDiffusionTheta::SolverDiffusionTheta(GeometryDiffusion* geometry) :
  SolverDiffusion(geometry){}


SolverDiffusionTheta::~SolverDiffusionTheta(){
}


void SolverDiffusionTheta::generateInitialMatrices(){

  int nx = _geometry_diffusion->getNumXShape();
  int ny = _geometry_diffusion->getNumYShape();
  int nz = _geometry_diffusion->getNumZShape();
  int ng = _num_energy_groups;
  double width  = _geometry->getWidth()  / nx;
  double height = _geometry->getHeight() / ny;
  double depth  = _geometry->getDepth()  / nz;
  double volume = width * height * depth;
  double val;

  _shape_A->clear();
  _shape_M->clear();

  for (int z=0; z < nz; z++){
    for (int y=0; y < ny; y++){
      for (int x=0; x < nx; x++){
        int cell = z*nx*ny + y*nx + x;
        
        double temp        = getTemperatureByValue(cell, CURRENT);
        Material* material = _geometry->getMaterial(cell);

        for (int e=0; e < ng; e++){
          
          double sig_a      = material->getSigmaAByGroup(e, CURRENT, temp);
          double dif_coef   = material->getDifCoefByGroup(e, CURRENT, temp);
          double chi        = material->getChiByGroup(e, CURRENT, temp);
          double sig_s, nu_sig_f;
          
          /* Absorption term on the diagonal */
          val = sig_a * volume;
          _shape_A->incrementValue(cell, e, cell, e, val);
          
          /* Buckling term on the diagonal */
          val = dif_coef * volume * _buckling;
          _shape_A->incrementValue(cell, e, cell, e, val);
          
          /* Outscattering term on diagonal */
          for (int g=0; g < ng; g++){
            if (e != g){
              sig_s = material->getSigmaSByGroup(e, g, CURRENT, temp);
              val = sig_s * volume;
              _shape_A->incrementValue(cell, e, cell, e, val);
            }
          }

          /* Fission terms */
          for (int g=0; g < ng; g++){
            nu_sig_f = material->getNuSigmaFByGroup(g, CURRENT, temp);
            val = chi * nu_sig_f * volume;
            _shape_M->incrementValue(cell, g, cell, e, val);
          }
          
          /* Inscattering term on off diagonals */
          for (int g=0; g < ng; g++){
            if (e != g){
              sig_s = material->getSigmaSByGroup(g, e, CURRENT, temp);
              val = - sig_s * volume;
              _shape_A->incrementValue(cell, g, cell, e, val);
            }
          }

          /* X_MIN SURFACE */
          
          /* Transport term on diagonal */
          dif_coef = getSurfaceDifCoef(cell, e, SURFACE_X_MIN, CURRENT);
          val = dif_coef * height * depth;
          _shape_A->incrementValue(cell, e, cell, e, val);
          
          /* Transport term on off diagonals */
          if (x != 0){
            dif_coef = getSurfaceDifCoef(cell, e, SURFACE_X_MIN, CURRENT);
            val = - dif_coef * height * depth;
            _shape_A->incrementValue(cell-1, e, cell, e, val);
          }          

          /* X_MAX SURFACE */
          
          /* Transport term on diagonal */
          dif_coef = getSurfaceDifCoef(cell, e, SURFACE_X_MAX, CURRENT);
          val = dif_coef * height * depth;
          _shape_A->incrementValue(cell, e, cell, e, val);
          
          /* Transport term on off diagonals */
          if (x != nx - 1){
            dif_coef = getSurfaceDifCoef(cell, e, SURFACE_X_MAX, CURRENT);
            val = - dif_coef * height * depth;
            _shape_A->incrementValue(cell+1, e, cell, e, val);
          }
          
          /* Y_MIN SURFACE */
          
          /* Transport term on diagonal */
          dif_coef = getSurfaceDifCoef(cell, e, SURFACE_Y_MIN, CURRENT);
          val = dif_coef * width * depth;
          _shape_A->incrementValue(cell, e, cell, e, val);

          /* Transport term on off diagonals */
          if (y != 0){
            dif_coef = getSurfaceDifCoef(cell, e, SURFACE_Y_MIN, CURRENT);
            val = - dif_coef * width * depth;
            _shape_A->incrementValue(cell-nx, e, cell, e, val);
          }          
          
          /* Y_MAX SURFACE */
          
          /* Transport term on diagonal */
          dif_coef = getSurfaceDifCoef(cell, e, SURFACE_Y_MAX, CURRENT);
          val = dif_coef * width * depth;
          _shape_A->incrementValue(cell, e, cell, e, val);
          
          /* Transport term on off diagonals */
          if (y != ny - 1){
            dif_coef = getSurfaceDifCoef(cell, e, SURFACE_Y_MAX, CURRENT);
            val = - dif_coef * width * depth;
            _shape_A->incrementValue(cell+nx, e, cell, e, val);
          }          

          /* Z_MIN SURFACE */
          
          /* Transport term on diagonal */
          dif_coef = getSurfaceDifCoef(cell, e, SURFACE_Z_MIN, CURRENT);
          val = dif_coef * width * height;
          _shape_A->incrementValue(cell, e, cell, e, val);
          
          /* Transport term on off diagonals */
          if (z != 0){
            dif_coef = getSurfaceDifCoef(cell, e, SURFACE_Z_MIN, CURRENT);
            val = - dif_coef * width * height;
            _shape_A->incrementValue(cell-nx*ny, e, cell, e, val);
          } 
          
          /* Z_MAX SURFACE */
          
          /* Transport term on diagonal */
          dif_coef = getSurfaceDifCoef(cell, e, SURFACE_Z_MAX, CURRENT);
          val = dif_coef * width * height;
          _shape_A->incrementValue(cell, e, cell, e, val);
          
          /* Transport term on off diagonals */
          if (z != nz - 1){
            dif_coef = getSurfaceDifCoef(cell, e, SURFACE_Z_MAX, CURRENT);
            val = - dif_coef * width * height;
            _shape_A->incrementValue(cell+nx*ny, e, cell, e, val);
          }
        }
      }
    }    
  }
}


double SolverDiffusionTheta::getSurfaceDifCoef(int cell, int group,
                                               int surface, int state){

  int nx = _geometry_diffusion->getNumXShape();
  int ny = _geometry_diffusion->getNumYShape();
  int nz = _geometry_diffusion->getNumZShape();
  double width  = _geometry->getWidth()  / nx;
  double height = _geometry->getHeight() / ny;
  double depth  = _geometry->getDepth()  / nz;

  int x = cell % (nx*ny) % nx;
  int y = cell % (nx*ny) / nx;
  int z = cell / (nx*ny);

  double length, length_perpen;
  double dif_linear, dif_coef_next, dif_coef;
  Material* material = _geometry->getMaterial(cell);
  Material* material_next;
  
  int cell_next = _geometry_diffusion->getNeighborShapeCell(x, y, z, surface);
    
  if (surface == SURFACE_X_MIN || surface == SURFACE_X_MAX){
    length = height * depth;
    length_perpen = width;
  }
  else if (surface == SURFACE_Y_MIN || surface == SURFACE_Y_MAX){
    length = width * depth;
    length_perpen = height;
  }
  else{
    length = width * height;
    length_perpen = depth;
  }
  
  dif_coef = material->getDifCoefByGroup(group, state);
            
  if (cell_next == -1){
    dif_linear = 2 * dif_coef / length_perpen / 
      (1 + 4 * dif_coef / length_perpen);
    dif_linear *= (1 - _geometry->getBoundary(surface));
  }
  else{
    material_next = _geometry->getMaterial(cell_next);
    dif_coef_next = material_next->getDifCoefByGroup(group, state);
    dif_linear = 2 * dif_coef * dif_coef_next / (length_perpen * dif_coef +
                                                 length_perpen * dif_coef_next);
  }

  return dif_linear;
}


void SolverDiffusionTheta::generateMatrices(){

  int state = FORWARD_OUT;
  int state_prev = PREVIOUS_OUT;
  int nx = _geometry_diffusion->getNumXShape();
  int ny = _geometry_diffusion->getNumYShape();
  int nz = _geometry_diffusion->getNumZShape();
  int ng = _num_energy_groups;
  double width  = _geometry->getWidth()  / nx;
  double height = _geometry->getHeight() / ny;
  double depth  = _geometry->getDepth()  / nz;
  double val;
  double dt = _clock->getTime(state) - _clock->getTime(state_prev);
  double volume = width * height * depth;
  
  _shape_b->clear();
  _shape_AM->clear();
  _shape_M->clear();
  _shape_A->clear();

  #pragma omp parallel for private(val) collapse(3)
  for (int z=0; z < nz; z++){
    for (int y=0; y < ny; y++){
      for (int x=0; x < nx; x++){

        int cell           = z*nx*ny + y*nx + x;
        Material* material = _geometry->getMaterial(cell);
        double temp        = getTemperatureByValue(cell, state);
        double volume      = _geometry->getVolume(cell);

        for (int e=0; e < ng; e++){

          double v           = material->getVelocityByGroup(e, state, temp);
          double flux_prev   = getFluxByValue(cell, e, state);

          double chi         = material->getChiByGroup(e, state, temp);
          double sig_a       = material->getSigmaAByGroup(e, state, temp); 
          double dif_coef    = material->getDifCoefByGroup(e, state, temp);
          double delay_tot   = material->getDelayedFractionTotal(state);
          double sig_s;
          double nu_sig_f;

          /* Time absorption term on the diagonal */
          val = 1.0 / v * volume / dt;
          _shape_AM->incrementValue(cell, e, cell, e, val);
          
          val = flux_prev / v * volume / dt;
          _shape_b->incrementValue(cell, e, val);
            
          /* Delayed neutron precursors */
          if (material->isFissionable()){
            for (int d=0; d < _num_delayed_groups; d++){
              double decay     = material->getDecayConstantByGroup(d);
              double conc      = material->getPrecursorConcByGroup(d, state);
              
              val = chi * decay * conc * volume;
              _shape_b->incrementValue(cell, e, val);
            }
          }
          
          /* Absorption term on the diagonal */
          val = sig_a * volume;
          _shape_AM->incrementValue(cell, e, cell, e, val);

          /* Buckling term on the diagonal */
          val = dif_coef * volume * _buckling;
          _shape_AM->incrementValue(cell, e, cell, e, val);
          
          /* Outscattering term on diagonal */
          for (int g=0; g < ng; g++){
            if (e != g){
              sig_s = material->getSigmaSByGroup(e, g, state, temp);
              val = sig_s * volume;
              _shape_AM->incrementValue(cell, e, cell, e, val);
            }
          }

          /* Fission terms */
          for (int g=0; g < ng; g++){
            nu_sig_f = material->getNuSigmaFByGroup(g, state, temp);
            
            val = - (1-delay_tot) * chi * nu_sig_f / _k_eff_0 * volume;
            _shape_AM->incrementValue(cell, g, cell, e, val);
            _shape_M->incrementValue(cell, g, cell, e, -val);
          }

          /* Inscattering term on off diagonals */
          for (int g=0; g < ng; g++){
            if (e != g){
              sig_s = material->getSigmaSByGroup(g, e, state, temp);
              val = - sig_s * volume;
              _shape_AM->incrementValue(cell, g, cell, e, val);
            }
          }

          /* X_MIN SURFACE */
          
          /* Transport term on diagonal */
          val = getSurfaceDifCoef(cell, e, SURFACE_X_MIN, state)
            * height * depth;
          _shape_AM->incrementValue(cell, e, cell, e, val);
          
          /* Transport term on off diagonals */
          if (x != 0){
            val = - getSurfaceDifCoef(cell, e, SURFACE_X_MIN, state)
              * height * depth;
            _shape_AM->incrementValue(cell-1, e, cell, e, val);
          }

          /* X_MAX SURFACE */
          
          /* Transport term on diagonal */
          val = getSurfaceDifCoef(cell, e, SURFACE_X_MAX, state)
            * height * depth;
          _shape_AM->incrementValue(cell, e, cell, e, val);
          
          /* Transport term on off diagonals */
          if (x != nx - 1){
            val = - getSurfaceDifCoef(cell, e, SURFACE_X_MAX, state)
              * height * depth;
            _shape_AM->incrementValue(cell+1, e, cell, e, val);
          }
          
          /* Y_MIN SURFACE */
          
          /* Transport term on diagonal */
          val = getSurfaceDifCoef(cell, e, SURFACE_Y_MIN, state)
            * width * depth;
          _shape_AM->incrementValue(cell, e, cell, e, val);
                    
          /* Transport term on off diagonals */
          if (y != 0){
            val = - getSurfaceDifCoef(cell, e, SURFACE_Y_MIN, state)
              * width * depth;
            _shape_AM->incrementValue(cell-nx, e, cell, e, val);
          }
          
          /* Y_MAX SURFACE */
          
          /* Transport term on diagonal */
          val = getSurfaceDifCoef(cell, e, SURFACE_Y_MAX, state)
            * width * depth;
          _shape_AM->incrementValue(cell, e, cell, e, val);
          
          /* Transport term on off diagonals */
          if (y != ny - 1){
            val = - getSurfaceDifCoef(cell, e, SURFACE_Y_MAX, state)
              * width * depth;
            _shape_AM->incrementValue(cell+nx, e, cell, e, val);
          }

          /* Z_MIN SURFACE */
          
          /* Transport term on diagonal */
          val = getSurfaceDifCoef(cell, e, SURFACE_Z_MIN, state)
            * width * height;
          _shape_AM->incrementValue(cell, e, cell, e, val);
          
          /* Transport term on off diagonals */
          if (z != 0){
            val = - getSurfaceDifCoef(cell, e, SURFACE_Z_MIN, state)
              * width * height;
            _shape_AM->incrementValue(cell-nx*ny, e, cell, e, val);
          }
          
          /* Z_MAX SURFACE */
          
          /* Transport term on diagonal */
          val = getSurfaceDifCoef(cell, e, SURFACE_Z_MAX, state)
            * width * height;
          _shape_AM->incrementValue(cell, e, cell, e, val);
          
          /* Transport term on off diagonals */
          if (z != nz - 1){
            val = - getSurfaceDifCoef(cell, e, SURFACE_Z_MAX, state)
              * width * height;
            _shape_AM->incrementValue(cell+nx*ny, e, cell, e, val);
          }
        }
      }
    }
  }
}


void SolverDiffusionTheta::computeInitialFlux(double tol){

  initializeClock();

  generateInitialMatrices();

  /* Solve the initial eigenvalue problem */
  _k_eff_0 = eigenvalueSolve(_shape_A, _shape_M, _flux[CURRENT], tol);

  log_printf(NORMAL, "Initial Eigenvalue: %.6f", _k_eff_0);

  /* Normalize flux to initial power level */
  normalizeFlux(_initial_power, CURRENT);

  /* Compute the initial precursor concentrations */
  computePrecursorConcentrations(CURRENT);

  /* Broadcast all field variables to all time steps */
  copyVariablesToAll(CURRENT);
}


void SolverDiffusionTheta::takeStep(){

  double residual = 1.e10;
  double tolerance = 1.e-2;
  
  /* Take step forward with clock */
  _clock->takeOuterStep();

  /* Integrate precursor concentrations */
  integratePrecursorConcentrations(PREVIOUS_OUT, FORWARD_OUT);

  /* Integrate the temperature to the forward time step */
  integrateTemperature(PREVIOUS_OUT, FORWARD_OUT);

  while(true){

    /* Generate the matrices */
    generateMatrices();

    /* Solve the linear problem to get the flux at the forward time step */
    linearSolve(_shape_AM, _shape_M, _flux[FORWARD_OUT], _shape_b, _flux_solve_tolerance);
    
    /* Reintegrate precursor concentrations */
    integratePrecursorConcentrations(PREVIOUS_OUT, FORWARD_OUT);
    
    /* Reintegrate the temperature to the forward time step */
    integrateTemperature(PREVIOUS_OUT, FORWARD_OUT);

    /* Check for convergence of the forward power */
    residual = computePowerRMSError(FORWARD_OUT, FORWARD_OUT_OLD);

    log_printf(NORMAL, "TIME = %1.4f, POWER = %.6e, RESIDUAL = %.6e",
               _clock->getTime(FORWARD_OUT), computeAveragePower(FORWARD_OUT), residual);
    
    if (residual < tolerance)
      break;
    else{
      copyVariables(FORWARD_OUT, FORWARD_OUT_OLD);
    }
  }

  /* Broadcast all field variables to all time steps */
  copyVariablesToAll(FORWARD_OUT);
}


void SolverDiffusionTheta::setStepSize(double time) {
  _clock->setDtOuter(time);
}
