#include "SolverDiffusion.h"

SolverDiffusion::SolverDiffusion(GeometryDiffusion* geometry) : Solver(geometry){

  /* Initialize variables */
  _geometry_diffusion = geometry;

  int num_x = _geometry_diffusion->getNumXShape();
  int num_y = _geometry_diffusion->getNumYShape();
  int num_z = _geometry_diffusion->getNumZShape();
  
  /* Initialize coarse mesh matrices */
  _shape_A_matrix  = new Matrix(num_x, num_y, num_z, _num_energy_groups);
  _shape_M_matrix  = new Matrix(num_x, num_y, num_z, _num_energy_groups);
  _shape_AM_matrix = new Matrix(num_x, num_y, num_z, _num_energy_groups);
  _shape_source    = new Vector(num_x, num_y, num_z, _num_energy_groups);

  /* Allocate memory for field variables */
  for (int t=0; t < NUM_STATES; t++){
    Vector* dif_linear  = new Vector(num_x, num_y, num_z, _num_energy_groups*NUM_SURFACES);
    Vector* temperature = new Vector(num_x, num_y, num_z, 1);
    Vector* flux        = new Vector(num_x, num_y, num_z, _num_energy_groups);
    Vector* shape       = new Vector(num_x, num_y, num_z, _num_energy_groups);
    Vector* power       = new Vector(num_x, num_y, num_z, 1);
    
    dif_linear->setAll(0.0);
    temperature->setAll(300.0);
    flux->setAll(1.0);
    shape->setAll(1.0);
    power->setAll(1.0);

    _dif_linear_fine[t] = dif_linear;
    _temperature[t] = temperature;
    _flux[t] = flux;
    _shape[t] = shape;
    _power[t] = power;
  }
}


SolverDiffusion::~SolverDiffusion(){
}


Matrix* SolverDiffusion::getShapeAMatrix(){
  return _shape_A_matrix;
}


Matrix* SolverDiffusion::getShapeMMatrix(){
  return _shape_M_matrix;
}


Matrix* SolverDiffusion::getShapeAMMatrix(){
  return _shape_AM_matrix;
}


Vector* SolverDiffusion::getShapeSource(){
  return _shape_source;
}

void SolverDiffusion::generateAdiabaticShapeMatrices(){

  int nx = _geometry_diffusion->getNumXShape();
  int ny = _geometry_diffusion->getNumYShape();
  int nz = _geometry_diffusion->getNumZShape();
  int ng = _num_energy_groups;
  double width  = _geometry->getWidth()  / nx;
  double height = _geometry->getHeight() / ny;
  double depth  = _geometry->getDepth()  / nz;
  double volume = width * height * depth;
  double val;

  _shape_A_matrix->clear();
  _shape_M_matrix->clear();

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
          _shape_A_matrix->incrementValue(cell, e, cell, e, val);
          
          /* Buckling term on the diagonal */
          val = dif_coef * volume * _buckling;
          _shape_A_matrix->incrementValue(cell, e, cell, e, val);
          
          /* Outscattering term on diagonal */
          for (int g=0; g < ng; g++){
            if (e != g){
              sig_s = material->getSigmaSByGroup(e, g, CURRENT, temp);
              val = sig_s * volume;
              _shape_A_matrix->incrementValue(cell, e, cell, e, val);
            }
          }

          /* Fission terms */
          for (int g=0; g < ng; g++){
            nu_sig_f = material->getNuSigmaFByGroup(g, CURRENT, temp);
            val = chi * nu_sig_f * volume;
            _shape_M_matrix->incrementValue(cell, g, cell, e, val);
          }
          
          /* Inscattering term on off diagonals */
          for (int g=0; g < ng; g++){
            if (e != g){
              sig_s = material->getSigmaSByGroup(g, e, CURRENT, temp);
              val = - sig_s * volume;
              _shape_A_matrix->incrementValue(cell, g, cell, e, val);
            }
          }

          /* X_MIN SURFACE */
          
          /* Transport term on diagonal */
          dif_coef = getDifLinearFineByValue(cell, e, SURFACE_X_MIN, CURRENT);
          val = dif_coef * height * depth;
          _shape_A_matrix->incrementValue(cell, e, cell, e, val);
          
          /* Transport term on off diagonals */
          if (x != 0){
            dif_coef = getDifLinearFineByValue(cell, e, SURFACE_X_MIN, CURRENT);
            val = - dif_coef * height * depth;
            _shape_A_matrix->incrementValue(cell-1, e, cell, e, val);
          }          

          /* X_MAX SURFACE */
          
          /* Transport term on diagonal */
          dif_coef = getDifLinearFineByValue(cell, e, SURFACE_X_MAX, CURRENT);
          val = dif_coef * height * depth;
          _shape_A_matrix->incrementValue(cell, e, cell, e, val);
          
          /* Transport term on off diagonals */
          if (x != nx - 1){
            dif_coef = getDifLinearFineByValue(cell, e, SURFACE_X_MAX, CURRENT);
            val = - dif_coef * height * depth;
            _shape_A_matrix->incrementValue(cell+1, e, cell, e, val);
          }
          
          /* Y_MIN SURFACE */
          
          /* Transport term on diagonal */
          dif_coef = getDifLinearFineByValue(cell, e, SURFACE_Y_MIN, CURRENT);
          val = dif_coef * width * depth;
          _shape_A_matrix->incrementValue(cell, e, cell, e, val);

          /* Transport term on off diagonals */
          if (y != 0){
            dif_coef = getDifLinearFineByValue(cell, e, SURFACE_Y_MIN, CURRENT);
            val = - dif_coef * width * depth;
            _shape_A_matrix->incrementValue(cell-nx, e, cell, e, val);
          }          
          
          /* Y_MAX SURFACE */
          
          /* Transport term on diagonal */
          dif_coef = getDifLinearFineByValue(cell, e, SURFACE_Y_MAX, CURRENT);
          val = dif_coef * width * depth;
          _shape_A_matrix->incrementValue(cell, e, cell, e, val);
          
          /* Transport term on off diagonals */
          if (y != ny - 1){
            dif_coef = getDifLinearFineByValue(cell, e, SURFACE_Y_MAX, CURRENT);
            val = - dif_coef * width * depth;
            _shape_A_matrix->incrementValue(cell+nx, e, cell, e, val);
          }          

          /* Z_MIN SURFACE */
          
          /* Transport term on diagonal */
          dif_coef = getDifLinearFineByValue(cell, e, SURFACE_Z_MIN, CURRENT);
          val = dif_coef * width * height;
          _shape_A_matrix->incrementValue(cell, e, cell, e, val);
          
          /* Transport term on off diagonals */
          if (z != 0){
            dif_coef = getDifLinearFineByValue(cell, e, SURFACE_Z_MIN, CURRENT);
            val = - dif_coef * width * height;
            _shape_A_matrix->incrementValue(cell-nx*ny, e, cell, e, val);
          } 
          
          /* Z_MAX SURFACE */
          
          /* Transport term on diagonal */
          dif_coef = getDifLinearFineByValue(cell, e, SURFACE_Z_MAX, CURRENT);
          val = dif_coef * width * height;
          _shape_A_matrix->incrementValue(cell, e, cell, e, val);
          
          /* Transport term on off diagonals */
          if (z != nz - 1){
            dif_coef = getDifLinearFineByValue(cell, e, SURFACE_Z_MAX, CURRENT);
            val = - dif_coef * width * height;
            _shape_A_matrix->incrementValue(cell+nx*ny, e, cell, e, val);
          }
        }
      }
    }    
  }
}


void SolverDiffusion::computeDiffusionCoefficientsFine(int state){

  int nx = _geometry_diffusion->getNumXShape();
  int ny = _geometry_diffusion->getNumYShape();
  int nz = _geometry_diffusion->getNumZShape();
  double width  = _geometry->getWidth()  / nx;
  double height = _geometry->getHeight() / ny;
  double depth  = _geometry->getDepth()  / nz;

  /* Compute the linear and nonlinear diffusion coefficients */
  #pragma omp parallel for collapse(3)
  for (int z=0; z < nz; z++){
    for (int y=0; y < ny; y++){
      for (int x=0; x < nx; x++){
        
        int cell = z*nx*ny + y*nx + x;

        double length, length_perpen;
        double dif_linear, dif_coef_next, dif_coef;
        Material* material = _geometry->getMaterial(cell);
        Material* material_next;
        
        for (int s=0; s < NUM_SURFACES; s++){
          
          int cell_next = _geometry_diffusion->getNeighborShapeCell(x, y, z, s);
          
          if (s == SURFACE_X_MIN || s == SURFACE_X_MAX){
            length = height * depth;
            length_perpen = width;
          }
          else if (s == SURFACE_Y_MIN || s == SURFACE_Y_MAX){
            length = width * depth;
            length_perpen = height;
          }
          else{
           length = width * height;
           length_perpen = depth;
          }
          
          for (int g=0; g < _num_energy_groups; g++){
            
            dif_coef = material->getDifCoefByGroup(g, state);
            
            if (cell_next == -1){
              dif_linear = 2 * dif_coef / length_perpen / 
                (1 + 4 * dif_coef / length_perpen);
              dif_linear *= (1 - _geometry->getBoundary(s));
            }
            else{
              material_next = _geometry->getMaterial(cell_next);
              dif_coef_next = material_next->getDifCoefByGroup(g, state);
              dif_linear = 2 * dif_coef * dif_coef_next / (length_perpen * dif_coef +
                                                           length_perpen * dif_coef_next);
            }
            
            setDifLinearFineByValue(dif_linear, cell, g, s, state);
          }
        }
      }
    }
  }
}


Vector* SolverDiffusion::getDifLinearFine(int state){
  return _dif_linear_fine[state];
}


double SolverDiffusion::getDifLinearFineByValue(int cell, int group, int side, int state){
  return _dif_linear_fine[state]->getValue(cell, side*_num_energy_groups+group);
}


void SolverDiffusion::setDifLinearFineByValue(double value, int cell, int group, 
                                              int side, int state){
  _dif_linear_fine[state]->setValue(cell, side*_num_energy_groups+group, value);
}


void SolverDiffusion::generateShapeMatrices(int state, int state_prev){

  int nx = _geometry_diffusion->getNumXShape();
  int ny = _geometry_diffusion->getNumYShape();
  int nz = _geometry_diffusion->getNumZShape();
  int ng = _num_energy_groups;
  double width  = _geometry->getWidth()  / nx;
  double height = _geometry->getHeight() / ny;
  double depth  = _geometry->getDepth()  / nz;
  double val;
  double dt = _clock->getTime(state) - _clock->getTime(state_prev);
  int* shape_to_amp = _geometry->getShapeToAmpMap();
  double volume = width * height * depth;
  
  _shape_source->clear();
  _shape_AM_matrix->clear();
  _shape_M_matrix->clear();

  #pragma omp parallel for private(val) collapse(3)
  for (int z=0; z < nz; z++){
    for (int y=0; y < ny; y++){
      for (int x=0; x < nx; x++){

        int cell_shape     = z*nx*ny + y*nx + x;
        int cell_amp       = shape_to_amp[cell_shape];
        Material* material = _geometry->getMaterial(cell_shape);
        double temp        = getTemperatureByValue(cell_shape, state);
        double volume      = _geometry->getVolume(cell_shape);

        for (int e=0; e < ng; e++){

          double v           = material->getVelocityByGroup(e, state, temp);
          double amp         = getAmplitudeByValue(cell_amp, e, state);
          double freq        = getFrequencyByValue(cell_amp, e, state);
          double flux_prev   = getFluxByValue(cell_shape, e, state_prev);
          double shape_prev  = getShapeByValue(cell_shape, e, state_prev);

          double chi         = material->getChiByGroup(e, state, temp);
          double sig_a       = material->getSigmaAByGroup(e, state, temp); 
          double dif_coef    = material->getDifCoefByGroup(e, state, temp);
          double delay_tot   = material->getDelayedFractionTotal(state);
          double sig_s;
          double nu_sig_f;

          /* Time absorption term on the diagonal */

          if (_method == IQS){
            val = freq / v * volume;
            _shape_AM_matrix->incrementValue(cell_shape, e, cell_shape, e, val);

            val = amp / v * volume / dt;
            _shape_AM_matrix->incrementValue(cell_shape, e, cell_shape, e, val);
            
            val = shape_prev * amp / v * volume / dt;
            _shape_source->incrementValue(cell_shape, e, val);
          }
          else if (_method == THETA) {
            val = freq / v * volume;
            _shape_AM_matrix->incrementValue(cell_shape, e, cell_shape, e, val);

            val = amp / v * volume / dt;
            _shape_AM_matrix->incrementValue(cell_shape, e, cell_shape, e, val);
            
            val = flux_prev / v * volume / dt;
            _shape_source->incrementValue(cell_shape, e, val);

          }
            
          /* Delayed neutron precursors */
          if (material->isFissionable()){
            for (int d=0; d < _num_delayed_groups; d++){
              double decay     = material->getDecayConstantByGroup(d);
              double conc      = material->getPrecursorConcByGroup(d, state);
              
              val = chi * decay * conc * volume;
              _shape_source->incrementValue(cell_shape, e, val);
            }
          }
          
          /* Absorption term on the diagonal */
          val = sig_a * volume * amp;
          _shape_AM_matrix->incrementValue(cell_shape, e, cell_shape, e, val);

          /* Buckling term on the diagonal */
          val = dif_coef * volume * _buckling * amp;
          _shape_AM_matrix->incrementValue(cell_shape, e, cell_shape, e, val);
          
          /* Outscattering term on diagonal */
          for (int g=0; g < ng; g++){
            if (e != g){
              sig_s = material->getSigmaSByGroup(e, g, state, temp);
              val = sig_s * volume * amp;
              _shape_AM_matrix->incrementValue(cell_shape, e, cell_shape, e, val);
            }
          }

          /* Fission terms */
          for (int g=0; g < ng; g++){
            nu_sig_f = material->getNuSigmaFByGroup(g, state, temp);
            
            val = - (1-delay_tot) * chi * nu_sig_f / _k_eff_0 * volume *
              getAmplitudeByValue(cell_amp, g, state);
            _shape_AM_matrix->incrementValue(cell_shape, g, cell_shape, e, val);
            _shape_M_matrix->incrementValue(cell_shape, g, cell_shape, e, -val);
          }

          /* Inscattering term on off diagonals */
          for (int g=0; g < ng; g++){
            if (e != g){
              sig_s = material->getSigmaSByGroup(g, e, state, temp);              
              val = - sig_s * volume * getAmplitudeByValue(cell_amp, g, state);
              _shape_AM_matrix->incrementValue(cell_shape, g, cell_shape, e, val);
            }
          }

          /* X_MIN SURFACE */
          
          /* Transport term on diagonal */
          val = getDifLinearFineByValue(cell_shape, e, SURFACE_X_MIN, state)
            * height * depth;
          _shape_AM_matrix->incrementValue(cell_shape, e, cell_shape, e, val);
          
          /* Transport term on off diagonals */
          if (x != 0){
            val = - getDifLinearFineByValue(cell_shape, e, SURFACE_X_MIN, state)
              * height * depth;
            _shape_AM_matrix->incrementValue(cell_shape-1, e, cell_shape, e, val);
          }

          /* X_MAX SURFACE */
          
          /* Transport term on diagonal */
          val = getDifLinearFineByValue(cell_shape, e, SURFACE_X_MAX, state)
            * height * depth;
          _shape_AM_matrix->incrementValue(cell_shape, e, cell_shape, e, val);
          
          /* Transport term on off diagonals */
          if (x != nx - 1){
            val = - getDifLinearFineByValue(cell_shape, e, SURFACE_X_MAX, state)
              * height * depth;
            _shape_AM_matrix->incrementValue(cell_shape+1, e, cell_shape, e, val);
          }
          
          /* Y_MIN SURFACE */
          
          /* Transport term on diagonal */
          val = getDifLinearFineByValue(cell_shape, e, SURFACE_Y_MIN, state)
            * width * depth;
          _shape_AM_matrix->incrementValue(cell_shape, e, cell_shape, e, val);
                    
          /* Transport term on off diagonals */
          if (y != 0){
            val = - getDifLinearFineByValue(cell_shape, e, SURFACE_Y_MIN, state)
              * width * depth;
            _shape_AM_matrix->incrementValue(cell_shape-nx, e, cell_shape, e, val);
          }
          
          /* Y_MAX SURFACE */
          
          /* Transport term on diagonal */
          val = getDifLinearFineByValue(cell_shape, e, SURFACE_Y_MAX, state)
            * width * depth;
          _shape_AM_matrix->incrementValue(cell_shape, e, cell_shape, e, val);
          
          /* Transport term on off diagonals */
          if (y != ny - 1){
            val = - getDifLinearFineByValue(cell_shape, e, SURFACE_Y_MAX, state)
              * width * depth;
            _shape_AM_matrix->incrementValue(cell_shape+nx, e, cell_shape, e, val);
          }

          /* Z_MIN SURFACE */
          
          /* Transport term on diagonal */
          val = getDifLinearFineByValue(cell_shape, e, SURFACE_Z_MIN, state)
            * width * height;
          _shape_AM_matrix->incrementValue(cell_shape, e, cell_shape, e, val);
          
          /* Transport term on off diagonals */
          if (z != 0){
            val = - getDifLinearFineByValue(cell_shape, e, SURFACE_Z_MIN, state)
              * width * height;
            _shape_AM_matrix->incrementValue(cell_shape-nx*ny, e, cell_shape, e, val);
          }
          
          /* Z_MAX SURFACE */
          
          /* Transport term on diagonal */
          val = getDifLinearFineByValue(cell_shape, e, SURFACE_Z_MAX, state)
            * width * height;
          _shape_AM_matrix->incrementValue(cell_shape, e, cell_shape, e, val);
          
          /* Transport term on off diagonals */
          if (z != nz - 1){
            val = - getDifLinearFineByValue(cell_shape, e, SURFACE_Z_MAX, state)
              * width * height;
            _shape_AM_matrix->incrementValue(cell_shape+nx*ny, e, cell_shape, e, val);
          }
        }
      }
    }
  }
}


void SolverDiffusion::generateAmpCurrent(int state){

  int nxa = _geometry->getNumXAmp();
  int nya = _geometry->getNumYAmp();
  int nza = _geometry->getNumZAmp();
  int nxs = _geometry_diffusion->getNumXShape();
  int nys = _geometry_diffusion->getNumYShape();
  int nzs = _geometry_diffusion->getNumZShape();
  double width  = _geometry->getWidth()  / nxs;
  double height = _geometry->getHeight() / nys;
  double depth  = _geometry->getDepth()  / nzs;
  double area;
  
  std::vector< std::vector<int> >* amp_to_shape = _geometry->getAmpToShapeMap();
  
  for (int z=0; z < nza; z++){
    for (int y=0; y < nya; y++){
      for (int x=0; x < nxa; x++){
        int cell_amp = z*nxa*nya + y*nxa + x;
        for (int e=0; e < _num_energy_groups; e++){
          for (int s=0; s < NUM_SURFACES; s++){

            double current = 0.0;

            std::vector<int>::iterator iter;
            for (iter = amp_to_shape->at(cell_amp).begin(); 
                 iter != amp_to_shape->at(cell_amp).end(); ++iter){

              int zs = (*iter)  / (nxs*nys);
              int ys = ((*iter) - zs*nxs*nys) / nxs;
              int xs = ((*iter) - zs*nxs*nys) % nxs;

              int neighbor_cell = _geometry_diffusion->getNeighborShapeCell(xs,ys,zs,s);
              double sense;
              
              if (s == SURFACE_X_MIN ||
                  s == SURFACE_Y_MIN ||
                  s == SURFACE_Z_MIN)
                sense = -1.0;
              else
                sense = 1.0;

              if (s == SURFACE_X_MIN || s == SURFACE_X_MAX)
                area = height * depth;
              else if (s == SURFACE_Y_MIN || s == SURFACE_Y_MAX)
                area = width * depth;
              else
                area = width * height;
              
              /* Check if shape surface coincides with amp surface */
              if ((s == SURFACE_X_MIN && xs % (nxs/nxa) == 0) || 
                  (s == SURFACE_Y_MIN && ys % (nys/nya) == 0) ||
                  (s == SURFACE_Z_MIN && zs % (nzs/nza) == 0) ||
                  (s == SURFACE_X_MAX && xs % (nxs/nxa) == nxs/nxa - 1) ||
                  (s == SURFACE_Y_MAX && ys % (nys/nya) == nys/nya - 1) ||
                  (s == SURFACE_Z_MAX && zs % (nzs/nza) == nzs/nza - 1)){
                if (neighbor_cell == -1)
                  current += sense * getDifLinearFineByValue(*iter, e, s, state) * 
                    area * getFluxByValue(*iter, e, state);
                else{
                  current += sense * getDifLinearFineByValue(*iter, e, s, state) * 
                    area * (getFluxByValue(*iter, e, state) -
                                      getFluxByValue(neighbor_cell, e, state));
                }
              }
            }

            setCurrentByValue(current, cell_amp, e, s, state);
          }
        }
      }
    }
  }
}


void SolverDiffusion::computeInitialShape(double tol){

  initializeClock();

  /* Generate the initial shape matrices */
  computeDiffusionCoefficientsFine(CURRENT);

  generateAdiabaticShapeMatrices();

  /* Solve the initial eigenvalue problem */
  _k_eff_0 = eigenvalueSolve(_shape_A_matrix, _shape_M_matrix, _flux[CURRENT], tol);

  log_printf(NORMAL, "Initial Eigenvalue: %.6f", _k_eff_0);

  /* Normalize flux to initial power level */
  normalizeFlux();

  /* Generate the amp mesh currents */
  generateAmpCurrent(CURRENT);

  /* Compute initial shape */
  computeShape(CURRENT, CURRENT, CURRENT);

  /* Compute the initial precursor concentrations */
  computeInitialPrecursorConcentrations();

  /* Broadcast all field variables to all time steps */
  broadcastToAll(CURRENT);
}


void SolverDiffusion::takeOuterStepOnly(){

  double tolerance = 1.e-2;
  double lin_solve_tol = 1.e-6;
  double residual = 1.e10;

  /* Take step forward with clock */
  _clock->takeOuterStep();

  /* Integrate precursor concentrations */
  integratePrecursorConcentrations(CURRENT, FORWARD_OUT);

  /* Integrate the temperature to the forward time step */
  integrateTemperature(CURRENT, FORWARD_OUT);
  
  computeDiffusionCoefficientsFine(PREVIOUS_OUT);

  while(true){

    /* Generate the shape matrices */
    computeDiffusionCoefficientsFine(FORWARD_OUT);
    generateShapeMatrices(FORWARD_OUT, PREVIOUS_OUT);

    /* Solve the linear problem to get the flux at the forward time step */
    linearSolve(_shape_AM_matrix, _shape_M_matrix, _flux[FORWARD_OUT],
                _shape_source, lin_solve_tol);
    
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
      copyFieldVariables(FORWARD_OUT, FORWARD_OUT_OLD);
    }
  }

  /* Broadcast all field variables to all time steps */
  broadcastToAll(FORWARD_OUT);
}


void SolverDiffusion::takeInnerStep(){
  return;
}


void SolverDiffusion::takeOuterStep(){

  double tolerance = 1.e-2;
  double lin_solve_tol = 1.e-6;
  double residual = 1.e10;

  /* Take step forward with clock */
  _clock->takeOuterStep();

  while(true){

    /* Prolongate flux */
    while (_clock->getTime(CURRENT) < _clock->getTime(FORWARD_OUT) - 1.e-6)
      takeInnerStep();
    
    /* Reconstruct flux at FORWARD_OUT */

    /* Generate the shape matrices */
    computeDiffusionCoefficientsFine(FORWARD_OUT);
    generateShapeMatrices(FORWARD_OUT, PREVIOUS_OUT);

    /* Solve the linear problem to get the shape at the forward time step */
    linearSolve(_shape_AM_matrix, _shape_M_matrix, _shape[FORWARD_OUT],
                _shape_source, lin_solve_tol);
    
    /* Reconstruct flux at FORWARD_OUT */

    /* Check for convergence of the forward power */
    residual = computePowerRMSError(FORWARD_OUT, FORWARD_OUT_OLD);

    log_printf(NORMAL, "TIME = %1.4f, POWER = %.6e, RESIDUAL = %.6e",
               _clock->getTime(FORWARD_OUT), computeAveragePower(FORWARD_OUT), residual);
    
    if (residual < tolerance)
      break;
    else{
      _clock->resetToPreviousOuterStep();
      copyFieldVariables(FORWARD_OUT, FORWARD_OUT_OLD);
    }
  }

  /* Broadcast all field variables to all time steps */
  broadcastToAll(FORWARD_OUT);

  return;
}


GeometryDiffusion* SolverDiffusion::getGeometryDiffusion(){
  return _geometry_diffusion;
}
