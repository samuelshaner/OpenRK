#include "SolverDiffusion.h"

SolverDiffusion::SolverDiffusion(GeometryDiffusion* geometry) : Solver(geometry){

  /* Initialize variables */
  _geometry_diffusion = geometry;

  /* Initialize coarse mesh matrices */
  _shape_A_matrix  = new double*[_num_shape_cells];
  _shape_M_matrix  = new double*[_num_shape_cells];
  _shape_AM_matrix = new double*[_num_shape_cells];
  _shape_source    = new double[_num_shape_cells*_num_energy_groups];
  
  for (int i=0; i < nc_amp; i++){
    _shape_A_matrix[i]  = new double[_num_energy_groups*(_num_energy_groups+6)];
    _shape_M_matrix[i]  = new double[_num_energy_groups*(_num_energy_groups)];
    _shape_AM_matrix[i] = new double[_num_energy_groups*(_num_energy_groups+6)];
  }

  /* Allocate memory for field variables */
  for (int t=0; t < 8; t++){
    double* dif_linear = new double[_num_shape_cells*_num_energy_groups*6];

    memset(dif_linear, 0.0, sizeof(double) * _num_shape_cells * _num_energy_groups * 6);

    _dif_linear_fine[t] = dif_linear;
  }
}


SolverDiffusion::~SolverDiffusion(){
}


double** SolverDiffusion::getShapeAMatrix(){
  return _shape_A_matrix;
}


double** SolverDiffusion::getShapeMMatrix(){
  return _shape_M_matrix;
}


double** SolverDiffusion::getShapeAMMatrix(){
  return _shape_AM_matrix;
}


double* SolverDiffusion::getShapeSource(){
  return _shape_source;
}

void SolverDiffusion::generateInitialShapeMatrices(){

  int nx = _geometry_diffusion->getNumXShape();
  int ny = _geometry_diffusion->getNumYShape();
  int nz = _geometry_diffusion->getNumZShape();
  int ng = _num_energy_groups;
  double width = _geometry->getWidth() / nx;
  double height = _geometry->getHeight() / ny;
  double depth = _geometry->getDepth() / nz;
  double volume = width * height * depth;

  #pragma omp parallel for
  for (int i = 0; i < nx*ny*nz; i++){
    for (int g=0; g < ng*ng; g++)
      _shape_M_matrix[i][g] = 0.0;

    for (int g=0; g < ng*(ng+6); g++){
      _shape_A_matrix[i][g] = 0.0;
    }
  }

  for (int z=0; z < nz; z++){
    #pragma omp parallel for
    for (int y=0; y < ny; y++){
      for (int x=0; x < nx; x++){
        int cell = z*nx*ny + y*nx+x;
        double temp = temps[cell];
        Material* material = _geometry->getMaterial(cell);
        
        for (int e=0; e < ng; e++){
          
          int diag = e * (ng+6) + e + 3;

          /* Absorption term on the diagonal */
          _shape_A_matrix[cell][diag] += 
            material->getSigmaAByGroup(e, CURRENT, temp) * volume;
          
          /* Buckling term on the diagonal */
          _shape_A_matrix[cell][diag] += 
            material->getDifCoefByGroup(e, CURRENT, temp) * volume * _buckling;
          
          /* Outscattering term on diagonal */
          for (int g=0; g < ng; g++){
            if (e != g){
              _shape_A_matrix[cell][diag] += 
                material->getSigmaSByGroup(e, g, CURRENT, temp) * volume;
            }
          }
          
          /* Fission terms */
          for (int g=0; g < ng; g++){
            _shape_M_matrix[cell][e * ng + g] += 
              material->getChiByGroup(e, CURRENT, temp) * 
              material->getNuSigmaFByGroup(g, CURRENT, temp) * volume;
          }
          
          /* Inscattering term on off diagonals */
          for (int g=0; g < ng; g++){
            if (e != g){
              _shape_A_matrix[cell][e * (ng+6) + g + 3] -= 
                material->getSigmaSByGroup(g, e, CURRENT, temp) * volume;
            }
          }
                    
          /* LEFT SURFACE */
          
          /* Transport term on diagonal */
          _shape_A_matrix[cell][diag] +=
            getDifLinearFineByValue(cell, e, 0, CURRENT) * height * depth;
          
          /* Transport term on off diagonals */
          if (x != 0){
            _shape_A_matrix[cell][e * (ng+6)] -= 
              getDifLinearFineByValue(cell, e, 0, CURRENT) * height * depth;
          }          

          /* RIGHT SURFACE */
          
          /* Transport term on diagonal */
          _shape_A_matrix[cell][diag] +=
            getDifLinearFineByValue(cell, e, 3, CURRENT) * height * depth;
          
          /* Transport term on off diagonals */
          if (x != nx - 1){
            _shape_A_matrix[cell][e * (ng+6) + ng + 3] -= 
              getDifLinearFineByValue(cell, e, 3, CURRENT) * height * depth;
          }
          
          /* BACK SURFACE */
          
          /* Transport term on diagonal */
          _shape_A_matrix[cell][diag] +=
            getDifLinearFineByValue(cell, e, 1, CURRENT) * width * depth;

          /* Transport term on off diagonals */
          if (y != 0){
            _shape_A_matrix[cell][e * (ng+6) + 1] -=
              getDifLinearFineByValue(cell, e, 1, CURRENT) * width * depth;
          }          
          
          /* FRONT SURFACE */
          
          /* Transport term on diagonal */
          _shape_A_matrix[cell][diag] +=
            getDifLinearFineByValue(cell, e, 4, CURRENT) * width * depth;
                    
          /* Transport term on off diagonals */
          if (y != ny - 1){
            _shape_A_matrix[cell][e * (ng+6) + ng + 4] -= 
              getDifLinearFineByValue(cell, e, 4, CURRENT) * width * depth;
          }          

          /* BOTTOM SURFACE */
          
          /* Transport term on diagonal */
          _shape_A_matrix[cell][diag] +=
            getDifLinearFineByValue(cell, e, 2, CURRENT) * width * height;
          
          /* Transport term on off diagonals */
          if (z != 0){
            _shape_A_matrix[cell][e * (ng+6) + 2] -=
              getDifLinearFineByValue(cell, e, 2, CURRENT) * width * height;
          }          
          
          /* TOP SURFACE */
          
          /* Transport term on diagonal */
          _shape_A_matrix[cell][diag] +=
            getDifLinearFineByValue(cell, e, 5, CURRENT) * width * height;
          
          /* Transport term on off diagonals */
          if (z != nz - 1){
            _shape_A_matrix[cell][e * (ng+6) + ng + 5] -= 
              getDifLinearFineByValue(cell, e, 5, CURRENT) * width * height;
          }
        }
      }
    }    
  }
}


void SolverDiffusion::computeDiffusionCoefficientsFine(int time){

  int nx = _geometry->getNumXShape();
  int ny = _geometry->getNumYShape();
  int nz = _geometry->getNumZShape();
  double width = _geometry->getWidth() / nx;
  double height = _geometry->getHeight() / ny;
  double depth = _geometry->getDepth() / nz;
  double* temps = getTemperature(time);
  std::vector< std::vector<int> > amp_to_shape = _geometry->getAmpToShapeMap();

  /* Compute the linear and nonlinear diffusion coefficients */
  for (int z=0; z < nz; z++){
#pragma omp parallel for
    for (int y=0; y < ny; y++){
      for (int x=0; x < nx; x++){
        
        int cell = z*nx*ny + y*nx + x;

        double length, length_perpen;
        double dif_linear, dif_coef_next, dif_coef;
        
        for (int s=0; s < 6; s++){
          
          int cell_next = _geometry_diffusion->getNeighborShapeCell(x, y, z, s);
          
          if (s == 0 || s == 3){
            length = height * depth;
            length_perpen = width;
          }
          else if (s == 1 || s == 4){
            length = width * depth;
            length_perpen = height;
          }
          else{
           length = width * height;
           length_perpen = depth;
          }
          
          for (int g=0; g < _num_energy_groups; g++){
            
            dif_coef = material->getDifCoefByValue(cell, g, time);
            
            if (cell_next == -1){
              dif_linear = 2 * dif_coef / length_perpen / 
                (1 + 4 * dif_coef / length_perpen);
              dif_linear *= _boundaries[s];
            }
            else{
              dif_coef_next = dif_coefs[cell_next*_num_energy_groups + g];
              dif_linear = 2 * dif_coef * dif_coef_next / (length_perpen * dif_coef +
                                                           length_perpen * dif_coef_next);
            }
            
            setDifLinearFineByValue(dif_linear, cell, g, s, time);
          }
        }
      }
    }
  }
}


double* SolverDiffusion::getDifLinearFine(int time){
  return _dif_linear_fine[time];
}


double SolverDiffusion::getDifLinearFineByValue(int cell, int group, int side, int time){
  return _dif_linear_fine[time][(cell*6+side) * _num_energy_groups+group];
}


void SolverDiffusion::setDifLinearFineByValue(double value, int cell, int group, 
                                              int side, int time){
  _dif_linear_fine[time][(cell*6+side) * _num_energy_groups+group] = value;
}


void SolverDiffusion::generateShapeMatrices(double wt){

  int nx = _geometry->getNumXShape();
  int ny = _geometry->getNumYShape();
  int nz = _geometry->getNumZShape();
  int ng = _num_energy_groups;
  double width = _geometry->getWidth() / nx;
  double height = _geometry->getHeight() / ny;
  double depth = _geometry->getDepth() / nz;

  double dt = _clock->getDtInner();

  double* temps = getTemperature(FORWARD_OUT);
  double* temps_prev = getTemperature(PREVIOUS_OUT);
  double* flux = getFlux(FORWARD_OUT);
  double* flux_prev = getShape(PREVIOUS_OUT);
  double* amp = getAmplitude(FORWARD_OUT);
  double* freq = getFrequency(FORWARD_OUT);

  int* shape_to_amp = _geometry->getAmpToShapeMap();
 
  #pragma omp parallel for
  for (int i = 0; i < _num_shape_cells; i++){
    for (int g=0; g < ng; g++)
      _shape_source[i*ng+g] = 0.0;
    
    for (int g=0; g < ng*(ng+6); g++)
      _shape_AM_matrix[i][g] = 0.0;
  }
  
  for (int z=0; z < nz; z++){
#pragma omp parallel for
    for (int y=0; y < ny; y++){
      for (int x=0; x < nx; x++){

        int cell_shape = z*nx*ny + y*nx+x;
        int cell_amp = shape_to_amp[cell_shape];
        Material* material = _geometry->getMaterial(cell_shape);
        double temp = temps[cell_shape];
        double temp_prev = temps_prev[cell_shape];
        double volume = _geometry->getVolume(cell_shape);
        
        for (int e=0; e < ng; e++){

          int diag = e * (ng+6) + e + 3;
          int row_amp = cell_amp*ng + e;
          int row_shape = cell_shape*ng + e;

          /* Time absorption term on the diagonal */
          if (method == IQS){
            _shape_AM_matrix[cell_shape][diag] += 1.0 / dt / 
              material->getVelocityByGroup(e, FORWARD_OUT, temp) * volume * flux[row_shape];

            _shape_AM_matrix[cell_shape][diag] += 1.0 / 
              material->getVelocityByGroup(e, FORWARD_OUT, temp)  / 
              amp[row_amp] * volume * freq[row_amp];
            
            _shape_source[row_shape] += flux_prev[row_shape] / dt / 
              material->getVelocityByGroup(e, PREVIOUS_OUT, temp_prev) * 
              volume * amp[row_amp] / amp_prev[row_amp];
          }
          else if (method == THETA){
            _shape_AM_matrix[cell_shape][diag] += 1.0 / dt / 
              material->getVelocityByGroup(e, FORWARD_OUT, temp) * volume;

            _shape_source[row_shape] += flux_prev[row_shape] / dt / 
              material->getVelocityByGroup(e, PREVIOUS_OUT, temp_prev) * volume;
          }
                    
            /* Delayed neutron precursors */
            for (int d=0; d < _num_delayed_groups; d++){
              _shape_source[row_shape] += wt * material->getChiByGroup(e, FORWARD_OUT, temp)
                * material->getDecayConstantByGroup(d) * 
                material->getPrecursorConcByGroup(d, FORWARD_OUT) * volume;
              
              _shape_source[row_shape] += (1-wt) * 
                material->getChiByGroup(e, PREVIOUS_OUT, temp_prev)
                * material->getDecayConstantByGroup(d) * 
                material->getPrecursorConcByGroup(d, PREVIOUS_OUT) * volume;
            }
          
            /* Absorption term on the diagonal */
            _shape_AM_matrix[cell_shape][diag] += wt * 
              material->getSigmaAByGroup(e, FORWARD_OUT, temp) * volume;

            _shape_source[row_shape] -= (1-wt) * 
              material->getSigmaAByGroup(e, PREVIOUS_OUT, temp_prev)
              * flux_prev[row_shape] * volume;
            
            /* Buckling term on the diagonal */
            _shape_AM_matrix[cell_shape][diag] += wt * 
              material->getDifCoefByGroup(e, FORWARD_OUT, temp) * volume * _buckling;
            
            _shape_source[row_shape] -= (1-wt) * 
              material->getDifCoefByGroup(e, PREVIOUS_OUT, temp_prev) * 
              flux_prev[row_shape] * volume * _buckling;
            
            /* Outscattering term on diagonal */
            for (int g=0; g < ng; g++){
              if (e != g){
                _shape_AM_matrix[cell_shape][diag] += wt * 
                  material->getSigmaSByGroup(e, g, FORWARD_OUT, temp) * volume;
                
                _shape_source[row_shape] -= (1-wt) * flux_prev[row_shape] * volume *
                  material->getSigmaSByGroup(e, g, PREVIOUS_OUT, temp_prev);
              }
            }
            
            /* Fission terms */
            for (int g=0; g < ng; g++){
              _shape_AM_matrix[cell_shape][e * (ng+6) + g + 3] -= wt * 
                (1-material->getDelayedFractionTotal(FORWARD_OUT)) * 
                material->getChiByGroup(e, FORWARD_OUT, temp) *
                material->getNuSigmaFByGroup(g, FORWARD_OUT, temp) / _k_eff_0 * volume;

              _shape_source[cell_shape*ng+e] += (1-wt) * 
                (1-material->getDelayedFractionTotal(PREVIOUS_OUT)) * 
                material->getChiByGroup(e, PREVIOUS_OUT, temp_prev) *
                material->getNuSigmaFByGroup(g, PREVIOUS_OUT, temp_prev) / _k_eff_0
                * flux_prev[celL_shape*ng+g] * volume;
            }
            
            /* Inscattering term on off diagonals */
            for (int g=0; g < ng; g++){
              if (e != g){
                _shape_AM_matrix[cell_shape][e * (ng+6) + g + 3] -= wt * 
                  material->getSigmaSByGroup(g, e, FORWARD_OUT, temp) * volume;
                _shape_source[cell_shape*ng+e] += (1-wt) * 
                  material->getSigmaSByGroup(g, e, PREVIOUS_OUT, temp_prev)
                  * flux_prev[cell_shape*ng+g] * volume;
              }
            }
          }          
           
          /* LEFT SURFACE */
          
          /* Transport term on diagonal */
          _shape_AM_matrix[cell_shape][diag] += wt *
            getDifLinearFineByValue(cell_shape, e, 0, FORWARD_OUT) * height * depth;

          _shape_source[row_shape] -= (1-wt) * flux_prev[row_shape] *
            getDifLinearFineByValue(cell_shape, e, 0, PREVIOUS_OUT) * height * depth;
          
          /* Transport term on off diagonals */
          if (x != 0){
            _shape_AM_matrix[cell_shape][e * (ng+6)] -= wt *
              getDifLinearFineByValue(cell_shape, e, 0, FORWARD_OUT) * height * depth;

            _shape_source[row_shape] += (1-wt) * flux_prev[(cell_shape-1)*ng+e] *
              getDifLinearFineByValue(cell_shape, e, 0, PREVIOUS_OUT) * height * depth;
          }

          /* RIGHT SURFACE */
          
          /* Transport term on diagonal */
          _shape_AM_matrix[cell_shape][diag] += wt *
            getDifLinearFineByValue(cell_shape, e, 3, FORWARD_OUT) * height * depth;

          _shape_source[row_shape] -= (1-wt) * flux_prev[row_shape] *
            getDifLinearFineByValue(cell_shape, e, 3, PREVIOUS_OUT) * height * depth;
          
          /* Transport term on off diagonals */
          if (x != nx - 1){
            _shape_AM_matrix[cell_shape][e * (ng+6) + ng + 3] -= wt * 
              getDifLinearFineByValue(cell_shape, e, 3, FORWARD_OUT) * height * depth;
            _shape_source[row_shape] += (1-wt) * flux_prev[(cell_shape+1)*ng+e] *
              getDifLinearFineByValue(cell_shape, e, 3, PREVIOUS_OUT) * height * depth;
          }
          
          /* BACK SURFACE */
          
          /* Transport term on diagonal */
          _shape_AM_matrix[cell_shape][diag] += wt *
            getDifLinearFineByValue(cell_shape, e, 1, FORWARD_OUT) * width * depth;

          _shape_source[row_shape] -= (1-wt) * flux_prev[row_shape] *
            getDifLinearFineByValue(cell_shape, e, 1, PREVIOUS_OUT) * width * depth;
                    
          /* Transport term on off diagonals */
          if (y != 0){
            _shape_AM_matrix[cell_shape][e * (ng+6) + 1] -= wt *
              getDifLinearFineByValue(cell_shape, e, 1, FORWARD_OUT) * width * depth;

            _shape_source[row_shape] += (1-wt) * flux_prev[(cell_shape-nx)*ng+e] *
              getDifLinearFineByValue(cell_shape, e, 1, PREVIOUS_OUT) * width * depth;
          }
          
          /* FRONT SURFACE */
          
          /* Transport term on diagonal */
          _shape_AM_matrix[cell_shape][diag] += wt *
            getDifLinearFineByValue(cell_shape, e, 4, FORWARD_OUT) * width * depth;

          _shape_source[row_shape] -= (1-wt) * flux_prev[row_shape] *
            getDifLinearFineByValue(cell_shape, e, 4, PREVIOUS_OUT) * width * depth;
          
          
          /* Transport term on off diagonals */
          if (y != ny - 1){
            _shape_AM_matrix[cell_shape][e * (ng+6) + ng + 4] -= wt *
              getDifLinearFineByValue(cell_shape, e, 4, FORWARD_OUT) * width * depth;

            _shape_source[row_shape] += (1-wt) * flux_prev[(cell_shape+nx)*ng+e] *
              getDifLinearFineByValue(cell_shape, e, 4, PREVIOUS_OUT) * width * depth;
          }

          /* BOTTOM SURFACE */
          
          /* Transport term on diagonal */
          _shape_AM_matrix[cell_shape][diag] += wt *
            getDifLinearFineByValue(cell_shape, e, 2, FORWARD_OUT) * width * height;

          _shape_source[row_shape] -= (1-wt) * flux_prev[row_shape] *
            getDifLinearFineByValue(cell_shape, e, 2, PREVIOUS_OUT) * width * height;
          
          
          /* Transport term on off diagonals */
          if (z != 0){
            _shape_AM_matrix[cell_shape][e * (ng+6) + 2] -= wt *
              getDifLinearFineByValue(cell_shape, e, 2, FORWARD_OUT) * width * height;

            _shape_source[row_shape] += (1-wt) * flux_prev[(cell_shape-nx)*ng+e] *
              getDifLinearFineByValue(cell_shape, e, 2, PREVIOUS_OUT) * width * height;
          }
          
          /* TOP SURFACE */
          
          /* Transport term on diagonal */
          _shape_AM_matrix[cell_shape][diag] += wt *
            getDifLinearFineByValue(cell_shape, e, 5, FORWARD_OUT) * width * height;

          _shape_source[row_shape] -= (1-wt) * flux_prev[row_shape] *
            getDifLinearFineByValue(cell_shape, e, 5, PREVIOUS_OUT) * width * height;
          
          
          /* Transport term on off diagonals */
          if (z != nz - 1){
            _shape_AM_matrix[cell_shape][e * (ng+6) + ng + 5] -= wt * 
              getDifLinearFineByValue(cell_shape, e, 5, FORWARD_OUT) * width * height;

            _shape_source[row_shape] += (1-wt) * flux_prev[(cell_shape+nx)*ng+e] *
              getDifLinearFineByValue(cell_shape, e, 5, PREVIOUS_OUT) * width * height;
          }
        }
      }
    }
  }
}


void SolverDiffusion::generateAmpCurrent(int time){

  int nxa = getNumXAmp();
  int nya = getNumYAmp();
  int nza = getNumZAmp();
  int nxs = getNumXShape();
  int nys = getNumYShape();
  int nzs = getNumZShape();

  std::vector< std::vector<int> > amp_to_shape = _geometry->getAmpToShapeMap();
  
  for (int z=0; x < nza; z++){
    for (int y=0; x < nya; y++){
      for (int x=0; x < nxa; x++){
        int cell_amp = z*nxa*nya + y*nxa + x;
        for (int e=0; e < _num_energy_groups; e++){
          for (int s=0; s < 6; s++){

            double current = 0.0;

            std::vector<int>::iterator iter;
            for (iter = amp_to_shape[cell_amp].begin(); 
                 iter != amp_to_shape[cell_amp].end(); ++iter){

              int zs = (*iter)  / (nxs*nys);
              int ys = ((*iter) - zs*nxs*nys) / nxs;
              int xs = ((*iter) - zs*nxs*nys) % nxs;

              int neighbor_cell = _geometry_diffusion->getNeighborShapeCell(xs,ys,zs,s);
              double sense;
              
              if (s == 0 || s == 1 || s == 2)
                sense = -1.0;
              else
                sense = 1.0;
              
              /* Check if shape surface coincides with amp surface */
              if ((s == 0 && xs % (nxs/nxa) == 0) || 
                  (s == 1 && ys % (nys/nya) == 0) ||
                  (s == 2 && zs % (nzs/nza) == 0) ||
                  (s == 3 && xs % (nxs/nxa) == nxs/nxa - 1) ||
                  (s == 4 && ys % (nys/nya) == nys/nya - 1) ||
                  (s == 5 && zs % (nzs/nza) == nzs/nza - 1)){
                if (neighbor_cell == -1)
                  current += sense * getDifLinearFineByValue(cell_shape, e, s, time) * 
                    height * depth * getFluxByValue(cell_shape, e, time);
                else{
                  current += sense * getDifLinearFineByValue(cell_shape, e, s, time) * 
                    height * depth * (getFluxByValue(cell_shape, e, time) -
                                      getFluxByValue(neighbor_cell, e, time));
                }
              }
            }

            setCurrentByValue(current, cell_amp, e, s, CURRENT);
          }
        }
      }
    }
  }
}


void SolverDiffusion::computeInitialShape(double tol){

  /* Generate the initial shape matrices */
  computeDiffusionCoefficientsFine(CURRENT);
  generateInitialShapeMatrices(CURRENT);

  int nx = _num_x_shape;
  int ny = _num_y_shape;
  int nz = _num_z_shape;
  int ng = _num_energy_groups;
  double* old_source = new double[nx*ny*nz*ng];
  double* flux_temp = new double[nx*ny*nz*ng];

  /* Solve the initial eigenvalue problem */
  _k_eff_0 = eigenvalueSolve2d(_shape_A_matrix, nx*ny, ng*(ng+6), _shape_M_matrix, 
                               nx*ny, ng*ng, getFlux(CURRENT), nx*ny*ng, _shape_source, 
                               nx*ny*ng, old_source, nx*ny*ng, flux_temp, nx*ny*ng, ng, 
                               nx, ny, nz, tol);

  log_printf(NORMAL, "done solving eigenvalue problem");

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
  
  delete [] old_source;
  delete [] flux_temp;
}


void SolverDiffusion::takeOuterStepOnly(){

  double tolerance = 1.e-4;
  double residual = 1.e10;

  /* Take step forward with clock */
  _clock->takeOuterStep();

  /* Integrate precursor concentrations */
  integratePrecursorConcentrations(CURRENT, FORWARD_OUT);

  /* Integrate the temperature to the forward time step */
  integrateTemperature(CURRENT, FORWARD_OUT);

  int nx = _num_x_shape;
  int ny = _num_y_shape;
  int nz = _num_z_shape;
  int ng = _num_energy_groups;
  double* old_source = new double[nx*ny*nz*ng];
  double* flux_temp = new double[nx*ny*nz*ng];

  while(true){

    /* Generate the shape matrices */
    computeDiffusionCoefficientsFine(CURRENT);
    computeDiffusionCoefficientsFine(FORWARD_OUT);
    generateShapeMatrices();

    /* Solve the linear problem to get the flux at the forward time step */
    linearSolve2d(_shape_AM_matrix, nx*ny*nz, ng*(ng+6), getFlux(CURRENT), nx*ny*ng, 
                  _shape_source, nx*ny*ng, 
                  flux_temp, nx*ny*ng, nx, ny, nz, ng, tol);
    
    /* Reintegrate precursor concentrations */
    integratePrecursorConcentrations(CURRENT, FORWARD_OUT);
    
    /* Reintegrate the temperature to the forward time step */
    integrateTemperature(CURRENT, FORWARD_OUT);

    /* Check for convergence of the forward power */
    residual = computePowerRMSError(FORWARD_OUT, FORWARD_OUT_OLD);

    if (residual < tolerance)
      break;
    else{
      copyFlux(FORWARD_OUT, FORWARD_OUT_OLD);
    }
  }
  
  /* Broadcast all field variables to all time steps */
  broadcastToAll(CURRENT);
  
  delete [] old_source;
  delete [] flux_temp;
}


void SolverDiffusion::takeInnerStep(){
  return;
}


void SolverDiffusion::takeOuterStep(){
  return;
}
