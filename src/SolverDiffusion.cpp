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
  for (int t=0; t < NUM_TIME_POINTS; t++){
    Vector* dif_linear  = new Vector(num_x, num_y, num_z, _num_energy_groups*6);
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

void SolverDiffusion::generateInitialShapeMatrices(){

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
          
          /* Absorption term on the diagonal */
          val = sig_a * volume;
          _shape_A_matrix->incrementValueByCell(cell, e, cell, e, val);
          
          /* Buckling term on the diagonal */
          val = dif_coef * volume * _buckling;
          _shape_A_matrix->incrementValueByCell(cell, e, cell, e, val);
          
          /* Outscattering term on diagonal */
          for (int g=0; g < ng; g++){
            if (e != g){
              val = material->getSigmaSByGroup(e, g, CURRENT, temp) * volume;
              _shape_A_matrix->incrementValueByCell(cell, e, cell, e, val);
            }
          }

          /* Fission terms */
          for (int g=0; g < ng; g++){
            val = chi * material->getNuSigmaFByGroup(g, CURRENT, temp) * volume;
            _shape_M_matrix->incrementValueByCell(cell, g, cell, e, val);
          }
          
          /* Inscattering term on off diagonals */
          for (int g=0; g < ng; g++){
            if (e != g){
              val = -material->getSigmaSByGroup(g, e, CURRENT, temp) * volume;
              _shape_A_matrix->incrementValueByCell(cell, g, cell, e, val);
            }
          }
                    
          /* X_MIN SURFACE */
          
          /* Transport term on diagonal */
          val = getDifLinearFineByValue(cell, e, SURFACE_X_MIN, CURRENT) * height * depth;
          _shape_A_matrix->incrementValueByCell(cell, e, cell, e, val);
          
          /* Transport term on off diagonals */
          if (x != 0){
            val = - getDifLinearFineByValue(cell, e, SURFACE_X_MIN, CURRENT)
              * height * depth;
            _shape_A_matrix->incrementValueByCell(cell-1, e, cell, e, val);
          }          

          /* X_MAX SURFACE */
          
          /* Transport term on diagonal */
          val = getDifLinearFineByValue(cell, e, SURFACE_X_MAX, CURRENT) * height * depth;
          _shape_A_matrix->incrementValueByCell(cell, e, cell, e, val);
          
          /* Transport term on off diagonals */
          if (x != nx - 1){
            val = - getDifLinearFineByValue(cell, e, SURFACE_X_MAX, CURRENT)
              * height * depth;
            _shape_A_matrix->incrementValueByCell(cell+1, e, cell, e, val);
          }
          
          /* Y_MIN SURFACE */
          
          /* Transport term on diagonal */
          val = getDifLinearFineByValue(cell, e, SURFACE_Y_MIN, CURRENT) * width * depth;
          _shape_A_matrix->incrementValueByCell(cell, e, cell, e, val);

          /* Transport term on off diagonals */
          if (y != 0){
            val = - getDifLinearFineByValue(cell, e, SURFACE_Y_MIN, CURRENT)
              * width * depth;
            _shape_A_matrix->incrementValueByCell(cell-nx, e, cell, e, val);
          }          
          
          /* Y_MAX SURFACE */
          
          /* Transport term on diagonal */
          val = getDifLinearFineByValue(cell, e, SURFACE_Y_MAX, CURRENT) * width * depth;
          _shape_A_matrix->incrementValueByCell(cell, e, cell, e, val);
          
          /* Transport term on off diagonals */
          if (y != ny - 1){
            val = - getDifLinearFineByValue(cell, e, SURFACE_Y_MAX, CURRENT)
              * width * depth;
            _shape_A_matrix->incrementValueByCell(cell+nx, e, cell, e, val);
          }          

          /* Z_MIN SURFACE */
          
          /* Transport term on diagonal */
          val = getDifLinearFineByValue(cell, e, SURFACE_Z_MIN, CURRENT) * width * height;
          _shape_A_matrix->incrementValueByCell(cell, e, cell, e, val);
          
          /* Transport term on off diagonals */
          if (z != 0){
            val = - getDifLinearFineByValue(cell, e, SURFACE_Z_MIN, CURRENT)
              * width * height;
            _shape_A_matrix->incrementValueByCell(cell-nx*ny, e, cell, e, val);
          } 
          
          /* Z_MAX SURFACE */
          
          /* Transport term on diagonal */
          val = getDifLinearFineByValue(cell, e, SURFACE_Z_MAX, CURRENT) * width * height;
          _shape_A_matrix->incrementValueByCell(cell, e, cell, e, val);
          
          /* Transport term on off diagonals */
          if (z != nz - 1){
            val = - getDifLinearFineByValue(cell, e, SURFACE_Z_MAX, CURRENT)
              * width * height;
            _shape_A_matrix->incrementValueByCell(cell+nx*ny, e, cell, e, val);
          }
        }
      }
    }    
  }
}


void SolverDiffusion::computeDiffusionCoefficientsFine(int time){

  int nx = _geometry_diffusion->getNumXShape();
  int ny = _geometry_diffusion->getNumYShape();
  int nz = _geometry_diffusion->getNumZShape();
  double width  = _geometry->getWidth()  / nx;
  double height = _geometry->getHeight() / ny;
  double depth  = _geometry->getDepth()  / nz;
  std::vector< std::vector<int> > amp_to_shape = _geometry->getAmpToShapeMap();

  /* Compute the linear and nonlinear diffusion coefficients */
  for (int z=0; z < nz; z++){
    #pragma omp parallel for
    for (int y=0; y < ny; y++){
      for (int x=0; x < nx; x++){
        
        int cell = z*nx*ny + y*nx + x;

        double length, length_perpen;
        double dif_linear, dif_coef_next, dif_coef;
        Material* material = _geometry->getMaterial(cell);
        Material* material_next;
        
        for (int s=0; s < 6; s++){
          
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
            
            dif_coef = material->getDifCoefByGroup(g, time);
            
            if (cell_next == -1){
              dif_linear = 2 * dif_coef / length_perpen / 
                (1 + 4 * dif_coef / length_perpen);
              dif_linear *= (1 - _geometry->getBoundary(s));
            }
            else{
              material_next = _geometry->getMaterial(cell_next);
              dif_coef_next = material_next->getDifCoefByGroup(g, time);
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


Vector* SolverDiffusion::getDifLinearFine(int time){
  return _dif_linear_fine[time];
}


double SolverDiffusion::getDifLinearFineByValue(int cell, int group, int side, int time){
  return _dif_linear_fine[time]->getValueByCell(cell, side*_num_energy_groups+group);
}


void SolverDiffusion::setDifLinearFineByValue(double value, int cell, int group, 
                                              int side, int time){
  _dif_linear_fine[time]->setValueByCell(cell, side*_num_energy_groups+group, value);
}


void SolverDiffusion::generateShapeMatrices(double wt){

  int nx = _geometry_diffusion->getNumXShape();
  int ny = _geometry_diffusion->getNumYShape();
  int nz = _geometry_diffusion->getNumZShape();
  int ng = _num_energy_groups;
  double width  = _geometry->getWidth()  / nx;
  double height = _geometry->getHeight() / ny;
  double depth  = _geometry->getDepth()  / nz;
  double val;
  double dt = _clock->getDtOuter();
  int* shape_to_amp = _geometry->getShapeToAmpMap();
  double volume = width * height * depth;
  
  _shape_source->zero();
  _shape_AM_matrix->clear();

  for (int z=0; z < nz; z++){
    #pragma omp parallel for private(val)
    for (int y=0; y < ny; y++){
      for (int x=0; x < nx; x++){

        int cell_shape     = z*nx*ny + y*nx + x;
        int cell_amp       = shape_to_amp[cell_shape];
        Material* material = _geometry->getMaterial(cell_shape);
        double temp        = getTemperatureByValue(cell_shape, FORWARD_OUT);
        double temp_prev   = getTemperatureByValue(cell_shape, PREVIOUS_OUT);
        double volume      = _geometry->getVolume(cell_shape);
        
        for (int e=0; e < ng; e++){

          double v           = material->getVelocityByGroup(e, FORWARD_OUT, temp);
          double v_prev      = material->getVelocityByGroup(e, PREVIOUS_OUT, temp);
          double flux        = getFluxByValue(cell_shape, e, FORWARD_OUT);
          double flux_prev   = getFluxByValue(cell_shape, e, PREVIOUS_OUT);
          double amp         = getAmplitudeByValue(cell_amp, e, FORWARD_OUT);
          double amp_prev    = getAmplitudeByValue(cell_amp, e, PREVIOUS_OUT);
          double freq        = getFrequencyByValue(cell_amp, e, FORWARD_OUT);          
          double chi         = material->getChiByGroup(e, FORWARD_OUT, temp);
          double chi_prev    = material->getChiByGroup(e, PREVIOUS_OUT, temp);          
          double sig_a       = material->getSigmaAByGroup(e, FORWARD_OUT, temp); 
          double sig_a_prev  = material->getSigmaAByGroup(e, PREVIOUS_OUT, temp);
          double dif_coef    = material->getDifCoefByGroup(e, FORWARD_OUT, temp);
          double dif_coef_prev = material->getDifCoefByGroup(e, PREVIOUS_OUT, temp);
          double delay       = material->getDelayedFractionTotal(FORWARD_OUT);
          double delay_prev  = material->getDelayedFractionTotal(PREVIOUS_OUT);
          double sig_s, sig_s_prev;
          double flux_prev_g;
          double nu_sig_f, nu_sig_f_prev;
          
          /* Time absorption term on the diagonal */
          if (_method == IQS){
            val = 1.0 / dt / v * volume * flux;
            _shape_AM_matrix->incrementValueByCell(cell_shape, e, cell_shape, e, val);
            
            val = 1.0 / v / amp * volume * freq;
            _shape_AM_matrix->incrementValueByCell(cell_shape, e, cell_shape, e, val);

            val = flux_prev / dt / v_prev * volume * amp / amp_prev;
            _shape_source->incrementValueByCell(cell_shape, e, val);
          }
          else if (_method == THETA){
            val = 1.0 / dt / v * volume;
            _shape_AM_matrix->incrementValueByCell(cell_shape, e, cell_shape, e, val);
            
            val = flux_prev / dt / v_prev * volume;
            _shape_source->incrementValueByCell(cell_shape, e, val);
          }
          
          /* Delayed neutron precursors */
          for (int d=0; d < _num_delayed_groups; d++){
            double decay     = material->getDecayConstantByGroup(d);
            double conc      = material->getPrecursorConcByGroup(d, FORWARD_OUT);
            double conc_prev = material->getPrecursorConcByGroup(d, PREVIOUS_OUT);
            
            val = wt * chi * decay * conc * volume;
            _shape_source->incrementValueByCell(cell_shape, e, val);
            
            val = (1-wt) * chi_prev * decay * conc_prev * volume;
            _shape_source->incrementValueByCell(cell_shape, e, val);
          }
          
          /* Absorption term on the diagonal */
          val = wt * sig_a * volume;
          _shape_AM_matrix->incrementValueByCell(cell_shape, e, cell_shape, e, val);
          
          val = - (1-wt) * sig_a_prev * flux_prev * volume;
          _shape_source->incrementValueByCell(cell_shape, e, val);
          
          /* Buckling term on the diagonal */
          val = wt * dif_coef * volume * _buckling;
          _shape_AM_matrix->incrementValueByCell(cell_shape, e, cell_shape, e, val);
          
          val = - (1-wt) * dif_coef_prev * volume * _buckling;
          _shape_source->incrementValueByCell(cell_shape, e, val);
          
          /* Outscattering term on diagonal */
          for (int g=0; g < ng; g++){
            if (e != g){
              sig_s      = material->getSigmaSByGroup(e, g, FORWARD_OUT, temp);
              sig_s_prev = material->getSigmaSByGroup(e, g, PREVIOUS_OUT, temp);
              
              val = wt * sig_s * volume;
              _shape_AM_matrix->incrementValueByCell(cell_shape, e, cell_shape, e, val);
              
              val = - (1-wt) * sig_s_prev * volume;
              _shape_source->incrementValueByCell(cell_shape, e, val);
            }
          }
          
          /* Fission terms */
          for (int g=0; g < ng; g++){
            nu_sig_f      = material->getNuSigmaFByGroup(g, FORWARD_OUT, temp_prev);
            nu_sig_f_prev = material->getNuSigmaFByGroup(g, PREVIOUS_OUT, temp_prev);
            flux_prev_g = getFluxByValue(cell_shape, g, PREVIOUS_OUT);
            
            val = - wt * (1-delay) * chi * nu_sig_f / _k_eff_0 * volume;
            _shape_AM_matrix->incrementValueByCell(cell_shape, g, cell_shape, e, val);
            
            val = (1-wt) * (1-delay_prev) * chi_prev * nu_sig_f_prev / _k_eff_0
              * flux_prev_g * volume;
            _shape_source->incrementValueByCell(cell_shape, e, val);
          }
          
          /* Inscattering term on off diagonals */
          for (int g=0; g < ng; g++){
            if (e != g){
              sig_s      = material->getSigmaSByGroup(g, e, FORWARD_OUT, temp);
              sig_s_prev = material->getSigmaSByGroup(g, e, PREVIOUS_OUT, temp);
              flux_prev_g = getFluxByValue(cell_shape, g, PREVIOUS_OUT);
              
              val = - wt * sig_s * volume;
              _shape_AM_matrix->incrementValueByCell(cell_shape, g, cell_shape, e, val);
              
              val = (1-wt) * sig_s_prev * flux_prev_g * volume;
                _shape_source->incrementValueByCell(cell_shape, e, val);
            }
          }
        
          /* X_MIN SURFACE */
          
          /* Transport term on diagonal */
          val = wt * getDifLinearFineByValue(cell_shape, e, SURFACE_X_MIN, FORWARD_OUT)
            * height * depth;
          _shape_AM_matrix->incrementValueByCell(cell_shape, e, cell_shape, e, val);

          val = - (1-wt) * flux_prev * height * depth * 
            getDifLinearFineByValue(cell_shape, e, SURFACE_X_MIN, PREVIOUS_OUT);
          _shape_source->incrementValueByCell(cell_shape, e, val);
          
          /* Transport term on off diagonals */
          if (x != 0){
            val = -wt * getDifLinearFineByValue(cell_shape, e, SURFACE_X_MIN, FORWARD_OUT)
              * height * depth;
            _shape_AM_matrix->incrementValueByCell(cell_shape-1, e, cell_shape, e, val);
            
            val = (1-wt) * getFluxByValue(cell_shape - 1, e, PREVIOUS_OUT) * height * depth
              * getDifLinearFineByValue(cell_shape, e, SURFACE_X_MIN, PREVIOUS_OUT);
            _shape_source->incrementValueByCell(cell_shape, e, val);
          }

          /* X_MAX SURFACE */
          
          /* Transport term on diagonal */
          val = wt * getDifLinearFineByValue(cell_shape, e, SURFACE_X_MAX, FORWARD_OUT)
            * height * depth;
          _shape_AM_matrix->incrementValueByCell(cell_shape, e, cell_shape, e, val);

          val = - (1-wt) * flux_prev * height * depth * 
            getDifLinearFineByValue(cell_shape, e, SURFACE_X_MAX, PREVIOUS_OUT);
          _shape_source->incrementValueByCell(cell_shape, e, val);
          
          /* Transport term on off diagonals */
          if (x != nx - 1){
            val = -wt * getDifLinearFineByValue(cell_shape, e, SURFACE_X_MAX, FORWARD_OUT)
              * height * depth;
            _shape_AM_matrix->incrementValueByCell(cell_shape+1, e, cell_shape, e, val);
            
            val = (1-wt) * getFluxByValue(cell_shape+1, e, PREVIOUS_OUT) * height * depth
              * getDifLinearFineByValue(cell_shape, e, SURFACE_X_MAX, PREVIOUS_OUT);
            _shape_source->incrementValueByCell(cell_shape, e, val);
          }
          
          /* Y_MIN SURFACE */
          
          /* Transport term on diagonal */
          val = wt * getDifLinearFineByValue(cell_shape, e, SURFACE_Y_MIN, FORWARD_OUT)
            * width * depth;
          _shape_AM_matrix->incrementValueByCell(cell_shape, e, cell_shape, e, val);

          val = - (1-wt) * flux_prev * width * depth * 
            getDifLinearFineByValue(cell_shape, e, SURFACE_Y_MIN, PREVIOUS_OUT);
          _shape_source->incrementValueByCell(cell_shape, e, val);
                    
          /* Transport term on off diagonals */
          if (y != 0){
            val = -wt * getDifLinearFineByValue(cell_shape, e, SURFACE_Y_MIN, FORWARD_OUT)
              * width * depth;
            _shape_AM_matrix->incrementValueByCell(cell_shape-nx, e, cell_shape, e, val);
            
            val = (1-wt) * getFluxByValue(cell_shape-nx, e, PREVIOUS_OUT) * width * depth
              * getDifLinearFineByValue(cell_shape, e, SURFACE_Y_MIN, PREVIOUS_OUT);
            _shape_source->incrementValueByCell(cell_shape, e, val);
          }
          
          /* Y_MAX SURFACE */
          
          /* Transport term on diagonal */
          val = wt * getDifLinearFineByValue(cell_shape, e, SURFACE_Y_MAX, FORWARD_OUT)
            * width * depth;
          _shape_AM_matrix->incrementValueByCell(cell_shape, e, cell_shape, e, val);

          val = - (1-wt) * flux_prev * width * depth * 
            getDifLinearFineByValue(cell_shape, e, SURFACE_Y_MAX, PREVIOUS_OUT);
          _shape_source->incrementValueByCell(cell_shape, e, val);
          
          /* Transport term on off diagonals */
          if (y != ny - 1){
            val = -wt * getDifLinearFineByValue(cell_shape, e, SURFACE_Y_MAX, FORWARD_OUT)
              * width * depth;
            _shape_AM_matrix->incrementValueByCell(cell_shape+nx, e, cell_shape, e, val);
            
            val = (1-wt) * getFluxByValue(cell_shape+nx, e, PREVIOUS_OUT) * width * depth
              * getDifLinearFineByValue(cell_shape, e, SURFACE_Y_MAX, PREVIOUS_OUT);
            _shape_source->incrementValueByCell(cell_shape, e, val);
          }

          /* Z_MIN SURFACE */
          
          /* Transport term on diagonal */
          val = wt * getDifLinearFineByValue(cell_shape, e, SURFACE_Z_MIN, FORWARD_OUT)
            * width * height;
          _shape_AM_matrix->incrementValueByCell(cell_shape, e, cell_shape, e, val);

          val = - (1-wt) * flux_prev * width * height * 
            getDifLinearFineByValue(cell_shape, e, SURFACE_Z_MIN, PREVIOUS_OUT);
          _shape_source->incrementValueByCell(cell_shape, e, val);
          
          
          /* Transport term on off diagonals */
          if (z != 0){
            val = -wt * getDifLinearFineByValue(cell_shape, e, SURFACE_Z_MIN, FORWARD_OUT)
              * width * height;
            _shape_AM_matrix->incrementValueByCell(cell_shape-nx*ny, e, cell_shape, e, val);
            
            val = (1-wt) * getFluxByValue(cell_shape-nx*ny, e, PREVIOUS_OUT) * width * height
              * getDifLinearFineByValue(cell_shape, e, SURFACE_Z_MIN, PREVIOUS_OUT);
            _shape_source->incrementValueByCell(cell_shape, e, val);
          }
          
          /* Z_MAX SURFACE */
          
          /* Transport term on diagonal */
          val = wt * getDifLinearFineByValue(cell_shape, e, SURFACE_Z_MAX, FORWARD_OUT)
            * width * height;
          _shape_AM_matrix->incrementValueByCell(cell_shape, e, cell_shape, e, val);

          val = - (1-wt) * flux_prev * width * height * 
            getDifLinearFineByValue(cell_shape, e, SURFACE_Z_MAX, PREVIOUS_OUT);
          _shape_source->incrementValueByCell(cell_shape, e, val);
          
          /* Transport term on off diagonals */
          if (z != nz - 1){
            val = -wt * getDifLinearFineByValue(cell_shape, e, SURFACE_Z_MAX, FORWARD_OUT)
              * width * height;
            _shape_AM_matrix->incrementValueByCell(cell_shape+nx*ny, e, cell_shape, e, val);
            
            val = (1-wt) * getFluxByValue(cell_shape+nx*ny, e, PREVIOUS_OUT) * width * height
              * getDifLinearFineByValue(cell_shape, e, SURFACE_Z_MAX, PREVIOUS_OUT);
            _shape_source->incrementValueByCell(cell_shape, e, val);
          }
        }
      }
    }
  }
}


void SolverDiffusion::generateAmpCurrent(int time){

  int nxa = _geometry->getNumXAmp();
  int nya = _geometry->getNumYAmp();
  int nza = _geometry->getNumZAmp();
  int nxs = _geometry_diffusion->getNumXShape();
  int nys = _geometry_diffusion->getNumYShape();
  int nzs = _geometry_diffusion->getNumZShape();
  double width  = _geometry->getWidth()  / nxa;
  double height = _geometry->getHeight() / nya;
  double depth  = _geometry->getDepth()  / nza;

  std::vector< std::vector<int> > amp_to_shape = _geometry->getAmpToShapeMap();
  
  for (int z=0; z < nza; z++){
    for (int y=0; y < nya; y++){
      for (int x=0; x < nxa; x++){
        int cell_amp = z*nxa*nya + y*nxa + x;
        for (int e=0; e < _num_energy_groups; e++){
          for (int s=0; s < NUM_SURFACES; s++){

            double current = 0.0;

            std::vector<int>::iterator iter;
            for (iter = amp_to_shape[cell_amp].begin(); 
                 iter != amp_to_shape[cell_amp].end(); ++iter){

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
              
              /* Check if shape surface coincides with amp surface */
              if ((s == SURFACE_X_MIN && xs % (nxs/nxa) == 0) || 
                  (s == SURFACE_Y_MIN && ys % (nys/nya) == 0) ||
                  (s == SURFACE_Z_MIN && zs % (nzs/nza) == 0) ||
                  (s == SURFACE_X_MAX && xs % (nxs/nxa) == nxs/nxa - 1) ||
                  (s == SURFACE_Y_MAX && ys % (nys/nya) == nys/nya - 1) ||
                  (s == SURFACE_Z_MAX && zs % (nzs/nza) == nzs/nza - 1)){
                if (neighbor_cell == -1)
                  current += sense * getDifLinearFineByValue(*iter, e, s, time) * 
                    height * depth * getFluxByValue(*iter, e, time);
                else{
                  current += sense * getDifLinearFineByValue(*iter, e, s, time) * 
                    height * depth * (getFluxByValue(*iter, e, time) -
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

  initializeClock();

  /* Generate the initial shape matrices */
  computeDiffusionCoefficientsFine(CURRENT);

  generateInitialShapeMatrices();

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

  double tolerance = 1.e-4;
  double residual = 1.e10;

  /* Take step forward with clock */
  _clock->takeOuterStep();

  /* Integrate precursor concentrations */
  integratePrecursorConcentrations(CURRENT, FORWARD_OUT);

  /* Integrate the temperature to the forward time step */
  integrateTemperature(CURRENT, FORWARD_OUT);

  while(true){

    /* Generate the shape matrices */
    computeDiffusionCoefficientsFine(CURRENT);
    computeDiffusionCoefficientsFine(FORWARD_OUT);
    generateShapeMatrices(0.0);

    /* Solve the linear problem to get the flux at the forward time step */
    linearSolve(_shape_AM_matrix, _flux[CURRENT], _shape_source, tolerance);
    
    /* Reintegrate precursor concentrations */
    integratePrecursorConcentrations(CURRENT, FORWARD_OUT);
    
    /* Reintegrate the temperature to the forward time step */
    integrateTemperature(CURRENT, FORWARD_OUT);

    /* Check for convergence of the forward power */
    residual = computePowerRMSError(FORWARD_OUT, FORWARD_OUT_OLD);

    if (residual < tolerance)
      break;
    else{
      _flux[FORWARD_OUT]->copyTo(_flux[FORWARD_OUT_OLD]);
    }
  }
  
  /* Broadcast all field variables to all time steps */
  broadcastToAll(CURRENT);
}


void SolverDiffusion::takeInnerStep(){
  return;
}


void SolverDiffusion::takeOuterStep(){
  return;
}
