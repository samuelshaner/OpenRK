#include "Solver.h"

Solver::Solver(Geometry* geometry){

  /* Initialize variables */
  _k_eff_0 = 0.0;
  _method = IQS;
  _geometry = geometry;
  _num_energy_groups = _geometry->getNumEnergyGroups();
  _num_delayed_groups = _geometry->getNumDelayedGroups();
  _num_shape_cells = _geometry->getNumShapeCells();
  _num_amp_cells = _geometry->getNumAmpCells();
  _clock = new Clock();
  _buckling = 0.0;
  _initial_power = 1.e-6;
  
  /* Set solver parameters */
  setNumThreads(1);

  /* Initialize coarse mesh matrices */
  int num_x = _geometry->getNumXAmp();
  int num_y = _geometry->getNumYAmp();
  int num_z = _geometry->getNumZAmp();
  _amp_matrix = new Matrix(num_x, num_y, num_z, _num_energy_groups);
  _amp_source = new Vector(num_x, num_y, num_z, _num_energy_groups);

  /* Allocate memory for coarse mesh field variables */
  for (int t=0; t < NUM_STATES; t++){

    Vector* amplitude     = new Vector(num_x, num_y, num_z, _num_energy_groups);
    Vector* current       = new Vector(num_x, num_y, num_z, _num_energy_groups*NUM_SURFACES);
    Vector* dif_linear    = new Vector(num_x, num_y, num_z, _num_energy_groups*NUM_SURFACES);
    Vector* dif_nonlinear = new Vector(num_x, num_y, num_z, _num_energy_groups*NUM_SURFACES);
    Vector* frequency     = new Vector(num_x, num_y, num_z, _num_energy_groups);

    amplitude->setAll(1.0);
    current->setAll(0.0);
    dif_linear->setAll(0.0);
    dif_nonlinear->setAll(0.0);
    frequency->setAll(0.0);

    _amplitude[t] = amplitude;
    _current[t] = current;
    _dif_linear[t] = dif_linear;
    _dif_nonlinear[t] = dif_nonlinear;
    _frequency[t] = frequency;
  }
}


Solver::~Solver(){

  if (_amp_matrix != NULL)
    delete _amp_matrix;

  if (_amp_source != NULL)
    delete _amp_source;
}


Matrix* Solver::getAmpMatrix(){
  return _amp_matrix;
}


Vector* Solver::getAmpSource(){
  return _amp_source;
}


void Solver::generateAmplitudeMatrix(){

  int nx = _geometry->getNumXAmp();
  int ny = _geometry->getNumYAmp();
  int nz = _geometry->getNumZAmp();
  int ng = _num_energy_groups;
  double width  = _geometry->getWidth()  / nx;
  double height = _geometry->getHeight() / ny;
  double depth  = _geometry->getDepth()  / nz;
  double dt = _clock->getDtInner();
  double val;

  std::vector< std::vector<int> >* amp_to_shape = _geometry->getAmpToShapeMap();

  _amp_source->clear();
  _amp_matrix->clear();

  #pragma omp parallel for private(val) collapse(3)
  for (int z=0; z < nz; z++){
    for (int y=0; y < ny; y++){
      for (int x=0; x < nx; x++){
        
        int cell_amp = z*nx*ny + y*nx + x;
        
        for (int e=0; e < ng; e++){

          double amp_prev = getAmplitudeByValue(cell_amp, e, PREVIOUS_IN);
          
          /* Loop over shape cells in this amp cell */
          std::vector<int>::iterator iter;
          for (iter = amp_to_shape->at(cell_amp).begin();
               iter != amp_to_shape->at(cell_amp).end(); ++iter){
            
            int cell_shape       = *iter;
            Material* material   = _geometry->getMaterial(cell_shape);
            double temp          = getTemperatureByValue(cell_shape, CURRENT);
            double volume        = _geometry->getVolume(cell_shape);
            double v             = material->getVelocityByGroup(e, CURRENT, temp);
            double chi           = material->getChiByGroup(e, CURRENT, temp);
            double sig_a         = material->getSigmaAByGroup(e, CURRENT, temp);
            double shape         = getShapeByValue(cell_shape, e, CURRENT);
            double shape_prev    = getShapeByValue(cell_shape, e, PREVIOUS_IN);
            double amp           = getAmplitudeByValue(cell_amp, e, CURRENT);
            double beta          = material->getDelayedFractionTotal(CURRENT);
            double dif_coef      = material->getDifCoefByGroup(e, CURRENT, temp);
            double sig_s;
            double nu_sig_f;
            
            /* Time absorption term on the diagonal */
            val = shape / v * volume / dt;
            _amp_matrix->incrementValue(cell_amp, e, cell_amp, e, val);
            
            val = shape * amp / v * volume / dt;
            _amp_source->incrementValue(cell_amp, e, val);

            if (_method == IQS){
              val = (shape - shape_prev) / v * volume / dt;
              _amp_matrix->incrementValue(cell_amp, e, cell_amp, e, val);
            }
                    
            /* Delayed neutron precursors */
            for (int d=0; d < _num_delayed_groups; d++){
              double decay = material->getDecayConstantByGroup(d);              
              
              val = chi * decay * volume *
                material->getPrecursorConcByGroup(d, CURRENT);
              _amp_source->incrementValue(cell_amp, e, val);
            }
          
            /* Absorption term on the diagonal */
            val = sig_a * volume * shape;
            _amp_matrix->incrementValue(cell_amp, e, cell_amp, e, val);
 
            /* Buckling term on the diagonal */
            val = dif_coef * volume * _buckling * shape;
            _amp_matrix->incrementValue(cell_amp, e, cell_amp, e, val);
            
            /* Outscattering term on diagonal */
            for (int g=0; g < ng; g++){
              if (e != g){
                sig_s = material->getSigmaSByGroup(e,g, CURRENT, temp);
                
                val = sig_s * volume * shape;
                _amp_matrix->incrementValue(cell_amp, e, cell_amp, e, val);
              }
            }
            
            /* Fission terms */
            for (int g=0; g < ng; g++){
              nu_sig_f = material->getNuSigmaFByGroup(g, CURRENT, temp);
                            
              val = - (1-beta) * chi * nu_sig_f / _k_eff_0 * volume *
                getShapeByValue(cell_shape, g, CURRENT);
              _amp_matrix->incrementValue(cell_amp, g, cell_amp, e, val);
            }
            
            /* Inscattering term on off diagonals */
            for (int g=0; g < ng; g++){
              if (e != g){
                sig_s = material->getSigmaSByGroup(g,e, CURRENT, temp);
                
                val = - sig_s * volume * getShapeByValue(cell_shape, g, CURRENT);
                _amp_matrix->incrementValue(cell_amp, g, cell_amp, e, val);
              }
            }
          }          
           
          /* X_MIN SURFACE */
          
          /* Transport term on diagonal */
          val = (getDifLinearByValue(cell_amp, e, SURFACE_X_MIN, CURRENT) +
                 getDifNonlinearByValue(cell_amp, e, SURFACE_X_MIN, CURRENT))
            * height * depth;

          _amp_matrix->incrementValue(cell_amp, e, cell_amp, e, val);
          
          /* Transport term on off diagonals */
          if (x != 0){
            val = - (getDifLinearByValue(cell_amp, e, SURFACE_X_MIN, CURRENT) -
                     getDifNonlinearByValue(cell_amp, e, SURFACE_X_MIN, CURRENT))
              * height * depth;
            _amp_matrix->incrementValue(cell_amp - 1, e, cell_amp, e, val);
          }          

          /* X_MAX SURFACE */
          
          /* Transport term on diagonal */
          val = (getDifLinearByValue(cell_amp, e, SURFACE_X_MAX, CURRENT) -
                 getDifNonlinearByValue(cell_amp, e, SURFACE_X_MAX, CURRENT))
            * height * depth;
          _amp_matrix->incrementValue(cell_amp, e, cell_amp, e, val);
          
          /* Transport term on off diagonals */
          if (x != nx - 1){
            val = - (getDifLinearByValue(cell_amp, e, SURFACE_X_MAX, CURRENT) +
                     getDifNonlinearByValue(cell_amp, e, SURFACE_X_MAX, CURRENT))
              * height * depth;
            _amp_matrix->incrementValue(cell_amp + 1, e, cell_amp, e, val);
          }
          
          /* Y_MIN SURFACE */
          
          /* Transport term on diagonal */
          val = (getDifLinearByValue(cell_amp, e, SURFACE_Y_MIN, CURRENT) +
                 getDifNonlinearByValue(cell_amp, e, SURFACE_Y_MIN, CURRENT))
            * width * depth;
          _amp_matrix->incrementValue(cell_amp, e, cell_amp, e, val);
          
          /* Transport term on off diagonals */
          if (y != 0){
            val = - (getDifLinearByValue(cell_amp, e, SURFACE_Y_MIN, CURRENT) -
                     getDifNonlinearByValue(cell_amp, e, SURFACE_Y_MIN, CURRENT))
              * width * depth;
            _amp_matrix->incrementValue(cell_amp - nx, e, cell_amp, e, val);
          }
          
          /* Y_MAX SURFACE */
          
          /* Transport term on diagonal */
          val = (getDifLinearByValue(cell_amp, e, SURFACE_Y_MAX, CURRENT) -
                 getDifNonlinearByValue(cell_amp, e, SURFACE_Y_MAX, CURRENT))
            * width * depth;
          _amp_matrix->incrementValue(cell_amp, e, cell_amp, e, val);

          /* Transport term on off diagonals */
          if (y != ny - 1){
            val = - (getDifLinearByValue(cell_amp, e, SURFACE_Y_MAX, CURRENT) +
                     getDifNonlinearByValue(cell_amp, e, SURFACE_Y_MAX, CURRENT))
              * width * depth;
            _amp_matrix->incrementValue(cell_amp + nx, e, cell_amp, e, val);
          }

          /* Z_MIN SURFACE */
          
          /* Transport term on diagonal */
          val = (getDifLinearByValue(cell_amp, e, SURFACE_Z_MIN, CURRENT) +
                 getDifNonlinearByValue(cell_amp, e, SURFACE_Z_MIN, CURRENT))
            * width * height;
          _amp_matrix->incrementValue(cell_amp, e, cell_amp, e, val);

          /* Transport term on off diagonals */
          if (z != 0){
            val = - (getDifLinearByValue(cell_amp, e, SURFACE_Z_MIN, CURRENT) -
                     getDifNonlinearByValue(cell_amp, e, SURFACE_Z_MIN, CURRENT))
              * width * height;
            _amp_matrix->incrementValue(cell_amp - nx*ny, e, cell_amp, e, val);
          }
          
          /* Z_MAX SURFACE */
          
          /* Transport term on diagonal */
          val = (getDifLinearByValue(cell_amp, e, SURFACE_Z_MAX, CURRENT) -
                 getDifNonlinearByValue(cell_amp, e, SURFACE_Z_MAX, CURRENT))
            * width * height;
          _amp_matrix->incrementValue(cell_amp, e, cell_amp, e, val);
          
          /* Transport term on off diagonals */
          if (z != nz - 1){
            val = - (getDifLinearByValue(cell_amp, e, SURFACE_Z_MAX, CURRENT) +
                     getDifNonlinearByValue(cell_amp, e, SURFACE_Z_MAX, CURRENT))
              * width * height;
            _amp_matrix->incrementValue(cell_amp + nx*ny, e, cell_amp, e, val);
          }
        }
      }
    }
  }
}


Vector* Solver::getTemperature(int state){
  return _temperature[state];
}


Vector* Solver::getFlux(int state){
  return _flux[state];
}


Vector* Solver::getShape(int state){
  return _shape[state];
}


Vector* Solver::getPower(int state){
  return _power[state];
}


Vector* Solver::getAmplitude(int state){
  return _amplitude[state];
}


Vector* Solver::getCurrent(int state){
  return _current[state];
}


Vector* Solver::getDifLinear(int state){
  return _dif_linear[state];
}


Vector* Solver::getDifNonlinear(int state){
  return _dif_nonlinear[state];
}


Vector* Solver::getFrequency(int state){
  return _frequency[state];
}


double Solver::getTemperatureByValue(int cell, int state){
  return _temperature[state]->getValue(cell, 0);
}


double Solver::getFluxByValue(int cell, int group, int state){
  return _flux[state]->getValue(cell, group);
}


double Solver::getShapeByValue(int cell, int group, int state){
  return _shape[state]->getValue(cell, group);
}


double Solver::getPowerByValue(int cell, int state){
  return _power[state]->getValue(cell, 0);
}


double Solver::getAmplitudeByValue(int cell, int group, int state){
  return _amplitude[state]->getValue(cell, group);
}


double Solver::getCurrentByValue(int cell, int group, int side, int state){
  return _current[state]->getValue(cell, side*_num_energy_groups+group);
}


double Solver::getDifLinearByValue(int cell, int group, int side, int state){
  return _dif_linear[state]->getValue(cell, side*_num_energy_groups+group);
}


double Solver::getDifNonlinearByValue(int cell, int group, int side, int state){
  return _dif_nonlinear[state]->getValue(cell, side*_num_energy_groups+group);
}


double Solver::getFrequencyByValue(int cell, int group, int state){
  return _frequency[state]->getValue(cell, group);
}


void Solver::copyPrecursors(int state_from, int state_to){

  for (int i=0; i < _num_shape_cells; i++){
    Material* material = _geometry->getMaterial(i);
    if (material->isFissionable()){
      for (int d=0; d < _num_delayed_groups; d++){
        double conc = material->getPrecursorConcByGroup(d, state_from);
        material->setPrecursorConcByGroup(conc, d, state_to);
      }
    }
  }
}


void Solver::copyFieldVariables(int state_from, int state_to){

  copyPrecursors(state_from, state_to);

  _temperature  [state_from]->copyTo(_temperature  [state_to]);
  _flux         [state_from]->copyTo(_flux         [state_to]);
  _shape        [state_from]->copyTo(_shape        [state_to]);
  _power        [state_from]->copyTo(_power        [state_to]);
  _amplitude    [state_from]->copyTo(_amplitude    [state_to]);
  _current      [state_from]->copyTo(_current      [state_to]);
  _dif_linear   [state_from]->copyTo(_dif_linear   [state_to]);
  _dif_nonlinear[state_from]->copyTo(_dif_nonlinear[state_to]);
  _frequency    [state_from]->copyTo(_frequency    [state_to]);
}


void Solver::broadcastToAll(int state_from){

  for (int t=0; t < NUM_STATES; t++)
    copyFieldVariables(state_from, t);

}


void Solver::setTemperatureByValue(double value, int cell, int state){
  _temperature[state]->setValue(cell, 0,  value);
}


void Solver::setFluxByValue(double value, int cell, int group, int state){
  _flux[state]->setValue(cell, group, value);
}


void Solver::setShapeByValue(double value, int cell, int group, int state){
  _shape[state]->setValue(cell, group, value);
}


void Solver::setPowerByValue(double value, int cell, int state){
  _power[state]->setValue(cell, 0, value);
}


void Solver::setAmplitudeByValue(double value, int cell, int group, int state){
  _amplitude[state]->setValue(cell, group, value);
}


void Solver::setCurrentByValue(double value, int cell, int group, int side, int state){
  _current[state]->setValue(cell, side*_num_energy_groups+group, value);
}


void Solver::setDifLinearByValue(double value, int cell, int group, int side, int state){
  _dif_linear[state]->setValue(cell, side*_num_energy_groups+group, value);
}


void Solver::setDifNonlinearByValue(double value, int cell, int group, int side, int state){
  _dif_nonlinear[state]->setValue(cell, side*_num_energy_groups+group, value);
}


void Solver::setFrequencyByValue(double value, int cell, int group, int state){
  _frequency[state]->setValue(cell, group, value);
}


void Solver::integratePrecursorConcentrations(int state_from, int state_to){

  double dt = _clock->getTime(state_to) - _clock->getTime(state_from);

  #pragma omp parallel for
  for (int i=0; i < _num_shape_cells; i++){

    double fission_rate_from = 0.0;
    double fission_rate_to   = 0.0;
    double temp_from         = getTemperatureByValue(i, state_from);
    double temp_to           = getTemperatureByValue(i, state_to);
    Material* material = _geometry->getMaterial(i);
    
    if (material->isFissionable()){
      for (int g=0; g < _num_energy_groups; g++){
        fission_rate_from += material->getNuSigmaFByGroup(g, state_from, temp_from) * 
          getFluxByValue(i, g, state_from);
        
        fission_rate_to   += material->getNuSigmaFByGroup(g, state_to, temp_to) * 
          getFluxByValue(i, g, state_to);
      }
      
      for (int d=0; d < _num_delayed_groups; d++){
        double delayed_fraction = material->getDelayedFractionByGroup(d, state_from);
        double decay_constant   = material->getDecayConstantByGroup(d);
        double k1 = exp(- decay_constant * dt);
        double k2 = 1.0 - (1.0 - k1) / (decay_constant * dt);
        double k3 = k1 + (k2 - 1.0);
        double conc = material->getPrecursorConcByGroup(d, state_from);
        double new_conc = k1 * conc + k2 * delayed_fraction / decay_constant / _k_eff_0
          * fission_rate_to - k3 * delayed_fraction / decay_constant / _k_eff_0 * 
          fission_rate_from;

        material->setPrecursorConcByGroup(new_conc, d, state_to);
      }
    }
  }
}


void Solver::computeInitialPrecursorConcentrations(){

  for (int i=0; i < _num_shape_cells; i++){

    double fission_rate = 0.0;
    Material* material = _geometry->getMaterial(i);
    double temp = getTemperatureByValue(i, CURRENT);

    
    if (material->isFissionable()){
      for (int g=0; g < _num_energy_groups; g++){
        fission_rate += material->getNuSigmaFByGroup(g, CURRENT, temp) * 
          getFluxByValue(i, g, CURRENT);
      }

      for (int d=0; d < _num_delayed_groups; d++){
        double delayed_fraction = material->getDelayedFractionByGroup(d, CURRENT);
        double decay_constant = material->getDecayConstantByGroup(d);
        double conc = delayed_fraction / decay_constant / _k_eff_0 * fission_rate;
        material->setPrecursorConcByGroup(conc, d, CURRENT);
      }
    }
  }
}


void Solver::computePower(int state){

  for (int i=0; i < _num_shape_cells; i++){

    double cell_power = 0.0;
    Material* material = _geometry->getMaterial(i);
    double temp = getTemperatureByValue(i, CURRENT);
    
    if (material->isFissionable()){
      for (int g=0; g < _num_energy_groups; g++){
        cell_power += material->getSigmaFByGroup(g, state, temp) * 
          getFluxByValue(i, g, state) * material->getEnergyPerFission()
          * _geometry->getVolume(i);
      }
    }

    setPowerByValue(cell_power, i, state);
  }
}


void Solver::normalizeFlux(){

  /* Compute the initial power */
  double average_power = computeAveragePower(CURRENT);

  for (int i=0; i < _num_shape_cells; i++){

    for (int e=0; e < _num_energy_groups; e++){
      double old_flux = getFluxByValue(i, e, CURRENT);
      setFluxByValue(old_flux * _initial_power / average_power, i, e, CURRENT);
    }
  }
  
  /* Recompute power with normalized flux */
  computePower(CURRENT);
}


double Solver::computeAveragePower(int state){
  
  computePower(state);
  double average_power = 0.0;
  double fuel_volume = 0.0;

  for (int i=0; i < _num_shape_cells; i++){

    Material* material = _geometry->getMaterial(i);
    average_power += getPowerByValue(i, state);

    if (material->isFissionable())
      fuel_volume += _geometry->getVolume(i);
  }

  /* Divide cumulative power by fuel volume */
  average_power /= fuel_volume;

  return average_power;
}


double Solver::computePowerRMSError(int state_1, int state_2){
  
  computePower(state_1);
  computePower(state_2);
  Vector* power_1 = getPower(state_1);
  Vector* power_2 = getPower(state_2);
  double power_rmse = computeRMSE(power_1, power_2, false);
  return power_rmse;
}


void Solver::integrateTemperature(int state_from, int state_to){
  
  double dt = _clock->getTime(state_to) - _clock->getTime(state_from);

  #pragma omp parallel for
  for (int i=0; i < _num_shape_cells; i++){

    double fission_rate_from = 0.0;
    double fission_rate_to   = 0.0;
    double temp_from         = getTemperatureByValue(i, state_from);
    double temp_to           = getTemperatureByValue(i, state_to);    
    Material* material = _geometry->getMaterial(i);
    
    if (material->isFissionable()){
      for (int g=0; g < _num_energy_groups; g++){
        fission_rate_from += material->getSigmaFByGroup(g, state_from, temp_from) *
          getFluxByValue(i, g, state_from);
        fission_rate_to   += material->getSigmaFByGroup(g, state_to  , temp_to  ) *
          getFluxByValue(i, g, state_to);
      }
      
      temp_to = getTemperatureByValue(i, state_from) + dt * 0.5 *
        (fission_rate_from + fission_rate_to) * material->getTemperatureConversionFactor();
      setTemperatureByValue(temp_to, i, state_to);
    }
  }
}


void Solver::computeShape(int state, int state_flux, int state_amp){

  std::vector< std::vector<int> >* amp_to_shape = _geometry->getAmpToShapeMap();

  #pragma omp parallel for
  for (int i=0; i < _num_amp_cells; i++){

    std::vector<int>::iterator iter;
    for (iter = amp_to_shape->at(i).begin(); iter != amp_to_shape->at(i).end(); ++iter){
      
      for (int g=0; g < _num_energy_groups; g++){
        double flux = getFluxByValue(*iter, g, state_flux);
        double amp = getAmplitudeByValue(i, g, state_amp);
        setShapeByValue(flux / amp, *iter, g, state);
      }
    }
  }
}


void Solver::interpolateShape(int state, int state_forward, int state_backward){

  double dt = _clock->getTime(state_forward) - _clock->getTime(state_backward);
  double wt_bwd = (_clock->getTime(state_forward) - _clock->getTime(state)) / dt;
  double wt_fwd = (_clock->getTime(state) - _clock->getTime(state_backward)) / dt;

  #pragma omp parallel for
  for (int i=0; i < _num_shape_cells; i++){
    for (int g=0; g < _num_energy_groups; g++){
      double shape_fwd = getShapeByValue(i, g, state_forward);
      double shape_bwd = getShapeByValue(i, g, state_backward);
      setShapeByValue(shape_bwd * wt_bwd + shape_fwd * wt_fwd, i, g, state);
    }
  }  
}


void Solver::interpolateDifNonlinear(int state, int state_forward, int state_backward){

  double dt = _clock->getTime(state_forward) - _clock->getTime(state_backward);
  double wt_bwd = (_clock->getTime(state_forward) - _clock->getTime(state)) / dt;
  double wt_fwd = (_clock->getTime(state) - _clock->getTime(state_backward)) / dt;

  #pragma omp parallel for
  for (int i=0; i < _num_shape_cells; i++){
    for (int s=0; s < 6; s++){
      for (int g=0; g < _num_energy_groups; g++){
        double dif_nonlinear_fwd = getDifNonlinearByValue(i, g, s, state_forward);
        double dif_nonlinear_bwd = getDifNonlinearByValue(i, g, s, state_backward);
        setDifNonlinearByValue(dif_nonlinear_bwd * wt_bwd + dif_nonlinear_fwd * wt_fwd,
                               i, g, s, state);
      }
    }  
  }
}


void Solver::reconstructFlux(int state, int state_shape, int state_amp){

  std::vector< std::vector<int> >* amp_to_shape = _geometry->getAmpToShapeMap();

  #pragma omp parallel for
  for (int i=0; i < _num_amp_cells; i++){

    std::vector<int>::iterator iter;
    for (iter = amp_to_shape->at(i).begin(); iter != amp_to_shape->at(i).end(); ++iter){
      
      for (int g=0; g < _num_energy_groups; g++){
        double shape = getShapeByValue(*iter, g, state_shape);
        double amp = getAmplitudeByValue(i, g, state_amp);
        setFluxByValue(shape * amp, *iter, g, state);
      }
    }
  }
}


void Solver::setEndTime(double time){
  _clock->setEndTime(time);
}


void Solver::setInnerTimeStepSize(double time){
  _clock->setDtInner(time);
}


void Solver::setOuterTimeStepSize(double time){
  _clock->setDtOuter(time);
}


void Solver::computeDiffusionCoefficients(int state){

  int nx = _geometry->getNumXAmp();
  int ny = _geometry->getNumYAmp();
  int nz = _geometry->getNumZAmp();
  int ng = _num_energy_groups;
  double width = _geometry->getWidth() / nx;
  double height = _geometry->getHeight() / ny;
  double depth = _geometry->getDepth() / nz;
  std::vector< std::vector<int> >* amp_to_shape = _geometry->getAmpToShapeMap();

  /* Reconstruct the fine mesh flux at t=state */
  reconstructFlux(state, state, state);

  /* Allocate array of diffusion coefficients */
  double* dif_coefs = new double[_num_amp_cells*ng];

  double trans_tally, rxn_tally, volume, dif_coef, flux, temp;
  Material* material;
  int cell_amp, cell_shape;

  /* Condense the diffusion coefficients */
  for (int i=0; i < _num_amp_cells; i++){
    for (int e=0; e < ng; e++){
      
      trans_tally = 0.0;
      rxn_tally = 0.0;
      
      /* Loop over shape cells in this amp cell */
      std::vector<int>::iterator iter;
      for (iter = amp_to_shape->at(i).begin(); iter != amp_to_shape->at(i).end(); ++iter){
        
        cell_shape = *iter;
        material = _geometry->getMaterial(cell_shape);
        temp = getTemperatureByValue(cell_shape, state);
        volume = _geometry->getVolume(cell_shape);
        dif_coef = material->getDifCoefByGroup(e, state, temp);
        flux = getFluxByValue(cell_shape, e, state);
        
        trans_tally += 1.0 / (3.0*dif_coef) * flux * volume;
        rxn_tally += flux * volume;
      }
      
      dif_coefs[i*ng+e] = 1.0 / (3.0 * (trans_tally / rxn_tally));
    }
  }
    
  /* Compute the linear and nonlinear diffusion coefficients */
  for (int z=0; z < nz; z++){
    for (int y=0; y < ny; y++){
      for (int x=0; x < nx; x++){
        
        int cell = z*nx*ny + y*nx + x;

        int sense;
        double length, length_perpen;
        double current, dif_linear, dif_nonlinear, flux_next, dif_coef_next;
        
        for (int s=0; s < 6; s++){
          
          int cell_next = _geometry->getNeighborAmpCell(x, y, z, s);
          
          if (s == SURFACE_X_MIN ||
              s == SURFACE_Y_MIN ||
              s == SURFACE_Z_MIN)
            sense = -1;
          else
            sense = 1;
          
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
          
          for (int g=0; g < ng; g++){
            
            dif_coef = dif_coefs[cell*ng+g];
            current = getCurrentByValue(cell, g, s, state);
            flux = getAmplitudeByValue(cell, g, state);
            
            if (cell_next == -1){
              dif_linear = 2 * dif_coef / length_perpen / 
                (1 + 4 * dif_coef / length_perpen);
              dif_nonlinear = (sense * dif_linear * flux - current / length) / flux;
              dif_linear *= (1.0 - _geometry->getBoundary(s));
              dif_nonlinear *= (1.0 - _geometry->getBoundary(s));
            }
            else{
              flux_next = getAmplitudeByValue(cell_next, g, state);
              dif_coef_next = dif_coefs[cell_next*ng + g];
              dif_linear = 2 * dif_coef * dif_coef_next / (length_perpen * dif_coef +
                                                           length_perpen * dif_coef_next);
              dif_nonlinear = - (sense * dif_linear * (flux_next - flux) + current / length)
                / (flux_next + flux);
              
              if (dif_nonlinear > dif_linear){
                if (sense == -1){
                  if (dif_nonlinear > 0.0){
                    dif_linear = - current / (2 * flux);
                    dif_nonlinear = - current / (2 * flux);
                  }
                  else{
                    dif_linear = current / (2 * flux_next);
                    dif_nonlinear = - current / (2 * flux_next);                  
                  }
                }
                else{
                  if (dif_nonlinear > 0.0){
                    dif_linear = - current / (2 * flux_next);
                    dif_nonlinear = - current / (2 * flux_next);
                  }
                  else{
                    dif_linear = current / (2 * flux);
                    dif_nonlinear = - current / (2 * flux);
                  }
                }
              }
            }
            
            setDifLinearByValue(dif_linear, cell, g, s, state);
            setDifNonlinearByValue(dif_nonlinear, cell, g, s, state);          
          }
        }
      }
    }
  }

  delete [] dif_coefs;
}


void Solver::computeFrequency(){

  int nx = _geometry->getNumXAmp();
  int ny = _geometry->getNumYAmp();
  int nz = _geometry->getNumZAmp();
  int ng = _num_energy_groups;
  double width = _geometry->getWidth() / nx;
  double height = _geometry->getHeight() / ny;
  double depth = _geometry->getDepth() / nz;
  double dt = _clock->getDtInner();
  double mult;
  
  std::vector< std::vector<int> >* amp_to_shape = _geometry->getAmpToShapeMap();
 
  reconstructFlux(CURRENT, CURRENT, CURRENT);

  _frequency[CURRENT]->clear();
  
  for (int z=0; z < nz; z++){
    for (int y=0; y < ny; y++){
      for (int x=0; x < nx; x++){

        int cell_amp = z*nx*ny + y*nx+x;
            
        for (int e=0; e < ng; e++){

          double freq = 0.0;
          double amp         = getAmplitudeByValue(cell_amp, e, CURRENT);
          double amp_prev    = getAmplitudeByValue(cell_amp, e, PREVIOUS_IN);
          
          /* Loop over shape cells in this amp cell */
          std::vector<int>::iterator iter;
          for (iter = amp_to_shape->at(cell_amp).begin(); iter != amp_to_shape->at(cell_amp).end();
               ++iter){
            
            int cell_shape = *iter;
            Material* material = _geometry->getMaterial(cell_shape);
            double temp        = getTemperatureByValue(cell_shape, CURRENT);
            double volume      = _geometry->getVolume(cell_shape);
            double v           = material->getVelocityByGroup(e, CURRENT, temp);
            double shape       = getShapeByValue(cell_shape, e, CURRENT);
            double shape_prev  = getShapeByValue(cell_shape, e, PREVIOUS_IN);
            double flux        = getFluxByValue(cell_shape, e, CURRENT);
            double chi         = material->getChiByGroup(e, CURRENT, temp);
            double sig_a       = material->getSigmaAByGroup(e, CURRENT, temp);
            double dif_coef    = material->getDifCoefByGroup(e, CURRENT, temp);
            double delay_total = material->getDelayedFractionTotal(CURRENT);
            double flux_g, sig_s, nu_sig_f;
            
            mult        = - v / (volume * flux);

            freq = - 1.0 / dt;
            
            /* Time absorption term on the diagonal */
            freq += mult / dt / v * volume * shape_prev * amp;
            
            /* Delayed neutron precursors */
            if (material->isFissionable()){
              for (int d=0; d < _num_delayed_groups; d++){
                double decay = material->getDecayConstantByGroup(d);
                double conc = material->getPrecursorConcByGroup(d, CURRENT);
                freq += mult * chi * volume * decay * conc;
              }
            }
          
            /* Absorption term on the diagonal */
            freq -= mult * sig_a * volume * flux;
            
            /* Buckling term on the diagonal */
            freq -= mult * dif_coef * volume * _buckling * flux;
            
            /* Outscattering term on diagonal */
            for (int g=0; g < ng; g++){
              if (e != g){
                freq -= mult * volume * flux
                  * material->getSigmaSByGroup(e, g, CURRENT, temp);
              }
            }
            
            /* Fission terms */
            for (int g=0; g < ng; g++){
              nu_sig_f = material->getNuSigmaFByGroup(g, CURRENT, temp);
              flux_g = getFluxByValue(cell_amp, g, CURRENT);              
              freq += mult * (1-delay_total) * chi / _k_eff_0 * volume * nu_sig_f * flux_g;
            }
            
            /* Inscattering term on off diagonals */
            for (int g=0; g < ng; g++){
              if (e != g){
                sig_s = material->getSigmaSByGroup(g, e, CURRENT, temp);
                flux_g = getFluxByValue(cell_amp, g, CURRENT);
                freq += mult * sig_s * flux_g * volume;
              }
            }
          }

          /* X_MIN SURFACE */
          
          /* Transport term on diagonal */
          freq -= mult * (getDifLinearByValue(cell_amp, e, SURFACE_X_MIN, CURRENT) +
                          getDifNonlinearByValue(cell_amp, e, SURFACE_X_MIN, CURRENT))
            * height * depth * amp;

          /* Transport term on off diagonals */
          if (x != 0){
            freq += mult * (getDifLinearByValue(cell_amp, e, SURFACE_X_MIN, CURRENT) -
                            getDifNonlinearByValue(cell_amp, e, SURFACE_X_MIN, CURRENT))
              * height * depth * getAmplitudeByValue(cell_amp-1, e, CURRENT);
          }

          /* X_MAX SURFACE */
          
          /* Transport term on diagonal */
          freq -= mult * (getDifLinearByValue(cell_amp, e, SURFACE_X_MAX, CURRENT) -
                          getDifNonlinearByValue(cell_amp, e, SURFACE_X_MAX, CURRENT))
            * height * depth * amp;
          
          /* Transport term on off diagonals */
          if (x != nx - 1){
            freq += mult * (getDifLinearByValue(cell_amp, e, SURFACE_X_MAX, CURRENT) +
                            getDifNonlinearByValue(cell_amp, e, SURFACE_X_MAX, CURRENT))
              * height * depth * getAmplitudeByValue(cell_amp + 1, e, CURRENT);
          }
          
          /* Y_MIN SURFACE */
          
          /* Transport term on diagonal */
          freq -= mult * (getDifLinearByValue(cell_amp, e, SURFACE_Y_MIN, CURRENT) +
                          getDifNonlinearByValue(cell_amp, e, SURFACE_Y_MIN, CURRENT))
            * width * depth * amp;          
          
          /* Transport term on off diagonals */
          if (y != 0){
            freq += mult * (getDifLinearByValue(cell_amp, e, SURFACE_Y_MIN, CURRENT) - 
                            getDifNonlinearByValue(cell_amp, e, SURFACE_Y_MIN, CURRENT))
              * width * depth * getAmplitudeByValue(cell_amp - nx, e, CURRENT);
          }
          
          /* Y_MAX SURFACE */
          
          /* Transport term on diagonal */
          freq -= mult * (getDifLinearByValue(cell_amp, e, SURFACE_Y_MAX, CURRENT) -
                          getDifNonlinearByValue(cell_amp, e, SURFACE_Y_MAX, CURRENT))
            * width * depth * amp;
          
          /* Transport term on off diagonals */
          if (y != ny - 1){
            freq += mult * (getDifLinearByValue(cell_amp, e, SURFACE_Y_MAX, CURRENT) +
                            getDifNonlinearByValue(cell_amp, e, SURFACE_Y_MAX, CURRENT))
              * width * depth * getAmplitudeByValue(cell_amp + nx, e, CURRENT);
          }

          /* Z_MIN SURFACE */
          
          /* Transport term on diagonal */
          freq -= mult * (getDifLinearByValue(cell_amp, e, SURFACE_Z_MIN, CURRENT) +
                          getDifNonlinearByValue(cell_amp, e, SURFACE_Z_MIN, CURRENT))
            * width * height * amp;
          
          /* Transport term on off diagonals */
          if (z != 0){
            freq += mult * (getDifLinearByValue(cell_amp, e, SURFACE_Z_MIN, CURRENT) -
                            getDifNonlinearByValue(cell_amp, e, SURFACE_Z_MIN, CURRENT))
              * width * height * getAmplitudeByValue(cell_amp - nx*ny, e, CURRENT);
          }
          
          /* Z_MAX SURFACE */
          
          /* Transport term on diagonal */
          freq -= mult * (getDifLinearByValue(cell_amp, e, SURFACE_Z_MAX, CURRENT) -
                          getDifNonlinearByValue(cell_amp, e, SURFACE_Z_MAX, CURRENT))
            * width * height * amp;
          
          /* Transport term on off diagonals */
          if (z != nz - 1){
            freq += mult * (getDifLinearByValue(cell_amp, e, SURFACE_Z_MAX, CURRENT) +
                            getDifNonlinearByValue(cell_amp, e, SURFACE_Z_MAX, CURRENT))
              * width * height * getAmplitudeByValue(cell_amp + nx*ny, e, CURRENT);
          }

          setFrequencyByValue(freq, cell_amp, e, CURRENT);
          log_printf(NORMAL, "cell: %i, group: %i, freq: %f", cell_amp, e, freq);
        }
      }
    }
  }
}


void Solver::setBuckling(double buckling){
  _buckling = buckling;
}


double Solver::getBuckling(){
  return _buckling;
}


double Solver::getKeff0(){
  return _k_eff_0;
}


void Solver::setInitialPower(double power){
  _initial_power = power;
}


int Solver::getMethod(){
  return _method;
}


void Solver::setMethod(int method){
  _method = method;
}


void Solver::initializeClock(){

  log_printf(NORMAL, "num shape cells: %i", _num_shape_cells);
  
  for (int i=0; i < _num_shape_cells; i++)
    _geometry->getMaterial(i)->setClock(_clock);
}


Geometry* Solver::getGeometry(){
  return _geometry;
}
