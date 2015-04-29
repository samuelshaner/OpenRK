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
  _amp_matrix = new double*[_num_amp_cells];
  _amp_source = new double[_num_amp_cells*_num_energy_groups];
  
  for (int i=0; i < _num_amp_cells; i++)
    _amp_matrix[i] = new double[_num_energy_groups*(_num_energy_groups+6)];

  /* Allocate memory for field variables */
  for (int t=0; t < 8; t++){
    double* temperature   = new double[_num_shape_cells];
    double* flux          = new double[_num_shape_cells*_num_energy_groups];
    double* shape         = new double[_num_shape_cells*_num_energy_groups];
    double* amplitude     = new double[_num_amp_cells*_num_energy_groups];
    double* power         = new double[_num_shape_cells];
    double* current       = new double[_num_amp_cells*_num_energy_groups*6];
    double* dif_linear    = new double[_num_amp_cells*_num_energy_groups*6];
    double* dif_nonlinear = new double[_num_amp_cells*_num_energy_groups*6];
    double* frequency     = new double[_num_amp_cells*_num_energy_groups];

    memset(temperature, 300.0, sizeof(double) * _num_shape_cells);
    memset(flux,          1.0, sizeof(double) * _num_shape_cells * _num_energy_groups);
    memset(shape,         0.0, sizeof(double) * _num_shape_cells * _num_energy_groups);
    memset(amplitude,     1.0, sizeof(double) * _num_amp_cells * _num_energy_groups);
    memset(power,         0.0, sizeof(double) * _num_shape_cells);
    memset(current,       0.0, sizeof(double) * _num_amp_cells * _num_energy_groups * 6);
    memset(dif_linear,    0.0, sizeof(double) * _num_amp_cells * _num_energy_groups * 6);
    memset(dif_nonlinear, 0.0, sizeof(double) * _num_amp_cells * _num_energy_groups * 6);
    memset(frequency,     0.0, sizeof(double) * _num_amp_cells * _num_energy_groups);

    _temperature[t] = temperature;
    _flux[t] = flux;
    _shape[t] = shape;
    _amplitude[t] = amplitude;
    _power[t] = power;
    _current[t] = current;
    _dif_linear[t] = dif_linear;
    _dif_nonlinear[t] = dif_nonlinear;
    _frequency[t] = frequency;
  }
}


Solver::~Solver(){

  if (_amp_matrix != NULL){
    for (int i=0; i < _num_amp_cells; i++)
      delete [] _amp_matrix[i];

    delete [] _amp_matrix;
  }    

  if (_amp_source != NULL){
    delete [] _amp_source;
  }    
}


double** Solver::getAmpMatrix(){
  return _amp_matrix;
}


double* Solver::getAmpShape(){
  return _amp_shape;
}


void Solver::generateAmplitudeMatrix(double wt){

  int nx = _geometry->getNumXAmp();
  int ny = _geometry->getNumYAmp();
  int nz = _geometry->getNumZAmp();
  int ng = _num_energy_groups;
  double width = _geometry->getWidth() / nx;
  double height = _geometry->getHeight() / ny;
  double depth = _geometry->getDepth() / nz;
  double dt = _clock->getDtInner();
  int row_amp, row_shape, diag;
  int ng = _num_energy_groups;

  double* shape = getShape(CURRENT);
  double* shape_prev = getShape(PREVIOUS_IN);
  double* amp_prev = getAmplitude(PREVIOUS_IN);

  std::vector< std::vector<int> > amp_to_shape = _geometry->getAmpToShapeMap();
 
  #pragma omp parallel for
  for (int i = 0; i < nx*ny*nz; i++){
    for (int g=0; g < ng; g++)
      _amp_source[i*ng+g] = 0.0;
    
    for (int g=0; g < ng*(ng+6); g++)
      _amp_matrix[i][g] = 0.0;
  }
  
  for (int z=0; z < nz; z++){
#pragma omp parallel for private(row_amp, row_shape, diag)
    for (int y=0; y < ny; y++){
      for (int x=0; x < nx; x++){
        int cell_amp = z*nx*ny + y*nx+x;
            
        /* Loop over energy groups */
        for (int e=0; e < ng; e++){

          diag = e * (ng+6) + e + 3;
          row_amp = cell_amp*ng + e;
          
          /* Loop over shape cells in this amp cell */
          std::vector<int>::iterator iter;
          for (iter = amp_to_shape.begin(); iter != amp_to_shape.end(); ++iter){
            
            int cell_shape = *iter;
            Material* material = _geometry->getMaterial(cell_shape);
            double temp = getTemperatureByValue(cell_shape, CURRENT);
            double temp_prev = getTemperatureByValue(cell_shape, PREVIOUS_IN);
            double volume = _geometry->getVolume(cell_shape);
            row_shape = cell_shape*ng + e;

            /* Time absorption term on the diagonal */
            if (method == IQS){
              
              _amp_matrix[cell_amp][diag] += 2.0 / dt / 
                _material->getVelocityByGroup(e, CURRENT, temp) * volume * 
                shape[row_shape];
              
              _amp_matrix[cell_amp][diag] -= 1.0 / dt / 
                material->getVelocityByGroup(e, CURRENT, temp) * volume * 
                shape_prev[row_shape];
          
              _amp_source[row_amp] += shape[row_shape] / dt / 
                material->getVelocityByGroup(e, PREVIOUS_IN, temp_prev) * 
                volume * amp_prev[row_amp];
            }
            else if (method == THETA){
              _amp_matrix[cell_amp][diag] += 1.0 / dt / 
                material->getVelocityByGroup(e, CURRENT, temp) * 
                volume * shape[row_shape];

              _amp_source[row_amp] += shape_prev[row_shape] / dt / 
                material->getVelocityByGroup(e, PREVIOUS_IN, temp_prev) * 
                volume * amp_prev[row_amp];
            }
                    
            /* Delayed neutron precursors */
            for (int d=0; d < _num_delayed_groups; d++){
              _amp_source[row_amp] += wt * material->getChiByGroup(e, CURRENT, temp)
                * material->getDecayConstantByGroup(d) * 
                material->getPrecursorConcByGroup(d, CURRENT) * volume;

              _amp_source[row_amp] += (1-wt) *
                material->getChiByGroup(e, PREVIOUS_IN, temp_prev)
                * material->getDecayConstantByGroup(d) * 
                material->getPrecursorConcByGroup(d, PREVIOUS_IN) * volume;
            }
          
            /* Absorption term on the diagonal */
            _amp_matrix[cell_amp][diag] += wt * 
              material->getSigmaAByGroup(e, CURRENT, temp) * volume * shape[row_shape];

            _amp_source[row_amp] -= (1-wt) * 
              material->getSigmaAByGroup(e, PREVIOUS_IN, temp_prev) *
              shape_prev[row_shape] * amp_prev[row_amp] * volume;
            
            /* Buckling term on the diagonal */
            _amp_matrix[cell_amp][diag] += wt * 
              material->getDifCoefByGroup(e, CURRENT, temp) * 
              volume * _buckling * shape[row_shape];

            _amp_source[row_amp] -= (1-wt) * 
              material->getDifCoefByGroup(e, PREVIOUS_IN, temp_prev) * 
              shape_prev[row_shape] * amp_prev[row_amp] * volume * _buckling;
            
            /* Outscattering term on diagonal */
            for (int g=0; g < ng; g++){
              if (e != g){
                _amp_matrix[cell_amp][diag] += wt * 
                  material->getSigmaSByGroup(e, g, CURRENT, temp) * volume * 
                  shape[row_shape];

                _amp_source[row_amp] -= (1-wt) * 
                  material->getSigmaSByGroup(e, g, PREVIOUS_IN, temp_prev)
                  * shape_prev[row_shape] * amp_prev[row_amp] * volume;
              }
            }
            
            /* Fission terms */
            for (int g=0; g < ng; g++){
              _amp_matrix[cell_amp][e * (ng+6) + g + 3] -= wt * 
                (1-material->getDelayedFractionTotal(CURRENT)) * 
                material->getChiByGroup(e, CURRENT, temp) *
                material->getNuSigmaFByGroup(g, CURRENT, temp) / 
                _k_eff_0 * volume * shape[cell_shape*ng+g];

              _amp_source[row_amp] += (1-wt) * 
                (1-material->getDelayedFractionTotal(PREVIOUS_IN)) * 
                material->getChiByGroup(e, PREVIOUS_IN, temp_prev) *
                material->getNuSigmaFByGroup(g, PREVIOUS_IN, temp_prev) / _k_eff_0 *
                shape_prev[cell_shape*ng+g] * 
                amp_prev[cell_amp*ng+g] * volume;
            }
            
            /* Inscattering term on off diagonals */
            for (int g=0; g < ng; g++){
              if (e != g){
                _amp_matrix[cell_amp][e * (ng+6) + g + 3] -= 
                  wt * material->getSigmaSByGroup(g, e, CURRENT, temp) * 
                  volume * shape[cell_shape*ng+g];

                _amp_source[row_amp] += (1-wt) * 
                  material->getSigmaSByGroup(g, e, PREVIOUS_IN, temp_prev) *
                  shape_prev[cell_shape*ng+g] * 
                  amp_prev[cell_amp*ng+g] * volume;
              }
            }
          }          
           
          /* LEFT SURFACE */
          
          /* Transport term on diagonal */
          _amp_matrix[cell_amp][diag] += wt *
            (getDifLinearByValue(cell_amp, e, 0, CURRENT) +
             getDifNonlinearByValue(cell_amp, e, 0, CURRENT)) * height * depth;

          _amp_source[row_amp] -= (1-wt) *
            (getDifLinearByValue(cell_amp, e, 0, PREVIOUS_IN) +
             getDifNonlinearByValue(cell_amp, e, 0, PREVIOUS_IN)) * height * depth * 
            amp_prev[row_amp];          
          
          /* Transport term on off diagonals */
          if (x != 0){
            _amp_matrix[cell_amp][e * (ng+6)] -= wt *
              (_amp_mesh->getDifLinearByValue(cell_amp, e, 0, CURRENT) -
               _amp_mesh->getDifNonlinearByValue(cell_amp, e, 0, CURRENT)) * height * depth;

            _amp_source[row_amp] += (1-wt) *
              (_amp_mesh->getDifLinearByValue(cell_amp, e, 0, PREVIOUS_IN) -
               _amp_mesh->getDifNonlinearByValue(cell_amp, e, 0, PREVIOUS_IN)) * height * 
              depth * amp_prev[(cell_amp-1)*ng+e];
          }          

          /* RIGHT SURFACE */
          
          /* Transport term on diagonal */
          _amp_matrix[cell_amp][diag] += wt *
            (getDifLinearByValue(cell_amp, e, 3, CURRENT) -
             getDifNonlinearByValue(cell_amp, e, 3, CURRENT)) * height * depth;

          _amp_source[row_amp] -= (1-wt) *
            (getDifLinearByValue(cell_amp, e, 3, PREVIOUS_IN) -
             getDifNonlinearByValue(cell_amp, e, 3, PREVIOUS_IN)) * height * depth * 
            amp_prev[row_amp];
          
          /* Transport term on off diagonals */
          if (x != nx - 1){
            _amp_matrix[cell_amp][e * (ng+6) + ng + 3] -= wt *
              (getDifLinearByValue(cell_amp, e, 3, CURRENT) +
               getDifNonlinearByValue(cell_amp, e, 3, CURRENT)) * height * depth;

            _amp_source[row_amp] += (1-wt) *
              (getDifLinearByValue(cell_amp, e, 3, PREVIOUS_IN) +
               getDifNonlinearByValue(cell_amp, e, 3, PREVIOUS_IN)) * height * depth * 
              amp_prev[(cell_amp+1)*ng+e];
          }
          
          /* BACK SURFACE */
          
          /* Transport term on diagonal */
          _amp_matrix[cell_amp][diag] += wt *
            (getDifLinearByValue(cell_amp, e, 1, CURRENT) +
             getDifNonlinearByValue(cell_amp, e, 1, CURRENT)) * width * depth;

          _amp_source[row_amp] -= (1-wt) *
            (getDifLinearByValue(cell_amp, e, 1, PREVIOUS_IN) +
             getDifNonlinearByValue(cell_amp, e, 1, PREVIOUS_IN)) * width * depth * 
            amp_prev[row_amp];
          
          
          /* Transport term on off diagonals */
          if (y != 0){
            _amp_matrix[cell_amp][e * (ng+6) + 1] -= wt *
              (getDifLinearByValue(cell_amp, e, 1, CURRENT) -
               getDifNonlinearByValue(cell_amp, e, 1, CURRENT)) * width * depth;

            _amp_source[row_amp] += (1-wt) *
              (getDifLinearByValue(cell_amp, e, 1, PREVIOUS_IN) -
               getDifNonlinearByValue(cell_amp, e, 1, PREVIOUS_IN)) * width * depth * 
              amp_prev[(cell_amp-nx)*ng+e];
          }
          
          /* FRONT SURFACE */
          
          /* Transport term on diagonal */
          _amp_matrix[cell_amp][diag] += wt *
            (getDifLinearByValue(cell_amp, e, 4, CURRENT) -
             getDifNonlinearByValue(cell_amp, e, 4, CURRENT)) * width * depth;

          _amp_source[row_amp] -= (1-wt) *
            (getDifLinearByValue(cell_amp, e, 4, PREVIOUS_IN) -
             getDifNonlinearByValue(cell_amp, e, 4, PREVIOUS_IN)) * width * depth * 
            amp_prev[row_amp];

          /* Transport term on off diagonals */
          if (y != ny - 1){
            _amp_matrix[cell_amp][e * (ng+6) + ng + 4] -= wt *
              (getDifLinearByValue(cell_amp, e, 4, CURRENT) +
               getDifNonlinearByValue(cell_amp, e, 4, CURRENT)) * width * depth;

            _amp_source[row_amp] += (1-wt) *
              (getDifLinearByValue(cell_amp, e, 4, PREVIOUS_IN) +
               getDifNonlinearByValue(cell_amp, e, 4, PREVIOUS_IN)) * width * depth * 
              amp_prev[(cell_amp+nx)*ng+e];
          }

          /* BOTTOM SURFACE */
          
          /* Transport term on diagonal */
          _amp_matrix[cell_amp][diag] += wt *
            (getDifLinearByValue(cell_amp, e, 2, CURRENT) +
             getDifNonlinearByValue(cell_amp, e, 2, CURRENT)) * width * height;

          _amp_source[row_amp] -= (1-wt) *
            (getDifLinearByValue(cell_amp, e, 2, PREVIOUS_IN) +
             getDifNonlinearByValue(cell_amp, e, 2, PREVIOUS_IN)) * width * height * 
            amp_prev[row_amp];

          /* Transport term on off diagonals */
          if (z != 0){
            _amp_matrix[cell_amp][e * (ng+6) + 2] -= wt *
              (getDifLinearByValue(cell_amp, e, 2, CURRENT) -
               getDifNonlinearByValue(cell_amp, e, 2, CURRENT)) * width * height;

            _amp_source[row_amp] += (1-wt) *
              (getDifLinearByValue(cell_amp, e, 2, PREVIOUS_IN) -
               getDifNonlinearByValue(cell_amp, e, 2, PREVIOUS_IN)) * width * height * 
              amp_prev[(cell_amp-nx)*ng+e];
          }
          
          /* TOP SURFACE */
          
          /* Transport term on diagonal */
          _amp_matrix[cell_amp][diag] += wt *
            (getDifLinearByValue(cell_amp, e, 5, CURRENT) -
             getDifNonlinearByValue(cell_amp, e, 5, CURRENT)) * width * height;

          _amp_source[row_amp] -= (1-wt) *
            (getDifLinearByValue(cell_amp, e, 5, PREVIOUS_IN) -
             getDifNonlinearByValue(cell_amp, e, 5, PREVIOUS_IN)) * width * height * 
            amp_prev[row_amp];

          /* Transport term on off diagonals */
          if (z != nz - 1){
            _amp_matrix[cell_amp][e * (ng+6) + ng + 5] -= wt *
              (getDifLinearByValue(cell_amp, e, 5, CURRENT) +
               getDifNonlinearByValue(cell_amp, e, 5, CURRENT)) * width * height;

            _amp_source[row_amp] += (1-wt) *
              (getDifLinearByValue(cell_amp, e, 5, PREVIOUS_IN) +
               getDifNonlinearByValue(cell_amp, e, 5, PREVIOUS_IN)) * width * height * 
              amp_prev[(cell_amp+nx)*ng+e];
          }
        }
      }
    }
  }
}


double* Solver::getTemperature(int time){
  return _temperature[time];
}


double* Solver::getFlux(int time){
  return _flux[time];
}


double* Solver::getShape(int time){
  return _shape[time];
}


double* Solver::getAmplitude(int time){
  return _amplitude[time];
}


double* Solver::getPower(int time){
  return _power[time];
}


double* Solver::getCurrent(int time){
  return _current[time];
}


double* Solver::getDifLinear(int time){
  return _dif_linear[time];
}


double* Solver::getDifNonlinear(int time){
  return _dif_nonlinear[time];
}


double Solver::getTemperatureByValue(int cell, int time){
  return _temperature[time][cell];
}


double Solver::getFluxByValue(int cell, int group, int time){
  return _flux[time][cell*_num_energy_groups+group];
}


double Solver::getShapeByValue(int cell, int group, int time){
  return _shape[time][cell*_num_energy_groups+group];
}


double Solver::getAmplitudeByValue(int cell, int group, int time){
  return _amplitude[time][cell*_num_energy_groups+group];
}


double Solver::getPowerByValue(int cell, int group, int time){
  return _power[time][cell];
}


double Solver::getCurrentByValue(int cell, int group, int side, int time){
  return _current[time][(cell*6+side) * _num_energy_groups+group];
}


double Solver::getDifLinearByValue(int cell, int group, int side, int time){
  return _dif_linear[time][(cell*6+side) * _num_energy_groups+group];
}


double Solver::getDifNonlinearByValue(int cell, int group, int side, int time){
  return _dif_nonlinear[time][(cell*6+side) * _num_energy_groups+group];
}


double Solver::getFrequencyByValue(int cell, int group, int time){
  return _frequency[time][cell * _num_energy_groups+group];
}


void Solver::copyTemperature(int time_from, int time_to){
  std::copy(_temperature[time_from], _temperature[time_from] + _num_shape_cells, _temperature[time_to]);
}


void Solver::copyFlux(int time_from, int time_to){
  std::copy(_flux[time_from], _flux[time_from] + _num_shape_cells*_num_energy_groups, _flux[time_to]);
}


void Solver::copyShape(int time_from, int time_to){
  std::copy(_shape[time_from], _shape[time_from] + _num_shape_cells*_num_energy_groups, _shape[time_to]);
}


void Solver::copyAmplitude(int time_from, int time_to){
  std::copy(_amplitude[time_from], _amplitude[time_from] + _num_amp_cells*_num_energy_groups, _amplitude[time_to]);
}


void Solver::copyPower(int time_from, int time_to){
  std::copy(_power[time_from], _power[time_from] + _num_shape_cells, _power[time_to]);
}


void Solver::copyCurrent(int time_from, int time_to){
  std::copy(_current[time_from], _current[time_from] + _num_amp_cells*_num_energy_groups*6, _current[time_to]);
}


void Solver::copyDifLinear(int time_from, int time_to){
  std::copy(_dif_linear[time_from], _dif_linear[time_from] + _num_amp_cells*_num_energy_groups*6, _dif_linear[time_to]);
}


void Solver::copyDifNonlinear(int time_from, int time_to){
  std::copy(_dif_nonlinear[time_from], _dif_nonlinear[time_from] + _num_amp_cells*_num_energy_groups*6, _dif_nonlinear[time_to]);
}


void Solver::copyFrequency(int time_from, int time_to){
  std::copy(_frequency[time_from], 
            _frequency[time_from] + _num_amp_cells*_num_energy_groups, _frequency[time_to]);
}


void Solver::copyPrecursors(int time_from, int time_to){
  for (int i=0; i < _num_shape_cells; i++){
    Material* material = _materials[i];
    for (int e=0; e < _num_energy_groups; e++){
      double conc = material->getPrecursorConcByGroup(e, time_from);
      material->setPrecursorConcByGroup(conc, e, time_to); 
    }
  }
}

void Solver::copyFieldVariables(int time_from, int time_to){

  copyTemperature(time_from, time_to);
  copyFlux(time_from, time_to);
  copyShape(time_from, time_to);
  copyAmplitude(time_from, time_to);
  copyPower(time_from, time_to);
  copyCurrent(time_from, time_to);
  copyDifLinear(time_from, time_to);
  copyDifNonlinear(time_from, time_to);
  copyFrequency(time_from, time_to);
  copyPrecursors(time_from, time_to);
}


void Solver::broadcastToAll(int time_from){

  for (int t=0; t < 8; t++)
    copyFieldVariables(time_from, t);

}


void Solver::setFrequencyByValue(double value, int cell, int group, int time){
  _frequency[time][cell*_num_energy_groups+group] = value;
}


void Solver::setPowerByValue(double value, int cell, int time){
  _power[time][cell] = value;
}


void Solver::setAmplitudeByValue(double value, int cell, int group, int time){
  _amplitude[time][cell*_num_energy_groups+group] = value;
}


void Solver::setShapeByValue(double value, int cell, int group, int time){
  _shape[time][cell*_num_energy_groups+group] = value;
}


void Solver::setTemperatureByValue(double value, int cell, int time){
  _temperature[time][cell] = value;
}


void Solver::setFluxByValue(double value, int cell, int group, int time){
  _flux[time][cell*_num_energy_groups+group] = value;
}


void Solver::setCurrentByValue(double value, int cell, int group, int side, int time){
  _current[time][(cell*6+side) * _num_energy_groups+group] = value;
}


void Solver::setDifLinearByValue(double value, int cell, int group, int side, int time){
  _dif_linear[time][(cell*6+side) * _num_energy_groups+group] = value;
}


void Solver::setDifNonlinearByValue(double value, int cell, int group, int side, int time){
  _dif_nonlinear[time][(cell*6+side) * _num_energy_groups+group] = value;
}


void Solver::integratePrecursorConcentrations(int time_from, int time_to){

  double* temps_from = getTemperature(time_from);
  double* temps_to = getTemperature(time_to);
  double dt = _clock->getTime(time_to) - _clock->getTime(time_from);

  #pragma omp parallel for
  for (int i=0; i < _num_shape_cells; i++){

    double fission_rate_from = 0.0;
    double fission_rate_to = 0.0;
    Material* material = _geometry->getMaterial(i);

    if (material->isFissionable()){
      for (int g=0; g < _num_energy_groups; g++){
        fission_rate_from += material->getNuSigmaFByGroup(g, time_from, temps_from[i]) * 
          getFluxByValue(i, g, time_from);
        
        fission_rate_to += material->getNuSigmaFByGroup(g, time_to, temps_to[i]) * 
          getFluxByValue(i, g, time_to);
      }
      
      for (int d=0; d < _num_delayed_groups; d++){
        double delayed_fraction = material->getDelayedFractionByGroup(d, time_from);
        double decay_constant = material->getDecayConstantByGroup(d, time_from);
        double k1 = exp(- decay_constant * dt);
        double k2 = 1.0 - (1.0 - k1) / (decay_constant * dt);
        double k3 = k1 + (k2 - 1.0);
        double conc = material->getPrecursorConcByGroup(d, time_from);
        double new_conc = k1 * conc + k2 * delayed_fraction / decay_constant / _k_eff_0
          * fission_rate_to - k3 * delayed_fraction / decay_constant / _k_eff_0 * 
          fission_rate_from;

        material->setPrecursorConcByGroup(new_conc, d, time_to);
      }
    }
  }
}


void Solver::computeInitialPrecursorConcentrations(){

  double* temps = getTemperature(CURRENT);

  #pragma omp parallel for
  for (int i=0; i < _num_shape_cells; i++){

    double fission_rate = 0.0;
    Material* material = _geometry->getMaterial(i);

    if (material->isFissionable()){
      for (int g=0; g < _num_energy_groups; g++){
        fission_rate += material->getNuSigmaFByGroup(g, CURRENT, temps[i]) * 
          getFluxByValue(i, g, CURRENT);
      }
      
      for (int d=0; d < _num_delayed_groups; d++){
        double delayed_fraction = material->getDelayedFractionByGroup(d, time_from);
        double decay_constant = material->getDecayConstantByGroup(d, time_from);
        double conc = delayed_fraction / decay_constant / _k_eff_0 * fission_rate;
        material->setPrecursorConcByGroup(conc, d, time);
      }
    }
  }
}


void Solver::computePower(int time){

  double* temps = getTemperature(time);
  double cell_power = 0.0;

  for (int i=0; i < _num_shape_cells; i++){

    double fission_rate = 0.0;
    Material* material = _geometry->getMaterial(i);

    if (material->isFissionable()){
      fuel_volume += _volumes[i];
      for (int g=0; g < _num_energy_groups; g++){
        cell_power += material->getSigmaFByGroup(g, CURRENT, temps[i]) * 
          getFluxByValue(i, g, CURRENT) * material->getEnergyPerFission() * _volumes[i];
      }
    }

    setPowerByValue(cell_power, i, time);
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


double Solver::computeAveragePower(int time){
  
  computerPower(time);
  double average_power = 0.0;
  double fuel_volume = 0.0;

  for (int i=0; i < _num_shape_cells; i++){

    Material* material = _geometry->getMaterial(i);
    average_power += getPowerByValue(i, time);

    if (material->isFissionable())
      fuel_volume += _volumes[i];
  }

  /* Divide cumulative power by fuel volume */
  average_power /= fuel_volume;

  return average_power;
}


void Solver::integrateTemperature(int time_from, int time_to){

  double* temps_from = getTemperature(time_from);
  double* temps_to = getTemperature(time_to);
  double dt = _clock->getTime(time_to) - _clock->getTime(time_from);

  #pragma omp parallel for
  for (int i=0; i < _num_shape_cells; i++){

    double fission_rate_from = 0.0;
    double fission_rate_to = 0.0;
    Material* material = _geometry->getMaterial(i);
    
    if (material->isFissionable()){
      for (int g=0; g < _num_energy_groups; g++){
        fission_rate_from += material->getSigmaFByGroup(g, time_from, temps_from[i]) *
          getFluxByValue(i, g, time_from);
        fission_rate_to += material->getSigmaFByGroup(g, time_to, temps_to[i]) *
          getFluxByValue(i, g, time_to);
      }
      
      
      temps_to[i] = temps_from[i] + dt * 0.5 * (fission_rate_from + fission_rate_to) * 
        material->getTemperatureConversionFactor();
    }
  }
}


void Solver::computeShape(int time, int time_flux, int time_amp){

  std::vector< std::vector<int> > amp_to_shape = _geometry->getAmpToShapeMap();
  double* shape = getShape(time);

  #pragma omp parallel for
  for (int i=0; i < _num_amp_cells; i++){

    std::vector<int>::iterator iter;
    for (iter = amp_to_shape.begin(); iter != amp_to_shape.end(); ++iter){
      
      for (int g=0; g < _num_energy_groups; g++){
        double flux = getFluxByValue(*iter, g, time_flux);
        double amp = getAmplitudeByValue(i, g, time_amp);
        shape[(*iter)*_num_energy_groups + g] = flux / amp;
      }
    }
  }
}


void Solver::interpolateShape(int time, int time_forward, int time_backward){

  double dt = _clock->getTime(time_forward) - _clock->getTime(time_backward);
  double wt_bwd = (_clock->getTime(time_forward) - _clock->getTime(time)) / dt;
  double wt_fwd = (_clock->getTime(time) - _clock->getTime(time_backward)) / dt;

  #pragma omp parallel for
  for (int i=0; i < _num_shape_cells; i++){
    for (int g=0; g < _num_energy_groups; g++){
      double shape_fwd = getShapeByValue(i, g, time_forward);
      double shape_bwd = getShapeByValue(i, g, time_backward);
      _shape[time][i*_num_energy_groups + g] = shape_bwd * wt_bwd + shape_fwd * wt_fwd;
    }
  }  
}


void Solver::interpolateDifNonlinear(int time, int time_forward, int time_backward){

  double dt = _clock->getTime(time_forward) - _clock->getTime(time_backward);
  double wt_bwd = (_clock->getTime(time_forward) - _clock->getTime(time)) / dt;
  double wt_fwd = (_clock->getTime(time) - _clock->getTime(time_backward)) / dt;

  #pragma omp parallel for
  for (int i=0; i < _num_shape_cells; i++){
    for (int s=0; s < 6; s++){
      for (int g=0; g < _num_energy_groups; g++){
        double dif_nonlinear_fwd = getDifNonlinearByValue(i, g, s, time_forward);
        double dif_nonlinear_bwd = getDifNonlinearByValue(i, g, s, time_backward);
        _dif_nonlinear[time][(i*6+s)*_num_energy_groups + g] = dif_nonlinear_bwd * wt_bwd 
          + dif_nonlinear_fwd * wt_fwd;
      }
    }  
  }
}


void Solver::reconstructFlux(int time, int time_shape, int time_amp){

  std::vector< std::vector<int> > amp_to_shape = _geometry->getAmpToShapeMap();
  double* flux = getFlux(time);

  #pragma omp parallel for
  for (int i=0; i < _num_amp_cells; i++){

    std::vector<int>::iterator iter;
    for (iter = amp_to_shape.begin(); iter != amp_to_shape.end(); ++iter){
      
      for (int g=0; g < _num_energy_groups; g++){
        double shape = getShapeByValue(*iter, g, time_shape);
        double amp = getAmplitudeByValue(i, g, time_amp);
        flux[(*iter)*_num_energy_groups + g] = shape * amp;
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


void Solver::computeDiffusionCoefficients(int time){

  int nx = _geometry->getNumXAmp();
  int ny = _geometry->getNumYAmp();
  int nz = _geometry->getNumZAmp();
  int ng = _num_energy_groups
  double width = _geometry->getWidth() / nx;
  double height = _geometry->getHeight() / ny;
  double depth = _geometry->getDepth() / nz;
  double* temps = getTemperature(time);
  std::vector< std::vector<int> > amp_to_shape = _geometry->getAmpToShapeMap();

  /* Reconstruct the fine mesh flux at t=time */
  reconstructFlux(time, time, time);

  /* Allocate array of diffusion coefficients */
  double* dif_coefs = new double[_num_amp_cells*ng];

  double trans_tally, rxn_tally, temp, volume, dif_coef, flux;
  Material* material;
  int cell_amp, cell_shape;
  double temp, volume;

  /* Condense the diffusion coefficients */
  for (int i=0; i < _num_amp_cells; i++){
    for (int e=0; e < ng; e++){
      
      trans_tally = 0.0;
      rxn_tally = 0.0;
      
      /* Loop over shape cells in this amp cell */
      std::vector<int>::iterator iter;
      for (iter = amp_to_shape.begin(); iter != amp_to_shape.end(); ++iter){
        
        cell_shape = *iter;
        material = _geometry->getMaterial(cell_shape);
        temp = temps[cell_shape];
        volume = _geometry->getVolume(cell_shape);
        dif_coef = material->getDifCoef(cell_shape, e, time);
        flux = getFluxByValue(cell_shape, e, time);
        
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
          
          if (s == 0 || s == 1 || s == 2)
            sense = -1;
          else
            sense = 1;
          
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
          
          for (int g=0; g < ng; g++){
            
            dif_coef = dif_coefs[cell*ng+g];
            current = getCurrentByValue(cell, g, s, time);
            flux = getAmplitudeByValue(cell, g, time);
            
            if (cell_next == -1){
              dif_linear = 2 * dif_coef / length_perpen / 
                (1 + 4 * dif_coef / length_perpen);
              dif_nonlinear = (sense * dif_linear * flux - current / length) / flux;
              dif_linear *= _boundaries[s];
              dif_nonlinear *= _boundaries[s];
            }
            else{
              flux_next = getAmplitudeByValue(cell_next, g, time);
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
            
            setDifLinearByValue(dif_linear, cell, g, s, time);
            setDifNonlinearByValue(dif_nonlinear, cell, g, s, time);          
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
  int ng = ng;
  double width = _geometry->getWidth() / nx;
  double height = _geometry->getHeight() / ny;
  double depth = _geometry->getDepth() / nz;
  double dt = _clock->getDtInner();

  double* temps = getTemperature(CURRENT);
  double* shape = getShape(CURRENT);
  double* flux = getFlux(CURRENT);
  double* amp = getAmplitude(CURRENT);
  double* frequency = getFrequency(CURRENT);
  double* shape_prev = getShape(PREVIOUS_IN);

  std::vector< std::vector<int> > amp_to_shape = _geometry->getAmpToShapeMap();
 
  reconstructFlux(CURRENT, CURRENT, CURRENT);

  #pragma omp parallel for
  for (int i = 0; i < nx*ny*nz; i++){
    for (int g=0; g < ng; g++)
      frequency[i*ng+g] = 0.0;
  }
  
  for (int z=0; z < nz; z++){
    for (int y=0; y < ny; y++){
      for (int x=0; x < nx; x++){
        int cell_amp = z*nx*ny + y*nx+x;
            
        for (int e=0; e < ng; e++){
          
          int row_amp = cell_amp*ng + e;

          /* Loop over shape cells in this amp cell */
          std::vector<int>::iterator iter;
          for (iter = amp_to_shape.begin(); iter != amp_to_shape.end(); ++iter){
            
            int cell_shape = *iter;
            Material* material = _geometry->getMaterial(cell_shape);
            double temp = temps[cell_shape];
            double volume = _geometry->getVolume(cell_shape);
            int row_shape = cell_shape*ng + e;

            double mult = - material->getVelocity(e, CURRENT, temp) / 
              (volume * shape[row_shape]);

            /* Time absorption term on the diagonal */
            frequency[row_amp] += mult / dt / 
              material->getVelocityByGroup(e, CURRENT, temp) * volume * flux[row_shape];

            frequency[row_amp] -= mult / dt / 
              material->getVelocityByGroup(e, CURRENT, temp) * volume * 
              shape_prev[row_shape] * amp[row_amp];

            /* Delayed neutron precursors */
            for (int d=0; d < _num_delayed_groups; d++){
              frequency[row_amp] -= mult * material->getChiByGroup(e, CURRENT, temp) * 
                material->getDecayConstantByGroup(d) * 
                material->getPrecursorConcByGroup(d, CURRENT) * volume;
            }
          
            /* Absorption term on the diagonal */
            frequency[row_amp] += mult * material->getSigmaAByGroup(e, CURRENT, temp) * 
              volume * flux[row_shape];
            
            /* Buckling term on the diagonal */
            frequency[row_amp] += mult * material->getDifCoefByGroup(e, CURRENT, temp) * 
              volume * _buckling * flux[row_shape];
            
            /* Outscattering term on diagonal */
            for (int g=0; g < ng; g++){
              if (e != g){
                frequency[row_amp] += mult * 
                  material->getSigmaSByGroup(e, g, CURRENT, temp) * volume * 
                  flux[row_shape];
              }
            }
            
            /* Fission terms */
            for (int g=0; g < ng; g++){
              frequency[row_amp] -= mult * 
                (1-material->getDelayedFractionTotal(CURRENT)) * 
                material->getChiByGroup(e, CURRENT, temp) * 
                material->getNuSigmaFByGroup(g, CURRENT, temp) / _k_eff_0 * volume * 
                flux[cell_shape*ng+g];
            }
            
            /* Inscattering term on off diagonals */
            for (int g=0; g < ng; g++){
              if (e != g){
                frequency[row_amp] -= mult * 
                  material->getSigmaSByGroup(g, e, CURRENT, temp) * volume * 
                  flux[cell_shape*ng+g];
              }
            }
          }          
           
          /* LEFT SURFACE */
          
          /* Transport term on diagonal */
          frequency[row_amp] += mult * 
            (getDifLinearByValue(cell_amp, e, 0, CURRENT) +
             getDifNonlinearByValue(cell_amp, e, 0, CURRENT)) * height * depth * 
            amp[row_amp];

          /* Transport term on off diagonals */
          if (x != 0){
            frequency[row_amp] -= mult * 
              (_amp_mesh->getDifLinearByValue(cell_amp, e, 0, CURRENT) -
               _amp_mesh->getDifNonlinearByValue(cell_amp, e, 0, CURRENT)) * height * 
              depth * amp[(cell_amp-1)*ng+e];
          }

          /* RIGHT SURFACE */
          
          /* Transport term on diagonal */
          frequency[row_amp] += mult * 
            (getDifLinearByValue(cell_amp, e, 3, CURRENT) -
             getDifNonlinearByValue(cell_amp, e, 3, CURRENT)) * height * depth * 
            amp[row_amp];
          
          /* Transport term on off diagonals */
          if (x != nx - 1){
            frequency[row_amp] -= mult * 
              (getDifLinearByValue(cell_amp, e, 3, CURRENT) +
               getDifNonlinearByValue(cell_amp, e, 3, CURRENT)) * height * depth * 
              amp[(cell_amp+1)*ng+e];
          }
          
          /* BACK SURFACE */
          
          /* Transport term on diagonal */
          frequency[row_amp] += mult * 
            (getDifLinearByValue(cell_amp, e, 1, CURRENT) +
             getDifNonlinearByValue(cell_amp, e, 1, CURRENT)) * width * depth * 
            amp[row_amp];          
          
          /* Transport term on off diagonals */
          if (y != 0){
            frequency[row_amp] -= mult * 
              (getDifLinearByValue(cell_amp, e, 1, CURRENT) - 
               getDifNonlinearByValue(cell_amp, e, 1, CURRENT)) * width * depth * 
              amp[(cell_amp-nx)*ng+e];
          }
          
          /* FRONT SURFACE */
          
          /* Transport term on diagonal */
          frequency[row_amp] += mult * 
            (getDifLinearByValue(cell_amp, e, 4, CURRENT) -
             getDifNonlinearByValue(cell_amp, e, 4, CURRENT)) * width * depth * 
            amp[row_amp];
          
          /* Transport term on off diagonals */
          if (y != ny - 1){
            frequency[row_amp] -= mult * 
              (getDifLinearByValue(cell_amp, e, 4, CURRENT) +
               getDifNonlinearByValue(cell_amp, e, 4, CURRENT)) * width * depth * 
              amp[(cell_amp+nx)*ng+e];
          }

          /* BOTTOM SURFACE */
          
          /* Transport term on diagonal */
          frequency[row_amp] += mult * 
            (getDifLinearByValue(cell_amp, e, 2, CURRENT) +
             getDifNonlinearByValue(cell_amp, e, 2, CURRENT)) * width * height * 
            amp[row_amp];
          
          /* Transport term on off diagonals */
          if (z != 0){
            frequency[row_amp] -= mult * 
              (getDifLinearByValue(cell_amp, e, 2, CURRENT) -
               getDifNonlinearByValue(cell_amp, e, 2, CURRENT)) * width * height * 
              amp[(cell_amp-nx)*ng+e];
          }
          
          /* TOP SURFACE */
          
          /* Transport term on diagonal */
          frequency[row_amp] += mult * 
            (getDifLinearByValue(cell_amp, e, 5, CURRENT) -
             getDifNonlinearByValue(cell_amp, e, 5, CURRENT)) * width * height * 
            amp[row_amp];
          
          /* Transport term on off diagonals */
          if (z != nz - 1){
            frequency[row_amp] -= mult * 
              (getDifLinearByValue(cell_amp, e, 5, CURRENT) +
               getDifNonlinearByValue(cell_amp, e, 5, CURRENT)) * width * height * 
              amp[(cell_amp+nx)*ng+e];
          }
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


double Solver::setInitialPower(double power){
  _initial_power = power;
}
