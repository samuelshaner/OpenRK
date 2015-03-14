#include "AmpMesh.h"


/**
 * @brief Constructor sets the ID and unique ID for the Material.
 * @param id the user-defined ID for the material
 */
AmpMesh::AmpMesh(double width, double height, double depth, 
                 int num_x, int num_y, int num_z) :
  Mesh(width, height, depth){

  setNumX(num_x);
  setNumY(num_y);
  setNumZ(num_z);

  _shape_mesh = NULL;
  _optically_thick = false;
  _group_indices = NULL;
  _mesh_type = AMPMESH;

  for (int c=0; c < 8; c++){
    _current[c] = NULL;
    _dif_linear[c] = NULL;
    _dif_nonlinear[c] = NULL;
  }

  
  return;
}


/**
 * @brief Destructor deletes all cross-section data structures from memory.
 */
AmpMesh::~AmpMesh() {

  for (int c=0; c < 8; c++){
    if (_current[c] != NULL)
      delete [] _current[c];

    if (_dif_linear[c] != NULL)
      delete [] _dif_linear[c];

    if (_dif_nonlinear[c] != NULL)
      delete [] _dif_nonlinear[c];

  }

  _current.clear();
  _dif_linear.clear();
  _dif_nonlinear.clear();

  std::vector< std::vector<int> >::iterator iter;
  for (iter=_shape_map.begin(); iter != _shape_map.end(); ++iter)
    *iter.clear();

  _shape_map.clear();

  if (_group_indices != NULL)
    delete [] _group_indices;
}


void AmpMesh::setOpticallyThick(bool optically_thick){
  _optically_thick = optically_thick;
}


void AmpMesh::setShapeMesh(ShapeMesh* mesh){
  _shape_mesh = mesh;

  for (int i=0; i < getNumCells(); i++){
    std::vector<int> *shape_cells = new std::vector<int>;
    _shape_map.push_back(*shape_cells);
  }  

  if (_shape_mesh->getMeshType() == STRUCTURED_SHAPE_MESH){
    int nx = static_cast<StructuredShapeMesh*>(_shape_mesh)->getNumX();
    int ny = static_cast<StructuredShapeMesh*>(_shape_mesh)->getNumY();
    int nz = static_cast<StructuredShapeMesh*>(_shape_mesh)->getNumZ();
    int num_refines_x = nx / _num_x;
    int num_refines_y = ny / _num_y;
    int num_refines_z = nz / _num_z;
    
    for (int z=0; z < nz; z++){
      int zz = z / num_refines_z;
      for (int y=0; y < ny; y++){
        int yy = y / num_refines_y;
        for (int x=0; x < nx; x++){
          int xx = x / num_refines_x;
          int amp_cell = zz*_num_x*_num_y + yy*_num_x + xx;
          int shape_cell = z*nx*ny + y*nx + x;
          _shape_map[amp_cell].push_back(shape_cell);
          _shape_mesh->setAmpCellContainingShapeCell(shape_cell, amp_cell);
        }
      }
    }
  }
}


void AmpMesh::setShapeMap(int** regions_to_cells, int* regions_per_cell, int num_cells){
  
  if (num_cells != getNumCells())
    log_printf(ERROR, "Unable to set the shape map with a regions to cells "
               "array of length %i since there are %i cells in the AmpMesh"
               , num_cells, getNumCells());

  if (_shape_mesh == NULL)
    log_printf(ERROR, "Unable to set the shape map since the shape mesh has not"
               " been set");
  
  for (int i=0; i < getNumCells(); i++){

    if (regions_per_cell[i] <= 0)
      log_printf(ERROR, "Unable to set a non-positive number of regions %i for"
                 " cell %i", regions_per_cell[i], i);

    for (int r=0; r < regions_per_cell[i]; r++){
      _shape_map[i].push_back(regions_to_cells[i][r]);
      _shape_mesh->setAmpCellContainingShapeCell(regions_to_cells[i][r], i);
    }
  }
}


void AmpMesh::setFluxByValue(double flux, int cell, int group, int position){

  if (cell >= getNumCells() || cell < 0)
    log_printf(ERROR, "Unable to set flux by value for Amp Mesh "
               "cell %i since there are only %i cells", cell, getNumCells());

  if (group >= _num_amp_energy_groups || group < 0)
    log_printf(ERROR, "Unable to set flux by value for Amp Mesh group "
               "%i since there are only %i groups", group, _num_amp_energy_groups);

  if (position >= 8  || position < 0)
    log_printf(ERROR, "Unable to set flux by value for Amp Mesh clock "
               "position %i since there are only 8 clock positions", position);

  if (flux < 0.0)
    log_printf(NORMAL, "Unable to set flux by value for Amp Mesh "
               "for cell %i and group %i to a negative value %f", cell, group, flux);

  _flux[position][cell*_num_amp_energy_groups + group] = flux;
}


void AmpMesh::setCurrentByValue(double current, int cell, int group, int side, int position){

  if (cell >= getNumCells() || cell < 0)
    log_printf(ERROR, "Unable to set current by value for Amp Mesh "
               "cell %i since there are only %i cells", cell, getNumCells());

  if (group >= _num_amp_energy_groups || group < 0)
    log_printf(ERROR, "Unable to set current by value for Amp Mesh group "
               "%i since there are only %i groups", group, _num_amp_energy_groups);
  
  if (position >= 8 || position < 0)
    log_printf(ERROR, "Unable to set current by value for Amp Mesh clock "
               "position %i since there are only 8 clock positions", position);

  if (side < 0 || side >= 6)
    log_printf(ERROR, "Unable to set current by value for Amp Mesh "
               "for side %i since there are only 6 sides", side);
  
  _current[position][(cell*6 + side)*_num_amp_energy_groups + group] = current;
}


void AmpMesh::setDifLinearByValue(double dif_linear, int cell, int group, int side, int position){

  if (cell >= getNumCells() || cell < 0)
    log_printf(ERROR, "Unable to set dif linear by value for Amp Mesh "
               "cell %i since there are only %i cells", cell, getNumCells());

  if (group >= _num_amp_energy_groups || group < 0)
    log_printf(ERROR, "Unable to set dif linear by value for Amp Mesh group "
               "%i since there are only %i groups", group, _num_amp_energy_groups);

  if (position >= 8 || position < 0)
    log_printf(ERROR, "Unable to set dif linear by value for Amp Mesh clock "
               "position %i since there are only 8 clock positions", position);

  if (side < 0 || side >= 6)
    log_printf(ERROR, "Unable to set dif linear by value for Amp Mesh "
               "for side %i since there are only 6 sides", side);

  _dif_linear[position][(cell*6 + side)*_num_amp_energy_groups + group] = dif_linear;
}


void AmpMesh::setDifNonlinearByValue(double dif_nonlinear, int cell, int group, int side, int position){

  if (cell >= getNumCells() || cell < 0)
    log_printf(ERROR, "Unable to set dif nonlinear by value for Amp Mesh "
               "cell %i since there are only %i cells", cell, getNumCells());

  if (group >= _num_amp_energy_groups || group < 0)
    log_printf(ERROR, "Unable to set dif nonlinear by value for Amp Mesh group "
               "%i since there are only %i groups", group, _num_amp_energy_groups);

  if (position >= 8 || position < 0)
    log_printf(ERROR, "Unable to set dif nonlinear by value for Amp Mesh clock "
               "position %i since there are only 8 clock positions", position);

  if (side < 0 || side >= 6)
    log_printf(ERROR, "Unable to set dif nonlinear by value for Amp Mesh "
               "for side %i since there are only 6 sides", side);

  _dif_nonlinear[position][(cell*6 + side)*_num_amp_energy_groups + group] 
    = dif_nonlinear;
}


double AmpMesh::getFluxByValue(int cell, int group, int position){

  if (cell >= getNumCells() || cell < 0)
    log_printf(ERROR, "Unable to get flux by value for Amp Mesh "
               "cell %i since there are only %i cells", cell, getNumCells());

  if (group >= _num_amp_energy_groups || group < 0)
    log_printf(ERROR, "Unable to get flux by value for Amp Mesh group "
               "%i since there are only %i groups", group, _num_amp_energy_groups);

  if (position >= 8 || position < 0)
    log_printf(ERROR, "Unable to get flux by value for Amp Mesh clock "
               "position %i since there are only 8 clock positions", position);

  return _flux[position][cell*_num_amp_energy_groups + group];
}


double AmpMesh::getCurrentByValue(int cell, int group, int side, int position){

  if (cell >= getNumCells() || cell < 0)
    log_printf(ERROR, "Unable to get current by value for Amp Mesh "
               "cell %i since there are only %i cells", cell, getNumCells());

  if (group >= _num_amp_energy_groups || group < 0)
    log_printf(ERROR, "Unable to get current by value for Amp Mesh group "
               "%i since there are only %i groups", group, _num_amp_energy_groups);
  
  if (position >= 8 || position < 0)
    log_printf(ERROR, "Unable to get current by value for Amp Mesh clock "
               "position %i since there are only 8 clock positions", position);

  if (side < 0 || side >= 6)
    log_printf(ERROR, "Unable to get current by value for Amp Mesh "
               "for side %i since there are only 6 sides", side);
  
  return _current[position][(cell*6 + side)*_num_amp_energy_groups + group];
}


double AmpMesh::getDifLinearByValue(int cell, int group, int side, int position){

  if (cell >= getNumCells() || cell < 0)
    log_printf(ERROR, "Unable to get dif linear by value for Amp Mesh "
               "cell %i since there are only %i cells", cell, getNumCells());

  if (group >= _num_amp_energy_groups || group < 0)
    log_printf(ERROR, "Unable to get dif linear by value for Amp Mesh group "
               "%i since there are only %i groups", group, _num_amp_energy_groups);

  if (position >= 8 || position < 0)
    log_printf(ERROR, "Unable to get dif linear by value for Amp Mesh clock "
               "position %i since there are only 8 clock positions", position);

  if (side < 0 || side >= 6)
    log_printf(ERROR, "Unable to get dif linear by value for Amp Mesh "
               "for side %i since there are only 6 sides", side);

  return _dif_linear[position][(cell*6 + side)*_num_amp_energy_groups + group];
}


double AmpMesh::getDifNonlinearByValue(int cell, int group, int side, int position){

  if (cell >= getNumCells() || cell < 0)
    log_printf(ERROR, "Unable to get dif nonlinear by value for Amp Mesh "
               "cell %i since there are only %i cells", cell, getNumCells());

  if (group >= _num_amp_energy_groups || group < 0)
    log_printf(ERROR, "Unable to get dif nonlinear by value for Amp Mesh group "
               "%i since there are only %i groups", group, _num_amp_energy_groups);

  if (position >= 8 || position < 0)
    log_printf(ERROR, "Unable to get dif nonlinear by value for Amp Mesh clock "
               "position %i since there are only 8 clock positions", position);

  if (side < 0 || side >= 6)
    log_printf(ERROR, "Unable to get dif nonlinear by value for Amp Mesh "
               "for side %i since there are only 6 sides", side);

  return _dif_nonlinear[position][(cell*6 + side)*_num_amp_energy_groups + group];
}


void AmpMesh::copyFlux(int position_from, int position_to){

  if (position_from >= 8 || position_from < 0)
    log_printf(ERROR, "Unable to copy flux from position %i since "
               "there are only 8 clock positions", position_from);

  if (position_to >= 8 || position_to < 0)
    log_printf(ERROR, "Unable to copy flux to position %i since "
               "there are only 8 clock positions", position_to);
    
  if (_flux[position_from] == NULL)
    log_printf(ERROR, "Unable to copy flux from position %i since "
               "the fluxs have not been initialized yet", position_from);

  if (_flux[position_to] == NULL)
    log_printf(ERROR, "Unable to copy flux to position %i since "
               "the fluxs have not been initialized yet", position_to);

  std::copy(_flux[position_from], 
            _flux[position_from] + getNumCells()*_num_amp_energy_groups, 
            _flux[position_to]);
}


void AmpMesh::copyCurrent(int position_from, int position_to){

  if (position_from >= 8 || position_from < 0)
    log_printf(ERROR, "Unable to copy current from position %i since "
               "there are only 8 clock positions", position_from);

  if (position_to >= 8 || position_to < 0)
    log_printf(ERROR, "Unable to copy current to position %i since "
               "there are only 8 clock positions", position_to);
    
  if (_current[position_from] == NULL)
    log_printf(ERROR, "Unable to copy current from position %i since "
               "the currents have not been initialized yet", position_from);

  if (_current[position_to] == NULL)
    log_printf(ERROR, "Unable to copy current to position %i since "
               "the currents have not been initialized yet", position_to);

  std::copy(_current[position_from], 
            _current[position_from] + getNumCells()*_num_amp_energy_groups*6, 
            _current[position_to]);
}


void AmpMesh::copyDifLinear(int position_from, int position_to){

  if (position_from >= 8 || position_from < 0)
    log_printf(ERROR, "Unable to copy dif linear from position %i since "
               "there are only 8 clock positions", position_from);

  if (position_to >= 8 || position_to < 0)
    log_printf(ERROR, "Unable to copy dif linear to position %i since "
               "there are only 8 clock positions", position_to);
    
  if (_dif_linear[position_from] == NULL)
    log_printf(ERROR, "Unable to copy dif linear from position %i since "
               "the dif linears have not been initialized yet", position_from);

  if (_dif_linear[position_to] == NULL)
    log_printf(ERROR, "Unable to copy dif linear to position %i since "
               "the dif linears have not been initialized yet", position_to);

  std::copy(_dif_linear[position_from], 
            _dif_linear[position_from] + getNumCells()*_num_amp_energy_groups*6, 
            _dif_linear[position_to]);
}


void AmpMesh::copyDifNonlinear(int position_from, int position_to){

  if (position_from >= 8 || position_from < 0)
    log_printf(ERROR, "Unable to copy dif nonlinear from position %i since "
               "there are only 8 clock positions", position_from);

  if (position_to >= 8 || position_to < 0)
    log_printf(ERROR, "Unable to copy dif nonlinear to position %i since "
               "there are only 8 clock positions", position_to);
    
  if (_dif_linear[position_from] == NULL)
    log_printf(ERROR, "Unable to copy dif nonlinear from position %i since "
               "the dif nonlinears have not been initialized yet", position_from);

  if (_dif_linear[position_to] == NULL)
    log_printf(ERROR, "Unable to copy dif nonlinear to position %i since "
               "the dif nonlinears have not been initialized yet", position_to);

  std::copy(_dif_nonlinear[position_from], 
            _dif_nonlinear[position_from] + getNumCells()*_num_amp_energy_groups*6, 
            _dif_nonlinear[position_to]);
}


void AmpMesh::condenseMaterials(int position, bool save_flux){

  if (_shape_mesh == NULL)
    log_printf(ERROR, "Unable to condense the shape mesh materials since "
               "the shape mesh has not been set");

  double shape_cell_volume;
  double amp_cell_volume = getCellVolume();
  double* temps = _shape_mesh->getTemperature(position);
  double chi[_num_amp_energy_groups];
  double sigma_s[_num_amp_energy_groups];
  Material* amp_mat;
  Material* shape_mat;

  #pragma omp parallel for private(chi, sigma_s, amp_mat, shape_mat, shape_cell_volume)
  for (int i=0; i < getNumCells(); i++){
    
    amp_mat = getMaterial(i);
    std::vector<int>::iterator iter;
    
    for (int g=0; g < _num_amp_energy_groups; g++){

      double sigma_a = 0.0;
      double sigma_t = 0.0;
      double sigma_f = 0.0;
      double nu_sigma_f = 0.0;
      double dif_coef = 0.0;
      double rxn = 0.0;
      double production = 0.0;
      double velocity = 0.0;
      double volume;
      double temp;
      
      for (int h=0; h < _num_amp_energy_groups; h++){
        sigma_s[h] = 0.0;
        chi[h] = 0.0;
      }
      
      for (iter = _shape_map[i].begin(); iter != _shape_map[i].end(); ++iter){

        shape_mat = _shape_mesh->getMaterial(*iter);
        shape_cell_volume = _shape_mesh->getCellVolume(*iter);
        temp = temps[*iter];

        for (int e=0; e < _num_amp_energy_groups; e++){
          double chi_tally = 0.0;

          for (int ee=_group_indices[e]; ee < _group_indices[e+1]; ee++)
            chi_tally += shape_mat->getChiByGroup(ee, position, temp);
          
          for (int h=0; h < _num_shape_energy_groups; h++){
            chi[e] += chi_tally * shape_mat->getNuSigmaFByGroup(h, position, temp) *
              _shape_mesh->getFluxByValue(*iter, h, position) * shape_cell_volume;
            production += chi_tally * shape_mat->getNuSigmaFByGroup(h, position, temp) *
              _shape_mesh->getFluxByValue(*iter, h, position) * shape_cell_volume;
          }
        }
      }

      for (int gg=_group_indices[g]; gg < _group_indices[g+1]; gg++){

        double rxn_group = 0.0;
        volume = 0.0;
        double sigma_t_group = 0.0;

      for (iter = _shape_map[i].begin(); iter != _shape_map[i].end(); ++iter){

          shape_mat = _shape_mesh->getMaterial(*iter);
          shape_cell_volume = _shape_mesh->getCellVolume(*iter);
          temp = temps[*iter];

          double shape = _shape_mesh->getFluxByValue(*iter, gg, position);
          sigma_a += shape_mat->getSigmaAByGroup(gg, position, temp) * shape * shape_cell_volume;
          sigma_t += shape_mat->getSigmaTByGroup(gg, position, temp) * shape * shape_cell_volume;
          sigma_f += shape_mat->getSigmaFByGroup(gg, position, temp) * shape * shape_cell_volume;
          nu_sigma_f += shape_mat->getNuSigmaFByGroup(gg, position, temp) * shape * shape_cell_volume;
          rxn += shape * shape_cell_volume;
          rxn_group += shape * shape_cell_volume;
          volume += shape_cell_volume;
          sigma_t_group += shape_mat->getSigmaTByGroup(gg, position, temp) 
            * shape * shape_cell_volume;
          velocity += 1.0 / shape_mat->getVelocityByGroup(gg, position, temp) 
            * shape * shape_cell_volume;

          for (int h=0; h < _num_shape_energy_groups; h++){
            sigma_s[_shape_mesh->getAmpGroup(h)] += 
              shape_mat->getSigmaSByGroup(gg, h, position, temp) *
              shape * shape_cell_volume;
          }
        }

        dif_coef += rxn_group / (3.0 * sigma_t_group / rxn_group);
      }

      amp_mat->setSigmaAByGroup(sigma_a / rxn, g, position);
      amp_mat->setSigmaTByGroup(sigma_t / rxn, g, position);
      amp_mat->setSigmaFByGroup(sigma_f / rxn, g, position);
      amp_mat->setNuSigmaFByGroup(nu_sigma_f / rxn, g, position);
      amp_mat->setDifCoefByGroup(dif_coef / rxn, g, position);
      amp_mat->setVelocityByGroup(1.0 / (velocity / rxn), g, position);

      if (save_flux)
        _flux[position][i*_num_amp_energy_groups + g] = rxn / volume;

      if (production != 0.0)
        amp_mat->setChiByGroup(chi[g] / production, g, position);
      else
        amp_mat->setChiByGroup(0.0, g, position);      

      for (int e=0; e < _num_amp_energy_groups; e++)
        amp_mat->setSigmaSByGroup(sigma_s[e] / rxn, g, e, position);
    }

    for (int d=0; d < _num_delayed_groups; d++){
      double precursor_conc = 0.0;

      for (iter = _shape_map[i].begin(); iter != _shape_map[i].end(); ++iter){
        if (_shape_mesh->getMaterial(*iter)->isFissionable()){
          precursor_conc += _shape_mesh->getMaterial(*iter)->getPrecursorConcByGroup(d, position) 
            * shape_cell_volume;
        }
      }
        
      amp_mat->setPrecursorConcByGroup(precursor_conc / amp_cell_volume, d, position);
    }
  }
}


void AmpMesh::initialize(){
  
  Mesh::initialize();

  for (int c=0; c < 8; c++){
    if (_flux[c] != NULL)
      delete [] _flux[c];

    if (_current[c] != NULL)
      delete [] _current[c];

    if (_dif_linear[c] != NULL)
      delete [] _dif_linear[c];

    if (_dif_nonlinear[c] != NULL)
      delete [] _dif_nonlinear[c];
  }

  for (int i=0; i < getNumCells(); i++){
    Material* material = new Material();
    material->setNumEnergyGroups(_num_amp_energy_groups);
    material->setNumDelayedGroups(_num_delayed_groups);
    _materials[i] = material;
  }

  for(int c=0; c < 8; c++){
    _dif_linear[c] = new double[getNumCells() * _num_amp_energy_groups * 6];
    _dif_nonlinear[c] = new double[getNumCells() * _num_amp_energy_groups * 6];
    _current[c] = new double[getNumCells() * _num_amp_energy_groups * 6];
    _flux[c] = new double[getNumCells() * _num_amp_energy_groups];

    memset(_dif_linear[c], 0.0, sizeof(double) * getNumCells() * _num_amp_energy_groups * 6);
    memset(_dif_nonlinear[c], 0.0, sizeof(double) * getNumCells() * _num_amp_energy_groups * 6);
    memset(_current[c], 0.0, sizeof(double) * getNumCells() * _num_amp_energy_groups * 6);
    memset(_flux[c], 1.0, sizeof(double) * getNumCells() * _num_amp_energy_groups);
  }  
}


void AmpMesh::computePower(int position){

  if (position >= 8 || position < 0)
    log_printf(ERROR, "Unable to compute the power for the AmpMesh for "
               "position %i since there are only 8 clock positions", position);

  if (_temperature[position] == NULL)
    log_printf(ERROR, "Unable to compute the power since the temperature has not"
               " been initialized");
  
  for (int i=0; i < getNumCells(); i++){
    double fission_rate = 0.0;
    double temp = _temperature[position][i];
    
    for (int g=0; g < _num_amp_energy_groups; g++)
      fission_rate += _materials[i]->getSigmaFByGroup(g, position, temp) * getFluxByValue(i, g, position);
    
    _power[position][i] = fission_rate;
  }
}


void AmpMesh::computeCurrent(int position){

  if (position >= 8 || position < 0)
    log_printf(ERROR, "Unable to compute the current for the AmpMesh for "
               "position %i since there are only 8 clock positions", position);

  if (_shape_mesh == NULL)
    log_printf(ERROR, "Unable to compute the current for the amp mesh "
               "since the shape mesh has not been set");

  if (_materials == NULL)
    log_printf(ERROR, "Unable to compute the current for the amp mesh "
               "since the materials have not been initialized");

  int sm_nx = _shape_mesh->getNumX();
  int sm_ny = _shape_mesh->getNumY();
  int sm_nz = _shape_mesh->getNumZ();
  int sm_cw = _shape_mesh->getCellWidth();
  int sm_ch = _shape_mesh->getCellHeight();
  int sm_cd = _shape_mesh->getCellDepth();
  int num_refines_x = sm_nx / _num_x;
  int num_refines_y = sm_ny / _num_y;
  int num_refines_z = sm_nz / _num_z;
  double* temps = _shape_mesh->getTemperature(position);

  #pragma omp parallel for 
  for (int i=0; i < getNumCells(); i++){

    std::vector<int>::iterator iter;
    double flux, dif_linear, d, d_next, temp, temp_next, flux_next;
    int cell_next;
    double length_perpen, length;
    Material *mat, *mat_next;

    for (int g=0; g < _num_amp_energy_groups; g++){
      for (int s=0; s < 6; s++){
        double current = 0.0;
        for (iter = _shape_map[i].begin(); iter != _shape_map[i].end(); ++iter){

          int x = (*iter % (sm_nx*sm_ny)) % sm_nx;
          int y = (*iter % (sm_nx*sm_ny)) / sm_nx;
          int z = (*iter / (sm_nx*sm_ny));
          mat = _shape_mesh->getMaterial(*iter);
          temp = temps[*iter];
          cell_next = _shape_mesh->getNeighborCell(x, y, z, s);
          temp_next = temps[cell_next];
          mat_next = _shape_mesh->getNeighborMaterial(x, y, z, s);

          if (s == 0 && x % num_refines_x == 0){
            for (int gg=_group_indices[g]; gg < _group_indices[g+1]; gg++){
              flux = _shape_mesh->getFluxByValue(*iter, gg, position);
              d = mat->getDifCoefByGroup(gg, position, temp);
              length_perpen = sm_cw;
              length = sm_ch * sm_cd;

              if (mat_next == NULL){
                dif_linear = 2 * d / length_perpen / (1 + 4 * d / length_perpen);
                dif_linear *= _shape_mesh->getBoundary(s);
                current += - dif_linear * flux * length;
              }
              else{
                d_next = mat_next->getDifCoefByGroup(gg, position, temp_next);
                dif_linear = 2 * d * d_next / (length_perpen * d + length_perpen * d_next);
                flux_next = _shape_mesh->getFluxByValue(cell_next, gg, position);
                current += - dif_linear * (flux - flux_next) * length;
              }
            }
          }

          else if (s == 1 && y % num_refines_y == 0){
            for (int gg=_group_indices[g]; gg < _group_indices[g+1]; gg++){
              flux = _shape_mesh->getFluxByValue(*iter, gg, position);
              d = mat->getDifCoefByGroup(gg, position, temp);
              length_perpen = sm_ch;
              length = sm_cw * sm_cd;

              if (mat_next == NULL){
                dif_linear = 2 * d / length_perpen / (1 + 4 * d / length_perpen);
                dif_linear *= _shape_mesh->getBoundary(s);
                current += - dif_linear * flux * length;
              }
              else{
                d_next = mat_next->getDifCoefByGroup(gg, position, temp_next);
                dif_linear = 2 * d * d_next / (length_perpen * d + length_perpen * d_next);
                flux_next = _shape_mesh->getFluxByValue(cell_next, gg, position);
                current += - dif_linear * (flux - flux_next) * length;
              }
            }            
          }

          else if (s == 2 && z % num_refines_z == 0){
            for (int gg=_group_indices[g]; gg < _group_indices[g+1]; gg++){
              flux = _shape_mesh->getFluxByValue(*iter, gg, position);
              d = mat->getDifCoefByGroup(gg, position, temp);
              length_perpen = sm_cd;
              length = sm_cw * sm_ch;

              if (mat_next == NULL){
                dif_linear = 2 * d / length_perpen / (1 + 4 * d / length_perpen);
                dif_linear *= _shape_mesh->getBoundary(s);
                current += - dif_linear * flux * length;
              }
              else{
                d_next = mat_next->getDifCoefByGroup(gg, position, temp_next);
                dif_linear = 2 * d * d_next / (length_perpen * d + length_perpen * d_next);
                flux_next = _shape_mesh->getFluxByValue(cell_next, gg, position);
                current += - dif_linear * (flux - flux_next) * length;
              }
            }            
          }

          else if (s == 3 && x % num_refines_x == num_refines_x - 1){
            for (int gg=_group_indices[g]; gg < _group_indices[g+1]; gg++){
              flux = _shape_mesh->getFluxByValue(*iter, gg, position);
              d = mat->getDifCoefByGroup(gg, position, temp);
              length_perpen = sm_cw;
              length = sm_ch * sm_cd;

              if (mat_next == NULL){
                dif_linear = 2 * d / length_perpen / (1 + 4 * d / length_perpen);
                dif_linear *= _shape_mesh->getBoundary(s);
                current += dif_linear * flux * length;
              }
              else{
                d_next = mat_next->getDifCoefByGroup(gg, position, temp_next);
                dif_linear = 2 * d * d_next / (length_perpen * d + length_perpen * d_next);
                flux_next = _shape_mesh->getFluxByValue(cell_next, gg, position);
                current += - dif_linear * (flux_next - flux) * length;
              }
            }            
          }

          else if (s == 4 && y % num_refines_y == num_refines_y - 1){
            for (int gg=_group_indices[g]; gg < _group_indices[g+1]; gg++){
              flux = _shape_mesh->getFluxByValue(*iter, gg, position);
              d = mat->getDifCoefByGroup(gg, position, temp);
              length_perpen = sm_ch;
              length = sm_cw * sm_cd;

              if (mat_next == NULL){
                dif_linear = 2 * d / length_perpen / (1 + 4 * d / length_perpen);
                dif_linear *= _shape_mesh->getBoundary(s);
                current += dif_linear * flux * length;
              }
              else{
                d_next = mat_next->getDifCoefByGroup(gg, position, temp_next);
                dif_linear = 2 * d * d_next / (length_perpen * d + length_perpen * d_next);
                flux_next = _shape_mesh->getFluxByValue(cell_next, gg, position);
                current += - dif_linear * (flux_next - flux) * length;
              }
            }            
          }
          else if (s == 5 && z % num_refines_z == num_refines_z - 1){
            for (int gg=_group_indices[g]; gg < _group_indices[g+1]; gg++){
              flux = _shape_mesh->getFluxByValue(*iter, gg, position);
              d = mat->getDifCoefByGroup(gg, position, temp);
              length_perpen = sm_cd;
              length = sm_cw * sm_ch;
              
              if (mat_next == NULL){
                dif_linear = 2 * d / length_perpen / (1 + 4 * d / length_perpen);
                dif_linear *= _shape_mesh->getBoundary(s);
                current += dif_linear * flux * length;
              }
              else{
                d_next = mat_next->getDifCoefByGroup(gg, position, temp_next);
                dif_linear = 2 * d * d_next / (length_perpen * d + length_perpen * d_next);
                flux_next = _shape_mesh->getFluxByValue(cell_next, gg, position);
                current += - dif_linear * (flux_next - flux) * length;
              }
            }            
          }
        }
        
        _current[position][(i*6+s) * _num_amp_energy_groups + g] = current;
      }
    }
  }
}


double AmpMesh::computeDifCorrect(double dif_coef, double length){

  double f = 0.0;
  
  if (_optically_thick){
    double mu = cos(asin(0.798184));
    double expon = exp(- length / (3 * dif_coef * mu));
    double alpha = (1 + expon) / (1 - expon) - 2 * 3 * dif_coef * mu / length;
    double rho = mu * alpha;
    f = 1 + length * rho / (2 * dif_coef);
  }
  else{
    f = 1.0;
  }
  
  return f;
}


void AmpMesh::computeDifCoefs(int position){

  if (position >= 8 || position < 0)
    log_printf(ERROR, "Unable to compute the dif coefs for the AmpMesh for "
               "position %i since there are only 8 clock positions", position);

  if (_materials == NULL)
    log_printf(ERROR, "Unable to compute the dif coefs for the amp mesh "
               "since the materials have not been initialized");

  int nx = _num_x;
  int ny = _num_y;
  int nz = _num_z;
  int ng = _num_amp_energy_groups;
  double width = getCellWidth();
  double height = getCellHeight();
  double depth = getCellDepth();
  double* temps = _temperature[position];
  

  for (int z=0; z < nz; z++){
    #pragma omp parallel for  
    for (int y=0; y < ny; y++){
      for (int x=0; x < nx; x++){


        int sense;
        double length, length_perpen;
        double dif_coef, current, flux, f, dif_linear, dif_nonlinear, flux_next, dif_coef_next, f_next;
        int cell = z*nx*ny+y*nx+x;
        double temp = temps[cell];
        
        for (int s=0; s < 6; s++){
          
          int cell_next = getNeighborCell(x, y, z, s);
          
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
            
            dif_coef = _materials[cell]->getDifCoefByGroup(g, position, temp);
            current = getCurrentByValue(cell, g, s, position);
            flux = getFluxByValue(cell, g, position);
            f = computeDifCorrect(dif_coef, length_perpen);
            
            if (cell_next == -1){
              dif_linear = 2 * dif_coef * f / length_perpen / (1 + 4 * dif_coef * f / length_perpen);
              dif_nonlinear = (sense * dif_linear * flux - current / length) / flux;
              dif_linear *= _boundaries[s];
              dif_nonlinear *= _boundaries[s];
            }
            else{
              flux_next = getFluxByValue(cell_next, g, position);
              dif_coef_next = _materials[cell_next]->getDifCoefByGroup(g, position, temps[cell_next]);
              f_next = computeDifCorrect(dif_coef_next, length_perpen);
              dif_linear = 2 * dif_coef * f * dif_coef_next * f_next / (length_perpen * dif_coef * f +
                                                                        length_perpen * dif_coef_next * f_next);
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
            
            setDifLinearByValue(dif_linear, cell, g, s, position);
            setDifNonlinearByValue(dif_nonlinear, cell, g, s, position);          
          }
        }
      }
    }
  }
}


AmpMesh* AmpMesh::clone(){

  if (_num_shape_energy_groups == 0)
    log_printf(ERROR, "Unable to clone the AmpMesh since the number "
               "of shape energy groups has not been set");

  if (_num_amp_energy_groups == 0)
    log_printf(ERROR, "Unable to clone the AmpMesh since the number "
               "of amp energy groups has not been set");

  if (_num_delayed_groups == 0)
    log_printf(ERROR, "Unable to clone the AmpMesh since the number "
               "of delayed groups has not been set");
  
  AmpMesh* mesh = new AmpMesh(getWidth(), getHeight(), getDepth(), _num_x, _num_y, _num_z);

  mesh->setNumShapeEnergyGroups(_num_shape_energy_groups);
  mesh->setNumAmpEnergyGroups(_num_amp_energy_groups);
  mesh->setNumDelayedGroups(_num_delayed_groups);
  mesh->setBuckling(_buckling);

  if (_k_eff_0 != 0.0)
    mesh->setKeff0(_k_eff_0);
 
  for (int s=0; s < 6; s++)
    mesh->setBoundary(s, _boundaries[s]);

  if (_clock != NULL)
    mesh->setClock(_clock);

  if (_decay_constants != NULL)
    mesh->setDecayConstants(_decay_constants, _num_delayed_groups);

  if (_delayed_fractions != NULL)
    mesh->setDelayedFractions(_delayed_fractions, _num_delayed_groups);  
}


void AmpMesh::setGroupStructure(int* group_indices, int length_group_indices){

  if (_group_indices != NULL)
    delete [] _group_indices;

  _group_indices = new int[length_group_indices];

  if (group_indices == NULL){
    for (int i = 0; i < length_group_indices; i++){
      _group_indices[i] = i;
    }
  }
  else{
    if (group_indices[0] != 0)
      log_printf(ERROR, "The first value in group indices must be 1!");

    /* Set first group indice to 0 */
    _group_indices[0] = 0;

    /* Set MOC group bounds for rest of CMFD energy groups */
    for (int i = 1; i < length_group_indices; i++){
      /* Check that the group indices are always increasing */
      if (group_indices[i] <= group_indices[i-1])
        log_printf(ERROR, "The group indices must be increasing!");

      _group_indices[i] = group_indices[i];
    }
  }
}


double AmpMesh::computeAveragePower(int position){

  if (position >= 8 || position < 0)
    log_printf(ERROR, "Unable to compute the average power for the AmpMesh for "
               "position %i since there are only 8 clock positions", position);

  if (_fuel_volume == 0.0)
    computeFuelVolume();   
  
  computePower(position);

  return pairwise_sum(_power[position], getNumCells()) * getCellVolume() / _fuel_volume;
}  


double AmpMesh::computePowerL2Norm(int position_1, int position_2){

  if (position_1 >= 8 || position_1 < 0)
    log_printf(ERROR, "Unable to compute the power l2 norm for the AmpMesh for "
               "position_1 %i since there are only 8 clock positions", position_1);

  if (position_2 >= 8 || position_2 < 0)
    log_printf(ERROR, "Unable to compute the power l2 norm for the AmpMesh for "
               "position_2 %i since there are only 8 clock positions", position_2);

  computePower(position_1);
  computePower(position_2);

  double* power_residual = new double[getNumCells()];
  memset(power_residual, 0.0, sizeof(double) * getNumCells());

  #pragma omp parallel for
  for (int i=0; i < getNumCells(); i++){
    if (_power[position_1][i] > 0.0)
      power_residual[i] = pow((_power[position_1][i] - _power[position_2][i]) 
                              / _power[position_1][i], 2);
  }

  double residual = sqrt(pairwise_sum(power_residual, getNumCells()));
  delete [] power_residual;

  return residual;
}


void AmpMesh::interpolateDifNonlinear(int position_begin, int position_end, int position){

  if (position_begin >= 8 || position_begin < 0)
    log_printf(ERROR, "Unable to interpolate the dif nonlinear for the AmpMesh for "
               "position begin %i since there are only 8 clock positions", position_begin);

  if (position_end >= 8 || position_end < 0)
    log_printf(ERROR, "Unable to interpolate the dif nonlinear for the AmpMesh for "
               "position end %i since there are only 8 clock positions", position_end);

  if (position >= 8 || position < 0)
    log_printf(ERROR, "Unable to interpolate the dif nonlinear for the AmpMesh for "
               "position %i since there are only 8 clock positions", position);

  double dt = _clock->getTime(position_end) - _clock->getTime(position_begin);
  double wt_begin = (_clock->getTime(position_end) - _clock->getTime(position)) / dt;
  double wt_end = (_clock->getTime(position) - _clock->getTime(position_begin)) / dt;

  #pragma omp parallel for
  for (int i=0; i < getNumCells()*_num_amp_energy_groups*6; i++){
    _dif_nonlinear[position][i] = _dif_nonlinear[position_begin][i] * wt_begin
      + _dif_nonlinear[position_end][i] * wt_end;
  }
}


void AmpMesh::setCurrent(double* current, int num_cells_time_groups){

  if (num_cells_times_groups != getNumCells() * _num_shape_energy_groups)
    log_printf(ERROR, "Unable to set current with %i values when the number of "
               "cells is %i and shape energy groups is %i", getNumCells(),
               _num_shape_energy_groups);

  for (int c=0; c < 8; c++){
    if (_current[c] == NULL)
      log_printf(ERROR, "Unable to set current for Amp Mesh since the"
                 " current has not been initialized");
    
    std::copy(current, current + num_cells_times_groups, _current[c]);
  }  
}


int AmpMesh::getNeighborCell(int x, int y, int z, int side){

  if (side < 0 || side > 5)
    log_printf(ERROR, "Unable to get neighbor cell for side %d as there are only"
               " 6 geometry sides", side);

  if (x < 0 || x >= _num_x || y < 0 || y >= _num_y || z < 0 || z > _num_z)
    log_printf(ERROR, "Unable to get the neighbor cell for invalid cell "
               " (%i, %i, %i) since the dimesions are (%i, %i, %i)"
               , x, y, z, _num_x, _num_y, _num_z);

  int neighbor_cell = -1;

  if (side == 0){
    if (x != 0)
      neighbor_cell = z*_num_x*_num_y + y*_num_x + x - 1;
  }
  else if (side == 1){
    if (y != 0)
      neighbor_cell = z*_num_x*_num_y + (y-1)*_num_x + x;
  }
  else if (side == 2){
    if (z != 0)
      neighbor_cell = (z-1)*_num_x*_num_y + y*_num_x + x;
  }
  else if (side == 3){
    if (x != _num_x-1)
      neighbor_cell = z*_num_x*_num_y + y*_num_x + x + 1;
  }
  else if (side == 4){
    if (y != _num_y-1)
      neighbor_cell = z*_num_x*_num_y + (y+1)*_num_x + x;
  }
  else if (side == 5){
    if (z != _num_z-1)
      neighbor_cell = (z+1)*_num_x*_num_y + y*_num_x + x;
  }

  return neighbor_cell;
}


Material* AmpMesh::getNeighborMaterial(int x, int y, int z, int side){

  if (x < 0 || x >= _num_x || y < 0 || y >= _num_y || z < 0 || z > _num_z)
    log_printf(ERROR, "Unable to get the neighbor material for invalid cell "
               " (%i, %i, %i) since the dimesions are (%i, %i, %i)"
               , x, y, z, _num_x, _num_y, _num_z);

  int neighbor_cell = getNeighborCell(x, y, z, side);

  if (neighbor_cell == -1)
    return NULL;
  else
    return _materials[neighbor_cell];
  
}


int AmpMesh::getNumX(){
  return _num_x;
}


int AmpMesh::getNumY(){
  return _num_y;
}


int AmpMesh::getNumZ(){
  return _num_z;
}


double AmpMesh::getCellWidth(){
  return _cell_width;
}


double AmpMesh::getCellHeight(){
  return _cell_height;
}

double AmpMesh::getCellDepth(){
  return _cell_depth;
}


double AmpMesh::getCellVolume(int cell){
  return _cell_width * _cell_height * _cell_depth;
}


void AmpMesh::setNumX(int num_x){

  if (num_x < 1)
    log_printf(ERROR, "Unable to set num x for Amp Mesh to non"
               " positive number: %i", num_x);

  _num_x = num_x;
  _cell_width = (_x_max - _x_min) / num_x;
}


void AmpMesh::setNumY(int num_y){

  if (num_y < 1)
    log_printf(ERROR, "Unable to set num y for Amp Mesh to non"
               " positive number: %i", num_y);

  _num_y = num_y;
  _cell_height = (_y_max - _y_min) / num_y;
}


void AmpMesh::setNumZ(int num_z){

  if (num_z < 1)
    log_printf(ERROR, "Unable to set num z for Amp Mesh to non"
               " positive number: %i", num_z);

  _num_z = num_z;
  _cell_depth = (_z_max - _z_min) / num_z;
}


int AmpMesh::findCell(double x, double y, double z){

  if (x > _x_max || x < _x_min)
    log_printf(ERROR, "Unable to find cell with x value %f outside "
               "the mesh x bounds (%f, %f)", x, _x_min, _x_max);

  if (y > _y_max || y < _y_min)
    log_printf(ERROR, "Unable to find cell with y value %f outside "
               "the mesh y bounds (%f, %f)", y, _y_min, _y_max);

  if (z > _z_max || z < _z_min)
    log_printf(ERROR, "Unable to find cell with z value %f outside "
               "the mesh z bounds (%f, %f)", z, _z_min, _z_max);

  int i = floor((x - _x_min) / _cell_width);
  int j = floor((y - _y_min) / _cell_height);
  int k = floor((z - _z_min) / _cell_depth);

  return k * _num_x * _num_y + j * _num_x + i;  
}


double* AmpMesh::getCurrent(int position){

  if (position >= 8  || position < 0)
    log_printf(ERROR, "Unable to get the current for position %i "
               "since there are only 8 clock positions", position);

  if (_current[position] == NULL)
    log_printf(ERROR, "Unable to get the current for position %i "
               "since the current has not been initialized", position);

  return _current[position];
}


double* AmpMesh::getDifLinear(int position){

  if (position >= 8  || position < 0)
    log_printf(ERROR, "Unable to get the dif linear for position %i "
               "since there are only 8 clock positions", position);

  if (_dif_linear[position] == NULL)
    log_printf(ERROR, "Unable to get the dif linear for position %i "
               "since the dif linear has not been initialized", position);

  return _dif_linear[position];
}


double* AmpMesh::getDifNonlinear(int position){

  if (position >= 8  || position < 0)
    log_printf(ERROR, "Unable to get the dif linear for position %i "
               "since there are only 8 clock positions", position);

  if (_dif_nonlinear[position] == NULL)
    log_printf(ERROR, "Unable to get the dif nonlinear for position %i "
               "since the dif nonlinear has not been initialized", position);

  return _dif_nonlinear[position];
}


int AmpMesh::getNumCells(){
  return _num_x*_num_y*_num_z;
}
