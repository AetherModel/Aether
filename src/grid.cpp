// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "aether.h"

// --------------------------------------------------------------------------
// Initialize Grid class
// --------------------------------------------------------------------------

Grid::Grid(std::string gridtype) {

  // At this point, we only need 2 ghostcells.  Hardcode this:
  nGCs = 2;

  Inputs::grid_input_struct grid_input = input.get_grid_inputs(gridtype);

  nX = grid_input.nX + nGCs * 2;
  nLons = nX;
  nY = grid_input.nY + nGCs * 2;
  nLats = nY;
  nZ = grid_input.nZ + nGCs * 2;
  nAlts = nZ;

  // No set all of the logicals to make the flow a bit easier:

  if (grid_input.nX == 1 &
      grid_input.nY == 1 &
      grid_input.nZ == 1)
    Is0D = true;
  else {
    if (grid_input.nY == 1 & grid_input.nZ == 1) Is1Dx = true;
    if (grid_input.nX == 1 & grid_input.nZ == 1) Is1Dy = true;
    if (grid_input.nX == 1 & grid_input.nY == 1) Is1Dz = true;
    if (!Is1Dx & !Is1Dy & !Is1Dz) {
      if (grid_input.nX == 1) Is2Dyz = true;
      if (grid_input.nY == 1) Is2Dxz = true;
      if (grid_input.nZ == 1) Is2Dxy = true;
      if (!Is2Dyz & !Is2Dxz & !Is2Dxy) Is3D = true;
    }
  }

  if (grid_input.nX == 1) HasXdim = false;
  if (grid_input.nY == 1) HasYdim = false;
  if (grid_input.nZ == 1) HasZdim = false;

  if (mklower(grid_input.shape) == "sphere") 
    iGridShape_ = iSphere_;
  if (mklower(grid_input.shape) == "cubesphere") 
    iGridShape_ = iCubesphere_;
  if (mklower(grid_input.shape) == "dipole") 
    iGridShape_ = iDipole_;

  geoLon_scgc.set_size(nX, nY, nZ);
  geoLat_scgc.set_size(nX, nY, nZ);
  geoAlt_scgc.set_size(nX, nY, nZ);
  geoLocalTime_scgc.set_size(nX, nY, nZ);

  refx_scgc.set_size(nX, nY, nZ);
  refy_scgc.set_size(nX, nY, nZ);
  refx_angle.set_size(nX, nY, nZ);
  refy_angle.set_size(nX, nY, nZ);

  A11_scgc.set_size(nX, nY, nZ);
  A12_scgc.set_size(nX, nY, nZ);
  A21_scgc.set_size(nX, nY, nZ);
  A22_scgc.set_size(nX, nY, nZ);
  A11_inv_scgc.set_size(nX, nY, nZ);
  A12_inv_scgc.set_size(nX, nY, nZ);
  A21_inv_scgc.set_size(nX, nY, nZ);
  A22_inv_scgc.set_size(nX, nY, nZ);
  g11_upper_scgc.set_size(nX, nY, nZ);
  g12_upper_scgc.set_size(nX, nY, nZ);
  g21_upper_scgc.set_size(nX, nY, nZ);
  g22_upper_scgc.set_size(nX, nY, nZ);
  sqrt_g_scgc.set_size(nX, nY, nZ);

  geoLon_Left.set_size(nX + 1, nY, nZ);
  geoLat_Left.set_size(nX + 1, nY, nZ);

  refx_Left.set_size(nX + 1, nY, nZ);
  refy_Left.set_size(nX + 1, nY, nZ);

  A11_Left.set_size(nX + 1, nY, nZ);
  A12_Left.set_size(nX + 1, nY, nZ);
  A21_Left.set_size(nX + 1, nY, nZ);
  A22_Left.set_size(nX + 1, nY, nZ);
  A11_inv_Left.set_size(nX + 1, nY, nZ);
  A12_inv_Left.set_size(nX + 1, nY, nZ);
  A21_inv_Left.set_size(nX + 1, nY, nZ);
  A22_inv_Left.set_size(nX + 1, nY, nZ);
  g11_upper_Left.set_size(nX + 1, nY, nZ);
  g12_upper_Left.set_size(nX + 1, nY, nZ);
  g21_upper_Left.set_size(nX + 1, nY, nZ);
  g22_upper_Left.set_size(nX + 1, nY, nZ);
  sqrt_g_Left.set_size(nX + 1, nY, nZ);

  geoLon_Down.set_size(nX, nY + 1, nZ);
  geoLat_Down.set_size(nX, nY + 1, nZ);

  refx_Down.set_size(nX, nY + 1, nZ);
  refy_Down.set_size(nX, nY + 1, nZ);

  A11_Down.set_size(nX, nY + 1, nZ);
  A12_Down.set_size(nX, nY + 1, nZ);
  A21_Down.set_size(nX, nY + 1, nZ);
  A22_Down.set_size(nX, nY + 1, nZ);
  A11_inv_Down.set_size(nX, nY + 1, nZ);
  A12_inv_Down.set_size(nX, nY + 1, nZ);
  A21_inv_Down.set_size(nX, nY + 1, nZ);
  A22_inv_Down.set_size(nX, nY + 1, nZ);
  g11_upper_Down.set_size(nX, nY + 1, nZ);
  g12_upper_Down.set_size(nX, nY + 1, nZ);
  g21_upper_Down.set_size(nX, nY + 1, nZ);
  g22_upper_Down.set_size(nX, nY + 1, nZ);
  sqrt_g_Down.set_size(nX, nY + 1, nZ);

  geoLon_Corner.set_size(nX + 1, nY + 1, nZ + 1);
  geoLat_Corner.set_size(nX + 1, nY + 1, nZ + 1);
  geoAlt_Corner.set_size(nX + 1, nY + 1, nZ + 1);

  refx_Corner.set_size(nX + 1, nY + 1, nZ + 1);
  refy_Corner.set_size(nX + 1, nY + 1, nZ + 1);

  geoAlt_Below.set_size(nX, nY, nZ + 1);

  geoX_scgc.set_size(nX, nY, nZ);
  geoY_scgc.set_size(nX, nY, nZ);
  geoZ_scgc.set_size(nX, nY, nZ);

  magLon_scgc.set_size(nX, nY, nZ);
  magLat_scgc.set_size(nX, nY, nZ);
  magAlt_scgc.set_size(nX, nY, nZ);

  magPhi_scgc.set_size(nX, nY, nZ);
  magP_scgc.set_size(nX, nY, nZ);
  magQ_scgc.set_size(nX, nY, nZ);

  magX_scgc.set_size(nX, nY, nZ);
  magY_scgc.set_size(nX, nY, nZ);
  magZ_scgc.set_size(nX, nY, nZ);

  magLocalTime_scgc.set_size(nX, nY, nZ);

  radius_scgc.set_size(nX, nY, nZ);
  radius2_scgc.set_size(nX, nY, nZ);
  radius2i_scgc.set_size(nX, nY, nZ);

  dalt_center_scgc.set_size(nX, nY, nZ);
  dalt_lower_scgc.set_size(nX, nY, nZ);
  dalt_ratio_scgc.set_size(nX, nY, nZ);
  dalt_ratio_sq_scgc.set_size(nX, nY, nZ);

  MeshCoef1s3rdp1.set_size(nX, nY, nGCs);
  MeshCoef1s3rdp2.set_size(nX, nY, nGCs);
  MeshCoef1s3rdp3.set_size(nX, nY, nGCs);
  MeshCoef1s3rdp4.set_size(nX, nY, nGCs);
  MeshCoef1s3rdp5.set_size(nX, nY, nGCs);

  dlat_center_scgc.set_size(nX, nY, nZ);
  dlat_center_dist_scgc.set_size(nX, nY, nZ);

  dlon_center_scgc.set_size(nX, nY, nZ);
  dlon_center_dist_scgc.set_size(nX, nY, nZ);

  sza_scgc.set_size(nX, nY, nZ);
  cos_sza_scgc.set_size(nX, nY, nZ);

  bfield_vcgc = make_cube_vector(nX, nY, nZ, 3);
  bfield_unit_vcgc = make_cube_vector(nX, nY, nZ, 3);
  bfield_mag_scgc.set_size(nX, nY, nZ);
  bfield_mag_scgc.zeros();

  GSE_XYZ_vcgc = make_cube_vector(nX, nY, nZ, 3);

  mag_pole_north_ll.set_size(2);
  mag_pole_south_ll.set_size(2);
  mag_pole_north_ll.zeros();
  mag_pole_south_ll.zeros();

  arma_cube tmp_col(1, 1, nZ);
  mag_pole_north_gse.push_back(tmp_col);
  mag_pole_north_gse.push_back(tmp_col);
  mag_pole_north_gse.push_back(tmp_col);

  mag_pole_south_gse.push_back(tmp_col);
  mag_pole_south_gse.push_back(tmp_col);
  mag_pole_south_gse.push_back(tmp_col);

  HasBField = 0;
  IsExperimental = false;

  cent_acc_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);
  for (int i = 0; i < 3; i++)
    cent_acc_vcgc[i].zeros();

}

// --------------------------------------------------------------------------
// Set Variable Sizes
// --------------------------------------------------------------------------

void Grid::set_variable_sizes() {

  return;
}


// --------------------------------------------------------------------------
// write restart out files for the grid
// --------------------------------------------------------------------------

bool Grid::write_restart(std::string dir) {
  bool DidWork = true;

  // All Ensemble member grids should be the same, so only need to write
  // out the 0th member
  if (iMember == 0) {
    try {
      OutputContainer RestartContainer;
      RestartContainer.set_directory(dir);
      RestartContainer.set_version(aether_version);
      RestartContainer.set_time(0.0);
      RestartContainer.set_nGhostCells(nGCs);

      // Output Cell Centers
      RestartContainer.set_filename("grid_" + cGrid);
      RestartContainer.store_variable(longitude_name,
                                      longitude_unit,
                                      geoLon_scgc);
      RestartContainer.store_variable(latitude_name,
                                      latitude_unit,
                                      geoLat_scgc);
      RestartContainer.store_variable(altitude_name,
                                      altitude_unit,
                                      geoAlt_scgc);
      RestartContainer.write();
      RestartContainer.clear_variables();

      // Output Corners
      RestartContainer.set_filename("grid_corners_" + cGrid);
      RestartContainer.store_variable(longitude_name + " Corners",
                                      longitude_unit,
                                      geoLon_Corner);
      RestartContainer.store_variable(latitude_name + " Corners",
                                      latitude_unit,
                                      geoLat_Corner);
      RestartContainer.store_variable(altitude_name + " Corners",
                                      altitude_unit,
                                      geoAlt_Corner);
      RestartContainer.write();
      RestartContainer.clear_variables();

      // Output Left Sides
      RestartContainer.set_filename("grid_left_" + cGrid);
      RestartContainer.store_variable(longitude_name + " Left",
                                      longitude_unit,
                                      geoLon_Left);
      RestartContainer.store_variable(latitude_name + " Left",
                                      latitude_unit,
                                      geoLat_Left);
      RestartContainer.write();
      RestartContainer.clear_variables();

      // Output Down Sides
      RestartContainer.set_filename("grid_down_" + cGrid);
      RestartContainer.store_variable(longitude_name + " Down",
                                      longitude_unit,
                                      geoLon_Down);
      RestartContainer.store_variable(latitude_name + " Down",
                                      latitude_unit,
                                      geoLat_Down);
      RestartContainer.write();
      RestartContainer.clear_variables();

      // Output Below
      RestartContainer.set_filename("grid_below_" + cGrid);
      RestartContainer.store_variable(altitude_name + " Below",
                                      altitude_unit,
                                      geoAlt_Below);

      RestartContainer.write();
      RestartContainer.clear_variables();

    } catch (...) {
      std::cout << "Error writing grid restart file!\n";
      DidWork = false;
    }
  }

  DidWork = sync_across_all_procs(DidWork);
  return DidWork;
}

// --------------------------------------------------------------------------
// read restart out files for the grid
// - Returns true if everything worked ok
// --------------------------------------------------------------------------

bool Grid::read_restart(std::string dir) {

  bool DidWork = true;
  int64_t iVar;

  // While only the 0th ensemble member writes, all ensemble members have
  // to read the grid.

  try {
    OutputContainer RestartContainer;
    RestartContainer.set_directory(dir);
    RestartContainer.set_version(aether_version);
    // Cell Centers:
    RestartContainer.set_filename("grid_" + cGrid);
    RestartContainer.read();
    geoLon_scgc = RestartContainer.get_element_value(longitude_name);
    geoLat_scgc = RestartContainer.get_element_value(latitude_name);
    geoAlt_scgc = RestartContainer.get_element_value(altitude_name);
    // Down Edges:
    RestartContainer.set_filename("grid_below_" + cGrid);
    RestartContainer.read();
    geoAlt_Below = RestartContainer.get_element_value(altitude_name +
                                                      " Below");
    // Cell Corners:
    RestartContainer.set_filename("grid_corners_" + cGrid);
    RestartContainer.read();
    geoLon_Corner = RestartContainer.get_element_value(longitude_name +
                                                       " Corners");
    geoLat_Corner = RestartContainer.get_element_value(latitude_name +
                                                       " Corners");
    geoAlt_Corner = RestartContainer.get_element_value(altitude_name +
                                                       " Corners");
    // Left Edges:
    RestartContainer.set_filename("grid_left_" + cGrid);
    RestartContainer.read();
    geoLon_Left = RestartContainer.get_element_value(longitude_name +
                                                     " Left");
    geoLat_Left = RestartContainer.get_element_value(latitude_name +
                                                     " Left");
    // Down Edges:
    RestartContainer.set_filename("grid_down_" + cGrid);
    RestartContainer.read();
    geoLon_Down = RestartContainer.get_element_value(longitude_name +
                                                     " Down");
    geoLat_Down = RestartContainer.get_element_value(latitude_name +
                                                     " Down");

  } catch (...) {
    std::cout << "Error reading grid restart file!\n";
    DidWork = false;
  }

  DidWork = sync_across_all_procs(DidWork);
  return DidWork;
}

// --------------------------------------------------------------------------
// Report Grid boundaries:
// --------------------------------------------------------------------------

void Grid::report_grid_boundaries() {
  std::cout << "---------------------------------------------------\n";
  std::cout << "Grid Boundaries (min / max):\n";
  std::cout << "Lon : "
            << geoLon_scgc.min() << " / "
            << geoLon_scgc.max() << "\n";
  std::cout << "Lat : "
            << geoLat_scgc.min() << " / "
            << geoLat_scgc.max() << "\n";
  std::cout << "Alt : "
            << geoAlt_scgc.min() << " / "
            << geoAlt_scgc.max() << "\n";
  std::cout << "---------------------------------------------------\n";
}

// --------------------------------------------------------------------------
// Get whether the grid is a geographic grid (or magnetic - return 0)
// --------------------------------------------------------------------------

bool Grid::get_IsGeoGrid() {
  return IsGeoGrid;
}

// --------------------------------------------------------------------------
// Get whether the grid is a experimental (return true for experimental)
// --------------------------------------------------------------------------

bool Grid::get_IsExperimental() {
  return IsExperimental;
}

// --------------------------------------------------------------------------
// Get whether the grid is a geographic grid (or magnetic - return 0)
// --------------------------------------------------------------------------

bool Grid::get_HasBField() {
  return HasBField;
}

// --------------------------------------------------------------------------
// Set whether the grid is a geographic grid (or magnetic - set to 0)
// --------------------------------------------------------------------------

void Grid::set_IsGeoGrid(bool value) {
  IsGeoGrid = value;
}

// --------------------------------------------------------------------------
// Set whether the grid is an experimental grid 
// --------------------------------------------------------------------------

void Grid::set_IsExperimental(bool value) {
  IsExperimental = value;
}

// --------------------------------------------------------------------------
// Set whether the grid is a dipole grid 
// --------------------------------------------------------------------------

void Grid::set_IsDipole(bool value) {
  IsDipole = value;
}

// --------------------------------------------------------------------------
// Get total number of grid points
// --------------------------------------------------------------------------

int64_t Grid::get_nPointsInGrid() {
  int64_t nPoints;
  nPoints = int64_t(nX) * int64_t(nY) * int64_t(nZ);
  return nPoints;
}

// --------------------------------------------------------------------------
// Get some grid definition things
// --------------------------------------------------------------------------

bool Grid::get_HasXdim() {
  return HasXdim;
}

bool Grid::get_HasYdim() {
  return HasYdim;
}

bool Grid::get_HasZdim() {
  return HasZdim;
}

bool Grid::get_Is0D() {
  return Is0D;
}

bool Grid::get_Is1Dx() {
  return Is1Dx;
}

bool Grid::get_Is1Dy() {
  return Is1Dy;
}

bool Grid::get_Is1Dz() {
  return Is1Dz;
}

// --------------------------------------------------------------------------
// Get a bunch of sizes within the grid
// --------------------------------------------------------------------------

int64_t Grid::get_nX() {
  return nX;
}
int64_t Grid::get_nY() {
  return nY;
}
int64_t Grid::get_nZ() {
  return nZ;
}

int64_t Grid::get_nX(bool includeGCs) {
  if (includeGCs)
    return nX;
  else
    return nX - 2*nGCs;
}
int64_t Grid::get_nY(bool includeGCs) {
  if (includeGCs)
    return nY;
  else
    return nY - 2*nGCs;
}
int64_t Grid::get_nZ(bool includeGCs) {
  if (includeGCs)
    return nZ;
  else
    return nZ - 2*nGCs;
}

int64_t Grid::get_nLons() {
  return nLons;
}
int64_t Grid::get_nLats() {
  return nLats;
}
int64_t Grid::get_nAlts() {
  return nAlts;
}

int64_t Grid::get_nLons(bool includeGCs) {
  if (includeGCs)
    return nLons;
  else
    return nLons - 2*nGCs;
}
int64_t Grid::get_nLats(bool includeGCs) {
  if (includeGCs)
    return nLats;
  else
    return nLats - 2*nGCs;
}
int64_t Grid::get_nAlts(bool includeGCs) {
  if (includeGCs)
    return nAlts;
  else
    return nAlts - 2*nGCs;
}

int64_t Grid::get_nGCs() {
  return nGCs;
}
