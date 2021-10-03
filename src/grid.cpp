// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "aether.h"

// --------------------------------------------------------------------------
// Initialize Grid class
// --------------------------------------------------------------------------

Grid::Grid(int nX_in, int nY_in, int nZ_in, int nGCs_in) {

  nX = nX_in + nGCs_in * 2;
  nLons = nX;
  nY = nY_in + nGCs_in * 2;
  nLats = nY;
  nZ = nZ_in + nGCs_in * 2;
  nAlts = nZ;
  nGCs = nGCs_in;

  geoLon_scgc.set_size(nX, nY, nZ);
  geoLat_scgc.set_size(nX, nY, nZ);
  geoAlt_scgc.set_size(nX, nY, nZ);
  geoLocalTime_scgc.set_size(nX, nY, nZ);

  geoX_scgc.set_size(nX, nY, nZ);
  geoY_scgc.set_size(nX, nY, nZ);
  geoZ_scgc.set_size(nX, nY, nZ);

  magLon_scgc.set_size(nX, nY, nZ);
  magLat_scgc.set_size(nX, nY, nZ);
  magAlt_scgc.set_size(nX, nY, nZ);

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
}

// --------------------------------------------------------------------------
// write restart out files for the grid
// --------------------------------------------------------------------------

bool Grid::write_restart(std::string dir) {
  std::string filename;
  bool DidWork = true;
  filename = dir+"/geolons.bin";
  if (DidWork) DidWork = geoLon_scgc.save(filename);
  filename = dir+"/geolats.bin";
  if (DidWork) DidWork = geoLat_scgc.save(filename);
  filename = dir+"/geoalts.bin";
  if (DidWork) DidWork = geoAlt_scgc.save(filename);
  return DidWork;
}

// --------------------------------------------------------------------------
// read restart out files for the grid
// - Returns true if everything worked ok
// --------------------------------------------------------------------------

bool Grid::read_restart(std::string dir) {
  std::string filename;
  bool DidWork = true;
  filename = dir+"/geolons.bin";
  if (DidWork) DidWork = geoLon_scgc.load(filename);
  filename = dir+"/geolats.bin";
  if (DidWork) DidWork = geoLat_scgc.load(filename);
  filename = dir+"/geoalts.bin";
  if (DidWork) DidWork = geoAlt_scgc.load(filename);
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

int Grid::get_IsGeoGrid() {
  return IsGeoGrid;
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

void Grid::set_IsGeoGrid(int value) {
  IsGeoGrid = value;
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

int64_t Grid::get_nLons() {
  return nLons;
}
int64_t Grid::get_nLats() {
  return nLats;
}
int64_t Grid::get_nAlts() {
  return nAlts;
}

int64_t Grid::get_nGCs() {
  return nGCs;
}
