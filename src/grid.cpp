// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "aether.h"

Grid::Grid(int nX_in, int nY_in, int nZ_in, int nGCs_in) {

  nX = nX_in + nGCs_in * 2; nLons = nX;
  nY = nY_in + nGCs_in * 2; nLats = nY;
  nZ = nZ_in + nGCs_in * 2; nAlts = nZ;
  nGCs = nGCs_in;

  int64_t nTotalPoints = int64_t(nX) * int64_t(nY) * int64_t(nZ);

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

  sza_scgc.set_size(nX, nY, nZ);
  cos_sza_scgc.set_size(nX, nY, nZ);

  fcube tmp(nX, nY, nZ);
  tmp.zeros();
  bfield_vcgc.push_back(tmp);  // x-component
  bfield_vcgc.push_back(tmp);  // y-component
  bfield_vcgc.push_back(tmp);  // z-component
  bfield_mag_scgc.set_size(nX, nY, nZ);
  bfield_mag_scgc.zeros();

  GSE_XYZ_vcgc.push_back(tmp);  // x-component
  GSE_XYZ_vcgc.push_back(tmp);  // y-component
  GSE_XYZ_vcgc.push_back(tmp);  // z-component

  mag_pole_north_ll.set_size(2);
  mag_pole_south_ll.set_size(2);
  mag_pole_north_ll.zeros();
  mag_pole_south_ll.zeros();

  fcube tmp_col(1,1,nZ);
  mag_pole_north_gse.push_back(tmp_col);
  mag_pole_north_gse.push_back(tmp_col);
  mag_pole_north_gse.push_back(tmp_col);

  mag_pole_south_gse.push_back(tmp_col);
  mag_pole_south_gse.push_back(tmp_col);
  mag_pole_south_gse.push_back(tmp_col);

}

int Grid::get_IsGeoGrid() {
  return IsGeoGrid;
}

void Grid::set_IsGeoGrid(int value) {
  IsGeoGrid = value;
}

int64_t Grid::get_nPointsInGrid() {
  int64_t nPoints;
  nPoints = int64_t(nX) * int64_t(nY) * int64_t(nZ);
  return nPoints;
}

int64_t Grid::get_nX() { return nX; }
int64_t Grid::get_nY() { return nY; }
int64_t Grid::get_nZ() { return nZ; }

int64_t Grid::get_nLons() { return nLons; }
int64_t Grid::get_nLats() { return nLats; }
int64_t Grid::get_nAlts() { return nAlts; }

int64_t Grid::get_nGCs() { return nGCs; }
