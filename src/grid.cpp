// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>
#include <armadillo>

#include "../include/inputs.h"
#include "../include/grid.h"
#include "../include/sizes.h"

using namespace arma;

Grid::Grid(int nX_in, int nY_in, int nZ_in, int nGCs_in) {

  nX = nX_in; nLons = nX;
  nY = nY_in; nLats = nY;
  nZ = nZ_in; nAlts = nZ;
  nGCs = nGCs_in;

  long nTotalPoints = long(nX) * long(nY) * long(nZ);

  geoLon_scgc.set_size(nX,nY,nZ);
  geoLat_scgc.set_size(nX,nY,nZ);
  geoAlt_scgc.set_size(nX,nY,nZ);

  geoX_scgc.set_size(nX,nY,nZ);
  geoY_scgc.set_size(nX,nY,nZ);
  geoZ_scgc.set_size(nX,nY,nZ);

  magLon_scgc.set_size(nX,nY,nZ);
  magLat_scgc.set_size(nX,nY,nZ);
  magAlt_scgc.set_size(nX,nY,nZ);

  magX_scgc.set_size(nX,nY,nZ);
  magY_scgc.set_size(nX,nY,nZ);
  magZ_scgc.set_size(nX,nY,nZ);

  magLocalTime_scgc.set_size(nX,nY,nZ);

  radius_scgc.set_size(nX,nY,nZ);
  radius2_scgc.set_size(nX,nY,nZ);
  radius2i_scgc.set_size(nX,nY,nZ);

  dalt_center_scgc.set_size(nX,nY,nZ);
  dalt_lower_scgc.set_size(nX,nY,nZ);

  sza_scgc.set_size(nX,nY,nZ);
  cos_sza_scgc.set_size(nX,nY,nZ);
  
  bfield_v3gc = (float*) malloc( 3 * nTotalPoints * sizeof(float) );
  bfield_mag_s3gc = (float*) malloc( nTotalPoints * sizeof(float) );

}

int Grid::get_IsGeoGrid() {
  return IsGeoGrid;
}

void Grid::set_IsGeoGrid(int value) {
  IsGeoGrid = value;
}

long Grid::get_nPointsInGrid() {
  long nPoints;
  nPoints = long(nX) * long(nY) * long(nZ);
  return nPoints;
}

long Grid::get_nX() { return nX; }
long Grid::get_nY() { return nY; }
long Grid::get_nZ() { return nZ; }

long Grid::get_nLons() { return nLons; }
long Grid::get_nLats() { return nLats; }
long Grid::get_nAlts() { return nAlts; }

long Grid::get_nGCs() { return nGCs; }
