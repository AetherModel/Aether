// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>
#include <armadillo>

#include "../include/inputs.h"
#include "../include/grid.h"
#include "../include/sizes.h"

Grid::Grid(int nX_in, int nY_in, int nZ_in, int nGCs_in) {

  nX = nX_in; nLons = nX;
  nY = nY_in; nLats = nY;
  nZ = nZ_in; nAlts = nZ;
  nGCs = nGCs_in;

  long nTotalPoints = long(nX) * long(nY) * long(nZ);

  geoLon_scgc.set_size(nX,nY,nZ);
  geoLat_scgc.set_size(nX,nY,nZ);
  geoAlt_scgc.set_size(nX,nY,nZ);

  radius_scgc.set_size(nX,nY,nZ);
  radius2_scgc.set_size(nX,nY,nZ);
  radius2i_scgc.set_size(nX,nY,nZ);

  dalt_center_scgc.set_size(nX,nY,nZ);
  dalt_lower_scgc.set_size(nX,nY,nZ);

  sza_scgc.set_size(nX,nY,nZ);
  cos_sza_scgc.set_size(nX,nY,nZ);
  
  geoLon_s3gc = (float*) malloc( nTotalPoints * sizeof(float) );
  geoLat_s3gc = (float*) malloc( nTotalPoints * sizeof(float) );
  geoAlt_s3gc = (float*) malloc( nTotalPoints * sizeof(float) );
  geoX_s3gc = (float*) malloc( nTotalPoints * sizeof(float) );
  geoY_s3gc = (float*) malloc( nTotalPoints * sizeof(float) );
  geoZ_s3gc = (float*) malloc( nTotalPoints * sizeof(float) );

  magLon_s3gc = (float*) malloc( nTotalPoints * sizeof(float) );
  magLat_s3gc = (float*) malloc( nTotalPoints * sizeof(float) );
  magAlt_s3gc = (float*) malloc( nTotalPoints * sizeof(float) );
  magX_s3gc = (float*) malloc( nTotalPoints * sizeof(float) );
  magY_s3gc = (float*) malloc( nTotalPoints * sizeof(float) );
  magZ_s3gc = (float*) malloc( nTotalPoints * sizeof(float) );

  magLocalTime_s3gc = (float*) malloc( nTotalPoints * sizeof(float) );

  radius_s3gc = (float*) malloc( nTotalPoints * sizeof(float) );
  radius_sq_s3gc = (float*) malloc( nTotalPoints * sizeof(float) );
  radius_inv_sq_s3gc = (float*) malloc( nTotalPoints * sizeof(float) );

  gravity_s3gc = (float*) malloc( nTotalPoints * sizeof(float) );

  sza_s3gc = (float*) malloc( nTotalPoints * sizeof(float) );
  cos_sza_s3gc = (float*) malloc( nTotalPoints * sizeof(float) );

  bfield_v3gc = (float*) malloc( 3 * nTotalPoints * sizeof(float) );
  bfield_mag_s3gc = (float*) malloc( nTotalPoints * sizeof(float) );

  dalt_center_s3gc = (float*) malloc( nTotalPoints * sizeof(float) );
  dalt_lower_s3gc = (float*) malloc( nTotalPoints * sizeof(float) );
  
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
