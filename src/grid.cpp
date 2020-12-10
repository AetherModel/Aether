// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "../include/inputs.h"
#include "../include/grid.h"
#include "../include/sizes.h"

Grid::Grid(int nX, int nY, int nZ) {

  long nTotalPoints = long(nX) * long(nY) * long(nZ);

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
  if (IsGeoGrid)
    nPoints = long(nGeoLons) * long(nGeoLats) * long(nGeoAlts);
  else
    nPoints = long(nMagLons) * long(nMagLats) * long(nMagAlts);
  return nPoints;
}
