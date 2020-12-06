// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "../include/inputs.h"
#include "../include/grid.h"
#include "../include/sizes.h"

Grid::Grid(int nX, int nY, int nZ) {

  long iTotal = nX * nY * nZ;

  geoLon_s3gc = (float*) malloc( iTotal * sizeof(float) );
  geoLat_s3gc = (float*) malloc( iTotal * sizeof(float) );
  geoAlt_s3gc = (float*) malloc( iTotal * sizeof(float) );
  geoX_s3gc = (float*) malloc( iTotal * sizeof(float) );
  geoY_s3gc = (float*) malloc( iTotal * sizeof(float) );
  geoZ_s3gc = (float*) malloc( iTotal * sizeof(float) );

  magLon_s3gc = (float*) malloc( iTotal * sizeof(float) );
  magLat_s3gc = (float*) malloc( iTotal * sizeof(float) );
  magAlt_s3gc = (float*) malloc( iTotal * sizeof(float) );
  magX_s3gc = (float*) malloc( iTotal * sizeof(float) );
  magY_s3gc = (float*) malloc( iTotal * sizeof(float) );
  magZ_s3gc = (float*) malloc( iTotal * sizeof(float) );

  magLocalTime_s3gc = (float*) malloc( iTotal * sizeof(float) );

  radius_s3gc = (float*) malloc( iTotal * sizeof(float) );
  radius_sq_s3gc = (float*) malloc( iTotal * sizeof(float) );
  radius_inv_sq_s3gc = (float*) malloc( iTotal * sizeof(float) );

  gravity_s3gc = (float*) malloc( iTotal * sizeof(float) );

  sza_s3gc = (float*) malloc( iTotal * sizeof(float) );
  cos_sza_s3gc = (float*) malloc( iTotal * sizeof(float) );

  dalt_center_s3gc = (float*) malloc( iTotal * sizeof(float) );
  dalt_lower_s3gc = (float*) malloc( iTotal * sizeof(float) );
  
}

int Grid::get_IsGeoGrid() {
  return IsGeoGrid;
}

void Grid::set_IsGeoGrid(int value) {
  IsGeoGrid = value;
}
