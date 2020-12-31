// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <math.h>
#include <vector>
#include <armadillo>

#include "../include/sizes.h"
#include "../include/grid.h"

using namespace arma;


// -----------------------------------------------------------------------
// copy from armidillo cube to 3d c-native array
// -----------------------------------------------------------------------

void copy_cube_to_array(fcube cube_in,
			float *array_out) {

  long nX = cube_in.n_rows;
  long nY = cube_in.n_cols;
  long nZ = cube_in.n_slices;
  long iX, iY, iZ, index;

  for (iX = 0; iX < nX; iX++) {
    for (iY = 0; iY < nY; iY++) {
      for (iZ = 0; iZ < nZ; iZ++) {
	index = iX*nY*nZ + iY*nZ + iZ;
	array_out[index] = cube_in(iX,iY,iZ);
      }
    }
  }
  
}



// -----------------------------------------------------------------------
// Transform Longitude, Latitude, Radius to X, Y, Z
// -----------------------------------------------------------------------

void transform_llr_to_xyz(float llr_in[3], float xyz_out[3]) {

  // llr_in[0] = longitude (in radians)
  // llr_in[1] = latitude (in radians)
  // llr_in[2] = radius
  
  xyz_out[0] = llr_in[2] * cos(llr_in[1]) * cos(llr_in[0]);
  xyz_out[1] = llr_in[2] * cos(llr_in[1]) * sin(llr_in[0]);
  xyz_out[2] = llr_in[2] * sin(llr_in[1]);
  
}

// -----------------------------------------------------------------------
// Rotate around the z-axis
//  - Angle needs to be in radians!!!
// -----------------------------------------------------------------------

void transform_rot_z(float xyz_in[3], float angle_in, float xyz_out[3]) {

  float ca = cos(angle_in);
  float sa = sin(angle_in);

  xyz_out[0] =  xyz_in[0] * ca + xyz_in[1] * sa;
  xyz_out[1] = -xyz_in[0] * sa + xyz_in[1] * ca;
  xyz_out[2] =  xyz_in[2];

  return;
  
}

// -----------------------------------------------------------------------
// Rotate around the y-axis
//  - Angle needs to be in radians!!!
// -----------------------------------------------------------------------

void transform_rot_y(float xyz_in[3], float angle_in, float xyz_out[3]) {

  float ca = cos(angle_in);
  float sa = sin(angle_in);

  xyz_out[0] = xyz_in[0] * ca - xyz_in[2] * sa;
  xyz_out[1] = xyz_in[1];
  xyz_out[2] = xyz_in[0] * sa + xyz_in[2] * ca;

  return;
  
}

// -----------------------------------------------------------------------
// Simply move data from a vector to a C-native array (float)
// -----------------------------------------------------------------------

void transform_float_vector_to_array(std::vector<float> input,
				     float output[3]) {

  for (int i=0; i<3; i++) output[i] = input[i];

  return;
}

// -----------------------------------------------------------------------
// Rotate a vector from XYZ to East, North, Vertical
// -----------------------------------------------------------------------

void transform_vector_xyz_to_env(float xyz_in[3],
				 float lon,
				 float lat,
				 float env_out[3]) {

  env_out[2] =   xyz_in[0] * cos(lat)*cos(lon) + xyz_in[1] * cos(lat)*sin(lon) + xyz_in[2]*sin(lat);
  env_out[1] = -(xyz_in[0] * sin(lat)*cos(lon) + xyz_in[1] * sin(lat)*sin(lon) - xyz_in[2]*cos(lat));
  env_out[0] = - xyz_in[0] * sin(lon)          + xyz_in[1] * cos(lon);
  
  return;

}

// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
// These are not realy transform functions.  I need a better place to
// put them.
// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
// Simple 3-element vector difference
// -----------------------------------------------------------------------

void vector_diff(float vect_in_1[3],
		 float vect_in_2[3],
		 float vect_out[3]) {
  
  for (int i=0; i<3; i++) vect_out[i] = vect_in_1[i] - vect_in_2[i];

  return;

}

// -----------------------------------------------------------------------
// grab one component of a vector
// -----------------------------------------------------------------------

void get_vector_component(float *vector_in_v3gc,
			  int iComponent,
			  int IsGeoGrid,
			  float *component_out_s3gc) {

  long nLons, nLats, nAlts, iLon, iLat, iAlt, index, indexv;

  if (IsGeoGrid) {
    nLons = nGeoLonsG;
    nLats = nGeoLatsG;
    nAlts = nGeoAltsG;
  } else {
    nLons = nMagLonsG;
    nLats = nMagLatsG;
    nAlts = nMagAltsG;
  }

  for (iLon = 0; iLon < nLons; iLon++) {
    for (iLat = 0; iLat < nLats; iLat++) {
      for (iAlt = 0; iAlt < nAlts; iAlt++) {

	if (IsGeoGrid) {
	  index = ijk_geo_s3gc(iLon,iLat,iAlt);
	  indexv = ijkl_geo_v3gc(iLon,iLat,iAlt,iComponent);
	} else {
	  index = ijk_mag_s3gc(iLon,iLat,iAlt);
	  indexv = ijkl_mag_v3gc(iLon,iLat,iAlt,iComponent);
	}

	component_out_s3gc[index] = vector_in_v3gc[indexv];

      }
    }
  }

  return;
  
}
