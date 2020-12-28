// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>
#include <math.h>

#include "../include/inputs.h"
#include "../include/constants.h"
#include "../include/report.h"
#include "../include/grid.h"
#include "../include/sizes.h"
#include "../include/planets.h"
#include "../include/transform.h"
#include "../include/bfield.h"

// -----------------------------------------------------------------------------
//  Fill in Solar Zenith Angle and cos(solar zenith angle)
// -----------------------------------------------------------------------------

void Grid::calc_sza(Planets planet, Times time, Report &report) {

  long iLon, iLat, iAlt, index;
  float local_time;

  std::string function = "Grid::calc_sza";
  static int iFunction = -1;
  report.enter(function, iFunction);  

  float lon_offset = planet.get_longitude_offset(time);
  float sin_dec = planet.get_sin_dec(time);
  float cos_dec = planet.get_cos_dec(time);

  fcube local_time3d(nLons,nLats,nAlts);
  local_time3d = geoLon_scgc + lon_offset;
  cos_sza_scgc = 
    sin_dec * sin(geoLat_scgc) +
    cos_dec * cos(geoLat_scgc) % cos(local_time3d-pi);
  sza_scgc = acos(cos_sza_scgc);
  
  for (iLon = 0; iLon < nGeoLonsG; iLon++) {
    for (iLat = 0; iLat < nGeoLatsG; iLat++) {
      for (iAlt = 0; iAlt < nGeoAltsG; iAlt++) {

	index = ijk_geo_s3gc(iLon,iLat,iAlt);

	// This is in radians:
	local_time = fmod(lon_offset + geoLon_s3gc[index] + twopi, twopi);
	
	cos_sza_s3gc[index] =
	  sin_dec * sin(geoLat_s3gc[index]) +
	  cos_dec * cos(geoLat_s3gc[index]) * cos(local_time-pi);
	
	sza_s3gc[index] = acos(cos_sza_s3gc[index]);

      }
    }
  }

  report.exit(function);

}

// -----------------------------------------------------------------------------
//  fill grid with magnetic field values
// -----------------------------------------------------------------------------

void Grid::fill_grid_bfield(Planets planet, Inputs input, Report &report) {

  std::string function = "Grid::fill_grid_bfield";
  static int iFunction = -1;
  report.enter(function, iFunction);  

  long nLons, nLats, nAlts, iLon, iLat, iAlt, iDim, index, indexv;
  float lon, lat, alt;
  bfield_info_type bfield_info;
  
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
	} else {
	  index = ijk_mag_s3gc(iLon,iLat,iAlt);
	}

	lon = geoLon_s3gc[index];
	lat = geoLat_s3gc[index];
	alt = geoAlt_s3gc[index];
	
	bfield_info = get_bfield(lon, lat, alt, planet, input, report);

	magLat_s3gc[index] = bfield_info.lat;
	magLon_s3gc[index] = bfield_info.lon;

	for (iDim = 0; iDim < 3; iDim++) {
	  if (IsGeoGrid) {
	    indexv = ijkl_geo_v3gc(iLon,iLat,iAlt,iDim);
	  } else {
	    indexv = ijkl_mag_v3gc(iLon,iLat,iAlt,iDim);
	  }
	  bfield_v3gc[indexv] = bfield_info.b[iDim];
	}
	
      }
    }
  }
  
  report.exit(function);
  return;
}

// -----------------------------------------------------------------------------
//  Fill in radius, radius^2, and 1/radius^2
// -----------------------------------------------------------------------------

void Grid::fill_grid_radius(Planets planet, Report &report) {

  long iLon, iLat, iAlt, index;

  float radius0;
  float mu = planet.get_mu();

  report.print(3, "starting fill_grid_radius");

  // Just in case we have a latitude-dependent planetary radius
  fvec radius0_1d(nLats);
  for (iLat = 0; iLat < nLats; iLat++) 
    radius0_1d(iLat) = planet.get_radius(geoLat_scgc(0,iLat,0));
  
  for (iLon=0; iLon < nLons; iLon++) 
    for (iAlt=0; iAlt < nAlts; iAlt++) 
      radius_scgc.subcube(iLon,0,iAlt,iLon,nLats-1,iAlt) = radius0_1d;

  radius2_scgc = radius_scgc % radius_scgc;
  radius2i_scgc = 1.0 / radius2_scgc;

  gravity_scgc = mu * radius2i_scgc;
  
  for (iLon = 0; iLon < nLons; iLon++) {
    for (iLat = 0; iLat < nLats; iLat++) {
      for (iAlt = 0; iAlt < nAlts; iAlt++) {

	if (IsGeoGrid) {
	  index = ijk_geo_s3gc(iLon,iLat,iAlt);
	} else {
	  index = ijk_mag_s3gc(iLon,iLat,iAlt);
	}

	// Set radius and radius^2:

	radius0 = planet.get_radius(geoLat_s3gc[index]);

	radius_s3gc[index] =
	  radius0 + geoAlt_s3gc[index];
	radius_sq_s3gc[index] =
	  radius_s3gc[index] * radius_s3gc[index];
	radius_inv_sq_s3gc[index] =
	  1.0/radius_sq_s3gc[index];

	gravity_s3gc[index] = mu * radius_inv_sq_s3gc[index];

      }
    }
  }
 
  report.print(3, "ending fill_grid_radius");

}

// -----------------------------------------------------------------------------
//  Fill in XYZ in geo and mag coordinates
// -----------------------------------------------------------------------------

void Grid::fill_grid(Planets planet, Report &report) {

  long iLon, iLat, iAlt, index;

  long indexp, indexm;
  
  report.print(3, "starting fill_grid");
  
  //if (IsGeoGrid) {
  //  nLons = nGeoLonsG;
  //  nLats = nGeoLatsG;
  //  nAlts = nGeoAltsG;
  //} else {
  //  nLons = nMagLonsG;
  //  nLats = nMagLatsG;
  //  nAlts = nMagAltsG;
  //}

  float xyz[3], llr[3];

  for (iAlt = 1; iAlt < nAlts-1; iAlt++) {
    dalt_center_scgc.slice(iAlt) =
      geoAlt_scgc.slice(iAlt+1) - geoAlt_scgc.slice(iAlt-1);
    dalt_lower_scgc.slice(iAlt) =
      geoAlt_scgc.slice(iAlt) - geoAlt_scgc.slice(iAlt-1);
  }
  dalt_center_scgc.slice(0) = dalt_center_scgc.slice(1); 
  dalt_center_scgc.slice(nAlts-1) = dalt_center_scgc.slice(nAlts-2);

  dalt_lower_scgc.slice(0) = dalt_lower_scgc.slice(1); 
  iAlt = nAlts-1;
  dalt_lower_scgc.slice(iAlt) =
    geoAlt_scgc.slice(iAlt) - geoAlt_scgc.slice(iAlt-1);

  for (iLon = 0; iLon < nLons; iLon++) {
    for (iLat = 0; iLat < nLats; iLat++) {
      for (iAlt = 0; iAlt < nAlts; iAlt++) {

	if (IsGeoGrid) {
	  index = ijk_geo_s3gc(iLon,iLat,iAlt);
	} else {
	  index = ijk_mag_s3gc(iLon,iLat,iAlt);
	}

	// Find XYZ coordinates (Geo):
	
	llr[0] = geoLon_s3gc[index];
	llr[1] = geoLat_s3gc[index];
	llr[2] = geoAlt_s3gc[index];

	transform_llr_to_xyz(llr, xyz);

	geoX_s3gc[index] = xyz[0];
	geoY_s3gc[index] = xyz[1];
	geoZ_s3gc[index] = xyz[2];

	//// Find XYZ coordinates (Mag):
	//
	//llr[0] = magLon_s3gc[index];
	//llr[1] = magLat_s3gc[index];
	//llr[2] = magAlt_s3gc[index];
	//
	//transform_llr_to_xyz(llr, xyz);
	//
	//magX_s3gc[index] = xyz[0];
	//magY_s3gc[index] = xyz[1];
	//magZ_s3gc[index] = xyz[2];

	if (iAlt > 0 && iAlt < nAlts-1) {
	  if (IsGeoGrid) {
	    indexp = ijk_geo_s3gc(iLon,iLat,iAlt+1);
	    indexm = ijk_geo_s3gc(iLon,iLat,iAlt-1);
	  } else {
	    indexp = ijk_mag_s3gc(iLon,iLat,iAlt+1);
	    indexm = ijk_mag_s3gc(iLon,iLat,iAlt-1);
	  }
	  
	  dalt_center_s3gc[index] =
	    (geoAlt_s3gc[indexp] - geoAlt_s3gc[indexm])/2;
	  dalt_lower_s3gc[index] =
	    geoAlt_s3gc[index] - geoAlt_s3gc[indexm];
	  
	}
	
      }
    }
  }

  report.print(4, "filling edges");

  for (iLon = 0; iLon < nLons; iLon++) {
    for (iLat = 0; iLat < nLats; iLat++) {

      iAlt = 0;

      if (IsGeoGrid) {
	index = ijk_geo_s3gc(iLon,iLat,iAlt);
	indexp = ijk_geo_s3gc(iLon,iLat,iAlt+1);
      } else {
	index = ijk_mag_s3gc(iLon,iLat,iAlt);
	indexp = ijk_mag_s3gc(iLon,iLat,iAlt+1);
      }

      dalt_center_s3gc[index] = dalt_center_s3gc[indexp];
      dalt_lower_s3gc[index] = dalt_lower_s3gc[indexp];

      iAlt = nAlts-1;

      if (IsGeoGrid) {
	index = ijk_geo_s3gc(iLon,iLat,iAlt);
	indexm = ijk_geo_s3gc(iLon,iLat,iAlt-1);
      } else {
	index = ijk_mag_s3gc(iLon,iLat,iAlt);
	indexm = ijk_mag_s3gc(iLon,iLat,iAlt-1);
      }

      dalt_center_s3gc[index] = dalt_center_s3gc[indexm];
      dalt_lower_s3gc[index] = dalt_lower_s3gc[indexm];
      
    }
  }

  report.print(3, "ending fill_grid");

  
}
