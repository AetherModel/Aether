// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
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

  std::string function = "Grid::calc_sza";
  static int iFunction = -1;
  report.enter(function, iFunction);

  float lon_offset = planet.get_longitude_offset(time);
  float declination = planet.get_declination(time);
  float sin_dec = planet.get_sin_dec(time);
  float cos_dec = planet.get_cos_dec(time);

  // Local time is in radians
  geoLocalTime_scgc = geoLon_scgc + lon_offset;
  cos_sza_scgc =
    sin_dec * sin(geoLat_scgc) +
    cos_dec * cos(geoLat_scgc) % cos(geoLocalTime_scgc-pi);
  sza_scgc = acos(cos_sza_scgc);

  // Compute GSE coordinates:
  // 1. use latitude / local time to derive XYZ
  // 2. rotate by declination to point x-axis towards the sun
  GSE_XYZ_vcgc[0] = radius_scgc % cos(geoLocalTime_scgc-pi) % cos(geoLat_scgc);
  GSE_XYZ_vcgc[1] = radius_scgc % sin(geoLocalTime_scgc-pi) % cos(geoLat_scgc);
  GSE_XYZ_vcgc[2] = radius_scgc % sin(geoLat_scgc);
  GSE_XYZ_vcgc = rotate_around_y_3d(GSE_XYZ_vcgc, -declination);
  
  report.exit(function);
}

// -----------------------------------------------------------------------------
//  fill grid with magnetic field values
// -----------------------------------------------------------------------------

void Grid::fill_grid_bfield(Planets planet, Inputs input, Report &report) {

  std::string function = "Grid::fill_grid_bfield";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t iLon, iLat, iAlt, iDim;
  float lon, lat, alt;
  bfield_info_type bfield_info;

  for (iLon = 0; iLon < nLons; iLon++) {
    for (iLat = 0; iLat < nLats; iLat++) {
      for (iAlt = 0; iAlt < nAlts; iAlt++) {

        lon = geoLon_scgc(iLon, iLat, iAlt);
        lat = geoLat_scgc(iLon, iLat, iAlt);
        alt = geoAlt_scgc(iLon, iLat, iAlt);

        bfield_info = get_bfield(lon, lat, alt, planet, input, report);

        magLat_scgc(iLon, iLat, iAlt) = bfield_info.lat;
        magLon_scgc(iLon, iLat, iAlt) = bfield_info.lon;

        for (iDim = 0; iDim < 3; iDim++)
          bfield_vcgc[iDim](iLon, iLat, iAlt) = bfield_info.b[iDim];
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

  int64_t iLon, iLat, iAlt;

  float radius0;
  float mu = planet.get_mu();

  report.print(3, "starting fill_grid_radius");

  // Just in case we have a latitude-dependent planetary radius
  fvec radius0_1d(nLats);
  for (iLat = 0; iLat < nLats; iLat++)
    radius0_1d(iLat) = planet.get_radius(geoLat_scgc(0, iLat, 0));

  for (iLon=0; iLon < nLons; iLon++)
    for (iAlt=0; iAlt < nAlts; iAlt++)
      radius_scgc.subcube(iLon, 0, iAlt, iLon, nLats-1, iAlt) = radius0_1d;
  radius_scgc = radius_scgc + geoAlt_scgc;

  radius2_scgc = radius_scgc % radius_scgc;
  radius2i_scgc = 1.0 / radius2_scgc;

  gravity_scgc = mu * radius2i_scgc;

  report.print(3, "ending fill_grid_radius");
}

// -----------------------------------------------------------------------------
//  Fill in XYZ in geo and mag coordinates
// -----------------------------------------------------------------------------

void Grid::fill_grid(Planets planet, Report &report) {

  int64_t iLon, iLat, iAlt;

  report.print(3, "starting fill_grid");

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

  transform_llr_to_xyz_3d(geoLon_scgc, geoLat_scgc, radius_scgc,
                          geoX_scgc, geoY_scgc, geoZ_scgc);
}
