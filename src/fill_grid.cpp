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

// ---------------------------------------------------------------------------
//  Fill in Solar Zenith Angle and cos(solar zenith angle)
// ---------------------------------------------------------------------------

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
    cos_dec * cos(geoLat_scgc) % cos(geoLocalTime_scgc);
  sza_scgc = acos(cos_sza_scgc);

  report.exit(function);
}

// ---------------------------------------------------------------------------
//  Fill in GSE Coordinates
// ---------------------------------------------------------------------------

void Grid::calc_gse(Planets planet, Times time, Report &report) {

  std::string function = "Grid::calc_sza";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // Compute GSE coordinates:
  // 1. use latitude / local time to derive XYZ
  std::vector<fcube> lon_lat_radius;
  // -pi is used here, since +X is pointed towards sun:
  lon_lat_radius.push_back(geoLocalTime_scgc-pi);
  lon_lat_radius.push_back(geoLat_scgc);
  lon_lat_radius.push_back(radius_scgc);
  GSE_XYZ_vcgc = transform_llr_to_xyz_3d(lon_lat_radius);
  // 2. rotate by declination to point x-axis towards the sun
  float declination = planet.get_declination(time);
  GSE_XYZ_vcgc = rotate_around_y_3d(GSE_XYZ_vcgc, -declination);

  // ---------------------------------------------------------------
  // Do the same thing for the magnetic poles:

  std::vector<fcube> lon_lat_radius_col;
  fcube tmp_col(1,1,nZ);

  // North:
  // Lon (converted to local time):

  float lon_offset = planet.get_longitude_offset(time);
  tmp_col.fill(mag_pole_north_ll[0] + lon_offset - pi);
  lon_lat_radius_col.push_back(tmp_col);
  // Lat:
  tmp_col.fill(mag_pole_north_ll[1]);
  lon_lat_radius_col.push_back(tmp_col);
  // Radius:
  tmp_col.tube(0,0) = radius_scgc.tube(nX/2,nY/2);
  lon_lat_radius_col.push_back(tmp_col);
  
  mag_pole_north_gse = transform_llr_to_xyz_3d(lon_lat_radius_col);
  mag_pole_north_gse = rotate_around_y_3d(mag_pole_north_gse, -declination);

  // South:
  // Lon:
  lon_lat_radius_col[0].fill(mag_pole_south_ll[0] + lon_offset - pi);
  // Lat:
  lon_lat_radius_col[1].fill(mag_pole_south_ll[1]);
  // Radius is already filled.
  mag_pole_south_gse = transform_llr_to_xyz_3d(lon_lat_radius_col);
  mag_pole_south_gse = rotate_around_y_3d(mag_pole_south_gse, -declination);

  report.exit(function);
}

// ---------------------------------------------------------------------------
//  Fill in Magnetic Local Time
//    - This assumes that the GE coordinates have been calculated
//      and filled in!
// ---------------------------------------------------------------------------

void Grid::calc_mlt(Report &report) {

  std::string function = "Grid::calc_mlt";
  static int iFunction = -1;
  report.enter(function, iFunction);

  fmat dx(nX, nY), dy(nX, nY);
  fmat dlat_north(nX, nY);
  fmat mlt(nX, nY);
  fmat x_blend(nX, nY), y_blend(nX, nY);

  // Need to blend north and south, so use the distance from the pole
  // To indicate which pole you should use for MLT:
  
  // If you look throught the Earth from the north, the magnetic poles
  // in the north and south are on different sides of the Earth. If
  // you calculate MLT based on one pole, you get the wrong answer
  // near the other pole, so you should calculate the MLT for each
  // pole independently, and blend them near the equator.  This is
  // hard.  So, the way that I solved the problem is that the
  // locations of the poles (GSE X and Y, not Z) are blended, such
  // that when you are near the north pole, the location is near the
  // north (and south near south), and when you are close to the
  // equator, the location of the pole is the average of the two
  // pole locations.
  // So, x_blend and y_blend are the X,Y locations of the poles, blended
  // Then the dx and dy are calculated between each grid point and
  // the blended pole location. Then the angle between x and y is
  // calculated and converted to an hour.

  for (int iZ = 0; iZ < nZ; iZ++) {
    dlat_north = 1.0 - (pi/2.0-magLat_scgc.slice(iZ))/pi;
    x_blend = dlat_north * mag_pole_north_gse[0](0, 0, iZ) +
      (1.0 - dlat_north) * mag_pole_south_gse[0](0, 0, iZ);
    y_blend = dlat_north * mag_pole_north_gse[1](0, 0, iZ) +
      (1.0 - dlat_north) * mag_pole_south_gse[1](0, 0, iZ);
    dx = GSE_XYZ_vcgc[0].slice(iZ) - x_blend;
    dy = GSE_XYZ_vcgc[1].slice(iZ) - y_blend;
    
    mlt = (atan2(dy, dx) + 2*pi)/pi*12.0;
    magLocalTime_scgc.slice(iZ) = mlt - 24.0 * floor(mlt/24.0);
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
  int IsNorth = 1, IsSouth = 0;
  mag_pole_north_ll = get_magnetic_pole(IsNorth, planet, input, report);
  mag_pole_south_ll = get_magnetic_pole(IsSouth, planet, input, report);

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

  std::vector<fcube> lon_lat_radius;
  lon_lat_radius.push_back(geoLon_scgc);
  lon_lat_radius.push_back(geoLat_scgc);
  lon_lat_radius.push_back(radius_scgc);
  std::vector<fcube> xyz;
  
  xyz = transform_llr_to_xyz_3d(lon_lat_radius);
  geoX_scgc = xyz[0];
  geoY_scgc = xyz[0];
  geoZ_scgc = xyz[0];
}
