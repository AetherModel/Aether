// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>
#include <math.h>

#include "aether.h"

// ---------------------------------------------------------------------------
//  Fill in Solar Zenith Angle and cos(solar zenith angle)
// ---------------------------------------------------------------------------

void Grid::calc_sza(Planets planet, Times time) {

  std::string function = "Grid::calc_sza";
  static int iFunction = -1;
  report.enter(function, iFunction);

  precision_t lon_offset = planet.get_longitude_offset(time);
  precision_t sin_dec = planet.get_sin_dec(time);
  precision_t cos_dec = planet.get_cos_dec(time);

  // Local time is in radians
  geoLocalTime_scgc = geoLon_scgc + lon_offset;
  geoLocalTime_scgc =
    geoLocalTime_scgc - cTWOPI * floor(geoLocalTime_scgc / (cTWOPI));
  cos_sza_scgc =
    sin_dec * sin(geoLat_scgc) +
    cos_dec * cos(geoLat_scgc) % cos(geoLocalTime_scgc - cPI);
  sza_scgc = acos(cos_sza_scgc);

  report.exit(function);
}

// ---------------------------------------------------------------------------
//  Fill in GSE Coordinates
// ---------------------------------------------------------------------------

void Grid::calc_gse(Planets planet, Times time) {

  std::string function = "Grid::calc_sza";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // Compute GSE coordinates:
  // 1. use latitude / local time to derive XYZ
  std::vector<arma_cube> lon_lat_radius;
  // -pi is used here, since +X is pointed towards sun:
  lon_lat_radius.push_back(geoLocalTime_scgc - cPI);
  lon_lat_radius.push_back(geoLat_scgc);
  lon_lat_radius.push_back(radius_scgc);
  GSE_XYZ_vcgc = transform_llr_to_xyz_3d(lon_lat_radius);
  // 2. rotate by declination to point x-axis towards the sun
  precision_t declination = planet.get_declination(time);
  GSE_XYZ_vcgc = rotate_around_y_3d(GSE_XYZ_vcgc, -declination);

  // ---------------------------------------------------------------
  // Do the same thing for the magnetic poles:

  std::vector<arma_cube> lon_lat_radius_col;
  arma_cube tmp_col(1, 1, nZ);

  // North:
  // Lon (converted to local time):

  precision_t lon_offset = planet.get_longitude_offset(time);
  tmp_col.fill(mag_pole_north_ll[0] + lon_offset - cPI);
  lon_lat_radius_col.push_back(tmp_col);
  // Lat:
  tmp_col.fill(mag_pole_north_ll[1]);
  lon_lat_radius_col.push_back(tmp_col);
  // Radius:
  tmp_col.tube(0, 0) = radius_scgc.tube(nX / 2, nY / 2);
  lon_lat_radius_col.push_back(tmp_col);

  mag_pole_north_gse = transform_llr_to_xyz_3d(lon_lat_radius_col);
  mag_pole_north_gse = rotate_around_y_3d(mag_pole_north_gse, -declination);

  // South:
  // Lon:
  lon_lat_radius_col[0].fill(mag_pole_south_ll[0] + lon_offset - cPI);
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

void Grid::calc_mlt() {

  std::string function = "Grid::calc_mlt";
  static int iFunction = -1;
  report.enter(function, iFunction);

  arma_mat dx(nX, nY), dy(nX, nY);
  arma_mat dlat_north(nX, nY);
  arma_mat mlt(nX, nY);
  arma_mat x_blend(nX, nY), y_blend(nX, nY);

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
    dlat_north = 1.0 - (cPI / 2.0 - magLat_scgc.slice(iZ)) / cPI;
    x_blend = dlat_north * mag_pole_north_gse[0](0, 0, iZ) +
              (1.0 - dlat_north) * mag_pole_south_gse[0](0, 0, iZ);
    y_blend = dlat_north * mag_pole_north_gse[1](0, 0, iZ) +
              (1.0 - dlat_north) * mag_pole_south_gse[1](0, 0, iZ);
    dx = GSE_XYZ_vcgc[0].slice(iZ) - x_blend;
    dy = GSE_XYZ_vcgc[1].slice(iZ) - y_blend;

    mlt = (atan2(dy, dx) + cTWOPI) / cPI * 12.0 + 12.0;
    magLocalTime_scgc.slice(iZ) = mlt - 24.0 * floor(mlt / 24.0);
  }

  report.exit(function);
}

// -----------------------------------------------------------------------------
//  fill grid with magnetic field values
// -----------------------------------------------------------------------------

void Grid::fill_grid_bfield(Planets planet) {

  std::string function = "Grid::fill_grid_bfield";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t iLon, iLat, iAlt, iDim;
  precision_t lon, lat, alt;
  bfield_info_type bfield_info;
  bool DoDebug = false;

  bfield_mag_scgc.zeros();

  if (abs(planet.get_dipole_strength()) > 0) {
    HasBField = 1;

    for (iLon = 0; iLon < nLons; iLon++) {
      for (iLat = 0; iLat < nLats; iLat++) {
        for (iAlt = 0; iAlt < nAlts; iAlt++) {

          lon = geoLon_scgc(iLon, iLat, iAlt);
          lat = geoLat_scgc(iLon, iLat, iAlt);
          alt = geoAlt_scgc(iLon, iLat, iAlt);

          bfield_info = get_bfield(lon, lat, alt, DoDebug,
                                   planet);

          magLat_scgc(iLon, iLat, iAlt) = bfield_info.lat;
          magLon_scgc(iLon, iLat, iAlt) = bfield_info.lon;

          bfield_mag_scgc(iLon, iLat, iAlt) = 0.0;

          for (iDim = 0; iDim < 3; iDim++) {
            bfield_vcgc[iDim](iLon, iLat, iAlt) = bfield_info.b[iDim] * cNTtoT;
            bfield_mag_scgc(iLon, iLat, iAlt) =
              bfield_mag_scgc(iLon, iLat, iAlt) +
              bfield_vcgc[iDim](iLon, iLat, iAlt) * bfield_vcgc[iDim](iLon, iLat, iAlt);
          }

          bfield_mag_scgc(iLon, iLat, iAlt) =
            sqrt(bfield_mag_scgc(iLon, iLat, iAlt));
        }
      }
    }

    for (iDim = 0; iDim < 3; iDim++)
      bfield_unit_vcgc[iDim] = bfield_vcgc[iDim] / (bfield_mag_scgc + 1e-32);

    int IsNorth = 1, IsSouth = 0;
    mag_pole_north_ll = get_magnetic_pole(IsNorth, planet);
    mag_pole_south_ll = get_magnetic_pole(IsSouth, planet);
  }

  report.exit(function);
  return;
}

// -----------------------------------------------------------------------------
//  Fill in radius, radius^2, and 1/radius^2
// -----------------------------------------------------------------------------

void Grid::fill_grid_radius(Planets planet) {

  std::string function = "Grid::fill_grid_radius";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t iLon, iLat, iAlt;

  // This generalizes things so that radius could be a function of all
  // three dimensions.  The Cubesphere has different latitudes in the first
  // and second dimensions.
  for (iLon = 0; iLon < nLons; iLon++)
    for (iLat = 0; iLat < nLats; iLat++)
      for (iAlt = 0; iAlt < nAlts; iAlt++)
        radius_scgc(iLon, iLat, iAlt) =
          planet.get_radius(geoLat_scgc(iLon, iLat, iAlt));

  radius_scgc = radius_scgc + geoAlt_scgc;
  radius2_scgc = radius_scgc % radius_scgc;
  radius2i_scgc = 1.0 / radius2_scgc;

  report.exit(function);
  return;
}

// -----------------------------------------------------------------------------
// Calculates radial unit vector - it does this by taking the gradient
// of the radius and makes a unit vector out of this.  With a sphere,
// this should be in the 3rd dimension, with an oblate spheriod, there
// will be a latitudinal component.
// -----------------------------------------------------------------------------

void Grid::calc_rad_unit(Planets planet) {

  std::string function = "Grid::calc_rad_unit";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // *this is the grid class....
  std::vector<arma_cube> gradient_vcgc = calc_gradient_vector(radius_scgc, *this);

  arma_cube mag_radius_gradient = sqrt( pow(gradient_vcgc[0], 2) +
                                        pow(gradient_vcgc[1], 2) +
                                        pow(gradient_vcgc[2], 2));
  arma_cube mag_radius_gradienti =  1.0 / mag_radius_gradient;

  rad_unit_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);

  for (int iV = 0; iV < 3; iV++)
    rad_unit_vcgc[iV].zeros();

  rad_unit_vcgc[0] = (gradient_vcgc[0] % mag_radius_gradienti);
  rad_unit_vcgc[1] = (gradient_vcgc[1] % mag_radius_gradienti);
  rad_unit_vcgc[2] = (gradient_vcgc[2] % mag_radius_gradienti);

  report.exit(function);
  return;
}

// -----------------------------------------------------------------------------
//  Calculates gravity, including J2 perturbation
// -----------------------------------------------------------------------------

void Grid::calc_gravity(Planets planet) {

  std::string function = "Grid::calc_gravity";
  static int iFunction = -1;
  report.enter(function, iFunction);

  precision_t mu = planet.get_mu();

  gravity_potential_scgc =
    - (mu / radius_scgc)
    + ((3 * (planet.get_J2() * mu)) /
       (2 * pow(radius_scgc, 3)) %
       ((sin(geoLat_scgc) % sin(geoLat_scgc)) - 1.0));

  // *this is the grid class....
  gravity_vcgc = calc_gradient_vector(gravity_potential_scgc, *this);

  gravity_vcgc[0] = - gravity_vcgc[0];
  gravity_vcgc[1] = - gravity_vcgc[1];
  gravity_vcgc[2] = - gravity_vcgc[2];
  gravity_mag_scgc = sqrt(
    gravity_vcgc[0] % gravity_vcgc[0] + 
    gravity_vcgc[1] % gravity_vcgc[1] + 
    gravity_vcgc[2] % gravity_vcgc[2]);

  report.exit(function);
  return;
}

// -----------------------------------------------------------------------------
//  Fill in XYZ in geo and mag coordinates
// -----------------------------------------------------------------------------

void Grid::calc_grid_spacing(Planets planet) {

  int64_t iLon, iLat, iAlt;

  report.print(3, "starting calc_grid_spacing");

  calc_alt_grid_spacing();
  calc_lat_grid_spacing();
  calc_long_grid_spacing();

  std::vector<arma_cube> lon_lat_radius;
  lon_lat_radius.push_back(geoLon_scgc);
  lon_lat_radius.push_back(geoLat_scgc);
  lon_lat_radius.push_back(radius_scgc);
  std::vector<arma_cube> xyz;

  xyz = transform_llr_to_xyz_3d(lon_lat_radius);
  geoX_scgc = xyz[0];
  geoY_scgc = xyz[0];
  geoZ_scgc = xyz[0];

  report.print(3, "ending calc_grid_spacing");
}

// ---------------------------------------
// Grid spacing for altitude:
// ---------------------------------------

void Grid::calc_alt_grid_spacing() {

  int64_t iAlt;

  for (iAlt = 1; iAlt < nAlts - 1; iAlt++) {
    dalt_center_scgc.slice(iAlt) =
      (geoAlt_scgc.slice(iAlt + 1) - geoAlt_scgc.slice(iAlt - 1)) / 2.0;
    dalt_lower_scgc.slice(iAlt) =
      geoAlt_scgc.slice(iAlt) - geoAlt_scgc.slice(iAlt - 1);
  }

  dalt_center_scgc.slice(0) = dalt_center_scgc.slice(1);
  dalt_center_scgc.slice(nAlts - 1) = dalt_center_scgc.slice(nAlts - 2);

  dalt_lower_scgc.slice(0) = dalt_lower_scgc.slice(1);
  iAlt = nAlts - 1;
  dalt_lower_scgc.slice(iAlt) =
    geoAlt_scgc.slice(iAlt) - geoAlt_scgc.slice(iAlt - 1);

  // For a stretched grid, calculate some useful quantities:
  // lower is defined for the current cell, which
  // means that upper(iAlt) is lower(iAlt+1)
  // ratio = upper / lower
  for (iAlt = 0; iAlt < nAlts - 1; iAlt++)
    dalt_ratio_scgc.slice(iAlt) =
      dalt_lower_scgc.slice(iAlt + 1) / dalt_lower_scgc.slice(iAlt);

  iAlt = nAlts - 1;
  dalt_ratio_scgc.slice(iAlt) = dalt_ratio_scgc.slice(iAlt - 1);

  // Need the square of the ratio:
  dalt_ratio_sq_scgc = dalt_ratio_scgc % dalt_ratio_scgc;
}

// ---------------------------------------
// Grid spacing for latitude:
// ---------------------------------------

void Grid::calc_lat_grid_spacing() {

  int64_t iLat;

  for (iLat = 1; iLat < nLats - 1; iLat++) {
    dlat_center_scgc.col(iLat) =
      (geoLat_scgc.col(iLat + 1) - geoLat_scgc.col(iLat - 1)) / 2.0;
  }

  // Bottom (one sided):
  iLat = 0;
  dlat_center_scgc.col(iLat) =
    geoLat_scgc.col(iLat + 1) - geoLat_scgc.col(iLat);
  // Top (one sided):
  iLat = nLats - 1;
  dlat_center_scgc.col(iLat) =
    geoLat_scgc.col(iLat) - geoLat_scgc.col(iLat - 1);

  // Make this into a distance:
  dlat_center_dist_scgc = dlat_center_scgc % radius_scgc;
}

// ---------------------------------------
// Grid spacing for longitude:
// ---------------------------------------

void Grid::calc_long_grid_spacing() {

  int64_t iLon;

  for (iLon = 1; iLon < nLons - 1; iLon++)
    dlon_center_scgc.row(iLon) =
      (geoLon_scgc.row(iLon + 1) - geoLon_scgc.row(iLon - 1)) / 2.0;

  // Bottom (one sided):
  iLon = 0;
  dlon_center_scgc.row(iLon) =
    geoLon_scgc.row(iLon + 1) - geoLon_scgc.row(iLon);
  // Top (one sided):
  iLon = nLons - 1;
  dlon_center_scgc.row(iLon) =
    geoLon_scgc.row(iLon) - geoLon_scgc.row(iLon - 1);

  // Make this into a distance:
  dlon_center_dist_scgc =
    dlon_center_scgc % radius_scgc % abs(cos(geoLat_scgc));
}

// -----------------------------------------------------------------------------
//  Calaculate centripetal acceleration
// -----------------------------------------------------------------------------

void Grid::calc_cent_acc(Planets planet) {
  // Longitudanal cent acc
  cent_acc_vcgc[0].zeros();

  float omega = planet.get_omega();
  float omega2 = omega * omega;

  // Latitudinal cent acc
  cent_acc_vcgc[1] = -1 * omega2 * radius_scgc  % cos(geoLat_scgc) % sin(
                       geoLat_scgc);

  // Radial cent acc
  cent_acc_vcgc[2] = omega2 * radius_scgc % cos(geoLat_scgc) % cos(geoLat_scgc);
}
