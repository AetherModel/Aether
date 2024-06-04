// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <cmath>
#include <iostream>

#include "aether.h"

// -----------------------------------------------------------------------------
// get the l-shell given latitude (in radians) and normalized radius
// -----------------------------------------------------------------------------

precision_t get_lshell(precision_t lat, precision_t rNorm) {
  precision_t cosLat = cos(lat);
  precision_t lshell = rNorm / (cosLat * cosLat);
  return lshell;
}

arma_vec get_lat_from_r_and_lshell(arma_vec r, precision_t lshell) {
  arma_vec cosLat = sqrt(r / lshell);
  cosLat.clamp(-1.0, 1.0);
  arma_vec lat = acos(cosLat);
  return lat;
}

precision_t get_lat_from_r_and_lshell(precision_t r, precision_t lshell) {
  precision_t cosLat = sqrt(r / lshell);
  if (cosLat < -1.0) cosLat = -1.0;
  if (cosLat > 1.0) cosLat = 1.0;
  precision_t lat = acos(cosLat);
  return lat;
}



// -----------------------------------------------------------------------------
// Calculate a tilted offset dipole field given the planetary
// characteristics
// -----------------------------------------------------------------------------

bfield_info_type get_dipole(precision_t lon,
                            precision_t lat,
                            precision_t alt,
                            bool DoDebug,
                            Planets planet) {

  std::string function = "dipole";
  static int iFunction = -1;

  if (DoDebug)
    report.enter(function, iFunction);

  bfield_info_type bfield_info;

  precision_t llr[3];
  precision_t xyz[3];
  precision_t radius = planet.get_radius(lat);

  llr[0] = lon;
  llr[1] = lat;
  llr[2] = alt + radius;
  transform_llr_to_xyz(llr, xyz);

  precision_t dipole_center[3];
  std::vector<float> temp_dipole_center = planet.get_dipole_center();
  transform_float_vector_to_array(temp_dipole_center, dipole_center);

  precision_t delta_pos_to_center[3];
  vector_diff(xyz, dipole_center, delta_pos_to_center);

  precision_t pos_rot_z[3];
  precision_t magnetic_pole_rotation = planet.get_dipole_rotation();
  transform_rot_z(delta_pos_to_center, -magnetic_pole_rotation, pos_rot_z);

  precision_t pos_rot_zy[3];
  precision_t magnetic_pole_tilt = planet.get_dipole_tilt();
  transform_rot_y(pos_rot_z, -magnetic_pole_tilt, pos_rot_zy);

  precision_t xypp2  = (pos_rot_zy[0] * pos_rot_zy[0] +
                        pos_rot_zy[1] * pos_rot_zy[1]);

  precision_t xyzpp2  = (pos_rot_zy[0] * pos_rot_zy[0] +
                         pos_rot_zy[1] * pos_rot_zy[1] +
                         pos_rot_zy[2] * pos_rot_zy[2]);

  precision_t xypp = sqrt(xypp2);
  precision_t xyzpp = sqrt(xyzpp2);

  precision_t normal_r = (radius / xyzpp);
  precision_t r3 = normal_r * normal_r * normal_r;

  // L-Shell is the equatorial distance to the magnetic field-line
  // normalized to some specific distance, typically the planet
  // radius.  In older models, this distance has been the planet
  // radius + distance to the bottom of the model.  This is so there
  // can be a 0 deg magnetic latitude, which can't really exist if
  // L-shell is defined with respect to the surface. But, to simplify
  // things to begin with (and make it planet agnostic), we use the
  // classic definition of L-Shell, which is with respect to the
  // planetary radius.

  precision_t cos_lat = xypp / xyzpp;
  precision_t lShell = 1.0 / normal_r / (cos_lat * cos_lat);

  precision_t mlat = acos(1.0 / sqrt(lShell));

  if (pos_rot_zy[2] < 0.0)
    mlat = -mlat;

  precision_t mlon = acos(pos_rot_zy[0] / xypp);

  if (pos_rot_zy[1] < 0.0)
    mlon = -mlon;

  precision_t dipole_strength = planet.get_dipole_strength();

  precision_t b[3];
  b[0] = dipole_strength * r3 * 3 * pos_rot_zy[0] * pos_rot_zy[2] / xyzpp2;
  b[1] = dipole_strength * r3 * 3 * pos_rot_zy[2] * pos_rot_zy[1] / xyzpp2;
  b[2] = dipole_strength * r3 / xyzpp2 *
         (2 * pos_rot_zy[2] * pos_rot_zy[2] - xypp2);

  precision_t b_rot_y[3];
  precision_t b_rot_yz[3];
  precision_t b_env[3];  // env = East, North, Vertical

  transform_rot_y(b, magnetic_pole_tilt, b_rot_y);
  transform_rot_z(b_rot_y, magnetic_pole_rotation, b_rot_yz);

  transform_vector_xyz_to_env(b_rot_yz, lon, lat, b_env);

  bfield_info.b[0] = b_env[0];
  bfield_info.b[1] = b_env[1];
  bfield_info.b[2] = b_env[2];
  bfield_info.lon = mlon;
  bfield_info.lat = mlat;

  if (DoDebug)
    report.exit(function);

  return bfield_info;
}
