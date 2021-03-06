// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <cmath>
#include <iostream>

#include "../include/aether.h"

bfield_info_type get_dipole(float lon,
                            float lat,
                            float alt,
                            Planets planet,
                            Inputs input,
                            Report &report) {

  std::string function = "dipole";
  static int iFunction = -1;
  report.enter(function, iFunction);

  bfield_info_type bfield_info;

  float llr[3];
  float xyz[3];
  float radius = planet.get_radius(lat);

  llr[0] = lon;
  llr[1] = lat;
  llr[2] = alt + radius;
  transform_llr_to_xyz(llr, xyz);

  float dipole_center[3];
  std::vector<float> temp_dipole_center = planet.get_dipole_center();
  transform_float_vector_to_array(temp_dipole_center, dipole_center);

  float delta_pos_to_center[3];
  vector_diff(xyz, dipole_center, delta_pos_to_center);

  float pos_rot_z[3];
  float magnetic_pole_rotation = planet.get_dipole_rotation();
  transform_rot_z(delta_pos_to_center, -magnetic_pole_rotation, pos_rot_z);

  float pos_rot_zy[3];
  float magnetic_pole_tilt = planet.get_dipole_tilt();
  transform_rot_y(pos_rot_z, -magnetic_pole_tilt, pos_rot_zy);

  float xypp2  = (pos_rot_zy[0]*pos_rot_zy[0] + pos_rot_zy[1]*pos_rot_zy[1]);

  float xyzpp2  = (pos_rot_zy[0]*pos_rot_zy[0] +
                   pos_rot_zy[1]*pos_rot_zy[1] +
                   pos_rot_zy[2]*pos_rot_zy[2]);

  float xypp = sqrt(xypp2);
  float xyzpp = sqrt(xyzpp2);

  float normal_r = (radius / xyzpp);
  float r3 = normal_r * normal_r * normal_r;

  // In GITM, L-Shell is defined with respect to the bottom of the
  // ionosphere.  This is so there can be a 0 deg magnetic latitude,
  // which can't really exist if L-shell is defined with respec to the
  // surface. But, to simplify things to begin with (and make it
  // planet agnostic), we use the classic definition of L-Shell, which
  // is with respect to the planetary radius.

  float cos_lat = xypp/xyzpp;
  float lShell = 1.0 / normal_r / (cos_lat * cos_lat);

  float mlat = acos(1.0 / sqrt(lShell));
  if (pos_rot_zy[2] < 0.0) mlat = -mlat;

  float mlon = acos(pos_rot_zy[0] / xypp);
  if (pos_rot_zy[1] < 0.0) mlon = -mlon;

  float dipole_strength = planet.get_dipole_strength();

  float b[3];
  b[0] = dipole_strength * r3 * 3 * pos_rot_zy[0] * pos_rot_zy[2] / xyzpp2;
  b[1] = dipole_strength * r3 * 3 * pos_rot_zy[2] * pos_rot_zy[1] / xyzpp2;
  b[2] = dipole_strength * r3 / xyzpp2 *
    (2 * pos_rot_zy[2]*pos_rot_zy[2] - xypp2);

  float b_rot_y[3];
  float b_rot_yz[3];
  float b_env[3];  // env = East, North, Vertical

  transform_rot_y(b, magnetic_pole_tilt, b_rot_y);
  transform_rot_z(b_rot_y, magnetic_pole_rotation, b_rot_yz);

  transform_vector_xyz_to_env(b_rot_yz, lon, lat, b_env);

  bfield_info.b[0] = b_env[0];
  bfield_info.b[1] = b_env[1];
  bfield_info.b[2] = b_env[2];
  bfield_info.lon = mlon;
  bfield_info.lat = mlat;

  report.exit(function);
  return bfield_info;
}
