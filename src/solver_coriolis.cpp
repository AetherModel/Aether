// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

// Added by S. Annamalai - June 12, 2023

// This function calculates the coriolis acceleration vector for an
// object in motion in a rotating frame of reference. It uses the
// coriolis terms in equations (14), (18), (19) in
// https://drive.google.com/file/d/1Q0cMzKhdd0IoXzBl3odTge9JhBwfAyk9

#include "../include/aether.h"

std::vector<arma_cube> coriolis(std::vector<arma_cube> velocity,
				precision_t rotation_rate,
				arma_cube lat_scgc) {
  std::vector<arma_cube> coriolis_vec(3);
  coriolis_vec[0] = -2 * rotation_rate * velocity[1] % sin(lat_scgc);
  coriolis_vec[1] =
    2 * rotation_rate * velocity[0] % sin(lat_scgc) -
    2 * rotation_rate * velocity[2] % cos(lat_scgc);
  coriolis_vec[2] = 2 * rotation_rate * cos(lat_scgc) % velocity[1];
  return coriolis_vec;
}
