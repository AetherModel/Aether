// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_BFIELD_H_
#define INCLUDE_BFIELD_H_

struct bfield_info_type {
  precision_t b[3];
  precision_t lon;
  precision_t lat;
};

precision_t get_lshell(precision_t lat, precision_t rNorm);
arma_vec get_lat_from_r_and_lshell(arma_vec r, precision_t lshell);
precision_t get_lat_from_r_and_lshell(precision_t r, precision_t lshell);

arma_vec get_magnetic_pole(int IsNorth,
		       Planets planet);

bfield_info_type get_bfield(precision_t lon,
                            precision_t lat,
                            precision_t alt,
			    bool DoDebug,
                            Planets planet);

bfield_info_type get_dipole(precision_t lon,
                            precision_t lat,
                            precision_t alt,
			    bool DoDebug,
                            Planets planet);

#endif // INCLUDE_BFIELD_H_
