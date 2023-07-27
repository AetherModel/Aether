// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_BFIELD_H_
#define INCLUDE_BFIELD_H_

struct bfield_info_type {
  precision_t b[3];
  precision_t lon;
  precision_t lat;
};

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
