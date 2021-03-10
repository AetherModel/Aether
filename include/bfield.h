// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_BFIELD_H_
#define INCLUDE_BFIELD_H_

struct bfield_info_type {
  float b[3];
  float lon;
  float lat;
};

fvec get_magnetic_pole(int IsNorth,
		       Planets planet,
		       Inputs input,
		       Report &report);

bfield_info_type get_bfield(float lon,
                            float lat,
                            float alt,
                            Planets planet,
                            Inputs input,
                            Report &report);

bfield_info_type get_dipole(float lon,
                            float lat,
                            float alt,
                            Planets planet,
                            Inputs input,
                            Report &report);

#endif // INCLUDE_BFIELD_H_
