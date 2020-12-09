// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef AETHER_INCLUDE_BFIELD_H_
#define AETHER_INCLUDE_BFIELD_H_

struct bfield_info_type {
  float b[3];
  float lon;
  float lat;
};

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

bfield_info_type get_dipole(float lon,
			    float lat,
			    float alt,
			    Planets planet,
			    Inputs input,
			    Report &report);

#endif // AETHER_INCLUDE_EUV_H_

