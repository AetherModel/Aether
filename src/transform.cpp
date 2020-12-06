// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <math.h>

// -----------------------------------------------------------------------
// Transform Longitude, Latitude, Radius to X, Y, Z
// -----------------------------------------------------------------------

void transform_llr_to_xyz(float llr[3], float xyz[3]) {

  // llr[0] = longitude (in radians)
  // llr[1] = latitude (in radians)
  // llr[2] = radius
  
  xyz[0] = llr[2] * cos(llr[1]) * cos(llr[0]);
  xyz[1] = llr[2] * cos(llr[1]) * sin(llr[0]);
  xyz[2] = llr[2] * sin(llr[1]);
  
}


