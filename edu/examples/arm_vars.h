// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef AETHER_EDU_ARM_VARS_H_
#define AETHER_EDU_ARM_VARS_H_

#include <armadillo>

using namespace arma;

class Grid {

 public:

  Grid(int nX, int nY, int nZ);

  void set_radius(float planet_radius, fcube alts3d);
  fcube get_radius();
  
 private:

  fcube radius3d;
  std::vector<fcube> gravity3dv;
  
};
  
#endif // AETHER_EDU_ARM_VARS_H_

