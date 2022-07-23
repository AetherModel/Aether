// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_SPHERE_H_
#define INCLUDE_SPHERE_H_

#include "aether.h"
#include <memory>

/*************************************************
 * \brief A namespace with all quad sphere grid logic.
 *************************************************/
namespace Sphere {

  /// The normalized origins of each face of the cube (i.e. corner)
  static const arma_mat ORIGINS = {
				   { 0.0, -0.5, 0.0}
  };

  /// Normalized right steps in cube
  static const arma_mat RIGHTS = {
				  {2.0, 0.0, 0.0}
  };

  /// Normalized right steps in cube
  static const arma_mat UPS = {
			       {0.0, 1.0, 0.0}
  };

};

#endif  // INCLUDE_SPHERE_H_
