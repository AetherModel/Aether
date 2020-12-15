// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef AETHER_INCLUDE_SOLVERS_H_
#define AETHER_INCLUDE_SOLVERS_H_

#include "sizes.h"

int solver_conduction(float value[nGeoAltsG],
		      float lambda[nGeoAltsG],
		      float front[nGeoAltsG],
		      float dt,
		      float dalt_lower[nGeoAltsG],
		      float *conduction);

float solver_chemistry(float old_density,
		       float source,
		       float loss,
		       float dt);

#endif // AETHER_INCLUDE_SOLVERS_H_
