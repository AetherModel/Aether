// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef AETHER_INCLUDE_SOLVERS_H_
#define AETHER_INCLUDE_SOLVERS_H_

#include <armadillo>
using namespace arma;
fvec solver_conduction_old(fvec value,
		       fvec lambda,
		       fvec front,
		       float dt,
		       fvec dx);

fvec solver_conduction_new(fvec value,
		       fvec lambda,
		       fvec front,
		       float dt,
		       fvec dx);

fvec solver_conduction(fvec value,
		       fvec lambda,
		       fvec front,
		       float dt,
		       fvec dx);


float solver_chemistry(float old_density,
		       float source,
		       float loss,
		       float dt);

#endif // AETHER_INCLUDE_SOLVERS_H_
