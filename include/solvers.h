// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_SOLVERS_H_
#define INCLUDE_SOLVERS_H_

// The armadillo library is to allow the use of 3d cubes and other
// array types, with array math built in. This eliminates loops!
#include <armadillo>
using namespace arma;

fvec solver_conduction(fvec value,
		       fvec lambda,
		       fvec front,
		       float dt,
		       fvec dx);


fcube solver_chemistry(fcube density,
		       fcube source,
		       fcube loss,
		       float dt);

#endif  // INCLUDE_SOLVERS_H_
