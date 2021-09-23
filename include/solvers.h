// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_SOLVERS_H_
#define INCLUDE_SOLVERS_H_

// The armadillo library is to allow the use of 3d cubes and other
// array types, with array math built in. This eliminates loops!
#include <armadillo>
using namespace arma;

arma_vec solver_conduction(arma_vec value,
		       arma_vec lambda,
		       arma_vec front,
		       precision_t dt,
		       arma_vec dx);


arma_cube solver_chemistry(arma_cube density,
		       arma_cube source,
		       arma_cube loss,
		       precision_t dt);

#endif  // INCLUDE_SOLVERS_H_
