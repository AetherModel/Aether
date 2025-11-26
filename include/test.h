// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_TEST_H_
#define INCLUDE_TEST_H_

#include "aether.h"


// Gradient tests
// Cubesphere is not done nor tested
bool test_gradient(Planets planet, Quadtree quadtree, json test_config,
                   Grid gGrid, Grid mGrid);
bool test_gradient_cubesphere(Planets planet, Quadtree quadtree, Grid grid);
bool test_gradient_ijk(Planets planet, Grid grid, bool debug);


#endif