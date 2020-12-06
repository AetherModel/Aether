// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef AETHER_INCLUDE_INIT_FILL_GRID_H_
#define AETHER_INCLUDE_INIT_FILL_GRID_H_

#include "inputs.h"
#include "planets.h"
#include "grid.h"
#include "times.h"

void fill_grid_radius(Grid &gGrid, Planets planet, Inputs input);
// void fill_grid_sza(Grid &gGrid, Planets planet, Times time, Inputs input);
// void fill_grid(Grid &gGrid, Planets planet, Inputs input);

#endif // AETHER_INCLUDE_INIT_FILL_GRID_H_
