// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_CALC_GRID_DERIVED_H_
#define INCLUDE_CALC_GRID_DERIVED_H_

#include <vector>

// ----------------------------------------------------------------------------
// 
// ----------------------------------------------------------------------------

std::vector<float> calc_bin_edges(std::vector<float> centers);
fvec calc_bin_edges(fvec centers);

std::vector<float> calc_bin_widths(std::vector<float> centers);
fvec calc_bin_widths(fvec centers);

#endif  // INCLUDE_CALC_GRID_DERIVED_H_
