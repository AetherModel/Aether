// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_CALC_GRID_DERIVED_H_
#define INCLUDE_CALC_GRID_DERIVED_H_

#include <vector>

// ----------------------------------------------------------------------------
// 
// ----------------------------------------------------------------------------

std::vector<precision_t> calc_bin_edges(std::vector<precision_t> centers);
arma_vec calc_bin_edges(arma_vec centers);

std::vector<precision_t> calc_bin_widths(std::vector<precision_t> centers);
arma_vec calc_bin_widths(arma_vec centers);

#endif  // INCLUDE_CALC_GRID_DERIVED_H_
