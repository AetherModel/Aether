// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_CALC_MOMENTUM_FRICTION_H_
#define INCLUDE_CALC_MOMENTUM_FRICTION_H_

arma_vec neutral_friction_one_cell(int64_t iLong, int64_t iLat, int64_t iAlt,
				   precision_t dt,
				   arma_vec &vels,
				   Neutrals &neutrals,
                                   Report &report);

void calc_neutral_friction(precision_t dt,
                           Grid &gGrid,
                           Neutrals &neutrals,
                           Report &report);

#endif  // INCLUDE_CALC_MOMENTUM_FRICTION_H_
