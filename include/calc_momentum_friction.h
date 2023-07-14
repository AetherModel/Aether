// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_CALC_MOMENTUM_FRICTION_H_
#define INCLUDE_CALC_MOMENTUM_FRICTION_H_

arma_vec neutral_friction_one_cell(int64_t iLong, int64_t iLat, int64_t iAlt,
				   arma_vec &vels,
				   Neutrals &neutrals,
                                   Report &report);

void calc_neutral_friction(Neutrals &neutrals,
                           Report &report);


/**********************************************************************
  \brief Calculate acceleration due to ion drag
  \param ions The ions with which we are calculating drag
  \param grid The grid to define the neutrals on
  \param dt The change in time
  \param report allow reporting to occur
**/
void calc_ion_collisions(Neutrals &neutrals,
			 Ions &ions,
			 Report &report);


#endif  // INCLUDE_CALC_MOMENTUM_FRICTION_H_
