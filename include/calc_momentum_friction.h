// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_CALC_MOMENTUM_FRICTION_H_
#define INCLUDE_CALC_MOMENTUM_FRICTION_H_

#include "../include/aether.h"


/**********************************************************************
  \brief Calculate acceleration due to ion drag
  \param ions The ions with which we are calculating drag
  \param grid The grid to define the neutrals on
  \param dt The change in time
  \param report allow reporting to occur
**/
void calc_ion_collisions(Neutrals &neutrals, Ions &ions);

#endif  // INCLUDE_CALC_MOMENTUM_FRICTION_H_
