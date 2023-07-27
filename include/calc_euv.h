// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_CALC_EUV_H_
#define INCLUDE_CALC_EUV_H_

#include <vector>
#include <string>

#include "inputs.h"
#include "times.h"
#include "indices.h"
#include "planets.h"
#include "grid.h"
#include "euv.h"
#include "neutrals.h"
#include "ions.h"

// -------------------------------------------------------------------------
//
// -------------------------------------------------------------------------

int calc_euv(Planets planet,
             Grid grid,
             Times time,
             Euv &euv,
             Neutrals &neutrals,
             Ions &ions,
             Indices indices);

void calc_ionization_heating(Euv euv,
			     Neutrals &neutrals,
			     Ions &ions);


#endif  // INCLUDE_CALC_EUV_H_
