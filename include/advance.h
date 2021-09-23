// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_ADVANCE_H_
#define INCLUDE_ADVANCE_H_

/**************************************************************
 * 
 * advance:
 *
 * - Function that advances the states in Aether by one time step
 *
 *   Pretty much all classes have to be passed into this function
 *   this function calls a bunch of functions that alters all
 *   of the states in the system.
 **************************************************************/

#include "../include/times.h"
#include "../include/inputs.h"
#include "../include/neutrals.h"
#include "../include/euv.h"
#include "../include/grid.h"
#include "../include/planets.h"
#include "../include/ions.h"


int advance(Planets &planet,
            Grid &gGrid,
            Times &time,
            Euv &euv,
            Neutrals &neutrals,
            Ions &ions,
            Chemistry &chemistry,
            Electrodynamics &electrodynamics,
            Indices &indices,
            Inputs &args,
            Report &report);

#endif // INCLUDE_ADVANCE_H_
