// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_ADVANCE_H_
#define INCLUDE_ADVANCE_H_

#include "../include/electrodynamics.h"
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
