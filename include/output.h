// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_OUTPUT_H_
#define INCLUDE_OUTPUT_H_

#include "../include/neutrals.h"
#include "../include/grid.h"
#include "../include/times.h"
#include "../include/planets.h"
#include "../include/inputs.h"
#include "../include/report.h"

int output(Neutrals neutrals,
           Ions ions,
           Grid grid,
           Times time,
           Planets planet,
           Inputs args,
           Report &report);

#endif  // INCLUDE_OUTPUT_H_
