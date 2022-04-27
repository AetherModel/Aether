// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"

// --------------------------------------------------------------------------
// Hack to just set the electon temperature to the neutral temperature
// --------------------------------------------------------------------------

void Ions::calc_electron_temperature(Neutrals neutrals, Grid grid,
                                     Report &report) {

  electron_temperature_scgc = neutrals.temperature_scgc;

}


