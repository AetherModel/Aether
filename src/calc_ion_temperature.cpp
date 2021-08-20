// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"

// --------------------------------------------------------------------------
// Hack to just set the ion temperature to the neutral temperature
// --------------------------------------------------------------------------

void Ions::calc_ion_temperature(Neutrals neutrals, Grid grid, Report &report) {

  int64_t iIon;

  for (iIon = 0; iIon < nIons; iIon++)
    species[iIon].temperature_scgc = neutrals.temperature_scgc;

  ion_temperature_scgc = neutrals.temperature_scgc;
  electron_temperature_scgc = neutrals.temperature_scgc;
  return;
}

