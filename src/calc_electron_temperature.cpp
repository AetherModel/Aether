// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"

// --------------------------------------------------------------------------
// TODO (#24): this currently just sets the electron temperature to the neutral temperature
// --------------------------------------------------------------------------

void Ions::calc_electron_temperature(Neutrals neutrals, Grid grid) {

  std::string function = "Ions::calc_electron_temperature";
  static int iFunction = -1;
  report.enter(function, iFunction);

  electron_temperature_scgc = neutrals.temperature_scgc;

  report.exit(function);
}

// --------------------------------------------------------------------------
// Heating terms:
//  - [ ] photoelectrons
//  - [ ] auroral ionization
//  - [ ] e- ion collisions
//  - [ ] e- neutral collisions (elastic & inelastic)
//  - [ ] e- chemistry (O2, V2 vibration; O fine structure, O exitation)
// --------------------------------------------------------------------------
