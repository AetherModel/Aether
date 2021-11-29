// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>
#include <string>

#include "../include/aether.h"

// -----------------------------------------------------------------------------
// Adds all of the sources to the states. Needs time to get dt.
// -----------------------------------------------------------------------------

void Neutrals::add_sources(Times time, Report &report) {

  std::string function = "add_sources";
  static int iFunction = -1;
  report.enter(function, iFunction);

  precision_t dt = time.get_dt();

  temperature_scgc =
    temperature_scgc +
    dt * (heating_euv_scgc +
          conduction_scgc);

  report.exit(function);
  return;
}
