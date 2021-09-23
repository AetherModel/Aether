// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "aether.h"

// -----------------------------------------------------------------------------
// main function to increment model states by one iteration. It needs
// so many inputs because it alters all of the states in the model.
// -----------------------------------------------------------------------------

int advance(Planets &planet,
            Grid &gGrid,
            Times &time,
            Euv &euv,
            Neutrals &neutrals,
            Ions &ions,
            Chemistry &chemistry,
            Indices &indices,
            Inputs &input,
            Report &report) {

  int iErr = 0;

  std::string function = "advance";
  static int iFunction = -1;
  report.enter(function, iFunction);

  if (time.check_time_gate(input.get_dt_report())) time.display();

  gGrid.calc_sza(planet, time, report);
  neutrals.calc_mass_density(report);
  neutrals.calc_specific_heat(report);
  time.calc_dt();

  iErr = calc_euv(planet,
                  gGrid,
                  time,
                  euv,
                  neutrals,
                  ions,
                  indices,
                  input,
                  report);

  neutrals.calc_conduction(gGrid, time, report);

  neutrals.add_sources(time, report);

  chemistry.calc_chemistry(neutrals, ions, time, gGrid, report);

  neutrals.set_bcs(report);
  neutrals.fill_with_hydrostatic(gGrid, report);

  time.increment_time();

  iErr = output(neutrals, ions, gGrid, time, planet, input, report);

  report.exit(function);
  return iErr;
}
