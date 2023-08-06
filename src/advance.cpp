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
            Electrodynamics &electrodynamics,
            Indices &indices,
            Logfile &logfile) {

  int iErr = 0;

  std::string function = "advance";
  static int iFunction = -1;
  report.enter(function, iFunction);

  if (time.check_time_gate(input.get_dt_report()) &&
      report.test_verbose(0))
    time.display();

  if (input.get_is_student())
    report.print(-1, "(1) What function is this " +
                 input.get_student_name() + "?");

  gGrid.calc_sza(planet, time);
  neutrals.calc_mass_density();
  neutrals.calc_specific_heat();
  neutrals.calc_concentration();
  neutrals.calc_mean_major_mass();
  neutrals.calc_pressure();
  neutrals.calc_bulk_velocity();
  neutrals.calc_kappa_eddy();
  time.calc_dt();

  iErr = calc_euv(planet,
                  gGrid,
                  time,
                  euv,
                  neutrals,
                  ions,
                  indices);

  iErr = electrodynamics.update(planet,
                                gGrid,
                                time,
                                ions);
  calc_ion_neutral_coll_freq(neutrals, ions);
  ions.calc_ion_drift(neutrals, gGrid, time.get_dt());

  calc_aurora(gGrid, neutrals, ions);

  // Calculate some neutral source terms:
  neutrals.calc_conduction(gGrid, time);
  chemistry.calc_chemistry(neutrals, ions, time, gGrid);
  if (input.get_O_cooling())
    neutrals.calc_O_cool();
  if (input.get_NO_cooling())
    neutrals.calc_NO_cool();
  neutrals.add_sources(time);

  neutrals.calc_conduction(gGrid, time);
  chemistry.calc_chemistry(neutrals, ions, time, gGrid);

  neutrals.vertical_momentum_eddy(gGrid);
  calc_ion_collisions(neutrals, ions);
  calc_neutral_friction(neutrals);

  neutrals.add_sources(time);
  ions.calc_ion_temperature(neutrals, gGrid, time);
  ions.calc_electron_temperature(neutrals, gGrid);

  neutrals.set_bcs(gGrid, time, indices);
  neutrals.calc_scale_height(gGrid);
  neutrals.fill_with_hydrostatic(gGrid);
  neutrals.exchange(gGrid);
  
  time.increment_time();

  if (time.check_time_gate(input.get_dt_write_restarts())) {
    report.print(3, "Writing restart files");
    neutrals.restart_file(input.get_restartout_dir(), DoWrite);
    ions.restart_file(input.get_restartout_dir(), DoWrite);
    time.restart_file(input.get_restartout_dir(), DoWrite);
  }

  if(input.get_nan_test()){
      neutrals.nan_test(input.get_nan_test_variable());
      ions.nan_test(input.get_nan_test_variable());
    }

  if(input.get_check_for_nans()){
    neutrals.check_for_nonfinites(report);
    ions.check_for_nonfinites(report);
  }

  iErr = output(neutrals, ions, gGrid, time, planet);

  logfile.write_logfile(indices, neutrals, ions, gGrid, time);

  report.exit(function);
  return iErr;
}
