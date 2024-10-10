// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "aether.h"

// -----------------------------------------------------------------------------
// main function to increment model states by one iteration. It needs
// so many inputs because it alters all of the states in the model.
// -----------------------------------------------------------------------------

bool advance(Planets &planet,
             Grid &gGrid,
             Grid &mGrid,
             Times &time,
             Euv &euv,
             Neutrals &neutrals,
             Neutrals &neutralsMag,
             Ions &ions,
             Ions &ionsMag,
             Chemistry &chemistry,
             Chemistry &chemistryMag,
             Electrodynamics &electrodynamics,
             Indices &indices,
             Logfile &logfile,
             Logfile &logfileMag) {

  bool didWork = true;

  std::string function = "advance";
  static int iFunction = -1;
  report.enter(function, iFunction);

  if (time.check_time_gate(input.get_dt_report()) &&
      report.test_verbose(0))
    time.display();

  if (input.get_is_student())
    report.print(-1, "(1) What function is this " +
                 input.get_student_name() + "?");

  if (didWork & input.get_check_for_nans()) {
    didWork = neutrals.check_for_nonfinites("Top of Advance - neu grid");
    didWork = neutralsMag.check_for_nonfinites("Top of Advance - ion grid");
  }

  gGrid.calc_sza(planet, time);
  mGrid.calc_sza(planet, time);

  neutrals.calc_mass_density();
  neutrals.calc_mean_major_mass();
  neutrals.calc_specific_heat();
  neutrals.calc_concentration();
  neutrals.calc_pressure();
  neutrals.calc_bulk_velocity();
  neutrals.calc_kappa_eddy();
  neutrals.calc_viscosity();
  //neutrals.calc_cMax();

  ions.fill_electrons();
  ions.calc_sound_speed();
  ions.calc_cMax();

  precision_t dtNeutral = calc_dt(gGrid, neutrals.cMax_vcgc);
  precision_t dtIon = calc_dt(gGrid, ions.cMax_vcgc);
  time.calc_dt(dtNeutral, dtIon);
  
  neutralsMag.calc_mass_density();
  neutralsMag.calc_mean_major_mass();
  neutralsMag.calc_specific_heat();
  neutralsMag.calc_concentration();
  neutralsMag.calc_pressure();
  neutralsMag.calc_bulk_velocity();
  neutralsMag.calc_kappa_eddy();
  neutralsMag.calc_cMax();

  // ------------------------------------
  // Do advection first :

  // Upper BCs requires the scale height to be calculated, so do that
  // first

  neutrals.calc_scale_height(gGrid);
  neutralsMag.calc_scale_height(mGrid);

  if (didWork)
    didWork = neutrals.set_bcs(gGrid, time, indices); 
  if (didWork)
    didWork = ions.set_bcs(gGrid, time, indices);
  if (didWork)
    didWork = neutralsMag.set_bcs(mGrid, time, indices);

  if (gGrid.get_nAlts(false) > 1)
    neutrals.advect_vertical(gGrid, time);

  neutrals.exchange_old(gGrid);
  //advect(gGrid, time, neutrals);

  if (didWork & input.get_check_for_nans())
    didWork = neutrals.check_for_nonfinites("After Horizontal Advection");
  
  // ------------------------------------
  // Calculate source terms next:

  if (didWork)
    didWork = calc_euv(planet,
                       gGrid,
                       time,
                       euv,
                       neutrals,
                       ions,
                       indices);
  if (didWork)
    didWork = calc_euv(planet,
                       mGrid,
                       time,
                       euv,
                       neutralsMag,
                       ionsMag,
                       indices);

  if (didWork)
    didWork = electrodynamics.update(planet,
                                     gGrid,
                                     time,
                                     indices,
                                     ions);
  if (didWork)
    didWork = electrodynamics.update(planet,
                                     mGrid,
                                     time,
                                     indices,
                                     ionsMag);


  if (didWork) {
    calc_ion_neutral_coll_freq(neutrals, ions);
    ions.calc_ion_drift(neutrals, gGrid, time.get_dt());

    calc_aurora(gGrid, neutrals, ions);
    calc_aurora(mGrid, neutralsMag, ionsMag);

    // Calculate chemistry on both grids:
    chemistry.calc_chemistry(neutrals, ions, time, gGrid);
    chemistryMag.calc_chemistry(neutralsMag, ionsMag, time, mGrid);

    if (input.get_O_cooling())
      neutrals.calc_O_cool();

    if (input.get_NO_cooling())
      neutrals.calc_NO_cool();

    calc_ion_collisions(neutrals, ions);

    neutrals.add_sources(time, planet, gGrid);
    if (didWork & input.get_check_for_nans())
      didWork = neutrals.check_for_nonfinites("After Add Sources");
    neutralsMag.add_sources(time, planet, mGrid);

    ions.calc_ion_temperature(neutrals, gGrid, time);
    ions.calc_electron_temperature(neutrals, gGrid);
    ionsMag.calc_ion_temperature(neutralsMag, mGrid, time);
    ionsMag.calc_electron_temperature(neutralsMag, mGrid);

    if (didWork & input.get_check_for_nans())
      didWork = neutrals.check_for_nonfinites("After Vertical Advection");

    neutrals.exchange_old(gGrid);

    time.increment_time();

    if (time.check_time_gate(input.get_dt_write_restarts())) {
      report.print(3, "Writing restart files");
      neutrals.restart_file(input.get_restartout_dir(), DoWrite);
      ions.restart_file(input.get_restartout_dir(), DoWrite);
      time.restart_file(input.get_restartout_dir(), DoWrite);
    }
  }

  if (didWork)
    didWork = output(neutrals, ions, gGrid, time, planet);
  if (didWork)
    didWork = output(neutralsMag, ionsMag, mGrid, time, planet);

  if (didWork)
    didWork = logfile.write_logfile(indices, neutrals, ions, gGrid, time);
  if (didWork)
    didWork = logfileMag.write_logfile(indices, neutralsMag, ionsMag, mGrid, time);

  if (!didWork)
    report.error("Error in Advance!");

  report.exit(function);
  return didWork;
}
