// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "aether.h"

// -----------------------------------------------------------------------
// Read input file - json format
// There are three potential input files:
//   - UA/inputs/defaults.json - sets the defaults for the simulation
//   - UA/restartIn/settings.json - sets the state for the old simulation
//   - aether.json - perturbs the state to what the user wants
// So, logic below is:
//   - read in defaults file
//   - read in user file
//   - if user file asks for a restart, then
//      - read in the restart settings file
//      - merge with the defaults file
//   - merge the user settings, so those overwrite anything being
//     asked for.
// -----------------------------------------------------------------------

bool Inputs::read_inputs_json(Times &time) {

  json defaults;
  json user_inputs;

  isOk = true;

  // Set the default values first:
  settings = read_json("UA/inputs/defaults.json");
  isOk = set_verbose(settings);

  try {

    // Then read in user perturbations on those defaults:
    user_inputs = read_json("aether.json");
    isOk = set_verbose(user_inputs);

    // Read in a restart file also if user specified it.
    //   - Here we merge the restart inputs with the defaults inputs
    //   - This is BEFORE the user inputs are merged!!!

    if (user_inputs.contains("Restart")) {
      if (user_inputs["Restart"].contains("do")) {
        if (user_inputs["Restart"]["do"]) {
          std::string restart_file = get_setting_str("Restart", "InDir");
          restart_file = restart_file + "/settings.json";
          json restart_inputs;
          restart_inputs = read_json(restart_file);
          // This forces the logfile to append.  User can override
          // if they really want:
          restart_inputs["Logfile"]["append"] = true;
          settings.merge_patch(restart_inputs);
        }
      }
    }

    // Merge the defaults/restart settings with the user provided
    // settings, with the default/restart settings being the default:
    settings.merge_patch(user_inputs);

    //change planet file to the one specified on aether.json:
    if (isOk)
      settings["PlanetSpeciesFile"] = get_setting_str("PlanetFile");

    std::string planet_filename = get_setting_str("PlanetSpeciesFile");
    report.print(1, "Using planet file : " + planet_filename);

    // Debug Stuff:
    if (isOk)
      report.set_verbose(get_setting_int("Debug", "iVerbose"));

    if (isOk)
      report.set_DefaultVerbose(get_setting_int("Debug", "iVerbose"));

    if (isOk)
      report.set_doInheritVerbose(get_setting_bool("Debug", "doInheritVerbose"));

    if (isOk)
      report.set_timing_depth(get_setting_int("Debug", "iTimingDepth"));

    if (isOk)
      report.set_timing_percent(get_setting_float("Debug", "TimingPercent"));

    if (isOk)
      report.set_iProc(get_setting_int("Debug", "iProc"));

    for (auto &item : settings["Debug"]["iFunctionVerbose"].items())
      report.set_FunctionVerbose(item.key(), item.value());

    // Capture time information:
    if (isOk)
      time.set_times(get_setting_timearr("StartTime"));

    if (isOk)
      time.set_end_time(get_setting_timearr("EndTime"));

  } catch (...) {
    report.error("Error in reading inputs!");
    isOk = false;
  }

  return isOk;
}
