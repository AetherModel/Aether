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

  bool DidWork = true;

  json defaults;
  json user_inputs;

  // Set the default values first:
  settings = read_json("UA/inputs/defaults.json");
  DidWork = set_verbose(settings);

  // Set the planet-specific file (user can change this in aether.in file!):
  settings["PlanetSpeciesFile"] = settings["Planet"]["file"];

  try {

    // Then read in user perturbations on those defaults:
    user_inputs = read_json("aether.json");
    DidWork = set_verbose(user_inputs);

    // Read in a restart file also if user specified it.
    //   - Here we merge the restart inputs with the defaults inputs
    //   - This is BEFORE the user inputs are merged!!!

    if (user_inputs.contains("Restart")) {
      if (user_inputs["Restart"].contains("do")) {
        if (user_inputs["Restart"]["do"]) {
          std::string restart_file = settings["Restart"]["InDir"];
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
    settings["PlanetSpeciesFile"] = settings["PlanetFile"];

    std::string planet_filename = settings["PlanetSpeciesFile"];
    report.print(1, "Using planet file : " + planet_filename);

    // Debug Stuff:
    report.set_verbose(get_settings("Debug", "iVerbose"));
    report.set_DefaultVerbose(get_settings("Debug", "iVerbose"));
    report.set_doInheritVerbose(settings["Debug"]["doInheritVerbose"]);
    report.set_timing_depth(get_settings("Debug", "iTimingDepth"));
    report.set_timing_percent(settings["Debug"]["TimingPercent"]);
    report.set_iProc(get_settings("Debug", "iProc"));

    for (auto &item : settings["Debug"]["iFunctionVerbose"].items())
      report.set_FunctionVerbose(item.key(), item.value());

    // Capture time information:
    std::vector<int> istart = get_settings_timearr("StartTime");
    time.set_times(istart);

    std::vector<int> iend = get_settings_timearr("EndTime");
    time.set_end_time(iend);
  } catch (...) {
    DidWork = false;
  }

  return DidWork;
}
