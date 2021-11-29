// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"

// -----------------------------------------------------------------------
// Read a generic json input file
// -----------------------------------------------------------------------
// this should be in tools now
//json Inputs::read_json(std::string json_file) {
//
//  int iErr = 0;
//
//  json json_inputs;
//  std::ifstream infile_ptr;
//  infile_ptr.open(json_file);
//
//  if (!infile_ptr.is_open())
//    std::cout << "Could not open input file: " << json_file << "!!!\n";
//
//  else
//    infile_ptr >> json_inputs;
//
//  return json_inputs;
//}

// -----------------------------------------------------------------------
// Read input file - json format
// -----------------------------------------------------------------------

int Inputs::read_inputs_json(Times &time, Report &report) {

  int iErr = 0;

  json defaults;
  json user_inputs;

  // Set the default values first:
  settings = read_json("UA/inputs/defaults.json");
  // Set the planet-specific file (user can change this in aether.in file!):
  settings["PlanetSpeciesFile"] =
    "UA/inputs/" + get_settings_str("Planet") + ".in";

  // Then read in user perturbations on those defaults:
  user_inputs = read_json("aether.json");

  // Merge the two, with settings being the default:
  settings.merge_patch(user_inputs);

  // Debug Stuff:
  report.set_verbose(settings["Debug"]["iVerbose"]);
  report.set_timing_depth(settings["Debug"]["iTimingDepth"]);

  std::vector<int> istart = get_settings_timearr("StartTime");
  time.set_times(istart);

  std::vector<int> iend = get_settings_timearr("EndTime");
  time.set_end_time(iend);

  return iErr;
}

