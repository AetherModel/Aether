// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "aether.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

Logfile::Logfile(std::string logfileNameIn,
                 precision_t dtIn,
                 bool doAppendIn,
                 Indices indices,
                 Inputs inputs,
                 Report &report) {

  std::string function = "Logfile::Logfile";
  static int iFunction = -1;
  report.enter(function, iFunction);

  isOk = true;

  doAppend = doAppendIn;
  logfileName = logfileNameIn;
  dt = dtIn;

  if (doAppend)
    logfilestream.open(logfileName, std::ofstream::app);
  else
    logfilestream.open(logfileName, std::ofstream::trunc);

  logfilestream << "year month day hour minute second milli ";
  std::vector<std::string> species_array = inputs.get_species_vector();

  for (int i = 0; i < indices.all_indices_array_size(); i++) {
    if (indices.get_nValues(i) > 0)
      logfilestream << indices.get_name(i) << " ";
  }

  logfilestream << "min_temp max_temp mean_temp ";
  std::string name;

  for (int i = 0; i < species_array.size(); i++) {
    name = species_array.at(i);
    logfilestream << name << "_min " << name << "_max " << name << "_mean ";
  }

  logfilestream << "specific_temp";
  logfilestream << "\n";

  if (doAppend)
    logfilestream.close();

  report.exit(function);
}


//-------------------------------------------------------------
// Add a new line to the log file
//-------------------------------------------------------------

bool Logfile::write_logfile(Times time,
                            Neutrals neutrals,
                            Ions ions,
                            Inputs inputs,
                            Indices indices,
                            Report &report) {

  std::string function = "Logfile::write_logfile";
  static int iFunction = -1;
  report.enter(function, iFunction);

  bool didWork = true;

  if (time.check_time_gate(dt)) {

    // If we openned in append mode, we have to re-open each time:
    if (doAppend)
      logfilestream.open(logfileName, std::ofstream::app);

    // Output time first:
    std::vector<int> itime = time.get_iCurrent();
    std::vector<std::string> species_array = inputs.get_species_vector();

    for (int i = 0; i <= 6; ++i)
      logfilestream << itime.at(i) << " ";

    // Output indices next:
    logfilestream.precision(4);

    for (int i = 0; i < indices.all_indices_array_size(); i++) {
      if (indices.get_nValues(i) > 0)
        logfilestream << indices.get_index(time.get_current(), i) << " ";
    }

    // Output Neutral Temperatures:
    std::vector<precision_t> temp_stats;
    temp_stats = get_min_mean_max(neutrals.temperature_scgc);

    for (int i = 0; i < 3; i++)
      logfilestream << temp_stats[i] << " ";

    // Output densities of requested species:
    for (int i = 0; i < species_array.size(); i++) {
      std::vector<precision_t> density_stats =
        get_min_mean_max_density(species_array[i], neutrals, ions, report);

      for (int i = 0; i < 3; i++)
        logfilestream << density_stats[i] << " ";
    }

    // Output closest temp given lat, lon, alt
    logfilestream << neutrals.temperature_scgc(lla.at(0), lla.at(1), lla.at(2));
    logfilestream << "\n";

    if (doAppend)
      logfilestream.close();
  } // time gate

  report.exit(function);
  return didWork;
}

void Logfile::close_logfile() {
  logfilestream.close();
}

