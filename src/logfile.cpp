// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "aether.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

//-------------------------------------------------------------
// Initialize the Logfile
//-------------------------------------------------------------

Logfile::Logfile(Indices &indices,
                 Inputs &input,
                 Report &report) {

    std::string function = "Logfile::Logfile";
    static int iFunction = -1;
    report.enter(function, iFunction);

    // Read the inputs
    logfileName = input.get_logfile();
    species = input.get_species_vector();
    dt = input.get_logfile_dt();
    doAppend = input.get_logfile_append();

    // Initialize other variables
    isOk = true;

    // Write the header
    // Open the file stream
    if (doAppend) {
        logfilestream.open(logfileName, std::ofstream::app);
    } else {
        logfilestream.open(logfileName, std::ofstream::trunc);
        logfilestream.precision(4);
    }
    // Report error if can not open the file stream
    if (!logfilestream.is_open()) {
        std::cout << "Could not open log file: " << logfileName << "!!!\n";
        isOk = false;
    }
    // The name of all indicies
    logfilestream << "year month day hour minute second milli ";
    for (int i = 0; i < indices.all_indices_array_size(); ++i) {
        logfilestream << indices.get_name(i) << ' ';
    }
    // The temperature
    logfilestream << "min_temp mean_temp max_temp ";
    // The specified variables
    for (auto &it : species) {
        logfilestream << it << "_min " << it << "_mean " << it << "_max ";
    }
    // The temp at the chosen point
    logfilestream << "specific_temp\n";

    // Close the file stream if append
    if (doAppend) {
        logfilestream.close();
    }

    report.exit(function);
}

//-------------------------------------------------------------
// Close the file stream if not append
//-------------------------------------------------------------

Logfile::~Logfile() {
    if (!doAppend) {
        logfilestream.close();
    }
}

//-------------------------------------------------------------
// Add a new line to the log file
//-------------------------------------------------------------

bool Logfile::write_logfile(Indices &indices,
                            Neutrals &neutrals,
                            Ions &ions,
                            Grid &gGrid,
                            Times &time,
                            Report &report) {

    std::string function = "Logfile::write_logfile";
    static int iFunction = -1;
    report.enter(function, iFunction);

    if (!time.check_time_gate(dt)) {
        // No work to be done
        return true;
    }

    // Open the file stream if append
    if (doAppend) {
        logfilestream.open(logfileName, std::ofstream::app);
        logfilestream.precision(4);
    }
    // Report error if the file stream is not open
    if (!logfilestream.is_open()) {
        isOk = false;
        return false;
    }

    // Get the current time
    std::vector<int> iCurrent = time.get_iCurrent();
    double current = time.get_current();

    // Display the year, month, day, hour, minute, second, and millisecond
    for (auto &it : iCurrent) {
        logfilestream << it << ' ';
    }

    // Display all indices
    for (int i = 0; i < indices.all_indices_array_size(); ++i) {
        logfilestream << indices.get_index(current, i) << ' ';
    }

    // Display temperatures
    std::vector<precision_t> min_mean_max;
    min_mean_max = get_min_mean_max(neutrals.temperature_scgc);
    for (auto &it : min_mean_max) {
        logfilestream << it << ' ';
    }

    // Display specified variables
    for (auto &it : species) {
        min_mean_max = get_min_mean_max_density(it, neutrals, ions, report);
        for (auto &it2 : min_mean_max) {
            logfilestream << it2 << ' ';
        }
    }

    // Display the temperature at the chosen point
    logfilestream << neutrals.temperature_scgc(lla[0], lla[1], lla[2]) << '\n';

    // Close the file stream if append
    if (doAppend) {
        logfilestream.close();
    }

    report.exit(function);
    return true;
}
