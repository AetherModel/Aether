// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "aether.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

// -------------------------------------------------------------------
// Construct the satellite class using file at csv_in
// -------------------------------------------------------------------

Satellite::Satellite(const std::string &csv_in, const precision_t dt_in) {
    // Set the time gate
    dt = dt_in;

    // Open the input file
    std::ifstream ifs(csv_in);
    if (!ifs.is_open()) {
        throw std::string("Can not open csv file");
    }

    // Read the first line
    const std::string sub = "Satellite: ";
    getline(ifs, name);
    // Erase the prefix "Satellite: " if its format is correct
    size_t pos_find = name.find(sub);
    if (pos_find == 0) {
        name.erase(0, sub.length());
    } else {
        throw std::string("Line 1 format error. It should be \"Satellite: <name>\"");
    }

    // The next line should be exactly the same as the following
    std::string line;
    const std::string line2 = "year, mon, day, hr, min, sec, lon (deg), lat (deg), alt (km), x (km), y (km), z (km), vx (km/s), vy (km/s), vz (km/s)";
    getline(ifs, line);
    if (line != line2) {
        throw std::string("Line 2 format error");
    }

    // Read time and position on each line
    std::vector<int> itime(7);
    precision_t lon_in, lat_in, alt_in;
    // CAUTION: THE SATELLITE INPUT FILE DO NOT HAS milliseconds
    // Fix the millisecond to be 0
    itime[6] = 0;    
    while (getline(ifs, line) && line.length() > 1) {
        std::istringstream iss(line);
        // ch is used to read comma
        char ch;
        if (iss >> itime[0] >> ch >> itime[1] >> ch >> itime[2] >> ch >> itime[3] >> ch
                >> itime[4] >> ch >> itime[5] >> ch >> lon_in >> ch >> lat_in >> ch
                >> alt_in) {
            double timereal = time_int_to_real(itime);
            // Exit if the time is not ascending
            if (!timereals.empty() && timereal <= timereals.back()) {
                display_itime(time_real_to_int(timereals.back()));
                display_itime(itime);
                throw std::string("Time is not strictly increasing");
            }
            timereals.push_back(timereal);
            // Convert degree to radians, kilometers to meters
            lons.push_back(lon_in * cDtoR);
            lats.push_back(lat_in * cDtoR);
            alts.push_back(alt_in * cKMtoM);
        } else {
            throw std::string("Data line format error");
        }
    }

    // Create the output file stream
    // output_dir fixed to be "UA/output/"
    std::string log_name = "UA/output/SAT_" + name + '_' + cGrid + "_log.txt";
    ofs.open(log_name);
    if (!ofs.is_open()) {
        throw std::string("Can not open the output logfile");
    }
}

// -------------------------------------------------------------------
// Close the output file stream
// -------------------------------------------------------------------

Satellite::~Satellite() {
    ofs.close();
}

// -------------------------------------------------------------------
// Return the position of the satellite at the given time
// -------------------------------------------------------------------

std::vector<precision_t> Satellite::get_position(const double time_in) const {
    // If the input is smaller than the first element or greater than
    // the last element, report input error
    if (time_in < timereals.front() || time_in > timereals.back()) {
        return std::vector<precision_t>();
    }

    // Use binary search to find the last element that is smaller than
    // or equal to the given value
    auto it = std::upper_bound(timereals.cbegin(), timereals.cend(), time_in);

    // Get the index
    size_t index = std::distance(timereals.cbegin(), it) - 1;
    // If it is the last element, decrement it by 1 so that
    // the point is promised to be inside [index, index+1]
    if (index == timereals.size() - 1) {
        --index;
    }

    // std::cout << "Using index = " << index << std::endl;

    // Get the two locations for interpolation
    const std::vector<precision_t> x0 = {lons[index], lats[index], alts[index]};
    std::vector<precision_t> x1 = {lons[index + 1], lats[index + 1], alts[index + 1]};

    // Calculate the interpolation coefficient
    precision_t ratio = (time_in - timereals[index]) / (timereals[index + 1] - timereals[index]);

    // std::cout << "Ratio = " << ratio << std::endl;

    // Do the interpolation
    // If it cross the prime meridian, first add 2pi to the latter item,
    // and then wrap the result
    if (x1[0] < x0[0]) {
        x1[0] += cTWOPI;
    }
    for (int i = 0; i < 3; ++i) {
        x1[i] = linear_interpolation(x0[i], x1[i], ratio);
    }
    if (x1[0] >= cTWOPI) {
        x1[0] -= cTWOPI;
    }

    return x1;
}

// -------------------------------------------------------------------
// Get the name of the satellite
// -------------------------------------------------------------------

std::string Satellite::get_name() const {
    return name;
}

// -------------------------------------------------------------------
// Get the time gate of the satellite
// -------------------------------------------------------------------

precision_t Satellite::get_dt() const {
    return dt;
}

// -------------------------------------------------------------------
// Satellite debug
// -------------------------------------------------------------------

void Satellite::print() {
    assert(timereals.size() == lons.size()
             && lons.size() == lats.size()
             && lats.size() == alts.size());
    for (size_t i = 0; i < timereals.size(); ++i) {
        std::cout << "At time = " << timereals[i] - timereals.front() << "\tpos = ("
                  << lons[i] << ", " << lats[i] << ", " << alts[i] << ")\n";
    }
}

//-------------------------------------------------------------
// Initialize the Logfile
//-------------------------------------------------------------

Logfile::Logfile(Indices &indices,
                 Inputs &input,
                 Report &report) {

    // Read the settings for general log file and satellites
    // Write the header to the general logfile
    // Write the header and the number of processors to the satellite combination helper
    // Initialize the satellite instances
    // Write the name of satellites to the satellite combination helper

    std::string function = "Logfile::Logfile";
    static int iFunction = -1;
    report.enter(function, iFunction);

    // Read the inputs
    logfileName = input.get_logfile();
    species = input.get_species_vector();
    dt = input.get_logfile_dt();
    doAppend = input.get_logfile_append();
    std::vector<std::string> sat_names = input.get_satellite_files();
    std::vector<precision_t> sat_dts = input.get_satellite_dt();

    // Initialize other variables
    isOk = true;

    // Open the general log file stream
    if (doAppend) {
        logfilestream.open(logfileName, std::ofstream::app);
    } else {
        logfilestream.open(logfileName, std::ofstream::trunc);
        logfilestream.precision(4);
    }
    // Open the satellite combination helper file stream
    std::ofstream sat_combine("UA/output/SAT_COMBINE.txt");
    // Report error if can not open the file stream
    if (!logfilestream.is_open() || !sat_combine.is_open()) {
        std::cout << "Could not open log file: " << logfileName << "!!!\n";
        isOk = false;
        // TRY TO EXIT GRACEFULLY HERE. ALL THE FOLLOWING CODE SHOULD NOT BE EXECUTED
        throw std::string("Could not open log file");
    }

    // Write the header to both gerneral log and combination helper
    std::string header;
    // The name of all indicies
    header = "year month day hour minute second milli ";
    for (int i = 0; i < indices.all_indices_array_size(); ++i) {
        header += indices.get_name(i) + ' ';
    }
    // The header for log and satellite make difference from here
    // The general log takes average, while the satellite takes the value at its point
    std::string header_sat = header;
    // The temperature
    header += "min_temp mean_temp max_temp ";
    header_sat += "temp ";
    // The specified variables
    for (auto &it : species) {
        header += it + "_min " + it + "_mean " + it + "_max ";
        header_sat += it + ' ';
    }
    // The temp at the chosen point only for general logfile
    header += "specific_temp\n";
    header_sat += '\n';

    // Write the header to log
    logfilestream << header;
    // Write the header and the total number of grids to sat_combine
    sat_combine << header_sat << nGrids << '\n';

    // Close the file stream if append
    if (doAppend) {
        logfilestream.close();
    }

    // Check that the size of files and dt are the same
    if (sat_names.size() != sat_dts.size()) {
        std::cout << "The size of files and dt for satellites are not the same!!!\n";
        isOk = false;
        // TRY TO EXIT GRACEFULLY HERE. ALL THE FOLLOWING CODE SHOULD NOT BE EXECUTED
        throw std::string("The size of files and dt for satellites are not the same");
    }
    // Initialize the instances of satellite
    for (size_t i = 0; i < sat_names.size(); ++i) {
        satellites.emplace_back(sat_names[i], sat_dts[i]);
    }
    // Write the name of the satellites to the sat_combine and close the file stream
    for (auto &it : satellites) {
        sat_combine << it.get_name() << ' ';
    }
    sat_combine << '\n';
    sat_combine.close();

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
