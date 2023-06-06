// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "aether.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

// -------------------------------------------------------------------
// Helper function to print vector
// Only used in this file and neither declared nor visible in any other file
// -------------------------------------------------------------------

void write_logfile_line(std::ofstream &ofs,
                        const std::vector<int> &iCurrent,
                        const std::vector<precision_t> &indices,
                        const std::vector<precision_t> &variables) {
    ofs << iCurrent[0];
    for (size_t i = 1; i < iCurrent.size(); ++i) {
        ofs << ' ' << iCurrent[i];
    }
    for (auto it : indices) {
        ofs << ' ' << it;
    }
    for (auto it : variables) {
        ofs << ' ' << it;
    }
    ofs << '\n';
}

// -------------------------------------------------------------------
// Construct the satellite class using file at csv_in
// -------------------------------------------------------------------

Satellite::Satellite(const std::string &csv_in, const precision_t dt_in) {
    // Set the time gate
    dt = dt_in;

    // Open the input file
    std::ifstream ifs(csv_in);
    if (!ifs.is_open()) {
        throw std::string("Satellite: Can not open csv file " + csv_in);
    }

    // Read the first line
    const std::string sub = "Satellite: ";
    getline(ifs, name);
    // Erase the prefix "Satellite: " if its format is correct
    size_t pos_find = name.find(sub);
    if (pos_find == 0) {
        name.erase(0, sub.length());
    } else {
        throw std::string("Satellite file " + csv_in + ": Line 1 format error. "
                          "It should be \"Satellite: <name>\"");
    }

    // The next line should be exactly the same as the following
    std::string line;
    const std::string line2 = "year, mon, day, hr, min, sec, lon (deg), lat (deg), alt (km), "
                              "x (km), y (km), z (km), vx (km/s), vy (km/s), vz (km/s)";
    getline(ifs, line);
    if (line != line2) {
        throw std::string("Satellite file " + csv_in + ": Line 2 format error");
    }

    // Read time and position on each line
    std::vector<int> itime(7);
    precision_t lon_in, lat_in, alt_in;
    // CAUTION: THE SATELLITE INPUT FILE DOES NOT HAVE milliseconds
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
                throw std::string("Satellite file " + csv_in + ": Time is not strictly increasing");
            }
            timereals.push_back(timereal);
            // Convert degree to radians, kilometers to meters
            lons.push_back(lon_in * cDtoR);
            lats.push_back(lat_in * cDtoR);
            alts.push_back(alt_in * cKMtoM);
        } else {
            throw std::string("Satellite file " + csv_in + ": Data line format error");
        }
    }

    // Create the output file stream
    // output_dir fixed to be "UA/output/"
    std::string log_name = "UA/output/SAT_" + name + '_' + cGrid + "_log.txt";
    ofs.open(log_name);
    if (!ofs.is_open()) {
        throw std::string("Satellite: Can not open the output logfile");
    }
    ofs.precision(4);
}

// -------------------------------------------------------------------
// Close the output file stream
// -------------------------------------------------------------------

Satellite::~Satellite() {
    ofs.close();
}

// -------------------------------------------------------------------
// Get the position of the satellite at given time
// -------------------------------------------------------------------

void Satellite::get_position(std::vector<precision_t> &lons_out,
                             std::vector<precision_t> &lats_out,
                             std::vector<precision_t> &alts_out,
                             const double time_in) const {
    // If the input is smaller than the first element or greater than
    // the last element, return false
    if (time_in < timereals.front() || time_in > timereals.back()) {
        throw std::string("Satellite:: The time_in for get_position is out of bound");
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

    // Calculate the interpolation coefficient
    precision_t ratio = (time_in - timereals[index]) / (timereals[index + 1] - timereals[index]);

    // std::cout << "Using index = " << index << std::endl;
    // std::cout << "Ratio = " << ratio << std::endl;

    // Do the interpolation, special treatment for longitude
    // If it cross the prime meridian, first add 2pi to the latter item, and then wrap the result
    if (lons[index + 1] < lons[index]) {
        precision_t lon_out = linear_interpolation(lons[index], lons[index + 1] + cTWOPI, ratio);
        if (lon_out >= cTWOPI) {
            lon_out -= cTWOPI;
        }
        lons_out.push_back(lon_out);
    } else {
        lons_out.push_back(linear_interpolation(lons[index], lons[index + 1], ratio));
    }
    lats_out.push_back(linear_interpolation(lats[index], lats[index + 1], ratio));
    alts_out.push_back(linear_interpolation(alts[index], alts[index + 1], ratio));
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
// Write the content to the satellite log
// -------------------------------------------------------------------

void Satellite::write_log(const std::vector<int> &iCurrent,
                          const std::vector<precision_t> &indices,
                          const std::vector<precision_t> &varialbes) {
    write_logfile_line(ofs, iCurrent, indices, varialbes);
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

    // Store the locations and indexes of satellites whose time gate is open
    std::vector<precision_t> lons;
    std::vector<precision_t> lats;
    std::vector<precision_t> alts;
    std::vector<size_t> idx;
    double current = time.get_current();
    for (size_t i = 0; i < satellites.size(); ++i) {
        if (time.check_time_gate(satellites[i].get_dt())) {
            satellites[i].get_position(lons, lats, alts, current);
            idx.push_back(i);
        }
    }

    for (size_t i = 0; i < idx.size(); ++i) {
        std::cout << lons[i] << '\t' << lats[i] << '\t' << alts[i] << '\n';
    }

    // If none of the satellites' time gate is open and (iProc != 0 or time gate
    // for general logfile is not open), then nor work to be done
    if (idx.empty() && (iProc != 0 || !time.check_time_gate(dt))) {
        return true;
    }

    // Generate the output of time and all indices first
    std::vector<precision_t> values;
    std::vector<int> iCurrent = time.get_iCurrent();
    for (int i = 0; i < indices.all_indices_array_size(); ++i) {
        values.push_back(indices.get_index(current, i));
    }

    // The contents for general log and satellite make difference from here

    // Store the log string for each satellite and the result of interpolation
    std::vector<std::vector<precision_t>> values_sat(idx.size());
    std::vector<precision_t> interp_result;

    // Set the interpolation coefficients using the location of satellites first
    if (!gGrid.set_interpolation_coefs(lons, lats, alts)) {
        std::cout << "Logfile: Can not set interpolation coefficients!\n";
        return false;
    }

    // Temperature
    interp_result = gGrid.get_interpolation_values(neutrals.temperature_scgc);
    for (size_t i = 0; i < interp_result.size(); ++i) {
        values_sat[i].push_back(interp_result[i]);
    }
    // Specified variables
    for (auto &it : species) {
        const arma_cube &density = find_species_density(it, neutrals, ions, report);
        interp_result = gGrid.get_interpolation_values(density);
        for (size_t i = 0; i < interp_result.size(); ++i) {
            values_sat[i].push_back(interp_result[i]);
        }
    }
    
    // Write both the common and special values to satellite log
    for (size_t i = 0; i < values_sat.size(); ++i) {
        satellites[idx[i]].write_log(iCurrent, values, values_sat[i]);
    }

    // Return here if there is only satellite work and no gerneral log
    if (iProc != 0 || !time.check_time_gate(dt)) {
        return true;
    }

    // Open the file stream if append
    if (doAppend) {
        logfilestream.open(logfileName, std::ofstream::app);
        logfilestream.precision(4);
    }
    // Report error if the file stream is not open
    if (!logfilestream.is_open()) {
        std::cout << "Logfile: Can not open output file stream!\n";
        return false;
    }

    std::vector<precision_t> values_log;
    std::vector<precision_t> min_mean_max;
    // Temperature
    min_mean_max = get_min_mean_max(neutrals.temperature_scgc);
    values_log.insert(values_log.end(), min_mean_max.begin(), min_mean_max.end());
    // Specified variables
    for (auto it : species) {
        min_mean_max = get_min_mean_max_density(it, neutrals, ions, report);
        values_log.insert(values_log.end(), min_mean_max.begin(), min_mean_max.end());
    }
    // Temperature at the chosen point
    values_log.push_back(neutrals.temperature_scgc(lla[0], lla[1], lla[2]));
    // Output
    write_logfile_line(logfilestream, iCurrent, values, values_log);

    // Close the file stream if append
    if (doAppend) {
        logfilestream.close();
    }

    report.exit(function);
    return true;
}
