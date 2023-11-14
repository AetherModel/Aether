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
                        const std::vector<precision_t> &variables) {
  // The first element is printed separately to achieve no trailing
  // white space
  ofs << iCurrent[0];

  for (size_t i = 1; i < iCurrent.size(); ++i)
    ofs << ' ' << iCurrent[i];

  for (auto it : variables)
    ofs << ' ' << it;

  ofs << '\n';
}

// -------------------------------------------------------------------
// Construct the satellite class using file at csv_in
// -------------------------------------------------------------------

Satellite::Satellite(const std::string &csv_in,
                     const std::string &log_name_in,
                     const std::string &sat_header,
                     const precision_t dt_in) {
  // Set the time gate
  dt = dt_in;

  // Open the input file
  std::ifstream ifs(csv_in);

  if (!ifs.is_open())
    throw std::string("Satellite: Can not open csv file " + csv_in);

  // Read the first line
  const std::string sub = "Satellite: ";
  getline(ifs, name);
  // Erase the prefix "Satellite: " if its format is correct
  size_t pos_find = name.find(sub);

  if (pos_find == 0)
    name.erase(0, sub.length());

  else {
    throw std::string("Satellite file " + csv_in + ": Line 1 format error. "
                      "It should be \"Satellite: <name>\"");
  }

  // The next line should be exactly the same as the following
  std::string line;
  const std::string line2 =
    "year, mon, day, hr, min, sec, lon (deg), lat (deg), alt (km), "
    "x (km), y (km), z (km), vx (km/s), vy (km/s), vz (km/s)";
  getline(ifs, line);

  if (line != line2) {
    throw std::string("Satellite file " + csv_in + ": Line 2 format error. "
                      "It should be \n" + line2);
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

    if (iss >> itime[0] >> ch
        >> itime[1] >> ch
        >> itime[2] >> ch
        >> itime[3] >> ch
        >> itime[4] >> ch
        >> itime[5] >> ch
        >> lon_in >> ch >> lat_in >> ch >> alt_in) {
      double timereal = time_int_to_real(itime);

      // Exit if the time is not ascending
      if (!timereals.empty() && timereal <= timereals.back()) {
        throw std::string("Satellite file " + csv_in +
                          ": Time is not strictly increasing");
      }

      timereals.push_back(timereal);
      // Convert degree to radians, kilometers to meters
      lons.push_back(lon_in * cDtoR);
      lats.push_back(lat_in * cDtoR);
      alts.push_back(alt_in * cKMtoM);
    } else {
      throw std::string("Satellite file " + csv_in +
                        ": Data line format error");
    }
  }

  // Store the log_name_in
  log_name = log_name_in;

  if (nMembers > 1)
    log_name += '_' + cMember;

  log_name += ".txt";

  // Grid 0 will initialize the output file
  if (iGrid == 0) {
    std::ofstream ofs(log_name);

    if (!ofs.is_open())
      throw std::string("Satellite: Can not open the output logfile");

    ofs << sat_header << '\n';
    ofs.close();
  }
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
  if (time_in < timereals.front() || time_in > timereals.back())
    throw std::string("Satellite:: The time_in for get_position is out of bound");

  // Use binary search to find the last element that is smaller than
  // or equal to the given value
  auto it = std::upper_bound(timereals.cbegin(), timereals.cend(), time_in);

  // Get the index
  size_t index = std::distance(timereals.cbegin(), it) - 1;

  // If it is the last element, decrement it by 1 so that
  // the point is promised to be inside [index, index+1]
  if (index == timereals.size() - 1)
    --index;

  // Calculate the interpolation coefficient
  precision_t ratio =
    (time_in - timereals[index]) /
    (timereals[index + 1] - timereals[index]);

  // Do the interpolation, special treatment for longitude: If it
  // crosses the prime meridian, first add 2pi to the latter item,
  // and then wrap the result
  if (lons[index + 1] < lons[index]) {
    precision_t lon_out =
      linear_interpolation(lons[index], lons[index + 1] + cTWOPI, ratio);

    if (lon_out >= cTWOPI)
      lon_out -= cTWOPI;

    lons_out.push_back(lon_out);
  } else {
    lons_out.push_back(linear_interpolation(lons[index],
                                            lons[index + 1],
                                            ratio));
  }

  lats_out.push_back(linear_interpolation(lats[index],
                                          lats[index + 1],
                                          ratio));
  alts_out.push_back(linear_interpolation(alts[index],
                                          alts[index + 1],
                                          ratio));
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
                          const std::vector<precision_t> &varialbes) {
  std::ofstream ofs(log_name, std::ios_base::app);
  ofs.precision(5);
  write_logfile_line(ofs, iCurrent, varialbes);
  ofs.close();
}

// -------------------------------------------------------------------
// Satellite debug
// -------------------------------------------------------------------

void Satellite::print() {
  assert(timereals.size() == lons.size()
         && lons.size() == lats.size()
         && lats.size() == alts.size());

  for (size_t i = 0; i < timereals.size(); ++i) {
    std::cout << "At time = "
              << timereals[i] - timereals.front()
              << "\tpos = ("
              << lons[i] << ", " << lats[i] << ", " << alts[i] << ")\n";
  }
}

//-------------------------------------------------------------
// Initialize the Logfile
//-------------------------------------------------------------

Logfile::Logfile(Indices &indices) {

  // Read the settings for general log file and satellites
  // Write the header to the general logfile
  // Initialize the satellite instances

  std::string function = "Logfile::Logfile";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // Read the inputs
  logfileName = input.get_logfile();
  species = input.get_species_vector();
  dt = input.get_logfile_dt();
  doAppend = input.get_logfile_append();
  std::vector<std::string> sat_files = input.get_satellite_files();
  std::vector<std::string> sat_names = input.get_satellite_names();
  std::vector<precision_t> sat_dts = input.get_satellite_dts();

  bool isRestart = input.get_do_restart();

  // Only the 0th processor in each member needs to open the logfile

  if (iGrid == 0) {
    // Open the general log file stream
    if (doAppend)
      logfilestream.open(logfileName, std::ofstream::app);
    else {
      logfilestream.open(logfileName, std::ofstream::trunc);
      logfilestream.precision(4);
    }

    // Report error if can not open the log file stream
    if (!logfilestream.is_open() & !isRestart) {
      // TRY TO EXIT GRACEFULLY HERE. ALL THE FOLLOWING CODE SHOULD
      // NOT BE EXECUTED
      throw std::string("Can not open log file");
    }
  }

  // The header of time. They are placed at the beginning of all values
  std::string header_time = "year month day hour minute second milli";
  // Satellite-specific and log-specific header
  std::string header_log, header_sat;
  // The general log takes average, while the satellite takes the
  // value at its point
  // The temperature
  header_log += " min_temp mean_temp max_temp specific_temp";
  header_sat += " lon(deg) lat(deg) alt(km) temp";

  // The specified variables
  for (auto &it : species) {
    header_log += ' ' + it + "_min " + it + "_mean " + it + "_max";
    header_sat += ' ' + it;
  }

  // The header of indicies. They are only used in general log file
  for (int i = 0; i < indices.all_indices_array_size(); ++i)
    header_log += ' ' + indices.get_name(i);

  // only the 0th processor writes to the logfile
  if (logfilestream.is_open()) {
    // Write the header to log
    if (!isRestart)
      logfilestream << header_time << header_log << '\n';

    // Close the file stream if append
    if (doAppend)
      logfilestream.close();
  }

  // Check whether the input settings contain Satellites
  if (!sat_files.empty()) {
    // Check that the size of files and dts are the same
    if (sat_files.size() != sat_names.size() ||
        sat_names.size() != sat_dts.size()) {
      // TRY TO EXIT GRACEFULLY HERE. ALL THE FOLLOWING CODE
      // SHOULD NOT BE EXECUTED
      throw std::string("The size of files and dts for satellites are not the same");
    }

    // Initialize the instances of satellite
    satellites.reserve(sat_files.size());

    for (size_t i = 0; i < sat_files.size(); ++i) {
      satellites.emplace_back(sat_files[i],
                              sat_names[i],
                              header_time + header_sat,
                              sat_dts[i]);
    }
  }

  // Synchronize here to make sure all satellite logfiles are initialized and
  // prevent any process writing data into logfile before initialization
  MPI_Barrier(aether_comm);

  report.exit(function);
  return;
}

//-------------------------------------------------------------
// Close the file stream if not append
//-------------------------------------------------------------

Logfile::~Logfile() {
  // only the 0th processor in each ensemble opened the file
  if (!doAppend & iGrid == 0)
    logfilestream.close();
}

//-------------------------------------------------------------
// Add a new line to the log file
//-------------------------------------------------------------

bool Logfile::write_logfile(Indices &indices,
                            Neutrals &neutrals,
                            Ions &ions,
                            Grid &gGrid,
                            Times &time) {

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

  // Time
  std::vector<int> iCurrent = time.get_iCurrent();

  // Check if there are open time gates for satellites
  if (!idx.empty()) {
    // Store the values to display for each satellite and the
    // result of interpolation
    std::vector<std::vector<precision_t>> values_sat(idx.size());
    std::vector<precision_t> interp_result;

    // The lon, lat, and alt of each satellite
    for (size_t i = 0; i < values_sat.size(); ++i) {
      values_sat[i].push_back(lons[i] * cRtoD);
      values_sat[i].push_back(lats[i] * cRtoD);
      values_sat[i].push_back(alts[i] * cMtoKM);
    }

    // Set the interpolation coefficients using the location of satellites
    if (!gGrid.set_interpolation_coefs(lons, lats, alts)) {
      std::cout << "Logfile: Can not set interpolation coefficients!\n";
      report.exit(function);
      return false;
    }

    // Temperature
    interp_result = gGrid.get_interpolation_values(neutrals.temperature_scgc);

    for (size_t i = 0; i < interp_result.size(); ++i)
      values_sat[i].push_back(interp_result[i]);

    // Specified variables
    for (auto &it : species) {
      const arma_cube &density = find_species_density(it, neutrals, ions);
      interp_result = gGrid.get_interpolation_values(density);

      for (size_t i = 0; i < interp_result.size(); ++i)
        values_sat[i].push_back(interp_result[i]);
    }

    // Write the values to satellite log
    for (size_t i = 0; i < values_sat.size(); ++i) {
      // Use the index 3 (temp) to check whether the satellite
      // is in grid or not
      if (values_sat[i][3] != cNinf)
        satellites[idx[i]].write_log(iCurrent, values_sat[i]);
    }
  }

  // Check if the time gate for general log file is open
  if (iGrid == 0 && time.check_time_gate(dt)) {
    // Open the file stream if append
    if (doAppend) {
      logfilestream.open(logfileName, std::ofstream::app);
      logfilestream.precision(4);
    }

    // Report error if the file stream is not open
    if (!logfilestream.is_open()) {
      report.error("Logfile: Can not open output file stream!");
      report.exit(function);
      return false;
    }

    std::vector<precision_t> values_log;
    std::vector<precision_t> min_mean_max;
    // Temperature
    min_mean_max = get_min_mean_max(neutrals.temperature_scgc);
    values_log.insert(values_log.end(),
                      min_mean_max.begin(),
                      min_mean_max.end());
    // Temperature at the chosen point
    values_log.push_back(neutrals.temperature_scgc(lla[0],
                                                   lla[1],
                                                   lla[2]));

    // Specified variables
    for (auto it : species) {
      min_mean_max = get_min_mean_max_density(it, neutrals, ions);
      values_log.insert(values_log.end(),
                        min_mean_max.begin(),
                        min_mean_max.end());
    }

    // Indices
    for (int i = 0; i < indices.all_indices_array_size(); ++i)
      values_log.push_back(indices.get_index(current, i));

    // Output
    write_logfile_line(logfilestream, iCurrent, values_log);

    // Close the file stream if append
    if (doAppend)
      logfilestream.close();
  }

  report.exit(function);
  return true;
}
