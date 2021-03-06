// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "../include/aether.h"

Inputs::Inputs(Times &time, Report &report) {

  // ------------------------------------------------
  // Set some defaults:

  iVerbose = 0;
  iTimingDepth = 3;
  euv_model = "euvac";
  planet = "Earth";

  // ------------------------------------------------
  // Grid Defaults:
  geo_grid_input.alt_file = "";
  geo_grid_input.IsUniformAlt = 1;
  geo_grid_input.alt_min = 100.0 * 1000.0;
  geo_grid_input.dalt = 5.0 * 1000.0;

  nLonsGeo = 12;
  nLatsGeo = 20;
  nAltsGeo = 40;

  if (nLonsGeo == 1) {
    geo_grid_input.lon_min = 0.0;
    geo_grid_input.lon_max = 0.0;
  } else {
    geo_grid_input.lon_min = 0.0;
    geo_grid_input.lon_max = 2.0*pi;
  }

  if (nLatsGeo == 1) {
    geo_grid_input.lat_min = 0.0;
    geo_grid_input.lat_max = 0.0;
  } else {
    geo_grid_input.lat_min = -pi/2;
    geo_grid_input.lat_max = pi/2;
  }

  euv_heating_eff_neutrals = 0.40;
  euv_heating_eff_electrons = 0.05;

  dt_output.push_back(300.0);
  type_output.push_back("states");
  dt_euv = 60.0;
  dt_report = 60.0;

  // ------------------------------------------------
  // Now read the input file:
  int iErr = read(time, report);
  if (iErr > 0) std::cout << "Error in reading input file!\n";
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

Inputs::grid_input_struct Inputs::get_grid_inputs() {
  return geo_grid_input;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

std::string Inputs::get_bfield_type() {
  return bfield;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

std::string Inputs::get_euv_model() {
  return euv_model;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

float Inputs::get_euv_heating_eff_neutrals() {
  return euv_heating_eff_neutrals;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

float Inputs::get_dt_euv() {
  return dt_euv;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

float Inputs::get_dt_report() {
  return dt_report;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

float Inputs::get_n_outputs() {
  return dt_output.size();
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

int Inputs::get_nLonsGeo() {
  return nLonsGeo;
}

int Inputs::get_nLatsGeo() {
  return nLatsGeo;
}

int Inputs::get_nAltsGeo() {
  return nAltsGeo;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

float Inputs::get_dt_output(int iOutput) {
  float value = 0.0;
  int iSize = dt_output.size();
  if (iOutput < iSize) value = dt_output[iOutput];
  return value;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

std::string Inputs::get_type_output(int iOutput) {
  std::string value = "";
  int iSize = dt_output.size();
  if (iOutput < iSize) value = type_output[iOutput];
  return value;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

std::string Inputs::get_euv_file() {
  return euv_file;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

std::string Inputs::get_chemistry_file() {
  return chemistry_file;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

int Inputs::get_number_of_omniweb_files() {
  return omniweb_files.size();
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

std::vector<std::string> Inputs::get_omniweb_files() {
  return omniweb_files;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

std::string Inputs::get_f107_file() {
  return f107_file;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

std::string Inputs::get_planet() {
  return planet;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

std::string Inputs::get_planetary_file() {
  return planetary_file;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

std::string Inputs::get_planet_species_file() {
  return planet_species_file;
}

// -----------------------------------------------------------------------
// Read input file!
// -----------------------------------------------------------------------

int Inputs::read(Times &time, Report &report) {

  int iErr;
  std::string line, hash;
  std::vector<int> itime(7, 0);

  iErr = 0;

  std::ifstream infile_ptr;
  infile_ptr.open(input_file);

  if (!infile_ptr.is_open()) {
    std::cout << "Could not open input file: " << input_file << "!!!\n";
    iErr = 1;
  } else {

    while (!infile_ptr.eof()) {

      // ---------------------------
      // Find the next hash:
      // ---------------------------

      hash = find_next_hash(infile_ptr);
      if (report.test_verbose(3))
        std::cout << "hash : -->" << hash << "<--\n";

      // ---------------------------
      // #debug or #verbose
      // ---------------------------

      if (hash == "#debug"  || hash == "#verbose") {
        iVerbose = read_int(infile_ptr, hash);
        report.set_verbose(iVerbose);
  // Need to read in the timing depth somewhere
        report.set_timing_depth(iTimingDepth);
      }

      // ---------------------------
      // #starttime
      // ---------------------------

      if (hash == "#starttime") {
        std::vector<int> istart = read_itime(infile_ptr, hash);
        if (istart[0] > 0) time.set_times(istart);
        if (report.test_verbose(3)) {
          std::cout << "Starttime : ";
          display_itime(istart);
        }
      }

      // ---------------------------
      // #endtime
      // ---------------------------

      if (hash == "#endtime") {
        std::vector<int> iend = read_itime(infile_ptr, hash);
        if (iend[0] > 0) time.set_end_time(iend);
        if (report.test_verbose(3)) {
          std::cout << "Endtime : ";
          display_itime(iend);
        }
      }

      // ---------------------------
      // #geoblocksize
      // ---------------------------

      if (hash == "#geoblocksize") {
        nLonsGeo = read_int(infile_ptr, hash);
        nLatsGeo = read_int(infile_ptr, hash);
        nAltsGeo = read_int(infile_ptr, hash);
        if (report.test_verbose(3)) {
          std::cout << "nLonsGeo : " << nLonsGeo << "\n";
          std::cout << "nLatsGeo : " << nLatsGeo << "\n";
          std::cout << "nAltsGeo : " << nAltsGeo << "\n";
        }
      }

      // ---------------------------
      // #geogrid
      // ---------------------------

      if (hash == "#geogrid") {
  geo_grid_input.lon_min = read_float(infile_ptr, hash) * dtor;
  geo_grid_input.lon_max = read_float(infile_ptr, hash) * dtor;
  geo_grid_input.lat_min = read_float(infile_ptr, hash) * dtor;
  geo_grid_input.lat_max = read_float(infile_ptr, hash) * dtor;
  geo_grid_input.alt_min = read_float(infile_ptr, hash) * 1000.0;
  geo_grid_input.dalt = read_float(infile_ptr, hash) * 1000.0;
  geo_grid_input.IsUniformAlt = 1;
        if (report.test_verbose(3)) {
          std::cout << "lon min/max : "
        << geo_grid_input.lon_min
        << geo_grid_input.lon_max << "\n";
          std::cout << "lat min/max : "
        << geo_grid_input.lat_min
        << geo_grid_input.lat_max << "\n";
          std::cout << "alt min + dalt: "
        << geo_grid_input.alt_min
        << geo_grid_input.dalt << "\n";
        }
      }

      // ---------------------------
      // #f107file
      // ---------------------------

      if (hash == "#f107file") {
        f107_file = read_string(infile_ptr, hash);
      }

      // ---------------------------
      // #bfield
      // ---------------------------

      if (hash == "#bfield") {
        bfield = read_string(infile_ptr, hash);
      }

      // ---------------------------
      // #chemistry
      // ---------------------------

      if (hash == "#chemistry") {
        chemistry_file = read_string(infile_ptr, hash);
      }

      // ---------------------------
      // #omniweb
      // This can actually be called multiple times:
      // ---------------------------

      if (hash == "#omniweb") {
        omniweb_files.push_back(read_string(infile_ptr, hash));
      }

      // ---------------------------
      // #planet
      // ---------------------------

      if (hash == "#planet") {
        planet = read_string(infile_ptr, hash);
        if (report.test_verbose(3))
          std::cout << "Setting planet to : " << planet << "\n";
        if (planet_species_file.length() <= 1)
          planet_species_file = "UA/inputs/"+planet+".in";
      }

      // ---------------------------
      // #output
      // ---------------------------

      if (hash == "#output") {
        std::vector<std::vector<std::string>> csv = read_csv(infile_ptr);
        // comma separated values, with type, then dt:
        int nOutputs = csv.size();
        int iOutput;
        if (nOutputs > 1) {
          type_output[0] = csv[0][0];
          dt_output[0] = stof(csv[0][1]);
          for (iOutput = 1; iOutput < nOutputs; iOutput++) {
            type_output.push_back(csv[iOutput][0]);
            dt_output.push_back(stof(csv[iOutput][1]));
          }
          // Allow users to enter 0 for dt, so they only get the
          // output at the beginning of the run:
          for (iOutput = 0; iOutput < nOutputs; iOutput++)
            if (dt_output[iOutput] <= 0.0) dt_output[iOutput] = 1.0e32;
        } else {
          std::cout << "Something wrong with #output. Need to report...\n";
        }
      }
    }
    infile_ptr.close();
  }
  return iErr;
}
