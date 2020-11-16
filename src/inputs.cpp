
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "../include/sizes.h"
#include "../include/constants.h"
#include "../include/times.h"
#include "../include/inputs.h"
#include "../include/file_input.h"
#include "../include/time_conversion.h"
#include "../include/report.h"

Inputs::Inputs(Times &time) {

  int iErr;

  iErr = 0;

  // ------------------------------------------------
  // Set some defaults:

  iVerbose = 3;
  euv_model="euvac";
  planet = "Earth";

  // ------------------------------------------------
  // Grid Defaults:
  IsUniformAlt = 1;
  alt_min = 100.0 * 1000.0;
  dalt = 5.0 * 1000.0;

  IsUniformAlt = 0;
  alt_min = 100.0 * 1000.0;
  dalt = 0.33;

  if (nGeoLons == 1) {
    lon_min = 0.0;
    lon_max = 0.0;
  } else {
    lon_min = 0.0;
    lon_max = 2.0*pi;
  }

  if (nGeoLats == 1) {
    lat_min = 0.0;
    lat_max = 0.0;
  } else {
    lat_min = -pi;
    lat_max = pi;
  }

  euv_heating_eff_neutrals = 0.40;
  euv_heating_eff_electrons = 0.05;

  dt_output.push_back(300.0);
  dt_euv = 60.0;

  // ------------------------------------------------
  // Now read the input file:
  iErr = read(time);
  
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

std::string Inputs::get_euv_file() {
  return euv_file;
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

int Inputs::read(Times &time) {

  int iErr;
  std::string line, hash;
  std::vector<int> itime(7,0);

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
      if (test_verbose(3, iVerbose))
	std::cout << "hash : -->" << hash << "<--\n";

      // ---------------------------
      // #debug or #verbose
      // ---------------------------

      if (hash == "#debug"  || hash == "#verbose") {
	iVerbose = read_int(infile_ptr, hash);
      }

      // ---------------------------
      // #starttime
      // ---------------------------

      if (hash == "#starttime") {
	std::vector<int> istart = read_itime(infile_ptr, hash);
	if (istart[0] > 0) time.set_times(istart);
	if (test_verbose(3, iVerbose)) {
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
	if (test_verbose(3, iVerbose)) {
	  std::cout << "Endtime : ";
	  display_itime(iend);
	}
      }

      // ---------------------------
      // #f107file
      // ---------------------------

      if (hash == "#f107file") {
	f107_file = read_string(infile_ptr, hash);
      }

      // ---------------------------
      // #planet
      // ---------------------------

      if (hash == "#planet") {
	planet = read_string(infile_ptr, hash);
      }

    }

    infile_ptr.close();

  }

  return iErr;


}
