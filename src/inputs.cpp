
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "../include/sizes.h"
#include "../include/constants.h"
#include "../include/times.h"
#include "../include/inputs.h"

Inputs::Inputs(Times &time) {

  int iErr;

  iErr = 0;

  // ------------------------------------------------
  // Set some defaults:

  iVerbose = 3;
  euv_file="euv.csv";
  planetary_file = "orbits.csv";
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

int Inputs::read(Times &time) {

  int iErr;
  std::string line;
  std::vector<int> itime(7,0);

  iErr = 0;

  std::ifstream myFile;
  myFile.open(input_file);

  if (!myFile.is_open()) {
    std::cout << "Could not open input file: " << input_file << "!!!\n";
    iErr = 1;
  } else {

    while (getline(myFile,line)) {

      iErr = 0;

      // ---------------------------
      // #starttime
      // ---------------------------

      if (line=="#starttime") {
	for (int i=0; i<6; i++) {
	  if (getline(myFile,line)) itime[i] = stoi(line);
	  else iErr = 1;
	}
	itime[6] = 0;
	if (iErr == 0) {
	  time.set_times(itime);
	} else {
	  std::cout << "Issue in read_inputs!\n";
	  std::cout << "Should be:\n";
	  std::cout << "#starttime\n";
	  std::cout << "year     (int)\n";
	  std::cout << "month    (int)\n";
	  std::cout << "day      (int)\n";
	  std::cout << "hour     (int)\n";
	  std::cout << "min      (int)\n";
	  std::cout << "sec      (int)\n";
	}

      }

      // ---------------------------
      // #endtime
      // ---------------------------

      if (line=="#endtime") {
	for (int i=0; i<6; i++) {
	  if (getline(myFile,line)) itime[i] = stoi(line);
	  else iErr = 1;
	}
	itime[6] = 0;
	if (iErr == 0) time.set_end_time(itime);
	else {
	  std::cout << "Issue in read_inputs!\n";
	  std::cout << "Should be:\n";
	  std::cout << "#endtime\n";
	  std::cout << "year     (int)\n";
	  std::cout << "month    (int)\n";
	  std::cout << "day      (int)\n";
	  std::cout << "hour     (int)\n";
	  std::cout << "min      (int)\n";
	  std::cout << "sec      (int)\n";
	}

      }

      // ---------------------------
      // #f107file
      // ---------------------------

      if (line=="#f107file") {
	if (getline(myFile,f107_file)) {
	  DoReadF107File = 1;
	} else {
	  std::cout << "Issue in read_inputs!\n";
	  std::cout << "Should be:\n";
	  std::cout << "#f107file\n";
	  std::cout << "f107_file     (string)\n";
	}
      }

      // ---------------------------
      // #f107file
      // ---------------------------

      if (line=="#planet") {
	getline(myFile,planet);
	// This will never happen....
	if (iErr > 0) {
	  std::cout << "Issue in read_inputs!\n";
	  std::cout << "Should be:\n";
	  std::cout << "#planet\n";
	  std::cout << "planet     (string)\n";
	}
      }

    }

    myFile.close();

  }

  return iErr;


}
