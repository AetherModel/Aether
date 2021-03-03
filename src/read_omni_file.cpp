// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iostream>

#include "aether.h"

index_file_output_struct read_omni_file(std::string omni_file,
					Indices indices,
					Report &report) {

  std::ifstream myFile;
  int iErr;

  index_file_output_struct file_contents;

  std::string function = "read_omni_file";
  static int iFunction = -1;
  report.enter(function, iFunction);

  file_contents.nTimes = 0;
  file_contents.nVars = 0;
  
  iErr = 0;

  myFile.open(omni_file);

  if (!myFile.is_open()) {
    std::cout << "Could not open input file: " << omni_file << "!!!\n";
    iErr = 1;
  } else {

    int IsFound = 0;
    std::string line;

    std::size_t test;
    
    // Read in the header. We can ignore this:
    while (getline(myFile, line)) {
      test = line.find("Selected parameters:");
      if (test != std::string::npos) {
        IsFound = 1;
        break;
      }
    }

    // Read in the variable names. Have to match them to indices:
    if (IsFound) {

      while (getline(myFile, line)) {
	if (line.length() < 2) break;
	// The first two characters are the variable number, third is space:
	file_contents.var_names.push_back(line.substr(3));
	file_contents.nVars++;
      }

      // Look up some characteristics of the variables by their names.
      // This is specific to the OMNIWEB data, so we have some special
      // functions that return the missing values for the variables
      // and the index id (within the Indices class).
      std::vector<float> values_dummy;
      float single_value;
      for (int iVar = 0; iVar < file_contents.nVars; iVar++) {
	std::string var = file_contents.var_names[iVar];
	int index_id = pair_omniweb_index(var, indices);
	float missing_value = get_omniweb_missing_value(var);
	file_contents.index_id.push_back(index_id);
	file_contents.missing_values.push_back(missing_value);
	file_contents.values.push_back(values_dummy);
      }

      std::vector<std::vector<std::string>> values_string = read_ssv(myFile);
      std::vector<int> iTime(7,0);
      
      // The last thing we read was a blank line.  The next line should be
      // a header for the numbers to come.  We need to read this and
      // figure out if we have hourly values or minute values.  The only
      // reason this matters is because we need to know which column is
      // the first variable, which is determined by whether there is a "MN"
      // in this line.

      int iMinuteData = 0;
      if (values_string[0][3] == "MN") iMinuteData = 1;

      int64_t nTimes_tmp = values_string.size();
      for (int64_t iLine = 1; iLine < nTimes_tmp; iLine++) {
	if (values_string[iLine].size() > 4) {
	  // Year:
	  iTime[0] = stoi(values_string[iLine][0]);
	  // This is Day of Year, which can be set by saying that it is in
	  // January:
	  iTime[1] = 1;
	  // Day of year:
	  iTime[2] = stoi(values_string[iLine][1]);
	  // Hour:
	  iTime[3] = stoi(values_string[iLine][2]);
	  if (iMinuteData) {
	    iTime[4] = stoi(values_string[iLine][3]); // minute
	  } else {
	    iTime[4] = 0;
	  }
	  iTime[5] = 0;
	  iTime[6] = 0;
	  file_contents.times.push_back(time_int_to_real(iTime));
	  file_contents.nTimes++;
	  for (int iVar = 0; iVar < file_contents.nVars; iVar++) {
	    single_value = stof(values_string[iLine][3+iMinuteData+iVar]);
	    file_contents.values[iVar].push_back(single_value);
	  }
	} else {
	  break;
	}
      }
    } // IsFound

    myFile.close();
  }

  report.exit(function);
  return file_contents;
}

// ----------------------------------------------------------------------
//
// ----------------------------------------------------------------------

int pair_omniweb_index(std::string var_name, Indices indices) {

  int index = -1;
  std::size_t test;

  test = var_name.find("BX");
  if (test != std::string::npos) {
    index = indices.get_imf_bx_index_id();
  }

  test = var_name.find("BY");
  if (test != std::string::npos) {
    index = indices.get_imf_by_index_id();
  }

  test = var_name.find("BZ");
  if (test != std::string::npos) {
    index = indices.get_imf_bz_index_id();
  }

  test = var_name.find("Vx");
  if (test != std::string::npos) {
    index = indices.get_sw_vx_index_id();
  }

  test = var_name.find("Vy");
  if (test != std::string::npos) {
    index = indices.get_sw_vy_index_id();
  }

  test = var_name.find("Vz");
  if (test != std::string::npos) {
    index = indices.get_sw_vz_index_id();
  }

  test = var_name.find("Proton Density");
  if (test != std::string::npos) {
    index = indices.get_sw_n_index_id();
  }

  test = var_name.find("Temperature");
  if (test != std::string::npos) {
    index = indices.get_sw_t_index_id();
  }

  test = var_name.find("AE-index");
  if (test != std::string::npos) {
    index = indices.get_ae_index_id();
  }

  test = var_name.find("AL-index");
  if (test != std::string::npos) {
    index = indices.get_al_index_id();
  }

  test = var_name.find("AU-index");
  if (test != std::string::npos) {
    index = indices.get_au_index_id();
  }

  return index;
  
}

// ----------------------------------------------------------------------
//
// ----------------------------------------------------------------------

float get_omniweb_missing_value(std::string var_name) {

  float missing = -1.0e32;
  std::size_t test;

  test = var_name.find("BX");
  if (test != std::string::npos) {
    missing = 9999.99;
  }

  test = var_name.find("BY");
  if (test != std::string::npos) {
    missing = 9999.99;
  }

  test = var_name.find("BZ");
  if (test != std::string::npos) {
    missing = 9999.99;
  }

  test = var_name.find("Vx");
  if (test != std::string::npos) {
    missing = 99999.9;
  }

  test = var_name.find("Vy");
  if (test != std::string::npos) {
    missing = 99999.9;
  }

  test = var_name.find("Vz");
  if (test != std::string::npos) {
    missing = 99999.9;
  }

  test = var_name.find("Proton Density");
  if (test != std::string::npos) {
    missing = 999.99;
  }

  test = var_name.find("Temperature");
  if (test != std::string::npos) {
    missing = 999999.0;
  }

  test = var_name.find("AE-index");
  if (test != std::string::npos) {
    missing = 9999.99;
  }

  test = var_name.find("AL-index");
  if (test != std::string::npos) {
    missing = 9999.99;
  }

  test = var_name.find("AU-index");
  if (test != std::string::npos) {
    missing = 9999.99;
  }

  return missing;
  
}
