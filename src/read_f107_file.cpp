// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iostream>

#include "aether.h"

// -----------------------------------------------------------------------------
// Read the F107 file produced by NOAA NDGC
// -----------------------------------------------------------------------------

index_file_output_struct read_f107_file(std::string f107_file,
                                        Indices indices,
                                        Report &report) {

  std::ifstream myFile;
  index_file_output_struct f107_contents;

  std::string function = "read_f107_file";
  static int iFunction = -1;
  report.enter(function, iFunction);

  f107_contents.nTimes = 0;
  f107_contents.nVars = 0;

  myFile.open(f107_file);

  if (!myFile.is_open())
    std::cout << "Could not open input file: " << f107_file << "!!!\n";

  else {

    int IsFound = 0;
    int IsAdjusted = 0;  // At some point need to take adjusted into account.
    std::string line;

    while (getline(myFile, line)) {
      if (line == "#Element: adjusted")
        IsAdjusted = 1;

      if (line == "#yyyy-MM-dd HH:mm value qualifier description") {
        IsFound = 1;
        break;
      }
    }

    // This means we found the line right before the data starts! Read
    // in the data!

    if (IsFound) {

      if (IsAdjusted && report.test_verbose(0))
        std::cout << "Need to NOT adjust F10.7, but that isn't included yet!!!"
                  << '\n';

      f107_contents.nVars = 1;
      f107_contents.var_names.push_back("F10.7");
      f107_contents.index_id.push_back(indices.get_f107_index_id());
      f107_contents.missing_values.push_back(1.0e32);

      std::string tmp;
      std::vector<int> itime(7, 0);
      std::vector<float> values;

      while (getline(myFile, line)) {
        std::stringstream ss(line);
        // year
        getline(ss, tmp, '-');
        itime[0] = stoi(tmp);
        // month
        getline(ss, tmp, '-');
        itime[1] = stoi(tmp);
        // day
        getline(ss, tmp, ' ');
        itime[2] = stoi(tmp);
        // hour
        getline(ss, tmp, ':');
        itime[3] = stoi(tmp);
        // minute
        getline(ss, tmp, ' ');
        itime[4] = stoi(tmp);
        itime[5] = 0;
        itime[6] = 0;
        f107_contents.times.push_back(time_int_to_real(itime));

        // f107
        getline(ss, tmp, '"');
        values.push_back(stof(tmp));
        f107_contents.nTimes++;
      }  // while

      // Push the vector into a vector of vectors:
      f107_contents.values.push_back(values);

    } else {
      std::cout << "Couldn't file line #yyyy-MM etc in file "
                << f107_file << "\n";
    }

    myFile.close();
  }

  report.exit(function);
  return f107_contents;
}
