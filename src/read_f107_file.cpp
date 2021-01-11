// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iostream>

#include "../include/time_conversion.h"

int read_f107_file(std::string f107_file,
                   std::vector<double> &time,
                   std::vector<float> &f107) {

  std::ifstream myFile;
  int iErr;

  iErr = 0;

  myFile.open(f107_file);

  if (!myFile.is_open()) {
    std::cout << "Could not open input file: " << f107_file << "!!!\n";
    iErr = 1;
  } else {

    int IsFound = 0;
    int IsAdjusted = 0;  // At some point need to take adjusted into account.
    std::string line;

    while (getline(myFile, line)) {
      if (line == "#Element: adjusted") IsAdjusted = 1;
      if (line == "#yyyy-MM-dd HH:mm value qualifier description") {
        IsFound = 1;
        break;
      }
    }

    // This means we found the line right before the data starts! Read
    // in the data!

    if (IsFound) {

      std::string tmp;
      std::vector<int> itime(7, 0);

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
        time.push_back(time_int_to_real(itime));

        // f107
        getline(ss, tmp, '"');
        f107.push_back(stof(tmp));
      }  // while

    } else {
      iErr = 1;
      std::cout << "Couldn't file line #yyyy-MM etc in file "
                << f107_file << "\n";
    }

    myFile.close();
  }

  return iErr;
}
