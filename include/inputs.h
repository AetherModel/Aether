// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef AETHER_INCLUDE_INPUTS_H_
#define AETHER_INCLUDE_INPUTS_H_

#include <vector>
#include <string>

#include "times.h"

class Inputs {

public:

  Inputs(Times &time);
  int read(Times &time);
  int get_verbose();
  
private:

  int iVerbose;

  std::string euv_file;
  std::string input_file="aether.in";
  std::string euv_model;
  std::string planetary_file;
  std::string planet;

  // ------------------------------
  // Grid inputs:
  int IsUniformAlt;
  float alt_min;
  float dalt;

  float lat_min;
  float lat_max;
  float lon_min;
  float lon_max;
  // ------------------------------

  std::string f107_file;
  int DoReadF107File = 0;

  float euv_heating_eff_neutrals;
  float euv_heating_eff_electrons;

  std::vector<float> dt_output;
  float dt_euv;

};

#endif // AETHER_INCLUDE_INPUTS_H_
