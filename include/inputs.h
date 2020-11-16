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
  std::string get_euv_model();
  std::string get_euv_file();
  std::string get_f107_file();
  std::string get_planet();
  std::string get_planetary_file();
  
  int iVerbose;

private:

  std::string euv_file = "UA/inputs/euv.csv";
  std::string input_file = "aether.in";
  std::string euv_model = "euvac";
  std::string planetary_file = "UA/inputs/orbits.csv";
  std::string planet = "earth";

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

  std::string f107_file = "";

  float euv_heating_eff_neutrals;
  float euv_heating_eff_electrons;

  std::vector<float> dt_output;
  float dt_euv;

};

#endif // AETHER_INCLUDE_INPUTS_H_
