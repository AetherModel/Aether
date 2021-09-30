// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_INPUTS_H_
#define INCLUDE_INPUTS_H_

#include <vector>
#include <string>

class Inputs {

public:

  Inputs(Times &time, Report &report);
  int read(Times &time, Report &report);
  json read_json(std::string json_file);
  int read_inputs_json(Times &time, Report &report);
  int get_verbose();
  precision_t get_dt_euv();
  precision_t get_dt_report();
  precision_t get_n_outputs();
  precision_t get_dt_output(int iOutput);
  std::string get_type_output(int iOutput);
  precision_t get_euv_heating_eff_neutrals();
  std::string get_euv_model();
  std::string get_euv_file();
  std::string get_aurora_file();
  std::string get_chemistry_file();
  std::vector<std::string> get_omniweb_files();
  int get_number_of_omniweb_files();
  std::string get_f107_file();
  std::string get_planet();
  std::string get_planetary_file();
  std::string get_planet_species_file();
  std::string get_collision_file();
  std::string get_bfield_type();
  std::string get_electrodynamics_file();

  // ------------------------------
  // Grid inputs:

  struct grid_input_struct {
    std::string alt_file;
    bool IsUniformAlt;
    precision_t alt_min;
    precision_t dalt;
    precision_t lat_min;
    precision_t lat_max;
    precision_t lon_min;
    precision_t lon_max;
  };

  grid_input_struct get_grid_inputs();

  int get_nLonsGeo();
  int get_nLatsGeo();
  int get_nAltsGeo();

  int iVerbose;
  int iTimingDepth;

  std::string get_settings_str(std::string key1);
  std::string get_settings_str(std::string key1, std::string key2);
  std::vector<int> get_settings_timearr(std::string key1);
  std::vector<int> get_settings_intarr(std::string key1);
  
private:

  json settings;
  
  std::string euv_file = "UA/inputs/euv.csv";
  std::string aurora_file = "UA/inputs/aurora_earth.csv";
  std::string chemistry_file = "UA/inputs/chemistry_earth.csv";
  std::string collision_file =
    "UA/inputs/ion_neutral_collision_frequencies.csv";
  std::string input_file = "aether.in";
  std::string euv_model = "euvac";
  std::string planetary_file = "UA/inputs/orbits.csv";
  std::vector<std::string> omniweb_files;
  std::string planet = "earth";
  std::string f107_file = "";
  std::string planet_species_file = "";
  std::string electrodynamics_file = "";

  std::string bfield = "none";

  grid_input_struct geo_grid_input;

  precision_t euv_heating_eff_neutrals;
  precision_t euv_heating_eff_electrons;

  std::vector<float> dt_output;
  std::vector<std::string> type_output;
  std::string output_directory = "UA/output";
  std::string restart_out_directory = "UA/restartOut";
  std::string restart_in_directory = "UA/restartIn";

  precision_t dt_euv;
  precision_t dt_report;

  int nLonsGeo;
  int nLatsGeo;
  int nAltsGeo;
};

#endif  // INCLUDE_INPUTS_H_
