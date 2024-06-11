// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_INPUTS_H_
#define INCLUDE_INPUTS_H_

#include <vector>
#include <string>

class Inputs {

public:

  Inputs() {}
  Inputs(Times &time);
  int read(Times &time);
  bool read_inputs_json(Times &time);
  bool set_verbose(json in);
  int get_verbose();
  int get_verbose_proc();
  precision_t get_dt_euv();
  bool get_include_photoelectrons();
  precision_t get_dt_report();
  precision_t get_n_outputs();
  precision_t get_dt_output(int iOutput);
  std::string get_type_output(int iOutput);
  std::string get_diffuse_auroral_model();
  std::string get_potential_model();
  std::string get_electrodynamics_dir();
  std::string get_electrodynamics_file();
  std::string get_electrodynamics_north_file();
  std::string get_electrodynamics_south_file();
  precision_t get_euv_heating_eff_neutrals();
  std::string get_euv_model();
  std::string get_euv_file();
  bool get_euv_douse();
  std::string get_aurora_file();
  std::string get_chemistry_file();
  std::string get_indices_lookup_file();
  std::vector<std::string> get_omniweb_files();
  int get_number_of_omniweb_files();
  std::string get_f107_file();
  std::string get_planet();
  std::string get_planetary_file();
  std::string get_planet_species_file();
  std::string get_collision_file();
  bool get_do_calc_bulk_ion_temp();
  precision_t get_eddy_coef();
  precision_t get_eddy_bottom();
  precision_t get_eddy_top();
  bool get_use_eddy_momentum();
  bool get_use_eddy_energy();
  std::string get_bfield_type();
  bool get_do_restart();
  std::string get_restartout_dir();
  std::string get_restartin_dir();
  precision_t get_dt_write_restarts();
  int get_original_seed();
  int get_updated_seed();
  void set_seed(int seed);
  bool write_restart();
  json get_perturb_values(); 
  bool get_do_lat_dependent_radius();
  bool get_do_J2();

  bool get_check_for_nans();
  bool get_nan_test();
  std::string get_nan_test_variable();
  
  bool get_is_cubesphere();

  bool get_NO_cooling();
  bool get_O_cooling();

  bool get_cent_acc();

  std::string get_student_name();
  bool get_is_student();
  
  json get_initial_condition_types();
  json get_boundary_condition_types();

  std::string get_advection_neutrals_vertical();

  // ------------------------------
  // Grid inputs:

  struct grid_input_struct {
    std::string alt_file;
    bool IsUniformAlt;
    bool IsDipole;
    precision_t alt_min;
    // Only really needed for Mag Field grid, since this could be overconstrained:
    precision_t alt_max;
    precision_t dalt;
    precision_t lat_min;
    precision_t lat_max;
    precision_t lon_min;
    precision_t lon_max;
    // Only needed for Mag Field grid:
    precision_t min_apex;
  };

  grid_input_struct get_grid_inputs(std::string gridtype);

  int get_nLons(std::string gridtype);
  int get_nLats(std::string gridtype);
  int get_nAlts(std::string gridtype);

  int get_nBlocksLonGeo();
  int get_nBlocksLatGeo();

  int get_nMembers();

  int iVerbose;
  int iVerboseProc;
  int iTimingDepth;

  std::string get_logfile();
  std::vector<std::string> get_species_vector();
  bool get_logfile_append();
  precision_t get_logfile_dt();

  // Satellites
  std::vector<std::string> get_satellite_files();
  std::vector<std::string> get_satellite_names();
  std::vector<precision_t> get_satellite_dts();
  
  // General get_setting functions with error checks:
  std::string get_setting_str(std::string key1);
  std::string get_setting_str(std::string key1, std::string key2);
  std::string get_setting_str(std::string key1,
                              std::string key2,
                              std::string key3);

  json get_setting_json(std::string key1);
  json get_setting_json(std::string key1, std::string key2);

  bool get_setting_bool(std::string key1);
  bool get_setting_bool(std::string key1, std::string key2);
  bool get_setting_bool(std::string key1, std::string key2, std::string key3);

  precision_t get_setting_float(std::string key1);
  precision_t get_setting_float(std::string key1, std::string key2);

  int64_t get_setting_int(std::string key1);
  int64_t get_setting_int(std::string key1, std::string key2);

  std::vector<int> get_setting_intarr(std::string key1);
  std::vector<int> get_setting_timearr(std::string key1);

  // Check settings functions:
  bool check_settings(std::string key1);
  bool check_settings(std::string key1, std::string key2);

  std::string check_settings_str(std::string key1);
  std::string check_settings_str(std::string key1, std::string key2);

  precision_t check_settings_pt(std::string key1, std::string key2);
  
  /**********************************************************************
     \brief Check to see if internal state of class is ok
   **/
  
  bool is_ok();
  
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

  bool DoRestart;
  
  precision_t dt_euv;
  precision_t dt_report;

  int nLonsGeo;
  int nLatsGeo;
  int nAltsGeo;

  int updated_seed;
  
  /// An internal variable to hold the state of the class
  bool isOk;

  std::vector<std::string> missing_settings;
};

extern Inputs input;

#endif  // INCLUDE_INPUTS_H_
