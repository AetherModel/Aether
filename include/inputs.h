// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_INPUTS_H_
#define INCLUDE_INPUTS_H_

#include <vector>
#include <string>

/*
  This is the input class for Aether.  It handles all of the specification of
  parameters and inputs into Aether. Some major features:
  - Aether handles a lot of inputs and outputs in json and CSV files
  - The settings variable is the main variable for controlling the settings
  - The settings variable is a json (like a python dictionary)
  - This class uses public "get" functions to access settings, which are private
  - The get functions have a naming standard that depends on the structure of the json
*/

class Inputs {

public:

  int iVerbose;
  int iVerboseProc;
  int iTimingDepth;

  // ------------------------------
  // Grid inputs:

  struct grid_input_struct {
    // Set the grid size (per block)
    // nX is nLonsPerBlock or nXsPerBlock
    // nY is nLatsPerBlock or nYsPerBlock
    // nZ is nAlts or nZs
    int64_t nX, nY, nZ;

    // min and max latitude to simulate:
    precision_t lat_min, lat_max;
    // min and max longitude to simulate:
    precision_t lon_min, lon_max;

    // The shape of the grid can be specified:
    // - "sphere", sets grid.iGridShape_ = iSphere_
    // - "cubesphere", sets grid.iGridShape_ = iCubesphere_
    // - "dipole", sets grid.iGridShape_ = iDipole_
    std::string shape;
    
    // Minimum altitude to simulate:
    precision_t alt_min;
    // Some grids allow the specification of the maximum altitude:
    precision_t alt_max;

    // Can specify the delta-altitude in either km or in scale-heights.
    // daltKm is used if there is a uniform grid
    // daltScale (0-1, typical) is used non-uniform
    precision_t daltKm, daltScale;
    // if IsUniformAlt is false, the code uses a temperature profile to
    // build an altitude profile of scale heights and uses these scale
    // heights to build the grid.
    bool IsUniformAlt;

    // Only needed for Mag Field grid:
    precision_t min_apex;

    // Some grid shapes allow specification of altitudes based on a file:
    std::string alt_file;
  };

  grid_input_struct get_grid_inputs(std::string gridtype);




  /**********************************************************************
     \brief 
     \param 
   **/
  Inputs() {}
  
  /**********************************************************************
     \brief 
     \param 
   **/
  Inputs(Times &time);

  /**********************************************************************
     \brief 
     \param 
   **/
  int read(Times &time);

  /**********************************************************************
     \brief 
     \param 
   **/
  bool read_inputs_json(Times &time);
  
  /**********************************************************************
     \brief 
     \param 
   **/
  bool set_verbose(json in);
  
  // --------------------------------------------------------------------
  // get functions:
  //  - These functions offer access to specific parts of the settings json.
  //  - They call the general functions that check whether the key(s) exists.
  //  - If the key does not exist, an error flag (in report) is set.
  //  - 
  // --------------------------------------------------------------------

  // ---------------------
  // Debugging information
  // ---------------------

  /**********************************************************************
     \brief returns settings["Debug"]["iVerbose"]
     \param none
   **/
  int get_verbose();
  
  /**********************************************************************
     \brief returns settings["Debug"]["iProc"]
     \param none
   **/
  int get_verbose_proc();
  
  /**********************************************************************
     \brief returns settings["Debug"]["dt"]
     \param none
   **/
  precision_t get_dt_report();

  // ---------------------
  // Output information
  // ---------------------


  /**********************************************************************
     \brief returns the size of settings["Outputs"]["type"]
     \param none
   **/
  precision_t get_n_outputs();
  
  /**********************************************************************
     \brief returns settings["Outputs"]["dt"][iOutput]
     \param iOutput int specifying which output file type to report on
   **/
  precision_t get_dt_output(int iOutput);
  
  /**********************************************************************
     \brief returns settings["Outputs"]["type"][iOutput]
     \param iOutput int specifying which output file type to report on
   **/
  std::string get_type_output(int iOutput);
  



  /**********************************************************************
     \brief returns settings["Euv"]["dt"]
     \param none
   **/
  precision_t get_dt_euv();
  
  /**********************************************************************
     \brief returns settings["Euv"]["IncludePhotoElectrons"]
     \param none
   **/
  bool get_include_photoelectrons();


  /**********************************************************************
     \brief returns settings["Electrodynamics"]["DiffuseAurora"]
     \param none
   **/
  std::string get_diffuse_auroral_model();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::string get_potential_model();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::string get_electrodynamics_dir();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::string get_electrodynamics_file();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::string get_electrodynamics_north_file();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::string get_electrodynamics_south_file();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  precision_t get_euv_heating_eff_neutrals();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::string get_euv_model();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::string get_euv_file();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  bool get_euv_douse();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::string get_aurora_file();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::string get_chemistry_file();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::string get_indices_lookup_file();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::vector<std::string> get_omniweb_files();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  int get_number_of_omniweb_files();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::string get_f107_file();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::string get_planet();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::string get_planetary_file();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::string get_planet_species_file();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::string get_collision_file();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  bool get_do_calc_bulk_ion_temp();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  precision_t get_eddy_coef();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  precision_t get_eddy_bottom();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  precision_t get_eddy_top();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  bool get_use_eddy_momentum();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  bool get_use_eddy_energy();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::string get_bfield_type();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  bool get_do_restart();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::string get_restartout_dir();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::string get_restartin_dir();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  precision_t get_dt_write_restarts();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  int get_original_seed();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  int get_updated_seed();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  void set_seed(int seed);
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  bool write_restart();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  json get_perturb_values(); 
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  bool get_do_lat_dependent_radius();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  bool get_do_J2();

  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  bool get_check_for_nans();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  bool get_nan_test();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::string get_nan_test_variable();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  bool get_is_cubesphere();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  bool get_NO_cooling();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  bool get_O_cooling();

  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  bool get_cent_acc();

  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::string get_student_name();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  bool get_is_student();
  
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  json get_initial_condition_types();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  json get_boundary_condition_types();

  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::string get_advection_neutrals_vertical();

  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  int get_nLons(std::string gridtype);
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  int get_nLats(std::string gridtype);
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  int get_nAlts(std::string gridtype);

  /**********************************************************************
     \brief returns settings[gridtype, "shape"]
     \param 
   **/
  std::string get_grid_shape(std::string gridtype);
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  int get_nMembers();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::string get_logfile();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::string get_logfile(int64_t iLog);
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::vector<std::string> get_species_vector();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  bool get_logfile_append();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  precision_t get_logfile_dt();

  // Satellites
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::vector<std::string> get_satellite_files();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::vector<std::string> get_satellite_names();
  
  /**********************************************************************
     \brief returns settings["
     \param 
   **/
  std::vector<precision_t> get_satellite_dts();
  
  // General get_setting functions with error checks:
  
  /**********************************************************************
     \brief 
     \param 
   **/
  std::string get_setting_str(std::string key1);
  
  /**********************************************************************
     \brief 
     \param 
   **/
  std::string get_setting_str(std::string key1, std::string key2);
  
  /**********************************************************************
     \brief 
     \param 
   **/
  std::string get_setting_str(std::string key1,
                              std::string key2,
                              std::string key3);

  
  /**********************************************************************
     \brief 
     \param 
   **/
  json get_setting_json(std::string key1);
  
  /**********************************************************************
     \brief 
     \param 
   **/
  json get_setting_json(std::string key1, std::string key2);

  
  /**********************************************************************
     \brief 
     \param 
   **/
  bool get_setting_bool(std::string key1);
  
  /**********************************************************************
     \brief 
     \param 
   **/
  bool get_setting_bool(std::string key1, std::string key2);
  
  /**********************************************************************
     \brief 
     \param 
   **/
  bool get_setting_bool(std::string key1, std::string key2, std::string key3);

  
  /**********************************************************************
     \brief 
     \param 
   **/
  precision_t get_setting_float(std::string key1);
  
  /**********************************************************************
     \brief 
     \param 
   **/
  precision_t get_setting_float(std::string key1, std::string key2);

  
  /**********************************************************************
     \brief 
     \param 
   **/
  int64_t get_setting_int(std::string key1);
  
  /**********************************************************************
     \brief 
     \param 
   **/
  int64_t get_setting_int(std::string key1, std::string key2);

  
  /**********************************************************************
     \brief 
     \param 
   **/
  std::vector<int> get_setting_intarr(std::string key1);
  
  /**********************************************************************
     \brief 
     \param 
   **/
  std::vector<int> get_setting_intarr(std::string key1, std::string key2);

  /**********************************************************************
     \brief 
     \param 
   **/
  std::vector<int> get_setting_timearr(std::string key1);

  // Check settings functions:
  
  /**********************************************************************
     \brief 
     \param 
   **/
  bool check_settings(std::string key1);
  
  /**********************************************************************
     \brief 
     \param 
   **/
  bool check_settings(std::string key1, std::string key2);

  
  /**********************************************************************
     \brief 
     \param 
   **/
  std::string check_settings_str(std::string key1);
  
  /**********************************************************************
     \brief 
     \param 
   **/
  std::string check_settings_str(std::string key1, std::string key2);

  
  /**********************************************************************
     \brief 
     \param 
   **/
  precision_t check_settings_pt(std::string key1, std::string key2);
  
  /**********************************************************************
     \brief Check to see if internal state of class is ok
   **/
  bool is_ok();
  
private:

  // This is the main variable that contains all of the settings in Aether:
  json settings;
  
  // These are a bunch of misc strings that should go away:
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

// This makes the input variable a global variable, so it can be used
// without passing it to every function.  This is done because the
// input class is needed by almost all functions within Aether
extern Inputs input;

#endif  // INCLUDE_INPUTS_H_
