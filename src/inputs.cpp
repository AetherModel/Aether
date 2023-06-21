// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "aether.h"

// -----------------------------------------------------------------------
// Initialize the Inputs class.  This also sets some initial values.
// The setting of initial values should probably be moved.
// -----------------------------------------------------------------------

Inputs::Inputs(Times &time, Report &report) {

  // ------------------------------------------------
  // Set some defaults:

  iVerbose = 0;
  iTimingDepth = 3;
  euv_model = "euvac";
  planet = "Earth";

  // ------------------------------------------------
  // Grid Defaults:
  geo_grid_input.alt_file = "";
  geo_grid_input.IsUniformAlt = true;
  geo_grid_input.alt_min = 100.0 * 1000.0;
  geo_grid_input.dalt = 5.0 * 1000.0;

  nLonsGeo = 12;
  nLatsGeo = 20;
  nAltsGeo = 40;

  if (nLonsGeo == 1) {
    geo_grid_input.lon_min = 0.0;
    geo_grid_input.lon_max = 0.0;
  } else {
    geo_grid_input.lon_min = 0.0;
    geo_grid_input.lon_max = 2.0 * cPI;
  }

  if (nLatsGeo == 1) {
    geo_grid_input.lat_min = 0.0;
    geo_grid_input.lat_max = 0.0;
  } else {
    geo_grid_input.lat_min = -cPI / 2;
    geo_grid_input.lat_max = cPI / 2;
  }

  euv_heating_eff_neutrals = 0.40;
  euv_heating_eff_electrons = 0.05;

  dt_output.push_back(300.0);
  type_output.push_back("states");
  dt_euv = 60.0;
  dt_report = 60.0;
  
  // ------------------------------------------------
  // Now read the input file:
  IsOk = read_inputs_json(time, report);

  if (!IsOk && iProc == 0)
    std::cout << "Error in reading input file!\n";
}

// -----------------------------------------------------------------------
// output the settings json to a file (for restart)
// -----------------------------------------------------------------------

bool Inputs::write_restart() {
  bool DidWork = true;

  if (iProc == 0) {
    std::string filename = check_settings_str("Restart", "OutDir");
    filename = filename + "/settings.json";
    DidWork = write_json(filename, settings);
  }

  DidWork = sync_across_all_procs(DidWork);
  return DidWork;
}

// -----------------------------------------------------------------------
// Settings for logfile
// -----------------------------------------------------------------------

std::string Inputs::get_logfile() {
  return check_settings_str("Logfile", "name");
}

precision_t Inputs::get_logfile_dt() {
  precision_t setting = 0;
  if(check_settings("Logfile", "dt")){
      setting = settings["Logfile", "dt"];
  }
  return setting;
}

bool Inputs::get_logfile_append() {
  bool setting = false;
  if(check_settings("Logfile", "append")){
      setting = settings["Logfile", "append"];
  }
  return setting;
}

int64_t Inputs::get_n_species() {
  int64_t setting = 0;
  if(check_settings("Logfile", "species")){
      setting = settings["Logfile", "species"].size();
  }
  return setting;
}

std::vector<std::string> Inputs::get_species_vector() {
  std::vector<std::string> species;
  std::string value = "unknown";

  for (int iOutput = 0; iOutput < get_n_species(); iOutput++) {
    value = settings.at("Logfile").at("species").at(iOutput);
    species.push_back(value);
  }

  return species;
}

// -----------------------------------------------------------------------
// Return value of a key in the json formatted inputs
// -----------------------------------------------------------------------


//set up dummy values for settings that aren't set

std::vector<int> Inputs::get_settings_intarr(std::string key1) {
  std::vector<int> value;

  if (settings.find(key1) != settings.end()) {
    int nPts = settings.at(key1).size();

    for (int i = 0; i < nPts; i++)
      value.push_back(settings.at(key1).at(i));
  } else {
    IsOk = false;
  }

  return value;
}

std::vector<int> Inputs::get_settings_timearr(std::string key1) {
  int nPtsTime = 7;
  std::vector<int> outarr(nPtsTime, 0);
  std::vector<int> timearr = get_settings_intarr(key1);

  int nPts = timearr.size();

  if (nPts > nPtsTime)
    nPts = nPtsTime;

  for (int i = 0; i < nPts; i++)
    outarr[i] = timearr[i];

  return outarr;
}

std::string Inputs::get_settings_str(std::string key1) {
  std::string value = "unknown";

  if (settings.find(key1) != settings.end())
    value = settings.at(key1);

  return value;
}

std::string Inputs::get_settings_str(std::string key1,
                                     std::string key2) {
  std::string value = "unknown";

  if (settings.find(key1) != settings.end())
    if (settings.at(key1).find(key2) != settings.at(key1).end())
      value = settings.at(key1).at(key2);

  return value;
}

// -----------------------------------------------------------------------
// Check for missing settings
// -----------------------------------------------------------------------

//dummy values, to use if the settings are not set
int dummy_int = -1;
float dummy_float = -1; 
std::string dummy_string = "unknown";

//check settings and throw invalid_argument error if the setting doesn't exist
bool Inputs::check_settings(std::string key1,
                            std::string key2) {
  //try to find the keys first
  if (settings.find(key1) != settings.end()) {
    if (settings.at(key1).find(key2) != settings.at(key1).end()) {
      return true;
    }
  }

  //if we haven't found the keys print a message & set IsOk to false
  IsOk = false;
  if(!IsOk) {
    std::cout << "Missing setting called! [" << key1 << ", " << key2 << "]\n";
  }
  return false;
}

std::string Inputs::check_settings_str(std::string key1,
                                       std::string key2) {
  if(check_settings(key1, key2))
    return settings.at(key1).at(key2);
  return dummy_string;
}

std::string Inputs::check_settings_str(std::string key) {
  if(get_settings_str(key) == dummy_string){
    IsOk = false;
    return dummy_string;
  }
  return settings[key];
}

precision_t Inputs::check_settings_pt(std::string key1,
                                      std::string key2) {
  if(check_settings(key1, key2)){
    return settings.at(key1).at(key2);
  }
  return dummy_float;
}

// -----------------------------------------------------------------------
// Return characteristics of the grid that are entered by the user
// -----------------------------------------------------------------------

Inputs::grid_input_struct Inputs::get_grid_inputs() {
  // First Get Values:
  geo_grid_input.alt_file = check_settings_str("GeoGrid", "AltFile");
  if(check_settings("GeoGrid", "IsUniformAlt"))
    geo_grid_input.IsUniformAlt = settings.at("GeoGrid").at("IsUniformAlt");
  else
    geo_grid_input.IsUniformAlt = false;
  geo_grid_input.alt_min = check_settings_pt("GeoGrid", "MinAlt");
  geo_grid_input.dalt = check_settings_pt("GeoGrid", "dAlt");
  geo_grid_input.lat_min = check_settings_pt("GeoGrid", "MinLat");
  geo_grid_input.lat_max = check_settings_pt("GeoGrid", "MaxLat");
  geo_grid_input.lon_min = check_settings_pt("GeoGrid", "MinLon");
  geo_grid_input.lon_max = check_settings_pt("GeoGrid", "MaxLon");

  // Second Change Units
  geo_grid_input.alt_min = geo_grid_input.alt_min * cKMtoM;
  geo_grid_input.lat_min = geo_grid_input.lat_min * cDtoR;
  geo_grid_input.lat_max = geo_grid_input.lat_max * cDtoR;
  geo_grid_input.lon_min = geo_grid_input.lon_min * cDtoR;
  geo_grid_input.lon_max = geo_grid_input.lon_max * cDtoR;

  // If the grid is uniform, dalt is in km, else it is in fractions of
  // scale height:
  if (geo_grid_input.IsUniformAlt)
    geo_grid_input.dalt = geo_grid_input.dalt * cKMtoM;

  return geo_grid_input;
}

// -----------------------------------------------------------------------
// Return whether user is student
// -----------------------------------------------------------------------

bool Inputs::get_is_student() {
  if(check_settings("Student", "is"))
    return settings.at("Student").at("is");
  return false;
}

// -----------------------------------------------------------------------
// Return student name
// -----------------------------------------------------------------------

std::string Inputs::get_student_name() {
  return check_settings_str("Student", "name");
}

// -----------------------------------------------------------------------
// Return whether grid is cubesphere or spherical
// -----------------------------------------------------------------------

bool Inputs::get_is_cubesphere() {
  if(check_settings("CubeSphere", "is"))
    return settings.at("CubeSphere").at("is");
  return false;
}

// -----------------------------------------------------------------------
// Return whether to restart or not
// -----------------------------------------------------------------------

bool Inputs::get_do_restart() {
  if(check_settings("Restart", "do"))
    return settings.at("Restart").at("do");
  return false;
}

// -----------------------------------------------------------------------
// Return restart OUT directory
// -----------------------------------------------------------------------

std::string Inputs::get_restartout_dir() {
  return check_settings_str("Restart", "OutDir");
}

// -----------------------------------------------------------------------
// Return restart OUT directory
// -----------------------------------------------------------------------

std::string Inputs::get_restartin_dir() {
  return check_settings_str("Restart", "InDir");
}

// -----------------------------------------------------------------------
// dt for writing restart files
// -----------------------------------------------------------------------

precision_t Inputs::get_dt_write_restarts() {
  return check_settings_pt("Restart", "dt");
}

// -----------------------------------------------------------------------
// Return magnetic field type (dipole and none defined now.)
// -----------------------------------------------------------------------

std::string Inputs::get_bfield_type() {
  return check_settings_str("BField");
}

// -----------------------------------------------------------------------
// Return the EUV model used (EUVAC only option now)
// -----------------------------------------------------------------------

std::string Inputs::get_euv_model() {
  return check_settings_str("Euv", "Model");
}

// -----------------------------------------------------------------------
// Return the heating efficiency of the neutrals for EUV
// -----------------------------------------------------------------------

precision_t Inputs::get_euv_heating_eff_neutrals() {
  return check_settings_pt("Euv", "HeatingEfficiency");
}

// -----------------------------------------------------------------------
// Return how often to calculate EUV energy deposition
// -----------------------------------------------------------------------

precision_t Inputs::get_dt_euv() {
  return check_settings_pt("Euv", "dt");
}

// -----------------------------------------------------------------------
// Return how often to report progress of simulation
// -----------------------------------------------------------------------

precision_t Inputs::get_dt_report() {
  return check_settings_pt("Debug", "dt");
}

// -----------------------------------------------------------------------
// Return number of output types
// -----------------------------------------------------------------------

precision_t Inputs::get_n_outputs() {
  return check_settings_str("Outputs", "type").size();
}

// -----------------------------------------------------------------------
// Return original random number seed
// -----------------------------------------------------------------------

int Inputs::get_original_seed() {
  if(settings.find("Seed") == settings.end()){
    IsOk = false;
    std::cout << "Error in getting seed!\n";
    return dummy_int;
  }
  return settings.at("Seed");
}

// -----------------------------------------------------------------------
// Set random number seed
// -----------------------------------------------------------------------

void Inputs::set_seed(int seed) {
  settings["Seed"] = seed;
  updated_seed = seed;
}

// -----------------------------------------------------------------------
// Return random number seed that has been updated
// -----------------------------------------------------------------------

int Inputs::get_updated_seed() {
  std::default_random_engine get_random(updated_seed);
  updated_seed = get_random();
  return updated_seed;
}

// -----------------------------------------------------------------------
// Return number of longitudes, latitudes, and altitudes in Geo grid
// -----------------------------------------------------------------------

int Inputs::get_nLonsGeo() {
  return check_settings_pt("GeoBlockSize", "nLons");
}

int Inputs::get_nLatsGeo() {
  return check_settings_pt("GeoBlockSize", "nLats");
}

int Inputs::get_nAltsGeo() {
  return check_settings_pt("GeoBlockSize", "nAlts");
}

// -----------------------------------------------------------------------
// Return number of Blocks of longitudes and latitudes in Geo grid
// -----------------------------------------------------------------------

int Inputs::get_nBlocksLonGeo() {
  return check_settings_pt("GeoBlockSize", "nBlocksLon");
}

int Inputs::get_nBlocksLatGeo() {
  return check_settings_pt("GeoBlockSize", "nBlocksLat");
}

// -----------------------------------------------------------------------
// Return number of ensemble members
// -----------------------------------------------------------------------

int Inputs::get_nMembers() {
  return check_settings_pt("Ensembles", "nMembers");
}

// -----------------------------------------------------------------------
// Return verbose variables
// -----------------------------------------------------------------------

int Inputs::get_verbose() {
  return check_settings_pt("Debug","iVerbose");
}

int Inputs::get_verbose_proc() {
  if(check_settings("Debug", "iProc"))
    return settings.at("Debug").at("iProc");
  return dummy_int;
}

// -----------------------------------------------------------------------
// Return how often to output a given output type
// -----------------------------------------------------------------------

precision_t Inputs::get_dt_output(int iOutput) {
  precision_t value = 0.0;
  int nOutputs = settings.at("Outputs").at("type").size();

  if (iOutput < nOutputs)
    value = settings.at("Outputs").at("dt").at(iOutput);

  return value;
}

// -----------------------------------------------------------------------
// Return the output type
// -----------------------------------------------------------------------

std::string Inputs::get_type_output(int iOutput) {
  std::string value = "unknown";
  int nOutputs = settings.at("Outputs").at("type").size();

  if (iOutput < nOutputs)
    value = settings.at("Outputs").at("type").at(iOutput);

  return value;
}

// -----------------------------------------------------------------------
// Return EUV file name
// -----------------------------------------------------------------------

std::string Inputs::get_euv_file() {
  return check_settings_str("Euv", "File");
}

// -----------------------------------------------------------------------
// Return aurora file name
// -----------------------------------------------------------------------

std::string Inputs::get_aurora_file() {
  return check_settings_str("AuroraFile");
}

// -----------------------------------------------------------------------
// Return Chemistry file name
// -----------------------------------------------------------------------

std::string Inputs::get_chemistry_file() {
  return check_settings_str("ChemistryFile");
}

// -----------------------------------------------------------------------
// Return Collision file name
// -----------------------------------------------------------------------

std::string Inputs::get_collision_file() {
  return check_settings_str("CollisionsFile");
}

// -----------------------------------------------------------------------
// Return Indices Lookup Filename
// -----------------------------------------------------------------------

std::string Inputs::get_indices_lookup_file() {
  return check_settings_str("IndicesLookupFile");
}

// -----------------------------------------------------------------------
// Return total number of OMNIWeb files to read
// -----------------------------------------------------------------------

int Inputs::get_number_of_omniweb_files() {
  if(settings.find("OmniwebFiles") != settings.end())
    return settings.at("OmniwebFiles").size();
  IsOk = false;
  return dummy_int;
}

// -----------------------------------------------------------------------
// Return OMNIWeb file names as a vector
// -----------------------------------------------------------------------

std::vector<std::string> Inputs::get_omniweb_files() {
  std::vector<std::string> omniweb_files;
  int nFiles = get_number_of_omniweb_files();

  for (int i = 0; i < nFiles; i++)
    omniweb_files.push_back(settings.at("OmniwebFiles").at(i));

  return omniweb_files;
}

// -----------------------------------------------------------------------
// Return F107 file to read
// -----------------------------------------------------------------------

std::string Inputs::get_f107_file() {
  return check_settings_str("F107File");
}

// -----------------------------------------------------------------------
// Return planet name
// -----------------------------------------------------------------------

std::string Inputs::get_planet() {
  return check_settings_str("Planet");
}

// -----------------------------------------------------------------------
// Return file that contains (all) planetary characteristics
// -----------------------------------------------------------------------

std::string Inputs::get_planetary_file() {
  return check_settings_str("PlanetCharacteristicsFile");
}

// -----------------------------------------------------------------------
// Return planetary file name that describes the species and such for
// a given planet
// -----------------------------------------------------------------------

std::string Inputs::get_planet_species_file() {
  return check_settings_str("PlanetSpeciesFile");
}

std::string Inputs::get_electrodynamics_file() {
  return check_settings_str("ElectrodynamicsFile");
}

// -----------------------------------------------------------------------
// Flag to do the bulk ion temperature calculation instead
// of individual ion specie temperature calculations
// -----------------------------------------------------------------------

bool Inputs::get_do_calc_bulk_ion_temp() {
  bool value = false;
  if(check_settings_str("DoCalcBulkIonTemp") != dummy_string)
    return check_settings_str("DoCalcBulkIonTemp") == "true";
  IsOk = false;
  return false;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

json Inputs::get_perturb_values() {
  json values;

  if (settings.contains("Perturb"))
    values = settings["Perturb"];
  else {
    values["Perturb"] = dummy_string;
    IsOk = false;
    std::cout << "Missing setting called! [Perturb]\n";
  }
  return values;
}

// --------------------------------------------------------------------------
// check to see if class is ok
// --------------------------------------------------------------------------

bool Inputs::is_ok() {
  return IsOk;
}