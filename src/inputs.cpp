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
    std::string filename = settings["Restart"]["OutDir"];
    filename = filename + "/settings.json";
    DidWork = write_json(filename, settings);
  }

  DidWork = sync_across_all_procs(DidWork);
  return DidWork;
}

// -----------------------------------------------------------------------
// Return value of a key in the json formatted inputs
// -----------------------------------------------------------------------

std::string Inputs::get_settings_str(std::string key1) {
  std::string value = "unknown";

  if (settings.contains(key1))
    value = settings.at(key1);

  return value;
}

std::string Inputs::get_settings_str(std::string key1,
                                     std::string key2) {
  std::string value = "unknown";

  if (settings.contains(key1))
    if (settings.at(key1).contains(key2))
      value = settings.at(key1).at(key2);

  return value;
}

std::vector<int> Inputs::get_settings_intarr(std::string key1) {
  std::vector<int> value;

  if (settings.contains(key1)) {
    int nPts = settings.at(key1).size();

    for (int i = 0; i < nPts; i++)
      value.push_back(settings.at(key1).at(i));
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

// -----------------------------------------------------------------------
// Return characteristics of the grid that are entered by the user
// -----------------------------------------------------------------------

Inputs::grid_input_struct Inputs::get_grid_inputs() {
  // First Get Values:
  geo_grid_input.alt_file = settings["GeoGrid"]["AltFile"];
  geo_grid_input.IsUniformAlt = settings["GeoGrid"]["IsUniformAlt"];
  geo_grid_input.alt_min = settings["GeoGrid"]["MinAlt"];
  geo_grid_input.dalt = settings["GeoGrid"]["dAlt"];
  geo_grid_input.lat_min = settings["GeoGrid"]["MinLat"];
  geo_grid_input.lat_max = settings["GeoGrid"]["MaxLat"];
  geo_grid_input.lon_min = settings["GeoGrid"]["MinLon"];
  geo_grid_input.lon_max = settings["GeoGrid"]["MaxLon"];

  // Second Change Units
  geo_grid_input.alt_min = geo_grid_input.alt_min * cKMtoM;
  geo_grid_input.dalt = geo_grid_input.dalt * cKMtoM;
  geo_grid_input.lat_min = geo_grid_input.lat_min * cDtoR;
  geo_grid_input.lat_max = geo_grid_input.lat_max * cDtoR;
  geo_grid_input.lon_min = geo_grid_input.lon_min * cDtoR;
  geo_grid_input.lon_max = geo_grid_input.lon_max * cDtoR;

  return geo_grid_input;
}

// -----------------------------------------------------------------------
// Return whether user is student
// -----------------------------------------------------------------------

bool Inputs::get_is_student() {
  return settings["Student"]["is"];
}

// -----------------------------------------------------------------------
// Return student name
// -----------------------------------------------------------------------

std::string Inputs::get_student_name() {
  return settings["Student"]["name"];
}

// -----------------------------------------------------------------------
// Return whether grid is cubesphere or spherical
// -----------------------------------------------------------------------

bool Inputs::get_is_cubesphere() {
  return settings["CubeSphere"]["is"];
}

// -----------------------------------------------------------------------
// Return whether to restart or not
// -----------------------------------------------------------------------

bool Inputs::get_do_restart() {
  return settings["Restart"]["do"];
}

// -----------------------------------------------------------------------
// Return restart OUT directory
// -----------------------------------------------------------------------

std::string Inputs::get_restartout_dir() {
  return settings["Restart"]["OutDir"];
}

// -----------------------------------------------------------------------
// Return restart OUT directory
// -----------------------------------------------------------------------

std::string Inputs::get_restartin_dir() {
  return settings["Restart"]["InDir"];
}

// -----------------------------------------------------------------------
// dt for writing restart files
// -----------------------------------------------------------------------

precision_t Inputs::get_dt_write_restarts() {
  return settings["Restart"]["dt"];
}

// -----------------------------------------------------------------------
// Return magnetic field type (dipole and none defined now.)
// -----------------------------------------------------------------------

std::string Inputs::get_bfield_type() {
  return settings["BField"];
}

// -----------------------------------------------------------------------
// Return the EUV model used (EUVAC only option now)
// -----------------------------------------------------------------------

std::string Inputs::get_euv_model() {
  return settings["Euv"]["Model"];
}

// -----------------------------------------------------------------------
// Return the heating efficiency of the neutrals for EUV
// -----------------------------------------------------------------------

precision_t Inputs::get_euv_heating_eff_neutrals() {
  return settings["Euv"]["HeatingEfficiency"];
}

// -----------------------------------------------------------------------
// Return how often to calculate EUV energy deposition
// -----------------------------------------------------------------------

precision_t Inputs::get_dt_euv() {
  return settings["Euv"]["dt"];
}

// -----------------------------------------------------------------------
// Return how often to report progress of simulation
// -----------------------------------------------------------------------

precision_t Inputs::get_dt_report() {
  return settings["Debug"]["dt"];
}

// -----------------------------------------------------------------------
// Return number of output types
// -----------------------------------------------------------------------

precision_t Inputs::get_n_outputs() {
  return settings["Outputs"]["type"].size();
}

// -----------------------------------------------------------------------
// Return original random number seed
// -----------------------------------------------------------------------

int Inputs::get_original_seed() {
  return settings["Seed"];
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
  return settings["GeoBlockSize"]["nLons"];
}

int Inputs::get_nLatsGeo() {
  return settings["GeoBlockSize"]["nLats"];;
}

int Inputs::get_nAltsGeo() {
  return settings["GeoBlockSize"]["nAlts"];
}

// -----------------------------------------------------------------------
// Return number of Blocks of longitudes and latitudes in Geo grid
// -----------------------------------------------------------------------

int Inputs::get_nBlocksLonGeo() {
  return settings["GeoBlockSize"]["nBlocksLon"];
}

int Inputs::get_nBlocksLatGeo() {
  return settings["GeoBlockSize"]["nBlocksLat"];;
}

// -----------------------------------------------------------------------
// Return number of ensemble members
// -----------------------------------------------------------------------

int Inputs::get_nMembers() {
  return settings["Ensembles"]["nMembers"];
}

// -----------------------------------------------------------------------
// Return verbose variables
// -----------------------------------------------------------------------

int Inputs::get_verbose() {
  return settings["Debug"]["iVerbose"];
}

int Inputs::get_verbose_proc() {
  return settings["Debug"]["iProc"];
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
  std::string value = "";
  int nOutputs = settings.at("Outputs").at("type").size();

  if (iOutput < nOutputs)
    value = settings.at("Outputs").at("type").at(iOutput);

  return value;
}

// -----------------------------------------------------------------------
// Return EUV file name
// -----------------------------------------------------------------------

std::string Inputs::get_euv_file() {
  return settings["Euv"]["File"];
}

// -----------------------------------------------------------------------
// Return aurora file name
// -----------------------------------------------------------------------

std::string Inputs::get_aurora_file() {
  return settings["AuroraFile"];
}

// -----------------------------------------------------------------------
// Return Chemistry file name
// -----------------------------------------------------------------------

std::string Inputs::get_chemistry_file() {
  return settings["ChemistryFile"];
}

// -----------------------------------------------------------------------
// Return Collision file name
// -----------------------------------------------------------------------

std::string Inputs::get_collision_file() {
  return settings["CollisionsFile"];
}

// -----------------------------------------------------------------------
// Return Indices Lookup Filename
// -----------------------------------------------------------------------

std::string Inputs::get_indices_lookup_file() {
  return settings["IndicesLookupFile"];
}

// -----------------------------------------------------------------------
// Return total number of OMNIWeb files to read
// -----------------------------------------------------------------------

int Inputs::get_number_of_omniweb_files() {
  return settings["OmniwebFiles"].size();
}

// -----------------------------------------------------------------------
// Return OMNIWeb file names as a vector
// -----------------------------------------------------------------------

std::vector<std::string> Inputs::get_omniweb_files() {
  std::vector<std::string> omniweb_files;
  int nFiles = settings["OmniwebFiles"].size();

  for (int i = 0; i < nFiles; i++)
    omniweb_files.push_back(settings.at("OmniwebFiles").at(i));

  return omniweb_files;
}

// -----------------------------------------------------------------------
// Return F107 file to read
// -----------------------------------------------------------------------

std::string Inputs::get_f107_file() {
  return settings["F107File"];
}

// -----------------------------------------------------------------------
// Return planet name
// -----------------------------------------------------------------------

std::string Inputs::get_planet() {
  return settings["Planet"];
}

// -----------------------------------------------------------------------
// Return file that contains (all) planetary characteristics
// -----------------------------------------------------------------------

std::string Inputs::get_planetary_file() {
  return settings["PlanetCharacteristicsFile"];
}

// -----------------------------------------------------------------------
// Return planetary file name that describes the species and such for
// a given planet
// -----------------------------------------------------------------------

std::string Inputs::get_planet_species_file() {
  return settings["PlanetSpeciesFile"];
}

std::string Inputs::get_electrodynamics_file() {
  return settings["ElectrodynamicsFile"];
}

// -----------------------------------------------------------------------
// Flag to do the bulk ion temperature calculation instead 
// of individual ion specie temperature calculations
// -----------------------------------------------------------------------

bool Inputs::get_do_calc_bulk_ion_temp() {
  return settings["DoCalcBulkIonTemp"];
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

json Inputs::get_perturb_values() {
  json values;

  if (settings.contains("Perturb"))
    values = settings["Perturb"];

  return values;
}

// --------------------------------------------------------------------------
// check to see if class is ok
// --------------------------------------------------------------------------

bool Inputs::is_ok() {
  return IsOk;
}
