// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "aether.h"

Inputs input;

// -----------------------------------------------------------------------
// Initialize the Inputs class.  This also sets some initial values.
// The setting of initial values should probably be moved.
// -----------------------------------------------------------------------

Inputs::Inputs(Times &time)
{

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
  geo_grid_input.alt_min = 100.0;
  geo_grid_input.daltKm = 5.0;

  nLonsGeo = 12;
  nLatsGeo = 20;
  nAltsGeo = 40;

  if (nLonsGeo == 1)
  {
    geo_grid_input.lon_min = 0.0;
    geo_grid_input.lon_max = 0.0;
  }
  else
  {
    geo_grid_input.lon_min = 0.0;
    geo_grid_input.lon_max = 2.0 * cPI;
  }

  if (nLatsGeo == 1)
  {
    geo_grid_input.lat_min = 0.0;
    geo_grid_input.lat_max = 0.0;
  }
  else
  {
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
  isOk = read_inputs_json(time);

  if (report.test_verbose(1))
  {
    std::cout << "Settings read in:\n";
    std::cout << std::setw(2) << settings;
  }

  if (!isOk && iProc == 0)
    std::cout << "Error in reading input file!\n";
}

// -----------------------------------------------------------------------
// output the settings json to a file (for restart)
// -----------------------------------------------------------------------

bool Inputs::write_restart()
{
  bool didWork = true;

  if (iProc == 0)
  {
    std::string filename = get_setting_str("Restart", "OutDir");
    filename = filename + "/settings.json";
    didWork = write_json(filename, settings);
  }

  didWork = sync_across_all_procs(didWork);
  return didWork;
}

// -----------------------------------------------------------------------
// General check functions to see if keys exist:
// check settings and throw invalid_argument error
// if the setting doesn't exist
// -----------------------------------------------------------------------

// dummy values, to use if the settings are not set
int dummy_int = -1;
float dummy_float = -1;
std::string dummy_string = "unknown";

// -----------------------------------------------------------------------
// 2 keys:

bool Inputs::check_settings(std::string key1,
                            std::string key2)
{
  if (report.test_verbose(2))
    std::cout << "checking setting : "
              << key1 << " and "
              << key2 << "\n";

  // try to find the keys first
  if (settings.find(key1) != settings.end())
  {
    if (settings.at(key1).find(key2) != settings.at(key1).end())
      isOk = true;
  }
  else
    // if we haven't found the keys print a message & set IsOk to false
    isOk = false;

  if (!isOk)
  {
    report.error("Error in setting : " + key1 + " : " + key2);
    std::cout << "Missing setting called! [" << key1 << ", " << key2 << "]\n";
  }

  return isOk;
}

// -----------------------------------------------------------------------
// 1 key:

bool Inputs::check_settings(std::string key1)
{
  if (report.test_verbose(2))
    std::cout << "checking setting : " << key1 << "\n";

  // try to find the keys first
  if (settings.find(key1) != settings.end())
    isOk = true;
  else
    // if we haven't found the key print a message & set IsOk to false
    isOk = false;

  // perturb is non-essential, otherwise print error message
  if (!isOk && key1 != "Perturb")
  {
    report.error("Error in setting : " + key1);
    std::cout << "Missing setting called! [" << key1 << "]\n";
  }

  return isOk;
}

// -----------------------------------------------------------------------
// Functions that check keys and return values:
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
// a general int vector

std::vector<int> Inputs::get_setting_intarr(std::string key1)
{
  std::vector<int> value;

  if (check_settings(key1))
  {
    int nPts = settings.at(key1).size();
    isOk = true;

    for (int i = 0; i < nPts; i++)
      value.push_back(settings.at(key1).at(i));
  }
  else
    isOk = false;

  return value;
}

std::vector<int> Inputs::get_setting_intarr(std::string key1,
                                            std::string key2)
{
  std::vector<int> value;

  if (check_settings(key1, key2))
  {
    int nPts = settings.at(key1).at(key2).size();
    isOk = true;

    for (int i = 0; i < nPts; i++)
      value.push_back(settings.at(key1).at(key2).at(i));
  }
  else
    isOk = false;

  return value;
}

// -----------------------------------------------------------------------
// A specific length int vector

std::vector<int> Inputs::get_setting_timearr(std::string key1)
{
  int nPtsTime = 7;
  std::vector<int> outarr(nPtsTime, 0);
  std::vector<int> timearr = get_setting_intarr(key1);

  if (isOk)
  {
    int nPts = timearr.size();

    if (nPts > nPtsTime)
      nPts = nPtsTime;

    for (int i = 0; i < nPts; i++)
      outarr[i] = timearr[i];
  }

  return outarr;
}

// -----------------------------------------------------------------------
// a string with 1 key:

std::string Inputs::get_setting_str(std::string key1)
{
  std::string value = "unknown";

  if (check_settings(key1))
    value = settings.at(key1);

  return value;
}

// -----------------------------------------------------------------------
// a string with 2 keys:

std::string Inputs::get_setting_str(std::string key1,
                                    std::string key2)
{
  std::string value = "unknown";

  if (check_settings(key1, key2))
    value = settings.at(key1).at(key2);

  return value;
}

// -----------------------------------------------------------------------
// a string with 3 keys:

std::string Inputs::get_setting_str(std::string key1,
                                    std::string key2,
                                    std::string key3)
{
  std::string value = "unknown";
  isOk = false;

  if (settings.find(key1) != settings.end())
    if (settings.at(key1).find(key2) != settings.at(key1).end())
      if (settings.at(key1).at(key2).find(key3) !=
          settings.at(key1).at(key2).end())
      {
        value = settings.at(key1).at(key2).at(key3);
        isOk = true;
      }

  if (!isOk)
    report.error("Error in setting : " + key1 + " : " + key2 + " : " + key3);

  return value;
}

// -----------------------------------------------------------------------
// an int with 1 key:

int64_t Inputs::get_setting_int(std::string key1)
{
  int64_t value = LONG_MIN;

  if (check_settings(key1))
    value = settings.at(key1);

  return value;
}

// -----------------------------------------------------------------------
// an int with 2 keys:

int64_t Inputs::get_setting_int(std::string key1,
                                std::string key2)
{
  int64_t value = LONG_MIN;

  if (check_settings(key1, key2))
    value = settings.at(key1).at(key2);

  return value;
}

// -----------------------------------------------------------------------
// a bool with 1 key:

bool Inputs::get_setting_bool(std::string key1)
{
  bool value = false;

  if (check_settings(key1))
    value = settings.at(key1);

  return value;
}

// -----------------------------------------------------------------------
// a bool with 2 keys:

bool Inputs::get_setting_bool(std::string key1,
                              std::string key2)
{
  bool value = false;

  if (check_settings(key1, key2))
    value = settings.at(key1).at(key2);

  return value;
}

// -----------------------------------------------------------------------
// a bool with 3 keys:

bool Inputs::get_setting_bool(std::string key1,
                              std::string key2,
                              std::string key3)
{
  bool value = false;
  isOk = false;

  if (settings.find(key1) != settings.end())
    if (settings.at(key1).find(key2) != settings.at(key1).end())
      if (settings.at(key1).at(key2).find(key3) !=
          settings.at(key1).at(key2).end())
      {
        value = settings.at(key1).at(key2).at(key3);
        isOk = true;
      }

  if (!isOk)
    report.error("Error in setting : " + key1 + " : " + key2 + " : " + key3);

  return value;
}

// -----------------------------------------------------------------------
// a float with 1 key:

precision_t Inputs::get_setting_float(std::string key1)
{
  precision_t value = std::numeric_limits<precision_t>::lowest();

  if (check_settings(key1))
    value = settings.at(key1);

  return value;
}

// -----------------------------------------------------------------------
// a float with 2 key:

precision_t Inputs::get_setting_float(std::string key1,
                                      std::string key2)
{
  precision_t value = std::numeric_limits<precision_t>::lowest();

  if (check_settings(key1, key2))
    value = settings.at(key1).at(key2);

  return value;
}

// -----------------------------------------------------------------------
// a json with 1 key:

json Inputs::get_setting_json(std::string key1)
{
  json value;

  if (settings.find(key1) != settings.end())
    value = settings.at(key1);
  else
  {
    isOk = false;
    report.error("Error in setting : " + key1);
  }

  return value;
}

// -----------------------------------------------------------------------
// a json with 2 keys:

json Inputs::get_setting_json(std::string key1,
                              std::string key2)
{
  json value;

  if (settings.find(key1) != settings.end())
    if (settings.at(key1).find(key2) != settings.at(key1).end())
      value = settings.at(key1).at(key2);
    else
    {
      isOk = false;
      report.error("Error in setting : " + key1 + " : " + key2);
    }
  else
  {
    isOk = false;
    report.error("Error in setting : " + key1);
  }

  return value;
}

// -----------------------------------------------------------------------
// a string with 2 keys:

std::string Inputs::check_settings_str(std::string key1,
                                       std::string key2)
{
  if (check_settings(key1, key2))
    return settings.at(key1).at(key2);

  return dummy_string;
}

// -----------------------------------------------------------------------
// a string with 1 key:

std::string Inputs::check_settings_str(std::string key)
{
  if (get_setting_str(key) == dummy_string)
  {
    isOk = false;
    return dummy_string;
  }

  return settings[key];
}

// -----------------------------------------------------------------------
// a float with 2 keys:

precision_t Inputs::check_settings_pt(std::string key1,
                                      std::string key2)
{
  if (check_settings(key1, key2))
    return settings.at(key1).at(key2);

  isOk = false;
  return dummy_float;
}

// -----------------------------------------------------------------------
// Return characteristics of the grid that are entered by the user
// gridtype needs to be "neuGrid" or "ionGrid"
// -----------------------------------------------------------------------

Inputs::grid_input_struct Inputs::get_grid_inputs(std::string gridtype)
{

  Inputs::grid_input_struct grid_specs;

  std::vector<int> min_max;

  grid_specs.shape = check_settings_str(gridtype, "Shape");
  grid_specs.nX = get_setting_int(gridtype, "nLonsPerBlock");
  grid_specs.nY = get_setting_int(gridtype, "nLatsPerBlock");
  grid_specs.nZ = get_setting_int(gridtype, "nAlts");

  min_max = get_setting_intarr(gridtype, "LonRange");
  grid_specs.lon_min = min_max[0] * cDtoR;
  grid_specs.lon_max = min_max[1] * cDtoR;

  grid_specs.alt_min = check_settings_pt(gridtype, "MinAlt");

  // The rest of the settings are different for mag/geo grids,
  // First take the magnetic options, then "else" should be (cube-)sphere
  if (grid_specs.shape == "dipole")
  {
    // min_apex MUST be more than min_alt:
    grid_specs.min_apex = check_settings_pt(gridtype, "MinApex");
    if (grid_specs.min_apex <= grid_specs.alt_min)
    {
      report.error("Error in Inputs! min_apex must be more than min_alt!");
    }
    grid_specs.LatStretch = check_settings_pt(gridtype, "LatStretch");
    grid_specs.max_lat_dipole = check_settings_pt(gridtype, "LatMax") * cDtoR;
    grid_specs.FieldLineStretch = check_settings_pt(gridtype, "LineSpacing");
  }
  else
  {
    min_max = get_setting_intarr(gridtype, "LatRange");
    grid_specs.lat_min = min_max[0] * cDtoR;
    grid_specs.lat_max = min_max[1] * cDtoR;
    grid_specs.alt_file = check_settings_str(gridtype, "AltFile");
    grid_specs.IsUniformAlt = get_setting_bool(gridtype, "IsUniformAlt");
    if (grid_specs.IsUniformAlt)
      grid_specs.daltKm = check_settings_pt(gridtype, "dAltkm");
    else
      grid_specs.daltScale = check_settings_pt(gridtype, "dAltScale");
  }

  return grid_specs;
}

// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
// More complicated get functions:
// -----------------------------------------------------------------------
// -----------------------------------------------------------------------

// This is needed, because we may want to check for verbose specifically
// in a given json and not the normal settings json:

bool Inputs::set_verbose(json in)
{
  bool didWork = true;
  int iVerbose = -1;

  // Want to set verbose level ASAP:
  if (in.contains("Debug"))
  {
    if (in.at("Debug").contains("iVerbose"))
    {
      iVerbose = in.at("Debug").at("iVerbose");

      if (in.at("Debug").contains("iProc"))
      {
        if (iProc != in.at("Debug").at("iProc"))
          iVerbose = -1;
      }
    }
    else
      didWork = false;
  }
  else
    didWork = false;

  if (iVerbose > 0)
  {
    std::cout << "Setting iVerbose : " << iVerbose << "\n";
    report.set_verbose(iVerbose);
  }

  return didWork;
}

// -----------------------------------------------------------------------
// Return total number of OMNIWeb files to read
// -----------------------------------------------------------------------

int Inputs::get_number_of_omniweb_files()
{
  if (settings.find("OmniwebFiles") != settings.end())
    return settings.at("OmniwebFiles").size();

  isOk = false;
  return dummy_int;
}

// -----------------------------------------------------------------------
// Return OMNIWeb file names as a vector
// -----------------------------------------------------------------------

std::vector<std::string> Inputs::get_omniweb_files()
{
  std::vector<std::string> omniweb_files;
  int nFiles = get_number_of_omniweb_files();

  for (int i = 0; i < nFiles; i++)
    omniweb_files.push_back(settings.at("OmniwebFiles").at(i));

  return omniweb_files;
}

// -----------------------------------------------------------------------
// Return how often to output a given output type
// -----------------------------------------------------------------------

precision_t Inputs::get_dt_output(int iOutput)
{
  precision_t value = 0.0;
  int nOutputs = settings.at("Outputs").at("type").size();

  if (iOutput < nOutputs)
    value = settings.at("Outputs").at("dt").at(iOutput);

  return value;
}

// -----------------------------------------------------------------------
// Return the output type
// -----------------------------------------------------------------------

std::string Inputs::get_type_output(int iOutput)
{
  std::string value = "unknown";
  int nOutputs = settings.at("Outputs").at("type").size();

  if (iOutput < nOutputs)
    value = settings.at("Outputs").at("type").at(iOutput);

  return value;
}

// -----------------------------------------------------------------------
// Set random number seed
// -----------------------------------------------------------------------

void Inputs::set_seed(int seed)
{
  settings["Seed"] = seed;
  updated_seed = seed;
}

// -----------------------------------------------------------------------
// Return random number seed that has been updated
// -----------------------------------------------------------------------

int Inputs::get_updated_seed()
{
  std::default_random_engine get_random(updated_seed);
  updated_seed = get_random();
  return updated_seed;
}

// -----------------------------------------------------------------------
// Return log file name
// -----------------------------------------------------------------------

std::string Inputs::get_logfile()
{
  std::string logfile = get_setting_str("Logfile", "name");

  if (nMembers > 1)
    logfile = add_cmember(logfile);

  return logfile;
}

// -----------------------------------------------------------------------
// Return log file name
// -----------------------------------------------------------------------

std::string Inputs::get_logfile(int64_t iLog)
{
  std::string logfile = "log.txt";
  if (check_settings("Logfile", "name"))
  {
    int64_t nLogs = settings.at("Logfile").at("name").size();
    if (nLogs == 1)
    {
      logfile = settings.at("Logfile").at("name").at(iLog);
      // logfile = get_setting_str("Logfile", "name");
    }
    else
    {
      if (iLog > nLogs - 1)
      {
        report.error("Error in getting logfile name!");
        logfile = settings.at("Logfile").at("name").at(nLogs - 1);
      }
      else
      {
        logfile = settings.at("Logfile").at("name").at(iLog);
      }
    }
  }
  if (nMembers > 1)
    logfile = add_cmember(logfile);

  return logfile;
}

// -----------------------------------------------------------------------
// Return the name of specified variables as a vector
// -----------------------------------------------------------------------

std::vector<std::string> Inputs::get_species_vector()
{
  std::vector<std::string> species;
  const json &json_species = get_setting_json("Logfile", "species");

  for (size_t iOutput = 0; iOutput < json_species.size(); iOutput++)
    species.push_back(json_species.at(iOutput));

  return species;
}

// -----------------------------------------------------------------------
// Return the name of satellite files as a vector
// -----------------------------------------------------------------------

std::vector<std::string> Inputs::get_satellite_files()
{
  std::vector<std::string> files;
  const json &json_files = get_setting_json("Satellites", "files");

  for (size_t i = 0; i < json_files.size(); ++i)
    files.push_back(json_files.at(i));

  return files;
}

// -----------------------------------------------------------------------
// Return the output file names of satellites as a vector
// -----------------------------------------------------------------------

std::vector<std::string> Inputs::get_satellite_names()
{
  std::vector<std::string> names;
  const json &json_names = get_setting_json("Satellites", "names");

  for (size_t i = 0; i < json_names.size(); ++i)
    names.push_back(json_names.at(i));

  return names;
}

// -----------------------------------------------------------------------
// Return how oftern to write log file for satellites as a vector
// -----------------------------------------------------------------------

std::vector<precision_t> Inputs::get_satellite_dts()
{
  std::vector<precision_t> dts;
  const json &json_dts = get_setting_json("Satellites", "dts");

  for (size_t i = 0; i < json_dts.size(); ++i)
    dts.push_back(json_dts.at(i));

  return dts;
}

// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
// Extremely simple get functions:
// -----------------------------------------------------------------------
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
// Return how oftern to write log file
// -----------------------------------------------------------------------

precision_t Inputs::get_logfile_dt()
{
  return get_setting_float("Logfile", "dt");
}

// -----------------------------------------------------------------------
// Return whether to append or rewrite
// -----------------------------------------------------------------------

bool Inputs::get_logfile_append()
{
  return get_setting_bool("Logfile", "append");
}

// -----------------------------------------------------------------------
// Return whether user is student
// -----------------------------------------------------------------------

bool Inputs::get_is_student()
{
  return get_setting_bool("Student", "is");
}

// -----------------------------------------------------------------------
// Return student name
// -----------------------------------------------------------------------

std::string Inputs::get_student_name()
{
  return check_settings_str("Student", "name");
}

// -----------------------------------------------------------------------
// Return whether grid is cubesphere or spherical
// -----------------------------------------------------------------------

bool Inputs::get_is_cubesphere()
{
  return get_setting_bool("CubeSphere", "is");
}

// -----------------------------------------------------------------------
// Return whether to restart or not
// -----------------------------------------------------------------------

bool Inputs::get_do_restart()
{
  return get_setting_bool("Restart", "do");
}

// -----------------------------------------------------------------------
// Return NO cooling
// -----------------------------------------------------------------------

bool Inputs::get_NO_cooling()
{
  return get_setting_bool("Sources", "Neutrals", "NO_cool");
}

// -----------------------------------------------------------------------
// Return O cooling
// -----------------------------------------------------------------------

bool Inputs::get_O_cooling()
{
  return get_setting_bool("Sources", "Neutrals", "O_cool");
}

// -----------------------------------------------------------------------
// Return centripetal acceleration
// -----------------------------------------------------------------------

bool Inputs::get_cent_acc()
{
  return get_setting_bool("Sources", "Grid", "Cent_acc");
}

// -----------------------------------------------------------------------
// Return restart OUT directory
// -----------------------------------------------------------------------

std::string Inputs::get_restartout_dir()
{
  return check_settings_str("Restart", "OutDir");
}

// -----------------------------------------------------------------------
// Return restart In directory
// -----------------------------------------------------------------------

std::string Inputs::get_restartin_dir()
{
  return check_settings_str("Restart", "InDir");
}

// -----------------------------------------------------------------------
// dt for writing restart files
// -----------------------------------------------------------------------

precision_t Inputs::get_dt_write_restarts()
{
  return check_settings_pt("Restart", "dt");
}

// -----------------------------------------------------------------------
// Return magnetic field type (dipole and none defined now.)
// -----------------------------------------------------------------------

std::string Inputs::get_bfield_type()
{
  return check_settings_str("BField");
}

// -----------------------------------------------------------------------
// Return whether to use EUV at all
// -----------------------------------------------------------------------

bool Inputs::get_euv_douse()
{
  return get_setting_bool("Euv", "doUse");
}

// -----------------------------------------------------------------------
// Return the Electrodynamics Dir - this is where all of the
//   files that are for the empirical models reside
// -----------------------------------------------------------------------

std::string Inputs::get_electrodynamics_north_file()
{
  return check_settings_str("Electrodynamics", "NorthFile");
}

// -----------------------------------------------------------------------
// Return the Electrodynamics Dir - this is where all of the
//   files that are for the empirical models reside
// -----------------------------------------------------------------------

std::string Inputs::get_electrodynamics_south_file()
{
  return check_settings_str("Electrodynamics", "SouthFile");
}

// -----------------------------------------------------------------------
// Return the Electrodynamics Dir - this is where all of the
//   files that are for the empirical models reside
// -----------------------------------------------------------------------

std::string Inputs::get_electrodynamics_file()
{
  return check_settings_str("Electrodynamics", "File");
}

// -----------------------------------------------------------------------
// Return the Electrodynamics Dir - this is where all of the
//   files that are for the empirical models reside
// -----------------------------------------------------------------------

std::string Inputs::get_electrodynamics_dir()
{
  return check_settings_str("Electrodynamics", "Dir");
}

// -----------------------------------------------------------------------
// Return the Electrodynamics Potential Model
// -----------------------------------------------------------------------

std::string Inputs::get_potential_model()
{
  return mklower(check_settings_str("Electrodynamics", "Potential"));
}

// -----------------------------------------------------------------------
// Return the Electrodynamics Diffuse Auroral Model
// -----------------------------------------------------------------------

std::string Inputs::get_diffuse_auroral_model()
{
  return mklower(check_settings_str("Electrodynamics", "DiffuseAurora"));
}

// -----------------------------------------------------------------------
// Return the EUV model used (EUVAC only option now)
// -----------------------------------------------------------------------

std::string Inputs::get_euv_model()
{
  return mklower(check_settings_str("Euv", "Model"));
}

// -----------------------------------------------------------------------
// Return the heating efficiency of the neutrals for EUV
// -----------------------------------------------------------------------

precision_t Inputs::get_euv_heating_eff_neutrals()
{
  return check_settings_pt("Euv", "HeatingEfficiency");
}

// -----------------------------------------------------------------------
// Return whether to include the photoelectron ionization
// -----------------------------------------------------------------------

bool Inputs::get_include_photoelectrons()
{
  return get_setting_bool("Euv", "IncludePhotoElectrons");
}

// -----------------------------------------------------------------------
// Return how often to calculate EUV energy deposition
// -----------------------------------------------------------------------

precision_t Inputs::get_dt_euv()
{
  return check_settings_pt("Euv", "dt");
}

// -----------------------------------------------------------------------
// Return how often to report progress of simulation
// -----------------------------------------------------------------------

precision_t Inputs::get_dt_report()
{
  return check_settings_pt("Debug", "dt");
}

// -----------------------------------------------------------------------
// Return number of output types
// -----------------------------------------------------------------------

precision_t Inputs::get_n_outputs()
{
  return settings.at("Outputs").at("type").size();
}

// -----------------------------------------------------------------------
// Return original random number seed
// -----------------------------------------------------------------------

int Inputs::get_original_seed()
{
  return get_setting_int("Seed");
}

// -----------------------------------------------------------------------
// Return number of longitudes, latitudes, and altitudes in grid
// -----------------------------------------------------------------------

int Inputs::get_nLons(std::string gridtype)
{
  return get_setting_int(gridtype, "nLonsPerBlock");
}

int Inputs::get_nLats(std::string gridtype)
{
  return get_setting_int(gridtype, "nLatsPerBlock");
}

int Inputs::get_nAlts(std::string gridtype)
{
  return get_setting_int(gridtype, "nAlts");
}

std::string Inputs::get_grid_shape(std::string gridtype)
{
  return mklower(get_setting_str(gridtype, "Shape"));
}

// -----------------------------------------------------------------------
// Return number of ensemble members
// -----------------------------------------------------------------------

int Inputs::get_nMembers()
{
  return check_settings_pt("Ensembles", "nMembers");
}

// -----------------------------------------------------------------------
// Return verbose variables
// -----------------------------------------------------------------------

int Inputs::get_verbose()
{
  return check_settings_pt("Debug", "iVerbose");
}

int Inputs::get_verbose_proc()
{
  return get_setting_int("Debug", "iProc");
}

// -----------------------------------------------------------------------
// Return EUV file name
// -----------------------------------------------------------------------

std::string Inputs::get_euv_file()
{
  return check_settings_str("Euv", "File");
}

// -----------------------------------------------------------------------
// Return aurora file name
// -----------------------------------------------------------------------

std::string Inputs::get_aurora_file()
{
  return check_settings_str("AuroraFile");
}

// -----------------------------------------------------------------------
// Return Chemistry file name
// -----------------------------------------------------------------------

std::string Inputs::get_chemistry_file()
{
  return settings.at("ChemistryFile");
}

// -----------------------------------------------------------------------
// Return Collision file name
// -----------------------------------------------------------------------

std::string Inputs::get_collision_file()
{
  return check_settings_str("CollisionsFile");
}

// -----------------------------------------------------------------------
// Return Indices Lookup Filename
// -----------------------------------------------------------------------

std::string Inputs::get_indices_lookup_file()
{
  return check_settings_str("IndicesLookupFile");
}

// -----------------------------------------------------------------------
// Return F107 file to read
// -----------------------------------------------------------------------

std::string Inputs::get_f107_file()
{
  return check_settings_str("F107File");
}

// -----------------------------------------------------------------------
// Return planet name
// -----------------------------------------------------------------------

std::string Inputs::get_planet()
{
  return get_setting_str("Planet", "name");
}

// -----------------------------------------------------------------------
// Return file that contains (all) planetary characteristics
// -----------------------------------------------------------------------

std::string Inputs::get_planetary_file()
{
  return check_settings_str("PlanetCharacteristicsFile");
}

// -----------------------------------------------------------------------
// Return planetary file name that describes the species and such for
// a given planet
// -----------------------------------------------------------------------

std::string Inputs::get_planet_species_file()
{
  return check_settings_str("PlanetSpeciesFile");
}

// -----------------------------------------------------------------------
// Flag to do the bulk ion temperature calculation instead
// of individual ion specie temperature calculations
// -----------------------------------------------------------------------

bool Inputs::get_do_calc_bulk_ion_temp()
{
  return get_setting_bool("DoCalcBulkIonTemp");
}

// -----------------------------------------------------------------------
// Return Eddy Coefficient
// -----------------------------------------------------------------------

precision_t Inputs::get_eddy_coef()
{
  return get_setting_float("Eddy", "Coefficient");
}

// -----------------------------------------------------------------------
// Return pressure where Eddy Diffusion starts to drop off
// -----------------------------------------------------------------------

precision_t Inputs::get_eddy_bottom()
{
  return get_setting_float("Eddy", "BottomPressure");
}

// -----------------------------------------------------------------------
// Return pressure where Eddy Diffusion becomes zero
// -----------------------------------------------------------------------

precision_t Inputs::get_eddy_top()
{
  return get_setting_float("Eddy", "TopPressure");
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

bool Inputs::get_use_eddy_momentum()
{
  return get_setting_bool("Eddy", "UseInMomentum");
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

bool Inputs::get_use_eddy_energy()
{
  return get_setting_bool("Eddy", "UseInEnergy");
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

json Inputs::get_perturb_values()
{
  return get_setting_json("Perturb");
}

// -----------------------------------------------------------------------
// Flag to check neutral and ions for nans and infinites
// -----------------------------------------------------------------------

bool Inputs::get_check_for_nans()
{
  return get_setting_bool("Debug", "check_for_nans");
}

// -----------------------------------------------------------------------
// Checks to see if nan_test is needed
// -----------------------------------------------------------------------

bool Inputs::get_nan_test()
{
  return get_setting_bool("Debug", "nan_test", "insert");
}

// -----------------------------------------------------------------------
// Returns which variable is being tested for nans
// -----------------------------------------------------------------------

std::string Inputs::get_nan_test_variable()
{
  return get_setting_str("Debug", "nan_test", "variable");
}

// -----------------------------------------------------------------------
// Flag to have a latitude dependent radius, and by extension gravity
// -----------------------------------------------------------------------

bool Inputs::get_do_lat_dependent_radius()
{
  return get_setting_bool("Oblate", "isOblate");
}

// -----------------------------------------------------------------------
// Flag to include J2 term in the gravity calculation
// -----------------------------------------------------------------------

bool Inputs::get_do_J2()
{
  return get_setting_bool("Oblate", "isJ2");
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

json Inputs::get_initial_condition_types()
{
  return get_setting_json("InitialConditions");
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------

json Inputs::get_boundary_condition_types()
{
  return get_setting_json("BoundaryConditions");
}

std::string Inputs::get_advection_neutrals_vertical()
{
  return get_setting_str("Advection", "Neutrals", "Vertical");
}

// --------------------------------------------------------------------------
// check to see if class is ok
// --------------------------------------------------------------------------

bool Inputs::is_ok()
{
  return isOk;
}
