// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <math.h>
#include <iostream>
#include <fstream>

#include "aether.h"

// -----------------------------------------------------------------------------
//  Create a single species by filling the species structure
// -----------------------------------------------------------------------------

Neutrals::species_chars Neutrals::create_species(Grid grid) {

  species_chars tmp;

  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();

  // Constants:
  tmp.DoAdvect = 0;
  tmp.thermal_cond = 0.0;
  tmp.thermal_exp = 0.0;
  tmp.lower_bc_density = -1.0;

  tmp.density_scgc.set_size(nLons, nLats, nAlts);
  tmp.newDensity_scgc.set_size(nLons, nLats, nAlts);
  tmp.velocity_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);
  tmp.newVelocity_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);

  for (int iDir = 0; iDir < 3; iDir++) {
    tmp.velocity_vcgc[iDir].zeros();
    tmp.newVelocity_vcgc[iDir].zeros();
  }

  tmp.chapman_scgc.set_size(nLons, nLats, nAlts);
  tmp.scale_height_scgc.set_size(nLons, nLats, nAlts);
  tmp.ionization_scgc.set_size(nLons, nLats, nAlts);

  tmp.acc_neutral_friction = make_cube_vector(nLons, nLats, nAlts, 3);
  tmp.acc_ion_drag = make_cube_vector(nLons, nLats, nAlts, 3);
  tmp.acc_eddy.set_size(nLons, nLats, nAlts);
  tmp.ionization_scgc.zeros();

  tmp.concentration_scgc.set_size(nLons, nLats, nAlts);

  tmp.density_scgc.ones();
  tmp.chapman_scgc.ones();
  tmp.scale_height_scgc.ones();
  tmp.rho_alt_int_scgc.zeros();

  tmp.sources_scgc.set_size(nLons, nLats, nAlts);
  tmp.sources_scgc.zeros();
  tmp.losses_scgc.set_size(nLons, nLats, nAlts);
  tmp.losses_scgc.zeros();

  tmp.nAuroraIonSpecies = 0;
  tmp.Aurora_Coef = -1.0;

  return tmp;
}

// -----------------------------------------------------------------------------
//  Initialize neutrals
// -----------------------------------------------------------------------------

Neutrals::Neutrals(Grid grid,
                   Planets planet,
                   Times time,
                   Indices indices) {

  int iErr;
  bool didWork = true;
  species_chars tmp;

  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();

  report.print(2, "Initializing Neutrals");

  json planet_neutrals = planet.get_neutrals();

  for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    tmp = create_species(grid);
    species.push_back(tmp);
  }

  velocity_name.push_back("velocity_east");
  velocity_name.push_back("velocity_north");
  velocity_name.push_back("velocity_up");

  // State variables:

  density_scgc.set_size(nLons, nLats, nAlts);
  density_scgc.ones();
  temperature_scgc.set_size(nLons, nLats, nAlts);
  temperature_scgc.ones();
  newTemperature_scgc.set_size(nLons, nLats, nAlts);
  newTemperature_scgc.ones();
  O_cool_scgc.set_size(nLons, nLats, nAlts);
  O_cool_scgc.zeros();
  NO_cool_scgc.set_size(nLons, nLats, nAlts);
  NO_cool_scgc.zeros();

  // Derived quantities:

  rho_scgc.set_size(nLons, nLats, nAlts);
  rho_scgc.ones();
  velocity_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);

  for (int iDir = 0; iDir < 3; iDir++)
    velocity_vcgc[iDir].zeros();

  cMax_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);
  mean_major_mass_scgc.set_size(nLons, nLats, nAlts);
  mean_major_mass_scgc.ones();
  pressure_scgc.set_size(nLons, nLats, nAlts);
  pressure_scgc.ones();
  sound_scgc.set_size(nLons, nLats, nAlts);
  sound_scgc.ones();

  // Heating and cooling parameters:
  Cv_scgc.set_size(nLons, nLats, nAlts);
  Cv_scgc.ones();
  gamma_scgc.set_size(nLons, nLats, nAlts);
  gamma_scgc.zeros();
  kappa_scgc.set_size(nLons, nLats, nAlts);
  kappa_scgc.zeros();
  kappa_eddy_scgc.set_size(nLons, nLats, nAlts);
  kappa_eddy_scgc.zeros();

  conduction_scgc.set_size(nLons, nLats, nAlts);
  heating_euv_scgc.set_size(nLons, nLats, nAlts);
  heating_chemical_scgc.set_size(nLons, nLats, nAlts);

  heating_efficiency = input.get_euv_heating_eff_neutrals();

  // This gets a bunch of the species-dependent characteristics:
  iErr = read_planet_file(planet);

  if (iErr > 0)
    report.error("Error reading planet file!");

  // This specifies the initial conditions for the neutrals:
  didWork = initial_conditions(grid, time, indices);

  if (!didWork)
    report.error("Error in setting neutral initial conditions!");

  return;
}

// -----------------------------------------------------------------------------
// Read in the planet file that describes the species - only neutrals
// -----------------------------------------------------------------------------

int Neutrals::read_planet_file(Planets planet) {

  int iErr = 0;
  std::string hash;
  std::ifstream infile_ptr;

  report.print(3, "In read_planet_file for Neutrals");

  json neutrals = planet.get_neutrals();

  nSpecies = neutrals["name"].size();
  nSpeciesAdvect = 0;

  for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    species[iSpecies].cName = neutrals["name"][iSpecies];
    double mass = neutrals["mass"][iSpecies];
    species[iSpecies].mass = mass * cAMU;
    species[iSpecies].vibe = neutrals["vibration"][iSpecies];
    species[iSpecies].thermal_cond = neutrals["thermal_cond"][iSpecies];
    species[iSpecies].thermal_exp = neutrals["thermal_exp"][iSpecies];
    species[iSpecies].DoAdvect = neutrals["advect"][iSpecies];
    species[iSpecies].lower_bc_density = neutrals["BC"][iSpecies];
  }

  // account for advected neutrals:
  for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    if (species[iSpecies].DoAdvect == 1) {
      nSpeciesAdvect++;
      species_to_advect.push_back(iSpecies);
    }
  }

  json temperatures = planet.get_temperatures();
  nInitial_temps = temperatures["alt"].size();

  for (int i = 0; i < nInitial_temps; i++) {
    initial_altitudes.push_back(double(temperatures["alt"][i]) * 1000.0);
    initial_temperatures.push_back(temperatures["temp"][i]);
  }

  return iErr;
}

//----------------------------------------------------------------------
// Fill With Hydrostatic Solution (all species)
//   - iEnd is NOT included (python style)!
//----------------------------------------------------------------------

void Neutrals::fill_with_hydrostatic(int64_t iStart,
                                     int64_t iEnd,
                                     Grid grid) {

  for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    // Integrate with hydrostatic equilibrium up:
    for (int iAlt = iStart; iAlt < iEnd; iAlt++) {
      species[iSpecies].density_scgc.slice(iAlt) =
        temperature_scgc.slice(iAlt - 1) /
        temperature_scgc.slice(iAlt) %
        species[iSpecies].density_scgc.slice(iAlt - 1) %
        exp(-grid.dalt_lower_scgc.slice(iAlt) /
            species[iSpecies].scale_height_scgc.slice(iAlt));
    }
  }

  calc_mass_density();
  return;
}

//----------------------------------------------------------------------
// Fill With Hydrostatic Solution (only one constituent)
//   - iEnd is NOT included (python style)!
//----------------------------------------------------------------------

void Neutrals::fill_with_hydrostatic(int64_t iSpecies,
                                     int64_t iStart,
                                     int64_t iEnd,
                                     Grid grid) {

  // Integrate with hydrostatic equilibrium up:
  for (int iAlt = iStart; iAlt < iEnd; iAlt++) {
    species[iSpecies].density_scgc.slice(iAlt) =
      temperature_scgc.slice(iAlt - 1) /
      temperature_scgc.slice(iAlt) %
      species[iSpecies].density_scgc.slice(iAlt - 1) %
      exp(-grid.dalt_lower_scgc.slice(iAlt) /
          species[iSpecies].scale_height_scgc.slice(iAlt));
  }

  calc_mass_density();
  return;
}

//----------------------------------------------------------------------
// Reports location of nans inserted into specified variable
//----------------------------------------------------------------------
void Neutrals::nan_test(std::string variable) {
  std::vector<int> locations;
  std::string message = ("For Neutrals " + variable + " ");

  if (variable == "temperature_scgc") {
    locations = insert_indefinites(temperature_scgc);
    message += print_nan_vector(locations, temperature_scgc);
  }

  if (variable == "density_scgc") {
    locations = insert_indefinites(density_scgc);
    message += print_nan_vector(locations, density_scgc);
  }

  if (variable == "velocity_vcgc") {
    locations = insert_indefinites(velocity_vcgc[0]);
    message +=
      "at the x loc " + print_nan_vector(locations, velocity_vcgc[0]);
    locations = insert_indefinites(velocity_vcgc[1]);
    message +=
      "at the y loc " + print_nan_vector(locations, velocity_vcgc[1]);
    locations = insert_indefinites(velocity_vcgc[2]);
    message +=
      "at the z loc " + print_nan_vector(locations, velocity_vcgc[2]);
  }

  std::cout << message;
}

//----------------------------------------------------------------------
// Checks for nans and +/- infinities in density, temp, and velocity
//----------------------------------------------------------------------

bool Neutrals::check_for_nonfinites() {
  bool isBad = false;
  bool didWork = true;

  isBad = !all_finite(density_scgc, "density_scgc");

  if (isBad) {
    report.error("non-finite found in neutral density!");
    didWork = false;
  }

  isBad = !all_finite(temperature_scgc, "temperature_scgc");

  if (isBad) {
    report.error("non-finite found in neutral temperature!");
    didWork = false;
  }

  isBad = !all_finite(velocity_vcgc, "velocity_vcgc");

  if (isBad) {
    report.error("non-finite found in neutral velocity!");
    didWork = false;
  }

  return didWork;
}

//----------------------------------------------------------------------
// return the index of the requested species
// This will return -1 if the species is not found or name is empty
//----------------------------------------------------------------------

int Neutrals::get_species_id(std::string name) {

  std::string function = "Neutrals::get_species_id";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int iSpecies;
  int id_ = -1;

  if (name.length() > 0) {
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      if (name == species[iSpecies].cName) {
        id_ = iSpecies;
        break;
      }
  }

  report.exit(function);
  return id_;
}

//----------------------------------------------------------------------
// Read/Write restart files for the neutrals
//----------------------------------------------------------------------

bool Neutrals::restart_file(std::string dir, bool DoRead) {

  std::string filename;
  bool DidWork = true;
  int64_t iVar;
  std::string cName;

  OutputContainer RestartContainer;
  RestartContainer.set_directory(dir);
  RestartContainer.set_filename("neutrals_" + cMember + "_" + cGrid);

  try {
    if (DoRead)
      RestartContainer.read();
    else {
      RestartContainer.set_version(0.1);
      RestartContainer.set_time(0.0);
    }

    for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      cName = species[iSpecies].cName;

      if (DoRead)
        species[iSpecies].density_scgc =
          RestartContainer.get_element_value(cName);
      else
        RestartContainer.store_variable(cName,
                                        density_unit,
                                        species[iSpecies].density_scgc);
    }

    cName = temperature_name;

    if (DoRead)
      temperature_scgc = RestartContainer.get_element_value(cName);
    else
      RestartContainer.store_variable(cName,
                                      temperature_unit,
                                      temperature_scgc);

    for (int iDir = 0; iDir < 3; iDir++) {
      cName = velocity_name[iDir];

      if (DoRead)
        velocity_vcgc[iDir] = RestartContainer.get_element_value(cName);
      else
        RestartContainer.store_variable(cName,
                                        velocity_unit,
                                        velocity_vcgc[iDir]);
    }

    if (!DoRead) {
      RestartContainer.write();
      RestartContainer.clear_variables();
    }
  } catch (...) {
    std::cout << "Error reading in neutral restart file!\n";
    DidWork = false;
  }

  return DidWork;
}

//----------------------------------------------------------------------
// Calculate value of NO Cooling
//----------------------------------------------------------------------

void Neutrals::calc_NO_cool() {
  // finds O & NO species
  int iO = get_species_id("O");
  int iNO = get_species_id("NO");

  if (iNO != -1) {
    // omega value using O density
    arma_cube omega = 3.6e-17 * species[iO].density_scgc / (3.6e-17 *
                                                            species[iO].density_scgc + 13.3);

    arma_cube v = -cH * cC / (5.3e-6 * cKB * temperature_scgc);

    // calculation for NO_cool_scgc
    arma_cube NO_cool_scgc_calc = cH * cC /
                                  5.3e-6 * omega * 13.3 % exp(v) % species[iNO].density_scgc;

    NO_cool_scgc = NO_cool_scgc_calc / (rho_scgc % Cv_scgc);
  }

  return;
}

//----------------------------------------------------------------------
// Calculate value of O Cooling
//----------------------------------------------------------------------

void Neutrals::calc_O_cool() {
  // find O species
  int iO = get_species_id("O");

  if (iO != -1) {
    arma_cube tmp2 = exp(-228 / temperature_scgc);
    arma_cube tmp3 = exp(-326 / temperature_scgc);

    // calculation for O_cool_scgc
    arma_cube O_cool_scgc_calc = (1.69e-18 * tmp2 + 4.59e-20 * tmp3) %
                                 (species[iO].density_scgc / 1.0e6) / (1.0 + 0.6 * tmp2 + 0.2 * tmp3);

    O_cool_scgc = O_cool_scgc_calc / (rho_scgc % Cv_scgc);
  }

  return;
}
