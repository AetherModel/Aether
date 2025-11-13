// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>
#include <fstream>

#include "aether.h"

// -----------------------------------------------------------------------------
// Initialize a single species for the ions
// -----------------------------------------------------------------------------

Ions::species_chars Ions::create_species(Grid &grid) {

  species_chars tmp;

  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();

  // Constants:
  tmp.DoAdvect = 0;
  tmp.vibe = 5;

  tmp.density_scgc.set_size(nLons, nLats, nAlts);
  tmp.density_scgc.fill(1e10);
  tmp.temperature_scgc.set_size(nLons, nLats, nAlts);
  tmp.temperature_scgc.fill(200.0);
  tmp.ionization_scgc.set_size(nLons, nLats, nAlts);
  tmp.ionization_scgc.zeros();

  tmp.sources_scgc.set_size(nLons, nLats, nAlts);
  tmp.sources_scgc.zeros();
  tmp.losses_scgc.set_size(nLons, nLats, nAlts);
  tmp.losses_scgc.zeros();

  tmp.heating_neutral_friction_scgc.set_size(nLons, nLats, nAlts);
  tmp.heating_neutral_friction_scgc.zeros();
  tmp.heating_electron_friction_scgc.set_size(nLons, nLats, nAlts);
  tmp.heating_electron_friction_scgc.zeros();
  tmp.heating_neutral_heat_transfer_scgc.set_size(nLons, nLats, nAlts);
  tmp.heating_neutral_heat_transfer_scgc.zeros();
  tmp.heating_electron_heat_transfer_scgc.set_size(nLons, nLats, nAlts);
  tmp.heating_electron_heat_transfer_scgc.zeros();
  tmp.heating_sources_total.set_size(nLons, nLats, nAlts);
  tmp.heating_sources_total.zeros();
  tmp.Cv_scgc.set_size(nLons, nLats, nAlts);
  tmp.Cv_scgc.zeros();
  tmp.lambda.set_size(nLons, nLats, nAlts);
  tmp.lambda.zeros();

  tmp.par_velocity_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);
  tmp.perp_velocity_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);
  tmp.velocity_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);

  for (int iDir = 0; iDir < 3; iDir++) {
    tmp.par_velocity_vcgc[iDir].zeros();
    tmp.perp_velocity_vcgc[iDir].zeros();
    tmp.velocity_vcgc[iDir].zeros();
  }

  // The collision frequencies need the neutrals, so those are
  // initialized in init_ion_temperature.

  return tmp;
}

// -----------------------------------------------------------------------------
//  Initialize Ions class
// -----------------------------------------------------------------------------

Ions::Ions(Grid &grid, Planets planet) {

  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();

  species_chars tmp;

  report.print(2, "Initializing Ions");

  for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    tmp = create_species(grid);
    species.push_back(tmp);
  }

  velocity_name.push_back("velocity_east");
  velocity_name.push_back("velocity_north");
  velocity_name.push_back("velocity_up");

  par_velocity_name.push_back("velocity_parallel_east");
  par_velocity_name.push_back("velocity_parallel_north");
  par_velocity_name.push_back("velocity_parallel_up");

  perp_velocity_name.push_back("velocity_perp_east");
  perp_velocity_name.push_back("velocity_perp_north");
  perp_velocity_name.push_back("velocity_perp_up");

  // Create one extra species for electrons
  tmp = create_species(grid);
  species.push_back(tmp);

  // State variables:

  density_scgc.set_size(nLons, nLats, nAlts);
  density_scgc.ones();
  velocity_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);
  cMax_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);

  // Some variables that will be used in calc_ion_v:
  gravity_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);
  wind_acc = make_cube_vector(nLons, nLats, nAlts, 3);
  total_acc = make_cube_vector(nLons, nLats, nAlts, 3);
  efield_acc = make_cube_vector(nLons, nLats, nAlts, 3);
  a_par = make_cube_vector(nLons, nLats, nAlts, 3);
  a_perp = make_cube_vector(nLons, nLats, nAlts, 3);
  a_x_b = make_cube_vector(nLons, nLats, nAlts, 3);
  grad_Pi_plus_Pe = make_cube_vector(nLons, nLats, nAlts, 3);
  rho.set_size(nLons, nLats, nAlts);
  nuin.set_size(nLons, nLats, nAlts);
  nuin_sum.set_size(nLons, nLats, nAlts);
  Nie.set_size(nLons, nLats, nAlts);
  sum_rho.set_size(nLons, nLats, nAlts);
  top.set_size(nLons, nLats, nAlts);
  bottom.set_size(nLons, nLats, nAlts);



  Cv_scgc.set_size(nLons, nLats, nAlts);
  Cv_scgc.zeros();
  lambda.set_size(nLons, nLats, nAlts);
  lambda.zeros();
  gamma_scgc.set_size(nLons, nLats, nAlts);
  gamma_scgc.ones();
  sound_scgc.set_size(nLons, nLats, nAlts);
  sound_scgc.ones();

  heating_neutral_friction_scgc.set_size(nLons, nLats, nAlts);
  heating_neutral_friction_scgc.zeros();
  heating_electron_friction_scgc.set_size(nLons, nLats, nAlts);
  heating_electron_friction_scgc.zeros();
  heating_neutral_heat_transfer_scgc.set_size(nLons, nLats, nAlts);
  heating_neutral_heat_transfer_scgc.zeros();
  heating_electron_heat_transfer_scgc.set_size(nLons, nLats, nAlts);
  heating_electron_heat_transfer_scgc.zeros();
  heating_sources_total.set_size(nLons, nLats, nAlts);
  heating_sources_total.zeros();

  for (int iDir = 0; iDir < 3; iDir++) {
    velocity_vcgc[iDir].zeros();
    cMax_vcgc[iDir].ones();
  }

  temperature_scgc.set_size(nLons, nLats, nAlts);
  temperature_scgc.fill(200.0);
  electron_temperature_scgc.set_size(nLons, nLats, nAlts);
  electron_temperature_scgc.fill(200);

  rho_scgc.set_size(nLons, nLats, nAlts);
  mean_major_mass_scgc.set_size(nLons, nLats, nAlts);
  mean_major_mass_scgc.ones();

  tmp.sources_scgc.set_size(nLons, nLats, nAlts);
  tmp.sources_scgc.zeros();
  tmp.losses_scgc.set_size(nLons, nLats, nAlts);
  tmp.losses_scgc.zeros();

  potential_scgc.set_size(nLons, nLats, nAlts);
  potential_scgc.zeros();

  conduction_scgc.set_size(nLons, nLats, nAlts);

  eflux.set_size(nLons, nLats);
  eflux.zeros();
  avee.set_size(nLons, nLats);
  avee.zeros();

  // Make vectors:
  efield_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);
  exb_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);

  // This gets a bunch of the species-dependent characteristics:
  int iErr = read_planet_file(planet);

  if (input.get_do_restart()) {
    report.print(1, "Restarting! Reading ion files!");
    bool DidWork = restart_file(input.get_restartin_dir(), grid.get_gridtype(),
                                DoRead);

    if (!DidWork)
      std::cout << "Reading Restart for Ions Failed!!!\n";
  }

  if (iErr > 0)
    std::cout << "Error in reading planet file!" << '\n';
}

// -----------------------------------------------------------------------------
// Read in the planet file that describes the species - only ions
// -----------------------------------------------------------------------------

int Ions::read_planet_file(Planets planet) {

  int iErr = 0;
  std::string hash;
  std::ifstream infile_ptr;
  int doAdvect;

  report.print(3, "In read_planet_file for Ions");

  json ions = planet.get_ions();

  nSpecies = ions["name"].size();

  for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    species[iSpecies].cName = ions["name"][iSpecies];
    double mass = ions["mass"][iSpecies];
    species[iSpecies].mass = mass * cAMU;
    species[iSpecies].charge = ions["charge"][iSpecies];
    doAdvect = ions["advect"][iSpecies];

    if (doAdvect == 0)
      species[iSpecies].DoAdvect = false;
    else
      species[iSpecies].DoAdvect = true;

    species[iSpecies].vibe = ions["vibration"][iSpecies];
  }

  // account for advected ions:
  for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    if (species[iSpecies].DoAdvect == 1) {
      nSpeciesAdvect++;
      species_to_advect.push_back(iSpecies);
    }
  }

  species[nSpecies].cName = "e-";
  species[nSpecies].mass = cME;
  species[nSpecies].charge = -1;
  species[nSpecies].DoAdvect = 0;

  return iErr;
}

//----------------------------------------------------------------------
// Reports location of nans inserted into specified variable
//----------------------------------------------------------------------
void Ions::nan_test(std::string variable) {
  std::vector<int> locations;
  std::string message = ("For Ions " + variable + " ");

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

bool Ions::check_for_nonfinites(std::string location) {
  bool didWork = true;

  if (!all_finite(density_scgc, "density_scgc") ||
      !all_finite(temperature_scgc, "temperature_scgc") ||
      !all_finite(velocity_vcgc, "velocity_vcgc"))
    didWork = false;

  if (!didWork)
    report.error("ions are nan from location : " + location);

  return didWork;
}

// -----------------------------------------------------------------------------
// Set a floor for ion densities
// -----------------------------------------------------------------------------

void Ions::set_floor() {

  int iSpecies;

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    species[iSpecies].density_scgc.clamp(1.0, 1e15);

  return;
}

// ----------------------------------------------------------------------
// Calculate a bunch of derived products:
//   - Gamma
//   - Speed of sound
// ----------------------------------------------------------------------

void Ions::calc_sound_speed() {

  int64_t iSpecies;

  std::string function = "Ion::calc_sound_speed";
  static int iFunction = -1;
  report.enter(function, iFunction);

  gamma_scgc.zeros();

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    gamma_scgc = gamma_scgc +
                 species[iSpecies].density_scgc / (species[iSpecies].vibe - 2);
  }

  gamma_scgc = gamma_scgc * 2.0 / density_scgc + 1.0;

  sound_scgc = sqrt(cKB *
                    gamma_scgc %
                    temperature_scgc /
                    mean_major_mass_scgc);

  if (report.test_verbose(2)) {
    std::cout << "max sound speed : " << sound_scgc.max() << "\n";
    std::cout << "max gamma : " << gamma_scgc.max() << "\n";

    if (!all_finite(sound_scgc, "sound speed")) {
      std::cout << "sound speed has nans!\n";
      report.report_errors();
    }
  }

  report.exit(function);
  return;
}

// ----------------------------------------------------------------------
// Calculate cMax, which is the sound speed + velocity in each
// direction
// ----------------------------------------------------------------------

void Ions::calc_cMax() {

  std::string function = "Ions::calc_cMax";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int iDir;

  // just take the bulk value for now:

  if (report.test_verbose(3)) {
    std::cout << "max sound speed : " << sound_scgc.max() << "\n";

    for (iDir = 0; iDir < 3; iDir++) {
      arma_cube tmp = abs(velocity_vcgc[iDir]);
      std::cout << "min velocity : " << tmp.min() << "\n";
      std::cout << "max velocity : " << tmp.max() << "\n";
    }
  }

  for (iDir = 0; iDir < 3; iDir++)
    cMax_vcgc[iDir] = sound_scgc + abs(velocity_vcgc[iDir]);

  report.exit(function);
  return;
}

// ----------------------------------------------------------------------
// Calculate thermal conduction (lambda)
// ----------------------------------------------------------------------

void Ions::calc_lambda() {

  int64_t iIon, jIon;

  std::string function = "Ions::calc_specific_heat";
  static int iFunction = -1;
  report.enter(function, iFunction);

  lambda.zeros();
  precision_t Mi, Mj;
  arma_cube density_ratio, ratios;

  // This sets the size:
  ratios = density_scgc;

  for (iIon = 0; iIon < nSpecies; iIon++) {

    Mi = species[iIon].mass / cAMU;
    ratios.zeros();

    for (jIon = 0; jIon < nSpecies; jIon++) {
      if (jIon != iIon) {
        Mj = species[jIon].mass / cAMU;
        density_ratio = species[jIon].density_scgc /
                        species[iIon].density_scgc;
        density_ratio.clamp(0.001, 1000.0);
        ratios = ratios + density_ratio *
                 (species[jIon].charge * species[jIon].charge /
                  species[iIon].charge / species[iIon].charge) *
                 sqrt(Mj / (Mi + Mj)) *
                 (3 * Mi * Mi + 1.6 * Mi * Mj + 1.3 * Mj * Mj) /
                 ((Mi + Mj) * (Mi + Mj));
      }

      species[iIon].lambda =
        3.1e6 / sqrt(Mi) / pow(species[iIon].charge, 4) *
        pow(species[iIon].temperature_scgc, 2.5) % (1 + 1.75 * ratios) * cE;
    }

    lambda = lambda + species[iIon].lambda % species[iIon].density_scgc;
  }

  lambda = lambda / density_scgc;

  //lambda1d = 25.0 * cKB * pow(temp1d, 2.5) * (cKB / species[iIon].mass)
  //           / species[iIon].nu_ion_ion[iIon] / 8.0;


  report.exit(function);
  return;
}


// ----------------------------------------------------------------------
// Calculate Specific Heat at Constant Volume (Cv)
// ----------------------------------------------------------------------

void Ions::calc_specific_heat() {

  int64_t iSpecies;

  std::string function = "Ions::calc_specific_heat";
  static int iFunction = -1;
  report.enter(function, iFunction);

  Cv_scgc.zeros();
  gamma_scgc.zeros();

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    // individual Cv for individual temperatures:
    species[iSpecies].Cv_scgc =
      (species[iSpecies].vibe - 2) * cKB / species[iSpecies].mass / 2;
    // Bulk Cv for the bulk temperature:
    Cv_scgc = Cv_scgc +
              (species[iSpecies].vibe - 2) *
              species[iSpecies].density_scgc *
              cKB / species[iSpecies].mass;
    gamma_scgc = gamma_scgc +
                 species[iSpecies].density_scgc / (species[iSpecies].vibe - 2);
  }

  // Bulk Cv and gamma are the density-weighted Cv and gamma
  Cv_scgc = Cv_scgc / (2 * density_scgc);
  gamma_scgc = gamma_scgc * 2.0 / density_scgc + 1.0;

  report.exit(function);
  return;
}


// -----------------------------------------------------------------------------
// Calculate the electron density from the sum of all ion species
// -----------------------------------------------------------------------------

void Ions::fill_electrons() {

  int iSpecies;

  std::string function = "Ions::fill_electrons";
  static int iFunction = -1;
  report.enter(function, iFunction);

  species[nSpecies].density_scgc.zeros();
  rho_scgc.zeros();

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    // Electron density is the last species:
    species[nSpecies].density_scgc =
      species[nSpecies].density_scgc + species[iSpecies].density_scgc;
    // While we are at it, calculate the mass density too:
    rho_scgc = rho_scgc +
               species[iSpecies].mass * species[iSpecies].density_scgc;
  }

  // The electron density is also stored as the bulk density:
  density_scgc = species[nSpecies].density_scgc;
  // We can now calculate the mean major mass for ions too:
  mean_major_mass_scgc = rho_scgc / density_scgc;

  report.exit(function);
  return;
}

//----------------------------------------------------------------------
// return the index of the requested species
// This will return -1 if the species is not found or name is empty
// Will return nSpecies for electrons
//----------------------------------------------------------------------

int Ions::get_species_id(const std::string &name) const {

  std::string function = "Ions::get_species_id";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int iSpecies;
  int id_ = -1;

  if (name.length() > 0) {
    for (iSpecies = 0; iSpecies <= nSpecies; iSpecies++)
      if (name == species[iSpecies].cName) {
        id_ = iSpecies;
        break;
      }
  }

  report.exit(function);
  return id_;
}

//----------------------------------------------------------------------
// Read/Write restart files for the ions
//----------------------------------------------------------------------

bool Ions::restart_file(std::string dir, std::string cGridtype, bool DoRead) {

  std::string filename;
  bool DidWork = true;
  int64_t iVar;

  OutputContainer RestartContainer;
  RestartContainer.set_directory(dir);
  RestartContainer.set_filename("ions_" + cMember + "_" + cGrid + "_" +
                                cGridtype);

  try {
    if (DoRead)
      RestartContainer.read();
    else {
      RestartContainer.set_version(aether_version);
      RestartContainer.set_time(0.0);
    }

    std::string cName;

    for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      // ----------------------------
      // Density and Temperature (per ion)
      // ----------------------------
      cName = species[iSpecies].cName;

      if (DoRead)
        species[iSpecies].density_scgc =
          RestartContainer.get_element_value(cName);
      else
        RestartContainer.store_variable(cName,
                                        density_unit,
                                        species[iSpecies].density_scgc);

      cName = temperature_name + " (" + species[iSpecies].cName + ")";

      if (DoRead)
        species[iSpecies].temperature_scgc =
          RestartContainer.get_element_value(cName);
      else
        RestartContainer.store_variable(cName,
                                        temperature_unit,
                                        species[iSpecies].temperature_scgc);

      // ----------------------------
      // Parallel Velocity (per ion)
      // ----------------------------
      for (int iDir = 0; iDir < 3; iDir++) {
        cName = par_velocity_name[iDir] + " (" + species[iSpecies].cName + ")";

        if (DoRead)
          species[iSpecies].par_velocity_vcgc[iDir] =
            RestartContainer.get_element_value(cName);
        else
          RestartContainer.store_variable(cName,
                                          velocity_unit,
                                          species[iSpecies].
                                          par_velocity_vcgc[iDir]);
      }

      // ----------------------------
      // Perpendicular Velocity (per ion)
      // ----------------------------
      for (int iDir = 0; iDir < 3; iDir++) {
        cName = perp_velocity_name[iDir] + " (" + species[iSpecies].cName + ")";

        if (DoRead)
          species[iSpecies].perp_velocity_vcgc[iDir] =
            RestartContainer.get_element_value(cName);
        else
          RestartContainer.store_variable(cName,
                                          velocity_unit,
                                          species[iSpecies].
                                          perp_velocity_vcgc[iDir]);
      }
    }

    // ----------------------------
    // Bulk ion and electron temperatures
    // ----------------------------
    cName = temperature_name + " (bulk ion)";

    if (DoRead)
      temperature_scgc = RestartContainer.get_element_value(cName);
    else
      RestartContainer.store_variable(cName,
                                      temperature_unit,
                                      temperature_scgc);

    cName = temperature_name + " (electron)";

    if (DoRead)
      electron_temperature_scgc = RestartContainer.get_element_value(cName);
    else
      RestartContainer.store_variable(temperature_name + " (electron)",
                                      temperature_unit,
                                      electron_temperature_scgc);

    if (!DoRead)
      RestartContainer.write();

    RestartContainer.clear_variables();
  } catch (...) {
    std::cout << "Error reading in ion restart file!\n";
    DidWork = false;
  }

  return DidWork;
}
