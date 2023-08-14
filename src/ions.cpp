// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>
#include <fstream>

#include "aether.h"

// -----------------------------------------------------------------------------
// Initialize a single species for the ions
// -----------------------------------------------------------------------------

Ions::species_chars Ions::create_species(Grid grid) {

  species_chars tmp;

  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();

  // Constants:
  tmp.DoAdvect = 0;

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

  tmp.par_velocity_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);
  tmp.perp_velocity_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);

  for (int iDir = 0; iDir < 3; iDir++) {
    tmp.par_velocity_vcgc[iDir].zeros();
    tmp.perp_velocity_vcgc[iDir].zeros();
  }

  // The collision frequencies need the neutrals, so those are
  // initialized in init_ion_temperature.

  return tmp;
}

// -----------------------------------------------------------------------------
//  Initialize Ions class
// -----------------------------------------------------------------------------

Ions::Ions(Grid grid, Planets planet) {

  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();

  species_chars tmp;

  report.print(2, "Initializing Ions");

  for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    tmp = create_species(grid);
    species.push_back(tmp);
  }

  velocity_name.push_back("Ion Velocity (Zonal)");
  velocity_name.push_back("Ion Velocity (Meridional)");
  velocity_name.push_back("Ion Velocity (Vertical)");

  par_velocity_name.push_back("Parallel Ion Velocity (Zonal)");
  par_velocity_name.push_back("Parallel Ion Velocity (Meridional)");
  par_velocity_name.push_back("Parallel Ion Velocity (Vertical)");

  perp_velocity_name.push_back("Perp. Ion Velocity (Zonal)");
  perp_velocity_name.push_back("Perp. Ion Velocity (Meridional)");
  perp_velocity_name.push_back("Perp. Ion Velocity (Vertical)");

  // Create one extra species for electrons
  tmp = create_species(grid);
  species.push_back(tmp);

  // State variables:

  density_scgc.set_size(nLons, nLats, nAlts);
  density_scgc.ones();
  velocity_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);

  for (int iDir = 0; iDir < 3; iDir++)
    velocity_vcgc[iDir].zeros();

  temperature_scgc.set_size(nLons, nLats, nAlts);
  temperature_scgc.fill(200.0);
  electron_temperature_scgc.set_size(nLons, nLats, nAlts);
  electron_temperature_scgc.fill(200);

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
    bool DidWork = restart_file(input.get_restartin_dir(), DoRead);

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

  report.print(3, "In read_planet_file for Ions");

  json ions = planet.get_ions();

  nSpecies = ions["name"].size();

  for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    species[iSpecies].cName = ions["name"][iSpecies];
    double mass = ions["mass"][iSpecies];
    species[iSpecies].mass = mass * cAMU;
    species[iSpecies].charge = ions["charge"][iSpecies];
    species[iSpecies].DoAdvect = ions["advect"][iSpecies];
  }

  species[nSpecies].cName = "e-";
  species[nSpecies].mass = cME;
  species[nSpecies].charge = -1;
  species[nSpecies].DoAdvect = 0;

  return iErr;
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

// -----------------------------------------------------------------------------
// Calculate the electron density from the sum of all ion species
// -----------------------------------------------------------------------------

void Ions::fill_electrons() {

  int iSpecies;

  std::string function = "Ions::fill_electrons";
  static int iFunction = -1;
  report.enter(function, iFunction);

  species[nSpecies].density_scgc.zeros();

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    species[nSpecies].density_scgc =
      species[nSpecies].density_scgc + species[iSpecies].density_scgc;

  density_scgc = species[nSpecies].density_scgc;

  report.exit(function);
  return;
}

//----------------------------------------------------------------------
// return the index of the requested species
// This will return -1 if the species is not found or name is empty
// Will return nSpecies for electrons
//----------------------------------------------------------------------

int Ions::get_species_id(std::string name) {

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

bool Ions::restart_file(std::string dir, bool DoRead) {

  std::string filename;
  bool DidWork = true;
  int64_t iVar;

  OutputContainer RestartContainer;
  RestartContainer.set_directory(dir);
  RestartContainer.set_filename("ions_" + cMember + "_" + cGrid);

  try {
    if (DoRead)
      RestartContainer.read_container_netcdf();
    else {
      RestartContainer.set_version(0.1);
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

