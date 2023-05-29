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
  tmp.chapman_scgc.set_size(nLons, nLats, nAlts);
  tmp.scale_height_scgc.set_size(nLons, nLats, nAlts);
  tmp.ionization_scgc.set_size(nLons, nLats, nAlts);

  tmp.density_scgc.ones();
  tmp.chapman_scgc.ones();
  tmp.scale_height_scgc.ones();
  tmp.rho_alt_int_scgc.zeros();

  tmp.ionization_scgc.zeros();

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
		   Indices indices,
		   Inputs input,
		   Report report) {

  int iErr;
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

  velocity_name.push_back("Zonal Wind");
  velocity_name.push_back("Meridional Wind");
  velocity_name.push_back("Vertical Wind");

  // State variables:

  density_scgc.set_size(nLons, nLats, nAlts);
  density_scgc.ones();
  temperature_scgc.set_size(nLons, nLats, nAlts);
  temperature_scgc.ones();

  // Derived quantities:

  rho_scgc.set_size(nLons, nLats, nAlts);
  rho_scgc.ones();
  velocity_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);
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

  conduction_scgc.set_size(nLons, nLats, nAlts);
  heating_euv_scgc.set_size(nLons, nLats, nAlts);

  heating_efficiency = input.get_euv_heating_eff_neutrals();

  // This gets a bunch of the species-dependent characteristics:
  iErr = read_planet_file(planet, input, report);

  if (iErr > 0)
    std::cout << "Error reading planet file!" << '\n';

  // This specifies the initial conditions for the neutrals:
  iErr = initial_conditions(grid, time, indices, input, report);

  if (iErr > 0)
    std::cout << "Error in setting neutral initial conditions!" << '\n';
}

// -----------------------------------------------------------------------------
// Read in the planet file that describes the species - only neutrals
// -----------------------------------------------------------------------------

int Neutrals::read_planet_file(Planets planet, Inputs input, Report report) {

  int iErr = 0;
  std::string hash;
  std::ifstream infile_ptr;

  report.print(3, "In read_planet_file for Neutrals");

  json neutrals = planet.get_neutrals();

  nSpecies = neutrals["name"].size();

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
//----------------------------------------------------------------------

void Neutrals::fill_with_hydrostatic(Grid grid, Report report) {

  int64_t nAlts = grid.get_nAlts();

  for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

    // Integrate with hydrostatic equilibrium up:
    species[iSpecies].scale_height_scgc =
      cKB * temperature_scgc / (species[iSpecies].mass * grid.gravity_scgc);
    
    for (int iAlt = 1; iAlt < nAlts; iAlt++) {
      species[iSpecies].density_scgc.slice(iAlt) =
        species[iSpecies].density_scgc.slice(iAlt - 1) %
        exp(-grid.dalt_lower_scgc.slice(iAlt) /
            species[iSpecies].scale_height_scgc.slice(iAlt));
    }
  }

  calc_mass_density(report);
}

//----------------------------------------------------------------------
// Fill With Hydrostatic Solution (only one constituent)
//----------------------------------------------------------------------

void Neutrals::fill_with_hydrostatic(int64_t iSpecies,
				     Grid grid, Report report) {

  int64_t nAlts = grid.get_nAlts();

  species[iSpecies].scale_height_scgc =
    cKB * temperature_scgc / (species[iSpecies].mass * grid.gravity_scgc);
    
  // Integrate with hydrostatic equilibrium up:
  for (int iAlt = 1; iAlt < nAlts; iAlt++) {
    species[iSpecies].density_scgc.slice(iAlt) =
      species[iSpecies].density_scgc.slice(iAlt - 1) %
      exp(-grid.dalt_lower_scgc.slice(iAlt) /
	  species[iSpecies].scale_height_scgc.slice(iAlt));
  }
  calc_mass_density(report);
}

//----------------------------------------------------------------------
// return the index of the requested species
// This will return -1 if the species is not found or name is empty
//----------------------------------------------------------------------

int Neutrals::get_species_id(std::string name, Report &report) {

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
  RestartContainer.set_netcdf();
  RestartContainer.set_directory(dir);
  RestartContainer.set_filename("neutrals_" + cMember + "_" + cGrid);

  try {
    if (DoRead)
      RestartContainer.read_container_netcdf();
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

