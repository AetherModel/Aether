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

Neutrals::Neutrals(Grid grid, Planets planet, Inputs input, Report report) {

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
  iErr = initial_conditions(grid, input, report);

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

// -----------------------------------------------------------------------------
//  This is a place holder for creating initial conditions
// -----------------------------------------------------------------------------

int Neutrals::initial_conditions(Grid grid, Inputs input, Report report) {

  int iErr = 0;
  int64_t iLon, iLat, iAlt, iA;
  precision_t alt, r;

  report.print(3, "Creating Neutrals initial_condition");

  if (input.get_do_restart()) {
    report.print(1, "Restarting! Reading neutral files!");
    bool DidWork = restart_file(input.get_restartin_dir(), DoRead);

    if (!DidWork)
      std::cout << "Reading Restart for Neutrals Failed!!!\n";
  } else {

    // ---------------------------------------------------------------------
    // This section assumes we want a hydrostatic solution given the
    // temperature profile in the planet.in file.
    // ---------------------------------------------------------------------

    int64_t nLons = grid.get_nLons();
    int64_t nLats = grid.get_nLats();
    int64_t nAlts = grid.get_nAlts();

    // Let's assume that the altitudes are not dependent on lat/lon:

    arma_vec alt1d(nAlts);
    arma_vec temp1d(nAlts);

    arma_mat H2d(nLons, nLats);

    alt1d = grid.geoAlt_scgc.tube(0, 0);

    if (nInitial_temps > 0) {
      for (iAlt = 0; iAlt < nAlts; iAlt++) {
        alt = alt1d(iAlt);

        // Find temperatures:
        if (alt <= initial_altitudes[0])
          temp1d[iAlt] = initial_temperatures[0];

        else {
          if (alt >= initial_altitudes[nInitial_temps - 1])
            temp1d[iAlt] = initial_temperatures[nInitial_temps - 1];

          else {
            // Linear interpolation!
            iA = 0;

            while (alt > initial_altitudes[iA])
              iA++;

            iA--;
            // alt will be between iA and iA+1:
            r = (alt - initial_altitudes[iA]) /
                (initial_altitudes[iA + 1] - initial_altitudes[iA]);
            temp1d[iAlt] =
              (1.0 - r) * initial_temperatures[iA] +
              (r) * initial_temperatures[iA + 1];
          }
        }
      }
    } else
      temp1d = 200.0;

    // spread the 1D temperature across the globe:
    for (iLon = 0; iLon < nLons; iLon++) {
      for (iLat = 0; iLat < nLats; iLat++)
        temperature_scgc.tube(iLon, iLat) = temp1d;
    }

    // Set the lower boundary condition:
    for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      species[iSpecies].density_scgc.slice(0).
      fill(species[iSpecies].lower_bc_density);
    }

    fill_with_hydrostatic(grid, report);
  }

  return iErr;
}

//----------------------------------------------------------------------
// Fill With Hydrostatic Solution
//----------------------------------------------------------------------

void Neutrals::fill_with_hydrostatic(Grid grid, Report report) {

  int64_t nAlts = grid.get_nAlts();

  for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

    // Integrate with hydrostatic equilibrium up:
    for (int iAlt = 1; iAlt < nAlts; iAlt++) {
      //CHANGED:
      species[iSpecies].scale_height_scgc.slice(iAlt) =
        cKB * temperature_scgc.slice(iAlt) /
        (species[iSpecies].mass * -grid.gravity_vcgc[2].slice(iAlt));
      species[iSpecies].density_scgc.slice(iAlt) =
        species[iSpecies].density_scgc.slice(iAlt - 1) %
        exp(-grid.dalt_lower_scgc.slice(iAlt) /
            species[iSpecies].scale_height_scgc.slice(iAlt));
    }
  }

  calc_mass_density(report);
}

//----------------------------------------------------------------------
// set_bcs
//----------------------------------------------------------------------

void Neutrals::set_bcs(Report &report) {

  std::string function = "Neutrals::set_bcs";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t nAlts = temperature_scgc.n_slices;

  // Set the lower boundary condition:
  for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    species[iSpecies].density_scgc.slice(0).
    fill(species[iSpecies].lower_bc_density);
  }

  temperature_scgc.slice(nAlts - 2) = temperature_scgc.slice(nAlts - 3);
  temperature_scgc.slice(nAlts - 1) = temperature_scgc.slice(nAlts - 2);

  report.exit(function);
}

//----------------------------------------------------------------------
// set_horizontal_bcs
//   iDir tells which direction to set:
//      iDir = 0 -> +x
//      iDir = 1 -> +y
//      iDir = 2 -> -x
//      iDir = 3 -> -y
//----------------------------------------------------------------------

void Neutrals::set_horizontal_bcs(int64_t iDir, Grid grid, Report &report) {

  std::string function = "Neutrals::set_horizontal_bcs";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t nX = grid.get_nX(), iX;
  int64_t nY = grid.get_nY(), iY;
  int64_t nAlts = grid.get_nAlts(), iAlt;
  int64_t nGCs = grid.get_nGCs();
  int64_t iV;

  // iDir = 0 is right BC:
  if (iDir == 0) {
    for (iX = nX - nGCs; iX < nX; iX++) {
      for (iY = 0; iY < nY; iY++) {
        // Constant Gradient for Temperature:
        temperature_scgc.tube(iX, iY) =
          2 * temperature_scgc.tube(iX - 1, iY) -
          temperature_scgc.tube(iX - 2, iY);

        // Constant Value for Velocity:
        for (iV = 0; iV < 3; iV++)
          velocity_vcgc[iV].tube(iX, iY) = velocity_vcgc[iV].tube(iX - 1, iY);

        // Constant Gradient for densities:
        for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++)
          species[iSpecies].density_scgc.tube(iX, iY) =
            2 * species[iSpecies].density_scgc.tube(iX - 1, iY) -
            species[iSpecies].density_scgc.tube(iX - 2, iY);
      }
    }
  }

  // iDir = 2 is left BC:
  if (iDir == 2) {
    for (iX = nGCs - 1; iX >= 0; iX--) {
      for (iY = 0; iY < nY; iY++) {
        // Constant Gradient for Temperature:
        temperature_scgc.tube(iX, iY) =
          2 * temperature_scgc.tube(iX + 1, iY) -
          temperature_scgc.tube(iX + 2, iY);

        // Constant Value for Velocity:
        for (iV = 0; iV < 3; iV++)
          velocity_vcgc[iV].tube(iX, iY) = velocity_vcgc[iV].tube(iX + 1, iY);

        // Constant Gradient for densities:
        for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++)
          species[iSpecies].density_scgc.tube(iX, iY) =
            2 * species[iSpecies].density_scgc.tube(iX + 1, iY) -
            species[iSpecies].density_scgc.tube(iX + 2, iY);
      }
    }
  }

  // iDir = 1 is upper BC:
  if (iDir == 1) {
    for (iX = 0; iX < nX; iX++) {
      for (iY = nX - nGCs; iY < nY; iY++) {
        // Constant Gradient for Temperature:
        temperature_scgc.tube(iX, iY) =
          2 * temperature_scgc.tube(iX, iY - 1) -
          temperature_scgc.tube(iX, iY - 2);

        // Constant Value for Velocity:
        for (iV = 0; iV < 3; iV++)
          velocity_vcgc[iV].tube(iX, iY) = velocity_vcgc[iV].tube(iX, iY - 1);

        // Constant Gradient for densities:
        for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++)
          species[iSpecies].density_scgc.tube(iX, iY) =
            2 * species[iSpecies].density_scgc.tube(iX, iY - 1) -
            species[iSpecies].density_scgc.tube(iX, iY - 2);
      }
    }
  }

  // iDir = 2 is left BC:
  if (iDir == 3) {
    for (iX = 0; iX < nX; iX++) {
      for (iY = nGCs - 1; iY >= 0; iY--) {
        // Constant Gradient for Temperature:
        temperature_scgc.tube(iX, iY) =
          2 * temperature_scgc.tube(iX, iY + 1) -
          temperature_scgc.tube(iX, iY + 2);

        // Constant Value for Velocity:
        for (iV = 0; iV < 3; iV++)
          velocity_vcgc[iV].tube(iX, iY) = velocity_vcgc[iV].tube(iX, iY + 1);

        // Constant Gradient for densities:
        for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++)
          species[iSpecies].density_scgc.tube(iX, iY) =
            2 * species[iSpecies].density_scgc.tube(iX, iY + 1) -
            species[iSpecies].density_scgc.tube(iX, iY + 2);
      }
    }
  }

  report.exit(function);
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

