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

Neutrals::Neutrals(Grid grid, Inputs input, Report report) {

  int iErr;
  species_chars tmp;

  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();

  report.print(2, "Initializing Neutrals");

  for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    tmp = create_species(grid);
    species.push_back(tmp);
  }

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

  initial_temperatures = NULL;
  initial_altitudes = NULL;

  // This gets a bunch of the species-dependent characteristics:
  iErr = read_planet_file(input, report);

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

int Neutrals::read_planet_file(Inputs input, Report report) {

  int iErr = 0;
  std::string hash;
  std::ifstream infile_ptr;

  report.print(3, "In read_planet_file for Neutrals");

  infile_ptr.open(input.get_planet_species_file());

  if (!infile_ptr.is_open()) {
    std::cout << "Could not open input file: "
              << input.get_planet_species_file() << "!!!\n";
    iErr = 1;
  } else {

    int IsDone = 0;

    while (!IsDone) {

      hash = find_next_hash(infile_ptr);

      if (report.test_verbose(4))
        std::cout << "hash : -->" << hash << "<--\n";

      if (hash == "#neutrals") {
        // Read in the characteristics as CSVs:
        report.print(4, "Found #neutrals!");

        std::vector<std::vector<std::string>> lines = read_csv(infile_ptr);

        // I should totally redo the initialization of the species,
        // since we could just do it here, but that is for the future.

        if (lines.size() - 1 != nSpecies) {
          std::cout << "number of neutrals species defined in planet.h file : "
                    << nSpecies << "\n";
          std::cout << "number of species defined in planet.in file : "
                    << lines.size() << "\n";
          std::cout << "These don't match!\n";
          iErr = 1;
        } else {
          // assume order of rows right now:
          // name, mass, vibration, thermal_cond, thermal_exp, advect, lower BC

          for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
            report.print(5, "setting neutral species " + lines[iSpecies + 1][0]);
            species[iSpecies].cName = lines[iSpecies + 1][0];
            species[iSpecies].mass = stof(lines[iSpecies + 1][1]) * cAMU;
            species[iSpecies].vibe = stof(lines[iSpecies + 1][2]);
            species[iSpecies].thermal_cond = stof(lines[iSpecies + 1][3]);
            species[iSpecies].thermal_exp = stof(lines[iSpecies + 1][4]);
            species[iSpecies].DoAdvect = stoi(lines[iSpecies + 1][5]);
            species[iSpecies].lower_bc_density = stof(lines[iSpecies + 1][6]);
          }  // iSpecies
        }  // else size
      }  // #neutrals

      if (hash == "#temperature") {

        report.print(4, "Found #temperatures!");

        std::vector<std::vector<std::string>> temps = read_csv(infile_ptr);

        int nTemps = temps.size() - 1;
        initial_temperatures =
          static_cast<float*>(malloc(nTemps * sizeof(float)));
        initial_altitudes =
          static_cast<float*>(malloc(nTemps * sizeof(float)));

        for (int iTemp = 0; iTemp < nTemps; iTemp++) {
          report.print(5, "reading initial temp alt " + temps[iTemp + 1][0]);
          // convert altitudes from km to m
          initial_altitudes[iTemp] = stof(temps[iTemp + 1][0]) * 1000;
          initial_temperatures[iTemp] = stof(temps[iTemp + 1][1]);
        }  // for iTemp

        nInitial_temps = nTemps;
      }  // #temperature

      if (infile_ptr.eof())
        IsDone = 1;
    }   // while !IsDone

    infile_ptr.close();
  }  // else (isopen)

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
      species[iSpecies].scale_height_scgc.slice(iAlt) =
        cKB * temperature_scgc.slice(iAlt) /
        (species[iSpecies].mass * grid.gravity_scgc.slice(iAlt));
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
  json description;
  
  // Output Densities
  for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    filename = dir + "/neu_s"+tostr(iSpecies,2)+"_n.bin";
    if (DidWork)
      if (DoRead)
	DidWork = species[iSpecies].density_scgc.load(filename);
      else {
	DidWork = species[iSpecies].density_scgc.save(filename);
	description["density"][species[iSpecies].cName] = filename;
      }
  }

  // Output Temperature
  filename = dir + "/neu_t.bin";
  if (DidWork)
    if (DoRead)
      DidWork = temperature_scgc.load(filename);
    else {
      DidWork = temperature_scgc.save(filename);
      description["temperature"]["bulk"] = filename;
    }

  // Output Velocity
  for (int iComp = 0; iComp < 3; iComp++) {
    filename = dir + "/neu_v"+tostr(iComp,1)+".bin";
    if (DidWork)
      if (DoRead)
	DidWork = velocity_vcgc[iComp].load(filename);
      else {
	DidWork = velocity_vcgc[iComp].save(filename);
	description["vel_comp"+tostr(iComp,1)]["bulk"] = filename;
      }
  }

  if (!DoRead && DidWork) {
    filename = dir + "/neutrals.json";
    DidWork = write_json(filename, description);
  }
  
  return DidWork;
}

