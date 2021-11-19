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
  tmp.temperature_scgc.ones();
  tmp.ionization_scgc.set_size(nLons, nLats, nAlts);
  tmp.ionization_scgc.zeros();

  tmp.sources_scgc.set_size(nLons, nLats, nAlts);
  tmp.sources_scgc.zeros();
  tmp.losses_scgc.set_size(nLons, nLats, nAlts);
  tmp.losses_scgc.zeros();

  tmp.par_velocity_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);
  tmp.perp_velocity_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);

  tmp.nu_ion_neutral_vcgc = make_cube_vector(nLons, nLats, nAlts, nSpecies);

  return tmp;
}

// -----------------------------------------------------------------------------
//  Initialize Ions class
// -----------------------------------------------------------------------------

Ions::Ions(Grid grid, Inputs input, Report report) {

  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();

  species_chars tmp;

  report.print(2, "Initializing Ions");

  for (int iSpecies = 0; iSpecies < nIons; iSpecies++) {
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
  ion_temperature_scgc.set_size(nLons, nLats, nAlts);
  ion_temperature_scgc.ones();
  electron_temperature_scgc.set_size(nLons, nLats, nAlts);
  electron_temperature_scgc.ones();

  tmp.sources_scgc.set_size(nLons, nLats, nAlts);
  tmp.sources_scgc.zeros();
  tmp.losses_scgc.set_size(nLons, nLats, nAlts);
  tmp.losses_scgc.zeros();

  potential_scgc.set_size(nLons, nLats, nAlts);
  potential_scgc.zeros();

  eflux.set_size(nLons, nLats);
  eflux.zeros();
  avee.set_size(nLons, nLats);
  avee.zeros();

  // Make vectors:
  efield_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);
  exb_vcgc = make_cube_vector(nLons, nLats, nAlts, 3);

  // This gets a bunch of the species-dependent characteristics:
  int iErr = read_planet_file(input, report);

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

int Ions::read_planet_file(Inputs input, Report report) {

  int iErr = 0;
  std::string hash;
  std::ifstream infile_ptr;

  report.print(3, "In read_planet_file for Ions");

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

      if (hash == "#ions") {

        // Read in the characteristics as CSVs:
        report.print(4, "Found #ions!");

        std::vector<std::vector<std::string>> lines = read_csv(infile_ptr);

        if (lines.size() - 1 != nIons) {
          std::cout << "num of ion species (nIons) defined in planet.h file : "
                    << nIons << "\n";
          std::cout << "number of ions defined in planet.in file : "
                    << lines.size() << "\n";
          std::cout << "These don't match!\n";
          iErr = 1;
        } else {
          // assume order of rows right now:
          // name, mass, charge, advect
          for (int iSpecies = 0; iSpecies < nIons; iSpecies++) {
            report.print(5, "setting ion species " + lines[iSpecies + 1][0]);
            species[iSpecies].cName = lines[iSpecies + 1][0];
            species[iSpecies].mass = stof(lines[iSpecies + 1][1]) * cAMU;
            species[iSpecies].charge = stoi(lines[iSpecies + 1][2]);
            species[iSpecies].DoAdvect = stoi(lines[iSpecies + 1][3]);
          }

          species[nIons].cName = "e-";
          species[nIons].mass = cME;
          species[nIons].charge = -1;
          species[nIons].DoAdvect = 0;
        }
      }

      if (infile_ptr.eof())
        IsDone = 1;
    }

    infile_ptr.close();
  }

  return iErr;
}

// -----------------------------------------------------------------------------
// Calculate the electron density from the sum of all ion species
// -----------------------------------------------------------------------------

void Ions::fill_electrons(Report &report) {

  int iSpecies;

  std::string function = "Ions::fill_electrons";
  static int iFunction = -1;
  report.enter(function, iFunction);

  species[nIons].density_scgc.zeros();

  for (iSpecies = 0; iSpecies < nIons; iSpecies++)
    species[nIons].density_scgc =
      species[nIons].density_scgc + species[iSpecies].density_scgc;

  density_scgc = species[nIons].density_scgc;

  report.exit(function);
  return;
}

//----------------------------------------------------------------------
// return the index of the requested species
// This will return -1 if the species is not found or name is empty
// Will return nIons for electrons
//----------------------------------------------------------------------

int Ions::get_species_id(std::string name, Report &report) {

  std::string function = "Ions::get_species_id";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int iSpecies;
  int id_ = -1;

  if (name.length() > 0) {
    for (iSpecies = 0; iSpecies <= nIons; iSpecies++)
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
  RestartContainer.set_netcdf();
  RestartContainer.set_directory(dir);
  RestartContainer.set_filename("ions");

  try {
    if (DoRead)
      RestartContainer.read_container_netcdf();
    else {
      RestartContainer.set_version(0.1);
      RestartContainer.set_time(0.0);
    }

    std::string cName;
    for (int iSpecies = 0; iSpecies < nIons; iSpecies++) { 
      // ----------------------------
      // Density and Temperature (per ion)
      // ----------------------------
      cName = species[iSpecies].cName;
      if (DoRead) {
	iVar = RestartContainer.find_variable(cName);
	species[iSpecies].density_scgc =
	  RestartContainer.get_element_value(iVar);
      } else {
	RestartContainer.store_variable(cName,
					density_unit,
					species[iSpecies].density_scgc);
      }
      cName = temperature_name + " ("+species[iSpecies].cName+")";
      if (DoRead) {
	iVar = RestartContainer.find_variable(cName);
	species[iSpecies].temperature_scgc =
	  RestartContainer.get_element_value(iVar);
      } else {
	RestartContainer.store_variable(cName,
					temperature_unit,
					species[iSpecies].temperature_scgc);
      }
      // ----------------------------
      // Parallel Velocity (per ion)
      // ----------------------------
      for (int iDir = 0; iDir < 3; iDir++) { 
	cName = par_velocity_name[iDir] + " ("+species[iSpecies].cName+")";
	if (DoRead) {
	  iVar = RestartContainer.find_variable(cName);
	  species[iSpecies].par_velocity_vcgc[iDir] =
	    RestartContainer.get_element_value(iVar);
	} else {
	  RestartContainer.store_variable(cName,
					  velocity_unit,
					  species[iSpecies].
					  par_velocity_vcgc[iDir]);
	}
      }
      // ----------------------------
      // Perpendicular Velocity (per ion)
      // ----------------------------
      for (int iDir = 0; iDir < 3; iDir++) {
	cName = perp_velocity_name[iDir] + " ("+species[iSpecies].cName+")";
	if (DoRead) {
	  iVar = RestartContainer.find_variable(cName);
	  species[iSpecies].perp_velocity_vcgc[iDir] =
	    RestartContainer.get_element_value(iVar);
	  iVar++;
	} else {
	  RestartContainer.store_variable(cName,
					  velocity_unit,
					  species[iSpecies].
					  perp_velocity_vcgc[iDir]);
	}
      }
    }
    // ----------------------------
    // Bulk ion and electron temperatures
    // ----------------------------
    cName = temperature_name + " (bulk ion)";
    if (DoRead) {
      iVar = RestartContainer.find_variable(cName);
      ion_temperature_scgc = RestartContainer.get_element_value(iVar);
    } else {
      RestartContainer.store_variable(cName,
				      temperature_unit,
				      ion_temperature_scgc);
    }
  
    cName = temperature_name + " (electron)";
    if (DoRead) {
      iVar = RestartContainer.find_variable(cName);
      electron_temperature_scgc = RestartContainer.get_element_value(iVar);
    } else {
      RestartContainer.store_variable(temperature_name + " (electron)",
				      temperature_unit,
				      electron_temperature_scgc);
      RestartContainer.write();
    }
    RestartContainer.clear_variables();
  } catch (...) {
    std::cout << "Error reading in ion restart file!\n";
    DidWork = false;
  }
  return DidWork;
}

