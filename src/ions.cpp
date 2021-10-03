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

  if (input.get_do_restart()) {
    report.print(1, "Restarting! Reading ion files!");
    bool DidWork = restart_file(input.get_restartin_dir(), DoRead);
    if (!DidWork)
      std::cout << "Reading Restart for Ions Failed!!!\n";
  }
  
  // This gets a bunch of the species-dependent characteristics:
  int iErr = read_planet_file(input, report);

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
  json description;

  for (int iSpecies = 0; iSpecies < nIons; iSpecies++) {

    // Output Densities
    filename = dir + "/ion_s"+tostr(iSpecies,2)+"_n.bin";
    if (DidWork)
      if (DoRead)
	DidWork = species[iSpecies].density_scgc.load(filename);
      else {
	DidWork = species[iSpecies].density_scgc.save(filename);	
	description["density"][species[iSpecies].cName] = filename;
      }

    // Output Temperature for each species
    filename = dir + "/ion_s"+tostr(iSpecies,2)+"_t.bin";
    if (DidWork)
      if (DoRead)
	DidWork = species[iSpecies].temperature_scgc.load(filename);
      else {
	DidWork = species[iSpecies].temperature_scgc.save(filename);
	description["temperature"][species[iSpecies].cName] = filename;
      }
    // Output Velocity (Parallel and Perp) for each species
    for (int iComp = 0; iComp < 3; iComp++) {
      filename = dir + "/ion_s"+tostr(iSpecies,2)+"_vpa"+tostr(iComp,1)+".bin";
      if (DidWork)
	if (DoRead)
	  DidWork = species[iSpecies].par_velocity_vcgc[iComp].load(filename);
	else {
	  DidWork = species[iSpecies].par_velocity_vcgc[iComp].save(filename);
	  description["vel_par_comp"+tostr(iComp,1)][species[iSpecies].cName] =
	    filename;
	}
      filename = dir + "/ion_s"+tostr(iSpecies,2)+"_vpe"+tostr(iComp,1)+".bin";
      if (DidWork)
	if (DoRead)
	  DidWork = species[iSpecies].perp_velocity_vcgc[iComp].load(filename);
	else {
	  DidWork = species[iSpecies].perp_velocity_vcgc[iComp].save(filename);
	  description["vel_perp_comp"+tostr(iComp,1)][species[iSpecies].cName] =
	    filename;
	}
    }
  }

  // Output bulk Temperature
  filename = dir + "/ion_t.bin";
  if (DidWork)
    if (DoRead)
      DidWork = ion_temperature_scgc.load(filename);
    else {
      DidWork = ion_temperature_scgc.save(filename);
      description["temperature"]["bulk"] = filename;
    }

  // Output electron Temperature
  filename = dir + "/ele_t.bin";
  if (DidWork)
    if (DoRead)
      DidWork = electron_temperature_scgc.load(filename);
    else {
      DidWork = electron_temperature_scgc.save(filename);
      description["temperature"]["electron"] = filename;
    }
  
  if (!DoRead && DidWork) {
    filename = dir + "/ions.json";
    DidWork = write_json(filename, description);
  }
  
  return DidWork;
}

