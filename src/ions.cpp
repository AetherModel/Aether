// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>
#include <fstream>

#include "../include/aether.h"

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
  tmp.density_scgc.ones();
  tmp.temperature_scgc.set_size(nLons, nLats, nAlts);
  tmp.temperature_scgc.ones();
  tmp.ionization_scgc.set_size(nLons, nLats, nAlts);
  tmp.ionization_scgc.zeros();

  tmp.sources_scgc.set_size(nLons, nLats, nAlts);
  tmp.sources_scgc.zeros();
  tmp.losses_scgc.set_size(nLons, nLats, nAlts);
  tmp.losses_scgc.zeros();

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

  for (int iSpecies=0; iSpecies < nIons; iSpecies++) {
    tmp = create_species(grid);
    species.push_back(tmp);
  }

  // Create one extra species for electrons
  tmp = create_species(grid);
  species.push_back(tmp);

  // State variables:

  density_scgc.set_size(nLons, nLats, nAlts);
  density_scgc.ones();
  ion_temperature_scgc.set_size(nLons, nLats, nAlts);
  ion_temperature_scgc.ones();
  electron_temperature_scgc.set_size(nLons, nLats, nAlts);
  electron_temperature_scgc.ones();

  tmp.sources_scgc.set_size(nLons, nLats, nAlts);
  tmp.sources_scgc.zeros();
  tmp.losses_scgc.set_size(nLons, nLats, nAlts);
  tmp.losses_scgc.zeros();

  // This gets a bunch of the species-dependent characteristics:
  int iErr = read_planet_file(input, report);
  if (iErr > 0) std::cout << "Error in reading planet file!" << '\n';
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
        if (lines.size()-1 != nIons) {
          std::cout << "num of ion species (nIons) defined in planet.h file : "
                    << nIons << "\n";
          std::cout << "number of ions defined in planet.in file : "
                    << lines.size() << "\n";
          std::cout << "These don't match!\n";
          iErr = 1;
        } else {
          // assume order of rows right now:
          // name, mass, charge, advect
          for (int iSpecies=0; iSpecies < nIons; iSpecies++) {
            report.print(5, "setting ion species " + lines[iSpecies+1][0]);
            species[iSpecies].cName = lines[iSpecies+1][0];
            species[iSpecies].mass = stof(lines[iSpecies+1][1])*amu;
            species[iSpecies].charge = stoi(lines[iSpecies+1][2]);
            species[iSpecies].DoAdvect = stoi(lines[iSpecies+1][3]);
          }
          species[nIons].cName = "e-";
          species[nIons].mass = mass_electron;
          species[nIons].charge = -1;
          species[nIons].DoAdvect = 0;
        }
      }
      if (infile_ptr.eof()) IsDone = 1;
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
  for (iSpecies=0; iSpecies < nIons; iSpecies++)
    species[nIons].density_scgc =
      species[nIons].density_scgc + species[iSpecies].density_scgc;
  density_scgc = species[nIons].density_scgc;

  report.exit(function);
  return;
}
