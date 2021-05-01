// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iostream>

#include "../include/aether.h"

// -----------------------------------------------------------------------------
// read collision frequencies and diffusion coefficient file
// -----------------------------------------------------------------------------

void read_collision_file(Neutrals &neutrals,
			 Ions &ions,
			 Inputs input,
			 Report &report) {

  std::string function = "read_collision_file";
  static int iFunction = -1;
  report.enter(function, iFunction);

  std::ifstream infile_ptr;
  std::string hash;
  int iErr = 0;

  report.print(1, "Reading Collision File : "+input.get_collision_file());

  infile_ptr.open(input.get_collision_file());

  if (!infile_ptr.is_open()) {
    std::cout << "Could not open collision file!\n";
    iErr = 1;
  } else {

    while (!infile_ptr.eof()) {

      // ---------------------------
      // Find the next hash:
      // ---------------------------

      hash = find_next_hash(infile_ptr);
      if (report.test_verbose(0))
        std::cout << "hash : -->" << hash << "<--\n";

      // ---------------------------
      // #nu_in
      // ---------------------------

      if (hash == "#nu_in") {
	std::cout << "Reading #nu_in\n";
        std::vector<std::vector<std::string>> csv = read_csv(infile_ptr);
        if (csv.size() > 1) {
	  parse_nu_in_table(csv, neutrals, ions, report);
	} else {
	  std::cout << "NU_in table is empty!!! Yikes!!!\n";
	}
      }

      // ---------------------------
      // #resonant_nu_in
      // ---------------------------

      if (hash == "#resonant_nu_in") {
	std::cout << "Reading #resonant_nu_in\n";
        std::vector<std::vector<std::string>> csv = read_csv(infile_ptr);
        if (csv.size() > 1) {
	  parse_resonant_nu_in_table(csv, neutrals, ions, report);
	} else {
	  std::cout << "Resonant_nu_in table is empty!!! Yikes!!!\n";
	}
      }



    }
    infile_ptr.close();
    // We should add something here to check the tables:
    // - Need to see if the nu_in exists for each ion species
    // - Need to make sure that if there is a resonant Nu, it exists

  }
  report.exit(function);
  return;
}

// -----------------------------------------------------------------------------
// parse nu_in table
// -----------------------------------------------------------------------------

void parse_nu_in_table(std::vector<std::vector<std::string>> csv,
		       Neutrals &neutrals,
		       Ions &ions,
		       Report &report) {

  std::string function = "parse_nu_in_table";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int nLines = csv.size();
  
  // 1. check to see that we have a coefficient
  // in the last line:
  float coef = stof(csv[nLines-1][0]);

  // 2. we figure out which neutrals we have (0th col):
  int nCols = csv[0].size();
  std::vector<int> iNeutralIds_;
  for (int iCol = 1; iCol < nCols; iCol++) {
    if (report.test_verbose(4))
      std::cout << "neutral : " << csv[0][iCol] << "\n";
    iNeutralIds_.push_back(neutrals.get_species_id(csv[0][iCol],
						   report));
    if (report.test_verbose(4))
      std::cout << "iCol : " << iCol << " "
		<< iNeutralIds_[iCol-1] << "\n";
  }

  // 3. figure out which ions we have (0th column), and match
  // to neutrals:
  int iIon;
  for (int iLine = 1; iLine < nLines-1; iLine++) {
    iIon = ions.get_species_id(csv[iLine][0], report);
    if (report.test_verbose(4))
      std::cout << "iLine : " << iLine
		<< " " << csv[iLine][0]
		<< " " << iIon << "\n";
    
    // This means that we have found the ion, now let's match
    // up to the neutrals:
    if (iIon > -1) {
      // Make the array the right size, filling with zeros,
      // and setting resonant to false:
      for (int iNeutral = 0; iNeutral < nSpecies; iNeutral++) {
	ions.species[iIon].nu_ion_neutral.push_back(0.0);
	ions.species[iIon].nu_is_resonant.push_back(false);
      }
      // Now go through all of the neutrals and see which we have:
      for (int iCol = 1; iCol < nCols; iCol++) {
	// Check to see if a neutral exists: 
	if (iNeutralIds_[iCol-1] > -1) {
	  // Check to see if it is supposed to be a Resonant collision:
	  if (csv[iLine][iCol] == "R") {
	    if (report.test_verbose(4))
	      std::cout << "resonant!!! " << csv[iLine][iCol] << "\n";
	    // If it is resonant, set it to -1. We will then check
	    // to make sure 
	    ions.species[iIon].nu_is_resonant[iNeutralIds_[iCol-1]] =
	      true;
	  } else {
	    if (report.test_verbose(4))
	      std::cout << "NONresonant!!! " << iIon << " "
			<< iNeutralIds_[iCol-1] << " "
			<< csv[iLine][iCol] << "\n";
	    ions.species[iIon].nu_ion_neutral[iNeutralIds_[iCol-1]] =
	      stof(csv[iLine][iCol]) * coef;
	  }
	}
      }
    }
  }

  // Report out the table, if verbose is high enough:
  
  if (report.test_verbose(0)) {
    std::cout << "nu_in table:\n";
    for (int iIon = 0; iIon < nIons; iIon++) {
      if (ions.species[iIon].nu_ion_neutral.size() > 0) {
	for (int iNeutral = 0; iNeutral < nSpecies; iNeutral++) {
	  std::cout << ions.species[iIon].cName << " -> ";
	  std::cout << neutrals.species[iNeutral].cName << " = ";
	  if (ions.species[iIon].nu_is_resonant[iNeutral]) {
	    std::cout << "Resonant!\n";
	  } else {
	    std::cout << ions.species[iIon].nu_ion_neutral[iNeutral] << "\n";
	  }
	}
      } else {
	std::cout << ions.species[iIon].cName
		  << " has no collision frequencies! \n";
      }
    }
  }
  
  report.exit(function);
  return;
}

// -----------------------------------------------------------------------------
// parse resonant_nu_in table
// -----------------------------------------------------------------------------

void parse_resonant_nu_in_table(std::vector<std::vector<std::string>> csv,
				Neutrals &neutrals,
				Ions &ions,
				Report &report) {

  std::string function = "parse_resonant_nu_in_table";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int nLines = csv.size();
  
  // Go row-by-row to:
  //   1. See which ion we are working with
  //   2a. Check to see if we have already created the needed arrays
  //       (This can happen if the ion has multiple resonant collision freq.)
  //   2b. If we haven't, create them
  //   3. Find the neutral
  //   4. If the neutral exists, fill in the values

  int iIon, iNeutral;
  for (int iLine = 1; iLine < nLines-1; iLine++) {
    iIon = ions.get_species_id(csv[iLine][0], report);
    if (report.test_verbose(4))
      std::cout << "iLine : " << iLine
		<< " " << csv[iLine][0]
		<< " " << iIon << "\n";
    // This means that we have found the ion, now let's match
    // up to the neutrals:
    if (iIon > -1) {
      if (report.test_verbose(4))
	std::cout << "Found Ion : " << ions.species[iIon].cName << "\n";
      // Make the array the right size, filling with zeros,
      // and setting resonant to false:
      if (ions.species[iIon].nu_in_res_temp_min.size() < nSpecies) {
	if (report.test_verbose(4))
	  std::cout << "Creating resonant arrays\n";
	for (int iNeutral = 0; iNeutral < nSpecies; iNeutral++) {
	  ions.species[iIon].nu_in_res_temp_min.push_back(0.0);
	  ions.species[iIon].nu_in_res_coef1.push_back(0.0);
	  ions.species[iIon].nu_in_res_coef2.push_back(0.0);
	  ions.species[iIon].nu_in_res_tn_frac.push_back(0.0);
	  ions.species[iIon].nu_in_res_ti_frac.push_back(0.0);
	}
      }

      iNeutral = neutrals.get_species_id(csv[iLine][1], report);
      if (iNeutral > -1) {
	if (report.test_verbose(4))
	  std::cout << "Found Neutral : " << iNeutral << " "
		    << neutrals.species[iNeutral].cName << "\n";	  
	ions.species[iIon].nu_in_res_temp_min[iNeutral] = stof(csv[iLine][2]);
	ions.species[iIon].nu_in_res_coef1[iNeutral] = stof(csv[iLine][3]);
	ions.species[iIon].nu_in_res_tn_frac[iNeutral] = stof(csv[iLine][4]);
	ions.species[iIon].nu_in_res_ti_frac[iNeutral] = stof(csv[iLine][5]);
	ions.species[iIon].nu_in_res_coef2[iNeutral] = stof(csv[iLine][6]);
      }
    }
  }
  return;
}
    
