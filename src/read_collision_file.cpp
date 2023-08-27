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
                         Ions &ions) {

  std::string function = "read_collision_file";
  static int iFunction = -1;
  report.enter(function, iFunction);

  std::ifstream infile_ptr;
  std::string hash;
  int iErr = 0;

  report.print(1, "Reading Collision File : " + input.get_collision_file());

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

      if (report.test_verbose(2))
        std::cout << "hash : -->" << hash << "<--\n";

      // ---------------------------
      // #nu_in
      // ---------------------------

      if (hash == "#nu_in") {
        std::vector<std::vector<std::string>> csv = read_csv(infile_ptr);

        if (csv.size() > 1)
          parse_nu_in_table(csv, neutrals, ions);

        else
          std::cout << "Nu_in table is empty!!! Yikes!!!\n";
      }

      // ---------------------------
      // #resonant_nu_in
      // ---------------------------

      if (hash == "#resonant_nu_in") {
        std::vector<std::vector<std::string>> csv = read_csv(infile_ptr);

        if (csv.size() > 1)
          parse_resonant_nu_in_table(csv, neutrals, ions);

        else
          std::cout << "Resonant_nu_in table is empty!!! Yikes!!!\n";
      }

      // ---------------------------
      // #Bst for ion-ion interactions
      // ---------------------------

      if (hash == "#bst") {
        std::vector<std::vector<std::string>> csv = read_csv(infile_ptr);

        if (csv.size() > 1)
          parse_bst_in_table(csv, neutrals, ions);

        else
          std::cout << "Bst table is empty!!! Yikes!!!\n";
      }

      // ---------------------------
      // #Diff0 neutral-neutral diffusion
      // ---------------------------

      if (hash == "#diff0") {
        std::vector<std::vector<std::string>> csv = read_csv(infile_ptr);

        if (csv.size() > 1)
          parse_diff0_in_table(csv, neutrals);

        else
          std::cout << "diff0 table is empty!!! Yikes!!!\n";
      }

      // ---------------------------
      // #Diffexp neutral-neutral diffusion
      // ---------------------------

      if (hash == "#diffexp") {
        std::vector<std::vector<std::string>> csv = read_csv(infile_ptr);

        if (csv.size() > 1)
          parse_diffexp_in_table(csv, neutrals);

        else
          std::cout << "diffexp table is empty!!! Yikes!!!\n";
      }

    }

    infile_ptr.close();
    check_collision_frequncies(ions, neutrals);
  }

  report.exit(function);
  return;
}

// -----------------------------------------------------------------------------
// Check the nu_in and resonant_nu_in arrays to see if they are consistent
// -----------------------------------------------------------------------------

void check_collision_frequncies(Ions ions,
                                Neutrals neutrals) {

  // Report out the table, if verbose is high enough:

  if (report.test_verbose(2)) {
    std::cout << "nu_in table:\n";

    for (int iIon = 0; iIon < ions.nSpecies; iIon++) {
      std::cout << "====> Looking at Ion : "
                << ions.species[iIon].cName
                << " " << iIon << " of " << ions.nSpecies << "\n";

      if (ions.species[iIon].nu_ion_neutral_coef.size() > 0) {
        for (int iNeutral = 0; iNeutral < neutrals.nSpecies; iNeutral++) {
          std::cout << ions.species[iIon].cName << " -> ";
          std::cout << neutrals.species[iNeutral].cName << " = ";

          if (ions.species[iIon].nu_is_resonant[iNeutral]) {
            std::cout << "Resonant! => Checking Resonant Arrays\n";

            if (ions.species[iIon].nu_in_res_coef1.size() < 1)
              std::cout << "  --> There are no resonants for this ion!!\n";

            else {
              if (ions.species[iIon].nu_in_res_coef1[iNeutral] == 0)
                std::cout << "  --> resonant coef is 0!!\n";

              else {
                std::cout << "  temp min : "
                          << ions.species[iIon].nu_in_res_temp_min[iNeutral]
                          << "\n";
                std::cout << "  coef 1 : "
                          << ions.species[iIon].nu_in_res_coef1[iNeutral]
                          << "\n";
                std::cout << "  coef 2 : "
                          << ions.species[iIon].nu_in_res_coef2[iNeutral]
                          << "\n";
                std::cout << "  Tn ratios (Tn & Ti): "
                          << ions.species[iIon].nu_in_res_tn_frac[iNeutral]
                          << " "
                          << ions.species[iIon].nu_in_res_ti_frac[iNeutral]
                          << "\n";
              }
            }
          } else {
            std::cout << ions.species[iIon].nu_ion_neutral_coef[iNeutral]
                      << "\n";
          }
        }
      } else {
        std::cout << ions.species[iIon].cName
                  << " has no collision frequencies! \n";
      }
    }

    std::cout << "Done with check_collision_frequncies\n";
  }

  return;
}

// -----------------------------------------------------------------------------
// parse nu_in table
// -----------------------------------------------------------------------------

void parse_nu_in_table(std::vector<std::vector<std::string>> csv,
                       Neutrals &neutrals,
                       Ions &ions) {

  std::string function = "parse_nu_in_table";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int nLines = csv.size();

  // 1. check to see that we have a coefficient
  // in the last line:
  float coef = str_to_num(csv[nLines - 1][0]);

  // 2. we figure out which neutrals we have (0th col):
  int nCols = csv[0].size();
  std::vector<int> iNeutralIds_;

  for (int iCol = 1; iCol < nCols; iCol++) {
    if (report.test_verbose(4))
      std::cout << "neutral : " << csv[0][iCol] << "\n";

    iNeutralIds_.push_back(neutrals.get_species_id(csv[0][iCol]));

    if (report.test_verbose(4))
      std::cout << "iCol : " << iCol << " "
                << iNeutralIds_[iCol - 1] << "\n";
  }

  // 3. figure out which ions we have (0th column), and match
  // to neutrals:
  int iIon;

  for (int iLine = 1; iLine < nLines - 1; iLine++) {
    iIon = ions.get_species_id(csv[iLine][0]);

    if (report.test_verbose(4))
      std::cout << "iLine : " << iLine
                << " " << csv[iLine][0]
                << " " << iIon << "\n";

    // This means that we have found the ion, now let's match
    // up to the neutrals:
    if (iIon > -1) {
      // Make the array the right size, filling with zeros,
      // and setting resonant to false:
      for (int iNeutral = 0; iNeutral < neutrals.nSpecies; iNeutral++) {
        ions.species[iIon].nu_ion_neutral_coef.push_back(0.0);
        ions.species[iIon].nu_is_resonant.push_back(false);
      }

      // Now go through all of the neutrals and see which we have:
      for (int iCol = 1; iCol < nCols; iCol++) {
        // Check to see if a neutral exists:
        if (iNeutralIds_[iCol - 1] > -1) {
          // Check to see if it is supposed to be a Resonant collision:
          if (csv[iLine][iCol] == "R") {
            if (report.test_verbose(4))
              std::cout << "resonant!!! " << csv[iLine][iCol] << "\n";

            // If it is resonant, set it to -1. We will then check
            // to make sure
            ions.species[iIon].nu_is_resonant[iNeutralIds_[iCol - 1]] =
              true;
          } else {
            if (report.test_verbose(4))
              std::cout << "NONresonant!!! " << iIon << " "
                        << iNeutralIds_[iCol - 1] << " "
                        << csv[iLine][iCol] << "\n";

            ions.species[iIon].nu_ion_neutral_coef[iNeutralIds_[iCol - 1]] =
              str_to_num(csv[iLine][iCol]) * coef;
          }
        }
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
                                Ions &ions) {

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

  for (int iLine = 1; iLine < nLines - 1; iLine++) {
    iIon = ions.get_species_id(csv[iLine][0]);

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
      if (ions.species[iIon].nu_in_res_temp_min.size() < neutrals.nSpecies) {
        if (report.test_verbose(4))
          std::cout << "Creating resonant arrays\n";

        for (int iNeutral = 0; iNeutral < neutrals.nSpecies; iNeutral++) {
          ions.species[iIon].nu_in_res_temp_min.push_back(0.0);
          ions.species[iIon].nu_in_res_coef1.push_back(0.0);
          ions.species[iIon].nu_in_res_coef2.push_back(0.0);
          ions.species[iIon].nu_in_res_tn_frac.push_back(0.0);
          ions.species[iIon].nu_in_res_ti_frac.push_back(0.0);
        }
      }

      iNeutral = neutrals.get_species_id(csv[iLine][1]);

      if (iNeutral > -1) {
        if (report.test_verbose(4))
          std::cout << "Found Neutral : " << iNeutral << " "
                    << neutrals.species[iNeutral].cName << "\n";

        ions.species[iIon].nu_in_res_temp_min[iNeutral] =
          str_to_num(csv[iLine][2]);
        ions.species[iIon].nu_in_res_coef1[iNeutral] =
          str_to_num(csv[iLine][3]);
        ions.species[iIon].nu_in_res_tn_frac[iNeutral] =
          str_to_num(csv[iLine][4]);
        ions.species[iIon].nu_in_res_ti_frac[iNeutral] =
          str_to_num(csv[iLine][5]);
        ions.species[iIon].nu_in_res_coef2[iNeutral] =
          str_to_num(csv[iLine][6]);
      }
    }
  }

  report.exit(function);
  return;
}

// -----------------------------------------------------------------------------
// parse Bst table - Coulomb collision frequency coefficients
//                 - Ionospheres Book, Table 4.3
// -----------------------------------------------------------------------------

void parse_bst_in_table(std::vector<std::vector<std::string>> csv,
                        Neutrals &neutrals,
                        Ions &ions) {

  std::string function = "parse_bst_in_table";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int nLines = csv.size();  // set number of lines in table
  int nCol = csv[0].size();  // set number of columns in table
  int iIonS, iIonT, iIonP, iIonD;
  int iIon, iCol, iLine;

  std::vector<int> iIonSIds_;

  // Look for a coefficient in the last line:
  float coef = str_to_num(csv[nLines - 1][0]);

  // Read ion specie names across first row of table, cell[0][0] is empty
  for (iCol = 1; iCol < nCol; iCol++)
    iIonSIds_.push_back(ions.get_species_id(csv[0][iCol]));

  // Set the array size and fill with zeros
  for (iIon = 0; iIon < ions.nSpecies; iIon++) {
    for (iIonS = 0; iIonS < ions.nSpecies; iIonS++)
      ions.species[iIon].nu_ion_ion.push_back(0.0);
  }

  //  Read ion specie names down first column of table
  for (iLine = 1; iLine < nLines - 1; iLine++) {
    iIonT = ions.get_species_id(csv[iLine][0]);

    // Found a used specie, time to extract Bst table data
    if (iIonT > -1) {

      // Cycle through all of the species
      for (iCol = 1; iCol < nCol; iCol++) {

        // If a matching specie exists extract Bst from table
        if (iIonSIds_[iCol - 1] > -1) {
          ions.species[iIonT].nu_ion_ion[iIonSIds_[iCol - 1]] =
            str_to_num(csv[iLine][iCol]) * coef;

          if (report.test_verbose(4)) {
            std::cout << "Species s vs t : "
                      << csv[iLine][0] << " and " << csv[0][iCol] << "\n";
            std::cout << "nu_ion_ion     : "
                      << ions.species[iIonT].nu_ion_ion[iIonSIds_[iCol - 1]]
                      << "\n";
          }
        }  // End iIonSIds_
      }  // End iCol
    }  // End iIonT
  }  // End iLine

  if (report.test_verbose(4)) {
    for (iIon = 0; iIon < ions.nSpecies; iIon++) {
      for (iIonS = 0; iIonS < ions.nSpecies; iIonS++) {
        std::cout << "Bst for : " << ions.species[iIon].cName << " and "
                  << ions.species[iIonS].cName << " is "
                  << ions.species[iIon].nu_ion_ion[iIonS] << "\n";
      }
    }
  }

  report.exit(function);
  return;
} // parse_bst_in_table

void parse_diff0_in_table(std::vector<std::vector<std::string>> csv,
                          Neutrals &neutrals) {

  std::string function = "parse_diff0_in_table";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int nLines = csv.size(); // number of lines in the table
  int nCol = csv[0].size(); // number of columns in the table

  std::vector<int> iNeutralSIds_; //lists every Neutrals ID

  // last line of the csv stores the coefficient:
  float coef = str_to_num(csv[nLines - 1][0]);

  int64_t iSpecies, jSpecies, iLine, iCol;

  // Read neutral species across first row of the table:
  for (iCol = 1; iCol < nCol; iCol++)
    iNeutralSIds_.push_back(neutrals.get_species_id(csv[0][iCol]));

  // set array size and fill with zeros
  for (iSpecies = 0; iSpecies < neutrals.nSpecies; iSpecies++)
    neutrals.species[iSpecies].diff0.resize(neutrals.nSpecies, 0);

  // If the user does not know diffusion coefficients, they may simply
  // put a 0.0 in for the coefficient.  We need to figure out a better
  // way to deal with this (physics!), but for now simply grab the
  // minimum diffusion coeffient and exp and assign that to the 0s.
  precision_t minDiff0 = 1e32;
  precision_t value;

  for (iLine = 1; iLine < nLines - 1; iLine++) {

    // Get species id of neutral on line:
    int neutral_id = neutrals.get_species_id(csv[iLine][0]);

    if (neutral_id > -1) {
      for (iCol = 1; iCol < nCol; iCol++) {

        if (iNeutralSIds_[iCol - 1] > -1) {
          value = str_to_num(csv[iLine][iCol]) * coef;

          if (value > 0.0 and value < minDiff0)
            minDiff0 = value;

          neutrals.species[neutral_id].diff0[iNeutralSIds_[iCol - 1]] =
            value;

          if (neutrals.species[neutral_id].diff0[iNeutralSIds_[iCol - 1]] > 0)
            if (report.test_verbose(4)) {
              std::cout << "Two diffusion species: "
                        << csv[iLine][0] << " and " << csv[0][iCol] << "\n";
              std::cout << "diff0     : "
                        << neutrals.species[neutral_id].diff0[iNeutralSIds_[iCol - 1]] << "\n";
            } //verbose

        } // if matching specie exist
      } // for iCol
    } // if neutral_id > -1
  } // for iLine

  // This takes all of the zeros and turns them into at least something
  // so that the code doesn't blow up
  for (iSpecies = 0; iSpecies < neutrals.nSpecies; iSpecies++)
    for (jSpecies = 0; jSpecies < neutrals.nSpecies; jSpecies++)
      if (neutrals.species[iSpecies].diff0[jSpecies] <= 0.0)
        neutrals.species[iSpecies].diff0[jSpecies] = minDiff0;

  if (report.test_verbose(4)) {
    for (int i = 0; i < neutrals.nSpecies; i++) {
      cout << "Current species: " << neutrals.species[i].cName << endl;

      for (int j = 0; j < neutrals.species[i].diff0.size(); j++)
        cout << neutrals.species[i].diff0[j] << ", ";

      cout << endl;
    }
  }

  report.exit(function);
  return;
} // parse_diff0_in_table

void parse_diffexp_in_table(std::vector<std::vector<std::string>> csv,
                            Neutrals &neutrals) {

  std::string function = "parse_diffexp_in_table";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int nLines = csv.size(); // number of lines in the table
  int nCol = csv[0].size(); // number of columns in the table

  std::vector<int> iNeutralSIds_; //lists every Neutrals ID

  int64_t iSpecies, jSpecies, iLine, iCol;

  // Read neutral species across first row of the table:
  for (iCol = 1; iCol < nCol; iCol++)
    iNeutralSIds_.push_back(neutrals.get_species_id(csv[0][iCol]));

  // set array size and fill with zeros
  for (iSpecies = 0; iSpecies < neutrals.nSpecies; iSpecies++)
    neutrals.species[iSpecies].diff_exp.resize(neutrals.nSpecies, 0);

  // If the user does not know diffusion coefficients, they may simply
  // put a 0.0 in for the coefficient.  We need to figure out a better
  // way to deal with this (physics!), but for now simply grab the
  // minimum diffusion coeffient and exp and assign that to the 0s.
  precision_t minDiffexp = 1e32;
  precision_t value;

  for (iLine = 1; iLine < nLines - 1; iLine++) {

    // Get species id of neutral on line:
    int neutral_id = neutrals.get_species_id(csv[iLine][0]);

    if (neutral_id > -1) {
      for (iCol = 1; iCol < nCol; iCol++) {

        if (iNeutralSIds_[iCol - 1] > -1) {
          value = str_to_num(csv[iLine][iCol]);

          if (value > 0.0 and value < minDiffexp)
            minDiffexp = value;

          neutrals.species[neutral_id].diff_exp[iNeutralSIds_[iCol - 1]] =
            value;

          if (neutrals.species[neutral_id].diff_exp[iNeutralSIds_[iCol - 1]] > 0)
            if (report.test_verbose(4)) {
              std::cout << "Two diffusion species: "
                        << csv[iLine][0] << " and " << csv[0][iCol] << "\n";
              std::cout << "diff_exp     : "
                        << neutrals.species[neutral_id].diff_exp[iNeutralSIds_[iCol - 1]] << "\n";
            } //verbose

        } // if matching specie exist
      } // for iCol
    } // if neutral_id > -1
  } // for iLine

  // This takes all of the zeros and turns them into at least something
  // so that the code doesn't blow up
  for (iSpecies = 0; iSpecies < neutrals.nSpecies; iSpecies++)
    for (jSpecies = 0; jSpecies < neutrals.nSpecies; jSpecies++)
      if (neutrals.species[iSpecies].diff_exp[jSpecies] <= 0.0)
        neutrals.species[iSpecies].diff_exp[jSpecies] = minDiffexp;

  if (report.test_verbose(4)) {
    for (int i = 0; i < neutrals.nSpecies; i++) {
      cout << "Current species: " << neutrals.species[i].cName << endl;

      for (int j = 0; j < neutrals.species[i].diff_exp.size(); j++)
        cout << neutrals.species[i].diff_exp[j] << ", ";

      cout << endl;
    }
  }

  report.exit(function);
  return;
} // parse_diffexp_in_table

