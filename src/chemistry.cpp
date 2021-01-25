// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <string>
#include <fstream>
#include <vector>
#include <iostream>

#include "../include/chemistry.h"
#include "../include/inputs.h"
#include "../include/report.h"
#include "../include/file_input.h"

// -----------------------------------------------------------------------------
// Initialize chemistry class
// -----------------------------------------------------------------------------

Chemistry::Chemistry(Neutrals neutrals,
                     Ions ions,
                     Inputs args,
                     Report &report) {

  std::string function = "Chemistry::Chemistry";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int iErr = 0;

  read_chemistry_file(neutrals, ions, args, report);

  report.exit(function);
  return;
}

// -----------------------------------------------------------------------------
// Read chemistry file
// -----------------------------------------------------------------------------

int Chemistry::read_chemistry_file(Neutrals neutrals,
                                   Ions ions,
                                   Inputs args,
                                   Report &report) {

  std::string function = "Chemistry::read_chemistry_file";
  static int iFunction = -1;
  report.enter(function, iFunction);

  std::ifstream infile_ptr;
  int iErr = 0;
  reaction_type reaction;

  report.print(1, "Reading Chemistry File : "+args.get_chemistry_file());

  infile_ptr.open(args.get_chemistry_file());

  if (!infile_ptr.is_open()) {
    std::cout << "Could not open chemistry file!\n";
    iErr = 1;
  } else {

    if (infile_ptr.good()) {

      std::vector<std::vector<std::string>> csv = read_csv(infile_ptr);

      int nLines = csv.size();

      if (nLines <= 2) {
        iErr = 1;
      } else {

        nReactions = 0;

        // Skip 2 lines of headers!
        for (int iLine = 2; iLine < nLines; iLine++) {
          // Some final rows can have comments in them, so we want to
          // skip anything where the length of the string in column 2
          // is == 0:
          if (csv[iLine][1].length() > 0) {
	    report.print(2, "interpreting chemistry line : "+csv[iLine][0]);
            reaction = interpret_reaction_line(neutrals, ions,
                                               csv[iLine], report);
            if (reaction.nLosses > 0 && reaction.nSources > 0) {
              reactions.push_back(reaction);
              nReactions++;
            }
          }
        }
      }
    }
  }
  report.exit(function);
  return iErr;
}

// -----------------------------------------------------------------------------
// Interpret a comma separated line of the chemical reaction file
// -----------------------------------------------------------------------------

Chemistry::reaction_type Chemistry::interpret_reaction_line(Neutrals neutrals,
                                                            Ions ions,
                                                            std::vector<std::string> line,
                                                            Report &report) {

  std::string function = "Chemistry::interpret_reaction_line";
  static int iFunction = -1;
  report.enter(function, iFunction);

  reaction_type reaction;

  int i;
  int id_, IsNeutral;

  // Losses (left side) first:
  reaction.nLosses = 0;
  for (i = 0; i < 3; i++) {
    find_species_id(line[i], neutrals, ions, id_, IsNeutral, report);
    if (id_ >= 0) {
      reaction.losses_names.push_back(line[i]);
      reaction.losses_ids.push_back(id_);
      reaction.losses_IsNeutral.push_back(IsNeutral);
      reaction.nLosses++;
    }
  }

  // Sources (right side) second:
  reaction.nSources = 0;
  for (i = 4; i < 7; i++) {
    find_species_id(line[i], neutrals, ions, id_, IsNeutral, report);
    if (id_ >= 0) {
      reaction.sources_names.push_back(line[i]);
      reaction.sources_ids.push_back(id_);
      reaction.sources_IsNeutral.push_back(IsNeutral);
      reaction.nSources++;
    }
  }

  // Reaction Rate:
  reaction.rate = stof(line[7]);

  // Branching Ratio:
  reaction.branching_ratio = stof(line[10]);

  // energy released as exo-thermic reaction:
  reaction.energy = stof(line[11]);

  report.exit(function);
  return reaction;
}

// -----------------------------------------------------------------------------
// Match a string to the neutral or ion species
// -----------------------------------------------------------------------------

void Chemistry::find_species_id(std::string name,
                                Neutrals neutrals,
                                Ions ions,
                                int &id_,
                                int &IsNeutral,
                                Report &report) {

  std::string function = "Chemistry::find_species_id";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int iSpecies;

  id_ = -1;
  IsNeutral = -1;
  if (name.length() > 0) {

    // Check Neutrals:
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      if (name == neutrals.neutrals[iSpecies].cName) {
        id_ = iSpecies;
        IsNeutral = 1;
        break;
      }

    if (id_ == -1) {
      // Check Ions:
      for (iSpecies = 0; iSpecies <= nIons; iSpecies++)
        if (name == ions.species[iSpecies].cName) {
          id_ = iSpecies;
          IsNeutral = 0;
          break;
        }
    }
  }

  report.exit(function);
  return;
}

// -----------------------------------------------------------------------------
// Display a reaction:
// -----------------------------------------------------------------------------

void Chemistry::display_reaction(Chemistry::reaction_type reaction) {

  int i;

  std::cout << "Number of Losses : " << reaction.nLosses << "\n";
  std::cout << "Number of Sources : " << reaction.nSources << "\n";

  for (i = 0; i < reaction.nLosses; i++)
    std::cout << reaction.losses_names[i] << " + ";
  std::cout << " -> ";
  for (i = 0; i < reaction.nSources; i++)
    std::cout << reaction.sources_names[i] << " + ";
  std::cout << " ( RR : " << reaction.rate << ")\n";

  for (i = 0; i < reaction.nLosses; i++)
    std::cout << reaction.losses_ids[i]
              << "(" << reaction.losses_IsNeutral[i] << ")" << " + ";
  std::cout << " -> ";
  for (i = 0; i < reaction.nSources; i++)
    std::cout << reaction.sources_ids[i]
              << "(" << reaction.sources_IsNeutral[i]
              << ")" << " + ";
  std::cout << " ( RR : " << reaction.rate << ")\n";
}
