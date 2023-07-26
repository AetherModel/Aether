// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_CHEMISTRY_H_
#define INCLUDE_CHEMISTRY_H_

#include <vector>
#include <string>

/**************************************************************
 * \class Chemistry
 *
 * \brief For keeping track of chemical reactions and conducting chemistry
 *
 * This class keeps track of the chemical reactions
 *
 *
 *
 * \author Aaron Ridley, Michael Rinaldi
 * \date 2022/07/05, 2023/03/10
 *
 **************************************************************/

class Chemistry {

 public:

  /// This structure keeps track of one reaction.

  struct reaction_type {
    // Reactions:
    // loss1 + loss2 + loss3 -> source1 + source2 + source3
    // Need to match the CSV file of the species or won't be recognized.
    std::vector<std::string> sources_names;
    std::vector<std::string> losses_names;

    /// These are sources species indices:
    std::vector<int> sources_ids;
    /// Determine whether source is a neutral or not:
    std::vector<bool> sources_IsNeutral;

    /// These are losses species indices:
    std::vector<int> losses_ids;
    /// Determine whether loss is a neutral or not:
    std::vector<bool> losses_IsNeutral;

    /// Keep track of the number of sources:
    int nSources;
    /// Keep track of the number of losses:
    int nLosses;

    /// Energy for exothermic reactions:
    precision_t energy;
    /// Reaction rate:
    precision_t rate;
    /// Branching ratio if needed:
    precision_t branching_ratio;

    // variables for temperature-dependent rate calculation
    /// numerator of temperature dependent rate:
    precision_t numerator;
    /// denominator of temperature dependent rate:
    std::string denominator;
    /// exponent of (num / den) ^ exp:
    precision_t exponent;
    /// which temperature to use w/ min/max:
    std::string piecewiseVar;
    /// temperature min and max within which to use the reaction:
    int min;
    /// temperature min and max within which to use the reaction:
    int max;
    /// type of formula to use for reaction rate:
    int type;
    /// name of the reaction
    std::string name;    

  };

  /// vector to keep track of all reactions:
  std::vector<reaction_type> reactions;
  /// Total number of reactions:
  int64_t nReactions;

  Chemistry(Neutrals neutrals,
            Ions ions);

  void calc_chemistry(Neutrals &neutrals,
                      Ions &ions,
                      Times time,
                      Grid grid);

  void calc_chemical_sources(Neutrals &neutrals,
                             Ions &ions);

 private:

  int read_chemistry_file(Neutrals neutrals,
                          Ions ions);

  reaction_type interpret_reaction_line(Neutrals neutrals,
                                        Ions ions,
                                        std::vector<std::string> line,
					json headers);

  void find_species_id(std::string name,
                       Neutrals neutrals,
                       Ions ions,
                       int &id_,
                       bool &IsNeutral);

  void display_reaction(reaction_type reaction);
};

#endif  // INCLUDE_CHEMISTRY_H_
