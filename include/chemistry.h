// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_CHEMISTRY_H_
#define INCLUDE_CHEMISTRY_H_

#include <vector>
#include <string>

class Chemistry {

 public:

  struct reaction_type {
    // Reactions:
    // loss1 + loss2 + loss3 -> source1 + source2 + source3
    std::vector<std::string> sources_names;
    std::vector<std::string> losses_names;

    std::vector<int> sources_ids;
    std::vector<int> sources_IsNeutral;

    std::vector<int> losses_ids;
    std::vector<int> losses_IsNeutral;

    int nSources;
    int nLosses;

    precision_t energy;
    precision_t rate;
    precision_t branching_ratio;
  };

  std::vector<reaction_type> reactions;
  int64_t nReactions;

  Chemistry(Neutrals neutrals,
            Ions ions,
            Inputs args,
            Report &report);

  void calc_chemistry(Neutrals &neutrals,
                      Ions &ions,
                      Times time,
                      Grid grid,
                      Report &report);

  void calc_chemical_sources(Neutrals &neutrals,
                             Ions &ions,
                             Report &report);

 private:

  struct sources_and_losses_type {

    precision_t neutral_sources[nSpecies];
    precision_t neutral_losses[nSpecies];
    precision_t ion_sources[nIons];
    precision_t ion_losses[nIons];

    precision_t heat_neutrals;
    precision_t heat_ions;
    precision_t heat_electrons;
  };

  sources_and_losses_type sources_and_losses;

  int read_chemistry_file(Neutrals neutrals,
                          Ions ions,
                          Inputs args,
                          Report &report);

  reaction_type interpret_reaction_line(Neutrals neutrals,
                                        Ions ions,
                                        std::vector<std::string> line,
                                        Report &report);

  void find_species_id(std::string name,
                       Neutrals neutrals,
                       Ions ions,
                       int &id_,
                       int &IsNeutral,
                       Report &report);

  void display_reaction(reaction_type reaction);
};

#endif  // INCLUDE_CHEMISTRY_H_
