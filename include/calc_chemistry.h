// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_CALC_CHEMISTRY_H_
#define INCLUDE_CALC_CHEMISTRY_H_

#include <vector>
#include <string>
#include "../include/ions.h"
#include "../include/neutrals.h"
#include "../include/times.h"
#include "../include/grid.h"
#include "../include/report.h"

class Chemistry {

 public:

  Chemistry();

  struct reaction_type {
    // Reactions:
    // loss1 + loss2 + loss3 -> source1 + source2 + source3

    std::vector<std::string> sources_names;
    std::vector<std::string> losses_names;

    std::vector<int> sources_ids;
    std::vector<int> sources_IsNeutrals;

    std::vector<int> losses_ids;
    std::vector<int> losses_IsNeutrals;

    precision_t rate;
  };

  std::vector<reaction_type> reactions;

  struct sources_and_losses_type {
    precision_t neutral_sources[nSpecies];
    precision_t neutral_losses[nSpecies];
    precision_t ion_sources[nIons];
    precision_t ion_losses[nIons];
    precision_t heat_neutrals;
    precision_t heat_ions;
    precision_t heat_electrons;
  };

#endif  // INCLUDE_CALC_CHEMISTRY_H_
