// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "../include/report.h"
#include "../include/chemistry.h"

void Chemistry::calc_chemical_sources(float *neutral_density,
				      float *ion_density,
				      float Tn,
				      float Ti,
				      float Te,
				      Report &report) {

  std::string function = "Chemistry::calc_chemical_sources";
  static int iFunction = -1;
  report.enter(function, iFunction);  

  long iReaction, iLoss, iSource;
  float change, rate;
  int IsNeutral, id_;
  
  for (iReaction = 0; iReaction < nReactions; iReaction++) {

    // std::cout << "iReaction : " << iReaction << " " << nReactions << "\n";
    // display_reaction(reactions[iReaction]);
    
    // First calculate reaction rate:

    rate = reactions[iReaction].rate;

    // Second calculate the amount of change:

    change = rate;
    for (iLoss = 0; iLoss < reactions[iReaction].nLosses; iLoss++) {

      IsNeutral = reactions[iReaction].losses_IsNeutral[iLoss];
      id_ = reactions[iReaction].losses_ids[iLoss];
      
      if (IsNeutral) change = change * neutral_density[id_];
      else change = change * ion_density[id_];
      
    }

    // Third add change to the different consituents:

    for (iLoss = 0; iLoss < reactions[iReaction].nLosses; iLoss++) {

      IsNeutral = reactions[iReaction].losses_IsNeutral[iLoss];
      id_ = reactions[iReaction].losses_ids[iLoss];

      if (IsNeutral)
	sources_and_losses.neutral_losses[id_] = 
	  sources_and_losses.neutral_losses[id_] + change;
      else
	sources_and_losses.ion_losses[id_] = 
	  sources_and_losses.ion_losses[id_] + change;

    }

    for (iSource = 0; iSource < reactions[iReaction].nSources; iSource++) {

      IsNeutral = reactions[iReaction].sources_IsNeutral[iSource];
      id_ = reactions[iReaction].sources_ids[iSource];
      
      if (IsNeutral)
	sources_and_losses.neutral_sources[id_] = 
	  sources_and_losses.neutral_sources[id_] + change;
      else
	sources_and_losses.ion_sources[id_] = 
	  sources_and_losses.ion_sources[id_] + change;

    }
    
  }
  
  report.exit(function);

  return;

}

