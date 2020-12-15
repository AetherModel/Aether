// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef AETHER_INCLUDE_CHEMISTRY_H_
#define AETHER_INCLUDE_CHEMISTRY_H_

#include <vector>
#include <string>

#include "../include/earth.h"
#include "../include/neutrals.h"
#include "../include/ions.h"
#include "../include/inputs.h"
#include "../include/report.h"


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
    
    float energy;
    float rate;
    float branching_ratio;
  
  };

  std::vector<reaction_type> reactions;
  long nReactions;
  
  Chemistry(Neutrals neutrals,
	    Ions ions,
	    Inputs args,
	    Report &report);

    
  void calc_chemistry(Neutrals &neutrals,
		      Ions &ions,
		      Times time,
		      Grid grid,
		      Report &report);

  void calc_chemical_sources(float *neutral_density,
			     float *ion_density,
			     float Tn,
			     float Ti,
			     float Te,
			     Report &report);

 private:

  struct sources_and_losses_type {

    float neutral_sources[nSpecies];
    float neutral_losses[nSpecies];
    float ion_sources[nIons];
    float ion_losses[nIons];
    
    float heat_neutrals;
    float heat_ions;
    float heat_electrons;

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

#endif // AETHER_INCLUDE_CHEMISTRY_H_
