// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef AETHER_INCLUDE_CALC_CHEMISTRY_H_
#define AETHER_INCLUDE_CALC_CHEMISTRY_H_

class Chemistry {

 public:

  Chemistry(Inputs args, Report report);
  
  struct reaction_type {

    // Reactions:
    // loss1 + loss2 + loss3 -> source1 + source2 + source3
  
    std::vector<std::string> sources_names;
    std::vector<std::string> losses_names;
    
    std::vector<int> sources_ids;
    std::vector<int> sources_IsNeutrals;
  
    std::vector<int> losses_ids;
    std::vector<int> losses_IsNeutrals;

    float rate;
  
  };

  std::vector<reaction_type> reactions;
  
  struct sources_and_losses_type {

    float neutral_sources[nSpecies];
    float neutral_losses[nSpecies];
    float ion_sources[nIons];
    float ion_losses[nIons];
    
    float heat_neutrals;
    float heat_ions;
    float heat_electrons;

  };

  void calc_chemistry(Neutrals &neutrals,
		      Ions &ions,
		      Times time,
		      Grid grid,
		      Report report);

  void calc_chemical_sources(float neutral_density,
			     float ion_density,
			     float Tn,
			     float Ti,
			     float Te,
			     sources_and_losses_type &sources_and_losses;
			     Report report);
  

#endif // AETHER_INCLUDE_CALC_CHEMISTRY_H_
