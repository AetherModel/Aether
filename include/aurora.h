// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_AURORA_H_
#define INCLUDE_AURORA_H_

/**********************************************************************
  * brief Read in a file containing information about splitting ionization
  *
  * param neutrals the class that contains all info about the neutrals
  * param ions the class that contains all info about the ions
  * param input info about how user has configured things
  * param report allow reporting to occur
 **/

void read_aurora(Neutrals &neutrals,
		 Ions &ions,
		 Inputs args,
		 Report &report);

arma_vec calculate_fang(float eflux,  // in ergs/cm2/s
                    float avee,   // in keV
                    float Ebin,   // eV
                    arma_vec rhoH,
                    std::vector<float> Ci,
                    float dE,     // eV
                    arma_vec H,
                    Report &report);

/**********************************************************************
  * brief Read in a file containing information about splitting ionization
  *
  * param grid the grid class to use
  * param neutrals the class that contains all info about the neutrals
  * param ions the class that contains all info about the ions
  * param input info about how user has configured things
  * param report allow reporting to occur
 **/

void calc_aurora(Grid grid,
		 Neutrals &neutrals,
		 Ions &ions,
		 Inputs args,
		 Report &report);

#endif  // INCLUDE_AURORA_H_
