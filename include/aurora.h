// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_AURORA_H_
#define INCLUDE_AURORA_H_

/**********************************************************************
  * brief Read in a file containing information about splitting ionization
  *
  * param neutrals the class that contains all info about the neutrals
  * param ions the class that contains all info about the ions
 **/

void read_aurora(Neutrals &neutrals,
		 Ions &ions);

arma_vec calculate_fang(float eflux,  // in ergs/cm2/s
			float avee,   // in keV
			float Ebin,   // eV
			arma_vec rhoH,
			std::vector<float> Ci,
			float dE,     // eV
			arma_vec H,
			bool DoDebug);

/**********************************************************************
  * brief Read in a file containing information about splitting ionization
  *
  * param grid the grid class to use
  * param neutrals the class that contains all info about the neutrals
  * param ions the class that contains all info about the ions
 **/

void calc_aurora(Grid grid,
		 Neutrals &neutrals,
		 Ions &ions);

#endif  // INCLUDE_AURORA_H_
