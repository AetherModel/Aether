// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "../include/sizes.h"

#include "../include/chemistry.h"
#include "../include/neutrals.h"
#include "../include/ions.h"
#include "../include/times.h"
#include "../include/grid.h"
#include "../include/report.h"
#include "../include/solvers.h"

void Chemistry::calc_chemistry(Neutrals &neutrals,
			       Ions &ions,
			       Times time,
			       Grid grid,
			       Report &report) {

  long nLons, nLats, nAlts, iLon, iLat, iAlt, index;
  int iSpecies;

  std::string function = "Chemistry::calc_chemistry";
  static int iFunction = -1;
  report.enter(function, iFunction);  

  nLons = grid.get_nLons();
  nLats = grid.get_nLats();
  nAlts = grid.get_nAlts();

  float neutral_density[nSpecies];
  float ion_density[nIons+1];  // add one for electron density
  float neutral_sources[nSpecies];
  float neutral_losses[nSpecies];
  float ion_sources[nIons];
  float ion_losses[nIons];

  float Tn, Te, Ti;

  float dt = time.get_dt();
  float old_density, source, loss, new_density;
  float heat_neutrals, heat_ions, heat_electrons;
  
  // ------------------------------------
  // Calculate electron densities
  // ------------------------------------

  ions.fill_electrons(report);

  for (iSpecies=0; iSpecies < nSpecies; iSpecies++) {
    neutrals.neutrals[iSpecies].losses_scgc = neutrals.neutrals[iSpecies].ionization_scgc;
    neutrals.neutrals[iSpecies].sources_scgc.zeros();
  }

  for (iSpecies=0; iSpecies < nIons; iSpecies++) {
    ions.species[iSpecies].losses_scgc.zeros();
    ions.species[iSpecies].sources_scgc = ions.species[iSpecies].ionization_scgc;
  }

  calc_chemical_sources(neutrals, ions, report);

  fcube norm_loss = neutrals.neutrals[0].losses_scgc;

  for (iSpecies=0; iSpecies < nSpecies; iSpecies++) {
    norm_loss = neutrals.neutrals[0].losses_scgc /
      (neutrals.neutrals[iSpecies].density_scgc + 1.0e-6);
    neutrals.neutrals[iSpecies].density_scgc =
      ( neutrals.neutrals[iSpecies].density_scgc + 
	dt * neutrals.neutrals[iSpecies].sources_scgc ) /
      (1.0 + dt * norm_loss);
  }

  for (iSpecies=0; iSpecies < nIons; iSpecies++) {
    norm_loss = ions.species[iSpecies].losses_scgc /
      (ions.species[iSpecies].density_scgc + 1.0e-6);
    ions.species[iSpecies].density_scgc =
      ( ions.species[iSpecies].density_scgc + 
    	dt * ions.species[iSpecies].sources_scgc ) /
      (1.0 + dt * norm_loss);
  }

  ions.fill_electrons(report);

  for (iLon = 0; iLon < nLons; iLon++) {
    for (iLat = 0; iLat < nLats; iLat++) {
      for (iAlt = 0; iAlt < nAlts; iAlt++) {
	
	index = ijk_geo_s3gc(iLon,iLat,iAlt);

	for (iSpecies=0; iSpecies < nSpecies; iSpecies++) {
	  neutrals.neutrals[iSpecies].density_s3gc[index] =
	    neutrals.neutrals[iSpecies].density_scgc(iLon,iLat,iAlt);
	}

	for (iSpecies=0; iSpecies <= nIons; iSpecies++) {
	  ions.species[iSpecies].density_s3gc[index] =
	    ions.species[iSpecies].density_scgc(iLon,iLat,iAlt);

	}
	ions.density_s3gc[index] = ions.density_scgc(iLon,iLat,iAlt);
	
      }
    }
  }

//  for (iAlt = 0; iAlt < nAlts; iAlt++) 
//    std::cout << iAlt << " " << ions.density_scgc(9,18,iAlt) << " " 
//	      << ions.species[0].ionization_scgc(9,18,iAlt) << " " 
//	      << "\n";
  
  report.exit(function);
  return;
}
