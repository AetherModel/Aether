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

  int iSpecies;

  std::string function = "Chemistry::calc_chemistry";
  static int iFunction = -1;
  report.enter(function, iFunction);  

  float dt = time.get_dt();
  
  // ------------------------------------
  // Calculate electron densities
  // ------------------------------------

  ions.fill_electrons(report);

  // ----------------------------------------------------------
  // Initialize the sources and losses with EUV stuff:
  // ----------------------------------------------------------

  for (iSpecies=0; iSpecies < nSpecies; iSpecies++) {
    neutrals.neutrals[iSpecies].losses_scgc = neutrals.neutrals[iSpecies].ionization_scgc;
    neutrals.neutrals[iSpecies].sources_scgc.zeros();
  }

  for (iSpecies=0; iSpecies < nIons; iSpecies++) {
    ions.species[iSpecies].losses_scgc.zeros();
    ions.species[iSpecies].sources_scgc = ions.species[iSpecies].ionization_scgc;
  }

  // ----------------------------------------------------
  // Calculate the chemical sources and losses
  // ----------------------------------------------------

  calc_chemical_sources(neutrals, ions, report);

  // ---------------------------------------------------------
  // Once sources and losses are done, solve for new densities
  // ---------------------------------------------------------

  for (iSpecies=0; iSpecies < nSpecies; iSpecies++) {
    neutrals.neutrals[iSpecies].density_scgc =
      solver_chemistry_new(neutrals.neutrals[iSpecies].density_scgc,
			   neutrals.neutrals[iSpecies].sources_scgc,
			   neutrals.neutrals[iSpecies].losses_scgc,
			   dt);
  }

  for (iSpecies=0; iSpecies < nIons; iSpecies++) {
    ions.species[iSpecies].density_scgc =
      solver_chemistry_new(ions.species[iSpecies].density_scgc,
			   ions.species[iSpecies].sources_scgc,
			   ions.species[iSpecies].losses_scgc,
			   dt);
  }

  // ---------------------------------------------------------
  // Recalculate electrons
  // ---------------------------------------------------------

  ions.fill_electrons(report);

  report.exit(function);
  return;
}
