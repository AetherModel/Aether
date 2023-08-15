// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <cmath>
#include "../include/aether.h"

// -----------------------------------------------------------------------------
// Calculates the ion-neutral collision frequencies given parameters
// read in from the CSV file.
// -----------------------------------------------------------------------------

void calc_ion_neutral_coll_freq(Neutrals &neutrals, Ions &ions) {

  std::string function = "calc_ion_neutral_coll_freq";
  static int iFunction = -1;
  report.enter(function, iFunction);
  arma_cube t, one_minus_log;

  for (int iIon = 0; iIon < ions.nSpecies; iIon++) {
    if (ions.species[iIon].nu_ion_neutral_coef.size() > 0) {
      for (int iNeutral = 0; iNeutral < neutrals.nSpecies; iNeutral++) {
        ions.species[iIon].nu_ion_neutral_vcgc[iNeutral].zeros();

        if (ions.species[iIon].nu_is_resonant[iNeutral]) {
          t = (ions.species[iIon].nu_in_res_tn_frac[iNeutral] *
               neutrals.temperature_scgc +
               ions.species[iIon].nu_in_res_tn_frac[iNeutral] *
               ions.temperature_scgc);
          one_minus_log =
            (1.0 - ions.species[iIon].nu_in_res_coef2[iNeutral] *
             log10(t));
          ions.species[iIon].nu_ion_neutral_vcgc[iNeutral] =
            ions.species[iIon].nu_in_res_coef1[iNeutral] *
            neutrals.species[iNeutral].density_scgc %
            sqrt(t) % one_minus_log % one_minus_log;
        } else {
          ions.species[iIon].nu_ion_neutral_vcgc[iNeutral] =
            ions.species[iIon].nu_ion_neutral_coef[iNeutral] *
            neutrals.species[iNeutral].density_scgc;
        } // is resonant
      } // for neutrals
    } // if size > 0
  } // for ions

  report.exit(function);
  return;
}
