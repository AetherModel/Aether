// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <cmath>
#include "../include/aether.h"

// -----------------------------------------------------------------------------
// Calculates the ion-neutral collision frequencies given parameters
// read in from the CSV file.
// -----------------------------------------------------------------------------

// removed this from the ions class // MB
void calc_ion_neutral_coll_freq(Neutrals &neutrals, Ions &ions, Report &report) {

  std::string function = "Ions::calc_ion_neutral_coll_freq";
  static int iFunction = -1;
  report.enter(function, iFunction);
  arma_cube t, one_minus_log;

  for (int iIon = 0; iIon < nIons; iIon++) {
    if (ions.species[iIon].nu_ion_neutral_coef.size() > 0) {
      for (int iNeutral = 0; iNeutral < nSpecies; iNeutral++) {
        species[iIon].nu_ion_neutral_vcgc[iNeutral].zeros();

        if (species[iIon].nu_is_resonant[iNeutral]) {
          t = (species[iIon].nu_in_res_tn_frac[iNeutral] *
               neutrals.temperature_scgc +
               species[iIon].nu_in_res_tn_frac[iNeutral] *
               ion_temperature_scgc);
          one_minus_log =
            (1.0 - species[iIon].nu_in_res_coef2[iNeutral] *
             log10(t));
          species[iIon].nu_ion_neutral_vcgc[iNeutral] =
            species[iIon].nu_in_res_coef1[iNeutral] *
            neutrals.species[iNeutral].density_scgc %
            sqrt(t) % one_minus_log % one_minus_log;
        } else {
          species[iIon].nu_ion_neutral_vcgc[iNeutral] =
            species[iIon].nu_ion_neutral_coef[iNeutral] *
            neutrals.species[iNeutral].density_scgc;
        }
      }
    }
  }

  report.exit(function);
  return;
}
