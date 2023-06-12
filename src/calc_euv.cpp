// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <string>
#include <vector>

#include "aether.h"

// -----------------------------------------------------------------------------
// Runs through the steps of calculating the EUV energy deposition:
//   - calculate chapman integrals
//   - calculate EUV spectrum
//   - calculate ionization and heating
// -----------------------------------------------------------------------------

int calc_euv(Planets planet,
             Grid grid,
             Times time,
             Euv &euv,
             Neutrals &neutrals,
             Ions &ions,
             Indices indices,
             Inputs input,
             Report &report) {

  int iErr = 0;

  if (time.check_time_gate(input.get_dt_euv())) {
    std::string function = "Euv::calc_euv";
    static int iFunction = -1;
    report.enter(function, iFunction);

    if (input.get_is_student())
      report.print(-1, "(2) What function is this " +
		   input.get_student_name() + "? Found it: def");
    
    // Chapman integrals for EUV energy deposition:
    neutrals.calc_chapman(grid, report);

    iErr = euv.euvac(time, indices, report);
    iErr = euv.scale_from_1au(planet, time, report);

    calc_ionization_heating(euv, neutrals, ions, report);

    report.exit(function);
  }

  return iErr;
}

// -----------------------------------------------------------------------------
// Calculate EUV driven ionization and heating rates
// -----------------------------------------------------------------------------

void calc_ionization_heating(Euv euv,
                             Neutrals &neutrals,
                             Ions &ions,
                             Report &report) {

  int64_t iAlt, iWave, iSpecies;
  int i_, iIon, iIonization;

  std::string function = "calc_ionization_heating";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // Zero out all source terms:

  neutrals.heating_euv_scgc.zeros();

  for (iSpecies = 0; iSpecies < neutrals.nSpecies; iSpecies++)
    neutrals.species[iSpecies].ionization_scgc.zeros();

  for (iSpecies = 0; iSpecies < ions.nSpecies; iSpecies++)
    ions.species[iSpecies].ionization_scgc.zeros();

  int64_t nAlts = neutrals.heating_euv_scgc.n_slices;

  arma_mat tau2d = neutrals.heating_euv_scgc.slice(0);
  arma_mat intensity2d = neutrals.heating_euv_scgc.slice(0);
  arma_mat ionization2d = neutrals.heating_euv_scgc.slice(0);

  for (iAlt = 2; iAlt < nAlts - 2; iAlt++) {
    for (iWave = 0; iWave < euv.nWavelengths; iWave++) {

      tau2d.zeros();

      for (iSpecies = 0; iSpecies < neutrals.nSpecies; iSpecies++) {
        if (neutrals.species[iSpecies].iEuvAbsId_ > -1) {
          i_ = neutrals.species[iSpecies].iEuvAbsId_;
          tau2d = tau2d +
                  euv.waveinfo[i_].values[iWave] *
                  neutrals.species[iSpecies].chapman_scgc.slice(iAlt);
        }
      }

      intensity2d = euv.wavelengths_intensity_top[iWave] * exp(-1.0 * tau2d);

      for (iSpecies = 0; iSpecies < neutrals.nSpecies; iSpecies++) {
        // Calculate Photo-Absorbtion for each species and add them up:
        // index of photo abs cross section
        i_ = neutrals.species[iSpecies].iEuvAbsId_;

        if (i_ > -1) {
          neutrals.heating_euv_scgc.slice(iAlt) =
            neutrals.heating_euv_scgc.slice(iAlt) +
            euv.wavelengths_energy[iWave] *
            euv.waveinfo[i_].values[iWave] *  // cross section
            (intensity2d %
             neutrals.species[iSpecies].density_scgc.slice(iAlt) );
        }

        for (iIonization = 0;
             iIonization < neutrals.species[iSpecies].nEuvIonSpecies;
             iIonization++) {

          i_ = neutrals.species[iSpecies].iEuvIonId_[iIonization];

          ionization2d =
            euv.waveinfo[i_].values[iWave] *  // cross section
            intensity2d %
            neutrals.species[iSpecies].density_scgc.slice(iAlt);

          neutrals.species[iSpecies].ionization_scgc.slice(iAlt) =
            neutrals.species[iSpecies].ionization_scgc(iAlt) + ionization2d;

          iIon = neutrals.species[iSpecies].iEuvIonSpecies_[iIonization];
          ions.species[iIon].ionization_scgc.slice(iAlt) =
            ions.species[iIon].ionization_scgc.slice(iAlt) + ionization2d;
        }  // iIonization
      }  // iSpecies
    }  // iWave
  }  // iAlt

  neutrals.heating_euv_scgc =
    neutrals.heating_efficiency *
    neutrals.heating_euv_scgc /
    neutrals.rho_scgc /
    neutrals.Cv_scgc;

  report.exit(function);
  return;
}
