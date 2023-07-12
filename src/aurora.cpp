// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iostream>
#include <cmath>
#include <iostream>

#include "../include/aether.h"

// -----------------------------------------------------------------------------
// Read in aurora information
//  - because we need both neutrals and ions, we can't do this while reading
//    in the planet file.
// -----------------------------------------------------------------------------

void read_aurora(Neutrals &neutrals,
                 Ions &ions,
                 Inputs args,
                 Report &report) {

  std::ifstream myFile;
  myFile.open(args.get_aurora_file());

  if (myFile.good()) {
    std::vector<std::vector<std::string>> csv = read_csv(myFile);

    int nLines = csv.size();

    for (int iLine = 0; iLine < nLines; iLine++) {
      //Set the auroral ion and neutral indices and coefficients
      int iNeutral_ = neutrals.get_species_id(csv[iLine][0], report);
      int iIon_ = ions.get_species_id(csv[iLine][1], report);
      neutrals.species[iNeutral_].iAuroraIonSpecies_.push_back(iIon_);
      neutrals.species[iNeutral_].nAuroraIonSpecies++;
      neutrals.species[iNeutral_].Aurora_Coef = stod(csv[iLine][2]);
    }

    myFile.close();
  }

}

// -----------------------------------------------------------------------------
// Calculate the Maxellian distribution (eqn 6 in Fang et al. [2010])
// -----------------------------------------------------------------------------

arma_vec calculate_maxwellian(precision_t eflux,  // in ergs/cm2/s
                              precision_t avee,   // in keV
                              arma_vec energies) {  // in keV

  // Change all units to be in eV and cm2:
  precision_t E0 = avee / 2; // characteristic energy in keV
  precision_t Q0 = eflux * 6.242e11 / 1000.0;  //  keV/cm2/s
  precision_t a = Q0 / 2 / (E0 * E0 * E0);  // cm2/s/keV2

  arma_vec diff_num_flux = a * energies % exp(-energies / E0);  //  keV/cm2/s
  return  diff_num_flux;
}

// -----------------------------------------------------------------------------
// Function to Calculate ionization rate for 1 Ebin for 1 alt profile
// -----------------------------------------------------------------------------

arma_vec calculate_fang_v2(precision_t energy_bin,
                           precision_t diff_energy_flux,
                           arma_vec rhoH,
                           std::vector<precision_t> Ci,
                           arma_vec scale_height,
                           bool DoDebug,
                           Report &report) {

  // Set up function reporting
  std::string function = "calc_fang";
  static int iFunction = -1;

  if (DoDebug)
    report.enter(function, iFunction);

  // rhoH is in kg/m2, but need it in g/cm2 (/10)
  // scale_height needs to be in cm
  // energy_bin needs to be in keV
  // diff_energy_flux needs to be in keV/cm2/s

  precision_t de = 0.035;  // keV
  precision_t E_mono = energy_bin;
  precision_t Q_mono = diff_energy_flux;
  dvec H = conv_to<dvec>::from(scale_height);

  // Eqn. 1 of Fang et al [2010]:
  dvec rhoHnorm = conv_to<dvec>::from(rhoH / 10.0 / 6e-6);
  dvec yE = (2.0 / E_mono) * pow( rhoHnorm, 0.7);

  // Eqn. 4 of Fang et al [2010]:
  dvec fyE =
    (Ci[0] * pow(yE, Ci[1])) % exp((-Ci[2] * pow(yE, Ci[3]))) +
    (Ci[4] * pow(yE, Ci[5])) % exp((-Ci[6] * pow(yE, Ci[7])));

  // Eqn. 3 of Fang et al [2010] (parenthesis):
  dvec fac = Q_mono / (de * H);

  // Eqn. 3 of Fang et al [2010] (solve for Qtot(z), ionization rate):
  arma_vec q_tot = conv_to<arma_vec>::from(fyE % fac);

  if (DoDebug)
    report.exit(function);

  return q_tot;
}

// -----------------------------------------------------------------------------
// Calculate aurora
// -----------------------------------------------------------------------------
void calc_aurora(Grid grid,
                 Neutrals &neutrals,
                 Ions &ions,
                 Inputs args,
                 Report &report) {

  std::string function = "calc_aurora";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // Set the grid, location, and species variables
  int64_t iSpecies;
  int64_t iAlt, iLon, iLat;
  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();

  // DENSITY INTEGRAL CALULATION ( done in calc_neutral_derived.cpp line
  // 170 rho_alt_int_scgc species[iSpecies].rho_alt_int_scgc =
  // integral3d * species[iSpecies].mass;

  // SET UP PIJ VALUES - these are directly from Fang et al. [2010]:
  static mat Pij = {
    {1.25, 1.45903, -2.42e-1, 5.95e-2},
    {2.24, -4.23e-7, 1.36e-2, 2.53e-3},
    {1.42, 1.45e-1, 1.70e-2, 6.40e-4},
    {0.248775, -1.51e-1, 6.31e-9, 1.24e-3},
    {-0.465119, -1.05e-1, -8.96e-2, 1.22e-2},
    {3.86e-1, 1.75e-3, -7.43e-4, 4.61e-4},
    {-6.45e-1, 8.50e-4, -4.29e-2, -2.99e-3},
    {9.49e-1, 1.97e-1, -2.51e-3, -2.07e-3}
  };

  static std::vector<std::vector<precision_t>> CiArray;
  static bool IsFirstTime = 1;

  // ENERGY BINS AND DE (E in eV)
  static precision_t min = 100;
  static precision_t max = 1000000;
  static precision_t Emin = log(min);
  static precision_t Emax = log(max);
  static int nBins = 101;
  static arma_vec auroral_energies(nBins);
  static arma_vec auroral_energy_widths(nBins);
  std::vector<precision_t> Ci;

  if (IsFirstTime) {
    // Initialize the aurora using the auroral csv file
    read_aurora(neutrals, ions, args, report);

    precision_t lnE;

    for (int64_t iBin = 0; iBin < nBins; iBin++) {
      precision_t energy = exp(Emin + iBin * (Emax - Emin) / (nBins - 1));
      // convert from eV -> keV
      auroral_energies(iBin) = energy / 1000.0;
    }

    auroral_energy_widths = calc_bin_widths(auroral_energies);

    for (int64_t iBin = 0; iBin < nBins; iBin++) {

      lnE = log(auroral_energies(iBin));

      // loop through Pij values to get vector of Ci values.  This is
      // directly from Fang et al., [2010]:
      for (int i = 0; i < 8; i++) {
        precision_t tot = 0;

        for (int j = 0; j < 4; j++)
          tot = tot +  Pij.at(i, j) * pow(lnE, j);

        Ci.push_back(exp(tot));
      }

      CiArray.push_back(Ci);
    }

    IsFirstTime = 0;
  }

  arma_vec rhoH1d;
  arma_cube scale_height;
  arma_vec ionization1d;
  arma_vec H;
  arma_vec yE;
  arma_vec rho_tube;
  arma_vec weighted_sum;
  precision_t coef;
  arma_vec neutral_density_tube;
  arma_vec ionization_tube, ionization_species;

  int iIon_;

  rhoH1d.set_size(nAlts);
  ionization1d.set_size(nAlts);
  weighted_sum.set_size(nAlts);

  scale_height = cKB * neutrals.temperature_scgc /
                 (neutrals.mean_major_mass_scgc % abs(grid.gravity_vcgc[2]));

  precision_t eflux;
  precision_t avee;
  arma_vec diff_num_flux;
  arma_vec diff_energy_flux;
  bool DoDebug = false;

  // loop through each altitude and calculate ionization
  for (iLon = 0; iLon < nLons ; iLon++) {
    for (iLat = 0; iLat < nLats ; iLat++) {

      eflux = ions.eflux(iLon, iLat);  // in ergs/cm2/s
      avee = ions.avee(iLon, iLat);  // in keV

      if (eflux > 0.1) {

        // Step 1: Calculate the height-integrated mass density:
        rhoH1d.zeros();

        for (iSpecies = 0; iSpecies < neutrals.nSpecies; iSpecies++) {
          rho_tube =
            neutrals.species[iSpecies].rho_alt_int_scgc.tube(iLon, iLat);
          rhoH1d = rhoH1d + rho_tube;
        }

        // Step 2: Calculate the distribution function:
        diff_num_flux = calculate_maxwellian(eflux,
                                             avee,
                                             auroral_energies);

        // Step 3: Calculate the differential energy flux:
        diff_energy_flux =
          diff_num_flux % auroral_energies % auroral_energy_widths;

        // Step 4: Calculate ionization rates from Fang (all energy bins):
        // in cm (from meters)
        H = scale_height.tube(iLon, iLat) * 100.0;
        arma_vec temp;

        ionization1d.zeros();

        for (int iBin = 0; iBin < nBins; iBin++) {
          Ci = CiArray[iBin];
          temp = calculate_fang_v2(auroral_energies(iBin),
                                   diff_energy_flux(iBin),
                                   rhoH1d,
                                   Ci,
                                   H,
                                   DoDebug,
                                   report);
          ionization1d = ionization1d + temp;
        }

        // /cm3 -> /m3
        ionization1d = ionization1d * pcm3topm3;

        // Step 5: Distribute ionization among neutrals:
        // Need to figure out which species get what percentage of the
        // ionization, so we compute a weighted average given the
        // weights (coef or Aurora_Coef) and the neutral density
        weighted_sum.zeros();

        for (int iSpecies = 0; iSpecies < neutrals.nSpecies; iSpecies++) {
          if (neutrals.species[iSpecies].nAuroraIonSpecies > 0) {
            neutral_density_tube =
              neutrals.species[iSpecies].density_scgc.tube(iLon, iLat);
            coef = neutrals.species[iSpecies].Aurora_Coef;
            weighted_sum = weighted_sum + (coef * neutral_density_tube);
          }
        }

        // Cycle through each species of neutrals that gets aurora,
        for (int iSpecies = 0; iSpecies < neutrals.nSpecies; iSpecies++) {

          if (neutrals.species[iSpecies].nAuroraIonSpecies > 0) {

            // Parse the ionization into the species-specific parts:
            neutral_density_tube =
              neutrals.species[iSpecies].density_scgc.tube(iLon, iLat);
            coef = neutrals.species[iSpecies].Aurora_Coef;
            ionization_species = (coef * ionization1d %
                                  neutral_density_tube / weighted_sum);

            // Add to neutrals:
            ionization_tube =
              neutrals.species[iSpecies].ionization_scgc.tube(iLon, iLat);
            neutrals.species[iSpecies].ionization_scgc.tube(iLon, iLat) =
              ionization_tube + ionization_species;

            // Add to ions:
            for (int iIon = 0;
                 iIon < neutrals.species[iSpecies].nAuroraIonSpecies;
                 iIon++) {
              iIon_ = neutrals.species[iSpecies].iAuroraIonSpecies_[iIon];
              ionization_tube =
                ions.species[iIon_].ionization_scgc.tube(iLon, iLat);
              ions.species[iIon_].ionization_scgc.tube(iLon, iLat) =
                ionization_tube + ionization_species;

            }  // nAuroraIonSpecies
          }  // if nAuroraIonSpecies > 0
        }  // nSpecies
      }  // eflux > 0.1
    }  // nLats
  }  // nLons

  report.exit(function);
}














