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
      //set aurora ion and neutral index and coef
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

fvec calculate_maxwellian(float eflux,  // in ergs/cm2/s
                          float avee,   // in keV
                          fvec energies) {  // in eV

  // Change all units to be in eV and cm2:
  float E0 = avee / 2 * 1000.0; // characteristic energy in eV
  float Q0 = eflux * 6.242e11;  // eV/cm2/s (ergs -> eV)
  float a = Q0 / 2 / (E0 * E0 * E0);  // cm2/s/eV2

  fvec diff_num_flux = a * energies % exp(-energies / E0);  //eV/cm2/s
  return  diff_num_flux;
}

// -----------------------------------------------------------------------------
// Function to Calculate ionization rate for 1 Ebin for 1 alt profile
// -----------------------------------------------------------------------------
fvec calculate_fang(float eflux,  // in ergs/cm2/s
                    float avee,   // in keV
                    float Ebin,   // eV
                    fvec rhoH,
                    std::vector<float> Ci,
                    float dE,     // eV
                    fvec H,
                    Report &report) {

  // report
  std::string function = "calc_fang";
  static int iFunction = -1;
  report.enter(function, iFunction);

  float de = 0.035;  // keV move to before any loops (beginning of function)
  float char_e = avee / 2;  // keV
  long double Q0 = eflux * 6.242e11;   // eV/cm2/s
  float E0 = char_e * 1000;  // eV

  long double a = Q0 / 2 / (E0 * E0 * E0);  // cm2/s/eV2
  long double flux_max = a * Ebin * exp(-Ebin / E0);  //  /cm2/s/eV
  long double QE = flux_max * Ebin * dE; //  eV/cm2/s

  long double E_keV = Ebin / 1000.0;   // Convert to keV
  long double Q_keV = QE / 1000.0;   // Convert to keV

  fvec yE = (2.0 / E_keV) * pow( (rhoH / (6e-6)), 0.7);
  fvec fyE =
    (Ci[0] * pow(yE, Ci[1])) % exp((-Ci[2] * pow(yE, Ci[3]))) +
    (Ci[4] * pow(yE, Ci[5])) % exp((-Ci[6] * pow(yE, Ci[7])));
  fvec fac = Q_keV / de / H;
  fvec qtot = fyE % fac;

  report.exit(function);
  return qtot;
}

// -----------------------------------------------------------------------------
// Function to Calculate ionization rate for 1 Ebin for 1 alt profile
// -----------------------------------------------------------------------------

fvec calculate_fang_v2(float energy_bin,
                       float diff_energy_flux,
                       fvec rhoH,
                       std::vector<float> Ci,
                       fvec scale_height,
                       Report &report) {

  // report
  std::string function = "calc_fang";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // rhoH is in kg/m2, but need it in g/cm2 (/10)
  // scale_height is in m, but need it in cm
  // energy_bin is in eV but need it in keV
  // diff_energy_flux is in eV/cm2/s but need it in keV/cm2/s

  float de = 0.035;  // keV
  float E_mono = energy_bin / 1000.0;   // Convert to keV
  float Q_mono = diff_energy_flux / 1000.0;   // Convert to keV
  dvec H = conv_to<dvec>::from(scale_height) * 100.0; // Convert to cm

  // Eqn. 1 of Fang et al [2010]:
  // This *100 below should be a /10, but we don't understand why it isn't:
  dvec rhoHnorm = conv_to<dvec>::from(rhoH * 100.0 / 6e-6);
  dvec yE = (2.0 / E_mono) * pow( rhoHnorm, 0.7);
  // Eqn. 4 of Fang et al [2010]:
  dvec fyE =
    (Ci[0] * pow(yE, Ci[1])) % exp((-Ci[2] * pow(yE, Ci[3]))) +
    (Ci[4] * pow(yE, Ci[5])) % exp((-Ci[6] * pow(yE, Ci[7])));
  // Eqn. 3 of Fang et al [2010] (parenthesis):
  dvec fac = Q_mono / de / H;
  // Eqn. 3 of Fang et al [2010] (solve for Qtot(z)):
  fvec q_tot = conv_to<fvec>::from(fyE % fac);
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

  // some variables
  int64_t iSpecies;
  int64_t iAlt, iLon, iLat;
  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();

  // DENSITY INTEGRAL CALULATION ( done in calc_neutral_derived.cpp line
  // 170 rho_alt_int_scgc species[iSpecies].rho_alt_int_scgc =
  // integral3d * species[iSpecies].mass;

  // SET UP PIJ VALUES
  static mat Pij = { {1.25, 1.45903, -2.42e-1, 5.95e-2},
    {2.24, -4.23e-7, 1.36e-2, 2.53e-3},
    {1.42, 1.45e-1, 1.70e-2, 6.40e-4},
    {0.248775, -1.51e-1, 6.31e-9, 1.24e-3},
    {-0.465119, -1.05e-1, -8.96e-2, 1.22e-2},
    {3.86e-1, 1.75e-3, -7.43e-4, 4.61e-4},
    {-6.45e-1, 8.50e-4, -4.29e-2, -2.99e-3},
    {9.49e-1, 1.97e-1, -2.51e-1, -2.07e-3}
  };

  static std::vector<std::vector<float>> CiArray;
  static bool IsFirstTime = 1;

  // ENERGY BINS AND DE (E in eV)
  static float min = 100;
  static float max = 1000000;
  static float Emin = log(min);
  static float Emax = log(max);
  static int nBins = 101;
  static fvec auroral_energies(nBins);
  static fvec auroral_energy_widths(nBins);
  std::vector<float> Ci;

  if (IsFirstTime) {
    // read aurora csv
    read_aurora(neutrals, ions, args, report);

    float lnE;

    for (int64_t iBin = 0; iBin < nBins; iBin++) {
      float energy = exp(Emin + iBin * (Emax - Emin) / (nBins - 1));
      auroral_energies(iBin) = energy;
    }

    auroral_energy_widths = calc_bin_widths(auroral_energies);

    for (int64_t iBin = 0; iBin < nBins; iBin++) {

      lnE = log(auroral_energies(iBin) / 1000.0); // eV -> keV

      // loop through Pij values to get vector of Ci values
      for (int i = 0; i < 8; i++) {
        float tot = 0;

        for (int j = 0; j < 4; j++)
          tot = tot +  Pij.at(i, j) * pow(lnE, j);

        Ci.push_back(exp(tot));
      }

      CiArray.push_back(Ci);
    }

    IsFirstTime = 0;
  }

  fvec rhoH1d;
  fcube scale_height;
  fvec ionization1d;
  fvec H;
  fvec yE;
  fvec rho_tube;
  fvec weighted_sum;
  float coef;
  fvec neutral_density_tube;
  fvec ionization_tube, ionization_species;

  int iIon_;

  rhoH1d.set_size(nAlts);
  ionization1d.set_size(nAlts);
  weighted_sum.set_size(nAlts);

  scale_height = cKB * neutrals.temperature_scgc /
                 (neutrals.mean_major_mass_scgc % abs(grid.gravity_scgc));

  float eflux;
  float avee;
  fvec diff_num_flux;
  fvec diff_energy_flux;

  // loop through each altitude and calculate ionization
  for (iLon = 0; iLon < nLons ; iLon++) {
    for (iLat = 0; iLat < nLats ; iLat++) {

      eflux = ions.eflux(iLon, iLat);  // in ergs/cm2/s
      avee = ions.avee(iLon, iLat);  // in keV

      if (eflux > 0.01) {

        // Step 1: Calculate the height-integrated mass density:

        rhoH1d.zeros();

        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
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

        H = scale_height.tube(iLon, iLat);
        fvec temp;

        ionization1d.zeros();

        for (int iBin = 0; iBin < nBins; iBin++) {
          Ci = CiArray[iBin];
          temp = calculate_fang_v2(auroral_energies(iBin),
                                   diff_energy_flux(iBin),
                                   rhoH1d,
                                   Ci,
                                   H,
                                   report);
          ionization1d = ionization1d + temp;
        }

        // This /1e4 shouldn't be there, but I don't know why it doesn't work
        // without it.
        ionization1d = ionization1d * pcm3topm3 / 1.0e4;  // /cm3 -> /m3

        // Step 3: Distribute ionization among neutrals:
        // Need to figure out which species get what percentage of the
        // ionization, so we compute a weighted average given the
        // weights (coef or Aurora_Coef) and the neutral density
        weighted_sum.zeros();

        for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          if (neutrals.species[iSpecies].nAuroraIonSpecies > 0) {
            neutral_density_tube =
              neutrals.species[iSpecies].density_scgc.tube(iLon, iLat);
            coef = neutrals.species[iSpecies].Aurora_Coef;
            weighted_sum = weighted_sum + (coef * neutral_density_tube);
          }
        }

        // for each species of neutrals that gets aurora,
        for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

          if (neutrals.species[iSpecies].nAuroraIonSpecies != 0) {

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
      }  // eflux > 0.01
    }  // nLats
  }  // nLons

  report.exit(function);
}














