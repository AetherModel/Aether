// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"

/// @brief Calculates photoelectron heating efficiency and heating rate
/// @details Based on two key references:
///   1. Swartz & Nisbet (1972)
///   2. Smithro & Solomon (2008)
///   
/// Uses equations 9-12 from (Zhu & Ridley, 2016)
///   https://doi.org/10.1016/j.jastp.2016.01.005
/// 
/// @param neutrals 
/// @param ions 
/// @return Qphe 
arma_cube calc_photoelectron_heating(Ions &ions, arma_cube epsilon);


// Since this is used a few times, calculate it separately & pass it to the functions.
// From (Smithro and Solomon, 2008)
arma_cube calc_epsilon(Neutrals &neutrals, Ions &ions) {

  // from Swartz & Nisbet (1972): we need neutral species id's for o2, n2 & o
  int64_t inO2 = neutrals.get_species_id("O2");
  int64_t inN2 = neutrals.get_species_id("N2");
  int64_t inO = neutrals.get_species_id("O");
  
  arma_cube epsilon, x, logx;
  
  x = ions.density_scgc / (neutrals.species[inO2].density_scgc 
                           + neutrals.species[inN2].density_scgc 
                           + neutrals.species[inO].density_scgc);
  // should rarely need to be used:
  x.clamp(2e-8, 10.0);

  logx = log(x);

  epsilon = exp(5.342 + 1.056*logx - pow(4.392e-2*logx, 2) 
                - pow(5.9e-2*logx,3) - 9.346e-3*pow(logx,4)
                - 5.755e-4*pow(logx,5) - 1.249e-5*pow(logx,6)
                )*1.6e-19;
  return exp(5.342 + 1.056*logx - pow(4.392e-2*logx, 2) 
                - pow(5.9e-2*logx,3) - 9.346e-3*pow(logx,4)
                - 5.755e-4*pow(logx,5) - 1.249e-5*pow(logx,6)
                )*1.6e-19;}

// --------------------------------------------------------------------------
// Heating terms:
//  - [x] photoelectrons
//  - [ ] auroral ionization
//  - [ ] e- ion collisions
//  - [ ] e- neutral collisions (elastic & inelastic)
//  - [ ] e- chemistry (O2, V2 vibration; O fine structure, O exitation)
// --------------------------------------------------------------------------




// --------------------------------------------------------------------------
// TODO (#24): this currently just sets the electron temperature to the neutral temperature
// --------------------------------------------------------------------------

void Ions::calc_electron_temperature(Neutrals neutrals, Grid grid) {

  std::string function = "Ions::calc_electron_temperature";
  static int iFunction = -1;
  report.enter(function, iFunction);

  arma_cube epsilon = calc_epsilon(neutrals, *this);

  arma_cube Qphe = calc_photoelectron_heating(*this, epsilon);

  electron_temperature_scgc = neutrals.temperature_scgc;

  report.exit(function);
}


// --------------------------------------------------------------------------
// Calculate photoelectron heating
// --------------------------------------------------------------------------
arma_cube calc_photoelectron_heating(Ions &ions,
                                     arma_cube epsilon) {

  int64_t nIons = ions.nSpecies;

  // Initialize Qphe & IonsIonizationRate (sum of ionization rates)
  // to the same size as epsilon and then zero them out:
  // (Qphe is the product of epsilon & IonsIonizationRate, so is not zeroed.)
  arma_cube Qphe = epsilon, IonsIonizationRate = epsilon;
  IonsIonizationRate.zeros();

  for (int64_t iIon = 0; iIon < nIons; iIon++) {
    IonsIonizationRate += ions.species[iIon].ionization_scgc;
  }

  Qphe = epsilon % IonsIonizationRate;

  return Qphe;
}
