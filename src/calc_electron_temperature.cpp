// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"

/// @brief Calculate epsilon
/// @details intermediate variable used in photoelectron & ionization heating
/// From (Smithro & Solomon, 2008).
/// @param neutrals 
/// @param ions 
/// @return epsilon 
arma_cube calc_epsilon(Neutrals &neutrals, Ions &ions);


/// @brief Calculates photoelectron heating
/// @details Based on (Swartz & Nisbet, 1972) & (Smithro & Solomon, 2008)
///   
/// Uses equations 9-12 from (Zhu & Ridley, 2016)
///   https://doi.org/10.1016/j.jastp.2016.01.005
/// 
/// @param ions 
/// @param epsilon 
/// @return Qphe 
arma_cube calc_photoelectron_heating(Ions &ions, arma_cube epsilon);


/// @brief Calculates auroral heating
/// @details NOTE: in GITM this is solved separately for ion precipitation & auroral 
/// ionization. In Aether these are both in ions.species[iIon].ionization_scgc...
/// @param ions 
/// @param epsilon 
/// @return Qaurora 
arma_cube calc_ionization_heating(Ions &ions, arma_cube epsilon);


/// @brief Calculates electron-ion (elastic) collisional heating
/// @details From Schunk and Nagy 2009, and Bei-Chen Zhang and Y. Kamide 2003
/// - This differs slightly from the GITM implementation, which assumes several ion species are present.
///   Instead, here we use each ion species for the sum.
/// - electon-ion collision frequency (from Schunk and Nagy 2009) = 5.45E-5
/// - This is capable of handling BOTH the bulk & individual ion temperatures
/// @param ions 
/// @return vector<Qeicp, Qeicm, Qeic_v>
std::vector<arma_cube> calc_electron_ion_collisions(Ions &ions);


/// @brief Calculates electron-neutral elastic collisional heating
/// @details From Schunk and Nagy 2009
/// @param ions
/// @param neutrals
/// @return vector<Qencp, Qencm, Qenc_v>
std::vector<arma_cube> calc_electron_neutral_collisions(Ions &ions, Neutrals &neutrals);


// --------------------------------------------------------------------------
// Heating terms:
//  - [x] photoelectrons
//  - [x] auroral ionization (from ion precipitation & auroral ionization)
//  - [x] e- ion collisions
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

  // Photoelectron, (all) ionization heating terms:
  arma_cube epsilon, Qphe, QIonization;

  // electron-ion collisions return a vector of cubes, one for Qe, one for Qi, one for friction:
  std::vector<arma_cube> Qeic;
  arma_cube Qeicm, Qeicp, Qeic_v;
  
  // (elastic) Electron-neutral collisions:
  std::vector<arma_cube> Qenc;
  arma_cube Qencm, Qencp, Qenc_v;

  // Initialize everything to zero!

  epsilon.set_size(grid.get_nLons(), grid.get_nLats(), grid.get_nAlts());
  epsilon.zeros();
  Qphe = epsilon;
  QIonization = epsilon;
  Qeicm = epsilon;
  Qeicp = epsilon;
  Qeic_v = epsilon;
  Qencm = epsilon;
  Qencp = epsilon;
  Qenc_v = epsilon;
  
  report.print(4, "Calculating epsilon");
  

  // Needed for both ionization & photoelectron heating:
  if (input.get_do_ionization_heating() || input.get_do_photoelectron_heating()) {
    epsilon = calc_epsilon(neutrals, *this);
  }

  report.print(4, "Calculating photoelectron heating");

  // Photoelectron heating
  if (input.get_do_photoelectron_heating()) {
    Qphe = calc_photoelectron_heating(*this, epsilon);
  }

  report.print(4, "Calculating ionization heating");

  // Ionization heating (includes all ionization sources)
  if (input.get_do_ionization_heating()) {
    QIonization = calc_ionization_heating(*this, epsilon);
  }

  report.print(4, "Calculating electron-ion collisions");

  // electron-ion collisions
  if (input.get_do_electron_ion_collisional_heating()) {
    Qeic = calc_electron_ion_collisions(*this);
    Qeicp = Qeic[0]; 
    Qeicm = Qeic[1]; 
    Qeic_v = Qeic[2]; // Friction
  }

  report.print(4, "Calculating electron-neutral collisions");

  // electron-neutral collisions
  if (input.get_do_electron_neutral_collisional_heating()) {
    Qenc = calc_electron_neutral_collisions(*this, neutrals);
    Qencp = Qenc[0]; 
    Qencm = Qenc[1]; 
    Qenc_v = Qenc[2]; // Friction
  }


  electron_temperature_scgc = neutrals.temperature_scgc;

  report.exit(function);
}



// Since this is used a few times, calculate it separately & pass it to the functions.
// From (Smithro and Solomon, 2008)
arma_cube calc_epsilon(Neutrals &neutrals, Ions &ions) {

  std::string function = "calc_epsilon";
  static int iFunction = -1;
  report.enter(function, iFunction);

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

  report.exit(function);
  return epsilon;
}



// --------------------------------------------------------------------------
// Calculate photoelectron heating
// --------------------------------------------------------------------------
arma_cube calc_photoelectron_heating(Ions &ions,
                                     arma_cube epsilon) {
  
  std::string function = "calc_photoelectron_heating";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t nIons = ions.nSpecies;

  // Initialize Qphe & IonsIonizationRate (sum of ionization rates)
  // to the same size as epsilon and then zero them out:
  // (Qphe is the product of epsilon & IonsIonizationRate, so is not zeroed.)
  arma_cube Qphe, IonsIonizationRate;
  IonsIonizationRate.set_size(epsilon.n_rows, epsilon.n_cols, epsilon.n_slices);
  IonsIonizationRate.zeros();

  for (int64_t iIon = 0; iIon < nIons; iIon++) {
    IonsIonizationRate += ions.species[iIon].ionization_scgc;
  }

  Qphe = epsilon % IonsIonizationRate;

  report.exit(function);
  return Qphe;
}


// --------------------------------------------------------------------------
// Calculate ionization heating
// --------------------------------------------------------------------------
arma_cube calc_ionization_heating(Ions &ions, arma_cube epsilon){

  std::string function = "calc_ionization_heating";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t nIons = ions.nSpecies;

  // auroral heating efficiency coefficient
  precision_t auroheat = 1.0;

  int64_t iO_3P = ions.get_species_id("O+");
  int64_t iO2P = ions.get_species_id("O2+");
  int64_t iN2P = ions.get_species_id("N2+");

  arma_cube QIonization = auroheat * epsilon % (ions.species[iO_3P].ionization_scgc 
                                                + ions.species[iO2P].ionization_scgc 
                                                + ions.species[iN2P].ionization_scgc);


  report.exit(function);
  return QIonization;
}

// --------------------------------------------------------------------------
// Calculate electron-ion collisions
// --------------------------------------------------------------------------
std::vector<arma_cube> calc_electron_ion_collisions(Ions &ions){

  std::string function = "calc_electron_ion_collisions";
  static int iFunction = -1;
  report.enter(function, iFunction);

  arma_cube Qeicp;
  Qeicp.set_size(ions.density_scgc.n_rows, ions.density_scgc.n_cols, ions.density_scgc.n_slices);
  Qeicp.zeros();
  arma_cube Qeicm = Qeicp;
  // Friction things
  arma_cube Qeic_v = Qeicp, dv2_ei = Qeicp;

  int64_t nSpecies = ions.nSpecies;

  if (input.get_do_calc_bulk_ion_temp()){
    // This is used when we calculate bulk ion temperature!
    report.print(3, "Using bulk ion temperature for electron-ion collisions");
    // Use all species, not just major species (different from GITM)
    for (int64_t iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      Qeicp += ions.species[iSpecies].density_scgc 
                / (cME + ions.species[iSpecies].mass);
    }

    Qeicp = Qeicp % ions.density_scgc * cME * 3.0 * cKB 
            % (ions.temperature_scgc - ions.electron_temperature_scgc)
            * 5.45e-5 / pow(ions.electron_temperature_scgc, 1.5);
  }
  else{
    // Individual ion temperatures:
    report.print(3, "Using individual ion temperatures for electron-ion collisions");
    // Use all species, not just major species (different from GITM)
    for (int64_t iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      Qeicp += ions.species[iSpecies].density_scgc 
               % (ions.species[iSpecies].temperature_scgc - ions.electron_temperature_scgc)
               / (cME + ions.species[iSpecies].mass);
    }

    Qeicp = Qeicp % ions.density_scgc * cME * 3.0 * cKB 
            * 5.45e-5 / pow(ions.electron_temperature_scgc, 1.5);
  }

  report.print(3, "Calculating frictional heating");

  // delta velocity **2 btwn e- & ions:
  // This uses the bulk ion velocity, not the individual ion velocity.
  // (Different from GITM): Uses all species' densities (so just ne), not just o+, o2+, n2+, no+, n+
  for (int64_t iDir = 0; iDir < 3; iDir++) {
    dv2_ei += pow(ions.velocity_vcgc[iDir] - ions.exb_vcgc[iDir], 2);
  }
  Qeic_v = ions.density_scgc * cME % dv2_ei * 5.45e-5 / pow(ions.electron_temperature_scgc, 1.5)
    % (ions.density_scgc);
    
  std::vector<arma_cube> Qeic = {Qeicp,
                                 Qeicp % ions.electron_temperature_scgc,
                                 Qeic_v};

  report.exit(function);
  
  return Qeic;
}

// --------------------------------------------------------------------------
// Calculate electron-neutral collisions
// --------------------------------------------------------------------------
std::vector<arma_cube> calc_electron_neutral_collisions(Ions &ions, Neutrals &neutrals){

  std::string function = "calc_electron_neutral_collisions";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // initialize & zero the quantities we need:
  arma_cube Qenc;
  Qenc.set_size(ions.density_scgc.n_rows, ions.density_scgc.n_cols, ions.density_scgc.n_slices);
  Qenc.zeros();
  arma_cube Qencp = Qenc;
  arma_cube Qencm = Qenc;
  // frictional things:
  arma_cube dv2_en = Qenc;
  arma_cube Qenc_v = Qenc;

  int64_t inO2 = neutrals.get_species_id("O2");
  int64_t inN2 = neutrals.get_species_id("N2");
  int64_t inO = neutrals.get_species_id("O");

  Qenc = ions.density_scgc * cME * 3.0 * cKB % (neutrals.temperature_scgc - ions.electron_temperature_scgc)
          % ((2.33e-11 * neutrals.species[inN2].density_scgc * 1.e-6 
              % (1 - 1.21e-4 * ions.electron_temperature_scgc) 
                % ions.electron_temperature_scgc / (cME + neutrals.species[inN2].mass))
            + (1.82e-10*neutrals.species[inO2].density_scgc*1.e-6
              % (1 + 3.60e-2 * pow(ions.electron_temperature_scgc, 0.5)) 
                % pow(ions.electron_temperature_scgc, 0.5)/(cME + neutrals.species[inO2].mass))
            + (8.90e-11*neutrals.species[inO].density_scgc*1.e-6
              % (1 + 5.70e-4 * ions.electron_temperature_scgc) 
                % pow(ions.electron_temperature_scgc, 0.5)/(cME + neutrals.species[inO].mass)) 
            );

  Qencp = ions.density_scgc * cME * 3.0 * cKB 
          % ((2.33e-11*neutrals.species[inN2].density_scgc*1.e-6
              % (1 - 1.21e-4*ions.electron_temperature_scgc) 
                % ions.electron_temperature_scgc / (cME + neutrals.species[inN2].mass))
            + (1.82e-10*neutrals.species[inO2].density_scgc*1.e-6
              % (1 + 3.60e-2 * pow(ions.electron_temperature_scgc, 0.5)) 
                % pow(ions.electron_temperature_scgc, 0.5)/(cME + neutrals.species[inO2].mass))
            + (8.90e-11*neutrals.species[inO].density_scgc*1.e-6
              % (1 + 5.70e-4*ions.electron_temperature_scgc) 
                % pow(ions.electron_temperature_scgc, 0.5)/(cME + neutrals.species[inO].mass))
          );

  // delta velocity **2 btwn e- & neutrals:
  for (int64_t iDir = 0; iDir < 3; iDir++) {
      dv2_en += pow(neutrals.velocity_vcgc[iDir] - ions.exb_vcgc[iDir], 2);
  }
  

  Qenc_v = ions.density_scgc * cME % dv2_en 
            %(2.33e-11 * neutrals.species[inN2].density_scgc * 1.e-6
                % (1 - 1.21e-4 * ions.electron_temperature_scgc) % ions.electron_temperature_scgc * neutrals.species[inN2].mass
                /(cME + neutrals.species[inN2].mass) 
              + 1.82e-10*neutrals.species[inO2].density_scgc*1.e-6%(1 + 3.60e-2*pow(ions.electron_temperature_scgc, 0.5))
                  % pow(ions.electron_temperature_scgc, 0.5) *neutrals.species[inO2].mass
                  /(cME + neutrals.species[inO2].mass) 
              + 8.90e-11*neutrals.species[inO].density_scgc*1.e-6%(1 + 5.70e-4*ions.electron_temperature_scgc)
                  %pow(ions.electron_temperature_scgc,0.5)*neutrals.species[inO2].mass/(cME + neutrals.species[inO2].mass) 
            );


  Qencm = Qencp % neutrals.temperature_scgc;

  report.exit(function);

  return std::vector<arma_cube> {Qencp, Qencm, Qenc_v};
  }