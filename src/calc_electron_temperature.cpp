// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include "../include/aether.h"


// --------------------------------------------------------------------------
// TODO (#24): this currently just sets the electron temperature to the neutral temperature
// --------------------------------------------------------------------------

void Ions::calc_electron_temperature(Neutrals neutrals, Grid &grid, Times time) {

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

  // Inelastic electron-neutral collisions:
  std::vector<arma_cube> Qenc_inelastic;
  arma_cube Qrotm, Qrotp, Qf, Qexc, Qvib_O2, Qvib_N2;

  // Thermoelectric Current
  arma_mat JParallel;

  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();
  int64_t nGCs = grid.get_nGCs();

  // Initialize everything to zero!
  epsilon.set_size(nLons, nLats, nAlts);
  epsilon.zeros();
  Qphe = epsilon;
  QIonization = epsilon;
  Qeicm = epsilon;
  Qeicp = epsilon;
  Qeic_v = epsilon;
  Qencm = epsilon;
  Qencp = epsilon;
  Qenc_v = epsilon;
  Qrotm = epsilon;
  Qrotp = epsilon;
  Qf = epsilon;
  Qexc = epsilon;
  Qvib_O2 = epsilon;
  Qvib_N2 = epsilon;

  report.print(4, "Calculating epsilon");

  // ============= ADD SOURCES =================


  // Needed for both ionization & photoelectron heating:
  if (input.get_do_ionization_heating() || input.get_do_photoelectron_heating())
    epsilon = calc_epsilon(neutrals);


  // Photoelectron heating
  if (input.get_do_photoelectron_heating()) {
    report.print(4, "Calculating photoelectron heating");
    Qphe = calc_photoelectron_heating(epsilon);
  }

  // Ionization heating (includes all ionization sources)
  if (input.get_do_ionization_heating()) {
    report.print(4, "Calculating ionization heating");
    QIonization = calc_ionization_heating(epsilon);
  }


  // electron-ion collisions
  if (input.get_do_electron_ion_collisional_heating()) {
    report.print(4, "Calculating electron-ion collisions");

    Qeic = calc_electron_ion_collisions();
    Qeicp = Qeic[0];
    Qeicm = Qeic[1];
    Qeic_v = Qeic[2]; // Friction
  }


  // electron-neutral Elastic collisions
  if (input.get_do_electron_neutral_elastic_collisional_heating()) {
    report.print(4, "Calculating electron-neutral elastic collisions");

    Qenc = calc_electron_neutral_elastic_collisions(neutrals);
    Qencp = Qenc[0];
    Qencm = Qenc[1];
    Qenc_v = Qenc[2]; // Friction
  }


  // electron-neutral inelastic collisions
  if (input.get_do_electron_neutral_inelastic_collisional_heating()) {
    report.print(4, "Calculating electron-neutral inelastic collisions");

    Qenc_inelastic = calc_electron_neutral_inelastic_collisions(neutrals);
    Qrotm = Qenc_inelastic[0];
    Qrotp = Qenc_inelastic[1];
    Qf = Qenc_inelastic[2];
    Qexc = Qenc_inelastic[3];
    Qvib_O2 = Qenc_inelastic[4];
    Qvib_N2 = Qenc_inelastic[5];
  }

  if (input.get_do_thermoelectric_heating()) {
    report.print(4, "Calculating thermoelectric heating");
    report.error("Thermoelectric heating is not working yet.");
    // Get JPar (Thermoelectric current)
    JParallel = calc_thermoelectric_current(grid);
  }

  report.print(4, "Calculating electron temperature");

  arma_cube sources = Qphe + Qencm + Qeicm + Qrotm + Qf + Qexc + Qvib_O2 + Qvib_N2
                      + QIonization + Qenc_v + Qeic_v;
  arma_cube sources_temp_dependent  = Qencp + Qeicp + Qrotp;


  // ============= CALCULATE TEMPERATURE =================
  arma_vec temp1d(nAlts);
  arma_vec lambda1d(nAlts);
  arma_vec front1d(nAlts);
  arma_vec dalt1d(nAlts);
  arma_vec conduction1d(nAlts);
  arma_vec sources1d(nAlts);
  arma_vec sources_tempdep_1d(nAlts);
  arma_vec ratios(nAlts);
  arma_vec density_ratio(nAlts);

  // Get the time step size
  precision_t dt = time.get_dt();


  for (int64_t iLon = nGCs; iLon < nLons - nGCs; iLon++) {
    for (int64_t iLat = nGCs; iLat < nLats - nGCs; iLat++) {
      temp1d = electron_temperature_scgc.tube(iLon, iLat);
      temp1d(0) = neutrals.temperature_scgc(iLon, iLat, 0);
      temp1d(1) = neutrals.temperature_scgc(iLon, iLat, 1);
      lambda1d = lambda.tube(iLon, iLat);
      lambda1d(1) = lambda1d(2);
      lambda1d(0) = lambda1d(2);
      front1d  = 3.0 / 2.0 * cKB * species[nSpecies].density_scgc.tube(iLon, iLat);
      dalt1d   = grid.dalt_lower_scgc.tube(iLon, iLat);
      sources1d = sources.tube(iLon, iLat);
      sources1d = sources1d / front1d;

      sources_tempdep_1d = sources_temp_dependent.tube(iLon, iLat);
      sources_tempdep_1d = sources_tempdep_1d / front1d;

      conduction1d.zeros();    // reset temp variable to zero
      conduction1d = solver_conduction(temp1d,
                                       lambda1d,
                                       front1d,
                                       sources1d,
                                       dalt1d,
                                       dt / 10.,
                                       nGCs,
                                       false,
                                       sources_tempdep_1d);

      // The conduction solver gives Tnew-Told, so divide by dt
      conduction1d.clamp(200, 6000);
      electron_temperature_scgc.tube(iLon, iLat) = conduction1d;
    }
  }

  // electron_temperature_scgc = neutrals.temperature_scgc;

  report.exit(function);
}



// Since this is used a few times, calculate it separately & pass it to the functions.
// From (Smithro and Solomon, 2008)
arma_cube Ions::calc_epsilon(Neutrals &neutrals) {

  std::string function = "Ions::calc_epsilon";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // from Swartz & Nisbet (1972): we need neutral species id's for o2, n2 & o
  int64_t inO2 = neutrals.get_species_id("O2");
  int64_t inN2 = neutrals.get_species_id("N2");
  int64_t inO = neutrals.get_species_id("O");

  if ((inO == -1) || (inN2 == -1) || (inO2 == -1))
    report.error("Could not find O, N2, or O2 in neutrals species list");

  arma_cube epsilon, x, logx;

  x = density_scgc / (neutrals.species[inO2].density_scgc
                      + neutrals.species[inN2].density_scgc
                      + neutrals.species[inO].density_scgc);
  // should rarely need to be used:
  x.clamp(2e-8, 10.0);

  logx = log(x);

  epsilon = exp(5.342 + 1.056 * logx - pow(4.392e-2 * logx, 2)
                - pow(5.9e-2 * logx, 3) - 9.346e-3 * pow(logx, 4)
                - 5.755e-4 * pow(logx, 5) - 1.249e-5 * pow(logx, 6)
               ) * 1.6e-19;

  report.exit(function);
  return epsilon;
}



// --------------------------------------------------------------------------
// Calculate photoelectron heating
// --------------------------------------------------------------------------
arma_cube Ions::calc_photoelectron_heating(arma_cube epsilon) {

  std::string function = "Ions::calc_photoelectron_heating";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t nIons = nSpecies;

  // Initialize Qphe & IonsIonizationRate (sum of ionization rates)
  // to the same size as epsilon and then zero them out:
  // (Qphe is the product of epsilon & IonsIonizationRate, so is not zeroed.)
  arma_cube Qphe, IonsIonizationRate;
  IonsIonizationRate.set_size(epsilon.n_rows, epsilon.n_cols, epsilon.n_slices);
  IonsIonizationRate.zeros();

  for (int64_t iIon = 0; iIon < nIons; iIon++)
    IonsIonizationRate += species[iIon].ionization_scgc;

  Qphe = epsilon % IonsIonizationRate;

  report.exit(function);
  return Qphe;
}


// --------------------------------------------------------------------------
// Calculate ionization heating
// --------------------------------------------------------------------------
arma_cube Ions::calc_ionization_heating(arma_cube epsilon) {

  std::string function = "Ions::calc_ionization_heating";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t nIons = nSpecies;

  // auroral heating efficiency coefficient
  precision_t auroheat = 1.0;

  int64_t iO_3P = get_species_id("O+");
  int64_t iO2P = get_species_id("O2+");
  int64_t iN2P = get_species_id("N2+");

  arma_cube QIonization = auroheat * epsilon % (species[iO_3P].ionization_scgc
                                                + species[iO2P].ionization_scgc
                                                + species[iN2P].ionization_scgc);


  report.exit(function);
  return QIonization;
}

// --------------------------------------------------------------------------
// Calculate electron-ion collisions
// --------------------------------------------------------------------------
std::vector<arma_cube> Ions::calc_electron_ion_collisions() {

  std::string function = "Ions::calc_electron_ion_collisions";
  static int iFunction = -1;
  report.enter(function, iFunction);

  arma_cube Qeicp;
  Qeicp.set_size(density_scgc.n_rows, density_scgc.n_cols, density_scgc.n_slices);
  Qeicp.zeros();
  arma_cube Qeicm = Qeicp;
  // Friction things
  arma_cube Qeic_v = Qeicp, dv2_ei = Qeicp;

  int64_t nSpecies = nSpecies;

  if (input.get_do_calc_bulk_ion_temp()) {
    // This is used when we calculate bulk ion temperature!
    report.print(3, "Using bulk ion temperature for electron-ion collisions");

    // Use all species, not just major species (different from GITM)
    for (int64_t iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      Qeicp += species[iSpecies].density_scgc
               / (cME + species[iSpecies].mass);
    }

    Qeicp = Qeicp % density_scgc * cME * 3.0 * cKB
            % (temperature_scgc - electron_temperature_scgc)
            * 5.45e-5 / pow(electron_temperature_scgc, 1.5);
  } else {
    // Individual ion temperatures:
    report.print(3,
                 "Using individual ion temperatures for electron-ion collisions");

    // Use all species, not just major species (different from GITM)
    for (int64_t iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      Qeicp += species[iSpecies].density_scgc
               % (species[iSpecies].temperature_scgc - electron_temperature_scgc)
               / (cME + species[iSpecies].mass);
    }

    Qeicp = Qeicp % density_scgc * cME * 3.0 * cKB
            * 5.45e-5 / pow(electron_temperature_scgc, 1.5);
  }

  report.print(3, "Calculating frictional heating");

  // delta velocity **2 btwn e- & ions:
  // This uses the bulk ion velocity, not the individual ion velocity.
  // (Different from GITM): Uses all species' densities (so just ne), not just o+, o2+, n2+, no+, n+
  for (int64_t iDir = 0; iDir < 3; iDir++)
    dv2_ei += pow(velocity_vcgc[iDir] - exb_vcgc[iDir], 2);

  Qeic_v = density_scgc * cME % dv2_ei * 5.45e-5 / pow(electron_temperature_scgc,
                                                       1.5)
           % (density_scgc);

  std::vector<arma_cube> Qeic = {Qeicp,
                                 Qeicp % electron_temperature_scgc,
                                 Qeic_v
                                };

  report.exit(function);

  return Qeic;
}

// --------------------------------------------------------------------------
// Calculate electron-neutral elastic collisions
// --------------------------------------------------------------------------
std::vector<arma_cube> Ions::calc_electron_neutral_elastic_collisions(
  Neutrals &neutrals) {

  std::string function = "Ions::calc_electron_neutral_elastic_collisions";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // initialize & zero the quantities we need:
  arma_cube Qenc;
  Qenc.set_size(density_scgc.n_rows, density_scgc.n_cols, density_scgc.n_slices);
  Qenc.zeros();
  arma_cube Qencp = Qenc;
  arma_cube Qencm = Qenc;
  // frictional things:
  arma_cube dv2_en = Qenc;
  arma_cube Qenc_v = Qenc;

  int64_t inO2 = neutrals.get_species_id("O2");
  int64_t inN2 = neutrals.get_species_id("N2");
  int64_t inO = neutrals.get_species_id("O");

  if ((inO == -1) || (inN2 == -1) || (inO2 == -1))
    report.error("Could not find O, N2, or O2 in neutrals species list");

  Qenc = density_scgc * cME * 3.0 * cKB % (neutrals.temperature_scgc -
                                           electron_temperature_scgc)
         % ((2.33e-11 * neutrals.species[inN2].density_scgc * 1.e-6
             % (1 - 1.21e-4 * electron_temperature_scgc)
             % electron_temperature_scgc / (cME + neutrals.species[inN2].mass))
            + (1.82e-10 * neutrals.species[inO2].density_scgc * 1.e-6
               % (1 + 3.60e-2 * pow(electron_temperature_scgc, 0.5))
               % pow(electron_temperature_scgc, 0.5) / (cME + neutrals.species[inO2].mass))
            + (8.90e-11 * neutrals.species[inO].density_scgc * 1.e-6
               % (1 + 5.70e-4 * electron_temperature_scgc)
               % pow(electron_temperature_scgc, 0.5) / (cME + neutrals.species[inO].mass))
           );
  report.print(6, "Qenc done");
  Qencp = density_scgc * cME * 3.0 * cKB
          % ((2.33e-11 * neutrals.species[inN2].density_scgc * 1.e-6
              % (1 - 1.21e-4 * electron_temperature_scgc)
              % electron_temperature_scgc / (cME + neutrals.species[inN2].mass))
             + (1.82e-10 * neutrals.species[inO2].density_scgc * 1.e-6
                % (1 + 3.60e-2 * pow(electron_temperature_scgc, 0.5))
                % pow(electron_temperature_scgc, 0.5) / (cME + neutrals.species[inO2].mass))
             + (8.90e-11 * neutrals.species[inO].density_scgc * 1.e-6
                % (1 + 5.70e-4 * electron_temperature_scgc)
                % pow(electron_temperature_scgc, 0.5) / (cME + neutrals.species[inO].mass))
            );
  report.print(6, "Qencp done");

  // delta velocity **2 btwn e- & neutrals:
  for (int64_t iDir = 0; iDir < 3; iDir++)
    dv2_en += pow(neutrals.velocity_vcgc[iDir] - exb_vcgc[iDir], 2);

  report.print(6, "dv2 done");

  Qenc_v = density_scgc * cME % dv2_en
           % (2.33e-11 * neutrals.species[inN2].density_scgc * 1.e-6
              % (1 - 1.21e-4 * electron_temperature_scgc) % electron_temperature_scgc *
              neutrals.species[inN2].mass
              / (cME + neutrals.species[inN2].mass)
              + 1.82e-10 * neutrals.species[inO2].density_scgc * 1.e-6 % (1 + 3.60e-2 * pow(
                    electron_temperature_scgc, 0.5))
              % pow(electron_temperature_scgc, 0.5) * neutrals.species[inO2].mass
              / (cME + neutrals.species[inO2].mass)
              + 8.90e-11 * neutrals.species[inO].density_scgc * 1.e-6 %
              (1 + 5.70e-4 * electron_temperature_scgc)
              % pow(electron_temperature_scgc,
                    0.5) * neutrals.species[inO2].mass / (cME + neutrals.species[inO2].mass)
             );

  report.print(6, "Qencv done");

  Qencm = Qencp % neutrals.temperature_scgc;

  report.exit(function);

  return std::vector<arma_cube> {Qencp, Qencm, Qenc_v};
}

// --------------------------------------------------------------------------
// Calculate electron-neutral inelasticcollisions
// --------------------------------------------------------------------------
std::vector<arma_cube> Ions::calc_electron_neutral_inelastic_collisions(
  Neutrals &neutrals) {

  std::string function = "Ions::calc_electron_neutral_inelastic_collisions";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // initialize & zero the quantities we need:
  arma_cube Qrot;
  Qrot.set_size(density_scgc.n_rows, density_scgc.n_cols, density_scgc.n_slices);
  Qrot.zeros(); // N2, O2 roration (Shunk & Nagy pp. 277)
  arma_cube Qrotp = Qrot;
  arma_cube Qrotm = Qrot;
  arma_cube Qf = Qrot; //fine structure heating rate (Shunk & Nagy pp. 282)
  arma_cube Qfp = Qrot,  Qfm = Qrot;
  arma_cube Qexc = Qrot; // O(1D) exitation
  arma_cube Qexcp = Qrot, Qexcm = Qrot;
  arma_cube Qvib_O2 = Qrot; // O2 vibration
  arma_cube logQ = Qrot, Qvib_O2p = Qrot, Qvib_O2m = Qrot;
  arma_cube Qvib_N2 = Qrot; // N2 vibration (Pavlov, 1998a)

  int64_t inO2 = neutrals.get_species_id("O2");
  int64_t inN2 = neutrals.get_species_id("N2");
  int64_t inO = neutrals.get_species_id("O");

  if ((inO == -1) || (inN2 == -1) || (inO2 == -1))
    report.error("Could not find O, N2, or O2 in neutrals species list");

  // some aliases for common quantities
  arma_cube ne = density_scgc;
  arma_cube Te = electron_temperature_scgc;
  arma_cube Ti = temperature_scgc;
  arma_cube Tn = neutrals.temperature_scgc;
  arma_cube no2 = neutrals.species[inO2].density_scgc;
  arma_cube nn2 = neutrals.species[inN2].density_scgc;
  arma_cube no = neutrals.species[inO].density_scgc;

  // Make sure (Ti && Te) > Tn everywhere (from GITM):
  arma::uvec Ti_mask = find(Ti <= Tn);
  arma::uvec Te_mask = find(Te <= Tn);
  Ti.elem(Ti_mask) = Tn.elem(Ti_mask) * 1.0001;
  Te.elem(Te_mask) = Tn.elem(Te_mask) * 1.0001;

  // N2, O2 rotation (Shunk & Nagy pp. 277)
  Qrot = 3.5e-14 * ne * 1.e-6 % nn2 * 1.e-6 % (Tn - Te) / (pow(Te, 0.5))
         + 5.2e-15 * ne * 1.e-6 % no2 * 1.e-6 % (Tn - Te) / (pow(Te, 0.5));
  Qrot = Qrot * 1.6e-13;   // eV cm-3 s -> J m-3

  Qrotp = 3.5e-14 * ne * 1.e-6 % nn2 * 1.e-6 / (pow(Te, 0.5)) // N2
          + 5.2e-15 * ne * 1.e-6 % no2 * 1.e-6 / (pow(Te, 0.5)); // O2
  Qrotp = Qrotp * 1.6e-13;
  Qrotm = Qrotp % Tn;

  // fine strcture heating rate  by Shunk and Nagy Page 282
  Qf = 0.0;
  arma_cube Dfine = 5.0 + exp(-326.6 / Tn) + 3.0 * exp(-227.7 / Tn);
  arma_cube s10 = 8.249e-16 * pow(Te, 0.6) % exp(-227.7 / Tn);
  Qf = ne % no * 1.e-12 / Dfine % (s10 % (1. - exp(98.9 * (1 / Te - 1 / Tn)))
                                   + 1.191e-11 * (1. - exp(326.6 * (1 / Te - 1 / Tn)))
                                   + 1.863e-11 * (1. - exp(227.7 * (1 / Te - 1 / Tn))));
  Qf = -Qf * 1.6e-13;      // eV cm-3 s-1 -> J m-3 s-1

  Qfp = 0.0;
  Qfm = Qf;

  // O(1D) excitation
  arma_cube te_exc = Te.clamp(0.0, 18000.0);
  arma_cube dexc = 2.4e4 + 0.3 * (te_exc - 1500.) - 1.947e-5 *
                   (te_exc - 1500.) % (te_exc - 4000.);
  Qexc = 1.57e-12 * ne % no * 1.e-12 % exp(dexc % (te_exc - 3000.) / 3000. /
                                           te_exc)
         % (exp(-22713 * (te_exc - Tn) / te_exc / Tn) - 1);
  Qexc = Qexc * 1.6e-13;      // eV cm-3 s-1 -> J m-3 s-1

  Qexcp = 0.0;
  Qexcm = Qexc;

  // O2 vibration:

  // Calculated differently from GITM, but the result is the same...
  // GITM clamped e- temp before calculating, but that makles things really hard here.
  // Instead we will use all Te's, and then limit the outputs after.

  logQ = (5.0148e-31 * pow(Te, 9) - 1.5346e-26 * pow(Te, 8)
          + 2.0127e-22 * pow(Te, 7) - 1.4791e-18 * pow(Te, 6)
          + 6.6865e-15 * pow(Te, 5) - 1.9228e-11 * pow(Te, 4)
          + 3.5187e-8 * pow(Te, 3) - 3.996e-5 * pow(Te, 2)
          + 0.0267 * Te - 19.9171);
  // GITM's Te_6000 was used when 300<T<6000, which corresponds to ~-15.9 & 198.3 for logQ
  // Mask the values outside of this range...
  logQ.clamp(-15.9, 198.273);
  logQ.elem( find(logQ <= -15.9) ).fill(-20.0);


  Qvib_O2 = ne % no2 * 1.e-12 % exp10(logQ) % (1 - exp(2239. *
                                                       (1. / Te - 1. / Tn)));
  Qvib_O2 = -Qvib_O2 * 1.6e-13;
  Qvib_O2p = 0.;
  Qvib_O2m = Qvib_O2;


  // N2 vibration from Pavlov 1998a

  // we'll use the same loop GITM uses:
  arma_vec Av0 = {2.025, -7.066, -8.211, -9.713, -10.353, -10.819, -10.183, -12.698, -14.710, -17.538};
  arma_vec Bv0 = {8.782e-4, 1.001e-2, 1.092e-2, 1.204e-2, 1.243e-2, 1.244e-2, 1.185e-2, 1.309e-2, 1.409e-2, 1.6e-2};
  arma_vec Cv0 = {2.954e-7, -3.066e-6, -3.369e-6, -3.732e-6, -3.850e-6, -3.771e-6, -3.570e-6, -3.952e-6, -4.249e-6, -4.916e-6};
  arma_vec Dv0 = {-9.562e-11, 4.436e-10, 4.891e-10, 5.431e-10, 5.6e-10, 5.385e-10, 5.086e-10, 5.636e-10, 6.058e-10, 7.128e-10};
  arma_vec Fv0 = {7.252e-15, -2.449e-14, -2.706e-14, -3.008e-14, -3.1e-14, -2.936e-14, -2.769e-14, -3.071e-14, -3.3e-14, -3.941e-14};

  precision_t Av0L = -6.462;
  precision_t Bv0L = 3.151e-2;
  precision_t Cv0L = -4.075e-5;
  precision_t Dv0L = 2.439e-8;
  precision_t Fv0L = -5.479e-12;

  arma_vec Av1 = {-3.413, -4.16, -5.193, -5.939, -8.261, -8.185, -10.823, -11.273};
  arma_vec Bv1 = {7.326e-3, 7.803e-3, 8.36e-3, 8.807e-3, 1.01e-2, 1.01e-2, 1.199e-2, 1.283e-2};
  arma_vec Cv1 = {-2.2e-6, -2.352e-6, -2.526e-6, -2.669e-6, -3.039e-6, -3.039e-6, -3.62e-6, -3.879e-6};
  arma_vec Dv1 = {3.128e-10, 3.352e-10, 3.606e-10, 3.806e-10, 4.318e-10, 4.318e-10, 5.159e-10, 5.534e-10};
  arma_vec Fv1 = {-1.702e-14, -1.828e-14, -1.968e-14, -2.073e-14, -2.347e-14, -2.347e-14, -2.81e-14, -3.016e-14};

  arma_vec logQv0, logQv1;
  logQv0.set_size(10);
  logQv1.set_size(8);

  double tte; // since we need the min of this & 6000 & std::min can't accept precision_t
  precision_t ttn, tte_6000;

  for (int64_t iLon = 0; iLon < density_scgc.n_rows; iLon++) { //nLons
    for (int64_t iLat = 0; iLat < density_scgc.n_cols; iLat++) { // nLats
      for (int64_t iAlt = 0; iAlt < density_scgc.n_slices; iAlt++) { // nAlts
        tte = Te(iLon, iLat, iAlt);
        ttn = Tn(iLon, iLat, iAlt);
        tte_6000 = std::min(tte, 6000.0);
        logQv0 = -20.0;
        logQv1 = -20.0;

        if (tte <= 1500 && tte > 300) {
          logQv0(0) = Av0L + Bv0L * tte + Cv0L * pow(tte, 2) + Dv0L * pow(tte,
                      3) + Fv0L * pow(tte, 4) - 16.0;
          Qvib_N2(iLon, iLat, iAlt) = (1.0 - exp(-3353.0 / ttn)) * exp10(logQv0(0))
                                      * (1.0 - exp(3353.0 * (1.0 / tte - 1.0 / ttn)));
        } else if (tte > 1500) {
          logQv0 = Av0 + Bv0 * tte + Cv0 * pow(tte, 2) + Dv0 * pow(tte,
                                                                   3) + Fv0 * pow(tte, 4) - 16.0;
          logQv1 = Av1 + Bv1 * tte + Cv1 * pow(tte, 2) + Dv1 * pow(tte,
                                                                   3) + Fv1 * pow(tte, 4) - 16.0;

          for (int iLevel = 0; iLevel < 10; iLevel++) {
            Qvib_N2(iLon, iLat, iAlt) += (1.0 - exp(-3353.0 / ttn)) * exp10(logQv0(iLevel))
                                         * (1.0 - exp((iLevel + 1) * 3353.0 * (1.0 / tte - 1.0 / ttn)));
          }

          for (int iLevel = 0; iLevel < 8; iLevel++) {
            // GITM did this loop from 2-9, which is out of bounds, so we correct for iLevel here)
            Qvib_N2(iLon, iLat, iAlt) += (1.0 - exp(-3353.0 / ttn)) * exp(-3353.0 / ttn)
                                         * exp10(logQv1(iLevel))
                                         * (1.0 - exp((iLevel + 2) * 3353.0 * (1.0 / tte - 1.0 / ttn)));
          }
        }
      }
    }
  }

  Qvib_N2 = -ne % nn2 * 1.e-12 % Qvib_N2 * 1.6e-13;

  report.exit(function);
  report.print(6, "Qenc done");

  return std::vector<arma_cube> {Qrotm, Qrotp, Qf, Qexc, Qvib_O2, Qvib_N2};
}

arma_mat Ions::calc_thermoelectric_current(Grid &grid) {

  std::string function = "Ions::calc_thermoelectric_current";
  static int iFunction = -1;
  report.enter(function, iFunction);

  std::vector <arma_cube> JTotal;

  arma_mat JParallel;
  JParallel.set_size(density_scgc.n_rows, density_scgc.n_cols);
  JParallel.zeros();

  // with the dipole, the field-aligned current is in the k^ direction
  // But we do not solve for e- velocity (and exb is 0 parallel to B), so we cannot do this:
  // if (grid.iGridShape_ == iDipole_){
  //   for (int64_t iAlt = 0; iAlt < ions.density_scgc.n_slices; iAlt++){
  //     JParallel += (ions.density_scgc.slice(iAlt) * cE % (ions.velocity_vcgc[2].slice(iAlt) - ions.exb_vcgc[2].slice(iAlt)))
  //                  * grid.dalt_center_scgc[iAlt];
  //   }
  // }

  // else{
  // We need to solve for J_Parallel with Equation (6) of (Zhu et al., 2016)
  // http://dx.doi.org/10.1016/j.jastp.2016.01.005 - using \nabla \dot  J
  JTotal.push_back(cE * density_scgc % (velocity_vcgc[0] - exb_vcgc[0]));
  JTotal.push_back(cE * density_scgc % (velocity_vcgc[1] - exb_vcgc[1]));
  JTotal.push_back(cE * density_scgc % (velocity_vcgc[2] - exb_vcgc[2]));

  arma_cube JuTotalDotB = dot_product(JTotal, grid.bfield_vcgc);
  arma_cube divJperp;

  divJperp.set_size(density_scgc.n_rows, density_scgc.n_cols,
                    density_scgc.n_slices);
  divJperp.zeros();

  for (int64_t iDir = 0; iDir < 3; iDir ++) {
    divJperp += calc_gradient_vector(JTotal[iDir], grid)[iDir];
    divJperp -= calc_gradient_vector(JuTotalDotB,
                                     grid)[iDir] % grid.bfield_unit_vcgc[iDir];
    divJperp -= calc_gradient_vector(JuTotalDotB, grid)[iDir] % JuTotalDotB;
  }

  for (int64_t iAlt = 0; iAlt < density_scgc.n_slices; iAlt++)
    JParallel -= divJperp.slice(iAlt) % grid.dalt_center_scgc.slice(iAlt);

  report.exit(function);
  return JParallel;
}