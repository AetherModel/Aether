// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <cmath>
#include <iostream>

#include "../include/constants.h"
#include "../include/neutrals.h"
#include "../include/earth.h"
#include "../include/report.h"
#include "../include/time.h"
#include "../include/solvers.h"

//----------------------------------------------------------------------
//
//----------------------------------------------------------------------

void Neutrals::calc_mass_density(Report &report) {

  int64_t iLon, iLat, iAlt, index, iSpecies;

  std::string function = "Neutrals::calc_mass_density";
  static int iFunction = -1;
  report.enter(function, iFunction);

  rho_scgc.zeros();
  density_scgc.zeros();

  for (iSpecies=0; iSpecies < nSpecies; iSpecies++) {
    rho_scgc = rho_scgc +
      neutrals[iSpecies].mass * neutrals[iSpecies].density_scgc;
    density_scgc = density_scgc + neutrals[iSpecies].density_scgc;
  }

  mean_major_mass_scgc = rho_scgc / density_scgc;
  pressure_scgc = boltzmanns_constant * density_scgc % temperature_scgc;

  report.exit(function);
  return;
}

//----------------------------------------------------------------------
//
//----------------------------------------------------------------------

void Neutrals::calc_specific_heat(Report &report) {

  int64_t iLon, iLat, iAlt, index, iSpecies;

  double t, p, r;

  std::string function = "Neutrals::calc_specific_heat";
  static int iFunction = -1;
  report.enter(function, iFunction);

  Cv_scgc.zeros();
  gamma_scgc.zeros();
  kappa_scgc.zeros();

  for (iSpecies=0; iSpecies < nSpecies; iSpecies++) {
    Cv_scgc = Cv_scgc +
      (neutrals[iSpecies].vibe - 2) *
      neutrals[iSpecies].density_scgc *
      boltzmanns_constant / neutrals[iSpecies].mass;
    gamma_scgc = gamma_scgc +
      neutrals[iSpecies].density_scgc / (neutrals[iSpecies].vibe-2);
    kappa_scgc = kappa_scgc +
      neutrals[iSpecies].thermal_cond *
      neutrals[iSpecies].density_scgc %
      pow(temperature_scgc, neutrals[iSpecies].thermal_exp);
  }

  Cv_scgc = Cv_scgc / (2*density_scgc);
  gamma_scgc = gamma_scgc * 2.0 / density_scgc + 1.0;
  kappa_scgc = kappa_scgc / density_scgc;

  sound_scgc = sqrt(boltzmanns_constant *
                    gamma_scgc %
                    temperature_scgc /
                    mean_major_mass_scgc);

  report.exit(function);
  return;
}

//----------------------------------------------------------------------
// Calculate the altitude integral of the different species for EUV
// calculations.  Scale these to be the slant path through the
// atmosphere.  This gets complicated behind the terminator.  All of
// this is taken from Smith and Smith, JGR 1972, vol. 77, page 3592
// ----------------------------------------------------------------------

void Neutrals::calc_chapman(Grid grid, Report &report) {

  int64_t iAlt, iLon, iLat, index, indexp;
  float H;  // scale height

  // This is all from Smith and Smith, JGR 1972, vol. 77, page 3592
  // "Numerical evaluation of chapman's grazing incidence integral ch(X,x)"
  // Xp is supposed to be R/H
  // JMB Update: 05/2017.  Corrected a small error in the y-calc for
  // erfc(y)
  //
  // Also Updated the Grazing Integral for SZA > 90.0
  // We now do log-linear interpolation for smoother transitions

  double a = 1.06069630;
  double b = 0.55643831;
  double c = 1.06198960;
  double d = 1.72456090;
  double f = 0.56498823;
  double g = 0.06651874;

  double y, dy;

  float Hp_up, Hp_dn, grad_hs, grad_xp, grad_in, Hg, Xg, in, int_g, int_p;
  int64_t index_bottom, iindex, iindexp, iiAlt;

  std::string function = "Neutrals::calc_chapman";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();
  int64_t nGCs = grid.get_nGCs();

  // New way of doing it with 3D arrays:

  fcube integral3d(nLons, nLats, nAlts);
  fcube log_int3d(nLons, nLats, nAlts);
  fcube xp3d(nLons, nLats, nAlts);
  fcube y3d(nLons, nLats, nAlts);
  fcube erfcy3d(nLons, nLats, nAlts);

  fvec integral1d(nAlts);
  fvec log_int1d(nAlts);
  fvec xp1d(nAlts);
  fvec y1d(nAlts);
  fvec erfcy1d(nAlts);
  fvec dAlt1d(nAlts);
  fvec sza1d(nAlts);
  fvec radius1d(nAlts);
  fvec H1d(nAlts);

  for (int iSpecies=0; iSpecies < nSpecies; iSpecies++) {

    neutrals[iSpecies].scale_height_scgc =
      boltzmanns_constant * temperature_scgc /
      (neutrals[iSpecies].mass * grid.gravity_scgc);

    xp3d = grid.radius_scgc / neutrals[iSpecies].scale_height_scgc;
    y3d = sqrt(0.5 * xp3d) % abs(grid.cos_sza_scgc);
    iAlt = nAlts-1;

    integral3d.slice(iAlt) =
      neutrals[iSpecies].density_scgc.slice(iAlt) %
      neutrals[iSpecies].scale_height_scgc.slice(iAlt);

    for (iAlt = nAlts-1; iAlt >= 0; iAlt--) {
      if (iAlt < nAlts-1) {
        integral3d.slice(iAlt) = integral3d.slice(iAlt+1) +
          neutrals[iSpecies].density_scgc.slice(iAlt) %
          grid.dalt_lower_scgc.slice(iAlt+1);
      }
    }

    erfcy3d = (a + b * y3d) / (c + d*y3d + y3d % y3d);
    for (iLon = 0; iLon < nLons ; iLon++)
      for (iLat = 0; iLat < nLats ; iLat++)
        for (iAlt = 0; iAlt < nAlts ; iAlt++)
          if (y3d(iLon, iLat, iAlt) >= 8.0)
            erfcy3d(iLon, iLat, iAlt) = f / (g + y3d(iLon, iLat, iAlt));

    log_int3d = log(integral3d);

    // Set chapman integrals to max in the lower ghostcells

    for (int iAlt=0; iAlt < nGCs; iAlt++)
      neutrals[iSpecies].chapman_scgc.slice(iAlt).fill(max_chapman);

    for (iLon = 0; iLon < nLons ; iLon++) {
      for (iLat = 0; iLat < nLats ; iLat++) {

        dAlt1d = grid.dalt_lower_scgc.tube(iLon, iLat);
        sza1d = grid.sza_scgc.tube(iLon, iLat);
        integral1d = integral3d.tube(iLon, iLat);
        log_int1d = log_int3d.tube(iLon, iLat);
        xp1d = xp3d.tube(iLon, iLat);
        y1d = y3d.tube(iLon, iLat);
        erfcy1d = erfcy3d.tube(iLon, iLat);
        radius1d = grid.radius_scgc.tube(iLon, iLat);
        H1d = neutrals[iSpecies].scale_height_scgc.tube(iLon, iLat);

        for (iAlt = nGCs; iAlt < nAlts; iAlt++) {
          // This is on the dayside:
          if (sza1d(iAlt) < pi/2 || sza1d(iAlt) > 3*pi/2) {
            neutrals[iSpecies].chapman_scgc(iLon, iLat, iAlt) =
              integral1d(iAlt) * sqrt(0.5 * pi * xp1d(iAlt)) * erfcy1d(iAlt);
          } else {
            // This is on the nghtside of the terminator:

            y = radius1d(iAlt) * abs(cos(sza1d(iAlt)-pi/2));

            // This sort of assumes that nGeoGhosts >= 2:
            if (y > radius1d(nGCs)) {

              iiAlt = iAlt;
              while (radius1d(iiAlt-1) > y) iiAlt--;
              iiAlt--;

              Hp_up = H1d(iiAlt+1);
              Hp_dn = H1d(iiAlt);

              // make sure to use the proper cell spacing (iiAlt+1 & lower):
              grad_hs = (Hp_up - Hp_dn) / dAlt1d(iiAlt+1);
              grad_xp = (xp1d(iiAlt+1) - xp1d(iiAlt)) / dAlt1d(iiAlt+1);
              grad_in = (log_int1d(iiAlt+1) - log_int1d(iiAlt)) / dAlt1d(iiAlt+1);

              // Linearly interpolate H and X:
              dy = y - radius1d(iiAlt);
              Hg = Hp_dn + grad_hs * dy;
              Xg = xp1d(iiAlt) + grad_xp * dy;
              in = log_int1d(iiAlt) + grad_in * dy;

              int_g = exp(in);
              int_p = integral1d(iAlt);
              // Equation (19) Smith & Smith
              neutrals[iSpecies].chapman_scgc(iLon, iLat, iAlt) =
                sqrt(0.5 * pi * Xg) * (2.0 * int_g - int_p * erfcy1d(iAlt));

              if (neutrals[iSpecies].chapman_scgc(iLon, iLat, iAlt) > max_chapman)
                neutrals[iSpecies].chapman_scgc(iLon, iLat, iAlt) = max_chapman;

            } else {
              // This says that we are in the shadow of the planet:

              neutrals[iSpecies].chapman_scgc(iLon, iLat, iAlt) = max_chapman;
            }
          }
        }  // iAlt
      }  // iLat
    }  // iLon
  }  // iSpecies

  report.exit(function);
  return;
}

// -----------------------------------------------------------------------------
// Calculate thermal conduction
// -----------------------------------------------------------------------------

void Neutrals::calc_conduction(Grid grid, Times time, Report &report) {

  float dt;

  int64_t iLon, iLat, iAlt, index;

  std::string function = "Neutrals::calc_conduction";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();

  fcube rhocvr23d(nLons, nLats, nAlts);
  fcube lambda3d(nLons, nLats, nAlts);
  fcube prandtl3d(nLons, nLats, nAlts);

  rhocvr23d = rho_scgc % Cv_scgc % grid.radius2_scgc;
  // Need to make this eddy * rho * cv:
  prandtl3d.zeros();
  lambda3d = (kappa_scgc + prandtl3d) % grid.radius2_scgc;

  fvec temp1d(nAlts);
  fvec lambda1d(nAlts);
  fvec rhocvr21d(nAlts);
  fvec dalt1d(nAlts);
  fvec conduction1d(nAlts);

  for (iLon = 0; iLon < nLons; iLon++) {
    for (iLat = 0; iLat < nLats; iLat++) {

      temp1d = temperature_scgc.tube(iLon, iLat);
      lambda1d = lambda3d.tube(iLon, iLat);
      rhocvr21d = rhocvr23d.tube(iLon, iLat);
      dalt1d = grid.dalt_lower_scgc.tube(iLon, iLat);
      conduction1d.zeros();

      dt = time.get_dt();

      conduction1d = solver_conduction(temp1d, lambda1d, rhocvr21d, dt, dalt1d);

      // We want the sources to be in terms of dT/dt, while the
      // conduction actually solves for Tnew-Told, so divide by dt

      conduction_scgc.tube(iLon, iLat) = conduction1d / dt;
    }  // lat
  }  // lon
  report.exit(function);
}



// -----------------------------------------------------------------------------
// Calculate EUV driven ionization and heating rates
// -----------------------------------------------------------------------------

void Neutrals::calc_ionization_heating(Euv euv, Ions &ions, Report &report) {

  int64_t iAlt, iLon, iLat, iWave, iSpecies, index, indexp;
  int i_, idion_, ideuv_, nIonizations, iIon, iIonization;
  float tau, intensity, photoion;

  float ionization;

  std::string function = "calc_ionization_heating";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // Zero out all source terms:

  heating_euv_scgc.zeros();
  for (iSpecies=0; iSpecies < nSpecies; iSpecies++)
    neutrals[iSpecies].ionization_scgc.zeros();

  int64_t nLons = heating_euv_scgc.n_rows;
  int64_t nLats = heating_euv_scgc.n_cols;
  int64_t nAlts = heating_euv_scgc.n_slices;

  fmat tau2d = heating_euv_scgc.slice(0);
  fmat intensity2d = heating_euv_scgc.slice(0);
  fmat ionization2d = heating_euv_scgc.slice(0);

  for (iAlt = 2; iAlt < nAlts-2; iAlt++) {
    for (iWave=0; iWave < euv.nWavelengths; iWave++) {

      tau2d.zeros();

      for (iSpecies=0; iSpecies < nSpecies; iSpecies++) {
        if (neutrals[iSpecies].iEuvAbsId_ > -1) {
          i_ = neutrals[iSpecies].iEuvAbsId_;
          tau2d = tau2d +
            euv.waveinfo[i_].values[iWave] *
            neutrals[iSpecies].chapman_scgc.slice(iAlt);
        }
      }

      intensity2d = euv.wavelengths_intensity_top[iWave] * exp(-1.0*tau2d);

      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        // Calculate Photo-Absorbtion for each species and add them up:
        i_ = neutrals[iSpecies].iEuvAbsId_;  // index of photo abs cross section
        if (i_ > -1) {
          heating_euv_scgc.slice(iAlt) = heating_euv_scgc.slice(iAlt) +
            euv.wavelengths_energy[iWave] *
            euv.waveinfo[i_].values[iWave] *  // cross section
            (intensity2d %
             neutrals[iSpecies].density_scgc.slice(iAlt) );
        }

        for (iIonization = 0;
             iIonization < neutrals[iSpecies].nEuvIonSpecies;
             iIonization++) {

          i_ = neutrals[iSpecies].iEuvIonId_[iIonization];

          ionization2d =
            euv.waveinfo[i_].values[iWave] *  // cross section
            intensity2d %
            neutrals[iSpecies].density_scgc.slice(iAlt);

          neutrals[iSpecies].ionization_scgc.slice(iAlt) =
            neutrals[iSpecies].ionization_scgc(iAlt) + ionization2d;

          iIon = neutrals[iSpecies].iEuvIonSpecies_[iIonization];
          ions.species[iIon].ionization_scgc.slice(iAlt) =
            ions.species[iIon].ionization_scgc.slice(iAlt) + ionization2d;
        }  // iIonization
      }  // iSpecies
    }  // iWave
  }  // iAlt

  heating_euv_scgc =
    heating_efficiency * heating_euv_scgc / rho_scgc / Cv_scgc;

  report.exit(function);
  return;
}
